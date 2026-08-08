#include "viewer_window.h"

#include "vector_glyphs.h"
#include "viewer_protocol.h"

#include <QAction>
#include <QComboBox>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QKeySequence>
#include <QLabel>
#include <QMatrix4x4>
#include <QMessageBox>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPainter>
#include <QStatusBar>
#include <QTimer>
#include <QToolBar>
#include <QVector3D>
#include <QWheelEvent>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace simwing::viewer {
namespace {

QString fromUtf8(const std::string& value) {
    return QString::fromUtf8(value.data(), static_cast<qsizetype>(value.size()));
}

QString abbreviated(const std::string& value) {
    const QString text = fromUtf8(value);
    return text.size() > 34 ? text.left(31) + QStringLiteral("...") : text;
}

struct RenderVertex {
    float position[3]{};
    float normal[3]{};
    float colour[3]{};
};

QVector3D colourMap(double value, double minimum, double maximum) {
    double t = 0.5;
    if (maximum > minimum) {
        t = (value - minimum) / (maximum - minimum);
    }
    t = std::clamp(t, 0.0, 1.0);

    // A compact blue-cyan-yellow-red map with useful contrast on both dark
    // and light surfaces. It is diagnostic rather than a categorical scale.
    const double red = std::clamp(1.5 - std::abs(4.0 * t - 3.0), 0.0, 1.0);
    const double green = std::clamp(1.5 - std::abs(4.0 * t - 2.0), 0.0, 1.0);
    const double blue = std::clamp(1.5 - std::abs(4.0 * t - 1.0), 0.0, 1.0);
    return QVector3D(static_cast<float>(red), static_cast<float>(green),
                     static_cast<float>(blue));
}

QVector3D lineRoleColour(std::uint32_t role) {
    static const QVector3D colours[] = {
        {0.92F, 0.93F, 0.96F}, {0.98F, 0.72F, 0.22F},
        {0.28F, 0.80F, 0.95F}, {0.78F, 0.46F, 0.94F},
        {0.44F, 0.90F, 0.48F}, {0.96F, 0.42F, 0.40F}};
    return colours[role % (sizeof(colours) / sizeof(colours[0]))];
}

void appendVertex(
    std::vector<RenderVertex>& output,
    const Vec3d& position,
    const QVector3D& normal,
    const QVector3D& colour) {
    RenderVertex vertex;
    vertex.position[0] = static_cast<float>(position.x);
    vertex.position[1] = static_cast<float>(position.y);
    vertex.position[2] = static_cast<float>(position.z);
    vertex.normal[0] = normal.x();
    vertex.normal[1] = normal.y();
    vertex.normal[2] = normal.z();
    vertex.colour[0] = colour.x();
    vertex.colour[1] = colour.y();
    vertex.colour[2] = colour.z();
    output.push_back(vertex);
}

class TraceStream final {
public:
    using HeaderCallback = std::function<void(TraceHeader)>;
    using FrameCallback =
        std::function<void(std::shared_ptr<const DiagnosticFrame>)>;
    using ErrorCallback = std::function<void(QString)>;
    using EndCallback = std::function<void()>;

    TraceStream(
        QString fileName,
        bool follow,
        HeaderCallback headerCallback,
        FrameCallback frameCallback,
        ErrorCallback errorCallback,
        EndCallback endCallback)
        : fileName_(std::move(fileName)),
          follow_(follow),
          headerCallback_(std::move(headerCallback)),
          frameCallback_(std::move(frameCallback)),
          errorCallback_(std::move(errorCallback)),
          endCallback_(std::move(endCallback)),
          thread_([this] { run(); }) {}

    ~TraceStream() {
        {
            std::lock_guard lock(mutex_);
            stopping_ = true;
        }
        condition_.notify_all();
        if (thread_.joinable()) {
            thread_.join();
        }
    }

    TraceStream(const TraceStream&) = delete;
    TraceStream& operator=(const TraceStream&) = delete;

    void requestNext() {
        {
            std::lock_guard lock(mutex_);
            if (stopping_ || ended_) {
                return;
            }
            requested_ = true;
        }
        condition_.notify_one();
    }

private:
    void run() {
        const QString absoluteFileName = QFileInfo(fileName_).absoluteFilePath();
#ifdef _WIN32
        const std::filesystem::path nativePath(absoluteFileName.toStdWString());
#else
        const QByteArray encodedPath = absoluteFileName.toUtf8();
        const std::filesystem::path nativePath(encodedPath.constData());
#endif
        std::ifstream input(nativePath, std::ios::binary);
        if (!input) {
            errorCallback_(QStringLiteral("Cannot open trace file: %1")
                               .arg(QFileInfo(fileName_).absoluteFilePath()));
            return;
        }

        TraceReader reader(
            input, follow_ ? TraceReadMode::Follow : TraceReadMode::Replay);
        TraceHeader header;
        if (!reader.readHeader(header)) {
            errorCallback_(QStringLiteral("Cannot read trace header: %1")
                               .arg(fromUtf8(reader.error().message)));
            return;
        }
        headerCallback_(std::move(header));

        for (;;) {
            {
                std::unique_lock lock(mutex_);
                condition_.wait(lock, [this] { return stopping_ || requested_; });
                if (stopping_) {
                    return;
                }
                requested_ = false;
            }

            DiagnosticFrame frame;
            const TraceReadStatus status = reader.readNext(frame);
            if (status == TraceReadStatus::Frame) {
                frameCallback_(
                    std::make_shared<const DiagnosticFrame>(std::move(frame)));
                continue;
            }

            if (status == TraceReadStatus::Pending) {
                // Follow mode deliberately treats both clean and partial
                // physical EOF as temporary. A condition-variable timeout
                // bounds polling and lets shutdown wake the reader at once.
                std::unique_lock lock(mutex_);
                if (condition_.wait_for(
                        lock, std::chrono::milliseconds(100),
                        [this] { return stopping_; })) {
                    return;
                }
                continue;
            }

            {
                std::lock_guard lock(mutex_);
                ended_ = true;
            }
            if (status == TraceReadStatus::End) {
                endCallback_();
            } else {
                errorCallback_(QStringLiteral("Cannot read trace frame: %1")
                                   .arg(fromUtf8(reader.error().message)));
            }
            return;
        }
    }

    QString fileName_;
    bool follow_ = false;
    HeaderCallback headerCallback_;
    FrameCallback frameCallback_;
    ErrorCallback errorCallback_;
    EndCallback endCallback_;
    std::mutex mutex_;
    std::condition_variable condition_;
    bool stopping_ = false;
    bool requested_ = true;
    bool ended_ = false;
    std::thread thread_;
};

class TraceView final : public QOpenGLWidget, protected QOpenGLFunctions {
public:
    explicit TraceView(QWidget* parent = nullptr) : QOpenGLWidget(parent) {
        setFocusPolicy(Qt::StrongFocus);
        setMinimumSize(640, 420);
    }

    void setFrame(std::shared_ptr<const DiagnosticFrame> frame) {
        // Keep a moving simulation in the same camera-relative position. A
        // free-flying wing must not simply leave the viewport after the first
        // frame was fitted; manual orbit, pan, and zoom remain relative to the
        // tracked surface centre.
        if (frame_ != nullptr && frame != nullptr
            && !frame_->vertices.empty() && !frame->vertices.empty()) {
            target_ += frameCentre(*frame) - frameCentre(*frame_);
        }
        frame_ = std::move(frame);
        geometryDirty_ = true;
        update();
    }

    void setField(QString name) {
        fieldName_ = std::move(name);
        geometryDirty_ = true;
        update();
    }

    void setVectorField(QString name) {
        vectorFieldName_ = std::move(name);
        geometryDirty_ = true;
        update();
    }

    void setMessage(QString message) {
        message_ = std::move(message);
        update();
    }

    [[nodiscard]] QString renderError() const { return glError_; }

    void fitView() {
        if (frame_ == nullptr || frame_->vertices.empty()) {
            target_ = {};
            distance_ = 5.0F;
            update();
            return;
        }
        QVector3D minimum(
            std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max());
        QVector3D maximum(
            std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest());
        for (const DiagnosticVertex& vertex : frame_->vertices) {
            const QVector3D point(
                static_cast<float>(vertex.positionMetres.x),
                static_cast<float>(vertex.positionMetres.y),
                static_cast<float>(vertex.positionMetres.z));
            minimum.setX(std::min(minimum.x(), point.x()));
            minimum.setY(std::min(minimum.y(), point.y()));
            minimum.setZ(std::min(minimum.z(), point.z()));
            maximum.setX(std::max(maximum.x(), point.x()));
            maximum.setY(std::max(maximum.y(), point.y()));
            maximum.setZ(std::max(maximum.z(), point.z()));
        }
        target_ = 0.5F * (minimum + maximum);
        const float radius = std::max(0.5F * (maximum - minimum).length(), 0.05F);
        // Leave enough perspective margin for modest deformation after a live
        // fit; the theoretical bounding-sphere minimum has no visual gutter
        // and clips as soon as a structural case changes shape.
        distance_ = std::max(radius * 3.8F, 0.2F);
        update();
    }

protected:
    void initializeGL() override {
        initializeOpenGLFunctions();
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glClearColor(0.045F, 0.055F, 0.075F, 1.0F);

        const bool coreProfile =
            context()->format().profile() == QSurfaceFormat::CoreProfile;
        const QString vertexShader = coreProfile
            ? QStringLiteral(
                  "#version 330 core\n"
                  "in vec3 position;\n"
                  "in vec3 normal;\n"
                  "in vec3 colour;\n"
                  "uniform mat4 mvp;\n"
                  "out vec3 vNormal;\n"
                  "out vec3 vColour;\n"
                  "void main() {\n"
                  "  vNormal = normal; vColour = colour;\n"
                  "  gl_Position = mvp * vec4(position, 1.0);\n"
                  "}\n")
            : QStringLiteral(
                  "attribute vec3 position;\n"
                  "attribute vec3 normal;\n"
                  "attribute vec3 colour;\n"
                  "uniform mat4 mvp;\n"
                  "varying vec3 vNormal;\n"
                  "varying vec3 vColour;\n"
                  "void main() {\n"
                  "  vNormal = normal; vColour = colour;\n"
                  "  gl_Position = mvp * vec4(position, 1.0);\n"
                  "}\n");
        const QString fragmentShader = coreProfile
            ? QStringLiteral(
                  "#version 330 core\n"
                  "uniform bool lit;\n"
                  "in vec3 vNormal;\n"
                  "in vec3 vColour;\n"
                  "out vec4 fragColor;\n"
                  "void main() {\n"
                  "  float d = lit ? 0.32 + 0.68 * abs(dot(normalize(vNormal), "
                  "normalize(vec3(0.35, -0.45, 0.82)))) : 1.0;\n"
                  "  fragColor = vec4(vColour * d, 1.0);\n"
                  "}\n")
            : QStringLiteral(
                  "uniform bool lit;\n"
                  "varying vec3 vNormal;\n"
                  "varying vec3 vColour;\n"
                  "void main() {\n"
                  "  float d = lit ? 0.32 + 0.68 * abs(dot(normalize(vNormal), "
                  "normalize(vec3(0.35, -0.45, 0.82)))) : 1.0;\n"
                  "  gl_FragColor = vec4(vColour * d, 1.0);\n"
                  "}\n");

        program_ = std::make_unique<QOpenGLShaderProgram>();
        if (!program_->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                               vertexShader)
            || !program_->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                                  fragmentShader)
            || !program_->link()) {
            glError_ = QStringLiteral("OpenGL shader error: %1")
                           .arg(program_->log().trimmed());
            program_.reset();
            return;
        }
        if (!vao_.create() || !buffer_.create()) {
            glError_ = QStringLiteral("OpenGL buffer initialization failed");
            program_.reset();
        }
    }

    void paintGL() override {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (program_ != nullptr) {
            if (geometryDirty_) {
                rebuildGeometry();
            }
            drawGeometry();
        }
        drawHud();
    }

    void mousePressEvent(QMouseEvent* event) override {
        lastMouse_ = event->position();
        event->accept();
    }

    void mouseMoveEvent(QMouseEvent* event) override {
        const QPointF delta = event->position() - lastMouse_;
        lastMouse_ = event->position();
        if (event->buttons().testFlag(Qt::LeftButton)) {
            yaw_ += static_cast<float>(delta.x()) * 0.45F;
            pitch_ = std::clamp(
                pitch_ + static_cast<float>(delta.y()) * 0.45F,
                -89.0F, 89.0F);
            update();
        } else if (event->buttons().testFlag(Qt::MiddleButton)
                   || event->buttons().testFlag(Qt::RightButton)) {
            const float scale = distance_ /
                                static_cast<float>(std::max(height(), 1));
            const QMatrix4x4 inverse = viewMatrix().inverted();
            const QVector3D right = inverse.mapVector(QVector3D(1, 0, 0));
            const QVector3D up = inverse.mapVector(QVector3D(0, 1, 0));
            target_ -= right * static_cast<float>(delta.x()) * scale;
            target_ += up * static_cast<float>(delta.y()) * scale;
            update();
        }
        event->accept();
    }

    void mouseDoubleClickEvent(QMouseEvent* event) override {
        fitView();
        event->accept();
    }

    void wheelEvent(QWheelEvent* event) override {
        const double steps = event->angleDelta().y() / 120.0;
        distance_ = std::clamp(
            distance_ * static_cast<float>(std::pow(0.82, steps)),
            0.01F, 1.0e7F);
        update();
        event->accept();
    }

private:
    [[nodiscard]] static QVector3D frameCentre(
        const DiagnosticFrame& frame) {
        QVector3D minimum(
            std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max());
        QVector3D maximum(
            std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest());
        for (const DiagnosticVertex& vertex : frame.vertices) {
            const QVector3D point(
                static_cast<float>(vertex.positionMetres.x),
                static_cast<float>(vertex.positionMetres.y),
                static_cast<float>(vertex.positionMetres.z));
            minimum.setX(std::min(minimum.x(), point.x()));
            minimum.setY(std::min(minimum.y(), point.y()));
            minimum.setZ(std::min(minimum.z(), point.z()));
            maximum.setX(std::max(maximum.x(), point.x()));
            maximum.setY(std::max(maximum.y(), point.y()));
            maximum.setZ(std::max(maximum.z(), point.z()));
        }
        return 0.5F * (minimum + maximum);
    }

    [[nodiscard]] QMatrix4x4 viewMatrix() const {
        QMatrix4x4 view;
        view.translate(0.0F, 0.0F, -distance_);
        view.rotate(pitch_, 1.0F, 0.0F, 0.0F);
        view.rotate(yaw_, 0.0F, 0.0F, 1.0F);
        view.translate(-target_);
        return view;
    }

    [[nodiscard]] QMatrix4x4 viewProjection() const {
        QMatrix4x4 projection;
        projection.perspective(
            42.0F,
            static_cast<float>(std::max(width(), 1)) /
                static_cast<float>(std::max(height(), 1)),
            std::max(distance_ * 0.001F, 0.0001F),
            std::max(distance_ * 100.0F, 100.0F));
        return projection * viewMatrix();
    }

    const ScalarField* selectedField(double& minimum, double& maximum) const {
        if (frame_ == nullptr || fieldName_.isEmpty()) {
            return nullptr;
        }
        const std::string requested = fieldName_.toUtf8().toStdString();
        const auto found = std::find_if(
            frame_->scalarFields.begin(), frame_->scalarFields.end(),
            [&](const ScalarField& field) { return field.name == requested; });
        if (found == frame_->scalarFields.end() || found->values.empty()) {
            return nullptr;
        }
        const auto range = std::minmax_element(found->values.begin(),
                                               found->values.end());
        minimum = *range.first;
        maximum = *range.second;
        return &*found;
    }

    const VectorField* selectedVectorField() const {
        if (frame_ == nullptr || vectorFieldName_.isEmpty()) {
            return nullptr;
        }
        const std::string requested =
            vectorFieldName_.toUtf8().toStdString();
        const auto found = std::find_if(
            frame_->vectorFields.begin(), frame_->vectorFields.end(),
            [&](const VectorField& field) {
                return field.name == requested
                    && field.association == FieldAssociation::Vertex;
            });
        return found == frame_->vectorFields.end() ? nullptr : &*found;
    }

    void rebuildGeometry() {
        geometryDirty_ = false;
        vertices_.clear();
        triangleCount_ = 0;
        lineStart_ = 0;
        lineCount_ = 0;
        pointStart_ = 0;
        pointCount_ = 0;
        fieldMinimum_ = 0.0;
        fieldMaximum_ = 0.0;
        haveFieldRange_ = false;
        vectorMaximum_ = 0.0;
        vectorGlyphCount_ = 0;
        haveVectorField_ = false;
        if (frame_ == nullptr) {
            return;
        }

        double minimum = 0.0;
        double maximum = 0.0;
        const ScalarField* field = selectedField(minimum, maximum);
        if (field != nullptr) {
            fieldMinimum_ = minimum;
            fieldMaximum_ = maximum;
            haveFieldRange_ = true;
        }
        VectorGlyphGeometry glyphs;
        if (const VectorField* vectorField = selectedVectorField()) {
            glyphs = buildVertexVectorGlyphs(
                frame_->vertices, *vectorField);
            vectorMaximum_ = glyphs.maximumVectorMagnitude;
            vectorGlyphCount_ = glyphs.glyphCount;
            haveVectorField_ = true;
        }

        vertices_.reserve(frame_->triangles.size() * 3
                          + frame_->lines.size() * 2
                          + glyphs.segments.size() * 2
                          + frame_->vertices.size()
                          + frame_->contacts.size() + frame_->sealing.size());
        std::vector<bool> referencedVertices(frame_->vertices.size(), false);
        for (std::size_t triangleIndex = 0;
             triangleIndex < frame_->triangles.size(); ++triangleIndex) {
            const DiagnosticTriangle& triangle = frame_->triangles[triangleIndex];
            referencedVertices[triangle.vertex0] = true;
            referencedVertices[triangle.vertex1] = true;
            referencedVertices[triangle.vertex2] = true;
            const Vec3d positions[] = {
                frame_->vertices[triangle.vertex0].positionMetres,
                frame_->vertices[triangle.vertex1].positionMetres,
                frame_->vertices[triangle.vertex2].positionMetres};
            const QVector3D a(
                static_cast<float>(positions[0].x),
                static_cast<float>(positions[0].y),
                static_cast<float>(positions[0].z));
            const QVector3D b(
                static_cast<float>(positions[1].x),
                static_cast<float>(positions[1].y),
                static_cast<float>(positions[1].z));
            const QVector3D c(
                static_cast<float>(positions[2].x),
                static_cast<float>(positions[2].y),
                static_cast<float>(positions[2].z));
            const QVector3D normal = QVector3D::crossProduct(b - a, c - a).normalized();
            const std::uint32_t indices[] = {
                triangle.vertex0, triangle.vertex1, triangle.vertex2};
            for (std::size_t corner = 0; corner < 3; ++corner) {
                QVector3D colour(0.40F, 0.68F, 0.92F);
                if (field != nullptr) {
                    if (field->association == FieldAssociation::Global) {
                        colour = colourMap(field->values[0], minimum, maximum);
                    } else if (field->association == FieldAssociation::Vertex) {
                        colour = colourMap(field->values[indices[corner]],
                                           minimum, maximum);
                    } else if (field->association == FieldAssociation::Triangle) {
                        colour = colourMap(field->values[triangleIndex],
                                           minimum, maximum);
                    }
                }
                appendVertex(vertices_, positions[corner], normal, colour);
            }
        }
        triangleCount_ = static_cast<GLsizei>(vertices_.size());

        lineStart_ = triangleCount_;
        for (std::size_t lineIndex = 0; lineIndex < frame_->lines.size(); ++lineIndex) {
            const DiagnosticLine& line = frame_->lines[lineIndex];
            referencedVertices[line.vertex0] = true;
            referencedVertices[line.vertex1] = true;
            QVector3D colour = lineRoleColour(line.role);
            if (field != nullptr) {
                if (field->association == FieldAssociation::Global) {
                    colour = colourMap(field->values[0], minimum, maximum);
                } else if (field->association == FieldAssociation::Line) {
                    colour = colourMap(field->values[lineIndex], minimum, maximum);
                }
            }
            appendVertex(vertices_, frame_->vertices[line.vertex0].positionMetres,
                         {}, colour);
            appendVertex(vertices_, frame_->vertices[line.vertex1].positionMetres,
                         {}, colour);
        }
        for (const VectorGlyphSegment& segment : glyphs.segments) {
            QVector3D colour(0.98F, 0.82F, 0.25F);
            if (field != nullptr) {
                if (field->association == FieldAssociation::Global) {
                    colour = colourMap(field->values[0], minimum, maximum);
                } else if (field->association == FieldAssociation::Vertex) {
                    colour = colourMap(
                        field->values[segment.sourceVertexIndex],
                        minimum, maximum);
                }
            }
            appendVertex(vertices_, segment.startMetres, {}, colour);
            appendVertex(vertices_, segment.endMetres, {}, colour);
        }
        lineCount_ = static_cast<GLsizei>(vertices_.size()) - lineStart_;

        pointStart_ = static_cast<GLsizei>(vertices_.size());
        for (std::size_t vertexIndex = 0;
             vertexIndex < frame_->vertices.size(); ++vertexIndex) {
            if (referencedVertices[vertexIndex]) {
                continue;
            }
            QVector3D colour(0.40F, 0.68F, 0.92F);
            if (field != nullptr) {
                if (field->association == FieldAssociation::Global) {
                    colour = colourMap(field->values[0], minimum, maximum);
                } else if (field->association == FieldAssociation::Vertex) {
                    colour = colourMap(
                        field->values[vertexIndex], minimum, maximum);
                }
            }
            appendVertex(
                vertices_, frame_->vertices[vertexIndex].positionMetres,
                {}, colour);
        }
        for (const ContactMarker& marker : frame_->contacts) {
            appendVertex(vertices_, marker.positionMetres, {},
                         QVector3D(1.0F, 0.22F, 0.18F));
        }
        for (const SealingMarker& marker : frame_->sealing) {
            appendVertex(vertices_, marker.positionMetres, {},
                         marker.sealed ? QVector3D(1.0F, 0.35F, 0.12F)
                                       : QVector3D(0.25F, 1.0F, 0.50F));
        }
        pointCount_ = static_cast<GLsizei>(vertices_.size()) - pointStart_;

        vao_.bind();
        buffer_.bind();
        buffer_.allocate(vertices_.data(),
                         static_cast<int>(vertices_.size() * sizeof(RenderVertex)));
        program_->bind();
        program_->enableAttributeArray("position");
        program_->setAttributeBuffer(
            "position", GL_FLOAT, offsetof(RenderVertex, position), 3,
            sizeof(RenderVertex));
        program_->enableAttributeArray("normal");
        program_->setAttributeBuffer(
            "normal", GL_FLOAT, offsetof(RenderVertex, normal), 3,
            sizeof(RenderVertex));
        program_->enableAttributeArray("colour");
        program_->setAttributeBuffer(
            "colour", GL_FLOAT, offsetof(RenderVertex, colour), 3,
            sizeof(RenderVertex));
        program_->release();
        buffer_.release();
        vao_.release();
    }

    void drawGeometry() {
        if (vertices_.empty()) {
            return;
        }
        program_->bind();
        vao_.bind();
        program_->setUniformValue("mvp", viewProjection());
        if (triangleCount_ > 0) {
            program_->setUniformValue("lit", true);
            glDrawArrays(GL_TRIANGLES, 0, triangleCount_);
        }
        if (lineCount_ > 0) {
            program_->setUniformValue("lit", false);
            glLineWidth(1.5F);
            glDrawArrays(GL_LINES, lineStart_, lineCount_);
        }
        if (pointCount_ > 0) {
            program_->setUniformValue("lit", false);
            glPointSize(7.0F);
            glDrawArrays(GL_POINTS, pointStart_, pointCount_);
        }
        vao_.release();
        program_->release();
    }

    void drawHud() {
        QPainter painter(this);
        painter.setRenderHint(QPainter::TextAntialiasing);
        const QColor panel(8, 12, 20, 205);
        painter.fillRect(QRect(12, 12, std::min(width() - 24, 610), 113), panel);
        painter.setPen(QColor(230, 235, 244));
        int y = 33;
        if (frame_ != nullptr) {
            painter.drawText(
                22, y,
                QStringLiteral("step %1   t %2 s   dt %3 s   coupling %4")
                    .arg(frame_->step)
                    .arg(frame_->simulationTimeSeconds, 0, 'g', 8)
                    .arg(frame_->timeStepSeconds, 0, 'g', 6)
                    .arg(frame_->couplingIteration));
            y += 20;
            painter.drawText(
                22, y,
                QStringLiteral("scene %1   solver %2")
                    .arg(abbreviated(frame_->sceneChecksum),
                         abbreviated(frame_->solverCommit)));
            y += 20;
            painter.drawText(
                22, y,
                QStringLiteral("residual: dx %1 m   traction %2 N   fluid %3")
                    .arg(frame_->couplingResiduals.displacementMetres, 0, 'g', 4)
                    .arg(frame_->couplingResiduals.tractionNewtons, 0, 'g', 4)
                    .arg(frame_->couplingResiduals.fluid, 0, 'g', 4));
            y += 20;
            QString field = fieldName_.isEmpty()
                ? QStringLiteral("plain surface") : fieldName_;
            if (haveFieldRange_) {
                field += QStringLiteral(" [%1, %2]")
                             .arg(fieldMinimum_, 0, 'g', 5)
                             .arg(fieldMaximum_, 0, 'g', 5);
            }
            if (haveVectorField_) {
                field += QStringLiteral("   vectors %1 (%2, max %3)")
                             .arg(vectorFieldName_)
                             .arg(vectorGlyphCount_)
                             .arg(vectorMaximum_, 0, 'g', 5);
            }
            painter.drawText(22, y, field);
        } else {
            painter.drawText(22, y, QStringLiteral("Waiting for a diagnostic frame..."));
        }
        painter.setPen(QColor(158, 171, 190));
        painter.drawText(22, 116, QStringLiteral(
            "Left drag: orbit   Middle/right drag: pan   Wheel: zoom   Double-click: fit"));

        const QString error = !glError_.isEmpty() ? glError_ : message_;
        if (!error.isEmpty()) {
            const QRect errorRect(12, height() - 62, width() - 24, 48);
            painter.fillRect(errorRect, QColor(105, 22, 22, 225));
            painter.setPen(Qt::white);
            painter.drawText(errorRect.adjusted(10, 5, -10, -5),
                             Qt::AlignLeft | Qt::AlignVCenter | Qt::TextWordWrap,
                             error);
        }
    }

    std::shared_ptr<const DiagnosticFrame> frame_;
    QString fieldName_;
    QString vectorFieldName_;
    QString message_;
    QString glError_;
    bool geometryDirty_ = true;
    std::vector<RenderVertex> vertices_;
    GLsizei triangleCount_ = 0;
    GLsizei lineStart_ = 0;
    GLsizei lineCount_ = 0;
    GLsizei pointStart_ = 0;
    GLsizei pointCount_ = 0;
    double fieldMinimum_ = 0.0;
    double fieldMaximum_ = 0.0;
    bool haveFieldRange_ = false;
    double vectorMaximum_ = 0.0;
    std::uint32_t vectorGlyphCount_ = 0;
    bool haveVectorField_ = false;

    std::unique_ptr<QOpenGLShaderProgram> program_;
    QOpenGLVertexArrayObject vao_;
    QOpenGLBuffer buffer_{QOpenGLBuffer::VertexBuffer};
    QVector3D target_;
    float distance_ = 5.0F;
    float yaw_ = 25.0F;
    float pitch_ = -55.0F;
    QPointF lastMouse_;
};

DiagnosticFrame smokeFrame() {
    DiagnosticFrame frame;
    frame.sceneChecksum = "smoke-scene";
    frame.solverCommit = "simwing-viewer-smoke";
    frame.step = 12;
    frame.simulationTimeSeconds = 0.2;
    frame.timeStepSeconds = 1.0 / 60.0;
    frame.couplingIteration = 2;
    frame.couplingResiduals = {1.2e-5, 0.018, 3.4e-6, 8.0e-7, 0.002};
    frame.vertices = {
        {1, {-2.4, 0.0, 0.0}}, {2, {-1.2, 0.45, 0.42}},
        {3, {0.0, 0.62, 0.58}}, {4, {1.2, 0.45, 0.42}},
        {5, {2.4, 0.0, 0.0}}, {6, {-1.2, -0.55, 0.18}},
        {7, {0.0, -0.72, 0.24}}, {8, {1.2, -0.55, 0.18}},
        {9, {0.0, 0.0, -2.7}}};
    frame.triangles = {
        {101, 0, 1, 5, 1, 2}, {102, 1, 6, 5, 1, 2},
        {103, 1, 2, 6, 1, 2}, {104, 2, 3, 7, 1, 2},
        {105, 2, 7, 6, 1, 2}, {106, 3, 4, 7, 1, 2}};
    frame.lines = {
        {201, 0, 8, 1}, {202, 2, 8, 2}, {203, 4, 8, 3}};
    frame.contacts = {
        {301, {EntityKind::Vertex, 2}, {EntityKind::Triangle, 3},
         {0.0, -0.1, 0.34}, {0.0, 0.0, 1.0}, -0.001}};
    frame.sealing = {
        {401, 1, 2, {0.0, -0.6, 0.22}, {0.0, -1.0, 0.0}, 0.0001, true}};
    frame.scalarFields = {
        {"strain", "1", FieldAssociation::Triangle,
         {0.002, 0.004, 0.008, 0.012, 0.009, 0.003}},
        {"line tension", "N", FieldAssociation::Line, {82.0, 118.0, 91.0}}};
    return frame;
}

} // namespace

class ViewerWindow::Impl {
public:
    explicit Impl(ViewerWindow* owner) : owner_(owner) {
        view_ = new TraceView(owner_);
        owner_->setCentralWidget(view_);
        owner_->resize(1100, 760);
        owner_->setWindowTitle(QStringLiteral("SimWing diagnostic viewer"));

        QToolBar* toolbar = owner_->addToolBar(QStringLiteral("Replay"));
        toolbar->setMovable(false);
        playAction_ = toolbar->addAction(QStringLiteral("Pause"));
        playAction_->setShortcut(QKeySequence(Qt::Key_Space));
        stepAction_ = toolbar->addAction(QStringLiteral("Step"));
        stepAction_->setShortcut(QKeySequence(Qt::Key_Right));
        restartAction_ = toolbar->addAction(QStringLiteral("Restart"));
        restartAction_->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_R));
        fitAction_ = toolbar->addAction(QStringLiteral("Fit"));
        fitAction_->setShortcut(QKeySequence(Qt::Key_F));
        toolbar->addSeparator();
        toolbar->addWidget(new QLabel(QStringLiteral(" Field: "), toolbar));
        fieldCombo_ = new QComboBox(toolbar);
        fieldCombo_->addItem(QStringLiteral("Plain"), QString());
        fieldCombo_->setMinimumContentsLength(20);
        toolbar->addWidget(fieldCombo_);
        toolbar->addWidget(new QLabel(QStringLiteral(" Vectors: "), toolbar));
        vectorCombo_ = new QComboBox(toolbar);
        vectorCombo_->addItem(QStringLiteral("None"), QString());
        vectorCombo_->setMinimumContentsLength(15);
        toolbar->addWidget(vectorCombo_);

        QObject::connect(playAction_, &QAction::triggered, owner_, [this] {
            playing_ = !playing_;
            updatePlayAction();
            playbackClock_.restart();
        });
        QObject::connect(stepAction_, &QAction::triggered, owner_, [this] {
            playing_ = false;
            updatePlayAction();
            if (pendingFrame_ != nullptr) {
                presentPending();
            } else if (!traceEnded_) {
                stepWhenReady_ = true;
                requestNext();
            }
        });
        QObject::connect(restartAction_, &QAction::triggered, owner_, [this] {
            if (!traceFile_.isEmpty()) {
                startTrace(traceFile_, follow_);
            }
        });
        QObject::connect(fitAction_, &QAction::triggered, view_, [this] {
            view_->fitView();
        });
        QObject::connect(fieldCombo_, &QComboBox::currentIndexChanged, owner_,
                         [this](int index) {
            desiredField_ = fieldCombo_->itemData(index).toString();
            view_->setField(desiredField_);
        });
        QObject::connect(
            vectorCombo_, &QComboBox::currentIndexChanged, owner_,
            [this](int index) {
                desiredVectorField_ =
                    vectorCombo_->itemData(index).toString();
                view_->setVectorField(desiredVectorField_);
            });

        timer_.setInterval(16);
        QObject::connect(&timer_, &QTimer::timeout, owner_, [this] { tick(); });
        timer_.start();
        owner_->statusBar()->showMessage(QStringLiteral("No trace loaded"));
    }

    ~Impl() { stream_.reset(); }

    bool loadTrace(
        const QString& fileName,
        QString* errorMessage,
        bool follow) {
        const QFileInfo info(fileName);
        if (!info.exists() || !info.isFile() || !info.isReadable()) {
            const QString message = QStringLiteral("Trace file is not readable: %1")
                                        .arg(info.absoluteFilePath());
            if (errorMessage != nullptr) {
                *errorMessage = message;
            }
            return false;
        }
        startTrace(info.absoluteFilePath(), follow);
        return true;
    }

    void showSmokeFrame() {
        ++generation_;
        stream_.reset();
        traceFile_.clear();
        traceEnded_ = true;
        playing_ = false;
        currentFrame_ = std::make_shared<const DiagnosticFrame>(smokeFrame());
        pendingFrame_.reset();
        view_->setMessage({});
        view_->setFrame(currentFrame_);
        refreshFields();
        view_->fitView();
        restartAction_->setEnabled(false);
        updatePlayAction();
        owner_->setWindowTitle(QStringLiteral("SimWing diagnostic viewer - smoke test"));
        owner_->statusBar()->showMessage(QStringLiteral("Built-in diagnostic frame"));
    }

    QString renderError() const { return view_->renderError(); }

private:
    void post(std::function<void()> callback) {
        QMetaObject::invokeMethod(owner_, std::move(callback),
                                  Qt::QueuedConnection);
    }

    void startTrace(const QString& fileName, bool follow) {
        ++generation_;
        const std::uint64_t generation = generation_;
        stream_.reset();
        traceFile_ = fileName;
        follow_ = follow;
        traceHeader_ = {};
        currentFrame_.reset();
        pendingFrame_.reset();
        traceEnded_ = false;
        readInFlight_ = true; // The stream requests its first frame at startup.
        stepWhenReady_ = false;
        playing_ = true;
        desiredField_.clear();
        view_->setField({});
        view_->setFrame({});
        view_->setMessage({});
        restartAction_->setEnabled(true);
        updatePlayAction();
        owner_->setWindowTitle(
            follow
                ? QStringLiteral("SimWing diagnostic viewer - %1 (following)")
                      .arg(QFileInfo(fileName).fileName())
                : QStringLiteral("SimWing diagnostic viewer - %1")
                      .arg(QFileInfo(fileName).fileName()));
        owner_->statusBar()->showMessage(
            follow ? QStringLiteral("Waiting for live trace frames...")
                   : QStringLiteral("Reading trace..."));

        auto headerCallback = [this, generation](TraceHeader header) {
            post([this, generation, header = std::move(header)]() mutable {
                if (generation == generation_) {
                    traceHeader_ = std::move(header);
                    owner_->statusBar()->showMessage(
                        QStringLiteral("Trace %1 / %2")
                            .arg(abbreviated(traceHeader_.sceneChecksum),
                                 abbreviated(traceHeader_.solverCommit)));
                }
            });
        };
        auto frameCallback =
            [this, generation](std::shared_ptr<const DiagnosticFrame> frame) {
                post([this, generation, frame = std::move(frame)]() mutable {
                    if (generation != generation_) {
                        return;
                    }
                    readInFlight_ = false;
                    if (currentFrame_ == nullptr) {
                        currentFrame_ = std::move(frame);
                        view_->setFrame(currentFrame_);
                        refreshFields();
                        view_->fitView();
                        playbackClock_.restart();
                        requestNext();
                    } else {
                        pendingFrame_ = std::move(frame);
                        if (stepWhenReady_) {
                            stepWhenReady_ = false;
                            presentPending();
                        }
                    }
                });
            };
        auto errorCallback = [this, generation](QString message) {
            post([this, generation, message = std::move(message)] {
                if (generation == generation_) {
                    readInFlight_ = false;
                    traceEnded_ = true;
                    playing_ = false;
                    updatePlayAction();
                    view_->setMessage(message);
                    owner_->statusBar()->showMessage(message);
                    QMessageBox::critical(owner_, QStringLiteral("Trace error"),
                                          message);
                }
            });
        };
        auto endCallback = [this, generation] {
            post([this, generation] {
                if (generation == generation_) {
                    readInFlight_ = false;
                    traceEnded_ = true;
                    if (pendingFrame_ == nullptr) {
                        playing_ = false;
                        updatePlayAction();
                    }
                    owner_->statusBar()->showMessage(
                        follow_ ? QStringLiteral("Live trace completed")
                                : QStringLiteral("Trace completed"));
                    if (follow_) {
                        owner_->setWindowTitle(
                            QStringLiteral("SimWing diagnostic viewer - %1 (complete)")
                                .arg(QFileInfo(traceFile_).fileName()));
                    }
                }
            });
        };
        stream_ = std::make_unique<TraceStream>(
            fileName, follow, std::move(headerCallback),
            std::move(frameCallback), std::move(errorCallback),
            std::move(endCallback));
    }

    void requestNext() {
        if (stream_ != nullptr && !readInFlight_ && !traceEnded_) {
            readInFlight_ = true;
            stream_->requestNext();
        }
    }

    void presentPending() {
        if (pendingFrame_ == nullptr) {
            return;
        }
        currentFrame_ = std::move(pendingFrame_);
        view_->setFrame(currentFrame_);
        refreshFields();
        playbackClock_.restart();
        requestNext();
    }

    void tick() {
        if (!playing_ || pendingFrame_ == nullptr || currentFrame_ == nullptr) {
            return;
        }
        double delaySeconds = pendingFrame_->simulationTimeSeconds
                              - currentFrame_->simulationTimeSeconds;
        if (!(delaySeconds > 0.0)) {
            delaySeconds = currentFrame_->timeStepSeconds;
        }
        delaySeconds = std::clamp(delaySeconds, 0.001, 1.0);
        if (!playbackClock_.isValid()
            || playbackClock_.elapsed() >= delaySeconds * 1000.0) {
            presentPending();
        }
    }

    void refreshFields() {
        if (currentFrame_ == nullptr) {
            return;
        }
        std::vector<const ScalarField*> fields;
        for (const ScalarField& field : currentFrame_->scalarFields) {
            if (field.association == FieldAssociation::Global
                || field.association == FieldAssociation::Vertex
                || field.association == FieldAssociation::Triangle
                || field.association == FieldAssociation::Line) {
                fields.push_back(&field);
            }
        }

        bool same = fieldCombo_->count() == static_cast<int>(fields.size() + 1);
        for (int i = 0; same && i < static_cast<int>(fields.size()); ++i) {
            same = fieldCombo_->itemData(i + 1).toString()
                   == fromUtf8(fields[i]->name);
        }
        if (!same) {
            fieldCombo_->blockSignals(true);
            fieldCombo_->clear();
            fieldCombo_->addItem(QStringLiteral("Plain"), QString());
            int desiredIndex = 0;
            for (std::size_t i = 0; i < fields.size(); ++i) {
                const ScalarField& field = *fields[i];
                const QString fieldName = fromUtf8(field.name);
                QString label = fieldName;
                if (!field.unit.empty()) {
                    label += QStringLiteral(" [%1]").arg(fromUtf8(field.unit));
                }
                fieldCombo_->addItem(label, fieldName);
                if (fieldName == desiredField_) {
                    desiredIndex = static_cast<int>(i + 1);
                }
            }
            fieldCombo_->setCurrentIndex(desiredIndex);
            fieldCombo_->blockSignals(false);
            if (desiredIndex == 0 && !desiredField_.isEmpty()) {
                desiredField_.clear();
                view_->setField({});
            }
        }

        std::vector<const VectorField*> vectorFields;
        for (const VectorField& field : currentFrame_->vectorFields) {
            if (field.association == FieldAssociation::Vertex) {
                vectorFields.push_back(&field);
            }
        }
        same = vectorCombo_->count()
            == static_cast<int>(vectorFields.size() + 1);
        for (int i = 0;
             same && i < static_cast<int>(vectorFields.size()); ++i) {
            same = vectorCombo_->itemData(i + 1).toString()
                == fromUtf8(vectorFields[i]->name);
        }
        if (!same) {
            vectorCombo_->blockSignals(true);
            vectorCombo_->clear();
            vectorCombo_->addItem(QStringLiteral("None"), QString());
            int desiredIndex = 0;
            for (std::size_t i = 0; i < vectorFields.size(); ++i) {
                const VectorField& field = *vectorFields[i];
                const QString fieldName = fromUtf8(field.name);
                QString label = fieldName;
                if (!field.unit.empty()) {
                    label += QStringLiteral(" [%1]").arg(fromUtf8(field.unit));
                }
                vectorCombo_->addItem(label, fieldName);
                if (fieldName == desiredVectorField_) {
                    desiredIndex = static_cast<int>(i + 1);
                }
            }
            vectorCombo_->setCurrentIndex(desiredIndex);
            vectorCombo_->blockSignals(false);
            if (desiredIndex == 0 && !desiredVectorField_.isEmpty()) {
                desiredVectorField_.clear();
                view_->setVectorField({});
            }
        }
    }

    void updatePlayAction() {
        playAction_->setText(playing_ ? QStringLiteral("Pause")
                                     : QStringLiteral("Play"));
    }

    ViewerWindow* owner_ = nullptr;
    TraceView* view_ = nullptr;
    QAction* playAction_ = nullptr;
    QAction* stepAction_ = nullptr;
    QAction* restartAction_ = nullptr;
    QAction* fitAction_ = nullptr;
    QComboBox* fieldCombo_ = nullptr;
    QComboBox* vectorCombo_ = nullptr;
    QTimer timer_;
    QElapsedTimer playbackClock_;
    QString traceFile_;
    bool follow_ = false;
    QString desiredField_;
    QString desiredVectorField_;
    TraceHeader traceHeader_;
    std::unique_ptr<TraceStream> stream_;
    std::shared_ptr<const DiagnosticFrame> currentFrame_;
    std::shared_ptr<const DiagnosticFrame> pendingFrame_;
    std::uint64_t generation_ = 0;
    bool playing_ = true;
    bool traceEnded_ = false;
    bool readInFlight_ = false;
    bool stepWhenReady_ = false;
};

ViewerWindow::ViewerWindow(QWidget* parent)
    : QMainWindow(parent), impl_(std::make_unique<Impl>(this)) {}

ViewerWindow::~ViewerWindow() = default;

bool ViewerWindow::loadTrace(
    const QString& fileName,
    QString* errorMessage,
    bool follow) {
    return impl_->loadTrace(fileName, errorMessage, follow);
}

void ViewerWindow::showSmokeFrame() {
    impl_->showSmokeFrame();
}

QString ViewerWindow::renderError() const {
    return impl_->renderError();
}

} // namespace simwing::viewer
