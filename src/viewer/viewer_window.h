#pragma once

#include <QMainWindow>
#include <QSurfaceFormat>
#include <QString>

#include <memory>

namespace simwing::viewer {

// The standalone entry point and the QOpenGLWidget both request this format.
// Keeping the per-widget request avoids silently losing the depth attachment
// when ViewerWindow is embedded by a test or another executable.
[[nodiscard]] QSurfaceFormat diagnosticViewerSurfaceFormat();

// Standalone, read-only diagnostic trace viewer. Trace decoding happens on a
// background thread and the renderer retains only the current and prefetched
// frames. It never exposes mutable frame data to OpenGL or to UI controls.
class ViewerWindow final : public QMainWindow {
public:
    explicit ViewerWindow(QWidget* parent = nullptr);
    ~ViewerWindow() override;

    [[nodiscard]] bool loadTrace(
        const QString& fileName,
        QString* errorMessage = nullptr,
        bool follow = false);
    void showSmokeFrame();

    [[nodiscard]] QString renderError() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simwing::viewer
