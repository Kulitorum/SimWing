#pragma once

#include <QMainWindow>
#include <QString>

#include <memory>

namespace simwing::viewer {

// Standalone, read-only diagnostic trace viewer. Trace decoding happens on a
// background thread and the renderer retains only the current and prefetched
// frames. It never exposes mutable frame data to OpenGL or to UI controls.
class ViewerWindow final : public QMainWindow {
public:
    explicit ViewerWindow(QWidget* parent = nullptr);
    ~ViewerWindow() override;

    [[nodiscard]] bool loadTrace(
        const QString& fileName,
        QString* errorMessage = nullptr);
    void showSmokeFrame();

    [[nodiscard]] QString renderError() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simwing::viewer
