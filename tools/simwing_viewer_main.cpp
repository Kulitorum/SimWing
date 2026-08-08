#include "viewer/viewer_window.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QTimer>

#include <cstdio>

int main(int argc, char* argv[]) {
    QSurfaceFormat::setDefaultFormat(
        simwing::viewer::diagnosticViewerSurfaceFormat());

    QApplication app(argc, argv);
    QApplication::setApplicationName(QStringLiteral("simwing-viewer"));
    QApplication::setApplicationVersion(QStringLiteral("0.1"));

    QCommandLineParser parser;
    parser.setApplicationDescription(QStringLiteral(
        "Replay immutable SimWing diagnostic trace files."));
    parser.addHelpOption();
    parser.addVersionOption();
    const QCommandLineOption smokeOption(
        QStringLiteral("smoke-test"),
        QStringLiteral("Render a built-in frame briefly and exit."));
    parser.addOption(smokeOption);
    const QCommandLineOption followOption(
        QStringLiteral("follow"),
        QStringLiteral(
            "Follow a growing trace until its explicit finish marker."));
    parser.addOption(followOption);
    parser.addPositionalArgument(
        QStringLiteral("trace"),
        QStringLiteral("SimWing viewer trace to replay."),
        QStringLiteral("[trace]"));
    parser.process(app);

    const QStringList positional = parser.positionalArguments();
    const bool smokeTest = parser.isSet(smokeOption);
    const bool follow = parser.isSet(followOption);
    if ((!smokeTest && positional.size() != 1)
        || (smokeTest && (positional.size() > 1 || follow))) {
        std::fprintf(stderr,
                     "Usage: simwing-viewer [--follow] <trace>\n"
                     "       simwing-viewer --smoke-test\n");
        return 2;
    }

    simwing::viewer::ViewerWindow window;
    if (smokeTest) {
        window.showSmokeFrame();
    } else {
        QString error;
        if (!window.loadTrace(positional.front(), &error, follow)) {
            std::fprintf(stderr, "%s\n", error.toLocal8Bit().constData());
            return 2;
        }
    }
    window.show();

    if (smokeTest) {
        QTimer::singleShot(650, &app, [&app, &window] {
            if (!window.renderError().isEmpty()) {
                std::fprintf(stderr, "%s\n",
                             window.renderError().toLocal8Bit().constData());
                app.exit(1);
            } else {
                app.exit(0);
            }
        });
    }
    return app.exec();
}
