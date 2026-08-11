#include "viewer_camera.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace simwing::viewer {

ViewerDepthRange viewerDepthRange(double cameraDistanceMetres,
                                  double sceneRadiusMetres) {
    if (!std::isfinite(cameraDistanceMetres)
        || !(cameraDistanceMetres > 0.0)
        || !std::isfinite(sceneRadiusMetres)
        || sceneRadiusMetres < 0.0) {
        throw std::invalid_argument(
            "Viewer depth range requires finite positive camera geometry");
    }

    const double padding = std::max(
        0.25 * sceneRadiusMetres, 0.01 * cameraDistanceMetres);
    const double coverage = sceneRadiusMetres + padding;
    const double minimumNear = std::max(
        1.0e-5, 1.0e-4 * cameraDistanceMetres);
    const double nearPlane = std::max(
        minimumNear, cameraDistanceMetres - coverage);
    const double minimumSpan = std::max(
        1.0e-4, 1.0e-3 * cameraDistanceMetres);
    const double farPlane = std::max(
        cameraDistanceMetres + coverage, nearPlane + minimumSpan);
    return {nearPlane, farPlane};
}

ViewerOrientationTriad viewerOrientationTriad(const double yawDegrees,
                                              const double pitchDegrees) {
    if (!std::isfinite(yawDegrees) || !std::isfinite(pitchDegrees)) {
        throw std::invalid_argument(
            "Viewer orientation requires finite orbit angles");
    }
    constexpr double radiansPerDegree =
        3.141592653589793238462643383279502884 / 180.0;
    const double yaw = yawDegrees * radiansPerDegree;
    const double pitch = pitchDegrees * radiansPerDegree;
    const double cy = std::cos(yaw);
    const double sy = std::sin(yaw);
    const double cp = std::cos(pitch);
    const double sp = std::sin(pitch);
    return {
        {cy, cp * sy, sp * sy},
        {-sy, cp * cy, sp * cy},
        {0.0, -sp, cp},
    };
}

} // namespace simwing::viewer
