#pragma once

namespace simwing::viewer {

struct ViewerDepthRange {
    double nearPlaneMetres = 0.0;
    double farPlaneMetres = 0.0;

    bool operator==(const ViewerDepthRange&) const = default;
};

// Returns a tight positive perspective depth interval around a scene-sized
// sphere centred at cameraDistanceMetres. Keeping this interval close to the
// actual geometry preserves depth-buffer precision across metre-scale and
// very small diagnostic cases.
[[nodiscard]] ViewerDepthRange viewerDepthRange(
    double cameraDistanceMetres,
    double sceneRadiusMetres);

} // namespace simwing::viewer
