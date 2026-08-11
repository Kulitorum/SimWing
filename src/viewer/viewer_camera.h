#pragma once

namespace simwing::viewer {

struct ViewerDepthRange {
    double nearPlaneMetres = 0.0;
    double farPlaneMetres = 0.0;

    bool operator==(const ViewerDepthRange&) const = default;
};

struct ViewerProjectedAxis {
    double horizontal = 0.0;
    double vertical = 0.0;
    double cameraDepth = 0.0;

    bool operator==(const ViewerProjectedAxis&) const = default;
};

struct ViewerOrientationTriad {
    ViewerProjectedAxis x;
    ViewerProjectedAxis y;
    ViewerProjectedAxis z;

    bool operator==(const ViewerOrientationTriad&) const = default;
};

// Returns a tight positive perspective depth interval around a scene-sized
// sphere centred at cameraDistanceMetres. Keeping this interval close to the
// actual geometry preserves depth-buffer precision across metre-scale and
// very small diagnostic cases.
[[nodiscard]] ViewerDepthRange viewerDepthRange(
    double cameraDistanceMetres,
    double sceneRadiusMetres);

// Projects the positive world axes through the viewer's orbit rotation.
// Horizontal and vertical use screen-space signs before the GUI's downward-Y
// pixel conversion; cameraDepth is retained so overlays can paint far axes
// before near axes. Translation, zoom, and perspective do not affect this
// orientation-only indicator.
[[nodiscard]] ViewerOrientationTriad viewerOrientationTriad(
    double yawDegrees,
    double pitchDegrees);

} // namespace simwing::viewer
