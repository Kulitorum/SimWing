#include "viewer_camera.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {

using namespace simwing::viewer;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(double actual,
               double expected,
               double tolerance,
               const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void testFittedSceneRange() {
    const ViewerDepthRange range = viewerDepthRange(3.8, 1.0);
    checkNear(range.nearPlaneMetres, 2.55, 1.0e-15,
              "depth range: fitted scene has a tight near plane");
    checkNear(range.farPlaneMetres, 5.05, 1.0e-15,
              "depth range: fitted scene has a tight far plane");
    check(range.farPlaneMetres / range.nearPlaneMetres < 2.0,
          "depth range: fitted scene uses available depth precision");
}

void testSmallAndCloseScenesRemainVisible() {
    const ViewerDepthRange point = viewerDepthRange(5.0, 0.0);
    checkNear(point.nearPlaneMetres, 4.95, 1.0e-15,
              "depth range: point scenes retain camera-relative padding");
    checkNear(point.farPlaneMetres, 5.05, 1.0e-15,
              "depth range: point scenes have a nonzero interval");

    const ViewerDepthRange inside = viewerDepthRange(0.25, 1.0);
    check(inside.nearPlaneMetres > 0.0
              && inside.farPlaneMetres >= 1.5,
          "depth range: close cameras retain positive near depth and full reach");
}

void testInvalidInputsAreRejected() {
    expectRejected(
        [] { static_cast<void>(viewerDepthRange(0.0, 1.0)); },
        "depth range: zero camera distance is rejected");
    expectRejected(
        [] { static_cast<void>(viewerDepthRange(1.0, -1.0)); },
        "depth range: negative scene radius is rejected");
    expectRejected(
        [] {
            static_cast<void>(viewerDepthRange(
                std::numeric_limits<double>::quiet_NaN(), 1.0));
        },
        "depth range: non-finite camera geometry is rejected");
}

void testOrientationTriadUsesWorldAxes() {
    const auto front = viewerOrientationTriad(0.0, 0.0);
    check(front.x == ViewerProjectedAxis{1.0, 0.0, 0.0}
              && front.y == ViewerProjectedAxis{0.0, 1.0, 0.0}
              && front.z == ViewerProjectedAxis{0.0, 0.0, 1.0},
          "orientation triad: unrotated world axes retain X/Y/Z identity");

    const auto yawed = viewerOrientationTriad(90.0, 0.0);
    checkNear(yawed.x.horizontal, 0.0, 1.0e-15,
              "orientation triad: yawed X has no horizontal component");
    checkNear(yawed.x.vertical, 1.0, 1.0e-15,
              "orientation triad: positive yaw maps X upward");
    checkNear(yawed.y.horizontal, -1.0, 1.0e-15,
              "orientation triad: positive yaw maps Y left");
    checkNear(yawed.z.cameraDepth, 1.0, 1.0e-15,
              "orientation triad: yaw leaves Z depth unchanged");

    const auto pitched = viewerOrientationTriad(0.0, 90.0);
    checkNear(pitched.y.cameraDepth, 1.0, 1.0e-15,
              "orientation triad: positive pitch maps Y toward camera");
    checkNear(pitched.z.vertical, -1.0, 1.0e-15,
              "orientation triad: positive pitch maps Z downward");

    expectRejected(
        [] {
            static_cast<void>(viewerOrientationTriad(
                std::numeric_limits<double>::infinity(), 0.0));
        },
        "orientation triad: non-finite orbit angles are rejected");
}

} // namespace

int main() {
    testFittedSceneRange();
    testSmallAndCloseScenesRemainVisible();
    testInvalidInputsAreRejected();
    testOrientationTriadUsesWorldAxes();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d viewer camera check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all viewer camera checks passed");
    return 0;
}
