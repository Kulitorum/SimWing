#include "open_piston_case.h"
#include "viewer_protocol.h"

#include <cmath>
#include <cstdio>
#include <ranges>
#include <string_view>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g, tolerance %.3g)\n",
                     message, actual, expected, tolerance);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "open piston: deterministic frame serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const std::string_view name) {
    const auto found = std::ranges::find(
        frame.scalarFields, name, &viewer::ScalarField::name);
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

void testVisibleOpenPistonAndDeterminism() {
    fsi::OpenPistonCase first;
    fsi::OpenPistonCase second;
    const auto firstInitialFrame = first.advance();
    const auto secondInitialFrame = second.advance();
    check(serialized(firstInitialFrame) == serialized(secondInitialFrame),
          "open piston: first accelerated frame replays bit-for-bit");

    const auto* initialPressure = scalarField(
        firstInitialFrame, "interface.pressure_traction");
    const auto* initialActuator = scalarField(
        firstInitialFrame, "actuator.step_impulse");
    const auto* initialGcl = scalarField(
        firstInitialFrame, "fluid.gcl_residual");
    check(initialPressure != nullptr
              && initialPressure->values.size() == 2
              && initialPressure->values.front() < 0.0
              && initialPressure->values[0] == initialPressure->values[1],
          "open piston: first frame exposes uniform resisting CFD pressure");
    check(initialActuator != nullptr
              && initialActuator->values.size() == 1
              && initialActuator->values.front() > 300.0,
          "open piston: explicit actuator overcomes inertia and fluid reaction");
    check(initialGcl != nullptr
              && std::abs(initialGcl->values.front()) < 1.0e-13,
          "open piston: first visible frame closes geometric conservation");

    viewer::DiagnosticFrame finalFrame = firstInitialFrame;
    constexpr std::uint64_t steps = 24;
    for (std::uint64_t step = 1; step < steps; ++step) {
        finalFrame = first.advance();
        const auto replay = second.advance();
        check(serialized(finalFrame) == serialized(replay),
              "open piston: continued accepted frames replay bit-for-bit");
    }

    const auto checkpoint = first.structure().checkpoint();
    check(checkpoint.acceptedStepCount == steps
              && finalFrame.step == steps
              && finalFrame.sceneChecksum == fsi::openPistonCaseChecksum
              && finalFrame.solverCommit == fsi::openPistonCaseSolverId,
          "open piston: trace provenance and accepted step remain aligned");
    check(finalFrame.triangles.size() == 2
              && finalFrame.triangles[0].negativeRegionId == 9
              && finalFrame.triangles[0].positiveRegionId == 9,
          "open piston: viewer exposes a nonseparating two-sided sheet");
    const double expectedOffset = 0.05
        * first.stepSettings().timeStepSeconds
        * static_cast<double>(steps);
    checkNear(first.surfaceOffsetMeters(), expectedOffset, 2.0e-17,
              "open piston: accepted plate displacement drives cut-cell offset");
    for (const auto& node : checkpoint.nodes) {
        checkNear(node.positionMeters.x, 3.0 + expectedOffset, 2.0e-14,
                  "open piston: prescribed actuator leaves visible translation");
        checkNear(node.velocityMetersPerSecond.x, 0.05, 2.0e-14,
                  "open piston: plate coasts at the prescribed verification speed");
    }
    checkNear(first.structure().diagnostics()
                  .linearMomentumKgMetersPerSecond.x,
              300.0, 2.0e-10,
              "open piston: structural momentum matches the driven plate");

    const auto& control = first.controlVolumeDiagnostics();
    check(control.accepted && control.finite,
          "open piston: final control-volume ledger remains accepted");
    checkNear(control.endVolumeCubicMeters,
              12.0 + 6.0 * expectedOffset, 2.0e-14,
              "open piston: visible motion accumulates analytic chamber volume");
    check(std::abs(control.surfaceGeometryResidualCubicMeters) < 2.0e-15
              && std::abs(control.continuityResidualCubicMeters) < 2.0e-15,
          "open piston: geometry sweep and opening transport remain closed");
    check(finalFrame.couplingResiduals.fluid < 2.0e-11
              && finalFrame.couplingResiduals.tractionNewtons < 1.0e-10
              && finalFrame.couplingResiduals.interfacePowerWatts < 1.0e-10,
          "open piston: visible numerical residuals meet their budgets");
}

} // namespace

int main() {
    testVisibleOpenPistonAndDeterminism();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing open piston case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing open piston case checks passed");
    return 0;
}
