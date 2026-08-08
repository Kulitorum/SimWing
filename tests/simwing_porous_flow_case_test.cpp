#include "porous_flow_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ranges>
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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "porous-flow worker frame serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [&](const auto& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

void testPressureDrivenWorker() {
    fsi::PorousFlowCase first;
    fsi::PorousFlowCase second;
    const auto firstFrame = first.advance();
    const auto secondFrame = second.advance();
    const auto firstBytes = serialized(firstFrame);
    check(firstBytes == serialized(secondFrame),
          "porous-flow workers publish deterministic first frames");
    const auto header = first.traceHeader();
    check(header.sceneChecksum == fsi::porousFlowCaseChecksum
              && header.solverCommit == fsi::porousFlowCaseSolverId,
          "porous-flow worker exposes its stable trace identity");
    check(first.acceptedStepCount() == 1
              && first.simulationTimeSeconds()
                  == first.stepSettings().timeStepSeconds
              && first.diagnostics().converged
              && first.flowDiagnostics().accepted
              && first.flowDiagnostics().projectedAdvection.firstProjection
                     .pressureJumpFaceCount == 24
              && first.flowDiagnostics().projectedAdvection.secondProjection
                     .pressureJumpFaceCount == 24
              && first.pressureJumps().faceCount() == 24,
          "porous-flow worker commits its plug and complete fluid step together");
    check(firstFrame.vertices.size() == 288
              && firstFrame.triangles.size() == 48,
          "porous-flow worker publishes cell samples and two interface planes");
    check(first.plugDiagnostics().accepted
              && std::abs(first.plugDiagnostics()
                              .momentumResidualNewtonSeconds) < 4.0e-15
              && std::abs(first.plugDiagnostics().energyResidualJoules)
                  < 4.0e-15
              && first.plugDiagnostics().porousDissipationJoules > 0.0,
          "porous-flow worker closes impulse, energy, and dissipation");

    const double speed = first.plugDiagnostics()
        .velocityAfterMetersPerSecond;
    check(speed > 0.0
              && std::ranges::all_of(
                  first.velocity().xFaces(),
                  [&](const double sample) {
                      return std::abs(sample - speed) < 2.0e-13;
                  })
              && std::ranges::all_of(
                  first.velocity().yFaces(),
                  [](const double sample) {
                      return std::abs(sample) < 2.0e-13;
                  })
              && std::ranges::all_of(
                  first.velocity().zFaces(),
                  [](const double sample) {
                      return std::abs(sample) < 2.0e-13;
                  }),
          "porous-flow grid retains the accepted uniform plug velocity");
    const double pressureDrop = first.plugDiagnostics()
        .endpointPressureDropPascals;
    const auto counts = first.grid().cellCounts();
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i >= 4 && i < 12
                    ? -0.5 * pressureDrop : 0.5 * pressureDrop;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(first.pressure().values()
                                 [first.grid().cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 2.0e-12,
          "porous-flow worker retains the analytic endpoint pressure loss");
    const auto* frameSpeed = scalarField(
        firstFrame, "porous normal velocity");
    const auto* frameDrop = scalarField(
        firstFrame, "porous endpoint pressure drop");
    const auto* frameDissipation = scalarField(
        firstFrame, "porous dissipation");
    const auto* frameDrive = scalarField(
        firstFrame, "driving pressure rise");
    check(frameSpeed != nullptr
              && frameSpeed->association
                  == viewer::FieldAssociation::Global
              && frameSpeed->values == std::vector<double>{speed}
              && frameDrop != nullptr
              && frameDrop->values == std::vector<double>{pressureDrop}
              && frameDissipation != nullptr
              && frameDissipation->values.front()
                  == first.plugDiagnostics().porousDissipationJoules
              && frameDrive != nullptr
              && frameDrive->values
                  == std::vector<double>{
                      first.plugSettings().drivingPressureRisePascals},
          "porous-flow frame exposes its constitutive state and loss");

    for (std::size_t step = 1; step < 240; ++step) {
        static_cast<void>(first.advance());
        static_cast<void>(second.advance());
    }
    const double steady = 0.5 * (std::sqrt(56.0) - 4.0);
    checkNear(first.plugDiagnostics().velocityAfterMetersPerSecond,
              steady, 2.0e-14,
              "porous-flow worker reaches the analytic steady velocity");
    checkNear(first.plugDiagnostics().endpointPressureDropPascals,
              first.plugSettings().drivingPressureRisePascals,
              3.0e-12,
              "porous-flow worker reaches the analytic steady pressure loss");
    check(serialized(first.advance()) == serialized(second.advance()),
          "porous-flow worker replay remains deterministic at steady state");
    check(serialized(firstFrame) == firstBytes,
          "porous-flow frames remain owning after later solver advances");
}

} // namespace

int main() {
    testPressureDrivenWorker();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing porous-flow case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing porous-flow case checks passed");
    return 0;
}
