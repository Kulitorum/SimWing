#include "pressure_jump_case.h"
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
          "pressure-jump worker frame serializes");
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

void testStaticProjectionAndDeterminism() {
    fsi::PressureJumpCase first;
    fsi::PressureJumpCase second;
    const auto firstFrame = first.advance();
    const auto secondFrame = second.advance();
    check(serialized(firstFrame) == serialized(secondFrame),
          "pressure-jump worker first accepted frame is deterministic");
    const auto header = first.traceHeader();
    check(header.sceneChecksum == fsi::pressureJumpCaseChecksum
              && header.solverCommit == fsi::pressureJumpCaseSolverId,
          "pressure-jump worker exposes its stable trace identity");
    check(first.acceptedStepCount() == 1
              && first.simulationTimeSeconds()
                  == first.stepSettings().timeStepSeconds
              && first.diagnostics().converged
              && first.diagnostics().pressureJumpFaceCount == 48,
          "pressure-jump worker commits one complete split-slab projection");
    check(firstFrame.vertices.size() == 384
              && firstFrame.triangles.size() == 96,
          "pressure-jump worker publishes cell samples and 48 layered quads");
    check(firstFrame.couplingResiduals.fluid < 1.0e-12,
          "pressure-jump worker publishes a divergence-free accepted state");

    double maximumVelocity = 0.0;
    for (std::size_t index = 0; index < first.grid().cellCount(); ++index) {
        maximumVelocity = std::max({
            maximumVelocity,
            std::abs(first.velocity().xFaces()[index]),
            std::abs(first.velocity().yFaces()[index]),
            std::abs(first.velocity().zFaces()[index]),
        });
    }
    check(maximumVelocity < 2.0e-13,
          "pressure-jump worker creates no spurious flow");
    const auto counts = first.grid().cellCounts();
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i >= 4 && i < 12
                    ? 125.0 : -125.0;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(first.pressure().values()
                                 [first.grid().cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 2.0e-12,
          "pressure-jump worker retains the exact two-level slab pressure");
    const auto* jumps = scalarField(firstFrame, "pressure jump");
    check(jumps != nullptr
              && jumps->association == viewer::FieldAssociation::Triangle
              && jumps->values.size() == firstFrame.triangles.size()
              && jumps->values[0] == 100.0
              && jumps->values[2] == 150.0,
          "pressure-jump worker exposes both subcell transition layers");

    const auto acceptedPressure = first.pressure();
    const auto acceptedVelocity = first.velocity();
    const auto continued = first.advance();
    const auto replay = second.advance();
    check(serialized(continued) == serialized(replay)
              && first.pressure() == acceptedPressure
              && first.velocity() == acceptedVelocity,
          "pressure-jump worker reprojects the static state bit-identically");
    checkNear(first.simulationTimeSeconds(),
              2.0 * first.stepSettings().timeStepSeconds, 0.0,
              "pressure-jump worker advances absolute diagnostic time");
}

} // namespace

int main() {
    testStaticProjectionAndDeterminism();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing pressure-jump case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing pressure-jump case checks passed");
    return 0;
}
