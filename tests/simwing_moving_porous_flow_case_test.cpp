#include "moving_porous_flow_case.h"
#include "viewer_protocol.h"

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

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "moving porous-flow frame serializes");
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

void testMovingPorousWorker() {
    fsi::MovingPorousFlowCase first;
    fsi::MovingPorousFlowCase second;
    const auto firstFrame = first.advance();
    const auto secondFrame = second.advance();
    const auto firstBytes = serialized(firstFrame);
    check(firstBytes == serialized(secondFrame),
          "moving porous-flow workers publish deterministic first frames");
    const auto header = first.traceHeader();
    check(header.sceneChecksum == fsi::movingPorousFlowCaseChecksum
              && header.solverCommit
                  == fsi::movingPorousFlowCaseSolverId,
          "moving porous-flow worker exposes stable trace identity");
    check(first.acceptedStepCount() == 1
              && std::abs(first.simulationTimeSeconds() - 0.1)
                  < 1.0e-15
              && std::abs(first.sheetPositionMeters() - 3.52)
                  < 1.0e-15
              && first.topologyRebaseCount() == 1
              && first.porousTopology().faceCoordinate == 0
              && first.porousTopology().periodicImage == 1,
          "moving porous-flow first macro-step crosses the positive wrap");
    const auto& diagnostics = first.diagnostics();
    check(diagnostics.accepted && diagnostics.finite
              && diagnostics.firstHalfSheet.topology.faceCoordinate == 3
              && diagnostics.firstHalfSheet.topology.periodicImage == 0
              && diagnostics.secondHalfSheet.topology.faceCoordinate == 0
              && diagnostics.secondHalfSheet.topology.periodicImage == 1
              && std::abs(diagnostics.kinematicResidualMeters)
                  <= diagnostics.kinematicToleranceMeters
              && diagnostics.flow.firstHalfPorous.accepted
              && diagnostics.flow.bulkFlow.accepted
              && diagnostics.flow.secondHalfPorous.accepted
              && diagnostics.flow.porousDissipationJoules > 0.0
              && diagnostics.flow.momentumResidualNormNewtonSeconds
                  < 2.0e-10,
          "moving porous-flow worker retains both accepted stages and ledgers");
    check(first.pressureJumps().faceCount() == 12
              && firstFrame.vertices.size() == 72
              && firstFrame.triangles.size() == 24,
          "moving porous-flow frame owns the fluid and both crossing planes");
    const auto* position = scalarField(firstFrame, "sheet position");
    const auto* secondFace = scalarField(
        firstFrame, "second porous face");
    const auto* secondImage = scalarField(
        firstFrame, "second porous image");
    const auto* residual = scalarField(
        firstFrame, "porous kinematic residual");
    check(position != nullptr
              && position->values
                  == std::vector<double>{first.sheetPositionMeters()}
              && secondFace != nullptr
              && secondFace->values == std::vector<double>{0.0}
              && secondImage != nullptr
              && secondImage->values == std::vector<double>{1.0}
              && residual != nullptr
              && residual->values
                  == std::vector<double>{
                      diagnostics.kinematicResidualMeters},
          "moving porous-flow frame exposes unwrapped topology and kinematics");

    for (std::size_t step = 1; step < 101; ++step) {
        static_cast<void>(first.advance());
        static_cast<void>(second.advance());
    }
    check(first.acceptedStepCount() == 101
              && first.topologyRebaseCount() == 5
              && first.porousTopology().faceCoordinate == 0
              && first.porousTopology().periodicImage == 2
              && std::abs(first.sheetPositionMeters() - 7.52)
                  < 2.0e-14,
          "moving porous-flow worker carries all epochs through a second wrap");
    check(serialized(first.advance()) == serialized(second.advance()),
          "moving porous-flow worker remains deterministic after two wraps");
    check(serialized(firstFrame) == firstBytes,
          "moving porous-flow frames remain owning after later advances");
}

} // namespace

int main() {
    testMovingPorousWorker();
    if (failures != 0) {
        std::fprintf(
            stderr,
            "%d SimWing moving porous-flow case check(s) failed\n",
            failures);
        return 1;
    }
    std::puts("all SimWing moving porous-flow case checks passed");
    return 0;
}
