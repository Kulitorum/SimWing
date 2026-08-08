#include "moving_porous_flow_case.h"
#include "viewer_protocol.h"

#include <cmath>
#include <cstdio>
#include <exception>
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

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected, message);
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

void testCheckpointReplayAndValidation() {
    fsi::MovingPorousFlowCase initialOwner;
    const auto initial = initialOwner.checkpoint();
    const auto expectedFirst = initialOwner.advance();
    fsi::MovingPorousFlowCase initialReplay;
    initialReplay.restore(initial);
    check(initial.version == fsi::movingPorousFlowCaseCheckpointVersion
              && initial.caseDefinitionFingerprint
                  == fsi::movingPorousFlowCaseDefinitionFingerprint
              && initial.cellCounts == initialReplay.grid().cellCounts()
              && initial.scalarSampleCount
                  == initialReplay.grid().cellCount()
              && initial.pressureJumpCount == 0
              && initial.acceptedStepCount == 0
              && serialized(expectedFirst)
                  == serialized(initialReplay.advance()),
          "moving porous-flow initial checkpoint replays its first wrap");

    fsi::MovingPorousFlowCase owner;
    constexpr std::uint64_t savedSteps = 37;
    for (std::uint64_t step = 0; step < savedSteps; ++step) {
        static_cast<void>(owner.advance());
    }
    const auto saved = owner.checkpoint();
    const auto savedVelocity = owner.velocity();
    const auto savedPressure = owner.pressure();
    const auto savedJumps = owner.pressureJumps();
    const auto savedDiagnostics = owner.diagnostics();
    const auto expectedNext = owner.advance();
    for (std::size_t step = 0; step < 5; ++step) {
        static_cast<void>(owner.advance());
    }
    owner.restore(saved);
    check(owner.velocity() == savedVelocity
              && owner.pressure() == savedPressure
              && owner.pressureJumps() == savedJumps
              && owner.diagnostics() == savedDiagnostics
              && owner.porousTopology() == saved.porousTopology
              && owner.topologyRebaseCount()
                  == saved.topologyRebaseCount
              && owner.sheetPositionMeters()
                  == saved.sheetPositionMeters
              && owner.sheetVelocityMetersPerSecond()
                  == saved.sheetVelocityMetersPerSecond
              && owner.acceptedStepCount() == savedSteps
              && owner.simulationTimeSeconds()
                  == saved.simulationTimeSeconds,
          "moving porous-flow checkpoint restores every accepted field and epoch");
    check(serialized(owner.advance()) == serialized(expectedNext),
          "moving porous-flow checkpoint continuation replays bit-for-bit");
    fsi::MovingPorousFlowCase rebuilt;
    rebuilt.restore(saved);
    check(serialized(rebuilt.advance()) == serialized(expectedNext),
          "rebuilt moving porous-flow worker resumes bit-for-bit");

    fsi::MovingPorousFlowCase wrappedOwner;
    for (std::uint64_t step = 0; step < 101; ++step) {
        static_cast<void>(wrappedOwner.advance());
    }
    const auto wrapped = wrappedOwner.checkpoint();
    const auto expectedWrappedNext = wrappedOwner.advance();
    fsi::MovingPorousFlowCase wrappedReplay;
    wrappedReplay.restore(wrapped);
    check(wrapped.topologyRebaseCount == 5
              && wrapped.porousTopology.faceCoordinate == 0
              && wrapped.porousTopology.periodicImage == 2
              && serialized(wrappedReplay.advance())
                  == serialized(expectedWrappedNext),
          "moving porous-flow checkpoint crosses the second wrap exactly");

    const auto beforeVelocity = owner.velocity();
    const auto beforePressure = owner.pressure();
    const auto beforeJumps = owner.pressureJumps();
    const auto beforeDiagnostics = owner.diagnostics();
    const auto beforeTopology = owner.porousTopology();
    const auto beforeStep = owner.acceptedStepCount();
    const auto beforeTime = owner.simulationTimeSeconds();
    const auto beforePosition = owner.sheetPositionMeters();
    const auto beforeRebases = owner.topologyRebaseCount();
    auto wrongVersion = saved;
    ++wrongVersion.version;
    expectRejected(
        [&] { owner.restore(wrongVersion); },
        "moving porous-flow checkpoint rejects unsupported versions");
    auto wrongDefinition = saved;
    ++wrongDefinition.caseDefinitionFingerprint;
    expectRejected(
        [&] { owner.restore(wrongDefinition); },
        "moving porous-flow checkpoint rejects foreign definitions");
    auto wrongGrid = saved;
    ++wrongGrid.cellCounts.x;
    expectRejected(
        [&] { owner.restore(wrongGrid); },
        "moving porous-flow checkpoint rejects edited grid metadata");
    auto wrongSamples = saved;
    ++wrongSamples.scalarSampleCount;
    expectRejected(
        [&] { owner.restore(wrongSamples); },
        "moving porous-flow checkpoint rejects edited sample counts");
    auto wrongJumps = saved;
    ++wrongJumps.pressureJumpCount;
    expectRejected(
        [&] { owner.restore(wrongJumps); },
        "moving porous-flow checkpoint rejects edited jump counts");
    auto wrongStep = saved;
    ++wrongStep.acceptedStepCount;
    expectRejected(
        [&] { owner.restore(wrongStep); },
        "moving porous-flow checkpoint rejects edited step metadata");
    auto wrongTime = saved;
    wrongTime.simulationTimeSeconds += 0.1;
    expectRejected(
        [&] { owner.restore(wrongTime); },
        "moving porous-flow checkpoint rejects edited time metadata");
    auto wrongPosition = saved;
    wrongPosition.sheetPositionMeters += 0.01;
    expectRejected(
        [&] { owner.restore(wrongPosition); },
        "moving porous-flow checkpoint rejects edited sheet position");
    auto wrongSheetVelocity = saved;
    wrongSheetVelocity.sheetVelocityMetersPerSecond = -0.4;
    expectRejected(
        [&] { owner.restore(wrongSheetVelocity); },
        "moving porous-flow checkpoint rejects edited sheet velocity");
    auto wrongRebaseCount = saved;
    ++wrongRebaseCount.topologyRebaseCount;
    expectRejected(
        [&] { owner.restore(wrongRebaseCount); },
        "moving porous-flow checkpoint rejects edited rebase count");
    auto wrongTopology = saved;
    wrongTopology.porousTopology.axis = fsi::fluid::GridFaceAxis::Y;
    expectRejected(
        [&] { owner.restore(wrongTopology); },
        "moving porous-flow checkpoint rejects edited topology");
    const fsi::MovingPorousFlowCaseCheckpoint empty;
    expectRejected(
        [&] { owner.restore(empty); },
        "moving porous-flow checkpoint rejects a missing private payload");
    check(owner.velocity() == beforeVelocity
              && owner.pressure() == beforePressure
              && owner.pressureJumps() == beforeJumps
              && owner.diagnostics() == beforeDiagnostics
              && owner.porousTopology() == beforeTopology
              && owner.acceptedStepCount() == beforeStep
              && owner.simulationTimeSeconds() == beforeTime
              && owner.sheetPositionMeters() == beforePosition
              && owner.topologyRebaseCount() == beforeRebases,
          "rejected moving porous-flow checkpoints leave state untouched");
}

} // namespace

int main() {
    testMovingPorousWorker();
    testCheckpointReplayAndValidation();
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
