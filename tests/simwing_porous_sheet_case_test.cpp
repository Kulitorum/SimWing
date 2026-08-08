#include "porous_sheet_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <ranges>
#include <sstream>
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
          "porous-sheet: deterministic frame serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const std::string_view name) {
    const auto found = std::ranges::find(
        frame.scalarFields, name, &viewer::ScalarField::name);
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

void testCoupledPorousSheet() {
    fsi::CoupledPorousSheetCase first;
    fsi::CoupledPorousSheetCase second;
    const auto firstFrame = first.advance();
    const auto secondFrame = second.advance();
    check(serialized(firstFrame) == serialized(secondFrame)
              && first.diagnostics() == second.diagnostics(),
          "porous-sheet: independent workers accept the same first state");
    const auto& firstDiagnostics = first.diagnostics();
    check(firstDiagnostics.accepted && firstDiagnostics.finite
              && firstDiagnostics.version
                  == fsi::coupledPorousSheetDiagnosticsVersion
              && firstDiagnostics.fluidProjection.accepted
              && firstDiagnostics.porousTraction.accepted
              && firstDiagnostics.bridge.accepted,
          "porous-sheet: every fluid, transfer, and structural boundary is accepted");
    check(firstFrame.sceneChecksum == fsi::coupledPorousSheetCaseChecksum
              && firstFrame.solverCommit
                  == fsi::coupledPorousSheetCaseSolverId
              && firstFrame.step == 1
              && firstFrame.vertices.size() == 4
              && firstFrame.triangles.size() == 2,
          "porous-sheet: first frame has stable provenance and complete geometry");
    check(firstFrame.triangles.front().negativeRegionId == 10
              && firstFrame.triangles.front().positiveRegionId == 11,
          "porous-sheet: viewer geometry retains two-sided fluid ownership");
    check(firstDiagnostics.fluidVelocityAfterMetersPerSecond
              > firstDiagnostics.sheetVelocityAfterMetersPerSecond
              && firstDiagnostics.sheetVelocityAfterMetersPerSecond > 0.0,
          "porous-sheet: pump accelerates fluid through and then loads the sheet");
    check(firstDiagnostics.momentumResidualNormNewtonSeconds < 3.0e-10
              && std::abs(firstDiagnostics.energyResidualJoules) < 2.0e-10
              && firstDiagnostics.porousDissipatedEnergyJoules > 0.0,
          "porous-sheet: momentum, pump work, and porous loss close independently");
    checkNear(
        firstDiagnostics.actualSheetImpulseNewtonSeconds.x,
        firstDiagnostics.porousImpulseOnSheetNewtonSeconds.x,
        2.0e-10,
        "porous-sheet: XPBD receives exactly the accepted sheet impulse");
    checkNear(
        firstDiagnostics.actualFluidImpulseNewtonSeconds.x,
        firstDiagnostics.pressureJumpImpulseOnFluidNewtonSeconds.x,
        2.0e-10,
        "porous-sheet: fluid momentum follows its complete jump ledger");
    const auto* jump = scalarField(firstFrame, "porous.pressure_jump");
    const auto* loss = scalarField(
        firstFrame, "porous.dissipated_energy");
    const auto* pumpWork = scalarField(firstFrame, "pump.step_work");
    check(jump != nullptr
              && jump->association == viewer::FieldAssociation::Triangle
              && jump->values.size() == 2
              && jump->values.front() < 0.0
              && loss != nullptr && loss->values.front() > 0.0
              && pumpWork != nullptr && pumpWork->values.front() > 0.0,
          "porous-sheet: visible fields separate pressure jump, loss, and pump work");

    viewer::DiagnosticFrame finalFrame = firstFrame;
    constexpr std::uint64_t steps = 120;
    for (std::uint64_t step = 1; step < steps; ++step) {
        finalFrame = first.advance();
        const auto replay = second.advance();
        check(serialized(finalFrame) == serialized(replay),
              "porous-sheet: coupled replay remains deterministic");
        check(first.diagnostics().accepted
                  && first.diagnostics().momentumResidualNormNewtonSeconds
                      < 3.0e-10
                  && std::abs(first.diagnostics().energyResidualJoules)
                      < 2.0e-10,
              "porous-sheet: every accepted step retains conservation closure");
    }
    const double elapsed = steps * first.stepSettings().timeStepSeconds;
    const double expectedTotalMomentum = 9.0 * elapsed;
    const auto structureDiagnostics = first.structure().diagnostics();
    const double fluidMass = 1.2
        * static_cast<double>(first.grid().cellCount())
        * first.grid().cellVolumeCubicMeters();
    checkNear(
        structureDiagnostics.linearMomentumKgMetersPerSecond.x
            + fluidMass
                * first.diagnostics().fluidVelocityAfterMetersPerSecond,
        expectedTotalMomentum, 2.0e-8,
        "porous-sheet: cumulative system momentum equals pump impulse");
    check(finalFrame.step == steps
              && first.diagnostics().sheetPositionAfterMeters > 1.45
              && first.diagnostics().sheetPositionAfterMeters < 1.75,
          "porous-sheet: visible sheet motion stays inside the declared topology epoch");
    check(serialized(firstFrame) == serialized(secondFrame),
          "porous-sheet: old immutable frames do not alias later solver state");
}

void testCompletedTrace() {
    fsi::CoupledPorousSheetCase simulation;
    std::stringstream trace(std::ios::in | std::ios::out | std::ios::binary);
    viewer::TraceWriter writer(trace);
    check(writer.writeHeader(simulation.traceHeader()),
          "porous-sheet trace: header writes");
    constexpr std::uint64_t steps = 5;
    for (std::uint64_t step = 0; step < steps; ++step) {
        check(writer.writeFrame(simulation.advance()),
              "porous-sheet trace: accepted frame writes");
    }
    check(writer.finish(), "porous-sheet trace: completion record writes");
    trace.seekg(0);
    viewer::TraceReader reader(trace);
    viewer::TraceHeader header;
    check(reader.readHeader(header)
              && header.sceneChecksum
                  == fsi::coupledPorousSheetCaseChecksum,
          "porous-sheet trace: provenance round-trips");
    std::uint64_t frameCount = 0;
    viewer::DiagnosticFrame frame;
    for (;;) {
        const auto status = reader.readNext(frame);
        if (status == viewer::TraceReadStatus::Frame) {
            ++frameCount;
            check(frame.step == frameCount,
                  "porous-sheet trace: accepted steps remain consecutive");
        } else {
            check(status == viewer::TraceReadStatus::End,
                  "porous-sheet trace: completed replay terminates cleanly");
            break;
        }
    }
    check(frameCount == steps,
          "porous-sheet trace: replay contains every accepted frame");
}

void testTopologyBoundaryRollback() {
    fsi::CoupledPorousSheetCase simulation;
    bool rejected = false;
    for (std::size_t attempt = 0; attempt < 1000; ++attempt) {
        const auto structureBefore = simulation.structure().checkpoint();
        const auto velocityBefore = simulation.velocity();
        const auto pressureBefore = simulation.pressure();
        const auto diagnosticsBefore = simulation.diagnostics();
        try {
            static_cast<void>(simulation.advance());
        } catch (const std::exception&) {
            const auto structureAfter = simulation.structure().checkpoint();
            check(structureAfter.acceptedStepCount
                      == structureBefore.acceptedStepCount
                      && structureAfter.simulationTimeSeconds
                          == structureBefore.simulationTimeSeconds
                      && structureAfter.nodes == structureBefore.nodes
                      && structureAfter.pendingExternalForcesNewtons
                          == structureBefore.pendingExternalForcesNewtons
                      && simulation.velocity() == velocityBefore
                      && simulation.pressure() == pressureBefore
                      && simulation.diagnostics() == diagnosticsBefore,
                  "porous-sheet topology: rejected segment crossing rolls back every owner");
            rejected = true;
            break;
        }
    }
    check(rejected,
          "porous-sheet topology: fixed epoch rejects before an implicit remap");
}

void expectRestoreRejected(
    fsi::CoupledPorousSheetCase& simulation,
    const fsi::CoupledPorousSheetCheckpoint& checkpoint,
    const char* message) {
    const auto before = simulation.checkpoint();
    const auto structureBefore = simulation.structure().checkpoint();
    const auto velocityBefore = simulation.velocity();
    const auto pressureBefore = simulation.pressure();
    const auto diagnosticsBefore = simulation.diagnostics();
    bool rejected = false;
    try {
        simulation.restore(checkpoint);
    } catch (const std::exception&) {
        rejected = true;
    }
    const auto structureAfter = simulation.structure().checkpoint();
    const auto after = simulation.checkpoint();
    check(rejected
              && after.acceptedStepCount == before.acceptedStepCount
              && after.simulationTimeSeconds == before.simulationTimeSeconds
              && structureAfter.acceptedStepCount
                  == structureBefore.acceptedStepCount
              && structureAfter.nodes == structureBefore.nodes
              && simulation.velocity() == velocityBefore
              && simulation.pressure() == pressureBefore
              && simulation.diagnostics() == diagnosticsBefore,
          message);
}

void testCheckpointReplayAndValidation() {
    fsi::CoupledPorousSheetCase initialOwner;
    const auto initial = initialOwner.checkpoint();
    fsi::CoupledPorousSheetCase initialRestored;
    static_cast<void>(initialRestored.advance());
    initialRestored.restore(initial);
    check(serialized(initialOwner.advance())
              == serialized(initialRestored.advance()),
          "porous-sheet checkpoint: initial state reproduces the first frame");

    fsi::CoupledPorousSheetCase owner;
    constexpr std::uint64_t checkpointStep = 60;
    for (std::uint64_t step = 0; step < checkpointStep; ++step) {
        static_cast<void>(owner.advance());
    }
    const auto checkpoint = owner.checkpoint();
    fsi::CoupledPorousSheetCase restored;
    restored.restore(checkpoint);
    check(checkpoint.version == fsi::coupledPorousSheetCheckpointVersion
              && checkpoint.caseFingerprint
                  == fsi::coupledPorousSheetCaseFingerprint
              && checkpoint.acceptedStepCount == checkpointStep
              && restored.diagnostics() == owner.diagnostics()
              && restored.velocity() == owner.velocity()
              && restored.pressure() == owner.pressure(),
          "porous-sheet checkpoint: accepted structure, fluid, and diagnostics restore exactly");
    checkNear(
        checkpoint.simulationTimeSeconds,
        checkpointStep * owner.stepSettings().timeStepSeconds,
        2.0e-15,
        "porous-sheet checkpoint: public time retains the accepted step epoch");
    const auto expected = owner.advance();
    const auto replay = restored.advance();
    check(serialized(expected) == serialized(replay)
              && owner.diagnostics() == restored.diagnostics(),
          "porous-sheet checkpoint: restored worker reproduces the exact next frame");
    restored.restore(checkpoint);
    const auto repeated = restored.advance();
    check(serialized(expected) == serialized(repeated),
          "porous-sheet checkpoint: repeated restore replays deterministically");

    auto invalid = checkpoint;
    invalid.version += 1;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: unsupported version is transactional");
    invalid = checkpoint;
    invalid.caseFingerprint ^= 1;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: foreign case is transactional");
    invalid = checkpoint;
    invalid.acceptedStepCount += 1;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: public step must match private state");
    invalid = checkpoint;
    invalid.simulationTimeSeconds += 0.5;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: public time must match private state");
    expectRestoreRejected(
        restored, {},
        "porous-sheet checkpoint validation: empty payload is transactional");
}

} // namespace

int main() {
    testCoupledPorousSheet();
    testCompletedTrace();
    testTopologyBoundaryRollback();
    testCheckpointReplayAndValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing coupled porous-sheet check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing coupled porous-sheet checks passed");
    return 0;
}
