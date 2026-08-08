#include "porous_sheet_case.h"
#include "porous_sheet_checkpoint_persistence.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <exception>
#include <ranges>
#include <sstream>
#include <string>
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

void testTopologyRebaseAndCollisionRollback() {
    fsi::CoupledPorousSheetCase simulation;
    fsi::CoupledPorousSheetCase replay;
    viewer::DiagnosticFrame rebaseFrame;
    viewer::DiagnosticFrame replayFrame;
    for (std::size_t attempt = 0; attempt < 400
         && simulation.topologyRebaseCount() == 0; ++attempt) {
        rebaseFrame = simulation.advance();
        replayFrame = replay.advance();
    }
    check(simulation.topologyRebaseCount() == 1
              && simulation.porousFaceCoordinate() == 4
              && simulation.diagnostics().topologyRebasedThisStep
              && simulation.diagnostics().topologyRebaseCount == 1
              && simulation.diagnostics().porousFaceCoordinate == 4
              && simulation.diagnostics()
                     .bridge.mapping.gridPlaneCoordinateMeters
                  == simulation.diagnostics()
                         .sheetPositionAtConstitutiveTimeMeters
              && simulation.diagnostics()
                     .bridge.mapping.physicalPlaneCoordinateMeters > 1.75,
          "porous-sheet topology: first dual-cell crossing commits an explicit rebase");
    check(serialized(rebaseFrame) == serialized(replayFrame),
          "porous-sheet topology: the accepted crossing frame replays bit-for-bit");

    const auto rebasedCheckpoint = simulation.checkpoint();
    fsi::CoupledPorousSheetCase restored;
    restored.restore(rebasedCheckpoint);
    check(rebasedCheckpoint.topologyRebaseCount == 1
              && rebasedCheckpoint.porousFaceCoordinate == 4
              && serialized(simulation.advance())
                  == serialized(restored.advance()),
          "porous-sheet topology: rebased checkpoint resumes the exact next frame");

    bool rejected = false;
    for (std::size_t attempt = 0; attempt < 1000; ++attempt) {
        const auto structureBefore = simulation.structure().checkpoint();
        const auto velocityBefore = simulation.velocity();
        const auto pressureBefore = simulation.pressure();
        const auto diagnosticsBefore = simulation.diagnostics();
        const auto checkpointBefore = simulation.checkpoint();
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
                      && simulation.diagnostics() == diagnosticsBefore
                      && simulation.topologyRebaseCount()
                          == checkpointBefore.topologyRebaseCount
                      && simulation.porousFaceCoordinate()
                          == checkpointBefore.porousFaceCoordinate,
                  "porous-sheet topology: pump collision rolls back every owner and epoch");
            rejected = true;
            break;
        }
    }
    check(rejected,
          "porous-sheet topology: later pump-surface collision is rejected explicitly");
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
              && after.topologyRebaseCount == before.topologyRebaseCount
              && after.porousFaceCoordinate == before.porousFaceCoordinate
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
    constexpr std::uint64_t checkpointStep = 350;
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
              && checkpoint.topologyRebaseCount == 1
              && checkpoint.porousFaceCoordinate == 4
              && restored.topologyRebaseCount() == 1
              && restored.porousFaceCoordinate() == 4
              && restored.diagnostics() == owner.diagnostics()
              && restored.velocity() == owner.velocity()
              && restored.pressure() == owner.pressure(),
          "porous-sheet checkpoint: accepted structure, fluid, and diagnostics restore exactly");
    checkNear(
        checkpoint.simulationTimeSeconds,
        checkpointStep * owner.stepSettings().timeStepSeconds,
        1.0e-14,
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
    invalid = checkpoint;
    invalid.topologyRebaseCount = 0;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: topology epoch must match state");
    invalid = checkpoint;
    invalid.porousFaceCoordinate = 3;
    expectRestoreRejected(
        restored, invalid,
        "porous-sheet checkpoint validation: topology face must match state");
    expectRestoreRejected(
        restored, {},
        "porous-sheet checkpoint validation: empty payload is transactional");
}

std::uint64_t readU64(
    const std::vector<std::uint8_t>& bytes,
    const std::size_t offset) {
    std::uint64_t result = 0;
    for (unsigned shift = 0; shift < 64; shift += 8) {
        result |= static_cast<std::uint64_t>(
            bytes[offset + shift / 8]) << shift;
    }
    return result;
}

void writeU64(std::vector<std::uint8_t>& bytes,
              const std::size_t offset,
              const std::uint64_t value) {
    for (unsigned shift = 0; shift < 64; shift += 8) {
        bytes[offset + shift / 8] =
            static_cast<std::uint8_t>(value >> shift);
    }
}

std::uint64_t payloadChecksum(
    const std::span<const std::uint8_t> bytes) {
    std::uint64_t result = 14695981039346656037ULL;
    for (const std::uint8_t byte : bytes) {
        result ^= byte;
        result *= 1099511628211ULL;
    }
    return result;
}

void recomputeEnvelopeChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t prefix = 16;
    constexpr std::size_t checksumBytes = 8;
    const std::size_t payloadSize = static_cast<std::size_t>(
        readU64(bytes, 8));
    check(bytes.size() == prefix + payloadSize + checksumBytes,
          "porous-sheet persistence fixture: envelope size is canonical");
    writeU64(
        bytes, prefix + payloadSize,
        payloadChecksum(std::span<const std::uint8_t>(
            bytes.data() + prefix, payloadSize)));
}

std::vector<std::uint8_t> encodeCheckpoint(
    const fsi::CoupledPorousSheetCase& owner,
    const fsi::CoupledPorousSheetCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    fsi::CoupledPorousSheetCheckpointPersistenceError error;
    check(fsi::serializeCoupledPorousSheetCheckpoint(
              owner, checkpoint, bytes, &error),
          "porous-sheet persistence: valid checkpoint serializes");
    return bytes;
}

void expectDecodeRejected(
    const std::vector<std::uint8_t>& bytes,
    const fsi::CoupledPorousSheetCase& owner,
    fsi::CoupledPorousSheetCheckpoint& destination,
    const fsi::CoupledPorousSheetCheckpointPersistenceErrorCode expectedCode,
    const char* message,
    const fsi::CoupledPorousSheetCheckpointPersistenceLimits& limits = {}) {
    const auto before = destination;
    fsi::CoupledPorousSheetCheckpointPersistenceError error;
    const bool decoded = fsi::deserializeCoupledPorousSheetCheckpoint(
        bytes, owner, destination, &error, limits);
    fsi::CoupledPorousSheetCase expected;
    fsi::CoupledPorousSheetCase actual;
    expected.restore(before);
    actual.restore(destination);
    check(!decoded && error.code == expectedCode
              && destination.version == before.version
              && destination.caseFingerprint == before.caseFingerprint
              && destination.acceptedStepCount == before.acceptedStepCount
              && destination.simulationTimeSeconds
                  == before.simulationTimeSeconds
              && destination.topologyRebaseCount
                  == before.topologyRebaseCount
              && destination.porousFaceCoordinate
                  == before.porousFaceCoordinate
              && serialized(expected.advance())
                  == serialized(actual.advance()),
          message);
}

void testPersistentCheckpoint() {
    fsi::CoupledPorousSheetCase initialOwner;
    const auto initialBytes = encodeCheckpoint(
        initialOwner, initialOwner.checkpoint());
    fsi::CoupledPorousSheetCheckpoint decodedInitial;
    fsi::CoupledPorousSheetCheckpointPersistenceError error;
    check(fsi::deserializeCoupledPorousSheetCheckpoint(
              initialBytes, initialOwner, decodedInitial, &error),
          "porous-sheet persistence: initial checkpoint decodes");
    fsi::CoupledPorousSheetCase initialReplay;
    initialReplay.restore(decodedInitial);
    check(serialized(initialOwner.advance())
              == serialized(initialReplay.advance()),
          "porous-sheet persistence: decoded initial state reproduces the first frame");

    fsi::CoupledPorousSheetCase owner;
    constexpr std::uint64_t stepCount = 350;
    for (std::uint64_t step = 0; step < stepCount; ++step) {
        static_cast<void>(owner.advance());
    }
    const auto checkpoint = owner.checkpoint();
    const auto first = encodeCheckpoint(owner, checkpoint);
    const auto second = encodeCheckpoint(owner, checkpoint);
    check(first == second && first.size() > 256,
          "porous-sheet persistence: repeated encoding is byte-deterministic");

    fsi::CoupledPorousSheetCase destinationOwner;
    for (std::size_t step = 0; step < 5; ++step) {
        static_cast<void>(destinationOwner.advance());
    }
    auto destination = destinationOwner.checkpoint();
    check(fsi::deserializeCoupledPorousSheetCheckpoint(
              first, destinationOwner, destination, &error)
              && destination.acceptedStepCount == stepCount
              && destination.topologyRebaseCount == 1
              && destination.porousFaceCoordinate == 4,
          "porous-sheet persistence: accepted checkpoint decodes transactionally");
    const auto reencoded = encodeCheckpoint(destinationOwner, destination);
    check(reencoded == first,
          "porous-sheet persistence: decode and re-encode preserve exact bytes");
    fsi::CoupledPorousSheetCase replay;
    replay.restore(destination);
    check(serialized(owner.advance()) == serialized(replay.advance()),
          "porous-sheet persistence: decoded accepted state reproduces the exact next frame");

    auto invalid = first;
    invalid[0] ^= 1;
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidMagic,
        "porous-sheet persistence validation: bad magic is transactional");
    invalid = first;
    invalid[4] += 1;
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::UnsupportedVersion,
        "porous-sheet persistence validation: protocol version is transactional");
    invalid = first;
    invalid[6] = 1;
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
        "porous-sheet persistence validation: reserved bits are transactional");
    invalid = first;
    invalid.pop_back();
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::Truncated,
        "porous-sheet persistence validation: truncation is transactional");
    invalid = first;
    invalid.push_back(0);
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::TrailingData,
        "porous-sheet persistence validation: trailing bytes are transactional");
    invalid = first;
    invalid[32] ^= 1;
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::ChecksumMismatch,
        "porous-sheet persistence validation: checksum damage is transactional");

    constexpr std::size_t payloadOffset = 16;
    constexpr std::size_t topologyRebaseOffsetInPayload = 28;
    constexpr std::size_t porousFaceOffsetInPayload = 36;
    invalid = first;
    invalid[payloadOffset + topologyRebaseOffsetInPayload] ^= 1;
    recomputeEnvelopeChecksum(invalid);
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
        "porous-sheet persistence validation: topology epoch corruption cannot evade replay");
    invalid = first;
    invalid[payloadOffset + porousFaceOffsetInPayload] ^= 1;
    recomputeEnvelopeChecksum(invalid);
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
        "porous-sheet persistence validation: topology face corruption cannot evade replay");

    fsi::CoupledPorousSheetCheckpointPersistenceLimits limits;
    limits.maximumEncodedBytes = first.size() - 1;
    expectDecodeRejected(
        first, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
        "porous-sheet persistence validation: byte limit is transactional",
        limits);
    limits = {};
    limits.maximumScalarSamples = 191;
    expectDecodeRejected(
        first, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
        "porous-sheet persistence validation: field limit is transactional",
        limits);
    limits = {};
    limits.maximumReplaySteps = stepCount - 1;
    expectDecodeRejected(
        first, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
        "porous-sheet persistence validation: replay limit is transactional",
        limits);
    limits = {};
    limits.maximumScalarSamples = 191;
    std::vector<std::uint8_t> preservedOutput{7, 8, 9};
    check(!fsi::serializeCoupledPorousSheetCheckpoint(
              owner, checkpoint, preservedOutput, &error, limits)
              && error.code
                  == fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded
              && preservedOutput == std::vector<std::uint8_t>({7, 8, 9}),
          "porous-sheet persistence validation: failed encoding preserves caller output");

    invalid = first;
    constexpr std::size_t structureSizeOffsetInPayload = 124;
    constexpr std::size_t structureOffsetInPayload = 132;
    const std::size_t structureSize = static_cast<std::size_t>(readU64(
        invalid, payloadOffset + structureSizeOffsetInPayload));
    check(structureSize > 32,
          "porous-sheet persistence fixture: nested Structure envelope is present");
    invalid[payloadOffset + structureOffsetInPayload + 24] ^= 1;
    recomputeEnvelopeChecksum(invalid);
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
        "porous-sheet persistence validation: nested Structure corruption is transactional");
    invalid = first;
    invalid[payloadOffset + structureOffsetInPayload + structureSize] ^= 1;
    recomputeEnvelopeChecksum(invalid);
    expectDecodeRejected(
        invalid, destinationOwner, destination,
        fsi::CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
        "porous-sheet persistence validation: recomputed field corruption cannot evade replay");
}

void testTerminalSafePointPersistence() {
    fsi::CoupledPorousSheetCase owner;
    std::string collisionError;
    for (std::uint64_t attempt = 0;
         attempt < 1000 && collisionError.empty(); ++attempt) {
        try {
            static_cast<void>(owner.advance());
        } catch (const std::exception& exception) {
            collisionError = exception.what();
        }
    }
    const auto terminal = owner.checkpoint();
    check(collisionError
              == "coupled porous sheet reached the pump-surface topology"
              && terminal.acceptedStepCount > 350
              && terminal.topologyRebaseCount == 1
              && terminal.porousFaceCoordinate == 4,
          "porous-sheet terminal checkpoint: collision leaves an accepted safe point");

    const auto bytes = encodeCheckpoint(owner, terminal);
    fsi::CoupledPorousSheetCheckpoint decoded;
    fsi::CoupledPorousSheetCheckpointPersistenceError error;
    check(fsi::deserializeCoupledPorousSheetCheckpoint(
              bytes, owner, decoded, &error)
              && decoded.acceptedStepCount == terminal.acceptedStepCount
              && decoded.simulationTimeSeconds
                  == terminal.simulationTimeSeconds
              && decoded.topologyRebaseCount == terminal.topologyRebaseCount
              && decoded.porousFaceCoordinate
                  == terminal.porousFaceCoordinate,
          "porous-sheet terminal checkpoint: accepted endpoint persists exactly");

    const auto verifyRepeatedCollision = [](const auto& checkpoint,
                                            const char* message) {
        fsi::CoupledPorousSheetCase replay;
        replay.restore(checkpoint);
        const auto structureBefore = replay.structure().checkpoint();
        const auto velocityBefore = replay.velocity();
        const auto pressureBefore = replay.pressure();
        const auto diagnosticsBefore = replay.diagnostics();
        bool rejected = false;
        try {
            static_cast<void>(replay.advance());
        } catch (const std::exception& exception) {
            rejected = std::string(exception.what())
                == "coupled porous sheet reached the pump-surface topology";
        }
        const auto structureAfter = replay.structure().checkpoint();
        check(rejected
                  && structureAfter.acceptedStepCount
                      == structureBefore.acceptedStepCount
                  && structureAfter.simulationTimeSeconds
                      == structureBefore.simulationTimeSeconds
                  && structureAfter.nodes == structureBefore.nodes
                  && replay.velocity() == velocityBefore
                  && replay.pressure() == pressureBefore
                  && replay.diagnostics() == diagnosticsBefore,
              message);
    };
    verifyRepeatedCollision(
        terminal,
        "porous-sheet terminal checkpoint: in-memory restore repeats collision transactionally");
    verifyRepeatedCollision(
        decoded,
        "porous-sheet terminal checkpoint: persisted restore repeats collision transactionally");
    check(encodeCheckpoint(owner, decoded) == bytes,
          "porous-sheet terminal checkpoint: decode and re-encode preserve bytes");
}

} // namespace

int main() {
    testCoupledPorousSheet();
    testCompletedTrace();
    testTopologyRebaseAndCollisionRollback();
    testCheckpointReplayAndValidation();
    testPersistentCheckpoint();
    testTerminalSafePointPersistence();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing coupled porous-sheet check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing coupled porous-sheet checks passed");
    return 0;
}
