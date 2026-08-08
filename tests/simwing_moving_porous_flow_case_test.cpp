#include "moving_porous_flow_case.h"
#include "moving_porous_flow_checkpoint_persistence.h"
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

std::vector<std::uint8_t> encodedCheckpoint(
    const fsi::MovingPorousFlowCaseCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    fsi::MovingPorousFlowCaseCheckpointPersistenceError error;
    check(fsi::serializeMovingPorousFlowCaseCheckpoint(
              checkpoint, bytes, &error),
          "moving porous-flow checkpoint serializes");
    check(!error,
          "successful moving porous-flow checkpoint encoding clears error");
    return bytes;
}

void refreshCheckpointChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 16;
    constexpr std::size_t checksumBytes = 8;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t checksum = offsetBasis;
    for (std::size_t index = envelopeBytes;
         index < bytes.size() - checksumBytes; ++index) {
        checksum ^= bytes[index];
        checksum *= prime;
    }
    for (unsigned shift = 0; shift < 64; shift += 8) {
        bytes[bytes.size() - checksumBytes + shift / 8] =
            static_cast<std::uint8_t>(checksum >> shift);
    }
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

void testPersistentCheckpoint() {
    fsi::MovingPorousFlowCase owner;
    constexpr std::uint64_t savedSteps = 37;
    for (std::uint64_t step = 0; step < savedSteps; ++step) {
        static_cast<void>(owner.advance());
    }
    const auto saved = owner.checkpoint();
    const auto firstEncoding = encodedCheckpoint(saved);
    check(!firstEncoding.empty()
              && firstEncoding == encodedCheckpoint(saved),
          "moving porous-flow checkpoint encoding is deterministic");

    fsi::MovingPorousFlowCaseCheckpointPersistenceError error;
    fsi::MovingPorousFlowCaseCheckpoint decoded;
    check(fsi::deserializeMovingPorousFlowCaseCheckpoint(
              firstEncoding, decoded, &error)
              && !error
              && encodedCheckpoint(decoded) == firstEncoding,
          "moving porous-flow checkpoint decode and re-encode preserves bytes");
    fsi::MovingPorousFlowCase replay;
    replay.restore(decoded);
    const auto expectedNext = owner.advance();
    check(decoded.acceptedStepCount == savedSteps
              && decoded.simulationTimeSeconds
                  == saved.simulationTimeSeconds
              && decoded.sheetPositionMeters == saved.sheetPositionMeters
              && decoded.topologyRebaseCount == saved.topologyRebaseCount
              && decoded.porousTopology == saved.porousTopology
              && serialized(replay.advance()) == serialized(expectedNext),
          "decoded moving porous-flow checkpoint continues bit-for-bit");

    fsi::MovingPorousFlowCase initialOwner;
    const auto initialBytes = encodedCheckpoint(initialOwner.checkpoint());
    fsi::MovingPorousFlowCaseCheckpoint decodedInitial;
    check(fsi::deserializeMovingPorousFlowCaseCheckpoint(
              initialBytes, decodedInitial, &error),
          "initial moving porous-flow checkpoint decodes");
    const auto expectedFirst = initialOwner.advance();
    fsi::MovingPorousFlowCase initialReplay;
    initialReplay.restore(decodedInitial);
    check(serialized(initialReplay.advance()) == serialized(expectedFirst),
          "persisted initial moving porous-flow state replays the first wrap");

    fsi::MovingPorousFlowCase wrappedOwner;
    for (std::uint64_t step = 0; step < 101; ++step) {
        static_cast<void>(wrappedOwner.advance());
    }
    const auto wrappedBytes = encodedCheckpoint(wrappedOwner.checkpoint());
    fsi::MovingPorousFlowCaseCheckpoint decodedWrapped;
    check(fsi::deserializeMovingPorousFlowCaseCheckpoint(
              wrappedBytes, decodedWrapped, &error),
          "second-wrap moving porous-flow checkpoint decodes");
    const auto expectedWrappedNext = wrappedOwner.advance();
    fsi::MovingPorousFlowCase wrappedReplay;
    wrappedReplay.restore(decodedWrapped);
    check(decodedWrapped.topologyRebaseCount == 5
              && decodedWrapped.porousTopology.faceCoordinate == 0
              && decodedWrapped.porousTopology.periodicImage == 2
              && serialized(wrappedReplay.advance())
                  == serialized(expectedWrappedNext),
          "persisted moving porous-flow second wrap continues bit-for-bit");

    fsi::MovingPorousFlowCase preservedOwner;
    static_cast<void>(preservedOwner.advance());
    fsi::MovingPorousFlowCaseCheckpoint output =
        preservedOwner.checkpoint();
    const auto preservedEncoding = encodedCheckpoint(output);
    const auto expectDecodeRejected = [&output, &error, &preservedEncoding](
        const std::vector<std::uint8_t>& candidate,
        const fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode code,
        const fsi::MovingPorousFlowCaseCheckpointPersistenceLimits& limits,
        const char* message) {
        error = {};
        check(!fsi::deserializeMovingPorousFlowCaseCheckpoint(
                  candidate, output, &error, limits)
                  && error.code == code
                  && encodedCheckpoint(output) == preservedEncoding,
              message);
    };

    auto corrupt = firstEncoding;
    corrupt[0] ^= 0xffU;
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidMagic,
        {}, "moving porous-flow checkpoint rejects invalid magic");
    corrupt = firstEncoding;
    ++corrupt[4];
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::UnsupportedVersion,
        {}, "moving porous-flow checkpoint rejects wire versions");
    corrupt = firstEncoding;
    corrupt[6] = 1;
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects reserved bits");
    corrupt = firstEncoding;
    corrupt[180] ^= 0x01U;
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::ChecksumMismatch,
        {}, "moving porous-flow checkpoint rejects checksum corruption");
    corrupt = firstEncoding;
    corrupt.pop_back();
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::Truncated,
        {}, "moving porous-flow checkpoint rejects truncation");
    corrupt = firstEncoding;
    corrupt.push_back(0);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::TrailingData,
        {}, "moving porous-flow checkpoint rejects trailing data");
    corrupt = firstEncoding;
    ++corrupt[16];
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::UnsupportedVersion,
        {}, "moving porous-flow checkpoint rejects state versions");
    corrupt = firstEncoding;
    corrupt[20] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects foreign definitions");
    corrupt = firstEncoding;
    corrupt[28] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects foreign grids");
    corrupt = firstEncoding;
    corrupt[116] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects inconsistent steps");
    corrupt = firstEncoding;
    corrupt[160] = 1;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects inconsistent topology");
    corrupt = firstEncoding;
    corrupt[180] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects inconsistent fluid fields");
    constexpr std::size_t firstJumpReservedOffset = 976;
    corrupt = firstEncoding;
    corrupt[firstJumpReservedOffset] = 1;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects crossing reserved bits");
    constexpr std::size_t firstJumpPressureOffset = 1004;
    corrupt = firstEncoding;
    corrupt[firstJumpPressureOffset] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
        {}, "moving porous-flow checkpoint rejects inconsistent crossings");

    auto limited = fsi::MovingPorousFlowCaseCheckpointPersistenceLimits{};
    limited.maximumEncodedBytes = firstEncoding.size() - 1;
    expectDecodeRejected(
        firstEncoding,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
        limited, "moving porous-flow checkpoint enforces byte limits");
    limited = {};
    limited.maximumScalarSamples = 95;
    expectDecodeRejected(
        firstEncoding,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
        limited, "moving porous-flow checkpoint enforces sample limits");
    limited = {};
    limited.maximumPressureJumpFaces = 11;
    expectDecodeRejected(
        firstEncoding,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
        limited, "moving porous-flow checkpoint enforces crossing limits");
    limited = {};
    limited.maximumReplaySteps = savedSteps - 1;
    expectDecodeRejected(
        firstEncoding,
        fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
        limited, "moving porous-flow checkpoint bounds deterministic replay");

    std::vector<std::uint8_t> preservedBytes{1, 2, 3};
    auto invalidCheckpoint = saved;
    ++invalidCheckpoint.version;
    error = {};
    check(!fsi::serializeMovingPorousFlowCaseCheckpoint(
              invalidCheckpoint, preservedBytes, &error)
              && error.code
                  == fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData
              && preservedBytes == std::vector<std::uint8_t>({1, 2, 3}),
          "moving porous-flow checkpoint failed encoding preserves output");
    preservedBytes = {4, 5, 6};
    limited = {};
    limited.maximumEncodedBytes = firstEncoding.size() - 1;
    error = {};
    check(!fsi::serializeMovingPorousFlowCaseCheckpoint(
              saved, preservedBytes, &error, limited)
              && error.code
                  == fsi::MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded
              && preservedBytes == std::vector<std::uint8_t>({4, 5, 6}),
          "moving porous-flow checkpoint encoding enforces bytes transactionally");
}

} // namespace

int main() {
    testMovingPorousWorker();
    testCheckpointReplayAndValidation();
    testPersistentCheckpoint();
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
