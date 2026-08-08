#include "fluid/checkpoint.h"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpoint;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointError;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointErrorCode;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointLimits;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::checkpointMovingInterfaceFluidState;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;
using simwing::fsi::fluid::restoreMovingInterfaceFluidState;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

std::vector<std::uint8_t> serializedCheckpoint(
    const MovingInterfaceFluidCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    MovingInterfaceFluidCheckpointError error;
    check(simwing::fsi::fluid::serializeMovingInterfaceFluidCheckpoint(
              checkpoint, bytes, &error)
              && !error && !bytes.empty(),
          "persistent checkpoint helper serializes accepted state");
    return bytes;
}

void refreshPersistentChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 24;
    constexpr std::size_t checksumOffset = 16;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t checksum = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        checksum ^= bytes[index];
        checksum *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(checksum); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            checksum >> (8U * byte));
    }
}

bool sameCheckpoint(
    const PeriodicCartesianGrid& grid,
    const MovingInterfaceFluidCheckpoint& first,
    const MovingInterfaceFluidCheckpoint& second) {
    const auto firstState = restoreMovingInterfaceFluidState(grid, first);
    const auto secondState = restoreMovingInterfaceFluidState(grid, second);
    return first.version == second.version
        && first.cellCounts == second.cellCounts
        && first.lowerMeters == second.lowerMeters
        && first.upperMeters == second.upperMeters
        && first.scalarSampleCount == second.scalarSampleCount
        && first.topologyFingerprint == second.topologyFingerprint
        && firstState.velocityMetersPerSecond
            == secondState.velocityMetersPerSecond
        && firstState.pressurePascals == secondState.pressurePascals
        && firstState.interfaces == secondState.interfaces
        && firstState.diagnostics == secondState.diagnostics;
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

PeriodicCartesianGrid makeGrid() {
    return PeriodicCartesianGrid({4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

FaceAlignedMovingInterface makeInterfaces(
    const PeriodicCartesianGrid& grid,
    const std::size_t plane,
    const double speedMetersPerSecond) {
    std::vector<GridFaceMovingInterface> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                300, 9, 9, GridFaceAxis::X,
                plane, j, k, speedMetersPerSecond,
            });
        }
    }
    return FaceAlignedMovingInterface(grid, std::move(faces));
}

MovingInterfaceProjectionSettings projectionSettings() {
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.1;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return settings;
}

void testAcceptedStateRoundTrip() {
    const auto grid = makeGrid();
    const auto interfaces = makeInterfaces(grid, 3, 0.125);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());
    const auto checkpoint = checkpointMovingInterfaceFluidState(
        grid, velocity, pressure, interfaces, diagnostics);
    const auto restored = restoreMovingInterfaceFluidState(grid, checkpoint);
    const auto replayed = restoreMovingInterfaceFluidState(grid, checkpoint);

    check(checkpoint.version
              == simwing::fsi::fluid::movingInterfaceFluidCheckpointVersion
              && checkpoint.cellCounts == grid.cellCounts()
              && checkpoint.lowerMeters == grid.lowerMeters()
              && checkpoint.upperMeters == grid.upperMeters()
              && checkpoint.scalarSampleCount == grid.cellCount()
              && checkpoint.topologyFingerprint != 0,
          "checkpoint: public grid and topology metadata are explicit");
    check(restored.velocityMetersPerSecond == velocity
              && restored.pressurePascals == pressure
              && restored.interfaces == interfaces
              && restored.diagnostics == diagnostics,
          "checkpoint: accepted fluid state restores bit-for-bit");
    check(replayed.velocityMetersPerSecond
              == restored.velocityMetersPerSecond
              && replayed.pressurePascals == restored.pressurePascals
              && replayed.interfaces == restored.interfaces
              && replayed.diagnostics == restored.diagnostics,
          "checkpoint: repeated restore replays bit-for-bit");

    velocity.xFaces().front() += 10.0;
    pressure.values().front() += 20.0;
    const auto afterSourceMutation = restoreMovingInterfaceFluidState(
        grid, checkpoint);
    check(afterSourceMutation.velocityMetersPerSecond
                  == restored.velocityMetersPerSecond
              && afterSourceMutation.pressurePascals
                  == restored.pressurePascals,
          "checkpoint: immutable payload does not alias mutable source fields");
}

void testTopologyEpochBinding() {
    const auto grid = makeGrid();
    const auto firstInterfaces = makeInterfaces(grid, 3, 0.125);
    MacVelocityField firstVelocity(grid);
    CellScalarField firstPressure(grid);
    const auto firstDiagnostics = projectVelocityWithMovingInterfaces(
        grid, firstVelocity, firstPressure,
        firstInterfaces, projectionSettings());
    const auto first = checkpointMovingInterfaceFluidState(
        grid, firstVelocity, firstPressure,
        firstInterfaces, firstDiagnostics);

    const auto rebasedInterfaces = makeInterfaces(grid, 0, 0.125);
    MacVelocityField rebasedVelocity = firstVelocity;
    CellScalarField rebasedPressure = firstPressure;
    const auto rebasedDiagnostics = projectVelocityWithMovingInterfaces(
        grid, rebasedVelocity, rebasedPressure,
        rebasedInterfaces, projectionSettings());
    const auto rebased = checkpointMovingInterfaceFluidState(
        grid, rebasedVelocity, rebasedPressure,
        rebasedInterfaces, rebasedDiagnostics);

    check(first.topologyFingerprint != rebased.topologyFingerprint,
          "checkpoint: a periodic face rebase starts a distinct topology epoch");
    const auto restoredFirst = restoreMovingInterfaceFluidState(grid, first);
    check(restoredFirst.interfaces == firstInterfaces
              && restoredFirst.velocityMetersPerSecond == firstVelocity
              && restoredFirst.pressurePascals == firstPressure,
          "checkpoint: an earlier topology epoch restores after a later solve");

    const auto sameTopologyDifferentSpeed = makeInterfaces(grid, 3, 0.25);
    MacVelocityField fasterVelocity(grid);
    CellScalarField fasterPressure(grid);
    const auto fasterDiagnostics = projectVelocityWithMovingInterfaces(
        grid, fasterVelocity, fasterPressure,
        sameTopologyDifferentSpeed, projectionSettings());
    const auto faster = checkpointMovingInterfaceFluidState(
        grid, fasterVelocity, fasterPressure,
        sameTopologyDifferentSpeed, fasterDiagnostics);
    check(faster.topologyFingerprint == first.topologyFingerprint,
          "checkpoint: prescribed speed changes state without changing topology identity");
}

void testPersistentRoundTrip() {
    const auto grid = makeGrid();
    const auto interfaces = makeInterfaces(grid, 3, 0.125);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());
    const auto checkpoint = checkpointMovingInterfaceFluidState(
        grid, velocity, pressure, interfaces, diagnostics);

    const auto first = serializedCheckpoint(checkpoint);
    const auto second = serializedCheckpoint(checkpoint);
    MovingInterfaceFluidCheckpoint decoded;
    MovingInterfaceFluidCheckpointError error{
        MovingInterfaceFluidCheckpointErrorCode::InvalidData, "old"};
    check(first == second
              && simwing::fsi::fluid::deserializeMovingInterfaceFluidCheckpoint(
                  first, decoded, &error)
              && !error
              && sameCheckpoint(grid, checkpoint, decoded)
              && serializedCheckpoint(decoded) == first,
          "persistent checkpoint round trip preserves every field and byte");

    const auto rebasedInterfaces = makeInterfaces(grid, 0, 0.125);
    MacVelocityField rebasedVelocity = velocity;
    CellScalarField rebasedPressure = pressure;
    const auto rebasedDiagnostics = projectVelocityWithMovingInterfaces(
        grid, rebasedVelocity, rebasedPressure,
        rebasedInterfaces, projectionSettings());
    const auto rebased = checkpointMovingInterfaceFluidState(
        grid, rebasedVelocity, rebasedPressure,
        rebasedInterfaces, rebasedDiagnostics);
    MovingInterfaceFluidCheckpoint decodedRebased;
    const auto rebasedBytes = serializedCheckpoint(rebased);
    check(simwing::fsi::fluid::deserializeMovingInterfaceFluidCheckpoint(
              rebasedBytes, decodedRebased, &error)
              && sameCheckpoint(grid, rebased, decodedRebased)
              && decodedRebased.topologyFingerprint
                  != decoded.topologyFingerprint,
          "persistent checkpoint retains a rebased topology epoch exactly");
}

void testPersistentCorruptionAndLimits() {
    const auto grid = makeGrid();
    const auto interfaces = makeInterfaces(grid, 3, 0.125);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());
    const auto checkpoint = checkpointMovingInterfaceFluidState(
        grid, velocity, pressure, interfaces, diagnostics);
    const auto valid = serializedCheckpoint(checkpoint);
    MovingInterfaceFluidCheckpoint output = checkpoint;
    MovingInterfaceFluidCheckpointError error;
    const auto expectRejected = [&]
        (const std::vector<std::uint8_t>& candidate,
         const MovingInterfaceFluidCheckpointErrorCode expected,
         const MovingInterfaceFluidCheckpointLimits& limits,
         const char* message) {
        output = checkpoint;
        error = {};
        check(!simwing::fsi::fluid::deserializeMovingInterfaceFluidCheckpoint(
                  candidate, output, &error, limits)
                  && error.code == expected
                  && sameCheckpoint(grid, checkpoint, output),
              message);
    };

    auto corrupt = valid;
    corrupt[0] ^= 0xffU;
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::InvalidMagic, {},
        "persistent checkpoint rejects bad magic transactionally");
    corrupt = valid;
    ++corrupt[4];
    expectRejected(
        corrupt,
        MovingInterfaceFluidCheckpointErrorCode::UnsupportedVersion, {},
        "persistent checkpoint rejects unsupported envelope versions");
    corrupt = valid;
    corrupt[6] = 1;
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::InvalidData, {},
        "persistent checkpoint rejects nonzero reserved bits");
    corrupt = valid;
    corrupt.pop_back();
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::Truncated, {},
        "persistent checkpoint rejects truncated payloads");
    corrupt = valid;
    corrupt.push_back(0);
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::TrailingData, {},
        "persistent checkpoint rejects trailing bytes");
    corrupt = valid;
    corrupt.back() ^= 0x80U;
    expectRejected(
        corrupt,
        MovingInterfaceFluidCheckpointErrorCode::ChecksumMismatch, {},
        "persistent checkpoint detects payload corruption");

    MovingInterfaceFluidCheckpointLimits byteLimit;
    byteLimit.maximumBytes = valid.size() - 1;
    expectRejected(
        valid, MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
        byteLimit, "persistent checkpoint enforces its byte limit");
    MovingInterfaceFluidCheckpointLimits sampleLimit;
    sampleLimit.maximumScalarSamples = grid.cellCount() - 1;
    expectRejected(
        valid, MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
        sampleLimit, "persistent checkpoint enforces its sample limit");
    MovingInterfaceFluidCheckpointLimits faceLimit;
    faceLimit.maximumInterfaceFaces = interfaces.faceCount() - 1;
    expectRejected(
        valid, MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
        faceLimit, "persistent checkpoint enforces its interface-face limit");
    std::vector<std::uint8_t> limitedBytes{1, 2, 3};
    check(!simwing::fsi::fluid::serializeMovingInterfaceFluidCheckpoint(
              checkpoint, limitedBytes, &error, faceLimit)
              && error.code
                  == MovingInterfaceFluidCheckpointErrorCode::LimitExceeded
              && limitedBytes.empty(),
          "persistent encoding enforces topology limits transactionally");
    MovingInterfaceFluidCheckpointLimits regionLimit;
    regionLimit.maximumFluidRegions = 0;
    expectRejected(
        valid, MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
        regionLimit, "persistent checkpoint enforces its region limit");
    MovingInterfaceFluidCheckpointLimits surfaceLimit;
    surfaceLimit.maximumDiagnosticSurfaces = 0;
    expectRejected(
        valid, MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
        surfaceLimit,
        "persistent checkpoint enforces its diagnostic-surface limit");

    constexpr std::size_t envelopeBytes = 24;
    constexpr std::size_t topologyFingerprintPayloadOffset = 84;
    constexpr std::size_t interfaceFaceCountPayloadOffset = 96;
    corrupt = valid;
    for (std::size_t byte = 0; byte < sizeof(std::uint64_t); ++byte) {
        corrupt[envelopeBytes + interfaceFaceCountPayloadOffset + byte] =
            static_cast<std::uint8_t>(1'000'000ULL >> (8U * byte));
    }
    refreshPersistentChecksum(corrupt);
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::Truncated, {},
        "persistent checkpoint checks record bytes before allocation");

    corrupt = valid;
    corrupt[envelopeBytes + topologyFingerprintPayloadOffset] ^= 1U;
    refreshPersistentChecksum(corrupt);
    expectRejected(
        corrupt, MovingInterfaceFluidCheckpointErrorCode::InvalidData, {},
        "persistent checkpoint recomputes its topology fingerprint");

    auto invalid = checkpoint;
    ++invalid.topologyFingerprint;
    std::vector<std::uint8_t> bytes{1, 2, 3};
    check(!simwing::fsi::fluid::serializeMovingInterfaceFluidCheckpoint(
              invalid, bytes, &error)
              && error.code
                  == MovingInterfaceFluidCheckpointErrorCode::InvalidData
              && bytes.empty(),
          "persistent serialization rejects invalid source metadata");
}

void testStrictValidation() {
    const auto grid = makeGrid();
    const auto interfaces = makeInterfaces(grid, 3, 0.125);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());
    const auto checkpoint = checkpointMovingInterfaceFluidState(
        grid, velocity, pressure, interfaces, diagnostics);

    auto wrongVersion = checkpoint;
    ++wrongVersion.version;
    expectRejected(
        [&] { static_cast<void>(restoreMovingInterfaceFluidState(
            grid, wrongVersion)); },
        "checkpoint validation: unsupported version is rejected");
    auto wrongFingerprint = checkpoint;
    ++wrongFingerprint.topologyFingerprint;
    expectRejected(
        [&] { static_cast<void>(restoreMovingInterfaceFluidState(
            grid, wrongFingerprint)); },
        "checkpoint validation: corrupted topology fingerprint is rejected");
    auto wrongGridMetadata = checkpoint;
    wrongGridMetadata.upperMeters.x += 1.0;
    expectRejected(
        [&] { static_cast<void>(restoreMovingInterfaceFluidState(
            grid, wrongGridMetadata)); },
        "checkpoint validation: corrupted public grid binding is rejected");
    auto wrongSampleCount = checkpoint;
    ++wrongSampleCount.scalarSampleCount;
    expectRejected(
        [&] { static_cast<void>(restoreMovingInterfaceFluidState(
            grid, wrongSampleCount)); },
        "checkpoint validation: corrupted field shape is rejected");

    const PeriodicCartesianGrid foreignGrid(
        {4, 4, 4}, {1.0, 0.0, 0.0}, {5.0, 4.0, 4.0});
    expectRejected(
        [&] { static_cast<void>(restoreMovingInterfaceFluidState(
            foreignGrid, checkpoint)); },
        "checkpoint validation: foreign grid geometry is rejected");

    auto failedDiagnostics = diagnostics;
    failedDiagnostics.projection.converged = false;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingInterfaceFluidState(
            grid, velocity, pressure, interfaces, failedDiagnostics)); },
        "checkpoint validation: failed projection cannot be captured");
    auto mismatchedDiagnostics = diagnostics;
    mismatchedDiagnostics.faces.front()
        .normalVelocityMetersPerSecond += 0.01;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingInterfaceFluidState(
            grid, velocity, pressure,
            interfaces, mismatchedDiagnostics)); },
        "checkpoint validation: diagnostics must bind exact interface kinematics");
    auto corruptedAggregate = diagnostics;
    corruptedAggregate.surfaces.front()
        .constraintReactionForceNewtons.x += 1.0;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingInterfaceFluidState(
            grid, velocity, pressure,
            interfaces, corruptedAggregate)); },
        "checkpoint validation: face and surface reaction ledgers must agree");
    auto incompatibleVelocity = velocity;
    incompatibleVelocity.xFaces()[grid.cellIndex(3, 0, 0)] += 0.01;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingInterfaceFluidState(
            grid, incompatibleVelocity, pressure,
            interfaces, diagnostics)); },
        "checkpoint validation: constrained MAC velocity must remain accepted");
    auto nonfiniteVelocity = velocity;
    nonfiniteVelocity.zFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(checkpointMovingInterfaceFluidState(
            grid, nonfiniteVelocity, pressure,
            interfaces, diagnostics)); },
        "checkpoint validation: non-finite fields cannot be captured");
}

} // namespace

int main() {
    testAcceptedStateRoundTrip();
    testTopologyEpochBinding();
    testPersistentRoundTrip();
    testPersistentCorruptionAndLimits();
    testStrictValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid checkpoint check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid checkpoint checks passed");
    return 0;
}
