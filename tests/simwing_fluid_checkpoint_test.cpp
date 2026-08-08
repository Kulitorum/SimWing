#include "fluid/checkpoint.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpoint;
using simwing::fsi::fluid::MovingPorousFluidCheckpoint;
using simwing::fsi::fluid::MovingPorousProjectionSettings;
using simwing::fsi::fluid::PorousConstitutiveEvaluation;
using simwing::fsi::fluid::PorousGridFaceCrossing;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointError;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointErrorCode;
using simwing::fsi::fluid::MovingInterfaceFluidCheckpointLimits;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::SharpPressureJumpField;
using simwing::fsi::fluid::checkpointMovingInterfaceFluidState;
using simwing::fsi::fluid::checkpointMovingPorousFluidState;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;
using simwing::fsi::fluid::projectVelocityWithMovingAndPorousInterfaces;
using simwing::fsi::fluid::restoreMovingInterfaceFluidState;
using simwing::fsi::fluid::restoreMovingPorousFluidState;

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

PeriodicCartesianGrid movingPorousGrid() {
    return PeriodicCartesianGrid({8, 2, 3}, {}, {4.0, 2.0, 3.0});
}

FaceAlignedMovingInterface movingPorousInterfaces(
    const PeriodicCartesianGrid& grid) {
    std::vector<GridFaceMovingInterface> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                100, 1, 2, GridFaceAxis::X, 2, j, k, 0.25});
            faces.push_back({
                200, 2, 1, GridFaceAxis::X, 6, j, k, 0.25});
        }
    }
    return FaceAlignedMovingInterface(grid, std::move(faces));
}

std::vector<PorousGridFaceCrossing> movingPorousCrossings(
    const PeriodicCartesianGrid& grid,
    const std::size_t plane = 3,
    const double surfaceVelocityMetersPerSecond = 0.1,
    const double linearResistance = 10.0) {
    std::vector<PorousGridFaceCrossing> crossings;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            crossings.push_back({
                300, 10, 11, GridFaceAxis::X, plane, j, k,
                0.4, surfaceVelocityMetersPerSecond,
                {linearResistance, 0.0}});
        }
    }
    return crossings;
}

SharpPressureJumpField movingPorousBalance(
    const PeriodicCartesianGrid& grid,
    const double pressureJumpPascals) {
    std::vector<simwing::fsi::fluid::GridFacePressureJump> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                400, 11, 10, GridFaceAxis::X, 5, j, k,
                pressureJumpPascals, 0.6});
        }
    }
    return SharpPressureJumpField(grid, std::move(faces));
}

MovingPorousProjectionSettings movingPorousSettings(
    const PorousConstitutiveEvaluation evaluation) {
    MovingPorousProjectionSettings settings;
    settings.movingProjection = projectionSettings();
    settings.movingProjection.projection.timeStepSeconds = 0.4;
    settings.movingProjection
        .absoluteRegionVolumeRateToleranceCubicMetersPerSecond = 1.0e-12;
    settings.iteration.constitutiveEvaluation = evaluation;
    settings.iteration.absoluteNormalVelocityToleranceMetersPerSecond =
        1.0e-12;
    settings.iteration.relativeNormalVelocityTolerance = 1.0e-12;
    settings.iteration.absolutePressureJumpTolerancePascals = 1.0e-11;
    settings.iteration.relativePressureJumpTolerance = 1.0e-12;
    settings.iteration.relaxation = 0.5;
    settings.iteration.maximumNonlinearIterations = 100;
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

void testMovingPorousAcceptedStateRoundTrip() {
    const auto grid = movingPorousGrid();
    const auto interfaces = movingPorousInterfaces(grid);
    const auto porous = movingPorousCrossings(grid);
    for (const auto evaluation : {
             PorousConstitutiveEvaluation::Endpoint,
             PorousConstitutiveEvaluation::Midpoint}) {
        const auto prescribed = movingPorousBalance(
            grid, evaluation == PorousConstitutiveEvaluation::Endpoint
                ? 1.5 : 0.25);
        const auto settings = movingPorousSettings(evaluation);
        MacVelocityField predicted(grid);
        MacVelocityField velocity = predicted;
        CellScalarField pressure(grid);
        const auto diagnostics =
            projectVelocityWithMovingAndPorousInterfaces(
                grid, velocity, pressure, interfaces,
                porous, prescribed, settings);
        const auto checkpoint = checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            porous, prescribed, diagnostics);
        const auto restored = restoreMovingPorousFluidState(
            grid, checkpoint);
        const auto replayed = restoreMovingPorousFluidState(
            grid, checkpoint);

        check(checkpoint.version
                  == simwing::fsi::fluid::movingPorousFluidCheckpointVersion
                  && checkpoint.cellCounts == grid.cellCounts()
                  && checkpoint.lowerMeters == grid.lowerMeters()
                  && checkpoint.upperMeters == grid.upperMeters()
                  && checkpoint.scalarSampleCount == grid.cellCount()
                  && checkpoint.topologyFingerprint != 0,
              "moving porous checkpoint: public epoch metadata are explicit");
        check(restored.predictedVelocityMetersPerSecond == predicted
                  && restored.velocityMetersPerSecond == velocity
                  && restored.pressurePascals == pressure
                  && restored.interfaces == interfaces
                  && restored.porousCrossings == porous
                  && restored.prescribedPressureJumps == prescribed
                  && restored.diagnostics == diagnostics,
              "moving porous checkpoint: accepted endpoint/midpoint state restores bit-for-bit");
        check(replayed.predictedVelocityMetersPerSecond
                      == restored.predictedVelocityMetersPerSecond
                  && replayed.velocityMetersPerSecond
                      == restored.velocityMetersPerSecond
                  && replayed.pressurePascals == restored.pressurePascals
                  && replayed.porousCrossings == restored.porousCrossings
                  && replayed.diagnostics == restored.diagnostics,
              "moving porous checkpoint: repeated restore replays bit-for-bit");

        predicted.xFaces().front() += 10.0;
        velocity.xFaces().front() += 20.0;
        pressure.values().front() += 30.0;
        const auto immutable = restoreMovingPorousFluidState(
            grid, checkpoint);
        check(immutable.predictedVelocityMetersPerSecond
                      == restored.predictedVelocityMetersPerSecond
                  && immutable.velocityMetersPerSecond
                      == restored.velocityMetersPerSecond
                  && immutable.pressurePascals == restored.pressurePascals,
              "moving porous checkpoint: immutable payload does not alias source fields");
    }
}

void testMovingPorousTopologyBinding() {
    const auto grid = movingPorousGrid();
    const auto interfaces = movingPorousInterfaces(grid);
    const auto settings = movingPorousSettings(
        PorousConstitutiveEvaluation::Endpoint);
    const auto capture = [&](const std::vector<PorousGridFaceCrossing>& porous,
                             const SharpPressureJumpField& prescribed) {
        MacVelocityField predicted(grid);
        MacVelocityField velocity = predicted;
        CellScalarField pressure(grid);
        const auto diagnostics =
            projectVelocityWithMovingAndPorousInterfaces(
                grid, velocity, pressure, interfaces,
                porous, prescribed, settings);
        return checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            porous, prescribed, diagnostics);
    };

    const auto porous = movingPorousCrossings(grid);
    const auto baseline = capture(
        porous, movingPorousBalance(grid, 1.5));
    auto reversed = porous;
    std::reverse(reversed.begin(), reversed.end());
    const auto reversedCheckpoint = capture(
        reversed, movingPorousBalance(grid, 1.5));
    check(baseline.topologyFingerprint
                  == reversedCheckpoint.topologyFingerprint
              && restoreMovingPorousFluidState(grid, baseline)
                      .porousCrossings
                  == restoreMovingPorousFluidState(
                      grid, reversedCheckpoint).porousCrossings,
          "moving porous checkpoint: authored order canonicalizes to one epoch");

    const auto moved = capture(
        movingPorousCrossings(grid, 4),
        movingPorousBalance(grid, 1.5));
    check(moved.topologyFingerprint != baseline.topologyFingerprint,
          "moving porous checkpoint: a changed porous face starts a new topology epoch");

    const auto changedState = capture(
        movingPorousCrossings(grid, 3, 0.05, 20.0),
        movingPorousBalance(grid, 4.0));
    check(changedState.topologyFingerprint == baseline.topologyFingerprint,
          "moving porous checkpoint: material, velocity, and jump values do not change topology identity");
}

void testMovingPorousStrictValidation() {
    const auto grid = movingPorousGrid();
    const auto interfaces = movingPorousInterfaces(grid);
    const auto porous = movingPorousCrossings(grid);
    const auto prescribed = movingPorousBalance(grid, 0.25);
    const auto settings = movingPorousSettings(
        PorousConstitutiveEvaluation::Midpoint);
    MacVelocityField predicted(grid);
    MacVelocityField velocity = predicted;
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithMovingAndPorousInterfaces(
        grid, velocity, pressure, interfaces,
        porous, prescribed, settings);
    const auto checkpoint = checkpointMovingPorousFluidState(
        grid, predicted, velocity, pressure, interfaces,
        porous, prescribed, diagnostics);

    auto wrongVersion = checkpoint;
    ++wrongVersion.version;
    expectRejected(
        [&] { static_cast<void>(restoreMovingPorousFluidState(
            grid, wrongVersion)); },
        "moving porous checkpoint validation: unsupported version is rejected");
    auto wrongFingerprint = checkpoint;
    ++wrongFingerprint.topologyFingerprint;
    expectRejected(
        [&] { static_cast<void>(restoreMovingPorousFluidState(
            grid, wrongFingerprint)); },
        "moving porous checkpoint validation: corrupted topology fingerprint is rejected");

    auto failed = diagnostics;
    failed.accepted = false;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            porous, prescribed, failed)); },
        "moving porous checkpoint validation: rejected outer solve cannot be captured");
    auto mismatchedProjection = diagnostics;
    ++mismatchedProjection.porous.projection.pressureJumpFaceCount;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            porous, prescribed, mismatchedProjection)); },
        "moving porous checkpoint validation: nested projections must match exactly");
    auto corruptedSample = diagnostics;
    corruptedSample.porous.samples.front()
        .relativeNormalVelocityMetersPerSecond += 0.01;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            porous, prescribed, corruptedSample)); },
        "moving porous checkpoint validation: samples must match the calibrated law");

    auto wrongProvenance = predicted;
    wrongProvenance.xFaces()[grid.cellIndex(5, 0, 0)] += 0.1;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, wrongProvenance, velocity, pressure, interfaces,
            porous, prescribed, diagnostics)); },
        "moving porous checkpoint validation: midpoint source power binds predicted velocity");
    auto changedResistance = porous;
    changedResistance.front().resistance.linearPascalSecondsPerMeter += 1.0;
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, predicted, velocity, pressure, interfaces,
            changedResistance, prescribed, diagnostics)); },
        "moving porous checkpoint validation: material definitions bind diagnostic samples");

    auto nonfinitePredicted = predicted;
    nonfinitePredicted.zFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(checkpointMovingPorousFluidState(
            grid, nonfinitePredicted, velocity, pressure, interfaces,
            porous, prescribed, diagnostics)); },
        "moving porous checkpoint validation: non-finite provenance cannot be captured");
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
    testMovingPorousAcceptedStateRoundTrip();
    testMovingPorousTopologyBinding();
    testMovingPorousStrictValidation();
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
