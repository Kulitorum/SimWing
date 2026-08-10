#include "scene_fluid_regional_opening_momentum_wall_structure_step_epoch.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    void bytes(const std::span<const std::uint8_t> values) {
        for (const std::uint8_t value : values) {
            integer(value);
        }
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void fingerprintVector(Fingerprint& fingerprint,
                       const StructureVector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening wall structural-step storage overflows");
    }
    return first + second;
}

std::uint64_t bytesFingerprint(
    const std::span<const std::uint8_t> bytes) {
    Fingerprint fingerprint;
    fingerprint.bytes(bytes);
    return fingerprint.value();
}

void fingerprintStepSettings(Fingerprint& fingerprint,
                             const StructureStepSettings& settings) {
    fingerprint.real(settings.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint32_t>(settings.substeps));
    fingerprint.integer(
        static_cast<std::uint32_t>(settings.constraintIterations));
    fingerprint.integer(
        static_cast<std::uint32_t>(settings.cableConstraintSweepPairs));
    fingerprintVector(
        fingerprint, settings.gravityMetersPerSecondSquared);
    fingerprint.real(settings.velocityDampingPerSecond);
    fingerprintVector(
        fingerprint, settings.dampingReferenceVelocityMetersPerSecond);
    fingerprint.integer(settings.workerThreads);
}

void fingerprintDiagnostics(Fingerprint& fingerprint,
                            const StructureDiagnostics& diagnostics) {
    for (const std::size_t value : {
             diagnostics.nodeCount,
             diagnostics.dynamicNodeCount,
             diagnostics.triangleCount,
             diagnostics.constraintCount,
             diagnostics.membraneCount,
             diagnostics.dihedralCount,
             diagnostics.contactPairCount,
             diagnostics.activeContactCount,
             diagnostics.suspensionSegmentCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(diagnostics.totalDynamicMassKg);
    fingerprintVector(fingerprint, diagnostics.centerOfMassMeters);
    fingerprintVector(
        fingerprint, diagnostics.linearMomentumKgMetersPerSecond);
    fingerprint.real(diagnostics.kineticEnergyJoules);
    fingerprint.real(diagnostics.maximumDistanceErrorMeters);
    fingerprint.real(diagnostics.maximumCableExtensionMeters);
    fingerprint.real(diagnostics.maximumAbsoluteMembraneStrain);
    fingerprint.real(diagnostics.maximumMembraneResidual);
    fingerprint.real(diagnostics.maximumContactPenetrationMeters);
    fingerprint.real(diagnostics.maximumSuspensionResidualMeters);
    fingerprintVector(
        fingerprint, diagnostics.pendingExternalForceNewtons);
    fingerprintVector(
        fingerprint, diagnostics.lastAppliedExternalForceNewtons);
    fingerprint.boolean(diagnostics.finite);
}

std::uint64_t settingsFingerprint(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings) {
    Fingerprint fingerprint;
    for (const double value : {
             settings.load.pressureState
                 .absolutePressureResidualTolerancePascals,
             settings.load.pressureState.relativePressureResidualTolerance,
             settings.load.pressureState
                 .absoluteForceResidualToleranceNewtons,
             settings.load.pressureState.relativeForceResidualTolerance,
             settings.load.pressureState
                 .absoluteWorkResidualToleranceJoules,
             settings.load.pressureState.relativeWorkResidualTolerance,
             settings.load.transfer.momentReferenceMeters.x,
             settings.load.transfer.momentReferenceMeters.y,
             settings.load.transfer.momentReferenceMeters.z,
             settings.load.transfer.minimumTriangleAreaSquareMeters,
             settings.load.transfer.minimumQuadratureAreaSquareMeters,
             settings.load.transfer.barycentricTolerance,
             settings.load
                 .absoluteWallActionReactionToleranceKilogramMetersPerSecond,
             settings.load.relativeWallActionReactionTolerance}) {
        fingerprint.real(value);
    }
    fingerprintStepSettings(fingerprint, settings.structure);
    return fingerprint.value();
}

bool sameStepSettings(const StructureStepSettings& first,
                      const StructureStepSettings& second) {
    return first.timeStepSeconds == second.timeStepSeconds
        && first.substeps == second.substeps
        && first.constraintIterations == second.constraintIterations
        && first.cableConstraintSweepPairs
            == second.cableConstraintSweepPairs
        && first.gravityMetersPerSecondSquared
            == second.gravityMetersPerSecondSquared
        && first.velocityDampingPerSecond
            == second.velocityDampingPerSecond
        && first.dampingReferenceVelocityMetersPerSecond
            == second.dampingReferenceVelocityMetersPerSecond
        && first.workerThreads == second.workerThreads;
}

void validateStructureSettings(const StructureStepSettings& settings) {
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || settings.substeps <= 0
        || settings.constraintIterations < 0
        || settings.cableConstraintSweepPairs < 0
        || !finite(settings.gravityMetersPerSecondSquared)
        || !std::isfinite(settings.velocityDampingPerSecond)
        || settings.velocityDampingPerSecond < 0.0
        || !finite(settings.dampingReferenceVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "regional opening wall structural-step settings are invalid");
    }
}

void validateDiagnostics(const StructureDiagnostics& diagnostics) {
    if (!diagnostics.finite
        || !std::isfinite(diagnostics.totalDynamicMassKg)
        || !finite(diagnostics.centerOfMassMeters)
        || !finite(diagnostics.linearMomentumKgMetersPerSecond)
        || !std::isfinite(diagnostics.kineticEnergyJoules)
        || !std::isfinite(diagnostics.maximumDistanceErrorMeters)
        || !std::isfinite(diagnostics.maximumCableExtensionMeters)
        || !std::isfinite(diagnostics.maximumAbsoluteMembraneStrain)
        || !std::isfinite(diagnostics.maximumMembraneResidual)
        || !std::isfinite(diagnostics.maximumContactPenetrationMeters)
        || !std::isfinite(diagnostics.maximumSuspensionResidualMeters)
        || !finite(diagnostics.pendingExternalForceNewtons)
        || !finite(diagnostics.lastAppliedExternalForceNewtons)) {
        throw std::invalid_argument(
            "regional opening wall structural-step diagnostics are invalid");
    }
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0
        || limits.checkpointPersistence.maximumEncodedBytes == 0
        || limits.checkpointPersistence.maximumNodes == 0) {
        throw std::invalid_argument(
            "regional opening wall structural-step limits are invalid");
    }
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch) {
    return checkedAdd(
        epoch.loadEpoch.ownedStorageBytes,
        checkedAdd(epoch.beforeStructureCheckpoint.size(),
                   epoch.afterStructureCheckpoint.size()));
}

std::size_t workingStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState) {
    return checkedAdd(
        cycleState.ownedStorageBytes,
        checkedAdd(epoch.loadEpoch.workingStorageBytes,
                   epoch.ownedStorageBytes));
}

std::uint64_t epochFingerprint(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch) {
    Fingerprint fingerprint;
    fingerprint.integer(epoch.version);
    for (const std::uint64_t value : {
             epoch.sourceCycleStateFingerprint,
             epoch.sourceLoadEpochFingerprint,
             epoch.targetDefinitionFingerprint,
             epoch.sourceSettingsFingerprint,
             epoch.beforeCheckpointFingerprint,
             epoch.afterCheckpointFingerprint,
             epoch.beforeAcceptedStepCount,
             epoch.afterAcceptedStepCount,
             epoch.loadEpoch.fingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.real(epoch.beforeSimulationTimeSeconds);
    fingerprint.real(epoch.afterSimulationTimeSeconds);
    fingerprintStepSettings(fingerprint, epoch.structureSettings);
    fingerprintDiagnostics(fingerprint, epoch.diagnostics);
    fingerprint.boolean(epoch.stepped);
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.workingStorageBytes));
    return fingerprint.value();
}

void serializeCheckpointOrThrow(
    const Structure& owner,
    const StructureCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    const StructureCheckpointPersistenceLimits& limits,
    const char* context) {
    StructureCheckpointPersistenceError error;
    if (!serializeStructureCheckpoint(
            owner, checkpoint, bytes, &error, limits)) {
        const std::string message = std::string(context) + ": "
            + error.message;
        if (error.code
            == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
            throw std::length_error(message);
        }
        throw std::invalid_argument(message);
    }
}

StructureCheckpoint deserializeCheckpointOrThrow(
    const std::span<const std::uint8_t> bytes,
    const Structure& owner,
    const StructureCheckpointPersistenceLimits& limits,
    const char* context) {
    StructureCheckpoint result;
    StructureCheckpointPersistenceError error;
    if (!deserializeStructureCheckpoint(
            bytes, owner, result, &error, limits)) {
        const std::string message = std::string(context) + ": "
            + error.message;
        if (error.code
            == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
            throw std::length_error(message);
        }
        throw std::invalid_argument(message);
    }
    return result;
}

StructureVector3 totalPendingForce(const StructureCheckpoint& checkpoint) {
    StructureVector3 result;
    for (const StructureVector3& force :
         checkpoint.pendingExternalForcesNewtons) {
        result.x += force.x;
        result.y += force.y;
        result.z += force.z;
    }
    return result;
}

bool zeroPendingForces(const StructureCheckpoint& checkpoint) {
    return std::ranges::all_of(
        checkpoint.pendingExternalForcesNewtons,
        [](const StructureVector3& value) {
            return value == StructureVector3{};
        });
}

void validateCheckpointEndpoints(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch,
    const StructureCheckpoint& before,
    const StructureCheckpoint& after) {
    if (before.definitionFingerprint != epoch.targetDefinitionFingerprint
        || after.definitionFingerprint != epoch.targetDefinitionFingerprint
        || before.acceptedStepCount != epoch.beforeAcceptedStepCount
        || after.acceptedStepCount != epoch.afterAcceptedStepCount
        || before.simulationTimeSeconds
            != epoch.beforeSimulationTimeSeconds
        || after.simulationTimeSeconds != epoch.afterSimulationTimeSeconds
        || before.nodes.size() != after.nodes.size()
        || totalPendingForce(before)
            != epoch.loadEpoch.priorPendingForceNewtons
        || !zeroPendingForces(after)
        || after.lastAppliedExternalForceNewtons
            != epoch.loadEpoch.resultingPendingForceNewtons
        || epoch.diagnostics.pendingExternalForceNewtons
            != StructureVector3{}
        || epoch.diagnostics.lastAppliedExternalForceNewtons
            != epoch.loadEpoch.resultingPendingForceNewtons) {
        throw std::invalid_argument(
            "regional opening wall structural-step checkpoints are inconsistent");
    }
}

void validateAggregateLimits(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits) {
    if (epoch.ownedStorageBytes > limits.maximumOwnedBytes
        || epoch.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening wall structural-step aggregate limit exceeded");
    }
}

} // namespace

SceneFluidRegionalOpeningMomentumWallStructureStepEpoch
advanceSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits) {
    validateLimits(limits);
    validateStructureSettings(settings.structure);
    if (settings.structure.timeStepSeconds
        != cycleState.adjustedMomentum.timeStepSeconds) {
        throw std::invalid_argument(
            "regional opening wall structural and fluid time steps differ");
    }

    const StructureCheckpoint before = target.checkpoint();
    try {
        SceneFluidRegionalOpeningMomentumWallStructureStepEpoch result;
        result.sourceCycleStateFingerprint = cycleState.fingerprint;
        result.targetDefinitionFingerprint = target.definitionFingerprint();
        result.sourceSettingsFingerprint = settingsFingerprint(settings);
        result.beforeAcceptedStepCount = before.acceptedStepCount;
        result.beforeSimulationTimeSeconds = before.simulationTimeSeconds;
        serializeCheckpointOrThrow(
            target, before, result.beforeStructureCheckpoint,
            limits.checkpointPersistence,
            "cannot retain pre-step Structure checkpoint");
        result.beforeCheckpointFingerprint =
            bytesFingerprint(result.beforeStructureCheckpoint);

        result.loadEpoch =
            applySceneFluidRegionalOpeningMomentumWallLoadEpoch(
                cycleState, transportMetric, acceptedPressureOperator,
                acceptedBasePressureOperator, grid, acceptedSweep,
                acceptedFragments, acceptedTopology, acceptedVolumeRates,
                acceptedOpeningDefinitions, acceptedOpenings,
                acceptedResistanceDefinitions, acceptedBaseMetric,
                acceptedMetric, surface, surfaceState, transfer,
                quadrature, target, settings.load, limits.load);
        result.sourceLoadEpochFingerprint = result.loadEpoch.fingerprint;
        result.structureSettings = settings.structure;
        result.diagnostics = target.step(settings.structure);
        validateDiagnostics(result.diagnostics);

        const StructureCheckpoint after = target.checkpoint();
        result.afterAcceptedStepCount = after.acceptedStepCount;
        result.afterSimulationTimeSeconds = after.simulationTimeSeconds;
        serializeCheckpointOrThrow(
            target, after, result.afterStructureCheckpoint,
            limits.checkpointPersistence,
            "cannot retain post-step Structure checkpoint");
        result.afterCheckpointFingerprint =
            bytesFingerprint(result.afterStructureCheckpoint);
        result.stepped = true;
        result.ownedStorageBytes = ownedStorageBytes(result);
        result.workingStorageBytes = workingStorageBytes(result, cycleState);
        validateCheckpointEndpoints(result, before, after);
        validateAggregateLimits(result, limits);
        result.fingerprint = epochFingerprint(result);
        validateSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
            result, cycleState, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedBaseMetric, acceptedMetric, surface, surfaceState,
            transfer, quadrature, target, settings, limits);
        return result;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch) {
    validateSceneFluidRegionalOpeningMomentumWallLoadEpochIntegrity(
        epoch.loadEpoch);
    validateStructureSettings(epoch.structureSettings);
    validateDiagnostics(epoch.diagnostics);
    if (epoch.version
            != sceneFluidRegionalOpeningMomentumWallStructureStepEpochVersion
        || epoch.fingerprint == 0
        || epoch.sourceCycleStateFingerprint == 0
        || epoch.sourceLoadEpochFingerprint == 0
        || epoch.targetDefinitionFingerprint == 0
        || epoch.sourceSettingsFingerprint == 0
        || epoch.beforeCheckpointFingerprint == 0
        || epoch.afterCheckpointFingerprint == 0
        || epoch.beforeStructureCheckpoint.empty()
        || epoch.afterStructureCheckpoint.empty()
        || !std::isfinite(epoch.beforeSimulationTimeSeconds)
        || !std::isfinite(epoch.afterSimulationTimeSeconds)
        || epoch.beforeSimulationTimeSeconds < 0.0
        || epoch.afterSimulationTimeSeconds < 0.0
        || epoch.beforeAcceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || epoch.afterAcceptedStepCount
            != epoch.beforeAcceptedStepCount + 1
        || epoch.afterSimulationTimeSeconds
            != epoch.beforeSimulationTimeSeconds
                + epoch.structureSettings.timeStepSeconds
        || epoch.sourceCycleStateFingerprint
            != epoch.loadEpoch.sourceCycleStateFingerprint
        || epoch.sourceLoadEpochFingerprint != epoch.loadEpoch.fingerprint
        || epoch.targetDefinitionFingerprint
            != epoch.loadEpoch.targetDefinitionFingerprint
        || epoch.beforeCheckpointFingerprint
            != bytesFingerprint(epoch.beforeStructureCheckpoint)
        || epoch.afterCheckpointFingerprint
            != bytesFingerprint(epoch.afterStructureCheckpoint)
        || !epoch.stepped
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "regional opening wall structural-step integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& target,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
        epoch);
    validateStructureSettings(settings.structure);
    if (settings.structure.timeStepSeconds
            != cycleState.adjustedMomentum.timeStepSeconds
        || !sameStepSettings(epoch.structureSettings, settings.structure)
        || epoch.sourceSettingsFingerprint != settingsFingerprint(settings)
        || epoch.sourceCycleStateFingerprint != cycleState.fingerprint
        || epoch.targetDefinitionFingerprint
            != target.definitionFingerprint()
        || epoch.workingStorageBytes
            != workingStorageBytes(epoch, cycleState)) {
        throw std::invalid_argument(
            "regional opening wall structural-step sources are foreign");
    }
    validateAggregateLimits(epoch, limits);
    validateSceneFluidRegionalOpeningMomentumWallLoadEpoch(
        epoch.loadEpoch, cycleState, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions,
        acceptedOpenings, acceptedResistanceDefinitions,
        acceptedBaseMetric, acceptedMetric, surface, surfaceState,
        transfer, quadrature, settings.load, limits.load);

    const StructureCheckpoint before = deserializeCheckpointOrThrow(
        epoch.beforeStructureCheckpoint, target,
        limits.checkpointPersistence,
        "cannot validate pre-step Structure checkpoint");
    const StructureCheckpoint after = deserializeCheckpointOrThrow(
        epoch.afterStructureCheckpoint, target,
        limits.checkpointPersistence,
        "cannot validate post-step Structure checkpoint");
    validateCheckpointEndpoints(epoch, before, after);

    Structure replay(target.definition());
    replay.restore(before);
    const auto replayLoad =
        applySceneFluidRegionalOpeningMomentumWallLoadEpoch(
            cycleState, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, surface, surfaceState, transfer,
            quadrature, replay, settings.load, limits.load);
    const StructureDiagnostics replayDiagnostics =
        replay.step(settings.structure);
    validateDiagnostics(replayDiagnostics);
    std::vector<std::uint8_t> replayAfterBytes;
    serializeCheckpointOrThrow(
        replay, replay.checkpoint(), replayAfterBytes,
        limits.checkpointPersistence,
        "cannot retain replayed post-step Structure checkpoint");
    if (replayLoad != epoch.loadEpoch
        || replayDiagnostics != epoch.diagnostics
        || replayAfterBytes != epoch.afterStructureCheckpoint) {
        throw std::invalid_argument(
            "regional opening wall structural-step replay differs");
    }
}

} // namespace simwing::fsi
