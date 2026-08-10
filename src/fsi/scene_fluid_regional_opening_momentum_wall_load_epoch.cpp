#include "scene_fluid_regional_opening_momentum_wall_load_epoch.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
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

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

double maximumAbsolute(const StructureVector3& value) {
    return std::max({std::abs(value.x),
                     std::abs(value.y),
                     std::abs(value.z)});
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening momentum wall-load epoch storage overflows");
    }
    return first + second;
}

void fingerprintVector(Fingerprint& fingerprint,
                       const StructureVector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening momentum wall-load epoch limits are invalid");
    }
}

SceneFluidRegionalOpeningLoadEpochSettings pressureSettings(
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings) {
    SceneFluidRegionalOpeningLoadEpochSettings result;
    result.pressureState = settings.pressureState;
    result.transfer = settings.transfer;
    return result;
}

SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings wallSettings(
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings) {
    SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings result;
    result.transfer = settings.transfer;
    result.absoluteActionReactionToleranceKilogramMetersPerSecond =
        settings
            .absoluteWallActionReactionToleranceKilogramMetersPerSecond;
    result.relativeActionReactionTolerance =
        settings.relativeWallActionReactionTolerance;
    return result;
}

SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits wallLimits(
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits) {
    auto result = limits.wallLoad;
    result.cycleState = limits.cycleState;
    return result;
}

std::uint64_t settingsFingerprint(
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings) {
    Fingerprint fingerprint;
    for (const double value : {
             settings.pressureState
                 .absolutePressureResidualTolerancePascals,
             settings.pressureState.relativePressureResidualTolerance,
             settings.pressureState
                 .absoluteForceResidualToleranceNewtons,
             settings.pressureState.relativeForceResidualTolerance,
             settings.pressureState
                 .absoluteWorkResidualToleranceJoules,
             settings.pressureState.relativeWorkResidualTolerance,
             settings.transfer.momentReferenceMeters.x,
             settings.transfer.momentReferenceMeters.y,
             settings.transfer.momentReferenceMeters.z,
             settings.transfer.minimumTriangleAreaSquareMeters,
             settings.transfer.minimumQuadratureAreaSquareMeters,
             settings.transfer.barycentricTolerance,
             settings
                 .absoluteWallActionReactionToleranceKilogramMetersPerSecond,
             settings.relativeWallActionReactionTolerance}) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch) {
    return checkedAdd(
        epoch.pressureLoad.ownedStorageBytes,
        epoch.wallLoad.ownedStorageBytes);
}

std::size_t workingStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState) {
    return checkedAdd(
        cycleState.ownedStorageBytes,
        checkedAdd(
            epoch.pressureLoad.workingStorageBytes,
            epoch.wallLoad.workingStorageBytes));
}

std::uint64_t epochFingerprint(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch) {
    Fingerprint fingerprint;
    fingerprint.integer(epoch.version);
    for (const std::uint64_t value : {
             epoch.sourceCycleStateFingerprint,
             epoch.sourceAdjustmentStateFingerprint,
             epoch.sourceAcceptedPressureFingerprint,
             epoch.sourceWallTractionFingerprint,
             epoch.transportMetricFingerprint,
             epoch.acceptedMetricFingerprint,
             epoch.surfaceStateFingerprint,
             epoch.quadratureFingerprint,
             epoch.couplingSurfaceFingerprint,
             epoch.targetDefinitionFingerprint,
             epoch.sourceSettingsFingerprint,
             epoch.acceptedStepCount,
             epoch.pressureLoad.fingerprint,
             epoch.wallLoad.fingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.real(epoch.simulationTimeSeconds);
    fingerprintVector(fingerprint, epoch.priorPendingForceNewtons);
    fingerprintVector(fingerprint, epoch.appliedPressureForceNewtons);
    fingerprintVector(fingerprint, epoch.appliedWallForceNewtons);
    fingerprintVector(fingerprint, epoch.combinedAppliedForceNewtons);
    fingerprintVector(fingerprint, epoch.resultingPendingForceNewtons);
    fingerprintVector(fingerprint, epoch.applicationResidualNewtons);
    fingerprint.boolean(epoch.applied);
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.workingStorageBytes));
    return fingerprint.value();
}

void validateAggregateLimits(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits) {
    if (epoch.ownedStorageBytes > limits.maximumOwnedBytes
        || epoch.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening momentum wall-load epoch aggregate limit exceeded");
    }
}

void validateSequentialHandoff(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch) {
    const auto& pressure = epoch.pressureLoad.application;
    const auto& wall = epoch.wallLoad;
    if (pressure.nodeLoads.size() != wall.nodeLoads.size()
        || pressure.resultingPendingForceNewtons
            != wall.priorPendingForceNewtons) {
        throw std::invalid_argument(
            "regional opening momentum wall-load aggregate handoff is invalid");
    }
    for (std::size_t index = 0;
         index < pressure.nodeLoads.size(); ++index) {
        const auto& pressureLoad = pressure.nodeLoads[index];
        const auto& wallLoad = wall.nodeLoads[index];
        if (pressureLoad.loadIndex != wallLoad.loadIndex
            || pressureLoad.stableId != wallLoad.stableId
            || pressureLoad.structureNode != wallLoad.structureNode
            || pressureLoad.resultingPendingForceNewtons
                != wallLoad.priorPendingForceNewtons) {
            throw std::invalid_argument(
                "regional opening momentum wall-load nodal handoff is invalid");
        }
    }
}

void validateSourceState(
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
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallCycleState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits);
}

} // namespace

SceneFluidRegionalOpeningMomentumWallLoadEpoch
applySceneFluidRegionalOpeningMomentumWallLoadEpoch(
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
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits) {
    validateLimits(limits);
    validateSourceState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits.cycleState);

    const StructureCheckpoint before = target.checkpoint();
    try {
        auto pressureLoad = applySceneFluidRegionalOpeningLoadEpoch(
            cycleState.acceptedPressure, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, surface, surfaceState, transfer,
            quadrature, target, pressureSettings(settings),
            limits.pressureLoad);
        auto wallLoad = applySceneFluidRegionalOpeningMomentumWallLoads(
            cycleState, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, surfaceState, transfer, quadrature, target,
            wallSettings(settings), wallLimits(limits));

        SceneFluidRegionalOpeningMomentumWallLoadEpoch result;
        result.sourceCycleStateFingerprint = cycleState.fingerprint;
        result.sourceAdjustmentStateFingerprint =
            cycleState.adjustedMomentum.fingerprint;
        result.sourceAcceptedPressureFingerprint =
            cycleState.acceptedPressure.fingerprint;
        result.sourceWallTractionFingerprint =
            cycleState.wallTractions.fingerprint;
        result.transportMetricFingerprint =
            cycleState.transportMetricFingerprint;
        result.acceptedMetricFingerprint =
            cycleState.acceptedMetricFingerprint;
        result.surfaceStateFingerprint = surfaceState.fingerprint;
        result.quadratureFingerprint = quadrature.fingerprint;
        result.couplingSurfaceFingerprint =
            transfer.couplingSurfaceFingerprint();
        result.targetDefinitionFingerprint =
            transfer.targetDefinitionFingerprint();
        result.sourceSettingsFingerprint = settingsFingerprint(settings);
        result.acceptedStepCount = surfaceState.acceptedStepCount;
        result.simulationTimeSeconds = surfaceState.simulationTimeSeconds;
        result.pressureLoad = std::move(pressureLoad);
        result.wallLoad = std::move(wallLoad);
        result.priorPendingForceNewtons =
            result.pressureLoad.application.priorPendingForceNewtons;
        result.appliedPressureForceNewtons =
            result.pressureLoad.application.appliedPressureForceNewtons;
        result.appliedWallForceNewtons =
            result.wallLoad.appliedWallForceNewtons;
        result.combinedAppliedForceNewtons = add(
            result.appliedPressureForceNewtons,
            result.appliedWallForceNewtons);
        result.resultingPendingForceNewtons =
            result.wallLoad.resultingPendingForceNewtons;
        result.applicationResidualNewtons = subtract(
            result.resultingPendingForceNewtons,
            add(result.priorPendingForceNewtons,
                result.combinedAppliedForceNewtons));
        result.applied = true;
        result.ownedStorageBytes = ownedStorageBytes(result);
        result.workingStorageBytes =
            workingStorageBytes(result, cycleState);
        validateSequentialHandoff(result);
        validateAggregateLimits(result, limits);
        result.fingerprint = epochFingerprint(result);
        validateSceneFluidRegionalOpeningMomentumWallLoadEpoch(
            result, cycleState, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedBaseMetric, acceptedMetric, surface, surfaceState,
            transfer, quadrature, settings, limits);
        return result;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalOpeningMomentumWallLoadEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch) {
    validateSceneFluidRegionalOpeningLoadEpochIntegrity(epoch.pressureLoad);
    validateSceneFluidRegionalOpeningMomentumWallLoadApplicationIntegrity(
        epoch.wallLoad);
    if (epoch.version
            != sceneFluidRegionalOpeningMomentumWallLoadEpochVersion
        || epoch.fingerprint == 0
        || epoch.sourceCycleStateFingerprint == 0
        || epoch.sourceAdjustmentStateFingerprint == 0
        || epoch.sourceAcceptedPressureFingerprint == 0
        || epoch.sourceWallTractionFingerprint == 0
        || epoch.transportMetricFingerprint == 0
        || epoch.acceptedMetricFingerprint == 0
        || epoch.surfaceStateFingerprint == 0
        || epoch.quadratureFingerprint == 0
        || epoch.couplingSurfaceFingerprint == 0
        || epoch.targetDefinitionFingerprint == 0
        || epoch.sourceSettingsFingerprint == 0
        || !std::isfinite(epoch.simulationTimeSeconds)
        || !finite(epoch.priorPendingForceNewtons)
        || !finite(epoch.appliedPressureForceNewtons)
        || !finite(epoch.appliedWallForceNewtons)
        || !finite(epoch.combinedAppliedForceNewtons)
        || !finite(epoch.resultingPendingForceNewtons)
        || !finite(epoch.applicationResidualNewtons)
        || epoch.sourceAcceptedPressureFingerprint
            != epoch.pressureLoad.sourceAcceptedStateFingerprint
        || epoch.sourceCycleStateFingerprint
            != epoch.wallLoad.sourceCycleStateFingerprint
        || epoch.sourceAdjustmentStateFingerprint
            != epoch.wallLoad.sourceAdjustmentStateFingerprint
        || epoch.sourceWallTractionFingerprint
            != epoch.wallLoad.sourceWallTractionFingerprint
        || epoch.transportMetricFingerprint
            != epoch.wallLoad.transportMetricFingerprint
        || epoch.acceptedMetricFingerprint
            != epoch.wallLoad.acceptedMetricFingerprint
        || epoch.surfaceStateFingerprint
            != epoch.pressureLoad.surfaceStateFingerprint
        || epoch.surfaceStateFingerprint
            != epoch.wallLoad.surfaceStateFingerprint
        || epoch.quadratureFingerprint
            != epoch.pressureLoad.quadratureFingerprint
        || epoch.quadratureFingerprint
            != epoch.wallLoad.quadratureFingerprint
        || epoch.couplingSurfaceFingerprint
            != epoch.pressureLoad.couplingSurfaceFingerprint
        || epoch.couplingSurfaceFingerprint
            != epoch.wallLoad.couplingSurfaceFingerprint
        || epoch.targetDefinitionFingerprint
            != epoch.pressureLoad.targetDefinitionFingerprint
        || epoch.targetDefinitionFingerprint
            != epoch.wallLoad.targetDefinitionFingerprint
        || epoch.acceptedStepCount != epoch.pressureLoad.acceptedStepCount
        || epoch.acceptedStepCount != epoch.wallLoad.acceptedStepCount
        || epoch.simulationTimeSeconds
            != epoch.pressureLoad.simulationTimeSeconds
        || epoch.simulationTimeSeconds
            != epoch.wallLoad.simulationTimeSeconds
        || epoch.priorPendingForceNewtons
            != epoch.pressureLoad.application.priorPendingForceNewtons
        || epoch.appliedPressureForceNewtons
            != epoch.pressureLoad.application.appliedPressureForceNewtons
        || epoch.appliedWallForceNewtons
            != epoch.wallLoad.appliedWallForceNewtons
        || epoch.combinedAppliedForceNewtons
            != add(epoch.appliedPressureForceNewtons,
                   epoch.appliedWallForceNewtons)
        || epoch.resultingPendingForceNewtons
            != epoch.wallLoad.resultingPendingForceNewtons
        || epoch.applicationResidualNewtons
            != subtract(
                epoch.resultingPendingForceNewtons,
                add(epoch.priorPendingForceNewtons,
                    epoch.combinedAppliedForceNewtons))
        || !epoch.applied || !epoch.pressureLoad.applied
        || !epoch.wallLoad.applied
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.workingStorageBytes
            < checkedAdd(
                epoch.pressureLoad.workingStorageBytes,
                epoch.wallLoad.workingStorageBytes)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "regional opening momentum wall-load epoch integrity is invalid");
    }
    validateSequentialHandoff(epoch);
    const double forceScale = std::max({
        maximumAbsolute(epoch.priorPendingForceNewtons),
        maximumAbsolute(epoch.combinedAppliedForceNewtons),
        maximumAbsolute(epoch.resultingPendingForceNewtons), 1.0});
    const double tolerance =
        16384.0 * std::numeric_limits<double>::epsilon() * forceScale;
    if (maximumAbsolute(epoch.applicationResidualNewtons) > tolerance) {
        throw std::invalid_argument(
            "regional opening momentum wall-load epoch closure is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallLoadEpoch(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch,
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
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallLoadEpochIntegrity(epoch);
    validateSourceState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits.cycleState);
    validateSceneFluidRegionalOpeningLoadEpoch(
        epoch.pressureLoad, cycleState.acceptedPressure,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions,
        acceptedOpenings, acceptedResistanceDefinitions, surface,
        surfaceState, transfer, quadrature, pressureSettings(settings),
        limits.pressureLoad);
    validateSceneFluidRegionalOpeningMomentumWallLoadApplication(
        epoch.wallLoad, cycleState, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions,
        acceptedOpenings, acceptedResistanceDefinitions,
        acceptedBaseMetric, acceptedMetric, surfaceState, transfer,
        quadrature, wallSettings(settings), wallLimits(limits));

    const std::size_t expectedWorking =
        workingStorageBytes(epoch, cycleState);
    if (epoch.sourceCycleStateFingerprint != cycleState.fingerprint
        || epoch.sourceAdjustmentStateFingerprint
            != cycleState.adjustedMomentum.fingerprint
        || epoch.sourceAcceptedPressureFingerprint
            != cycleState.acceptedPressure.fingerprint
        || epoch.sourceWallTractionFingerprint
            != cycleState.wallTractions.fingerprint
        || epoch.transportMetricFingerprint != transportMetric.fingerprint
        || epoch.acceptedMetricFingerprint != acceptedMetric.fingerprint
        || epoch.surfaceStateFingerprint != surfaceState.fingerprint
        || epoch.quadratureFingerprint != quadrature.fingerprint
        || epoch.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || epoch.targetDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || epoch.sourceSettingsFingerprint != settingsFingerprint(settings)
        || epoch.acceptedStepCount != surfaceState.acceptedStepCount
        || epoch.simulationTimeSeconds != surfaceState.simulationTimeSeconds
        || epoch.workingStorageBytes != expectedWorking) {
        throw std::invalid_argument(
            "regional opening momentum wall-load epoch sources are foreign");
    }
    validateAggregateLimits(epoch, limits);
}

} // namespace simwing::fsi
