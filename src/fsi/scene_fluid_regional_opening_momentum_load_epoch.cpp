#include "scene_fluid_regional_opening_momentum_load_epoch.h"

#include <limits>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    void integer(std::uint64_t value) {
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }
    void boolean(const bool value) { integer(value ? 1U : 0U); }
    [[nodiscard]] std::uint64_t value() const noexcept { return value_; }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening momentum-load storage overflows");
    }
    return first + second;
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening momentum-load limits are invalid");
    }
}

std::uint64_t epochFingerprint(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& epoch) {
    Fingerprint fingerprint;
    for (const std::uint64_t value : {
             static_cast<std::uint64_t>(epoch.version),
             epoch.sourceCycleStateFingerprint,
             epoch.sourceTransportFingerprint,
             epoch.sourceAcceptedStateFingerprint,
             epoch.sourceTransportMetricFingerprint,
             epoch.sourceAcceptedMetricFingerprint,
             epoch.loadEpoch.fingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.boolean(epoch.applied);
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.workingStorageBytes));
    return fingerprint.value();
}

void validateAggregateLimits(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits) {
    if (epoch.ownedStorageBytes > limits.maximumOwnedBytes
        || epoch.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening momentum-load aggregate limit exceeded");
    }
}

} // namespace

SceneFluidRegionalOpeningMomentumLoadEpoch
applySceneFluidRegionalOpeningMomentumLoadEpoch(
    const fluid::PlanarPressureRegionFragmentOpeningMomentumCycleState&
        cycleState,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        transportVolumeRates,
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
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits) {
    validateLimits(limits);
    const StructureCheckpoint before = target.checkpoint();
    try {
        fluid::validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
            cycleState, transportVolumeRates, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedMetric, limits.cycleState);
        auto loadEpoch = applySceneFluidRegionalOpeningLoadEpoch(
            cycleState.acceptedState, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, surface, surfaceState, transfer,
            quadrature, target, settings, limits.loadEpoch);

        SceneFluidRegionalOpeningMomentumLoadEpoch result;
        result.sourceCycleStateFingerprint = cycleState.fingerprint;
        result.sourceTransportFingerprint = cycleState.transport.fingerprint;
        result.sourceAcceptedStateFingerprint =
            cycleState.acceptedState.fingerprint;
        result.sourceTransportMetricFingerprint =
            cycleState.transportMetricFingerprint;
        result.sourceAcceptedMetricFingerprint =
            cycleState.acceptedMetricFingerprint;
        result.loadEpoch = std::move(loadEpoch);
        result.applied = true;
        result.ownedStorageBytes = result.loadEpoch.ownedStorageBytes;
        result.workingStorageBytes = checkedAdd(
            cycleState.ownedStorageBytes,
            result.loadEpoch.workingStorageBytes);
        validateAggregateLimits(result, limits);
        result.fingerprint = epochFingerprint(result);
        validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(result);
        return result;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& epoch) {
    validateSceneFluidRegionalOpeningLoadEpochIntegrity(epoch.loadEpoch);
    if (epoch.version
            != sceneFluidRegionalOpeningMomentumLoadEpochVersion
        || epoch.fingerprint == 0
        || epoch.sourceCycleStateFingerprint == 0
        || epoch.sourceTransportFingerprint == 0
        || epoch.sourceAcceptedStateFingerprint == 0
        || epoch.sourceTransportMetricFingerprint == 0
        || epoch.sourceAcceptedMetricFingerprint == 0
        || epoch.sourceAcceptedStateFingerprint
            != epoch.loadEpoch.sourceAcceptedStateFingerprint
        || !epoch.applied || !epoch.loadEpoch.applied
        || epoch.ownedStorageBytes != epoch.loadEpoch.ownedStorageBytes
        || epoch.workingStorageBytes
            < epoch.loadEpoch.workingStorageBytes
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "regional opening momentum-load epoch integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumLoadEpoch(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& epoch,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumCycleState&
        cycleState,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        transportVolumeRates,
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
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(epoch);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
        cycleState, transportVolumeRates, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedMetric, limits.cycleState);
    validateSceneFluidRegionalOpeningLoadEpoch(
        epoch.loadEpoch, cycleState.acceptedState,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, surface, surfaceState, transfer,
        quadrature, settings, limits.loadEpoch);
    const std::size_t expectedWorking = checkedAdd(
        cycleState.ownedStorageBytes,
        epoch.loadEpoch.workingStorageBytes);
    if (epoch.sourceCycleStateFingerprint != cycleState.fingerprint
        || epoch.sourceTransportFingerprint
            != cycleState.transport.fingerprint
        || epoch.sourceAcceptedStateFingerprint
            != cycleState.acceptedState.fingerprint
        || epoch.sourceTransportMetricFingerprint
            != transportMetric.fingerprint
        || epoch.sourceAcceptedMetricFingerprint
            != acceptedMetric.fingerprint
        || epoch.workingStorageBytes != expectedWorking) {
        throw std::invalid_argument(
            "regional opening momentum-load epoch sources are foreign");
    }
    validateAggregateLimits(epoch, limits);
}

} // namespace simwing::fsi
