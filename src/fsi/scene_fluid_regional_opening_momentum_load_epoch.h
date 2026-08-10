#pragma once

#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"
#include "scene_fluid_regional_opening_load_epoch.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumLoadEpochVersion = 1;

struct SceneFluidRegionalOpeningMomentumLoadEpochLimits {
    fluid::PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits
        cycleState;
    SceneFluidRegionalOpeningLoadEpochLimits loadEpoch;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// Transport-provenance wrapper around the established regional opening load
// application. The complete staggered cycle state is validated against its
// current transport and consecutive accepted-pressure epochs before the
// retained-solid pressure load can reach Structure. The nested load epoch owns
// the conservative scene mapping and nodal application.
//
// A Structure checkpoint encloses the whole outer transaction, so even a late
// aggregate/provenance failure after nested application restores the exact
// pending-load state. This endpoint does not step Structure, commit a cycle
// owner, add wall shear, rebase topology, or select a production worker.
struct SceneFluidRegionalOpeningMomentumLoadEpoch {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumLoadEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceCycleStateFingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourceTransportMetricFingerprint = 0;
    std::uint64_t sourceAcceptedMetricFingerprint = 0;
    SceneFluidRegionalOpeningLoadEpoch loadEpoch;
    bool applied = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumLoadEpoch&) const = default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumLoadEpoch
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits = {});

void validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& epoch);

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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits = {});

} // namespace simwing::fsi
