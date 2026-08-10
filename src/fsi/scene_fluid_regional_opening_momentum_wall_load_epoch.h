#pragma once

#include "scene_fluid_regional_opening_load_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_load_application.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallLoadEpochVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallLoadEpochSettings {
    fluid::PlanarPressureRegionFragmentOpeningPressureStateSettings
        pressureState;
    ConservativeTransferSettings transfer;
    double absoluteWallActionReactionToleranceKilogramMetersPerSecond =
        1.0e-12;
    double relativeWallActionReactionTolerance = 1.0e-11;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings&) const =
        default;
};

struct SceneFluidRegionalOpeningMomentumWallLoadEpochLimits {
    SceneFluidRegionalOpeningMomentumWallCycleStateLimits cycleState;
    SceneFluidRegionalOpeningLoadEpochLimits pressureLoad;
    SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits wallLoad;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One explicit accepted wall-cycle endpoint to pending Structure-load
// transaction. The established retained-solid pressure epoch is applied first;
// the conservative wall application then starts from that exact per-node
// pending-load endpoint. One outer Structure checkpoint encloses both nested
// mutations and all late aggregate/lineage validation.
//
// This endpoint does not step Structure, consume pending loads, commit the
// cycle owner, advance fluid momentum, or select a production worker.
struct SceneFluidRegionalOpeningMomentumWallLoadEpoch {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallLoadEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceCycleStateFingerprint = 0;
    std::uint64_t sourceAdjustmentStateFingerprint = 0;
    std::uint64_t sourceAcceptedPressureFingerprint = 0;
    std::uint64_t sourceWallTractionFingerprint = 0;
    std::uint64_t transportMetricFingerprint = 0;
    std::uint64_t acceptedMetricFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t targetDefinitionFingerprint = 0;
    std::uint64_t sourceSettingsFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidRegionalOpeningLoadEpoch pressureLoad;
    SceneFluidRegionalOpeningMomentumWallLoadApplication wallLoad;
    StructureVector3 priorPendingForceNewtons;
    StructureVector3 appliedPressureForceNewtons;
    StructureVector3 appliedWallForceNewtons;
    StructureVector3 combinedAppliedForceNewtons;
    StructureVector3 resultingPendingForceNewtons;
    StructureVector3 applicationResidualNewtons;
    bool applied = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallLoadEpoch&) const =
        default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallLoadEpoch
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
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
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings =
        {},
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits = {});

void validateSceneFluidRegionalOpeningMomentumWallLoadEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch);

// Rebuilds both nested applications, proves the exact per-node pressure-to-wall
// pending-load handoff, and validates the combined aggregate ledger without
// mutating Structure.
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
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
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings =
        {},
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits = {});

} // namespace simwing::fsi
