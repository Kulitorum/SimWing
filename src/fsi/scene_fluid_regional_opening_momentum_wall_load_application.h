#pragma once

#include "scene_fluid_regional_opening_momentum_wall_cycle_state.h"
#include "scene_fluid_surface_transfer.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallLoadApplicationVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings {
    ConservativeTransferSettings transfer;
    double absoluteActionReactionToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeActionReactionTolerance = 1.0e-11;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&)
        const = default;
};

struct SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits {
    SceneFluidRegionalOpeningMomentumWallCycleStateLimits cycleState;
    std::size_t maximumNodeLoads = 20'000'000;
    std::size_t maximumStructureNodes = 20'000'000;
    std::size_t maximumOwnedBytes = 2048ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionalOpeningMomentumWallAppliedNodeLoad {
    std::size_t loadIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;
    StructureVector3 priorPendingForceNewtons;
    StructureVector3 appliedWallForceNewtons;
    StructureVector3 resultingPendingForceNewtons;
    StructureVector3 applicationResidualNewtons;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallAppliedNodeLoad&) const =
        default;
};

// Immutable receipt for the opt-in wall-cycle endpoint's conservative
// tangential Structure load. The exact adjusted fluid impulse must close
// against integrated wall traction over the same step before any pending load
// changes. Every nodal mutation is then recorded against the prior Structure
// checkpoint and replayable from the source state, scene quadrature, and
// transfer topology.
//
// This transaction applies wall traction only. It does not apply the retained
// pressure endpoint, step Structure, consume pending loads, or select a worker.
struct SceneFluidRegionalOpeningMomentumWallLoadApplication {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallLoadApplicationVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceCycleStateFingerprint = 0;
    std::uint64_t sourceAdjustmentStateFingerprint = 0;
    std::uint64_t sourceWallTractionFingerprint = 0;
    std::uint64_t sourceWallExchangeFingerprint = 0;
    std::uint64_t transportMetricFingerprint = 0;
    std::uint64_t acceptedMetricFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t targetDefinitionFingerprint = 0;
    std::uint64_t sourceSettingsFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double timeStepSeconds = 0.0;
    std::size_t structureNodeCount = 0;
    std::vector<SceneFluidRegionalOpeningMomentumWallAppliedNodeLoad>
        nodeLoads;
    StructureVector3 fluidAdjustmentImpulseKilogramMetersPerSecond;
    StructureVector3 transferredWallForceNewtons;
    StructureVector3 structureWallImpulseKilogramMetersPerSecond;
    StructureVector3 actionReactionResidualKilogramMetersPerSecond;
    StructureVector3 priorPendingForceNewtons;
    StructureVector3 appliedWallForceNewtons;
    StructureVector3 resultingPendingForceNewtons;
    StructureVector3 applicationResidualNewtons;
    bool applied = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallLoadApplication&) const =
        default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallLoadApplication
applySceneFluidRegionalOpeningMomentumWallLoads(
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
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings = {},
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits =
        {});

void validateSceneFluidRegionalOpeningMomentumWallLoadApplicationIntegrity(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application);

// Re-evaluates the compact cycle state, exact quadrature transfer,
// action/reaction ledger, and every applied node load. Historical pending
// loads remain receipt-owned data and no Structure is mutated.
void validateSceneFluidRegionalOpeningMomentumWallLoadApplication(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application,
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
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings = {},
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits =
        {});

} // namespace simwing::fsi
