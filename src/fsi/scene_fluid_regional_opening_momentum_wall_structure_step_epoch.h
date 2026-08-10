#pragma once

#include "scene_fluid_regional_opening_momentum_wall_load_epoch.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallStructureStepEpochVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings {
    SceneFluidRegionalOpeningMomentumWallLoadEpochSettings load;
    StructureStepSettings structure;
};

struct SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits {
    SceneFluidRegionalOpeningMomentumWallLoadEpochLimits load;
    StructureCheckpointPersistenceLimits checkpointPersistence;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One opt-in structural macro-step that applies the accepted regional pressure
// and conservative wall traction transaction, then lets XPBD consume the
// resulting pending loads exactly once. The structural and fluid time steps
// must match exactly. Complete deterministic Structure checkpoint encodings
// retain the before/after endpoints for independent replay.
//
// One outer checkpoint encloses load application, stepping, persistence, and
// all late receipt validation. This endpoint does not commit the fluid cycle
// owner, advance fluid momentum, rebase topology, or select a production
// worker.
struct SceneFluidRegionalOpeningMomentumWallStructureStepEpoch {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallStructureStepEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceCycleStateFingerprint = 0;
    std::uint64_t sourceLoadEpochFingerprint = 0;
    std::uint64_t targetDefinitionFingerprint = 0;
    std::uint64_t sourceSettingsFingerprint = 0;
    std::uint64_t beforeCheckpointFingerprint = 0;
    std::uint64_t afterCheckpointFingerprint = 0;
    std::uint64_t beforeAcceptedStepCount = 0;
    std::uint64_t afterAcceptedStepCount = 0;
    double beforeSimulationTimeSeconds = 0.0;
    double afterSimulationTimeSeconds = 0.0;
    SceneFluidRegionalOpeningMomentumWallLoadEpoch loadEpoch;
    StructureStepSettings structureSettings;
    StructureDiagnostics diagnostics;
    std::vector<std::uint8_t> beforeStructureCheckpoint;
    std::vector<std::uint8_t> afterStructureCheckpoint;
    bool stepped = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallStructureStepEpoch
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
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits = {});

void validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch);

// Rebuilds a Structure from the retained pre-step checkpoint, reapplies both
// load paths, advances XPBD once, and requires byte-identical post-step state
// and diagnostics. The supplied target contributes only its trusted
// definition and is not mutated.
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
    const Structure& target,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits = {});

} // namespace simwing::fsi
