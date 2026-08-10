#pragma once

#include "scene_fluid_grid_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_coupled_state.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallPostStepGeometryVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits {
    SceneFluidGridEpochLimits gridEpoch;
    StructureCheckpointPersistenceLimits structureCheckpoint;
    std::size_t maximumSurfaceVertices = 10'000'000;
    std::size_t maximumOwnedBytes = 4ULL * 1024ULL * 1024ULL * 1024ULL;
};

// Immutable authoritative geometry handoff immediately after one accepted
// coupled structural step. It proves that the supplied live Structure is the
// exact post-step endpoint retained by the coupled state, captures the next
// scene surface epoch, and rebuilds every established grid-geometry stage.
//
// This receipt deliberately does not classify pressure topology as stable,
// build control volumes/links, rebase fluid state, run another solve, or select
// a production worker.
struct SceneFluidRegionalOpeningMomentumWallPostStepGeometry {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallPostStepGeometryVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceCoupledStateFingerprint = 0;
    std::uint64_t sourceStructureStepFingerprint = 0;
    std::uint64_t sourcePostStepCheckpointFingerprint = 0;
    std::uint64_t sourceSurfaceStateFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    SceneFluidSurfaceState surfaceState;
    SceneFluidGridEpoch gridEpoch;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallPostStepGeometry&) const =
        default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallPostStepGeometry
buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits& limits =
        {});

void validateSceneFluidRegionalOpeningMomentumWallPostStepGeometryIntegrity(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry);

void validateSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry,
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits& limits =
        {});

} // namespace simwing::fsi
