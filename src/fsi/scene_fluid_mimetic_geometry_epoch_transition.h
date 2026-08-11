#pragma once

#include "scene_fluid_mimetic_geometry_epoch.h"
#include "scene_fluid_pressure_topology_transition.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticGeometryEpochTransitionVersion = 1;

struct SceneFluidMimeticGeometryEpochTransitionLimits {
    SceneFluidMimeticGeometryEpochLimits geometryEpoch;
    SceneFluidPressureTopologyTransitionLimits topologyTransition;
    std::size_t maximumOwnedBytes = 12ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One consecutive graph-free mixed-hybrid geometry transaction. The complete
// current epoch is rebuilt from authoritative Structure motion and must equal
// the independently accepted current grid epoch. Stable pressure-control IDs
// are then classified through the shared topology-transition policy. No flow
// state is rebased and no pressure solve or Structure load occurs here.
struct SceneFluidMimeticGeometryEpochTransition {
    std::uint32_t version =
        sceneFluidMimeticGeometryEpochTransitionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t previousGeometryEpochFingerprint = 0;
    std::uint64_t acceptedCurrentGridEpochFingerprint = 0;
    std::uint64_t currentGeometryEpochFingerprint = 0;
    std::uint64_t topologyTransitionFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousSurfaceStateFingerprint = 0;
    std::uint64_t currentSurfaceStateFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    bool controlVolumeTopologyStable = false;
    SceneFluidMimeticGeometryEpoch currentGeometryEpoch;
    SceneFluidPressureTopologyTransition topologyTransition;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const SceneFluidMimeticGeometryEpochTransition&) const = default;
};

[[nodiscard]] SceneFluidMimeticGeometryEpochTransition
buildSceneFluidMimeticGeometryEpochTransition(
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidMimeticGeometryEpochSettings& settings = {},
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits = {});

void validateSceneFluidMimeticGeometryEpochTransitionIntegrity(
    const SceneFluidMimeticGeometryEpochTransition& transition);

void validateSceneFluidMimeticGeometryEpochTransition(
    const SceneFluidMimeticGeometryEpochTransition& transition,
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidMimeticGeometryEpochSettings& settings = {},
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits = {});

} // namespace simwing::fsi
