#pragma once

#include "fluid/scene_surface_face_partition.h"
#include "scene_fluid_quadrature.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidGridEpochVersion = 1;

struct SceneFluidGridEpochSettings {
    fluid::SceneFluidGridCandidateSettings candidates;
    fluid::SceneFluidGridIntersectionSettings intersections;
    fluid::SceneFluidFaceGraphSettings faceGraph;
    fluid::SceneFluidFaceLoopSettings faceLoops;
    fluid::SceneFluidFacePartitionSettings facePartitions;

    bool operator==(const SceneFluidGridEpochSettings&) const = default;
};

struct SceneFluidGridEpochLimits {
    fluid::SceneFluidGridCandidateLimits candidates;
    fluid::SceneFluidGridIntersectionLimits intersections;
    fluid::SceneFluidGridPatchLimits patches;
    fluid::SceneFluidPatchOwnershipLimits ownership;
    fluid::SceneFluidFaceCrossingLimits crossings;
    fluid::SceneFluidFaceTopologyLimits faceTopology;
    fluid::SceneFluidFaceGraphLimits faceGraph;
    fluid::SceneFluidFaceChainLimits faceChains;
    fluid::SceneFluidFaceLoopLimits faceLoops;
    fluid::SceneFluidFacePartitionLimits facePartitions;
    // Sum of vector payload storage owned by the completed epoch. Individual
    // stages retain their tighter count/byte limits as well.
    std::size_t maximumEpochBytes = 2ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One immutable, accepted Structure-surface remap onto a Cartesian grid. It
// composes every currently implemented geometry stage and unique conservative
// quadrature under one epoch identity. Face partitions may remain explicitly
// unresolved; no cut-cell volume, fluid-region reconstruction, pressure, or
// moving-boundary solve is invented here.
struct SceneFluidGridEpoch {
    std::uint32_t version = sceneFluidGridEpochVersion;
    std::uint64_t fingerprint = 0;
    std::size_t ownedStorageBytes = 0;
    fluid::SceneFluidGridCandidateSet candidates;
    fluid::SceneFluidGridIntersectionSet intersections;
    fluid::SceneFluidGridPatchSet patches;
    fluid::SceneFluidPatchOwnership ownership;
    fluid::SceneFluidFaceCrossingSet crossings;
    fluid::SceneFluidFaceTopology faceTopology;
    fluid::SceneFluidFaceGraph faceGraph;
    fluid::SceneFluidFaceChainSet faceChains;
    fluid::SceneFluidFaceLoopSet faceLoops;
    fluid::SceneFluidFacePartitionSet facePartitions;
    SceneFluidQuadratureDefinition quadrature;

    bool operator==(const SceneFluidGridEpoch&) const = default;
};

// Every stage is built into a local candidate and fully validated before the
// epoch is returned. Failure therefore cannot publish a partially mixed epoch.
[[nodiscard]] SceneFluidGridEpoch buildSceneFluidGridEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpochSettings& settings = {},
    const SceneFluidGridEpochLimits& limits = {});

void validateSceneFluidGridEpoch(
    const SceneFluidGridEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer);

} // namespace simwing::fsi
