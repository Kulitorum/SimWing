#pragma once

#include "fluid/scene_surface_face_chains.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFaceLoopVersion = 1;

struct SceneFluidFaceLoopSettings {
    double minimumAbsoluteAreaSquareMeters = 1.0e-18;
    double intersectionToleranceMeters = 1.0e-12;

    bool operator==(const SceneFluidFaceLoopSettings&) const = default;
};

struct SceneFluidFaceLoop {
    std::uint64_t stableId = 0;
    std::size_t chainIndex = 0;
    std::size_t activeFaceIndex = 0;
    double signedAreaSquareMeters = 0.0;
    double areaSquareMeters = 0.0;
    Vec3 centroidMeters;
    StableId enclosedRegionId = invalidStableId;
    StableId exteriorRegionId = invalidStableId;
    bool positiveSideIsInterior = false;

    bool operator==(const SceneFluidFaceLoop& other) const {
        return stableId == other.stableId
            && chainIndex == other.chainIndex
            && activeFaceIndex == other.activeFaceIndex
            && signedAreaSquareMeters == other.signedAreaSquareMeters
            && areaSquareMeters == other.areaSquareMeters
            && centroidMeters.x == other.centroidMeters.x
            && centroidMeters.y == other.centroidMeters.y
            && centroidMeters.z == other.centroidMeters.z
            && enclosedRegionId == other.enclosedRegionId
            && exteriorRegionId == other.exteriorRegionId
            && positiveSideIsInterior == other.positiveSideIsInterior;
    }
};

struct SceneFluidFaceLoopLimits {
    std::size_t maximumLoops = 5'000'000;
    std::size_t maximumLoopBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFaceLoopSet {
    std::uint32_t version = sceneFluidFaceLoopVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t faceChainFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidFaceLoopSettings settings;
    std::size_t unresolvedOpenChainCount = 0;
    double summedLoopAreaSquareMeters = 0.0;
    std::vector<SceneFluidFaceLoop> loops;

    bool operator==(const SceneFluidFaceLoopSet&) const = default;
};

// Measures simple closed face chains in a right-handed local face chart.
// Winding makes the positive authored side the left side of every directed
// edge, so signed area determines which authored region is enclosed. Open
// chains remain an explicit unresolved count. Separate nested loops remain
// separate; summed area is not a union or a final face partition.
[[nodiscard]] SceneFluidFaceLoopSet buildSceneFluidFaceLoops(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSettings& settings = {},
    const SceneFluidFaceLoopLimits& limits = {});

void validateSceneFluidFaceLoops(
    const SceneFluidFaceLoopSet& loops,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains);

} // namespace simwing::fsi::fluid
