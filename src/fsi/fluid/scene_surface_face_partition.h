#pragma once

#include "fluid/scene_surface_face_loops.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFacePartitionVersion = 3;
inline constexpr std::size_t noParentFaceLoop = static_cast<std::size_t>(-1);

enum class SceneFluidFacePartitionKind : std::uint8_t {
    ClosedLoops = 1,
    BoundaryOpenChain = 2,
    BoundaryChainArrangement = 3,
};

struct SceneFluidFacePartitionSettings {
    double geometryToleranceMeters = 1.0e-12;
    double areaClosureToleranceSquareMeters = 1.0e-12;

    bool operator==(const SceneFluidFacePartitionSettings&) const = default;
};

struct SceneFluidFaceLoopContainment {
    std::size_t loopIndex = 0;
    std::size_t parentLoopIndex = noParentFaceLoop;
    std::size_t depth = 0;

    bool operator==(const SceneFluidFaceLoopContainment&) const = default;
};

struct SceneFluidFaceRegionArea {
    StableId regionId = invalidStableId;
    double areaSquareMeters = 0.0;

    bool operator==(const SceneFluidFaceRegionArea&) const = default;
};

struct SceneFluidFacePartition {
    std::uint64_t stableId = 0;
    std::size_t activeFaceIndex = 0;
    SceneFluidFacePartitionKind kind =
        SceneFluidFacePartitionKind::ClosedLoops;
    StableId rootExteriorRegionId = invalidStableId;
    std::size_t firstLoopReference = 0;
    std::size_t loopReferenceCount = 0;
    std::size_t firstOpenChainReference = 0;
    std::size_t openChainReferenceCount = 0;
    std::size_t firstRegionArea = 0;
    std::size_t regionAreaCount = 0;
    double faceAreaSquareMeters = 0.0;
    double assignedAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;

    bool operator==(const SceneFluidFacePartition&) const = default;
};

struct SceneFluidFacePartitionLimits {
    std::size_t maximumPartitions = 5'000'000;
    std::size_t maximumReferences = 15'000'000;
    std::size_t maximumSegmentPairTests = 50'000'000;
    std::size_t maximumPartitionBytes = 768ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFacePartitionSet {
    std::uint32_t version = sceneFluidFacePartitionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t faceLoopFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidFacePartitionSettings settings;
    std::size_t unresolvedActiveFaceCount = 0;
    std::size_t segmentPairTestCount = 0;
    std::vector<SceneFluidFaceLoopContainment> loopContainment;
    std::vector<SceneFluidFacePartition> partitions;
    std::vector<std::size_t> loopReferences;
    std::vector<std::size_t> openChainReferences;
    std::vector<SceneFluidFaceRegionArea> regionAreas;

    bool operator==(const SceneFluidFacePartitionSet&) const = default;
};

// Builds containment and exact region-area accounting for active faces whose
// interface consists solely of non-touching simple closed loops, or of one
// simple directed open-chain arrangement whose leaves all lie on the rectangular
// face boundary. Parent/child authored regions must form a continuous nesting
// chain. Opening-ended chains, coplanar sheet area, boundary-touching closed
// loops, or absent interfaces stay explicitly unresolved. This is a face
// partition, not a cut-cell volume.
[[nodiscard]] SceneFluidFacePartitionSet buildSceneFluidFacePartitions(
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
    const SceneFluidFaceLoopSet& loops,
    const SceneFluidFacePartitionSettings& settings = {},
    const SceneFluidFacePartitionLimits& limits = {});

void validateSceneFluidFacePartitions(
    const SceneFluidFacePartitionSet& partitions,
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
    const SceneFluidFaceLoopSet& loops);

} // namespace simwing::fsi::fluid
