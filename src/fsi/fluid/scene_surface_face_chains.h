#pragma once

#include "fluid/scene_surface_face_graph.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFaceChainVersion = 1;

enum class SceneFluidFaceChainKind : std::uint8_t {
    Open = 1,
    Closed = 2,
};

struct SceneFluidFaceChain {
    std::uint64_t stableId = 0;
    std::size_t activeFaceIndex = 0;
    SceneFluidFaceChainKind kind = SceneFluidFaceChainKind::Open;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    std::size_t firstNodeReference = 0;
    std::size_t nodeReferenceCount = 0;
    std::size_t firstSegmentReference = 0;
    std::size_t segmentReferenceCount = 0;
    double lengthMeters = 0.0;
    std::array<std::uint8_t, 2> endpointFaceBoundaryMasks{};
    std::array<bool, 2> endpointOnAuthoredOpening{};

    bool operator==(const SceneFluidFaceChain&) const = default;
};

struct SceneFluidFaceChainLimits {
    std::size_t maximumChains = 5'000'000;
    std::size_t maximumNodeReferences = 10'000'000;
    std::size_t maximumSegmentReferences = 5'000'000;
    std::size_t maximumChainBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFaceChainSet {
    std::uint32_t version = sceneFluidFaceChainVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t faceGraphFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t openChainCount = 0;
    std::size_t closedChainCount = 0;
    std::size_t openingEndpointCount = 0;
    std::size_t gridBoundaryEndpointCount = 0;
    std::vector<SceneFluidFaceChain> chains;
    // Directed node/segment indices into SceneFluidFaceGraph. Closed chains
    // contain one unique node per segment and close implicitly to node zero.
    std::vector<std::size_t> nodeReferences;
    std::vector<std::size_t> segmentReferences;

    bool operator==(const SceneFluidFaceChainSet&) const = default;
};

// Converts each degree-at-most-two face graph into consistently directed open
// chains and closed loops. Direction follows triangle winding and therefore the
// authored negative-to-positive region convention. Branched nodes, conflicting
// directions, or a region-pair change inside one stitched chain are rejected
// rather than silently choosing a fluid partition.
[[nodiscard]] SceneFluidFaceChainSet buildSceneFluidFaceChains(
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
    const SceneFluidFaceChainLimits& limits = {});

void validateSceneFluidFaceChains(
    const SceneFluidFaceChainSet& chains,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph);

} // namespace simwing::fsi::fluid
