#pragma once

#include "fluid/scene_surface_face_topology.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFaceGraphVersion = 1;

enum class SceneFluidFaceNodeKind : std::uint8_t {
    SurfaceVertex = 1,
    SurfaceEdge = 2,
    GridEdge = 3,
};

enum SceneFluidFaceBoundary : std::uint8_t {
    FaceBoundaryNone = 0,
    FaceBoundaryUMinus = 1U << 0U,
    FaceBoundaryUPlus = 1U << 1U,
    FaceBoundaryVMinus = 1U << 2U,
    FaceBoundaryVPlus = 1U << 3U,
};

struct SceneFluidFaceGraphSettings {
    double endpointToleranceMeters = 1.0e-12;
    double barycentricZeroTolerance = 1.0e-12;

    bool operator==(const SceneFluidFaceGraphSettings&) const = default;
};

struct SceneFluidFaceGraphNode {
    std::uint64_t stableId = 0;
    std::size_t activeFaceIndex = 0;
    SceneFluidFaceNodeKind kind = SceneFluidFaceNodeKind::SurfaceEdge;
    Vec3 positionMeters;
    StableId sourceVertexId = invalidStableId;
    std::array<StableId, 2> sourceEdgeVertexIds{};
    StableId sourceTriangleId = invalidStableId;
    std::uint8_t faceBoundaryMask = FaceBoundaryNone;
    bool authoredOpeningBoundary = false;
    std::size_t firstIncidentSegmentReference = 0;
    std::size_t incidentSegmentReferenceCount = 0;

    bool operator==(const SceneFluidFaceGraphNode& other) const {
        return stableId == other.stableId
            && activeFaceIndex == other.activeFaceIndex
            && kind == other.kind
            && positionMeters.x == other.positionMeters.x
            && positionMeters.y == other.positionMeters.y
            && positionMeters.z == other.positionMeters.z
            && sourceVertexId == other.sourceVertexId
            && sourceEdgeVertexIds == other.sourceEdgeVertexIds
            && sourceTriangleId == other.sourceTriangleId
            && faceBoundaryMask == other.faceBoundaryMask
            && authoredOpeningBoundary == other.authoredOpeningBoundary
            && firstIncidentSegmentReference
                == other.firstIncidentSegmentReference
            && incidentSegmentReferenceCount
                == other.incidentSegmentReferenceCount;
    }
};

struct SceneFluidFaceGraphSegment {
    std::uint64_t stableId = 0;
    std::size_t activeFaceIndex = 0;
    std::size_t crossingIndex = 0;
    std::array<std::size_t, 2> nodeIndices{};

    bool operator==(const SceneFluidFaceGraphSegment&) const = default;
};

struct SceneFluidFaceGraphRange {
    std::size_t activeFaceIndex = 0;
    std::size_t firstNode = 0;
    std::size_t nodeCount = 0;
    std::size_t firstSegment = 0;
    std::size_t segmentCount = 0;

    bool operator==(const SceneFluidFaceGraphRange&) const = default;
};

struct SceneFluidFaceGraphLimits {
    std::size_t maximumNodes = 10'000'000;
    std::size_t maximumSegments = 5'000'000;
    std::size_t maximumIncidentReferences = 10'000'000;
    std::size_t maximumGraphBytes = 768ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFaceGraph {
    std::uint32_t version = sceneFluidFaceGraphVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t faceTopologyFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidFaceGraphSettings settings;
    std::size_t degreeOneNodeCount = 0;
    std::size_t degreeTwoNodeCount = 0;
    std::size_t higherDegreeNodeCount = 0;
    std::size_t authoredOpeningNodeCount = 0;
    std::vector<SceneFluidFaceGraphRange> faceRanges;
    std::vector<SceneFluidFaceGraphNode> nodes;
    std::vector<SceneFluidFaceGraphSegment> segments;
    std::vector<std::size_t> incidentSegmentReferences;

    bool operator==(const SceneFluidFaceGraph&) const = default;
};

// Stitches triangle-local transverse segments using topological provenance:
// stable scene vertices and edges are shared nodes, while endpoints clipped by
// the rectangular MAC-face boundary remain triangle/grid-edge nodes, with
// positions canonicalized from the authored triangle and exact grid planes.
// This does not infer fluid regions or require closed chains; degree
// diagnostics expose
// mesh boundaries, openings, and unresolved junctions for the next stage.
[[nodiscard]] SceneFluidFaceGraph buildSceneFluidFaceGraph(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraphSettings& settings = {},
    const SceneFluidFaceGraphLimits& limits = {});

void validateSceneFluidFaceGraph(
    const SceneFluidFaceGraph& graph,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology);

} // namespace simwing::fsi::fluid
