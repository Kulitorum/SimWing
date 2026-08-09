#pragma once

#include "scene_fluid_surface.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidOpeningCapVersion = 1;

struct SceneFluidOpeningCapSettings {
    double planarityToleranceMeters = 1.0e-10;
    double minimumTriangleAreaSquareMeters = 1.0e-18;
    // Compatibility name retained from the convex-only v1 subset. This is the
    // scale-relative orientation tolerance for simple-loop validation and
    // deterministic triangulation as well as convexity classification.
    double convexityTolerance = 1.0e-12;
    // A directed unresolved material-boundary loop may remain cap-free only
    // while both reference and accepted geometry have area no larger than
    // this fraction of the loop's squared extent. This represents a closed
    // collapsed wingtip, not a finite aperture.
    double collapsedBoundaryRelativeAreaTolerance = 1.0e-6;

    bool operator==(const SceneFluidOpeningCapSettings&) const = default;
};

struct SceneFluidOpeningCapLimits {
    std::size_t maximumBoundaryEdges = 20'000'000;
    std::size_t maximumCaps = 1'000'000;
    std::size_t maximumCapTriangles = 10'000'000;
    std::size_t maximumCapBytes = 512ULL * 1024ULL * 1024ULL;
    std::size_t maximumBoundaryIntersectionTests = 100'000'000;
    std::size_t maximumTriangulationPointTests = 100'000'000;
};

struct SceneFluidOpeningCapTriangle {
    std::size_t openingIndex = 0;
    std::size_t triangleOrdinal = 0;
    std::array<std::size_t, 3> vertexIndices{};
    Vec3 unitNormalNegativeToPositive;
    Vec3 centroidMeters;
    double areaSquareMeters = 0.0;

    bool operator==(const SceneFluidOpeningCapTriangle& other) const {
        return openingIndex == other.openingIndex
            && triangleOrdinal == other.triangleOrdinal
            && vertexIndices == other.vertexIndices
            && unitNormalNegativeToPositive.x
                == other.unitNormalNegativeToPositive.x
            && unitNormalNegativeToPositive.y
                == other.unitNormalNegativeToPositive.y
            && unitNormalNegativeToPositive.z
                == other.unitNormalNegativeToPositive.z
            && centroidMeters.x == other.centroidMeters.x
            && centroidMeters.y == other.centroidMeters.y
            && centroidMeters.z == other.centroidMeters.z
            && areaSquareMeters == other.areaSquareMeters;
    }
};

struct SceneFluidOpeningCap {
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    std::size_t negativeSideRegionIndex = 0;
    std::size_t positiveSideRegionIndex = 0;
    OpeningRole role = OpeningRole::Intake;
    std::size_t firstTriangle = 0;
    std::size_t triangleCount = 0;
    Vec3 unitNormalNegativeToPositive;
    Vec3 centroidMeters;
    double areaSquareMeters = 0.0;

    bool operator==(const SceneFluidOpeningCap& other) const {
        return openingIndex == other.openingIndex
            && openingId == other.openingId
            && negativeSideRegionIndex == other.negativeSideRegionIndex
            && positiveSideRegionIndex == other.positiveSideRegionIndex
            && role == other.role
            && firstTriangle == other.firstTriangle
            && triangleCount == other.triangleCount
            && unitNormalNegativeToPositive.x
                == other.unitNormalNegativeToPositive.x
            && unitNormalNegativeToPositive.y
                == other.unitNormalNegativeToPositive.y
            && unitNormalNegativeToPositive.z
                == other.unitNormalNegativeToPositive.z
            && centroidMeters.x == other.centroidMeters.x
            && centroidMeters.y == other.centroidMeters.y
            && centroidMeters.z == other.centroidMeters.z
            && areaSquareMeters == other.areaSquareMeters;
    }
};

// A virtual cap is fluid topology, not fabric. Material and cap incidences
// must form one oriented closed region cycle around every finite-area edge;
// this includes consistently wound three-region sheet/cap junctions. A closed
// material boundary may remain cap-free only while it is collapsed in both
// reference and accepted geometry. Planar loops may omit cap triangles:
// convex loops then retain the exact fan and concave loops receive
// deterministic reference-geometry ear clipping. A nonplanar loop must author
// one oriented boundary-vertex disk; its individual triangle normals survive
// accepted motion. Cap winding is derived from an adjacent fabric edge with
// the same region pair. No cap enters Structure or traction transfer.
struct SceneFluidOpeningCapSet {
    std::uint32_t version = sceneFluidOpeningCapVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidOpeningCapSettings settings;
    std::size_t separatingBoundaryEdgeCount = 0;
    double totalAreaSquareMeters = 0.0;
    std::vector<SceneFluidOpeningCap> caps;
    std::vector<SceneFluidOpeningCapTriangle> triangles;

    bool operator==(const SceneFluidOpeningCapSet&) const = default;
};

[[nodiscard]] SceneFluidOpeningCapSet buildSceneFluidOpeningCaps(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSettings& settings = {},
    const SceneFluidOpeningCapLimits& limits = {});

void validateSceneFluidOpeningCaps(
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state);

} // namespace simwing::fsi
