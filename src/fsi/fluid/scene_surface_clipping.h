#pragma once

#include "fluid/scene_surface_intersection.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidGridPatchVersion = 2;

enum class SceneFluidPatchDimension : std::uint8_t {
    Point = 1,
    Segment = 2,
    Area = 3,
};

enum SceneFluidCellBoundaryPlane : std::uint8_t {
    CellBoundaryNone = 0,
    CellBoundaryXMinus = 1U << 0U,
    CellBoundaryXPlus = 1U << 1U,
    CellBoundaryYMinus = 1U << 2U,
    CellBoundaryYPlus = 1U << 3U,
    CellBoundaryZMinus = 1U << 4U,
    CellBoundaryZPlus = 1U << 5U,
};

struct SceneFluidClippedVertex {
    Vec3 positionMeters;
    std::array<double, 3> barycentricCoordinates{};

    bool operator==(const SceneFluidClippedVertex& other) const {
        return positionMeters.x == other.positionMeters.x
            && positionMeters.y == other.positionMeters.y
            && positionMeters.z == other.positionMeters.z
            && barycentricCoordinates == other.barycentricCoordinates;
    }
};

// Exact convex clipping primitive shared by material-surface and virtual
// opening geometry. It carries source-triangle barycentric coordinates so a
// caller can integrate accepted piecewise-linear kinematics without inventing
// a second interpolation rule.
struct SceneFluidTriangleBoxClip {
    SceneFluidPatchDimension dimension = SceneFluidPatchDimension::Point;
    std::uint8_t coincidentBoundaryPlanes = CellBoundaryNone;
    double areaSquareMeters = 0.0;
    Vec3 centroidMeters;
    std::array<double, 3> centroidBarycentricCoordinates{};
    std::vector<SceneFluidClippedVertex> vertices;
};

[[nodiscard]] std::optional<SceneFluidTriangleBoxClip>
clipSceneFluidTriangleToAxisAlignedBox(
    const std::array<Vec3, 3>& triangle,
    const Vec3& lowerMeters,
    const Vec3& upperMeters);

struct SceneFluidCellPatch {
    std::size_t cellIndex = 0;
    GridCellCoordinate cell;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    std::size_t firstVertex = 0;
    std::size_t vertexCount = 0;
    SceneFluidPatchDimension dimension = SceneFluidPatchDimension::Point;
    std::uint8_t coincidentBoundaryPlanes = CellBoundaryNone;
    double areaSquareMeters = 0.0;
    Vec3 centroidMeters;
    std::array<double, 3> centroidBarycentricCoordinates{};

    bool operator==(const SceneFluidCellPatch& other) const {
        return cellIndex == other.cellIndex
            && cell == other.cell
            && triangleIndex == other.triangleIndex
            && triangleId == other.triangleId
            && firstVertex == other.firstVertex
            && vertexCount == other.vertexCount
            && dimension == other.dimension
            && coincidentBoundaryPlanes == other.coincidentBoundaryPlanes
            && areaSquareMeters == other.areaSquareMeters
            && centroidMeters.x == other.centroidMeters.x
            && centroidMeters.y == other.centroidMeters.y
            && centroidMeters.z == other.centroidMeters.z
            && centroidBarycentricCoordinates
                == other.centroidBarycentricCoordinates;
    }
};

struct SceneFluidGridPatchLimits {
    std::size_t maximumPatches = 5'000'000;
    std::size_t maximumVertices = 45'000'000;
    std::size_t maximumPatchBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidGridPatchSet {
    std::uint32_t version = sceneFluidGridPatchVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t intersectionSetFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<SceneFluidCellPatch> patches;
    std::vector<SceneFluidClippedVertex> vertices;

    [[nodiscard]] std::span<const SceneFluidCellPatch>
    patchesForCell(std::size_t cellIndex) const noexcept;
    [[nodiscard]] std::span<const SceneFluidClippedVertex>
    verticesForPatch(const SceneFluidCellPatch& patch) const;

    bool operator==(const SceneFluidGridPatchSet&) const = default;
};

// Clips every exact triangle/cell intersection against the six cell planes.
// Point and segment contacts remain explicit with zero area. Area patches that
// are coincident with a shared cell plane remain duplicated and flagged; this
// stage deliberately does not choose the unique fluid-face owner required for
// pressure integration.
[[nodiscard]] SceneFluidGridPatchSet clipSceneFluidSurfaceToCells(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchLimits& limits = {});

void validateSceneFluidGridPatches(
    const SceneFluidGridPatchSet& patches,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections);

} // namespace simwing::fsi::fluid
