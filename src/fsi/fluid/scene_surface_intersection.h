#pragma once

#include "fluid/scene_surface_geometry.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidGridIntersectionVersion = 1;

struct SceneFluidGridIntersectionSettings {
    // Contact within this Euclidean separation is retained. The broad-phase
    // padding must be at least this large so the narrow phase cannot invent a
    // pair that was never considered.
    double separationToleranceMeters = 0.0;

    bool operator==(const SceneFluidGridIntersectionSettings&) const = default;
};

struct SceneFluidGridIntersectionLimits {
    std::size_t maximumIntersections = 5'000'000;
    std::size_t maximumIntersectionBytes = 256ULL * 1024ULL * 1024ULL;
};

struct SceneFluidCellIntersection {
    std::size_t cellIndex = 0;
    GridCellCoordinate cell;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;

    bool operator==(const SceneFluidCellIntersection&) const = default;
};

struct SceneFluidGridIntersectionSet {
    std::uint32_t version = sceneFluidGridIntersectionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t candidateSetFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidGridIntersectionSettings settings;
    std::vector<std::size_t> triangleIntersectionCounts;
    // Canonical cell-major, then triangle-major order inherited from the
    // validated candidate set.
    std::vector<SceneFluidCellIntersection> intersections;

    [[nodiscard]] std::span<const SceneFluidCellIntersection>
    intersectionsForCell(std::size_t cellIndex) const noexcept;

    bool operator==(const SceneFluidGridIntersectionSet&) const = default;
};

// Applies a normalized separating-axis triangle/AABB test to every validated
// broad-phase pair. Touching counts as intersection. The result still does not
// assign regions, construct clipped polygons, or estimate cut-cell volume.
[[nodiscard]] SceneFluidGridIntersectionSet
intersectSceneFluidSurfaceWithGrid(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSettings& settings = {},
    const SceneFluidGridIntersectionLimits& limits = {});

void validateSceneFluidGridIntersections(
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates);

} // namespace simwing::fsi::fluid
