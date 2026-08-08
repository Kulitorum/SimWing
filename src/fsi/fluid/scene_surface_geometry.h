#pragma once

#include "fluid/grid.h"
#include "scene_fluid_surface.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidGridCandidateVersion = 1;

struct SceneFluidGridCandidateSettings {
    double boundingPaddingMeters = 0.0;
    double domainBoundaryToleranceMeters = 1.0e-12;

    bool operator==(const SceneFluidGridCandidateSettings&) const = default;
};

struct SceneFluidGridCandidateLimits {
    std::size_t maximumCandidates = 5'000'000;
    std::size_t maximumCandidateBytes = 256ULL * 1024ULL * 1024ULL;
};

struct GridCellCoordinate {
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;

    bool operator==(const GridCellCoordinate&) const = default;
};

// One conservative broad-phase pair. It says only that the current triangle's
// padded axis-aligned bounding box overlaps the cell; exact intersection,
// volume fractions, region labels, and face crossings remain unclassified.
struct SceneFluidCellCandidate {
    std::size_t cellIndex = 0;
    GridCellCoordinate cell;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;

    bool operator==(const SceneFluidCellCandidate&) const = default;
};

struct SceneFluidTriangleCandidateBounds {
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    GridCellCoordinate firstCell;
    GridCellCoordinate lastCell;
    std::size_t candidateCount = 0;

    bool operator==(const SceneFluidTriangleCandidateBounds&) const = default;
};

struct SceneFluidGridCandidateSet {
    std::uint32_t version = sceneFluidGridCandidateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    SceneFluidGridCandidateSettings settings;
    std::vector<SceneFluidTriangleCandidateBounds> triangleBounds;
    // Canonical cell-major, then triangle-major order.
    std::vector<SceneFluidCellCandidate> candidates;

    [[nodiscard]] std::span<const SceneFluidCellCandidate>
    candidatesForCell(std::size_t cellIndex) const noexcept;

    bool operator==(const SceneFluidGridCandidateSet&) const = default;
};

// The entire state surface must lie inside the grid bounds within the explicit
// tolerance. Periodic image selection and exact cut-cell geometry are later
// topology decisions and are not silently inferred here.
[[nodiscard]] SceneFluidGridCandidateSet buildSceneFluidGridCandidates(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSettings& settings = {},
    const SceneFluidGridCandidateLimits& limits = {});

void validateSceneFluidGridCandidates(
    const SceneFluidGridCandidateSet& candidateSet,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid);

} // namespace simwing::fsi::fluid
