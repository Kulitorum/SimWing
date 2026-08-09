#pragma once

#include "fluid/scene_surface_clipping.h"
#include "scene_fluid_opening_quadrature.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidOpeningGridPatchVersion = 2;

enum class SceneFluidOpeningPatchOwnerKind : std::uint8_t {
    Cell = 1,
    Face = 2,
};

enum class SceneFluidOpeningPatchFaceAxis : std::uint8_t {
    X = 0,
    Y = 1,
    Z = 2,
};

struct SceneFluidOpeningGridPatchSettings {
    double boundingPaddingMeters = 0.0;
    double domainBoundaryToleranceMeters = 1.0e-12;
    double absoluteAreaToleranceSquareMeters = 1.0e-12;
    double relativeAreaTolerance = 1.0e-10;
    double absoluteSweepRateToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeSweepRateTolerance = 1.0e-10;

    bool operator==(const SceneFluidOpeningGridPatchSettings&) const = default;
};

struct SceneFluidOpeningGridPatchLimits {
    std::size_t maximumCandidateCells = 20'000'000;
    std::size_t maximumPatches = 10'000'000;
    std::size_t maximumVertices = 90'000'000;
    std::size_t maximumPatchBytes = 1ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidOpeningGridPatch {
    std::uint64_t stableId = 0;
    std::uint64_t sourcePointStableId = 0;
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    std::size_t triangleOrdinal = 0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    SceneFluidOpeningPatchOwnerKind ownerKind =
        SceneFluidOpeningPatchOwnerKind::Cell;
    std::size_t cellIndex = 0;
    fluid::GridCellCoordinate cell;
    SceneFluidOpeningPatchFaceAxis faceAxis =
        SceneFluidOpeningPatchFaceAxis::X;
    std::size_t faceI = 0;
    std::size_t faceJ = 0;
    std::size_t faceK = 0;
    std::size_t firstVertex = 0;
    std::size_t vertexCount = 0;
    Vec3 centroidMeters;
    std::array<double, 3> centroidBarycentricCoordinates{};
    Vec3 velocityMetersPerSecond;
    Vec3 unitNormalNegativeToPositive;
    double areaSquareMeters = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningGridPatch& other) const;
};

struct SceneFluidOpeningGridPointRange {
    std::uint64_t sourcePointStableId = 0;
    std::size_t firstPatch = 0;
    std::size_t patchCount = 0;
    double areaSquareMeters = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningGridPointRange&) const = default;
};

// Exact cap-triangle partitions own positive area once, either inside one cell
// or on one non-periodic grid face. Point/segment contacts carry no opening
// flux and are omitted. The result exposes geometry and accepted cap velocity,
// but deliberately does not interpolate an Eulerian velocity or claim a mass
// flow rate.
struct SceneFluidOpeningGridPatchSet {
    std::uint32_t version = sceneFluidOpeningGridPatchVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t openingCapFingerprint = 0;
    std::uint64_t openingQuadratureFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidOpeningGridPatchSettings settings;
    std::size_t candidateCellCount = 0;
    std::size_t ownedStorageBytes = 0;
    double totalAreaSquareMeters = 0.0;
    double totalSurfaceSweepRateCubicMetersPerSecond = 0.0;
    std::vector<SceneFluidOpeningGridPointRange> pointRanges;
    std::vector<SceneFluidOpeningGridPatch> patches;
    std::vector<fluid::SceneFluidClippedVertex> vertices;

    [[nodiscard]] std::span<const SceneFluidOpeningGridPatch>
    patchesForPoint(const SceneFluidOpeningGridPointRange& range) const;
    [[nodiscard]] std::span<const fluid::SceneFluidClippedVertex>
    verticesForPatch(const SceneFluidOpeningGridPatch& patch) const;

    bool operator==(const SceneFluidOpeningGridPatchSet&) const = default;
};

[[nodiscard]] SceneFluidOpeningGridPatchSet
buildSceneFluidOpeningGridPatches(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningGridPatchSettings& settings = {},
    const SceneFluidOpeningGridPatchLimits& limits = {});

void validateSceneFluidOpeningGridPatches(
    const SceneFluidOpeningGridPatchSet& patches,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const fluid::PeriodicCartesianGrid& grid);

} // namespace simwing::fsi
