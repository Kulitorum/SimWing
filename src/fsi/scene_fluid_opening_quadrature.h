#pragma once

#include "scene_fluid_opening_cap.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidOpeningQuadratureVersion = 1;

struct SceneFluidOpeningQuadratureLimits {
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumPoints = 10'000'000;
    std::size_t maximumQuadratureBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidOpeningQuadraturePoint {
    std::uint64_t stableId = 0;
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    std::size_t triangleOrdinal = 0;
    std::array<std::size_t, 3> vertexIndices{};
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    Vec3 positionMeters;
    Vec3 velocityMetersPerSecond;
    Vec3 unitNormalNegativeToPositive;
    double areaSquareMeters = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningQuadraturePoint& other) const;
};

struct SceneFluidOpeningQuadrature {
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    std::size_t firstPoint = 0;
    std::size_t pointCount = 0;
    Vec3 unitNormalNegativeToPositive;
    Vec3 centroidMeters;
    double areaSquareMeters = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningQuadrature& other) const;
};

// One-point triangle quadrature exactly integrates the piecewise-linear cap
// velocity used by the geometric conservation law. These samples describe a
// virtual fluid boundary only: they are not Structure nodes, material patches,
// or aerodynamic traction targets. Grid-resolved fluid velocity and mass flux
// remain separate later operators.
struct SceneFluidOpeningQuadratureSet {
    std::uint32_t version = sceneFluidOpeningQuadratureVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t openingCapFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    double totalAreaSquareMeters = 0.0;
    double totalSurfaceSweepRateCubicMetersPerSecond = 0.0;
    std::vector<SceneFluidOpeningQuadrature> openings;
    std::vector<SceneFluidOpeningQuadraturePoint> points;

    bool operator==(const SceneFluidOpeningQuadratureSet&) const = default;
};

[[nodiscard]] SceneFluidOpeningQuadratureSet
buildSceneFluidOpeningQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureLimits& limits = {});

void validateSceneFluidOpeningQuadrature(
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps);

} // namespace simwing::fsi
