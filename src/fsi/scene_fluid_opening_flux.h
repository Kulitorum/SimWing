#pragma once

#include "scene_fluid_opening_patch.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidOpeningFluxVersion = 2;

struct SceneFluidOpeningFluxLimits {
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumPatchSamples = 10'000'000;
    std::size_t maximumVelocityEvaluations = 100'000'000;
    std::size_t maximumFluxBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidOpeningRegionFlux {
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    double outwardFluidVolumeFlowRateCubicMetersPerSecond = 0.0;
    double outwardSurfaceSweepRateCubicMetersPerSecond = 0.0;
    double outwardRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningRegionFlux&) const = default;
};

struct SceneFluidOpeningFluxSample {
    std::uint64_t patchStableId = 0;
    std::uint64_t sourcePointStableId = 0;
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    double areaSquareMeters = 0.0;
    double fluidNormalVelocityMetersPerSecond = 0.0;
    double surfaceNormalVelocityMetersPerSecond = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;
    double fluidVolumeFlowRateCubicMetersPerSecond = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningFluxSample&) const = default;
};

struct SceneFluidOpeningFlux {
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    std::size_t firstSample = 0;
    std::size_t sampleCount = 0;
    double areaSquareMeters = 0.0;
    double fluidVolumeFlowRateCubicMetersPerSecond = 0.0;
    double surfaceSweepRateCubicMetersPerSecond = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidOpeningFlux&) const = default;
};

// Positive flow travels from the authored negative region to the positive
// region. Face-owned patches read the exact owning MAC normal degree of
// freedom. Cell-owned polygons use deterministic cubic triangle quadrature of
// periodic staggered trilinear velocity, then subtract accepted cap motion.
// This is an immutable diagnostic ledger; it does not alter projection,
// pressure, region connectivity, or Structure loads.
struct SceneFluidOpeningFluxSet {
    std::uint32_t version = sceneFluidOpeningFluxVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t velocityFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t velocityEvaluationCount = 0;
    std::size_t ownedStorageBytes = 0;
    double totalAreaSquareMeters = 0.0;
    double totalFluidVolumeFlowRateCubicMetersPerSecond = 0.0;
    double totalSurfaceSweepRateCubicMetersPerSecond = 0.0;
    double totalRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double globalFluidRegionBalanceResidualCubicMetersPerSecond = 0.0;
    double globalSurfaceRegionBalanceResidualCubicMetersPerSecond = 0.0;
    double globalRelativeRegionBalanceResidualCubicMetersPerSecond = 0.0;
    std::vector<SceneFluidOpeningRegionFlux> regions;
    std::vector<SceneFluidOpeningFlux> openings;
    std::vector<SceneFluidOpeningFluxSample> samples;

    bool operator==(const SceneFluidOpeningFluxSet&) const = default;
};

[[nodiscard]] SceneFluidOpeningFluxSet evaluateSceneFluidOpeningFlux(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocityMetersPerSecond,
    const SceneFluidOpeningFluxLimits& limits = {});

void validateSceneFluidOpeningFlux(
    const SceneFluidOpeningFluxSet& flux,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocityMetersPerSecond);

} // namespace simwing::fsi
