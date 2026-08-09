#pragma once

#include "scene_fluid_cell_volume.h"
#include "scene_fluid_opening_flux.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionContinuityVersion = 1;

struct SceneFluidRegionContinuitySettings {
    double absoluteVolumeToleranceCubicMeters = 1.0e-11;
    double relativeVolumeTolerance = 1.0e-9;

    bool operator==(
        const SceneFluidRegionContinuitySettings&) const = default;
};

struct SceneFluidRegionContinuityLimits {
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumContinuityBytes = 256ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionContinuity {
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double previousOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double currentOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double integratedOutwardRelativeVolumeCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double toleranceCubicMeters = 0.0;
    bool withinTolerance = false;

    bool operator==(const SceneFluidRegionContinuity&) const = default;
};

// Two consecutive accepted geometry/flow epochs form one immutable
// macro-step audit. For each authored region, incompressible moving-control-
// volume continuity is
//
//   volume change + integral(outward relative opening flow) = 0.
//
// Endpoint opening flow is trapezoidally integrated. This adapter diagnoses
// compatibility only; it does not modify projection, pressure, connectivity,
// or Structure state.
struct SceneFluidRegionContinuitySet {
    std::uint32_t version = sceneFluidRegionContinuityVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousSurfaceStateFingerprint = 0;
    std::uint64_t currentSurfaceStateFingerprint = 0;
    std::uint64_t previousCellVolumeFingerprint = 0;
    std::uint64_t currentCellVolumeFingerprint = 0;
    std::uint64_t previousOpeningFluxFingerprint = 0;
    std::uint64_t currentOpeningFluxFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    double durationSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidRegionContinuitySettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t failedRegionCount = 0;
    double maximumAbsoluteContinuityResidualCubicMeters = 0.0;
    double globalGeometryVolumeChangeCubicMeters = 0.0;
    double globalIntegratedOutwardRelativeVolumeCubicMeters = 0.0;
    double globalContinuityResidualCubicMeters = 0.0;
    bool allRegionsWithinTolerance = false;
    std::vector<SceneFluidRegionContinuity> regions;

    bool operator==(
        const SceneFluidRegionContinuitySet&) const = default;
};

[[nodiscard]] SceneFluidRegionContinuitySet
auditSceneFluidRegionContinuity(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionContinuitySettings& settings = {},
    const SceneFluidRegionContinuityLimits& limits = {});

void validateSceneFluidRegionContinuity(
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux);

} // namespace simwing::fsi
