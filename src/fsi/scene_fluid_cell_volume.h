#pragma once

#include "scene_fluid_grid_epoch.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidCellVolumeVersion = 1;

struct SceneFluidCellVolumeSettings {
    double absoluteVolumeToleranceCubicMeters = 1.0e-12;
    double relativeVolumeTolerance = 1.0e-10;

    bool operator==(const SceneFluidCellVolumeSettings&) const = default;
};

struct SceneFluidCellVolumeLimits {
    std::size_t maximumCells = 10'000'000;
    std::size_t maximumContributionEvents = 100'000'000;
    std::size_t maximumCellRegionVolumes = 50'000'000;
    std::size_t maximumVolumeBytes = 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidCellRegionVolume {
    StableId regionId = invalidStableId;
    double volumeCubicMeters = 0.0;
    double volumeFraction = 0.0;

    bool operator==(const SceneFluidCellRegionVolume&) const = default;
};

struct SceneFluidCellVolume {
    std::size_t cellIndex = 0;
    fluid::GridCellCoordinate cell;
    std::size_t firstRegionVolume = 0;
    std::size_t regionVolumeCount = 0;
    double assignedVolumeCubicMeters = 0.0;
    double volumeResidualCubicMeters = 0.0;

    bool operator==(const SceneFluidCellVolume&) const = default;
};

struct SceneFluidRegionVolume {
    StableId regionId = invalidStableId;
    double summedCellVolumeCubicMeters = 0.0;
    double wholeSurfaceVolumeCubicMeters = 0.0;
    double volumeResidualCubicMeters = 0.0;

    bool operator==(const SceneFluidRegionVolume&) const = default;
};

// First bounded cut-cell volume subset. All interface components must be
// closed, consistently wound two-sided triangle manifolds; the scene must have
// one Outside root, no authored opening, no coplanar face-owned area, and no
// unresolved active MAC face. Sparse positive region volumes are published for
// every cell. Whole-surface divergence volume must equal their sum, preventing
// an enclosed full cell from being silently labelled as Outside. General open
// intakes, junctions, boundary contact, and unresolved face topology reject.
struct SceneFluidCellVolumeSet {
    std::uint32_t version = sceneFluidCellVolumeVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidCellVolumeSettings settings;
    StableId outsideRegionId = invalidStableId;
    double cellVolumeCubicMeters = 0.0;
    double maximumCellVolumeResidualCubicMeters = 0.0;
    double maximumRegionVolumeResidualCubicMeters = 0.0;
    std::vector<SceneFluidCellVolume> cells;
    std::vector<SceneFluidCellRegionVolume> cellRegionVolumes;
    std::vector<SceneFluidRegionVolume> regionVolumes;

    bool operator==(const SceneFluidCellVolumeSet&) const = default;
};

[[nodiscard]] SceneFluidCellVolumeSet buildSceneFluidCellVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidCellVolumeSettings& settings = {},
    const SceneFluidCellVolumeLimits& limits = {});

void validateSceneFluidCellVolumes(
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch);

} // namespace simwing::fsi
