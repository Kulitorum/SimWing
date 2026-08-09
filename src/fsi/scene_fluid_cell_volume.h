#pragma once

#include "scene_fluid_grid_epoch.h"
#include "scene_fluid_opening_cap.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidCellVolumeVersion = 3;

struct SceneFluidCellVolumeSettings {
    double absoluteVolumeToleranceCubicMeters = 1.0e-12;
    double relativeVolumeTolerance = 1.0e-10;
    SceneFluidOpeningCapSettings openingCaps;

    bool operator==(const SceneFluidCellVolumeSettings&) const = default;
};

struct SceneFluidCellVolumeLimits {
    SceneFluidOpeningCapLimits openingCaps;
    std::size_t maximumCells = 10'000'000;
    std::size_t maximumContributionEvents = 100'000'000;
    std::size_t maximumTetrahedronCellClips = 100'000'000;
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

// First bounded cut-cell volume subset. Separating fabric plus optional planar,
// strictly convex authored opening caps must form consistently wound two-sided
// triangle manifolds around one Outside root. Sparse positive region volumes
// are published for every cell. Each oriented material/cap triangle defines
// one signed tetrahedron against the grid origin; exact convex clipping
// distributes that chain into intersected cells, including cells wholly inside
// a region. Caps remain topology only and never enter Structure or traction.
struct SceneFluidCellVolumeSet {
    std::uint32_t version = sceneFluidCellVolumeVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t openingCapFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidCellVolumeSettings settings;
    StableId outsideRegionId = invalidStableId;
    double cellVolumeCubicMeters = 0.0;
    std::size_t openingCapCount = 0;
    double openingCapAreaSquareMeters = 0.0;
    std::size_t tetrahedronCellClipCount = 0;
    std::size_t nonzeroTetrahedronCellClipCount = 0;
    double maximumTetrahedronVolumeResidualCubicMeters = 0.0;
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
