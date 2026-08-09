#pragma once

#include "scene_fluid_region_connectivity.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureControlVolumeVersion = 2;

struct SceneFluidPressureControlVolumeLimits {
    std::size_t maximumCells = 10'000'000;
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumControlVolumeBytes =
        1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureCell {
    std::size_t cellIndex = 0;
    fluid::GridCellCoordinate cell;
    std::size_t firstControlVolume = 0;
    std::size_t controlVolumeCount = 0;
    double assignedVolumeCubicMeters = 0.0;
    double volumeResidualCubicMeters = 0.0;

    bool operator==(const SceneFluidPressureCell&) const = default;
};

struct SceneFluidPressureControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    std::size_t componentIndex = 0;
    double volumeCubicMeters = 0.0;
    double volumeFraction = 0.0;
    Vec3 centroidMeters;
    bool belongsToGaugeRegion = false;

    bool operator==(
        const SceneFluidPressureControlVolume& other) const {
        return controlVolumeIndex == other.controlVolumeIndex
            && stableId == other.stableId
            && cellIndex == other.cellIndex
            && regionIndex == other.regionIndex
            && regionId == other.regionId
            && kind == other.kind
            && componentIndex == other.componentIndex
            && volumeCubicMeters == other.volumeCubicMeters
            && volumeFraction == other.volumeFraction
            && centroidMeters.x == other.centroidMeters.x
            && centroidMeters.y == other.centroidMeters.y
            && centroidMeters.z == other.centroidMeters.z
            && belongsToGaugeRegion == other.belongsToGaugeRegion;
    }
};

struct SceneFluidPressureRegion {
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    std::size_t componentIndex = 0;
    std::size_t firstControlVolumeMember = 0;
    std::size_t controlVolumeCount = 0;
    double summedControlVolumeCubicMeters = 0.0;
    double sourceVolumeResidualCubicMeters = 0.0;

    bool operator==(const SceneFluidPressureRegion&) const = default;
};

struct SceneFluidPressureComponent {
    std::size_t componentIndex = 0;
    StableId gaugeRegionId = invalidStableId;
    std::size_t gaugeControlVolumeIndex = 0;
    std::size_t firstControlVolumeMember = 0;
    std::size_t controlVolumeCount = 0;
    double totalVolumeCubicMeters = 0.0;

    bool operator==(const SceneFluidPressureComponent&) const = default;
};

// Each positive sparse cell-region volume becomes one immutable pressure
// unknown. Stable IDs derive from the Cartesian cell index and authored region
// ID, so accepted motion preserves identity while that pair remains occupied.
// Each pressure point preserves its exact source cell-region centroid.
// Components and gauge regions come only from authored opening connectivity.
// This stage owns no face conductance and performs no pressure solve.
struct SceneFluidPressureControlVolumeSet {
    std::uint32_t version = sceneFluidPressureControlVolumeVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t cellVolumeFingerprint = 0;
    std::uint64_t regionConnectivityFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    double cellVolumeCubicMeters = 0.0;
    double domainVolumeCubicMeters = 0.0;
    double maximumCellVolumeResidualCubicMeters = 0.0;
    double maximumRegionVolumeResidualCubicMeters = 0.0;
    double domainVolumeResidualCubicMeters = 0.0;
    std::vector<SceneFluidPressureCell> cells;
    std::vector<SceneFluidPressureControlVolume> controlVolumes;
    std::vector<SceneFluidPressureRegion> regions;
    std::vector<SceneFluidPressureComponent> components;
    std::vector<std::size_t> regionControlVolumeIndices;
    std::vector<std::size_t> componentControlVolumeIndices;

    bool operator==(
        const SceneFluidPressureControlVolumeSet&) const = default;
};

[[nodiscard]] SceneFluidPressureControlVolumeSet
buildSceneFluidPressureControlVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeLimits& limits = {});

// Lightweight immutable-product check for consecutive-epoch adapters that
// already received an accepted pressure-volume topology and need to reject
// accidental mutation without rebuilding the complete geometry chain.
void validateSceneFluidPressureControlVolumeIntegrity(
    const SceneFluidPressureControlVolumeSet& pressureVolumes);

void validateSceneFluidPressureControlVolumes(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity);

} // namespace simwing::fsi
