#pragma once

#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_pressure_face_link.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticGeometryEpochVersion = 1;

struct SceneFluidMimeticGeometryEpochSettings {
    SceneFluidGridEpochSettings gridEpoch;
    SceneFluidCellVolumeSettings cellVolumes;
    SceneFluidOpeningGridPatchSettings openingPatches;
    SceneFluidCappedFacePartitionSettings cappedFacePartitions;
    SceneFluidPressureFaceLinkSettings faceLinks;

    bool operator==(
        const SceneFluidMimeticGeometryEpochSettings&) const = default;
};

struct SceneFluidMimeticGeometryEpochLimits {
    SceneFluidGridEpochLimits gridEpoch;
    SceneFluidCellVolumeLimits cellVolumes;
    SceneFluidOpeningQuadratureLimits openingQuadrature;
    SceneFluidOpeningGridPatchLimits openingPatches;
    SceneFluidOpeningFaceCrossingLimits openingFaceCrossings;
    SceneFluidCappedFacePartitionLimits cappedFacePartitions;
    SceneFluidPressureControlVolumeLimits pressureControlVolumes;
    SceneFluidPressureFaceLinkLimits faceLinks;
    std::size_t maximumEpochBytes =
        8ULL * 1024ULL * 1024ULL * 1024ULL;
};

// Complete graph-free geometry/topology epoch for the mixed-hybrid pressure
// path. It deliberately stops before building the older graph Laplacian, so
// authored intake traces that are valid mixed-hybrid faces do not need to be
// representable by that reference operator. Every constituent remains bound
// to one accepted Structure surface state and one Cartesian grid.
struct SceneFluidMimeticGeometryEpoch {
    std::uint32_t version = sceneFluidMimeticGeometryEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t regionConnectivityFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    SceneFluidGridEpoch gridEpoch;
    SceneFluidOpeningCapSet openingCaps;
    SceneFluidOpeningQuadratureSet openingQuadrature;
    SceneFluidOpeningGridPatchSet openingPatches;
    SceneFluidOpeningFaceCrossingSet openingFaceCrossings;
    SceneFluidCappedFacePartitionSet cappedFacePartitions;
    SceneFluidCellVolumeSet cellVolumes;
    SceneFluidPressureControlVolumeSet pressureControlVolumes;
    SceneFluidPressureFaceLinkSet pressureFaceLinks;

    bool operator==(const SceneFluidMimeticGeometryEpoch&) const = default;
};

[[nodiscard]] SceneFluidMimeticGeometryEpoch
buildSceneFluidMimeticGeometryEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidMimeticGeometryEpochSettings& settings = {},
    const SceneFluidMimeticGeometryEpochLimits& limits = {});

void validateSceneFluidMimeticGeometryEpochIntegrity(
    const SceneFluidMimeticGeometryEpoch& epoch);

void validateSceneFluidMimeticGeometryEpoch(
    const SceneFluidMimeticGeometryEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity);

} // namespace simwing::fsi
