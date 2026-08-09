#pragma once

#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_pressure_operator.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureEpochVersion = 3;

struct SceneFluidPressureEpochSettings {
    SceneFluidGridEpochSettings gridEpoch;
    SceneFluidCellVolumeSettings cellVolumes;
    SceneFluidOpeningGridPatchSettings openingPatches;
    SceneFluidCappedFacePartitionSettings cappedFacePartitions;
    SceneFluidPressureFaceLinkSettings faceLinks;

    bool operator==(const SceneFluidPressureEpochSettings&) const = default;
};

struct SceneFluidPressureEpochLimits {
    SceneFluidGridEpochLimits gridEpoch;
    SceneFluidCellVolumeLimits cellVolumes;
    SceneFluidOpeningQuadratureLimits openingQuadrature;
    SceneFluidOpeningGridPatchLimits openingPatches;
    SceneFluidOpeningFaceCrossingLimits openingFaceCrossings;
    SceneFluidCappedFacePartitionLimits cappedFacePartitions;
    SceneFluidPressureControlVolumeLimits pressureControlVolumes;
    SceneFluidPressureFaceLinkLimits faceLinks;
    SceneFluidPressureOperatorLimits pressureOperator;
    // Sum of vector payload storage owned by the completed pressure epoch.
    // Every constituent retains its tighter stage-specific limits.
    std::size_t maximumEpochBytes =
        8ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One immutable, fully resolved pressure geometry/operator epoch for an
// accepted Structure surface state. It composes the grid epoch, authored
// opening topology, exact opening grid patches and transverse cap crossings,
// sparse cell/region volumes, pressure unknowns, conservative face links, and
// the ungauged graph
// Laplacian. It samples no velocity, solves no pressure, and applies no load.
struct SceneFluidPressureEpoch {
    std::uint32_t version = sceneFluidPressureEpochVersion;
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
    SceneFluidPressureOperator pressureOperator;

    bool operator==(const SceneFluidPressureEpoch&) const = default;
};

[[nodiscard]] SceneFluidPressureEpoch buildSceneFluidPressureEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureEpochSettings& settings = {},
    const SceneFluidPressureEpochLimits& limits = {});

void validateSceneFluidPressureEpoch(
    const SceneFluidPressureEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity);

} // namespace simwing::fsi
