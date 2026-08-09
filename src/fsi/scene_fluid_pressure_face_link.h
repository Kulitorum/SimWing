#pragma once

#include "scene_fluid_pressure_control_volume.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureFaceLinkVersion = 1;

struct SceneFluidPressureFaceLinkSettings {
    double areaToleranceSquareMeters = 1.0e-12;

    bool operator==(const SceneFluidPressureFaceLinkSettings&) const = default;
};

struct SceneFluidPressureFaceLinkLimits {
    std::size_t maximumFaces = 30'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumLinkBytes = 2048ULL * 1024ULL * 1024ULL;
};

enum class SceneFluidPressureFaceStatus : std::uint8_t {
    ResolvedFull = 1,
    ResolvedPartition = 2,
    UnresolvedActive = 3,
    UnresolvedAmbiguous = 4,
    UnresolvedOpening = 5,
};

struct SceneFluidPressureFace {
    std::size_t faceIndex = 0;
    std::uint64_t stableId = 0;
    fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t minusCellIndex = 0;
    std::size_t plusCellIndex = 0;
    SceneFluidPressureFaceStatus status =
        SceneFluidPressureFaceStatus::UnresolvedAmbiguous;
    std::size_t firstLink = 0;
    std::size_t linkCount = 0;
    double faceAreaSquareMeters = 0.0;
    double linkedAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;

    bool operator==(const SceneFluidPressureFace&) const = default;
};

struct SceneFluidPressureFaceLink {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t faceIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t regionIndex = 0;
    std::size_t componentIndex = 0;
    std::size_t minusControlVolumeIndex = 0;
    std::size_t plusControlVolumeIndex = 0;
    double areaSquareMeters = 0.0;
    double centerDistanceMeters = 0.0;
    double geometryWeightMeters = 0.0;

    bool operator==(const SceneFluidPressureFaceLink&) const = default;
};

// Conservative first pressure-face subset. Exact active-face partitions
// create one same-region link per positive region area. An untouched face is
// linked over its full area only when the two adjacent sparse cells have one
// unambiguous common region. Open chains, coplanar sheets, authored-opening
// boundaries, and ambiguous inactive faces remain explicit unresolved faces;
// no dominant-cell or smeared-interface fallback is permitted.
struct SceneFluidPressureFaceLinkSet {
    std::uint32_t version = sceneFluidPressureFaceLinkVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidPressureFaceLinkSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t resolvedFullFaceCount = 0;
    std::size_t resolvedPartitionFaceCount = 0;
    std::size_t unresolvedActiveFaceCount = 0;
    std::size_t unresolvedAmbiguousFaceCount = 0;
    std::size_t unresolvedOpeningFaceCount = 0;
    double totalFaceAreaSquareMeters = 0.0;
    double totalLinkedAreaSquareMeters = 0.0;
    double maximumResolvedAreaResidualSquareMeters = 0.0;
    std::vector<SceneFluidPressureFace> faces;
    std::vector<SceneFluidPressureFaceLink> links;

    bool operator==(const SceneFluidPressureFaceLinkSet&) const = default;
};

[[nodiscard]] SceneFluidPressureFaceLinkSet
buildSceneFluidPressureFaceLinks(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSettings& settings = {},
    const SceneFluidPressureFaceLinkLimits& limits = {});

void validateSceneFluidPressureFaceLinks(
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes);

} // namespace simwing::fsi
