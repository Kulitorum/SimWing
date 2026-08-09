#pragma once

#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_pressure_control_volume.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureFaceLinkVersion = 8;
inline constexpr std::size_t invalidSceneFluidPressureFaceIndex =
    std::numeric_limits<std::size_t>::max();

struct SceneFluidPressureFaceLinkSettings {
    double areaToleranceSquareMeters = 1.0e-12;
    double minimumCenterDistanceMeters = 1.0e-10;

    bool operator==(const SceneFluidPressureFaceLinkSettings&) const = default;
};

struct SceneFluidPressureFaceLinkLimits {
    std::size_t maximumFaces = 30'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumEmbeddedOpeningRejections = 10'000'000;
    std::size_t maximumEmbeddedOpeningOneRingSupports = 120'000'000;
    std::size_t maximumLinkBytes = 2048ULL * 1024ULL * 1024ULL;
};

enum class SceneFluidPressureFaceStatus : std::uint8_t {
    ResolvedFull = 1,
    ResolvedPartition = 2,
    UnresolvedActive = 3,
    UnresolvedAmbiguous = 4,
    UnresolvedOpening = 5,
    ResolvedOpening = 6,
    UnresolvedCapped = 7,
};

enum class SceneFluidPressureFaceLinkKind : std::uint8_t {
    SameRegion = 1,
    AuthoredOpening = 2,
};

enum class SceneFluidPressureLinkGeometryKind : std::uint8_t {
    CartesianFace = 1,
    EmbeddedOpening = 2,
};

enum class SceneFluidEmbeddedOpeningRejectionStatus : std::uint8_t {
    NonPositiveProjectedDistance = 1,
    BelowMinimumProjectedDistance = 2,
};

enum class SceneFluidEmbeddedOpeningOneRingSide : std::uint8_t {
    NegativeRegion = 1,
    PositiveRegion = 2,
};

enum class SceneFluidEmbeddedOpeningOneRingStatus : std::uint8_t {
    NeitherSide = 1,
    NegativeSideOnly = 2,
    PositiveSideOnly = 3,
    BothSides = 4,
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
    std::size_t faceIndex = invalidSceneFluidPressureFaceIndex;
    SceneFluidPressureFaceLinkKind kind =
        SceneFluidPressureFaceLinkKind::SameRegion;
    SceneFluidPressureLinkGeometryKind geometryKind =
        SceneFluidPressureLinkGeometryKind::CartesianFace;
    StableId minusRegionId = invalidStableId;
    StableId plusRegionId = invalidStableId;
    std::size_t minusRegionIndex = 0;
    std::size_t plusRegionIndex = 0;
    std::size_t componentIndex = 0;
    std::size_t minusControlVolumeIndex = 0;
    std::size_t plusControlVolumeIndex = 0;
    StableId openingId = invalidStableId;
    std::uint64_t openingPatchStableId = 0;
    double areaSquareMeters = 0.0;
    fluid::Vector3 faceCentroidMeters;
    double centerDistanceMeters = 0.0;
    double geometryWeightMeters = 0.0;
    fluid::Vector3 unitNormalMinusToPlus;

    bool operator==(const SceneFluidPressureFaceLink&) const = default;
};

struct SceneFluidEmbeddedOpeningRejection {
    std::size_t rejectionIndex = 0;
    std::uint64_t openingPatchStableId = 0;
    StableId openingId = invalidStableId;
    std::size_t cellIndex = 0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    std::size_t negativeControlVolumeIndex = 0;
    std::size_t positiveControlVolumeIndex = 0;
    SceneFluidEmbeddedOpeningRejectionStatus status =
        SceneFluidEmbeddedOpeningRejectionStatus::
            NonPositiveProjectedDistance;
    SceneFluidEmbeddedOpeningOneRingStatus oneRingStatus =
        SceneFluidEmbeddedOpeningOneRingStatus::NeitherSide;
    std::size_t firstOneRingSupport = 0;
    std::size_t oneRingSupportCount = 0;
    std::size_t negativeOneRingSupportCount = 0;
    std::size_t positiveOneRingSupportCount = 0;
    std::size_t negativeAdmissibleOneRingSupportCount = 0;
    std::size_t positiveAdmissibleOneRingSupportCount = 0;
    double areaSquareMeters = 0.0;
    double projectedCenterDistanceMeters = 0.0;
    double negativeCentroidSignedDistanceMeters = 0.0;
    double positiveCentroidSignedDistanceMeters = 0.0;

    bool operator==(
        const SceneFluidEmbeddedOpeningRejection&) const = default;
};

struct SceneFluidEmbeddedOpeningOneRingSupport {
    std::size_t supportIndex = 0;
    std::size_t rejectionIndex = 0;
    SceneFluidEmbeddedOpeningOneRingSide side =
        SceneFluidEmbeddedOpeningOneRingSide::NegativeRegion;
    std::size_t cartesianFaceLinkIndex = 0;
    std::uint64_t cartesianFaceLinkStableId = 0;
    std::size_t rootControlVolumeIndex = 0;
    std::size_t donorControlVolumeIndex = 0;
    fluid::Vector3 donorOffsetFromOpeningCentroidMeters;
    double donorProjectedDistanceMeters = 0.0;
    bool isCorrectlySided = false;

    bool operator==(
        const SceneFluidEmbeddedOpeningOneRingSupport&) const = default;
};

// Conservative pressure-face subset. Exact active-face partitions create one
// same-region link per positive region area, including the exact subface
// centroid. On a face crossed transversely by
// a virtual opening cap, the capped material-plus-opening partition supersedes
// that material-only result. An untouched face is linked over its full area
// only when the two adjacent sparse cells have one unambiguous common region.
// A face-aligned authored opening instead contributes oriented cross-region
// links over its exact patch area and, when unambiguous, one complementary
// same-region link. A cell-owned opening patch instead connects its two
// same-cell pressure controls along the authored normal, using their projected
// centroid separation. A patch without admissible positive projected
// separation remains an explicit typed rejection with source identity and
// signed centroid-to-patch distances, and publishes no link. Such a rejection
// also retains every same-region Cartesian one-ring neighbor of both side
// controls, with periodic-image geometry and explicit sidedness. This is
// bounded reconstruction evidence only; it does not reroute aperture flux or
// claim that a conservative symmetric multipoint stencil exists.
// Material/open-chain/coplanar overlap, unresolved capped arrangements, and
// ambiguous faces likewise remain explicit and unresolved; no dominant-cell or
// smeared-interface fallback is permitted.
struct SceneFluidPressureFaceLinkSet {
    std::uint32_t version = sceneFluidPressureFaceLinkVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t cappedFacePartitionFingerprint = 0;
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
    std::size_t resolvedOpeningFaceCount = 0;
    std::size_t embeddedOpeningLinkCount = 0;
    std::size_t unresolvedEmbeddedOpeningPatchCount = 0;
    std::size_t embeddedOpeningBothSideOneRingCount = 0;
    std::size_t embeddedOpeningSingleSideOneRingCount = 0;
    std::size_t embeddedOpeningNoSideOneRingCount = 0;
    std::size_t unresolvedActiveFaceCount = 0;
    std::size_t unresolvedCappedFaceCount = 0;
    std::size_t unresolvedAmbiguousFaceCount = 0;
    std::size_t unresolvedOpeningFaceCount = 0;
    double totalFaceAreaSquareMeters = 0.0;
    double totalLinkedAreaSquareMeters = 0.0;
    double totalEmbeddedOpeningAreaSquareMeters = 0.0;
    double unresolvedEmbeddedOpeningAreaSquareMeters = 0.0;
    double maximumResolvedAreaResidualSquareMeters = 0.0;
    std::vector<SceneFluidPressureFace> faces;
    std::vector<SceneFluidPressureFaceLink> links;
    std::vector<SceneFluidEmbeddedOpeningRejection>
        embeddedOpeningRejections;
    std::vector<SceneFluidEmbeddedOpeningOneRingSupport>
        embeddedOpeningOneRingSupports;

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
    const SceneFluidOpeningFaceCrossingSet& openingFaceCrossings,
    const SceneFluidCappedFacePartitionSet& cappedFacePartitions,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSettings& settings = {},
    const SceneFluidPressureFaceLinkLimits& limits = {});

// Lightweight immutable-product check for downstream current-topology
// adapters that already received an accepted face-link set and must reject
// accidental nested mutation without rebuilding the complete geometry chain.
void validateSceneFluidPressureFaceLinkIntegrity(
    const SceneFluidPressureFaceLinkSet& faceLinks);

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
    const SceneFluidOpeningFaceCrossingSet& openingFaceCrossings,
    const SceneFluidCappedFacePartitionSet& cappedFacePartitions,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes);

} // namespace simwing::fsi
