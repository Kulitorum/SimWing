#pragma once

#include "scene_fluid_grid_epoch.h"
#include "scene_fluid_opening_face_crossing.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidCappedFacePartitionVersion = 4;
inline constexpr std::size_t invalidSceneFluidActiveFaceIndex =
    std::numeric_limits<std::size_t>::max();
inline constexpr std::size_t invalidSceneFluidCappedFacePartitionIndex =
    std::numeric_limits<std::size_t>::max();

struct SceneFluidCappedFacePartitionSettings {
    double geometryToleranceMeters = 1.0e-12;
    double areaClosureToleranceSquareMeters = 1.0e-12;

    bool operator==(
        const SceneFluidCappedFacePartitionSettings&) const = default;
};

struct SceneFluidCappedFacePartitionLimits {
    std::size_t maximumTouchedFaces = 10'000'000;
    std::size_t maximumPartitions = 10'000'000;
    std::size_t maximumReferences = 30'000'000;
    std::size_t maximumSegmentPairTests = 100'000'000;
    std::size_t maximumPartitionBytes = 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidCappedFacePartition {
    std::uint64_t stableId = 0;
    fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t activeFaceIndex = invalidSceneFluidActiveFaceIndex;
    std::size_t firstMaterialChainReference = 0;
    std::size_t materialChainReferenceCount = 0;
    std::size_t firstOpeningCrossingReference = 0;
    std::size_t openingCrossingReferenceCount = 0;
    std::size_t firstRegionArea = 0;
    std::size_t regionAreaCount = 0;
    double faceAreaSquareMeters = 0.0;
    double assignedAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;

    bool operator==(const SceneFluidCappedFacePartition&) const = default;
};

enum class SceneFluidCappedFaceStatus : std::uint8_t {
    Resolved = 0,
    FaceOwnedOpeningArea = 1,
    CoplanarMaterial = 2,
    InvalidSourceGeometry = 3,
    UnpairedMaterialEndpoint = 4,
    UnpairedOpeningEndpoint = 5,
    UnstitchedIntersection = 6,
    ConflictingWinding = 7,
    AmbiguousRegionOwnership = 8,
    AreaClosureFailure = 9,
    InvalidArrangementTopology = 10,
};

struct SceneFluidCappedFace {
    fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    SceneFluidCappedFaceStatus status =
        SceneFluidCappedFaceStatus::InvalidSourceGeometry;
    // The unique material-chain or opening-crossing stable ID responsible for
    // a classified failure, when one source can be named without guessing.
    StableId failureSourceStableId = invalidStableId;
    std::size_t partitionIndex =
        invalidSceneFluidCappedFacePartitionIndex;

    bool operator==(const SceneFluidCappedFace&) const = default;
};

// Exact same-region areas on Cartesian faces crossed transversely by one or
// more virtual opening-cap triangles. Directed material-chain segments and cap
// segments share one bounded planar arrangement. Face-owned aperture patches
// remain excluded because pressure links already own their cross-region area.
// `faces` contains exactly one deterministic record per touched face. Status
// records why unresolved geometry was rejected; only Resolved has a valid
// partition index. Unsupported coplanar material, unpaired endpoints,
// unstitched crossings, or conflicting authored winding remain explicit. Each
// resolved region area retains its exact global first moment and face centroid.
struct SceneFluidCappedFacePartitionSet {
    std::uint32_t version = sceneFluidCappedFacePartitionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t openingFaceCrossingFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidCappedFacePartitionSettings settings;
    std::size_t touchedFaceCount = 0;
    std::size_t unresolvedTouchedFaceCount = 0;
    std::size_t segmentPairTestCount = 0;
    std::size_t ownedStorageBytes = 0;
    std::vector<SceneFluidCappedFace> faces;
    std::vector<SceneFluidCappedFacePartition> partitions;
    std::vector<std::size_t> materialChainReferences;
    std::vector<std::size_t> openingCrossingReferences;
    std::vector<fluid::SceneFluidFaceRegionArea> regionAreas;

    bool operator==(
        const SceneFluidCappedFacePartitionSet&) const = default;
};

[[nodiscard]] SceneFluidCappedFacePartitionSet
buildSceneFluidCappedFacePartitions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings,
    const SceneFluidCappedFacePartitionSettings& settings = {},
    const SceneFluidCappedFacePartitionLimits& limits = {});

void validateSceneFluidCappedFacePartitions(
    const SceneFluidCappedFacePartitionSet& partitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings);

} // namespace simwing::fsi
