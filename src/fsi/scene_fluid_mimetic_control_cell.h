#pragma once

#include "scene_fluid_pressure_face_link.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticControlCellVersion = 1;
inline constexpr std::size_t invalidSceneFluidMimeticControlVolumeIndex =
    std::numeric_limits<std::size_t>::max();

enum class SceneFluidMimeticHalfFaceKind : std::uint8_t {
    CartesianTrace = 1,
    MaterialWall = 2,
    AuthoredOpeningTrace = 3,
};

enum class SceneFluidMimeticHalfFaceSide : std::uint8_t {
    MinusOrNegative = 1,
    PlusOrPositive = 2,
};

struct SceneFluidMimeticControlCellSettings {
    double absoluteAreaClosureToleranceSquareMeters = 1.0e-12;
    double absoluteDivergenceTheoremToleranceCubicMeters = 1.0e-12;
    double relativeGeometryTolerance = 1.0e-10;
    double unitNormalTolerance = 1.0e-10;

    bool operator==(
        const SceneFluidMimeticControlCellSettings&) const = default;
};

struct SceneFluidMimeticControlCellLimits {
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumHalfFaces = 200'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticHalfFace {
    std::size_t halfFaceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    SceneFluidMimeticHalfFaceSide side =
        SceneFluidMimeticHalfFaceSide::MinusOrNegative;
    std::size_t sourceIndex = 0;
    std::uint64_t sourceStableId = 0;
    std::uint64_t traceStableId = 0;
    std::size_t controlVolumeIndex = 0;
    std::size_t otherControlVolumeIndex =
        invalidSceneFluidMimeticControlVolumeIndex;
    double areaSquareMeters = 0.0;
    fluid::Vector3 centroidMeters;
    fluid::Vector3 periodicImageOffsetMeters;
    fluid::Vector3 outwardUnitNormal;

    bool operator==(const SceneFluidMimeticHalfFace&) const = default;
};

struct SceneFluidMimeticControlCell {
    std::size_t controlCellIndex = 0;
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    std::size_t firstHalfFace = 0;
    std::size_t halfFaceCount = 0;
    std::size_t cartesianHalfFaceCount = 0;
    std::size_t materialWallHalfFaceCount = 0;
    std::size_t openingHalfFaceCount = 0;
    std::size_t unresolvedCartesianIncidentCount = 0;
    std::size_t missingOtherControlCount = 0;
    double volumeCubicMeters = 0.0;
    fluid::Vector3 centroidMeters;
    fluid::Vector3 areaClosureVectorSquareMeters;
    std::array<double, 9> normalMomentMatrixCubicMeters{};
    double maximumAreaClosureErrorSquareMeters = 0.0;
    double maximumDivergenceTheoremErrorCubicMeters = 0.0;
    bool topologyComplete = false;
    bool areaVectorClosed = false;
    bool divergenceTheoremClosed = false;
    bool readyForMimeticOperator = false;

    bool operator==(const SceneFluidMimeticControlCell&) const = default;
};

// Audit-only scene adapter for the mixed-hybrid local kernel. Every positive
// cell/region pressure volume owns one control cell. Resolved Cartesian region
// subfaces and authored openings are shared traces; material quadrature is an
// impermeable boundary on each authored side. Centroids are unwrapped into the
// nearest periodic image of their owning control centroid. A material side
// whose cell/region volume is exactly absent is counted and omitted; a missing
// authored-opening side remains incomplete because shared flux continuity
// would have no owner. The product retains incomplete topology and geometric
// nonclosure per control instead of fabricating a coefficient or changing the
// production pressure graph.
struct SceneFluidMimeticControlCellSet {
    std::uint32_t version = sceneFluidMimeticControlCellVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t gridEpochFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidMimeticControlCellSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t readyControlCellCount = 0;
    std::size_t incompleteTopologyControlCellCount = 0;
    std::size_t nonclosingControlCellCount = 0;
    std::size_t unresolvedCartesianFaceCount = 0;
    std::size_t omittedZeroVolumeMaterialSideCount = 0;
    std::size_t missingOpeningControlSideCount = 0;
    std::size_t cartesianHalfFaceCount = 0;
    std::size_t materialWallHalfFaceCount = 0;
    std::size_t openingHalfFaceCount = 0;
    std::size_t maximumHalfFaceCountPerControl = 0;
    double maximumAreaClosureErrorSquareMeters = 0.0;
    double maximumDivergenceTheoremErrorCubicMeters = 0.0;
    std::vector<SceneFluidMimeticControlCell> controlCells;
    std::vector<SceneFluidMimeticHalfFace> halfFaces;

    bool operator==(
        const SceneFluidMimeticControlCellSet&) const = default;
};

[[nodiscard]] SceneFluidMimeticControlCellSet
buildSceneFluidMimeticControlCells(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticControlCellSettings& settings = {},
    const SceneFluidMimeticControlCellLimits& limits = {});

void validateSceneFluidMimeticControlCellIntegrity(
    const SceneFluidMimeticControlCellSet& controlCells);

void validateSceneFluidMimeticControlCells(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks);

} // namespace simwing::fsi
