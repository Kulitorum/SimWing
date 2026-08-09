#pragma once

#include "scene_fluid_pressure_face_link.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureOperatorVersion = 1;

struct SceneFluidPressureOperatorLimits {
    std::size_t maximumRows = 50'000'000;
    std::size_t maximumEntries = 200'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumOperatorBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureOperatorRow {
    std::size_t rowIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t controlVolumeIndex = 0;
    std::size_t componentIndex = 0;
    bool isGauge = false;
    std::size_t firstEntry = 0;
    std::size_t entryCount = 0;
    double diagonalGeometryWeightMeters = 0.0;

    bool operator==(const SceneFluidPressureOperatorRow&) const = default;
};

struct SceneFluidPressureOperatorEntry {
    std::size_t entryIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::size_t columnControlVolumeIndex = 0;
    double geometryWeightMeters = 0.0;

    bool operator==(const SceneFluidPressureOperatorEntry&) const = default;
};

struct SceneFluidPressureOperatorComponent {
    std::size_t componentIndex = 0;
    StableId gaugeRegionId = invalidStableId;
    std::size_t gaugeControlVolumeIndex = 0;
    std::size_t firstControlVolumeMember = 0;
    std::size_t controlVolumeCount = 0;
    double totalVolumeCubicMeters = 0.0;
    double totalGeometryWeightMeters = 0.0;

    bool operator==(
        const SceneFluidPressureOperatorComponent&) const = default;
};

// Integrated finite-volume graph Laplacian. Every conservative face link
// contributes +w to both endpoint diagonals and -w to both directed row
// entries, where w is area / centre distance. The operator is deliberately
// ungauged and therefore has one constant null mode per pressure component;
// the existing deterministic gauge control volume is retained as metadata for
// a future solve. Construction requires every Cartesian face to be resolved
// and the link graph to contain exactly one connected graph per authored
// pressure component. Missing opening/interface ownership therefore rejects
// instead of silently becoming a no-flow boundary.
struct SceneFluidPressureOperator {
    std::uint32_t version = sceneFluidPressureOperatorVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    double totalGeometryWeightMeters = 0.0;
    std::vector<SceneFluidPressureOperatorRow> rows;
    std::vector<SceneFluidPressureOperatorEntry> entries;
    std::vector<SceneFluidPressureOperatorComponent> components;
    std::vector<std::size_t> componentControlVolumeIndices;

    bool operator==(const SceneFluidPressureOperator&) const = default;
};

[[nodiscard]] SceneFluidPressureOperator buildSceneFluidPressureOperator(
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
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperatorLimits& limits = {});

void validateSceneFluidPressureOperator(
    const SceneFluidPressureOperator& pressureOperator,
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
    const SceneFluidPressureFaceLinkSet& faceLinks);

[[nodiscard]] std::vector<double> applySceneFluidPressureOperator(
    const SceneFluidPressureOperator& pressureOperator,
    std::span<const double> pressureValues);

} // namespace simwing::fsi
