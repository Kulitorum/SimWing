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

struct SceneFluidPressureSolveSettings {
    double absoluteResidualTolerancePascalsMeters = 1.0e-12;
    double relativeResidualTolerance = 1.0e-10;
    double absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-12;
    std::size_t maximumIterations = 4000;

    bool operator==(const SceneFluidPressureSolveSettings&) const = default;
};

struct SceneFluidPressureSolveComponentDiagnostics {
    std::size_t componentIndex = 0;
    std::size_t controlVolumeCount = 0;
    std::size_t gaugeControlVolumeIndex = 0;
    double rightHandSideSumPascalsMeters = 0.0;
    double compatibilityCorrectionPascalsMeters = 0.0;
    double pressureGaugeBeforePascals = 0.0;
    double pressureGaugeAfterPascals = 0.0;

    bool operator==(
        const SceneFluidPressureSolveComponentDiagnostics&) const = default;
};

struct SceneFluidPressureSolveDiagnostics {
    bool compatible = false;
    bool converged = false;
    bool finite = false;
    std::uint64_t pressureOperatorFingerprint = 0;
    std::size_t rowCount = 0;
    std::size_t componentCount = 0;
    std::size_t iterationCount = 0;
    double maximumAbsoluteComponentCompatibilityPascalsMeters = 0.0;
    double initialResidualL2PascalsMeters = 0.0;
    double finalResidualL2PascalsMeters = 0.0;
    double finalResidualMaximumPascalsMeters = 0.0;
    std::vector<SceneFluidPressureSolveComponentDiagnostics> components;

    bool operator==(const SceneFluidPressureSolveDiagnostics&) const = default;
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
// a future solve. Construction requires every Cartesian face and embedded
// opening patch to be resolved, and the link graph to contain exactly one
// connected graph per authored pressure component. Missing opening/interface
// ownership therefore rejects instead of silently becoming a no-flow boundary.
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

// Solves A*p=b in the operator's integrated units. Each component RHS must
// sum to zero within the declared absolute tolerance; only admitted roundoff
// is removed. Pressure is committed only after an explicitly recomputed
// residual converges, then shifted so every retained gauge control volume is
// exactly zero. Incompatibility or non-convergence leaves the caller's warm
// start bit-for-bit unchanged.
[[nodiscard]] SceneFluidPressureSolveDiagnostics
solveSceneFluidPressureSystem(
    const SceneFluidPressureOperator& pressureOperator,
    std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& pressurePascals,
    const SceneFluidPressureSolveSettings& settings = {});

} // namespace simwing::fsi
