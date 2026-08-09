#pragma once

#include "scene_fluid_opening_flux.h"
#include "scene_fluid_pressure_operator.h"
#include "scene_fluid_pressure_volume_rate.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureProjectionVersion = 6;

struct SceneFluidPressureLinkFlowContinuation;
struct SceneFluidRegionLinkFlowPrediction;

struct SceneFluidPressureProjectionSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    double absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond =
        1.0e-12;
    double relativeCorrectedVolumeRateTolerance = 1.0e-10;
    SceneFluidPressureSolveSettings pressureSolve;

    bool operator==(
        const SceneFluidPressureProjectionSettings&) const = default;
};

struct SceneFluidPressureProjectionLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumProjectionBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureProjectedLink {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t faceIndex = 0;
    SceneFluidPressureFaceLinkKind kind =
        SceneFluidPressureFaceLinkKind::SameRegion;
    std::size_t minusControlVolumeIndex = 0;
    std::size_t plusControlVolumeIndex = 0;
    std::uint64_t openingPatchStableId = 0;
    double predictedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double pressureCorrectionVolumeFlowRateCubicMetersPerSecond = 0.0;
    double correctedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(const SceneFluidPressureProjectedLink&) const = default;
};

struct SceneFluidPressureProjectedControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t componentIndex = 0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;
    double predictedNetOutwardVolumeFlowRateCubicMetersPerSecond = 0.0;
    double predictedContinuityResidualCubicMetersPerSecond = 0.0;
    double integratedRightHandSidePascalsMeters = 0.0;
    double correctedNetOutwardVolumeFlowRateCubicMetersPerSecond = 0.0;
    double correctedContinuityResidualCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidPressureProjectedControlVolume&) const = default;
};

struct SceneFluidPressureProjectionDiagnostics {
    bool accepted = false;
    bool finite = false;
    std::size_t controlVolumeCount = 0;
    std::size_t linkCount = 0;
    std::size_t authoredOpeningLinkCount = 0;
    bool usesMovingVolumeRates = false;
    double maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond = 0.0;
    double predictedNetOutwardVolumeRateL2CubicMetersPerSecond = 0.0;
    double predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond = 0.0;
    double maximumPredictedComponentBalanceResidualCubicMetersPerSecond =
        0.0;
    double predictedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double predictedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double maximumPredictedComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double correctedNetOutwardVolumeRateL2CubicMetersPerSecond = 0.0;
    double correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond = 0.0;
    double maximumCorrectedComponentBalanceResidualCubicMetersPerSecond =
        0.0;
    double correctedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double correctedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double maximumCorrectedComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    SceneFluidPressureSolveDiagnostics pressureSolve;

    bool operator==(
        const SceneFluidPressureProjectionDiagnostics&) const = default;
};

// Link-resolved finite-volume projection. Same-region links read
// their owning MAC normal velocity. Face-owned and embedded authored-opening
// links instead reuse the accepted negative-to-positive relative opening flux
// (fluid minus cap sweep) and orient it from the link's minus control to its
// plus control. The integrated
// pressure solve corrects one flow per link; a partitioned Cartesian face is
// deliberately not collapsed back to one MAC value.
//
// The result is an immutable attempt. Pressure and corrected flows are present
// only when the pressure solve and the explicit post-correction continuity
// check both succeed. The base overload holds accepted geometry fixed. The
// consecutive-epoch overload additionally enforces dV/dt + net flow = 0 from
// a topology-stable pressure-volume-rate product. Cell-crossing topology
// rebases remain outside this contract.
struct SceneFluidPressureProjection {
    std::uint32_t version = sceneFluidPressureProjectionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t pressureOperatorFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t pressureVolumeRateFingerprint = 0;
    std::uint64_t openingFluxFingerprint = 0;
    std::uint64_t velocityFingerprint = 0;
    std::uint64_t linkFlowContinuationFingerprint = 0;
    std::uint64_t regionLinkFlowPredictionFingerprint = 0;
    std::uint64_t regionWallExchangeFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidPressureProjectionSettings settings;
    std::size_t ownedStorageBytes = 0;
    SceneFluidPressureProjectionDiagnostics diagnostics;
    std::vector<double> pressurePascals;
    std::vector<SceneFluidPressureProjectedControlVolume> controlVolumes;
    std::vector<SceneFluidPressureProjectedLink> links;

    bool operator==(const SceneFluidPressureProjection&) const = default;
};

[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

// Uses the moving-epoch link predictor reconstructed from an accepted region
// momentum transport candidate. The supplied opening-flux/MAC chain remains
// the immutable geometry and fallback identity; its bulk link values are
// overridden by the region prediction.
[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidRegionLinkFlowPrediction& regionLinkFlowPrediction,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidRegionLinkFlowPrediction& regionLinkFlowPrediction,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    const SceneFluidPressureVolumeRateSet& volumeRates,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

// Uses an explicitly continued relative flow for every current face link.
// The continuation remains bound to the same bulk MAC/opening-flux input and
// changes only the link-resolved predictor used to assemble the pressure RHS.
[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureLinkFlowContinuation& linkFlowContinuation,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    const SceneFluidPressureVolumeRateSet& volumeRates,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

[[nodiscard]] SceneFluidPressureProjection
projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureLinkFlowContinuation& linkFlowContinuation,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    const SceneFluidPressureVolumeRateSet& volumeRates,
    std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings = {},
    const SceneFluidPressureProjectionLimits& limits = {});

void validateSceneFluidPressureProjectionIntegrity(
    const SceneFluidPressureProjection& projection);

} // namespace simwing::fsi
