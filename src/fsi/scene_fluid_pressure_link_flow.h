#pragma once

#include "scene_fluid_pressure_projection.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureLinkFlowContinuationVersion = 1;

struct SceneFluidPressureLinkFlowContinuationLimits {
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumContinuationBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureContinuedLinkFlow {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t faceIndex = 0;
    SceneFluidPressureFaceLinkKind kind =
        SceneFluidPressureFaceLinkKind::SameRegion;
    std::uint64_t openingPatchStableId = 0;
    double carriedAbsoluteVelocityDeviationMetersPerSecond = 0.0;
    double predictedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidPressureContinuedLinkFlow&) const = default;
};

struct SceneFluidPressureLinkFlowContinuationDiagnostics {
    std::size_t faceCount = 0;
    std::size_t linkCount = 0;
    std::size_t openingLinkCount = 0;
    std::size_t multiLinkFaceCount = 0;
    double maximumCarriedAbsoluteVelocityDeviationMetersPerSecond = 0.0;
    double maximumAreaRenormalizationMetersPerSecond = 0.0;
    double maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteFaceFlowClosureCubicMetersPerSecond = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureLinkFlowContinuationDiagnostics&) const =
        default;
};

// Topology-stable memory for link-resolved pressure flow. The accepted
// absolute velocity deviation of each previous subface from its Cartesian
// face mean is carried to the matching current link. Current-area weighted
// deviations are recentered to zero on every face, so the supplied bulk MAC
// field remains the exact face-total flow owner. Authored-opening relative
// flow uses the current cap sweep from openingFlux. This is a continuation
// model only: it neither advects the subface residual nor applies wall
// viscosity.
struct SceneFluidPressureLinkFlowContinuation {
    std::uint32_t version =
        sceneFluidPressureLinkFlowContinuationVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t previousPressureProjectionFingerprint = 0;
    std::uint64_t previousPressureFaceLinkFingerprint = 0;
    std::uint64_t previousOpeningPatchFingerprint = 0;
    std::uint64_t currentPressureFaceLinkFingerprint = 0;
    std::uint64_t currentOpeningFluxFingerprint = 0;
    std::uint64_t currentVelocityFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    SceneFluidPressureLinkFlowContinuationDiagnostics diagnostics;
    std::vector<SceneFluidPressureContinuedLinkFlow> links;

    bool operator==(
        const SceneFluidPressureLinkFlowContinuation&) const = default;
};

[[nodiscard]] SceneFluidPressureLinkFlowContinuation
continueSceneFluidPressureLinkFlows(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidOpeningGridPatchSet& previousOpeningPatches,
    const SceneFluidPressureProjection& previousProjection,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux,
    const fluid::MacVelocityField& currentBulkVelocityMetersPerSecond,
    const SceneFluidPressureLinkFlowContinuationLimits& limits = {});

void validateSceneFluidPressureLinkFlowContinuationIntegrity(
    const SceneFluidPressureLinkFlowContinuation& continuation);

void validateSceneFluidPressureLinkFlowContinuation(
    const SceneFluidPressureLinkFlowContinuation& continuation,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux);

} // namespace simwing::fsi
