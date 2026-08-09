#pragma once

#include "scene_fluid_region_transport.h"
#include "scene_fluid_region_wall.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionLinkFlowVersion = 2;

struct SceneFluidRegionLinkFlowLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumLinkFlowBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionPredictedLinkFlow {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t faceIndex = 0;
    SceneFluidPressureFaceLinkKind kind =
        SceneFluidPressureFaceLinkKind::SameRegion;
    std::uint64_t openingPatchStableId = 0;
    double predictedAbsoluteVelocityMetersPerSecond = 0.0;
    double predictedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidRegionPredictedLinkFlow&) const = default;
};

struct SceneFluidRegionLinkFlowDiagnostics {
    std::size_t controlVolumeCount = 0;
    std::size_t faceCount = 0;
    std::size_t linkCount = 0;
    std::size_t openingLinkCount = 0;
    std::size_t multiLinkFaceCount = 0;
    fluid::Vector3 sourceMomentumKilogramMetersPerSecond;
    fluid::Vector3 remappedMomentumKilogramMetersPerSecond;
    fluid::Vector3 geometricMomentumChangeKilogramMetersPerSecond;
    double maximumAbsoluteVolumeChangeCubicMeters = 0.0;
    double maximumAbsolutePredictedVelocityMetersPerSecond = 0.0;
    double maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidRegionLinkFlowDiagnostics&) const = default;
};

// Maps one accepted fixed-epoch region transport candidate onto a consecutive
// topology-stable moving pressure epoch. Stable cell/region velocities are
// retained while momentum is recomputed from the current physical volumes.
// Each pressure link receives the arithmetic mean endpoint normal velocity;
// an authored opening additionally subtracts its exact current oriented cap
// sweep. The wall-source overload uses already-adjusted current control-volume
// velocities and records that source explicitly. This is a first-order GCL
// remap/predictor, not a pressure correction or topology rebase.
struct SceneFluidRegionLinkFlowPrediction {
    std::uint32_t version = sceneFluidRegionLinkFlowVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceWallExchangeFingerprint = 0;
    std::uint64_t currentPressureControlVolumeFingerprint = 0;
    std::uint64_t currentPressureFaceLinkFingerprint = 0;
    std::uint64_t currentOpeningFluxFingerprint = 0;
    std::uint64_t currentVelocityFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionLinkFlowDiagnostics diagnostics;
    std::vector<SceneFluidRegionPredictedLinkFlow> links;

    bool operator==(
        const SceneFluidRegionLinkFlowPrediction&) const = default;
};

[[nodiscard]] SceneFluidRegionLinkFlowPrediction
predictSceneFluidRegionLinkFlows(
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux,
    const SceneFluidRegionLinkFlowLimits& limits = {});

[[nodiscard]] SceneFluidRegionLinkFlowPrediction
predictSceneFluidRegionLinkFlows(
    const SceneFluidRegionWallExchange& wallExchange,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux,
    const SceneFluidRegionLinkFlowLimits& limits = {});

void validateSceneFluidRegionLinkFlowPredictionIntegrity(
    const SceneFluidRegionLinkFlowPrediction& prediction);

void validateSceneFluidRegionLinkFlowPrediction(
    const SceneFluidRegionLinkFlowPrediction& prediction,
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux);

void validateSceneFluidRegionLinkFlowPrediction(
    const SceneFluidRegionLinkFlowPrediction& prediction,
    const SceneFluidRegionWallExchange& wallExchange,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux);

} // namespace simwing::fsi
