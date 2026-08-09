#pragma once

#include "scene_fluid_pressure_control_volume.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureVolumeRateVersion = 1;

struct SceneFluidPressureVolumeRateLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumVolumeRateBytes =
        1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureControlVolumeRate {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double volumeChangeCubicMeters = 0.0;
    double volumeChangeRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidPressureControlVolumeRate&) const = default;
};

struct SceneFluidPressureComponentVolumeRate {
    std::size_t componentIndex = 0;
    std::size_t controlVolumeCount = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double volumeChangeCubicMeters = 0.0;
    double volumeChangeRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidPressureComponentVolumeRate&) const = default;
};

// Consecutive accepted sparse pressure volumes are reduced to one exact
// geometry-volume rate per retained cell/region unknown. The cell-major stable
// IDs, component ownership, and gauge topology must be identical at both
// endpoints. Appearance/disappearance of a positive cut-cell region is a
// topology rebase and rejects here until a conservative remap owner exists.
// This immutable product supplies only dV/dt; it does not sample velocity,
// construct a pressure RHS, or modify either endpoint.
struct SceneFluidPressureVolumeRateSet {
    std::uint32_t version = sceneFluidPressureVolumeRateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousSurfaceStateFingerprint = 0;
    std::uint64_t currentSurfaceStateFingerprint = 0;
    std::uint64_t previousCellVolumeFingerprint = 0;
    std::uint64_t currentCellVolumeFingerprint = 0;
    std::uint64_t currentPressureControlVolumeFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    double durationSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    double maximumAbsoluteControlVolumeChangeCubicMeters = 0.0;
    double maximumAbsoluteControlVolumeRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteComponentVolumeRateCubicMetersPerSecond = 0.0;
    double globalVolumeChangeCubicMeters = 0.0;
    double globalVolumeChangeRateCubicMetersPerSecond = 0.0;
    std::vector<SceneFluidPressureControlVolumeRate> controlVolumes;
    std::vector<SceneFluidPressureComponentVolumeRate> components;

    bool operator==(
        const SceneFluidPressureVolumeRateSet&) const = default;
};

[[nodiscard]] SceneFluidPressureVolumeRateSet
buildSceneFluidPressureVolumeRates(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureVolumeRateLimits& limits = {});

void validateSceneFluidPressureVolumeRateIntegrity(
    const SceneFluidPressureVolumeRateSet& rates);

void validateSceneFluidPressureVolumeRates(
    const SceneFluidPressureVolumeRateSet& rates,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes);

} // namespace simwing::fsi
