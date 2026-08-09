#pragma once

#include "scene_fluid_region_transport.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionRebaseVersion = 1;

struct SceneFluidRegionRebaseLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumRebaseBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionRebaseControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    bool appearedThisEpoch = false;
    std::size_t sourceControlVolumeIndex =
        std::numeric_limits<std::size_t>::max();
    std::size_t donorControlVolumeCount = 0;
    double donorLinkAreaSquareMeters = 0.0;
    double volumeCubicMeters = 0.0;
    fluid::Vector3 velocityMetersPerSecond;
    fluid::Vector3 momentumKilogramMetersPerSecond;

    bool operator==(
        const SceneFluidRegionRebaseControlVolume&) const = default;
};

struct SceneFluidRegionRebaseDiagnostics {
    std::size_t previousControlVolumeCount = 0;
    std::size_t currentControlVolumeCount = 0;
    std::size_t retainedControlVolumeCount = 0;
    std::size_t appearedControlVolumeCount = 0;
    std::size_t maximumDonorControlVolumeCount = 0;
    fluid::Vector3 sourceMomentumKilogramMetersPerSecond;
    fluid::Vector3 rebasedMomentumKilogramMetersPerSecond;
    fluid::Vector3 geometricMomentumChangeKilogramMetersPerSecond;
    double sourceKineticEnergyJoules = 0.0;
    double rebasedKineticEnergyJoules = 0.0;
    double maximumAbsoluteVolumeChangeCubicMeters = 0.0;
    double maximumAppearedVolumeCubicMeters = 0.0;
    double maximumAbsoluteVelocityMetersPerSecond = 0.0;
    bool finite = false;

    bool operator==(const SceneFluidRegionRebaseDiagnostics&) const =
        default;
};

// Transactionally maps one accepted transported region state onto a current
// pressure topology in which old controls are retained and one or more new
// positive controls may appear. Retained rows keep their transported velocity.
// Each appeared row receives the area-weighted velocity of directly linked,
// retained controls in the same authored region; no cross-material or opening
// donor is admitted. Momentum is then recomputed from current physical volume,
// matching the topology-stable geometric remap. A disappeared old row or an
// appeared row without a retained one-ring donor rejects before publication.
struct SceneFluidRegionRebase {
    std::uint32_t version = sceneFluidRegionRebaseVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t previousPressureControlVolumeFingerprint = 0;
    std::uint64_t currentPressureControlVolumeFingerprint = 0;
    std::uint64_t currentPressureFaceLinkFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionRebaseDiagnostics diagnostics;
    std::vector<SceneFluidRegionRebaseControlVolume> controlVolumes;

    bool operator==(const SceneFluidRegionRebase&) const = default;
};

[[nodiscard]] SceneFluidRegionRebase rebaseSceneFluidRegionTransport(
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidRegionRebaseLimits& limits = {});

void validateSceneFluidRegionRebaseIntegrity(
    const SceneFluidRegionRebase& rebase);

void validateSceneFluidRegionRebase(
    const SceneFluidRegionRebase& rebase,
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks);

// Maps an accepted pressure warm start onto the same supported retained-plus-
// appearance subset. Retained stable IDs keep their exact pressure; an
// appeared control receives the area-weighted pressure of its retained,
// same-region one-ring neighbours. This is solver initialization only and
// does not publish or alter an accepted pressure state.
[[nodiscard]] std::vector<double> rebaseSceneFluidPressureWarmStart(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    std::span<const double> previousPressurePascals);

} // namespace simwing::fsi
