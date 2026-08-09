#pragma once

#include "scene_fluid_pressure_face_link.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureTopologyTransitionVersion = 1;

struct SceneFluidPressureTopologyTransitionLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumMappings = 100'000'000;
    std::size_t maximumTransitionBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureRetainedControl {
    std::size_t previousControlVolumeIndex = 0;
    std::size_t currentControlVolumeIndex = 0;
    std::uint64_t stableId = 0;

    bool operator==(
        const SceneFluidPressureRetainedControl&) const = default;
};

struct SceneFluidPressureAppearanceDonor {
    std::size_t appearedCurrentControlVolumeIndex = 0;
    std::size_t retainedPreviousControlVolumeIndex = 0;
    double linkAreaSquareMeters = 0.0;
    double normalizedWeight = 0.0;

    bool operator==(
        const SceneFluidPressureAppearanceDonor&) const = default;
};

struct SceneFluidPressureRetirementRecipient {
    std::size_t disappearedPreviousControlVolumeIndex = 0;
    std::size_t retainedCurrentControlVolumeIndex = 0;
    double linkAreaSquareMeters = 0.0;
    double normalizedWeight = 0.0;

    bool operator==(
        const SceneFluidPressureRetirementRecipient&) const = default;
};

// One immutable consecutive-epoch topology transaction shared by pressure
// volume rates, transported region state, and pressure warm starts. Stable
// rows are paired explicitly. An appeared row receives current-topology
// same-region donors that are retained across the transition. A disappeared
// row currently admits exactly one unique retained same-region recipient from
// the previous topology. The explicit mapping is the common policy boundary;
// it does not claim a general swept-volume remap.
struct SceneFluidPressureTopologyTransition {
    std::uint32_t version =
        sceneFluidPressureTopologyTransitionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t previousPressureControlVolumeFingerprint = 0;
    std::uint64_t previousPressureFaceLinkFingerprint = 0;
    std::uint64_t currentPressureControlVolumeFingerprint = 0;
    std::uint64_t currentPressureFaceLinkFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousCellVolumeFingerprint = 0;
    std::uint64_t currentCellVolumeFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    std::size_t previousControlVolumeCount = 0;
    std::size_t currentControlVolumeCount = 0;
    std::size_t retainedControlVolumeCount = 0;
    std::size_t appearedControlVolumeCount = 0;
    std::size_t disappearedControlVolumeCount = 0;
    std::size_t maximumAppearanceDonorCount = 0;
    std::size_t maximumRetirementRecipientCount = 0;
    double totalAppearanceDonorLinkAreaSquareMeters = 0.0;
    double totalRetirementRecipientLinkAreaSquareMeters = 0.0;
    std::vector<SceneFluidPressureRetainedControl> retainedControls;
    std::vector<SceneFluidPressureAppearanceDonor> appearanceDonors;
    std::vector<SceneFluidPressureRetirementRecipient>
        retirementRecipients;

    bool operator==(
        const SceneFluidPressureTopologyTransition&) const = default;
};

[[nodiscard]] SceneFluidPressureTopologyTransition
buildSceneFluidPressureTopologyTransition(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidPressureTopologyTransitionLimits& limits = {});

void validateSceneFluidPressureTopologyTransitionIntegrity(
    const SceneFluidPressureTopologyTransition& transition);

void validateSceneFluidPressureTopologyTransition(
    const SceneFluidPressureTopologyTransition& transition,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks);

} // namespace simwing::fsi
