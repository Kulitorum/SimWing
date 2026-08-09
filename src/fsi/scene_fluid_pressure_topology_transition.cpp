#include "scene_fluid_pressure_topology_transition.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

template<typename Value>
bool checkedStorageAdd(const std::size_t count, std::size_t& total) {
    if (count > std::numeric_limits<std::size_t>::max() / sizeof(Value)) {
        return false;
    }
    std::size_t next = 0;
    return checkedAdd(total, count * sizeof(Value), next)
        ? (total = next, true) : false;
}

std::size_t storageBytesForCounts(
    const std::size_t retainedCount,
    const std::size_t appearanceDonorCount,
    const std::size_t retirementRecipientCount) {
    std::size_t total = 0;
    if (!checkedStorageAdd<SceneFluidPressureRetainedControl>(
            retainedCount, total)
        || !checkedStorageAdd<SceneFluidPressureAppearanceDonor>(
            appearanceDonorCount, total)
        || !checkedStorageAdd<SceneFluidPressureRetirementRecipient>(
            retirementRecipientCount, total)) {
        throw std::length_error(
            "scene fluid pressure topology-transition storage overflows");
    }
    return total;
}

std::size_t storageBytes(
    const SceneFluidPressureTopologyTransition& transition) {
    return storageBytesForCounts(
        transition.retainedControls.size(),
        transition.appearanceDonors.size(),
        transition.retirementRecipients.size());
}

std::uint64_t transitionFingerprint(
    const SceneFluidPressureTopologyTransition& transition) {
    Fingerprint fingerprint;
    fingerprint.integer(transition.version);
    for (const std::uint64_t value : {
             transition.previousPressureControlVolumeFingerprint,
             transition.previousPressureFaceLinkFingerprint,
             transition.currentPressureControlVolumeFingerprint,
             transition.currentPressureFaceLinkFingerprint,
             transition.surfaceDefinitionFingerprint,
             transition.structureDefinitionFingerprint,
             transition.previousCellVolumeFingerprint,
             transition.currentCellVolumeFingerprint,
             transition.previousAcceptedStepCount,
             transition.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    for (const double value : {
             transition.previousSimulationTimeSeconds,
             transition.currentSimulationTimeSeconds,
             transition.lowerMeters.x,
             transition.lowerMeters.y,
             transition.lowerMeters.z,
             transition.upperMeters.x,
             transition.upperMeters.y,
             transition.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.cellCounts.z));
    for (const std::size_t value : {
             transition.ownedStorageBytes,
             transition.previousControlVolumeCount,
             transition.currentControlVolumeCount,
             transition.retainedControlVolumeCount,
             transition.appearedControlVolumeCount,
             transition.disappearedControlVolumeCount,
             transition.maximumAppearanceDonorCount,
             transition.maximumRetirementRecipientCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(
        transition.totalAppearanceDonorLinkAreaSquareMeters);
    fingerprint.real(
        transition.totalRetirementRecipientLinkAreaSquareMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.retainedControls.size()));
    for (const auto& retained : transition.retainedControls) {
        fingerprint.integer(static_cast<std::uint64_t>(
            retained.previousControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            retained.currentControlVolumeIndex));
        fingerprint.integer(retained.stableId);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.appearanceDonors.size()));
    for (const auto& donor : transition.appearanceDonors) {
        fingerprint.integer(static_cast<std::uint64_t>(
            donor.appearedCurrentControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            donor.retainedPreviousControlVolumeIndex));
        fingerprint.real(donor.linkAreaSquareMeters);
        fingerprint.real(donor.normalizedWeight);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.retirementRecipients.size()));
    for (const auto& recipient : transition.retirementRecipients) {
        fingerprint.integer(static_cast<std::uint64_t>(
            recipient.disappearedPreviousControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            recipient.retainedCurrentControlVolumeIndex));
        fingerprint.real(recipient.linkAreaSquareMeters);
        fingerprint.real(recipient.normalizedWeight);
    }
    return fingerprint.value();
}

void validateSources(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks) {
    validateSceneFluidPressureControlVolumeIntegrity(previousPressureVolumes);
    validateSceneFluidPressureControlVolumeIntegrity(currentPressureVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(previousFaceLinks);
    validateSceneFluidPressureFaceLinkIntegrity(currentFaceLinks);
    if (previousPressureVolumes.surfaceDefinitionFingerprint == 0
        || previousPressureVolumes.surfaceDefinitionFingerprint
            != currentPressureVolumes.surfaceDefinitionFingerprint
        || previousPressureVolumes.structureDefinitionFingerprint
            != currentPressureVolumes.structureDefinitionFingerprint
        || previousPressureVolumes.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentPressureVolumes.acceptedStepCount
            != previousPressureVolumes.acceptedStepCount + 1
        || !(currentPressureVolumes.simulationTimeSeconds
            > previousPressureVolumes.simulationTimeSeconds)
        || previousPressureVolumes.cellCounts
            != currentPressureVolumes.cellCounts
        || previousPressureVolumes.lowerMeters
            != currentPressureVolumes.lowerMeters
        || previousPressureVolumes.upperMeters
            != currentPressureVolumes.upperMeters
        || previousFaceLinks.pressureControlVolumeFingerprint
            != previousPressureVolumes.fingerprint
        || previousFaceLinks.acceptedStepCount
            != previousPressureVolumes.acceptedStepCount
        || previousFaceLinks.simulationTimeSeconds
            != previousPressureVolumes.simulationTimeSeconds
        || previousFaceLinks.cellCounts != previousPressureVolumes.cellCounts
        || previousFaceLinks.lowerMeters != previousPressureVolumes.lowerMeters
        || previousFaceLinks.upperMeters != previousPressureVolumes.upperMeters
        || currentFaceLinks.pressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || currentFaceLinks.acceptedStepCount
            != currentPressureVolumes.acceptedStepCount
        || currentFaceLinks.simulationTimeSeconds
            != currentPressureVolumes.simulationTimeSeconds
        || currentFaceLinks.cellCounts != currentPressureVolumes.cellCounts
        || currentFaceLinks.lowerMeters != currentPressureVolumes.lowerMeters
        || currentFaceLinks.upperMeters != currentPressureVolumes.upperMeters) {
        throw std::invalid_argument(
            "scene fluid pressure topology-transition source identity is invalid");
    }
}

SceneFluidPressureTopologyTransition buildTransition(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidPressureTopologyTransitionLimits& limits) {
    const std::size_t previousCount =
        previousPressureVolumes.controlVolumes.size();
    const std::size_t currentCount =
        currentPressureVolumes.controlVolumes.size();
    if (previousCount > limits.maximumControlVolumes
        || currentCount > limits.maximumControlVolumes
        || previousFaceLinks.links.size() > limits.maximumLinks
        || currentFaceLinks.links.size()
            > limits.maximumLinks - previousFaceLinks.links.size()) {
        throw std::length_error(
            "scene fluid pressure topology transition exceeds its count limit");
    }

    SceneFluidPressureTopologyTransition result;
    result.previousPressureControlVolumeFingerprint =
        previousPressureVolumes.fingerprint;
    result.previousPressureFaceLinkFingerprint = previousFaceLinks.fingerprint;
    result.currentPressureControlVolumeFingerprint =
        currentPressureVolumes.fingerprint;
    result.currentPressureFaceLinkFingerprint = currentFaceLinks.fingerprint;
    result.surfaceDefinitionFingerprint =
        currentPressureVolumes.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        currentPressureVolumes.structureDefinitionFingerprint;
    result.previousCellVolumeFingerprint =
        previousPressureVolumes.cellVolumeFingerprint;
    result.currentCellVolumeFingerprint =
        currentPressureVolumes.cellVolumeFingerprint;
    result.previousAcceptedStepCount =
        previousPressureVolumes.acceptedStepCount;
    result.currentAcceptedStepCount = currentPressureVolumes.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousPressureVolumes.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentPressureVolumes.simulationTimeSeconds;
    result.cellCounts = currentPressureVolumes.cellCounts;
    result.lowerMeters = currentPressureVolumes.lowerMeters;
    result.upperMeters = currentPressureVolumes.upperMeters;
    result.previousControlVolumeCount = previousCount;
    result.currentControlVolumeCount = currentCount;

    std::map<std::uint64_t, std::size_t> previousByStableId;
    for (const auto& previous : previousPressureVolumes.controlVolumes) {
        previousByStableId.emplace(
            previous.stableId, previous.controlVolumeIndex);
    }
    std::vector<std::size_t> previousIndexByCurrent(
        currentCount, std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> currentIndexByPrevious(
        previousCount, std::numeric_limits<std::size_t>::max());
    for (const auto& current : currentPressureVolumes.controlVolumes) {
        const auto found = previousByStableId.find(current.stableId);
        if (found == previousByStableId.end()) {
            ++result.appearedControlVolumeCount;
            continue;
        }
        const auto& previous =
            previousPressureVolumes.controlVolumes[found->second];
        if (previous.cellIndex != current.cellIndex
            || previous.regionId != current.regionId
            || previous.kind != current.kind
            || previous.componentIndex != current.componentIndex) {
            throw std::invalid_argument(
                "scene fluid pressure topology-transition retained identity changed");
        }
        previousIndexByCurrent[current.controlVolumeIndex] = found->second;
        currentIndexByPrevious[found->second] = current.controlVolumeIndex;
        result.retainedControls.push_back({
            found->second, current.controlVolumeIndex, current.stableId,
        });
    }
    result.retainedControlVolumeCount = result.retainedControls.size();

    std::map<std::size_t, std::map<std::size_t, double>>
        appearanceAreaByCurrent;
    const auto addAppearanceDonor = [&](
        const std::size_t appearedCurrentIndex,
        const std::size_t donorCurrentIndex,
        const double areaSquareMeters) {
        if (previousIndexByCurrent[appearedCurrentIndex]
                != std::numeric_limits<std::size_t>::max()
            || previousIndexByCurrent[donorCurrentIndex]
                == std::numeric_limits<std::size_t>::max()) {
            return;
        }
        const auto& appeared = currentPressureVolumes.controlVolumes[
            appearedCurrentIndex];
        const auto& donor = currentPressureVolumes.controlVolumes[
            donorCurrentIndex];
        if (donor.regionId != appeared.regionId) {
            return;
        }
        appearanceAreaByCurrent[appearedCurrentIndex][
            previousIndexByCurrent[donorCurrentIndex]] += areaSquareMeters;
    };
    for (const auto& link : currentFaceLinks.links) {
        if (link.kind != SceneFluidPressureFaceLinkKind::SameRegion) {
            continue;
        }
        addAppearanceDonor(
            link.minusControlVolumeIndex,
            link.plusControlVolumeIndex,
            link.areaSquareMeters);
        addAppearanceDonor(
            link.plusControlVolumeIndex,
            link.minusControlVolumeIndex,
            link.areaSquareMeters);
    }
    for (const auto& current : currentPressureVolumes.controlVolumes) {
        if (previousIndexByCurrent[current.controlVolumeIndex]
            != std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        const auto foundArea = appearanceAreaByCurrent.find(
            current.controlVolumeIndex);
        if (foundArea == appearanceAreaByCurrent.end()) {
            throw std::invalid_argument(
                "scene fluid pressure topology-transition appeared control "
                "lacks a retained same-region donor");
        }
        const auto& areaByPreviousIndex = foundArea->second;
        double totalArea = 0.0;
        for (const auto& [donor, area] : areaByPreviousIndex) {
            static_cast<void>(donor);
            totalArea += area;
        }
        if (areaByPreviousIndex.empty() || !(totalArea > 0.0)
            || !std::isfinite(totalArea)) {
            throw std::invalid_argument(
                "scene fluid pressure topology-transition appeared control "
                "lacks a retained same-region donor");
        }
        result.maximumAppearanceDonorCount = std::max(
            result.maximumAppearanceDonorCount, areaByPreviousIndex.size());
        result.totalAppearanceDonorLinkAreaSquareMeters += totalArea;
        for (const auto& [donor, area] : areaByPreviousIndex) {
            if (result.appearanceDonors.size()
                == limits.maximumMappings) {
                throw std::length_error(
                    "scene fluid pressure topology transition exceeds its "
                    "mapping limit");
            }
            result.appearanceDonors.push_back({
                current.controlVolumeIndex, donor, area, area / totalArea,
            });
        }
    }

    std::map<std::size_t, std::map<std::size_t, double>>
        retirementAreaByPrevious;
    const auto addRetirementRecipient = [&](
        const std::size_t disappearedPreviousIndex,
        const std::size_t neighbourPreviousIndex,
        const double areaSquareMeters) {
        if (currentIndexByPrevious[disappearedPreviousIndex]
                != std::numeric_limits<std::size_t>::max()
            || currentIndexByPrevious[neighbourPreviousIndex]
                == std::numeric_limits<std::size_t>::max()) {
            return;
        }
        const auto& disappeared = previousPressureVolumes.controlVolumes[
            disappearedPreviousIndex];
        const auto& neighbour = previousPressureVolumes.controlVolumes[
            neighbourPreviousIndex];
        const std::size_t recipientCurrentIndex =
            currentIndexByPrevious[neighbourPreviousIndex];
        const auto& recipient = currentPressureVolumes.controlVolumes[
            recipientCurrentIndex];
        if (neighbour.regionId != disappeared.regionId
            || recipient.componentIndex != disappeared.componentIndex) {
            return;
        }
        retirementAreaByPrevious[disappearedPreviousIndex][
            recipientCurrentIndex] += areaSquareMeters;
    };
    for (const auto& link : previousFaceLinks.links) {
        if (link.kind != SceneFluidPressureFaceLinkKind::SameRegion) {
            continue;
        }
        addRetirementRecipient(
            link.minusControlVolumeIndex,
            link.plusControlVolumeIndex,
            link.areaSquareMeters);
        addRetirementRecipient(
            link.plusControlVolumeIndex,
            link.minusControlVolumeIndex,
            link.areaSquareMeters);
    }
    for (const auto& previous : previousPressureVolumes.controlVolumes) {
        if (currentIndexByPrevious[previous.controlVolumeIndex]
            != std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        ++result.disappearedControlVolumeCount;
        const auto foundArea = retirementAreaByPrevious.find(
            previous.controlVolumeIndex);
        if (foundArea == retirementAreaByPrevious.end()) {
            throw std::invalid_argument(
                "scene fluid pressure topology-transition disappeared "
                "control lacks one unique retained same-region recipient");
        }
        const auto& areaByCurrentIndex = foundArea->second;
        double totalArea = 0.0;
        for (const auto& [recipient, area] : areaByCurrentIndex) {
            static_cast<void>(recipient);
            totalArea += area;
        }
        if (areaByCurrentIndex.size() != 1 || !(totalArea > 0.0)
            || !std::isfinite(totalArea)) {
            throw std::invalid_argument(
                "scene fluid pressure topology-transition disappeared "
                "control lacks one unique retained same-region recipient");
        }
        result.maximumRetirementRecipientCount = std::max(
            result.maximumRetirementRecipientCount,
            areaByCurrentIndex.size());
        result.totalRetirementRecipientLinkAreaSquareMeters += totalArea;
        for (const auto& [recipient, area] : areaByCurrentIndex) {
            if (result.appearanceDonors.size()
                    > limits.maximumMappings
                || result.retirementRecipients.size()
                    == limits.maximumMappings
                        - result.appearanceDonors.size()) {
                throw std::length_error(
                    "scene fluid pressure topology transition exceeds its "
                    "mapping limit");
            }
            result.retirementRecipients.push_back({
                previous.controlVolumeIndex, recipient, area,
                area / totalArea,
            });
        }
    }

    const std::size_t mappingCount = result.appearanceDonors.size()
        + result.retirementRecipients.size();
    if (mappingCount < result.appearanceDonors.size()
        || mappingCount > limits.maximumMappings) {
        throw std::length_error(
            "scene fluid pressure topology transition exceeds its mapping limit");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumTransitionBytes) {
        throw std::length_error(
            "scene fluid pressure topology transition exceeds its byte limit");
    }
    if (!std::isfinite(
            result.totalAppearanceDonorLinkAreaSquareMeters)
        || !std::isfinite(
            result.totalRetirementRecipientLinkAreaSquareMeters)) {
        throw std::overflow_error(
            "scene fluid pressure topology transition is non-finite");
    }
    result.fingerprint = transitionFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureTopologyTransition
buildSceneFluidPressureTopologyTransition(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidPressureTopologyTransitionLimits& limits) {
    validateSources(
        previousPressureVolumes, previousFaceLinks,
        currentPressureVolumes, currentFaceLinks);
    auto result = buildTransition(
        previousPressureVolumes, previousFaceLinks,
        currentPressureVolumes, currentFaceLinks, limits);
    validateSceneFluidPressureTopologyTransition(
        result, previousPressureVolumes, previousFaceLinks,
        currentPressureVolumes, currentFaceLinks);
    return result;
}

void validateSceneFluidPressureTopologyTransitionIntegrity(
    const SceneFluidPressureTopologyTransition& transition) {
    if (transition.version
            != sceneFluidPressureTopologyTransitionVersion
        || transition.fingerprint == 0
        || transition.previousPressureControlVolumeFingerprint == 0
        || transition.previousPressureFaceLinkFingerprint == 0
        || transition.currentPressureControlVolumeFingerprint == 0
        || transition.currentPressureFaceLinkFingerprint == 0
        || transition.previousControlVolumeCount == 0
        || transition.currentControlVolumeCount == 0
        || transition.retainedControlVolumeCount
                + transition.disappearedControlVolumeCount
            != transition.previousControlVolumeCount
        || transition.retainedControlVolumeCount
                + transition.appearedControlVolumeCount
            != transition.currentControlVolumeCount
        || transition.retainedControls.size()
            != transition.retainedControlVolumeCount
        || transition.retirementRecipients.size()
            != transition.disappearedControlVolumeCount
        || (transition.appearedControlVolumeCount != 0
            && transition.appearanceDonors.empty())
        || transition.ownedStorageBytes != storageBytes(transition)
        || !std::isfinite(
            transition.totalAppearanceDonorLinkAreaSquareMeters)
        || !std::isfinite(
            transition.totalRetirementRecipientLinkAreaSquareMeters)
        || transition.fingerprint != transitionFingerprint(transition)) {
        throw std::invalid_argument(
            "scene fluid pressure topology-transition integrity is invalid");
    }
}

void validateSceneFluidPressureTopologyTransition(
    const SceneFluidPressureTopologyTransition& transition,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks) {
    validateSources(
        previousPressureVolumes, previousFaceLinks,
        currentPressureVolumes, currentFaceLinks);
    validateSceneFluidPressureTopologyTransitionIntegrity(transition);
    const SceneFluidPressureTopologyTransitionLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildTransition(
        previousPressureVolumes, previousFaceLinks,
        currentPressureVolumes, currentFaceLinks, unlimited);
    if (transition != expected) {
        throw std::invalid_argument(
            "scene fluid pressure topology-transition payload is invalid");
    }
}

} // namespace simwing::fsi
