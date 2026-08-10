#include "fluid/planar_region_fragment_opening_continuation.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
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
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening continuation storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening continuation storage overflows");
    }
    return first + second;
}

template<typename... Values>
std::size_t storageSum(const Values... values) {
    std::size_t result = 0;
    ((result = checkedAdd(result, values)), ...);
    return result;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits) {
    if (limits.maximumTopologyLinkVelocities == 0
        || limits.maximumOpeningSamples == 0
        || limits.maximumPressureCorrections == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening continuation limits are invalid");
    }
}

bool finite(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation) {
    return storageSum(
        checkedMultiply(
            continuation.componentGaugeShiftsPascals.size(),
            sizeof(double)),
        checkedMultiply(
            continuation.orientedTopologyLinkVelocityMetersPerSecond.size(),
            sizeof(double)),
        checkedMultiply(
            continuation.openingVelocitySamples.size(),
            sizeof(PlanarPressureRegionFragmentOpeningVelocitySample)),
        continuation.openingFlux.ownedStorageBytes,
        checkedMultiply(
            continuation.pressureCorrectionPascals.size(),
            sizeof(double)));
}

std::size_t workingStorageBytes(const std::size_t previousLinkCount,
                                const std::size_t currentLinkCount,
                                const std::size_t previousPatchCount,
                                const std::size_t currentPatchCount,
                                const std::size_t previousRowCount,
                                const std::size_t currentRowCount,
                                const std::size_t previousComponentCount,
                                const std::size_t currentComponentCount) {
    return checkedMultiply(
        storageSum(previousLinkCount, currentLinkCount,
                   previousPatchCount, currentPatchCount,
                   previousRowCount, currentRowCount,
                   previousComponentCount, currentComponentCount),
        sizeof(std::size_t));
}

std::uint64_t continuationFingerprint(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation) {
    Fingerprint fingerprint;
    fingerprint.integer(continuation.version);
    for (const std::uint64_t value : {
             continuation.sourceAcceptedStateFingerprint,
             continuation.previousPressureOperatorFingerprint,
             continuation.previousBasePressureOperatorFingerprint,
             continuation.previousOpeningFingerprint,
             continuation.previousFragmentFingerprint,
             continuation.previousTopologyFingerprint,
             continuation.previousVolumeRateFingerprint,
             continuation.currentPressureOperatorFingerprint,
             continuation.currentBasePressureOperatorFingerprint,
             continuation.currentOpeningFingerprint,
             continuation.currentFragmentFingerprint,
             continuation.currentTopologyFingerprint,
             continuation.currentVolumeRateFingerprint,
             continuation.currentOpeningFluxFingerprint}) {
        fingerprint.integer(value);
    }
    for (const std::size_t value : {
             continuation.topologyLinkCount,
             continuation.openingPatchCount,
             continuation.pressureCorrectionCount,
             continuation.connectedComponentCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(continuation.maximumAbsoluteGaugeShiftPascals);
    fingerprint.real(
        continuation.maximumAbsolutePressureCorrectionPascals);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.componentGaugeShiftsPascals.size()));
    for (const double value : continuation.componentGaugeShiftsPascals)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.orientedTopologyLinkVelocityMetersPerSecond.size()));
    for (const double value
         : continuation.orientedTopologyLinkVelocityMetersPerSecond)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.openingVelocitySamples.size()));
    for (const auto& sample : continuation.openingVelocitySamples) {
        fingerprint.integer(sample.patchStableId);
        fingerprint.real(sample.relativeNormalVelocityMetersPerSecond);
    }
    fingerprint.integer(continuation.openingFlux.fingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.pressureCorrectionPascals.size()));
    for (const double value : continuation.pressureCorrectionPascals)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.workingStorageBytes));
    return fingerprint.value();
}

template<typename Values, typename StableId>
std::vector<std::size_t> stableOrder(const Values& values,
                                     StableId stableId) {
    std::vector<std::size_t> result(values.size());
    std::iota(result.begin(), result.end(), 0);
    std::ranges::sort(result, [&](const std::size_t first,
                                 const std::size_t second) {
        return stableId(values[first]) < stableId(values[second]);
    });
    return result;
}

bool sameLinkIdentity(
    const PlanarPressureRegionFragmentFaceLink& previous,
    const PlanarPressureRegionFragmentFaceLink& current) {
    return previous.stableId == current.stableId
        && previous.kind == current.kind
        && previous.axis == current.axis
        && previous.i == current.i
        && previous.j == current.j
        && previous.k == current.k
        && previous.facePeriodicImage == current.facePeriodicImage
        && previous.surfaceStableId == current.surfaceStableId
        && previous.minusFragmentStableId
            == current.minusFragmentStableId
        && previous.plusFragmentStableId
            == current.plusFragmentStableId
        && previous.minusRegionStableId == current.minusRegionStableId
        && previous.plusRegionStableId == current.plusRegionStableId
        && previous.pressureJumpPascals == current.pressureJumpPascals
        && previous.unitNormalMinusToPlus
            == current.unitNormalMinusToPlus
        && previous.crossesPeriodicBoundary
            == current.crossesPeriodicBoundary;
}

bool samePatchIdentity(
    const PlanarPressureRegionFragmentOpeningPatch& previous,
    const PlanarPressureRegionFragmentOpeningPatch& current) {
    return previous.patchStableId == current.patchStableId
        && previous.openingStableId == current.openingStableId
        && previous.surfaceStableId == current.surfaceStableId
        && previous.axis == current.axis
        && previous.i == current.i
        && previous.j == current.j
        && previous.k == current.k
        && previous.sourceFaceLinkStableId
            == current.sourceFaceLinkStableId
        && previous.minusFragmentStableId
            == current.minusFragmentStableId
        && previous.plusFragmentStableId
            == current.plusFragmentStableId
        && previous.negativeSideRegionStableId
            == current.negativeSideRegionStableId
        && previous.positiveSideRegionStableId
            == current.positiveSideRegionStableId
        && previous.areaSquareMeters == current.areaSquareMeters
        && previous.sourceWallAreaSquareMeters
            == current.sourceWallAreaSquareMeters
        && previous.sourceWallAreaFraction
            == current.sourceWallAreaFraction
        && previous.usesAuthoredCentroid
            == current.usesAuthoredCentroid
        && previous.unitNormalNegativeToPositive
            == current.unitNormalNegativeToPositive;
}

void validateSources(
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, limits.previousAcceptedStateLimits);
    validatePlanarPressureRegionFragmentOpeningPressureOperator(
        currentPressureOperator, currentBasePressureOperator, grid,
        currentSweep, currentFragments, currentTopology,
        currentOpeningDefinitions, currentOpenings,
        limits.currentPressureOperatorLimits);
    validatePlanarPressureRegionFragmentVolumeRates(
        currentVolumeRates, grid, currentSweep, currentFragments,
        currentTopology, limits.currentVolumeRateLimits);
    if (previousSweep.currentProfile != currentSweep.previousProfile) {
        throw std::invalid_argument(
            "opening continuation epochs are not consecutive");
    }
}

PlanarPressureRegionFragmentOpeningContinuation buildContinuation(
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits) {
    validateSources(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);

    if (currentTopology.links.size()
            > limits.maximumTopologyLinkVelocities
        || currentOpenings.patches.size() > limits.maximumOpeningSamples
        || currentPressureOperator.rows.size()
            > limits.maximumPressureCorrections
        || currentPressureOperator.components.size()
            > limits.maximumComponents) {
        throw std::length_error(
            "opening continuation count limit exceeded");
    }
    const std::size_t workingBytes = workingStorageBytes(
        previousTopology.links.size(), currentTopology.links.size(),
        previousOpenings.patches.size(), currentOpenings.patches.size(),
        previousPressureOperator.rows.size(),
        currentPressureOperator.rows.size(),
        previousPressureOperator.components.size(),
        currentPressureOperator.components.size());
    if (workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening continuation working-storage limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningContinuation result;
    result.sourceAcceptedStateFingerprint = previousState.fingerprint;
    result.previousPressureOperatorFingerprint =
        previousPressureOperator.fingerprint;
    result.previousBasePressureOperatorFingerprint =
        previousBasePressureOperator.fingerprint;
    result.previousOpeningFingerprint = previousOpenings.fingerprint;
    result.previousFragmentFingerprint = previousFragments.fingerprint;
    result.previousTopologyFingerprint = previousTopology.fingerprint;
    result.previousVolumeRateFingerprint = previousVolumeRates.fingerprint;
    result.currentPressureOperatorFingerprint =
        currentPressureOperator.fingerprint;
    result.currentBasePressureOperatorFingerprint =
        currentBasePressureOperator.fingerprint;
    result.currentOpeningFingerprint = currentOpenings.fingerprint;
    result.currentFragmentFingerprint = currentFragments.fingerprint;
    result.currentTopologyFingerprint = currentTopology.fingerprint;
    result.currentVolumeRateFingerprint = currentVolumeRates.fingerprint;
    result.topologyLinkCount = currentTopology.links.size();
    result.openingPatchCount = currentOpenings.patches.size();
    result.pressureCorrectionCount = currentPressureOperator.rows.size();
    result.connectedComponentCount =
        currentPressureOperator.components.size();
    result.workingStorageBytes = workingBytes;
    result.orientedTopologyLinkVelocityMetersPerSecond.assign(
        result.topologyLinkCount, 0.0);
    result.pressureCorrectionPascals.assign(
        result.pressureCorrectionCount, 0.0);

    const auto previousLinkOrder = stableOrder(
        previousTopology.links,
        [](const auto& link) { return link.stableId; });
    const auto currentLinkOrder = stableOrder(
        currentTopology.links,
        [](const auto& link) { return link.stableId; });
    if (previousLinkOrder.size() != currentLinkOrder.size()) {
        throw std::invalid_argument(
            "opening continuation link coverage changed");
    }
    for (std::size_t offset = 0; offset < currentLinkOrder.size(); ++offset) {
        const auto& previous =
            previousTopology.links[previousLinkOrder[offset]];
        const auto& current = currentTopology.links[currentLinkOrder[offset]];
        if (!sameLinkIdentity(previous, current)) {
            throw std::invalid_argument(
                "opening continuation link identity changed");
        }
        result.orientedTopologyLinkVelocityMetersPerSecond[
            current.linkIndex] = previousState
                .orientedTopologyLinkVelocityMetersPerSecond[
                    previous.linkIndex];
    }

    const auto previousPatchOrder = stableOrder(
        previousOpenings.patches,
        [](const auto& patch) { return patch.patchStableId; });
    const auto currentPatchOrder = stableOrder(
        currentOpenings.patches,
        [](const auto& patch) { return patch.patchStableId; });
    if (previousPatchOrder.size() != currentPatchOrder.size()) {
        throw std::invalid_argument(
            "opening continuation aperture coverage changed");
    }
    result.openingVelocitySamples.reserve(currentPatchOrder.size());
    for (std::size_t offset = 0; offset < currentPatchOrder.size(); ++offset) {
        const auto& previous =
            previousOpenings.patches[previousPatchOrder[offset]];
        const auto& current = currentOpenings.patches[currentPatchOrder[offset]];
        if (!samePatchIdentity(previous, current)
            || previousState.openingVelocitySamples[offset].patchStableId
                != previous.patchStableId) {
            throw std::invalid_argument(
                "opening continuation aperture identity changed");
        }
        result.openingVelocitySamples.push_back({
            current.patchStableId,
            previousState.openingVelocitySamples[offset]
                .relativeNormalVelocityMetersPerSecond,
        });
    }

    const auto previousRowOrder = stableOrder(
        previousPressureOperator.rows,
        [](const auto& row) { return row.fragmentStableId; });
    const auto currentRowOrder = stableOrder(
        currentPressureOperator.rows,
        [](const auto& row) { return row.fragmentStableId; });
    if (previousRowOrder.size() != currentRowOrder.size()) {
        throw std::invalid_argument(
            "opening continuation fragment coverage changed");
    }
    for (std::size_t offset = 0; offset < currentRowOrder.size(); ++offset) {
        const auto& previous =
            previousPressureOperator.rows[previousRowOrder[offset]];
        const auto& current =
            currentPressureOperator.rows[currentRowOrder[offset]];
        const auto& previousComponent = previousPressureOperator.components[
            previous.connectedComponentIndex];
        const auto& currentComponent = currentPressureOperator.components[
            current.connectedComponentIndex];
        if (previous.fragmentStableId != current.fragmentStableId
            || previous.stableId != current.stableId
            || previousTopology.fragments[previous.fragmentIndex]
                    .regionStableId
                != currentTopology.fragments[current.fragmentIndex]
                       .regionStableId
            || previousComponent.stableId != currentComponent.stableId) {
            throw std::invalid_argument(
                "opening continuation pressure identity changed");
        }
        result.pressureCorrectionPascals[current.rowIndex] =
            previousState.pressureCorrectionPascals[previous.rowIndex];
    }

    const auto previousComponentOrder = stableOrder(
        previousPressureOperator.components,
        [](const auto& component) { return component.stableId; });
    const auto currentComponentOrder = stableOrder(
        currentPressureOperator.components,
        [](const auto& component) { return component.stableId; });
    if (previousComponentOrder.size() != currentComponentOrder.size()) {
        throw std::invalid_argument(
            "opening continuation pressure components changed");
    }
    for (std::size_t offset = 0;
         offset < currentComponentOrder.size(); ++offset) {
        const auto& previous = previousPressureOperator.components[
            previousComponentOrder[offset]];
        const auto& current = currentPressureOperator.components[
            currentComponentOrder[offset]];
        if (previous.stableId != current.stableId
            || previous.baseComponentCount != current.baseComponentCount
            || previous.fragmentCount != current.fragmentCount) {
            throw std::invalid_argument(
                "opening continuation pressure component identity changed");
        }
    }

    result.componentGaugeShiftsPascals.assign(
        result.connectedComponentCount, 0.0);
    for (const auto& component : currentPressureOperator.components) {
        double pressureMoment = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            const std::size_t fragmentIndex =
                currentPressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset];
            pressureMoment +=
                currentFragments.fragments[fragmentIndex].volumeCubicMeters
                * result.pressureCorrectionPascals[fragmentIndex];
        }
        const double shift = pressureMoment
            / component.totalVolumeCubicMeters;
        if (!std::isfinite(shift)) {
            throw std::overflow_error(
                "opening continuation pressure gauge is non-finite");
        }
        result.componentGaugeShiftsPascals[component.componentIndex] = shift;
        result.maximumAbsoluteGaugeShiftPascals = std::max(
            result.maximumAbsoluteGaugeShiftPascals, std::abs(shift));
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            const std::size_t fragmentIndex =
                currentPressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset];
            result.pressureCorrectionPascals[fragmentIndex] -= shift;
        }
    }
    for (const double pressure : result.pressureCorrectionPascals) {
        if (!std::isfinite(pressure)) {
            throw std::overflow_error(
                "opening continuation pressure is non-finite");
        }
        result.maximumAbsolutePressureCorrectionPascals = std::max(
            result.maximumAbsolutePressureCorrectionPascals,
            std::abs(pressure));
    }

    result.openingFlux =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            grid, currentSweep, currentFragments, currentTopology,
            currentOpeningDefinitions, currentOpenings,
            result.openingVelocitySamples, limits.currentOpeningFluxLimits);
    result.currentOpeningFluxFingerprint = result.openingFlux.fingerprint;
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening continuation owned-storage limit exceeded");
    }
    result.fingerprint = continuationFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningContinuation
buildPlanarPressureRegionFragmentOpeningContinuation(
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits) {
    return buildContinuation(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);
}

void validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation) {
    validatePlanarPressureRegionFragmentOpeningFluxStateIntegrity(
        continuation.openingFlux);
    double maximumGaugeShift = 0.0;
    for (const double value : continuation.componentGaugeShiftsPascals)
        maximumGaugeShift = std::max(maximumGaugeShift, std::abs(value));
    double maximumPressure = 0.0;
    for (const double value : continuation.pressureCorrectionPascals)
        maximumPressure = std::max(maximumPressure, std::abs(value));
    bool sampleIdentity = continuation.openingVelocitySamples.size()
        == continuation.openingFlux.patches.size();
    for (std::size_t index = 0;
         sampleIdentity
             && index < continuation.openingVelocitySamples.size();
         ++index) {
        const auto& sample = continuation.openingVelocitySamples[index];
        const auto& flux = continuation.openingFlux.patches[index];
        sampleIdentity = sample.patchStableId != 0
            && std::isfinite(
                sample.relativeNormalVelocityMetersPerSecond)
            && (index == 0
                || continuation.openingVelocitySamples[index - 1]
                       .patchStableId < sample.patchStableId)
            && sample.patchStableId == flux.patchStableId
            && sample.relativeNormalVelocityMetersPerSecond
                == flux.relativeNormalVelocityMetersPerSecond;
    }
    if (continuation.version
            != planarPressureRegionFragmentOpeningContinuationVersion
        || continuation.fingerprint == 0
        || continuation.sourceAcceptedStateFingerprint == 0
        || continuation.previousPressureOperatorFingerprint == 0
        || continuation.previousBasePressureOperatorFingerprint == 0
        || continuation.previousOpeningFingerprint == 0
        || continuation.previousFragmentFingerprint == 0
        || continuation.previousTopologyFingerprint == 0
        || continuation.previousVolumeRateFingerprint == 0
        || continuation.currentPressureOperatorFingerprint == 0
        || continuation.currentBasePressureOperatorFingerprint == 0
        || continuation.currentOpeningFingerprint == 0
        || continuation.currentFragmentFingerprint == 0
        || continuation.currentTopologyFingerprint == 0
        || continuation.currentVolumeRateFingerprint == 0
        || continuation.currentOpeningFluxFingerprint == 0
        || continuation.topologyLinkCount
            != continuation
                   .orientedTopologyLinkVelocityMetersPerSecond.size()
        || continuation.openingPatchCount
            != continuation.openingVelocitySamples.size()
        || continuation.pressureCorrectionCount
            != continuation.pressureCorrectionPascals.size()
        || continuation.connectedComponentCount
            != continuation.componentGaugeShiftsPascals.size()
        || continuation.currentOpeningFluxFingerprint
            != continuation.openingFlux.fingerprint
        || continuation.currentOpeningFingerprint
            != continuation.openingFlux.sourceOpeningFingerprint
        || continuation.currentFragmentFingerprint
            != continuation.openingFlux.sourceFragmentFingerprint
        || continuation.currentTopologyFingerprint
            != continuation.openingFlux.sourceTopologyFingerprint
        || !sampleIdentity
        || !finite(continuation.componentGaugeShiftsPascals)
        || !finite(
            continuation.orientedTopologyLinkVelocityMetersPerSecond)
        || !finite(continuation.pressureCorrectionPascals)
        || !std::isfinite(
            continuation.maximumAbsoluteGaugeShiftPascals)
        || !std::isfinite(
            continuation.maximumAbsolutePressureCorrectionPascals)
        || continuation.maximumAbsoluteGaugeShiftPascals
            != maximumGaugeShift
        || continuation.maximumAbsolutePressureCorrectionPascals
            != maximumPressure
        || continuation.ownedStorageBytes
            != ownedStorageBytes(continuation)
        || continuation.fingerprint
            != continuationFingerprint(continuation)) {
        throw std::invalid_argument(
            "opening continuation integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningContinuation(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation,
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(
        continuation);
    const auto expected = buildContinuation(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);
    if (continuation != expected) {
        throw std::invalid_argument(
            "opening continuation payload is invalid");
    }
}

} // namespace simwing::fsi::fluid
