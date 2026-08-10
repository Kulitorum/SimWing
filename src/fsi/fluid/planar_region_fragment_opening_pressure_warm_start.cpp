#include "fluid/planar_region_fragment_opening_pressure_warm_start.h"

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

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening pressure warm-start storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening pressure warm-start storage overflows");
    }
    return first * second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits) {
    if (limits.maximumPressureCorrections == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening pressure warm-start limits are invalid");
    }
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart) {
    return checkedAdd(
        checkedMultiply(
            warmStart.componentGaugeShiftsPascals.size(), sizeof(double)),
        checkedMultiply(
            warmStart.pressureCorrectionPascals.size(), sizeof(double)));
}

std::uint64_t warmStartFingerprint(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart) {
    Fingerprint fingerprint;
    fingerprint.integer(warmStart.version);
    for (const std::uint64_t value : {
             warmStart.sourceAcceptedStateFingerprint,
             warmStart.previousPressureOperatorFingerprint,
             warmStart.previousBasePressureOperatorFingerprint,
             warmStart.previousOpeningFingerprint,
             warmStart.previousFragmentFingerprint,
             warmStart.previousTopologyFingerprint,
             warmStart.previousVolumeRateFingerprint,
             warmStart.currentPressureOperatorFingerprint,
             warmStart.currentBasePressureOperatorFingerprint,
             warmStart.currentOpeningFingerprint,
             warmStart.currentFragmentFingerprint,
             warmStart.currentTopologyFingerprint,
             warmStart.currentVolumeRateFingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.pressureCorrectionCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.connectedComponentCount));
    fingerprint.real(warmStart.maximumAbsoluteGaugeShiftPascals);
    fingerprint.real(
        warmStart.maximumAbsolutePressureCorrectionPascals);
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.componentGaugeShiftsPascals.size()));
    for (const double value : warmStart.componentGaugeShiftsPascals)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.pressureCorrectionPascals.size()));
    for (const double value : warmStart.pressureCorrectionPascals)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.workingStorageBytes));
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

PlanarPressureRegionFragmentOpeningPressureWarmStart buildWarmStart(
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
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions,
        limits.previousAcceptedStateLimits);
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
            "opening pressure warm-start epochs are not consecutive");
    }

    const std::size_t correctionCount =
        currentPressureOperator.rows.size();
    const std::size_t componentCount =
        currentPressureOperator.components.size();
    if (correctionCount == 0 || componentCount == 0
        || correctionCount > limits.maximumPressureCorrections
        || componentCount > limits.maximumComponents) {
        throw std::length_error(
            "opening pressure warm-start count limit exceeded");
    }
    std::size_t workingBytes = checkedMultiply(
        checkedAdd(
            previousPressureOperator.rows.size(), correctionCount),
        sizeof(std::size_t));
    workingBytes = checkedAdd(
        workingBytes,
        checkedMultiply(
            checkedAdd(
                previousPressureOperator.components.size(),
                componentCount),
            sizeof(std::size_t)));
    if (workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening pressure warm-start working-storage limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningPressureWarmStart result;
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
    result.pressureCorrectionCount = correctionCount;
    result.connectedComponentCount = componentCount;
    result.workingStorageBytes = workingBytes;
    result.pressureCorrectionPascals.assign(correctionCount, 0.0);

    const auto previousRowOrder = stableOrder(
        previousPressureOperator.rows,
        [](const auto& row) { return row.fragmentStableId; });
    const auto currentRowOrder = stableOrder(
        currentPressureOperator.rows,
        [](const auto& row) { return row.fragmentStableId; });
    if (previousRowOrder.size() != currentRowOrder.size()) {
        throw std::invalid_argument(
            "opening pressure warm-start fragment coverage changed");
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
                "opening pressure warm-start fragment identity changed");
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
            "opening pressure warm-start component coverage changed");
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
                "opening pressure warm-start component identity changed");
        }
    }

    result.componentGaugeShiftsPascals.assign(componentCount, 0.0);
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
                "opening pressure warm-start gauge is non-finite");
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
                "opening pressure warm-start correction is non-finite");
        }
        result.maximumAbsolutePressureCorrectionPascals = std::max(
            result.maximumAbsolutePressureCorrectionPascals,
            std::abs(pressure));
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening pressure warm-start owned-storage limit exceeded");
    }
    result.fingerprint = warmStartFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureWarmStart
buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
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
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits) {
    return buildWarmStart(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);
}

void validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart) {
    bool valuesFinite = true;
    double maximumGaugeShift = 0.0;
    for (const double value : warmStart.componentGaugeShiftsPascals) {
        valuesFinite = valuesFinite && std::isfinite(value);
        maximumGaugeShift = std::max(maximumGaugeShift, std::abs(value));
    }
    double maximumPressureCorrection = 0.0;
    for (const double value : warmStart.pressureCorrectionPascals) {
        valuesFinite = valuesFinite && std::isfinite(value);
        maximumPressureCorrection = std::max(
            maximumPressureCorrection, std::abs(value));
    }
    if (warmStart.version
            != planarPressureRegionFragmentOpeningPressureWarmStartVersion
        || warmStart.fingerprint == 0
        || warmStart.sourceAcceptedStateFingerprint == 0
        || warmStart.previousPressureOperatorFingerprint == 0
        || warmStart.previousBasePressureOperatorFingerprint == 0
        || warmStart.previousOpeningFingerprint == 0
        || warmStart.previousFragmentFingerprint == 0
        || warmStart.previousTopologyFingerprint == 0
        || warmStart.previousVolumeRateFingerprint == 0
        || warmStart.currentPressureOperatorFingerprint == 0
        || warmStart.currentBasePressureOperatorFingerprint == 0
        || warmStart.currentOpeningFingerprint == 0
        || warmStart.currentFragmentFingerprint == 0
        || warmStart.currentTopologyFingerprint == 0
        || warmStart.currentVolumeRateFingerprint == 0
        || warmStart.pressureCorrectionCount == 0
        || warmStart.connectedComponentCount == 0
        || warmStart.pressureCorrectionPascals.size()
            != warmStart.pressureCorrectionCount
        || warmStart.componentGaugeShiftsPascals.size()
            != warmStart.connectedComponentCount
        || !std::isfinite(warmStart.maximumAbsoluteGaugeShiftPascals)
        || warmStart.maximumAbsoluteGaugeShiftPascals < 0.0
        || !std::isfinite(
            warmStart.maximumAbsolutePressureCorrectionPascals)
        || warmStart.maximumAbsolutePressureCorrectionPascals < 0.0
        || !valuesFinite
        || warmStart.maximumAbsoluteGaugeShiftPascals
            != maximumGaugeShift
        || warmStart.maximumAbsolutePressureCorrectionPascals
            != maximumPressureCorrection
        || warmStart.ownedStorageBytes != ownedStorageBytes(warmStart)
        || warmStart.workingStorageBytes == 0
        || warmStart.fingerprint != warmStartFingerprint(warmStart)) {
        throw std::invalid_argument(
            "opening pressure warm-start integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningPressureWarmStart(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart,
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
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
        warmStart);
    if (warmStart.pressureCorrectionCount
            > limits.maximumPressureCorrections
        || warmStart.connectedComponentCount > limits.maximumComponents
        || warmStart.ownedStorageBytes > limits.maximumOwnedBytes
        || warmStart.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening pressure warm-start validation limit exceeded");
    }
    if (warmStart != buildWarmStart(
            previousState, previousPressureOperator,
            previousBasePressureOperator, grid, previousSweep,
            previousFragments, previousTopology, previousVolumeRates,
            previousOpeningDefinitions, previousOpenings,
            previousResistanceDefinitions, currentPressureOperator,
            currentBasePressureOperator, currentSweep, currentFragments,
            currentTopology, currentVolumeRates,
            currentOpeningDefinitions, currentOpenings, limits)) {
        throw std::invalid_argument(
            "opening pressure warm-start is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
