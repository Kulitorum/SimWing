#include "fluid/planar_region_fragment_opening_continuation_momentum_audit.h"

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

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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
            "opening continuation momentum-audit storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening continuation momentum-audit storage overflows");
    }
    return first * second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits) {
    if (limits.maximumSamples == 0 || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening continuation momentum-audit limits are invalid");
    }
}

double& coordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return value.x;
    case GridFaceAxis::Y:
        return value.y;
    case GridFaceAxis::Z:
        return value.z;
    }
    throw std::invalid_argument(
        "opening continuation momentum-audit axis is invalid");
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit&
        audit) {
    return checkedMultiply(
        audit.samples.size(),
        sizeof(
            PlanarPressureRegionFragmentOpeningContinuationMomentumSample));
}

std::uint64_t auditFingerprint(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit&
        audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    for (const std::uint64_t value : {
             audit.sourceContinuationFingerprint,
             audit.sourceAcceptedStateFingerprint,
             audit.previousTopologyFingerprint,
             audit.previousOpeningFingerprint,
             audit.currentTopologyFingerprint,
             audit.currentOpeningFingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.real(audit.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(audit.samples.size()));
    for (const auto& sample : audit.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(static_cast<std::uint8_t>(sample.kind));
        fingerprint.integer(sample.stableId);
        fingerprint.integer(static_cast<std::uint8_t>(sample.axis));
        fingerprint.integer(
            static_cast<std::uint64_t>(sample.previousSourceIndex));
        fingerprint.integer(
            static_cast<std::uint64_t>(sample.currentSourceIndex));
        for (const double value : {
                 sample.previousDualVolumeCubicMeters,
                 sample.currentDualVolumeCubicMeters,
                 sample.carriedVelocityMetersPerSecond,
                 sample.previousMassKilograms,
                 sample.currentMassKilograms,
                 sample.previousMomentumKilogramMetersPerSecond,
                 sample.carriedMomentumKilogramMetersPerSecond,
                 sample.momentumChangeKilogramMetersPerSecond,
                 sample.previousKineticEnergyJoules,
                 sample.carriedKineticEnergyJoules,
                 sample.kineticEnergyChangeJoules}) {
            fingerprint.real(value);
        }
        fingerprint.boolean(sample.metricChanged);
        fingerprint.boolean(sample.momentumChanged);
    }
    for (const std::size_t value : {
             audit.sameRegionGridSampleCount,
             audit.openingPatchSampleCount,
             audit.metricChangedSampleCount,
             audit.momentumChangedSampleCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprintVector(fingerprint, audit.previousMassByAxisKilograms);
    fingerprintVector(fingerprint, audit.currentMassByAxisKilograms);
    fingerprintVector(
        fingerprint, audit.previousMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, audit.carriedMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, audit.momentumChangeKilogramMetersPerSecond);
    for (const double value : {
             audit.previousKineticEnergyJoules,
             audit.carriedKineticEnergyJoules,
             audit.kineticEnergyChangeJoules,
             audit.maximumAbsoluteDualVolumeChangeCubicMeters,
             audit.maximumRelativeDualVolumeChange,
             audit.maximumAbsoluteMomentumChangeKilogramMetersPerSecond,
             audit.maximumAbsoluteKineticEnergyChangeJoules}) {
        fingerprint.real(value);
    }
    fingerprint.boolean(audit.metricChanged);
    fingerprint.boolean(audit.warmCarryChangesMomentum);
    fingerprint.boolean(audit.audited);
    fingerprint.integer(
        static_cast<std::uint64_t>(audit.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(audit.workingStorageBytes));
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

PlanarPressureRegionFragmentOpeningContinuationMomentumSample makeSample(
    const std::size_t sampleIndex,
    const PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind kind,
    const std::uint64_t stableId,
    const GridFaceAxis axis,
    const std::size_t previousSourceIndex,
    const std::size_t currentSourceIndex,
    const double previousDualVolumeCubicMeters,
    const double currentDualVolumeCubicMeters,
    const double carriedVelocityMetersPerSecond,
    const double densityKgPerCubicMeter) {
    PlanarPressureRegionFragmentOpeningContinuationMomentumSample sample;
    sample.sampleIndex = sampleIndex;
    sample.kind = kind;
    sample.stableId = stableId;
    sample.axis = axis;
    sample.previousSourceIndex = previousSourceIndex;
    sample.currentSourceIndex = currentSourceIndex;
    sample.previousDualVolumeCubicMeters =
        previousDualVolumeCubicMeters;
    sample.currentDualVolumeCubicMeters = currentDualVolumeCubicMeters;
    sample.carriedVelocityMetersPerSecond =
        carriedVelocityMetersPerSecond;
    sample.previousMassKilograms =
        densityKgPerCubicMeter * previousDualVolumeCubicMeters;
    sample.currentMassKilograms =
        densityKgPerCubicMeter * currentDualVolumeCubicMeters;
    sample.previousMomentumKilogramMetersPerSecond =
        sample.previousMassKilograms * carriedVelocityMetersPerSecond;
    sample.carriedMomentumKilogramMetersPerSecond =
        sample.currentMassKilograms * carriedVelocityMetersPerSecond;
    sample.momentumChangeKilogramMetersPerSecond =
        sample.carriedMomentumKilogramMetersPerSecond
        - sample.previousMomentumKilogramMetersPerSecond;
    const double velocitySquared = carriedVelocityMetersPerSecond
        * carriedVelocityMetersPerSecond;
    sample.previousKineticEnergyJoules =
        0.5 * sample.previousMassKilograms * velocitySquared;
    sample.carriedKineticEnergyJoules =
        0.5 * sample.currentMassKilograms * velocitySquared;
    sample.kineticEnergyChangeJoules =
        sample.carriedKineticEnergyJoules
        - sample.previousKineticEnergyJoules;
    sample.metricChanged = previousDualVolumeCubicMeters
        != currentDualVolumeCubicMeters;
    sample.momentumChanged =
        sample.momentumChangeKilogramMetersPerSecond != 0.0;
    return sample;
}

void clearSummary(
    PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit) {
    audit.sameRegionGridSampleCount = 0;
    audit.openingPatchSampleCount = 0;
    audit.metricChangedSampleCount = 0;
    audit.momentumChangedSampleCount = 0;
    audit.previousMassByAxisKilograms = {};
    audit.currentMassByAxisKilograms = {};
    audit.previousMomentumKilogramMetersPerSecond = {};
    audit.carriedMomentumKilogramMetersPerSecond = {};
    audit.momentumChangeKilogramMetersPerSecond = {};
    audit.previousKineticEnergyJoules = 0.0;
    audit.carriedKineticEnergyJoules = 0.0;
    audit.kineticEnergyChangeJoules = 0.0;
    audit.maximumAbsoluteDualVolumeChangeCubicMeters = 0.0;
    audit.maximumRelativeDualVolumeChange = 0.0;
    audit.maximumAbsoluteMomentumChangeKilogramMetersPerSecond = 0.0;
    audit.maximumAbsoluteKineticEnergyChangeJoules = 0.0;
    audit.metricChanged = false;
    audit.warmCarryChangesMomentum = false;
    audit.audited = false;
    audit.ownedStorageBytes = 0;
    audit.fingerprint = 0;
}

void summarize(
    PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit) {
    const auto workingStorageBytes = audit.workingStorageBytes;
    clearSummary(audit);
    audit.workingStorageBytes = workingStorageBytes;
    std::uint8_t previousKind = 0;
    std::uint64_t previousStableId = 0;
    for (std::size_t index = 0; index < audit.samples.size(); ++index) {
        const auto& sample = audit.samples[index];
        const auto kind = static_cast<std::uint8_t>(sample.kind);
        const double expectedPreviousMass = audit.densityKgPerCubicMeter
            * sample.previousDualVolumeCubicMeters;
        const double expectedCurrentMass = audit.densityKgPerCubicMeter
            * sample.currentDualVolumeCubicMeters;
        const double expectedPreviousMomentum = expectedPreviousMass
            * sample.carriedVelocityMetersPerSecond;
        const double expectedCurrentMomentum = expectedCurrentMass
            * sample.carriedVelocityMetersPerSecond;
        const double expectedMomentumChange =
            expectedCurrentMomentum - expectedPreviousMomentum;
        const double velocitySquared =
            sample.carriedVelocityMetersPerSecond
            * sample.carriedVelocityMetersPerSecond;
        const double expectedPreviousEnergy =
            0.5 * expectedPreviousMass * velocitySquared;
        const double expectedCurrentEnergy =
            0.5 * expectedCurrentMass * velocitySquared;
        const double expectedEnergyChange =
            expectedCurrentEnergy - expectedPreviousEnergy;
        const bool expectedMetricChanged =
            sample.previousDualVolumeCubicMeters
            != sample.currentDualVolumeCubicMeters;
        const bool expectedMomentumChanged =
            expectedMomentumChange != 0.0;
        if (sample.sampleIndex != index || sample.stableId == 0
            || (sample.kind
                    != PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
                        SameRegionGrid
                && sample.kind
                    != PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
                        OpeningPatch)
            || (index != 0
                && (kind < previousKind
                    || (kind == previousKind
                        && sample.stableId <= previousStableId)))
            || !std::isfinite(sample.previousDualVolumeCubicMeters)
            || !(sample.previousDualVolumeCubicMeters > 0.0)
            || !std::isfinite(sample.currentDualVolumeCubicMeters)
            || !(sample.currentDualVolumeCubicMeters > 0.0)
            || !std::isfinite(sample.carriedVelocityMetersPerSecond)
            || !std::isfinite(expectedPreviousMass)
            || !(expectedPreviousMass > 0.0)
            || !std::isfinite(expectedCurrentMass)
            || !(expectedCurrentMass > 0.0)
            || !std::isfinite(expectedPreviousMomentum)
            || !std::isfinite(expectedCurrentMomentum)
            || !std::isfinite(expectedMomentumChange)
            || !std::isfinite(expectedPreviousEnergy)
            || !std::isfinite(expectedCurrentEnergy)
            || !std::isfinite(expectedEnergyChange)
            || sample.previousMassKilograms != expectedPreviousMass
            || sample.currentMassKilograms != expectedCurrentMass
            || sample.previousMomentumKilogramMetersPerSecond
                != expectedPreviousMomentum
            || sample.carriedMomentumKilogramMetersPerSecond
                != expectedCurrentMomentum
            || sample.momentumChangeKilogramMetersPerSecond
                != expectedMomentumChange
            || sample.previousKineticEnergyJoules
                != expectedPreviousEnergy
            || sample.carriedKineticEnergyJoules
                != expectedCurrentEnergy
            || sample.kineticEnergyChangeJoules != expectedEnergyChange
            || sample.metricChanged != expectedMetricChanged
            || sample.momentumChanged != expectedMomentumChanged) {
            throw std::invalid_argument(
                "opening continuation momentum sample is invalid");
        }
        previousKind = kind;
        previousStableId = sample.stableId;
        if (sample.kind
            == PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
                SameRegionGrid) {
            ++audit.sameRegionGridSampleCount;
        } else {
            ++audit.openingPatchSampleCount;
        }
        audit.metricChangedSampleCount += sample.metricChanged ? 1 : 0;
        audit.momentumChangedSampleCount +=
            sample.momentumChanged ? 1 : 0;
        coordinate(audit.previousMassByAxisKilograms, sample.axis) +=
            sample.previousMassKilograms;
        coordinate(audit.currentMassByAxisKilograms, sample.axis) +=
            sample.currentMassKilograms;
        coordinate(
            audit.previousMomentumKilogramMetersPerSecond, sample.axis) +=
            sample.previousMomentumKilogramMetersPerSecond;
        coordinate(
            audit.carriedMomentumKilogramMetersPerSecond, sample.axis) +=
            sample.carriedMomentumKilogramMetersPerSecond;
        audit.previousKineticEnergyJoules +=
            sample.previousKineticEnergyJoules;
        audit.carriedKineticEnergyJoules +=
            sample.carriedKineticEnergyJoules;
        const double volumeChange = std::abs(
            sample.currentDualVolumeCubicMeters
            - sample.previousDualVolumeCubicMeters);
        const double volumeScale = std::max(
            sample.previousDualVolumeCubicMeters,
            sample.currentDualVolumeCubicMeters);
        audit.maximumAbsoluteDualVolumeChangeCubicMeters = std::max(
            audit.maximumAbsoluteDualVolumeChangeCubicMeters,
            volumeChange);
        audit.maximumRelativeDualVolumeChange = std::max(
            audit.maximumRelativeDualVolumeChange,
            volumeChange / volumeScale);
        audit.maximumAbsoluteMomentumChangeKilogramMetersPerSecond =
            std::max(
                audit.maximumAbsoluteMomentumChangeKilogramMetersPerSecond,
                std::abs(sample.momentumChangeKilogramMetersPerSecond));
        audit.maximumAbsoluteKineticEnergyChangeJoules = std::max(
            audit.maximumAbsoluteKineticEnergyChangeJoules,
            std::abs(sample.kineticEnergyChangeJoules));
    }
    audit.momentumChangeKilogramMetersPerSecond = {
        audit.carriedMomentumKilogramMetersPerSecond.x
            - audit.previousMomentumKilogramMetersPerSecond.x,
        audit.carriedMomentumKilogramMetersPerSecond.y
            - audit.previousMomentumKilogramMetersPerSecond.y,
        audit.carriedMomentumKilogramMetersPerSecond.z
            - audit.previousMomentumKilogramMetersPerSecond.z,
    };
    audit.kineticEnergyChangeJoules =
        audit.carriedKineticEnergyJoules
        - audit.previousKineticEnergyJoules;
    audit.metricChanged = audit.metricChangedSampleCount != 0;
    audit.warmCarryChangesMomentum =
        audit.momentumChangedSampleCount != 0;
    audit.audited = true;
    audit.ownedStorageBytes = ownedStorageBytes(audit);
    audit.fingerprint = auditFingerprint(audit);
}

PlanarPressureRegionFragmentOpeningContinuationMomentumAudit buildAudit(
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
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningContinuation(
        continuation, previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits.continuation);

    const std::size_t sampleCount = checkedAdd(
        previousTopology.sameRegionGridLinkCount,
        previousOpenings.patches.size());
    if (sampleCount > limits.maximumSamples) {
        throw std::length_error(
            "opening continuation momentum-audit sample limit exceeded");
    }
    const std::size_t mappingCount = checkedAdd(
        checkedAdd(previousTopology.links.size(), currentTopology.links.size()),
        checkedAdd(previousOpenings.patches.size(), currentOpenings.patches.size()));
    const std::size_t workingBytes =
        checkedMultiply(mappingCount, sizeof(std::size_t));
    const std::size_t ownedBytes = checkedMultiply(
        sampleCount,
        sizeof(
            PlanarPressureRegionFragmentOpeningContinuationMomentumSample));
    if (workingBytes > limits.maximumWorkingBytes
        || ownedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening continuation momentum-audit byte limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningContinuationMomentumAudit result;
    result.sourceContinuationFingerprint = continuation.fingerprint;
    result.sourceAcceptedStateFingerprint = previousState.fingerprint;
    result.previousTopologyFingerprint = previousTopology.fingerprint;
    result.previousOpeningFingerprint = previousOpenings.fingerprint;
    result.currentTopologyFingerprint = currentTopology.fingerprint;
    result.currentOpeningFingerprint = currentOpenings.fingerprint;
    result.densityKgPerCubicMeter =
        previousState.settings.projection.densityKgPerCubicMeter;
    result.workingStorageBytes = workingBytes;
    result.samples.reserve(sampleCount);

    const auto previousLinkOrder = stableOrder(
        previousTopology.links,
        [](const auto& link) { return link.stableId; });
    const auto currentLinkOrder = stableOrder(
        currentTopology.links,
        [](const auto& link) { return link.stableId; });
    for (std::size_t offset = 0; offset < currentLinkOrder.size(); ++offset) {
        const auto& previous =
            previousTopology.links[previousLinkOrder[offset]];
        const auto& current = currentTopology.links[currentLinkOrder[offset]];
        if (current.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        result.samples.push_back(makeSample(
            result.samples.size(),
            PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
                SameRegionGrid,
            current.stableId, current.axis, previous.linkIndex,
            current.linkIndex,
            previous.areaSquareMeters * previous.centerDistanceMeters,
            current.areaSquareMeters * current.centerDistanceMeters,
            continuation.orientedTopologyLinkVelocityMetersPerSecond[
                current.linkIndex],
            result.densityKgPerCubicMeter));
    }

    const auto previousPatchOrder = stableOrder(
        previousOpenings.patches,
        [](const auto& patch) { return patch.patchStableId; });
    const auto currentPatchOrder = stableOrder(
        currentOpenings.patches,
        [](const auto& patch) { return patch.patchStableId; });
    for (std::size_t offset = 0; offset < currentPatchOrder.size(); ++offset) {
        const auto& previous =
            previousOpenings.patches[previousPatchOrder[offset]];
        const auto& current =
            currentOpenings.patches[currentPatchOrder[offset]];
        result.samples.push_back(makeSample(
            result.samples.size(),
            PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
                OpeningPatch,
            current.patchStableId, current.axis, previous.patchIndex,
            current.patchIndex,
            previous.areaSquareMeters * previous.centerDistanceMeters,
            current.areaSquareMeters * current.centerDistanceMeters,
            continuation.openingVelocitySamples[offset]
                .relativeNormalVelocityMetersPerSecond,
            result.densityKgPerCubicMeter));
    }
    summarize(result);
    if (result.ownedStorageBytes != ownedBytes) {
        throw std::logic_error(
            "opening continuation momentum-audit storage changed");
    }
    validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningContinuationMomentumAudit
auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
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
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits) {
    return buildAudit(
        continuation, previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);
}

void validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit) {
    if (audit.version
            != planarPressureRegionFragmentOpeningContinuationMomentumAuditVersion
        || audit.fingerprint == 0
        || audit.sourceContinuationFingerprint == 0
        || audit.sourceAcceptedStateFingerprint == 0
        || audit.previousTopologyFingerprint == 0
        || audit.previousOpeningFingerprint == 0
        || audit.currentTopologyFingerprint == 0
        || audit.currentOpeningFingerprint == 0
        || !std::isfinite(audit.densityKgPerCubicMeter)
        || !(audit.densityKgPerCubicMeter > 0.0)
        || !finiteVector(audit.previousMassByAxisKilograms)
        || !finiteVector(audit.currentMassByAxisKilograms)
        || !finiteVector(audit.previousMomentumKilogramMetersPerSecond)
        || !finiteVector(audit.carriedMomentumKilogramMetersPerSecond)
        || !finiteVector(audit.momentumChangeKilogramMetersPerSecond)
        || !std::isfinite(audit.previousKineticEnergyJoules)
        || !std::isfinite(audit.carriedKineticEnergyJoules)
        || !std::isfinite(audit.kineticEnergyChangeJoules)
        || !std::isfinite(
            audit.maximumAbsoluteDualVolumeChangeCubicMeters)
        || !std::isfinite(audit.maximumRelativeDualVolumeChange)
        || !std::isfinite(
            audit.maximumAbsoluteMomentumChangeKilogramMetersPerSecond)
        || !std::isfinite(
            audit.maximumAbsoluteKineticEnergyChangeJoules)
        || !audit.audited || audit.workingStorageBytes == 0) {
        throw std::invalid_argument(
            "opening continuation momentum-audit integrity is invalid");
    }
    auto expected = audit;
    summarize(expected);
    if (expected != audit) {
        throw std::invalid_argument(
            "opening continuation momentum-audit summary is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit,
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
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
        audit);
    if (audit.samples.size() > limits.maximumSamples
        || audit.ownedStorageBytes > limits.maximumOwnedBytes
        || audit.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening continuation momentum-audit validation limit exceeded");
    }
    const auto expected = buildAudit(
        continuation, previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, limits);
    if (audit != expected) {
        throw std::invalid_argument(
            "opening continuation momentum audit is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
