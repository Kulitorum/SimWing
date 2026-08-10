#include "fluid/planar_region_fragment_opening_pressure_epoch.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

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
            "opening pressure-epoch storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening pressure-epoch storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningPressureEpochResult& result) {
    return checkedAdd(
        checkedAdd(
            checkedMultiply(
                result.diagnostics.pressureStep.resistance.patches.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningResistancePatchDiagnostics)),
            checkedMultiply(
                result.diagnostics.pressureStep.projection.pressureSolve
                    .components.size(),
                sizeof(
                    PlanarPressureRegionFragmentPressureSolveComponentDiagnostics))),
        result.acceptedState.ownedStorageBytes);
}

std::uint64_t resistanceDefinitionFingerprint(
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    if (definitions.size() != openings.patches.size()) {
        throw std::invalid_argument(
            "opening pressure-epoch resistance coverage is incomplete");
    }
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        canonical(definitions.begin(), definitions.end());
    std::ranges::sort(canonical, {},
        &PlanarPressureRegionFragmentOpeningResistanceDefinition::
            patchStableId);
    std::vector<std::uint64_t> patchIds;
    patchIds.reserve(openings.patches.size());
    for (const auto& patch : openings.patches)
        patchIds.push_back(patch.patchStableId);
    std::ranges::sort(patchIds);

    Fingerprint fingerprint;
    fingerprint.integer(static_cast<std::uint64_t>(canonical.size()));
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto& definition = canonical[index];
        const auto& resistance = definition.resistance;
        if (definition.patchStableId == 0
            || definition.patchStableId != patchIds[index]
            || (index != 0
                && canonical[index - 1].patchStableId
                    == definition.patchStableId)
            || !std::isfinite(resistance.linearPascalSecondsPerMeter)
            || resistance.linearPascalSecondsPerMeter < 0.0
            || !std::isfinite(
                resistance.quadraticPascalSecondsSquaredPerSquareMeter)
            || resistance.quadraticPascalSecondsSquaredPerSquareMeter
                < 0.0) {
            throw std::invalid_argument(
                "opening pressure-epoch resistance identity is invalid");
        }
        fingerprint.integer(definition.patchStableId);
        fingerprint.real(resistance.linearPascalSecondsPerMeter);
        fingerprint.real(
            resistance.quadraticPascalSecondsSquaredPerSquareMeter);
    }
    return fingerprint.value();
}

PlanarPressureRegionFragmentOpeningPressureEpochFailureStage failureStage(
    const PlanarPressureRegionFragmentOpeningPressureStepDiagnostics& step) {
    if (step.accepted) {
        return PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
            None;
    }
    if (!step.resistance.accepted) {
        return PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
            Resistance;
    }
    if (!step.projection.accepted) {
        return PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
            PressureProjection;
    }
    return PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
        AggregateEnergy;
}

bool emptyAcceptedState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state) {
    return state.fingerprint == 0 && !state.accepted
        && state.orientedTopologyLinkVelocityMetersPerSecond.empty()
        && state.openingVelocitySamples.empty()
        && state.openingFlux.fingerprint == 0
        && state.openingFlux.patches.empty()
        && state.openingFlux.fragments.empty()
        && state.openingFlux.baseComponents.empty()
        && state.openingFlux.connectedComponents.empty()
        && state.pressureCorrectionPascals.empty();
}

PlanarPressureRegionFragmentOpeningPressureEpochResult acceptEpoch(
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureEpochLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "opening pressure-epoch limits are invalid");
    }
    const auto continuation =
        buildPlanarPressureRegionFragmentOpeningContinuation(
            previousState, previousPressureOperator,
            previousBasePressureOperator, grid, previousSweep,
            previousFragments, previousTopology, previousVolumeRates,
            previousOpeningDefinitions, previousOpenings,
            previousResistanceDefinitions, currentPressureOperator,
            currentBasePressureOperator, currentSweep, currentFragments,
            currentTopology, currentVolumeRates,
            currentOpeningDefinitions, currentOpenings,
            limits.continuation);

    PlanarPressureRegionFragmentOpeningPressureEpochResult result;
    result.sourceAcceptedStateFingerprint = previousState.fingerprint;
    result.continuationFingerprint = continuation.fingerprint;
    result.continuationOpeningFluxFingerprint =
        continuation.currentOpeningFluxFingerprint;
    result.currentPressureOperatorFingerprint =
        currentPressureOperator.fingerprint;
    result.currentBasePressureOperatorFingerprint =
        currentBasePressureOperator.fingerprint;
    result.currentOpeningFingerprint = currentOpenings.fingerprint;
    result.currentFragmentFingerprint = currentFragments.fingerprint;
    result.currentTopologyFingerprint = currentTopology.fingerprint;
    result.currentVolumeRateFingerprint = currentVolumeRates.fingerprint;
    result.currentResistanceDefinitionFingerprint =
        resistanceDefinitionFingerprint(
            currentResistanceDefinitions, currentOpenings);
    result.settings = settings;
    result.diagnostics.usedConsecutiveContinuation = true;

    auto topologyVelocity =
        continuation.orientedTopologyLinkVelocityMetersPerSecond;
    auto openingSamples = continuation.openingVelocitySamples;
    auto openingFlux = continuation.openingFlux;
    auto pressureCorrection = continuation.pressureCorrectionPascals;
    result.diagnostics.pressureStep =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            currentPressureOperator, currentBasePressureOperator, grid,
            currentSweep, currentFragments, currentTopology,
            currentVolumeRates, currentOpeningDefinitions, currentOpenings,
            currentResistanceDefinitions, topologyVelocity, openingSamples,
            openingFlux, pressureCorrection, settings,
            limits.pressureStep);
    result.diagnostics.accepted =
        result.diagnostics.pressureStep.accepted;
    result.diagnostics.failureStage = failureStage(
        result.diagnostics.pressureStep);
    if (result.diagnostics.accepted) {
        result.acceptedState =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                currentPressureOperator, currentBasePressureOperator, grid,
                currentSweep, currentFragments, currentTopology,
                currentVolumeRates, currentOpeningDefinitions,
                currentOpenings, currentResistanceDefinitions,
                result.diagnostics.pressureStep, topologyVelocity,
                openingSamples, openingFlux, pressureCorrection, settings,
                limits.acceptedState);
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening pressure-epoch owned-storage limit exceeded");
    }
    validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureEpochResult
acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureEpochLimits& limits) {
    return acceptEpoch(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, currentResistanceDefinitions, settings, limits);
}

void validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureEpochResult& result) {
    const auto expectedFailureStage = failureStage(
        result.diagnostics.pressureStep);
    const bool acceptedPayload = result.acceptedState.fingerprint != 0
        && result.acceptedState.accepted;
    const bool emptyPayload = emptyAcceptedState(result.acceptedState);
    if (result.version
            != planarPressureRegionFragmentOpeningPressureEpochVersion
        || result.sourceAcceptedStateFingerprint == 0
        || result.continuationFingerprint == 0
        || result.continuationOpeningFluxFingerprint == 0
        || result.currentPressureOperatorFingerprint == 0
        || result.currentBasePressureOperatorFingerprint == 0
        || result.currentOpeningFingerprint == 0
        || result.currentFragmentFingerprint == 0
        || result.currentTopologyFingerprint == 0
        || result.currentVolumeRateFingerprint == 0
        || result.currentResistanceDefinitionFingerprint == 0
        || !result.diagnostics.usedConsecutiveContinuation
        || result.diagnostics.accepted
            != result.diagnostics.pressureStep.accepted
        || result.diagnostics.failureStage != expectedFailureStage
        || result.diagnostics.pressureStep.sourceOpeningFluxFingerprint
            != result.continuationOpeningFluxFingerprint
        || result.diagnostics.pressureStep.resistance
                .sourceOpeningFingerprint
            != result.currentOpeningFingerprint
        || result.diagnostics.pressureStep.resistance
                .sourceOpeningFluxFingerprint
            != result.continuationOpeningFluxFingerprint
        || result.diagnostics.pressureStep.resistance
                .resistanceDefinitionFingerprint
            != result.currentResistanceDefinitionFingerprint
        || (result.diagnostics.accepted
                ? (!acceptedPayload
                    || result.diagnostics.failureStage
                        != PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
                            None)
                : (!emptyPayload
                    || result.diagnostics.failureStage
                        == PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
                            None))
        || result.ownedStorageBytes != ownedStorageBytes(result)) {
        throw std::invalid_argument(
            "opening pressure-epoch result integrity is invalid");
    }
    if (result.diagnostics.accepted) {
        validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
            result.acceptedState);
        if (result.acceptedState.sourcePressureOperatorFingerprint
                != result.currentPressureOperatorFingerprint
            || result.acceptedState.sourceBasePressureOperatorFingerprint
                != result.currentBasePressureOperatorFingerprint
            || result.acceptedState.sourceOpeningFingerprint
                != result.currentOpeningFingerprint
            || result.acceptedState.sourceFragmentFingerprint
                != result.currentFragmentFingerprint
            || result.acceptedState.sourceTopologyFingerprint
                != result.currentTopologyFingerprint
            || result.acceptedState.sourceVolumeRateFingerprint
                != result.currentVolumeRateFingerprint
            || result.acceptedState.sourceOpeningFluxFingerprint
                != result.continuationOpeningFluxFingerprint
            || result.acceptedState.resistanceDefinitionFingerprint
                != result.currentResistanceDefinitionFingerprint
            || result.acceptedState.settings != result.settings) {
            throw std::invalid_argument(
                "opening pressure-epoch accepted endpoint is invalid");
        }
    }
}

void validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
    const PlanarPressureRegionFragmentOpeningPressureEpochResult& result,
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureEpochLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
        result);
    const auto expected = acceptEpoch(
        previousState, previousPressureOperator,
        previousBasePressureOperator, grid, previousSweep,
        previousFragments, previousTopology, previousVolumeRates,
        previousOpeningDefinitions, previousOpenings,
        previousResistanceDefinitions, currentPressureOperator,
        currentBasePressureOperator, currentSweep, currentFragments,
        currentTopology, currentVolumeRates, currentOpeningDefinitions,
        currentOpenings, currentResistanceDefinitions, settings, limits);
    if (result != expected) {
        throw std::invalid_argument(
            "opening pressure-epoch result payload is invalid");
    }
}

} // namespace simwing::fsi::fluid
