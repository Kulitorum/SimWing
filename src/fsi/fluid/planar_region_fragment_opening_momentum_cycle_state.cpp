#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    void integer(std::uint64_t value) {
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    [[nodiscard]] std::uint64_t value() const noexcept { return value_; }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening momentum-cycle state storage overflows");
    }
    return first + second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.transport.maximumFragments == 0
        || limits.transport.maximumDofs == 0
        || limits.transport.maximumOwnedBytes == 0
        || limits.transport.maximumWorkingBytes == 0
        || limits.acceptedState.maximumTopologyLinkVelocities == 0
        || limits.acceptedState.maximumOpeningSamples == 0
        || limits.acceptedState.maximumPressureCorrections == 0
        || limits.acceptedState.maximumOwnedBytes == 0
        || limits.acceptedState.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening momentum-cycle state limits are invalid");
    }
}

void validateStorageLimits(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    if (state.transport.controls.size()
            > limits.transport.maximumFragments
        || state.transport.diagnostics.transportDofCount
            > limits.transport.maximumDofs
        || state.transport.ownedStorageBytes
            > limits.transport.maximumOwnedBytes
        || state.transport.workingStorageBytes
            > limits.transport.maximumWorkingBytes
        || state.acceptedState
                   .orientedTopologyLinkVelocityMetersPerSecond.size()
            > limits.acceptedState.maximumTopologyLinkVelocities
        || state.acceptedState.openingVelocitySamples.size()
            > limits.acceptedState.maximumOpeningSamples
        || state.acceptedState.pressureCorrectionPascals.size()
            > limits.acceptedState.maximumPressureCorrections
        || state.acceptedState.ownedStorageBytes
            > limits.acceptedState.maximumOwnedBytes
        || state.acceptedState.workingStorageBytes
            > limits.acceptedState.maximumWorkingBytes
        || state.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening momentum-cycle state validation limit exceeded");
    }
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state) {
    return checkedAdd(
        state.transport.ownedStorageBytes,
        state.acceptedState.ownedStorageBytes);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.transportMetricFingerprint);
    fingerprint.integer(state.acceptedMetricFingerprint);
    fingerprint.integer(state.currentAcceptedFlowStateFingerprint);
    fingerprint.integer(state.predictionFingerprint);
    fingerprint.integer(state.pressureWarmStartFingerprint);
    fingerprint.integer(state.predictedOpeningFluxFingerprint);
    fingerprint.integer(state.transport.fingerprint);
    fingerprint.integer(state.acceptedState.fingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    return fingerprint.value();
}

bool equalWithinRoundoff(const double first, const double second) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second)
        <= 128.0 * std::numeric_limits<double>::epsilon() * scale;
}

void validateTransportEndpoint(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric) {
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(metric);
    if (state.transportMetricFingerprint != metric.fingerprint
        || state.transport.targetMetricFingerprint != metric.fingerprint
        || state.transport.targetVolumeRateFingerprint
            != volumeRates.fingerprint
        || volumeRates.sourceFragmentFingerprint
            != metric.sourceFragmentFingerprint
        || volumeRates.sourceTopologyFingerprint
            != metric.sourceTopologyFingerprint
        || state.transport.controls.size() != volumeRates.fragments.size()
        || state.transport.controls.size() != metric.fragments.size()) {
        throw std::invalid_argument(
            "opening momentum-cycle transport restart sources are foreign");
    }
    for (std::size_t index = 0;
         index < state.transport.controls.size(); ++index) {
        const auto& control = state.transport.controls[index];
        const auto& rate = volumeRates.fragments[index];
        const auto& fragment = metric.fragments[index];
        if (control.fragmentIndex != index
            || rate.fragmentIndex != index
            || fragment.fragmentIndex != index
            || control.stableId != rate.stableId
            || control.stableId != fragment.stableId
            || control.regionStableId != rate.regionStableId
            || control.regionStableId != fragment.regionStableId
            || control.connectedComponentIndex
                != fragment.connectedComponentIndex
            || !equalWithinRoundoff(
                control.volumeCubicMeters,
                rate.currentVolumeCubicMeters)
            || !equalWithinRoundoff(
                control.volumeCubicMeters,
                fragment.sourceVolumeCubicMeters)) {
            throw std::invalid_argument(
                "opening momentum-cycle transport controls are foreign");
        }
    }
}

} // namespace

PlanarPressureRegionFragmentOpeningMomentumCycleState
capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
        result);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        transportMetric);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        acceptedMetric);
    if (!result.diagnostics.accepted
        || result.currentMetricFingerprint != transportMetric.fingerprint
        || result.nextMetricFingerprint != acceptedMetric.fingerprint) {
        throw std::invalid_argument(
            "opening momentum-cycle state requires one accepted matching cycle");
    }

    PlanarPressureRegionFragmentOpeningMomentumCycleState state;
    state.transportMetricFingerprint = transportMetric.fingerprint;
    state.acceptedMetricFingerprint = acceptedMetric.fingerprint;
    state.currentAcceptedFlowStateFingerprint =
        result.currentAcceptedFlowStateFingerprint;
    state.predictionFingerprint = result.predictionFingerprint;
    state.pressureWarmStartFingerprint =
        result.pressureWarmStartFingerprint;
    state.predictedOpeningFluxFingerprint =
        result.predictedOpeningFluxFingerprint;
    state.transport = result.transport;
    state.acceptedState = result.acceptedState;
    state.ownedStorageBytes = ownedStorageBytes(state);
    state.fingerprint = stateFingerprint(state);
    validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
        state);
    validateStorageLimits(state, limits);
    return state;
}

void validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state) {
    validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        state.transport);
    validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
        state.acceptedState);
    if (state.version
            != planarPressureRegionFragmentOpeningMomentumCycleStateVersion
        || state.fingerprint == 0
        || state.transportMetricFingerprint == 0
        || state.acceptedMetricFingerprint == 0
        || state.currentAcceptedFlowStateFingerprint == 0
        || state.predictionFingerprint == 0
        || state.pressureWarmStartFingerprint == 0
        || state.predictedOpeningFluxFingerprint == 0
        || !state.transport.diagnostics.accepted
        || !state.acceptedState.accepted
        || state.transport.targetFlowStateFingerprint
            != state.currentAcceptedFlowStateFingerprint
        || state.transport.targetMetricFingerprint
            != state.transportMetricFingerprint
        || state.acceptedState.sourceOpeningFluxFingerprint
            != state.predictedOpeningFluxFingerprint
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening momentum-cycle state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    const PlanarPressureRegionFragmentVolumeRateSet& transportVolumeRates,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& acceptedSweep,
    const PlanarPressureRegionFragmentSet& acceptedFragments,
    const PlanarPressureRegionFragmentTopology& acceptedTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
        state);
    validateStorageLimits(state, limits);
    validateTransportEndpoint(
        state, transportVolumeRates, transportMetric);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        acceptedMetric);
    if (state.acceptedMetricFingerprint != acceptedMetric.fingerprint
        || acceptedMetric.sourceFragmentFingerprint
            != state.acceptedState.sourceFragmentFingerprint
        || acceptedMetric.sourceTopologyFingerprint
            != state.acceptedState.sourceTopologyFingerprint
        || acceptedMetric.sourceOpeningFingerprint
            != state.acceptedState.sourceOpeningFingerprint) {
        throw std::invalid_argument(
            "opening momentum-cycle accepted restart metric is foreign");
    }
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        state.acceptedState, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, limits.acceptedState);
}

} // namespace simwing::fsi::fluid
