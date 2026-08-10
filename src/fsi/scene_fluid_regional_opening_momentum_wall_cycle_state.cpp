#include "scene_fluid_regional_opening_momentum_wall_cycle_state.h"

#include <limits>
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
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
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
            "regional opening wall-cycle state storage overflows");
    }
    return first + second;
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    if (limits.maximumWallTractions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "regional opening wall-cycle state limits are invalid");
    }
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state) {
    return checkedAdd(
        state.adjustedMomentum.ownedStorageBytes,
        checkedAdd(
            state.acceptedPressure.ownedStorageBytes,
            state.wallTractions.ownedStorageBytes));
}

std::uint64_t stateFingerprint(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state) {
    Fingerprint fingerprint;
    for (const std::uint64_t value : {
             static_cast<std::uint64_t>(state.version),
             state.sourceWallPressureEpochFingerprint,
             state.sourceWallExchangeFingerprint,
             state.transportMetricFingerprint,
             state.acceptedMetricFingerprint,
             state.predictionFingerprint,
             state.pressureWarmStartFingerprint,
             state.predictedOpeningFluxFingerprint,
             state.adjustedMomentum.fingerprint,
             state.acceptedPressure.fingerprint,
             state.wallTractions.fingerprint,
             static_cast<std::uint64_t>(state.ownedStorageBytes)}) {
        fingerprint.integer(value);
    }
    return fingerprint.value();
}

} // namespace

SceneFluidRegionalOpeningMomentumWallCycleState
captureSceneFluidRegionalOpeningMomentumWallCycleState(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& pressureEpoch,
    const SceneFluidRegionalOpeningMomentumWallExchange& wallExchange,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
        pressureEpoch);
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(
        wallExchange);
    fluid::validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        transportMetric);
    fluid::validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        acceptedMetric);
    validateSceneFluidQuadratureDefinition(quadrature);
    if (!pressureEpoch.pressureEpoch.diagnostics.accepted
        || pressureEpoch.sourceWallExchangeFingerprint
            != wallExchange.fingerprint
        || pressureEpoch.sourceTransportMetricFingerprint
            != transportMetric.fingerprint
        || pressureEpoch.currentMetricFingerprint
            != acceptedMetric.fingerprint) {
        throw std::invalid_argument(
            "regional opening wall-cycle state source is invalid");
    }

    SceneFluidRegionalOpeningMomentumWallCycleState state;
    state.sourceWallPressureEpochFingerprint = pressureEpoch.fingerprint;
    state.sourceWallExchangeFingerprint = wallExchange.fingerprint;
    state.transportMetricFingerprint = transportMetric.fingerprint;
    state.acceptedMetricFingerprint = acceptedMetric.fingerprint;
    state.predictionFingerprint = pressureEpoch.prediction.fingerprint;
    state.pressureWarmStartFingerprint =
        pressureEpoch.sourcePressureWarmStartFingerprint;
    state.predictedOpeningFluxFingerprint =
        pressureEpoch.pressureEpoch
            .sourcePredictedOpeningFluxFingerprint;
    state.adjustedMomentum = pressureEpoch.adjustmentState;
    state.acceptedPressure = pressureEpoch.pressureEpoch.acceptedState;
    state.wallTractions = captureSceneFluidAcceptedWallTractions(
        wallExchange);
    validateSceneFluidAcceptedWallTractions(
        state.wallTractions, quadrature, wallExchange.fingerprint);
    state.ownedStorageBytes = ownedStorageBytes(state);
    if (state.adjustedMomentum.controls.size()
            > limits.adjustmentLimits.maximumFragments
        || state.adjustedMomentum.ownedStorageBytes
            > limits.adjustmentLimits.maximumOwnedBytes
        || state.acceptedPressure.ownedStorageBytes
            > limits.acceptedStateLimits.maximumOwnedBytes
        || state.wallTractions.tractions.size()
            > limits.maximumWallTractions
        || state.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "regional opening wall-cycle state limit exceeded");
    }
    state.fingerprint = stateFingerprint(state);
    validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(state);
    return state;
}

void validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state) {
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
        state.adjustedMomentum);
    fluid::validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
        state.acceptedPressure);
    validateSceneFluidAcceptedWallTractionSetIntegrity(state.wallTractions);
    if (state.version
            != sceneFluidRegionalOpeningMomentumWallCycleStateVersion
        || state.fingerprint == 0
        || state.sourceWallPressureEpochFingerprint == 0
        || state.sourceWallExchangeFingerprint == 0
        || state.transportMetricFingerprint == 0
        || state.acceptedMetricFingerprint == 0
        || state.predictionFingerprint == 0
        || state.predictedOpeningFluxFingerprint == 0
        || state.sourceWallExchangeFingerprint
            != state.adjustedMomentum.sourceAdjustmentFingerprint
        || state.sourceWallExchangeFingerprint
            != state.wallTractions.wallExchangeFingerprint
        || state.transportMetricFingerprint
            != state.adjustedMomentum.sourceMetricFingerprint
        || state.predictedOpeningFluxFingerprint
            != state.acceptedPressure.sourceOpeningFluxFingerprint
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "regional opening wall-cycle state integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallCycleState(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(state);
    fluid::validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        transportMetric);
    fluid::validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        acceptedMetric, grid, acceptedSweep, acceptedFragments,
        acceptedTopology, acceptedBaseMetric, acceptedOpeningDefinitions,
        acceptedOpenings, limits.adjustmentLimits.metricLimits);
    fluid::validatePlanarPressureRegionFragmentOpeningAcceptedState(
        state.acceptedPressure, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, limits.acceptedStateLimits);
    validateSceneFluidAcceptedWallTractions(
        state.wallTractions, quadrature,
        state.sourceWallExchangeFingerprint);
    if (state.transportMetricFingerprint != transportMetric.fingerprint
        || state.acceptedMetricFingerprint != acceptedMetric.fingerprint) {
        throw std::invalid_argument(
            "regional opening wall-cycle state sources are foreign");
    }
    if (state.adjustedMomentum.controls.size()
            > limits.adjustmentLimits.maximumFragments
        || state.adjustedMomentum.ownedStorageBytes
            > limits.adjustmentLimits.maximumOwnedBytes
        || state.wallTractions.tractions.size()
            > limits.maximumWallTractions
        || state.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "regional opening wall-cycle state validation limit exceeded");
    }
}

} // namespace simwing::fsi
