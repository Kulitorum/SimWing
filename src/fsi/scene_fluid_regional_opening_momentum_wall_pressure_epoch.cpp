#include "scene_fluid_regional_opening_momentum_wall_pressure_epoch.h"

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        integer(static_cast<std::make_unsigned_t<Underlying>>(value));
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
            "regional opening wall pressure-epoch storage overflows");
    }
    return first + second;
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening wall pressure-epoch limits are invalid");
    }
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result) {
    return checkedAdd(
        result.adjustmentState.ownedStorageBytes,
        checkedAdd(
            result.prediction.ownedStorageBytes,
            result.pressureEpoch.ownedStorageBytes));
}

std::size_t workingStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result) {
    std::size_t bytes = result.ownedStorageBytes;
    bytes = checkedAdd(bytes, result.adjustmentState.ownedStorageBytes);
    bytes = checkedAdd(bytes, result.prediction.workingStorageBytes);
    bytes = checkedAdd(bytes, result.pressureEpoch.workingStorageBytes);
    return bytes;
}

std::uint64_t resultFingerprint(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result) {
    Fingerprint fingerprint;
    for (const std::uint64_t value : {
             static_cast<std::uint64_t>(result.version),
             result.sourceWallExchangeFingerprint,
             result.sourceTransportFingerprint,
             result.sourceTransportMetricFingerprint,
             result.currentMetricFingerprint,
             result.sourceAdjustmentStateFingerprint,
             result.sourcePredictionFingerprint,
             result.sourcePressureWarmStartFingerprint,
             result.adjustmentState.fingerprint,
             result.prediction.fingerprint,
             result.pressureEpoch.sourcePredictionFingerprint,
             result.pressureEpoch.sourcePredictedVelocityStateFingerprint,
             result.pressureEpoch.sourcePredictedOpeningFluxFingerprint,
             result.pressureEpoch.sourcePressureWarmStartFingerprint,
             result.pressureEpoch.currentPressureOperatorFingerprint,
             result.pressureEpoch.currentBasePressureOperatorFingerprint,
             result.pressureEpoch.currentOpeningFingerprint,
             result.pressureEpoch.currentFragmentFingerprint,
             result.pressureEpoch.currentTopologyFingerprint,
             result.pressureEpoch.currentVolumeRateFingerprint,
             result.pressureEpoch.currentResistanceDefinitionFingerprint,
             result.pressureEpoch.acceptedState.fingerprint,
             static_cast<std::uint64_t>(result.ownedStorageBytes),
             static_cast<std::uint64_t>(result.workingStorageBytes)}) {
        fingerprint.integer(value);
    }
    fingerprint.enumeration(result.pressureEpoch.diagnostics.failureStage);
    fingerprint.integer(static_cast<std::uint8_t>(
        result.pressureEpoch.diagnostics.accepted ? 1U : 0U));
    return fingerprint.value();
}

SceneFluidRegionalOpeningMomentumWallPressureEpoch buildResult(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureWarmStart*
        warmStart,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        adjustmentSettings,
    const fluid::PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    validateLimits(limits);
    SceneFluidRegionalOpeningMomentumWallPressureEpoch result;
    result.adjustmentState =
        captureSceneFluidRegionalOpeningMomentumAdjustmentState(
            exchange, transport, transportTargetMetric, adjustmentSettings,
            limits.adjustmentLimits);
    result.prediction =
        fluid::predictPlanarPressureRegionFragmentOpeningMomentum(
            result.adjustmentState, transportTargetMetric, grid, sweep,
            fragments, topology, volumeRates, openingDefinitions, openings,
            baseMetric, metric, limits.predictionLimits);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
        result.prediction, result.adjustmentState, transportTargetMetric,
        grid, sweep, fragments, topology, volumeRates, openingDefinitions,
        openings, baseMetric, metric, limits.predictionLimits);
    if (warmStart == nullptr) {
        result.pressureEpoch =
            fluid::acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                result.prediction, pressureOperator, basePressureOperator,
                grid, sweep, fragments, topology, volumeRates,
                openingDefinitions, openings, resistanceDefinitions,
                baseMetric, metric, pressureSettings,
                limits.pressureEpochLimits);
        fluid::validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
            result.pressureEpoch, result.prediction, pressureOperator,
            basePressureOperator, grid, sweep, fragments, topology,
            volumeRates, openingDefinitions, openings,
            resistanceDefinitions, baseMetric, metric, pressureSettings,
            limits.pressureEpochLimits);
    } else {
        result.pressureEpoch =
            fluid::acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                result.prediction, *warmStart, pressureOperator,
                basePressureOperator, grid, sweep, fragments, topology,
                volumeRates, openingDefinitions, openings,
                resistanceDefinitions, baseMetric, metric, pressureSettings,
                limits.pressureEpochLimits);
        fluid::validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
            result.pressureEpoch, result.prediction, *warmStart,
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, baseMetric, metric, pressureSettings,
            limits.pressureEpochLimits);
    }

    result.sourceWallExchangeFingerprint = exchange.fingerprint;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.sourceTransportMetricFingerprint =
        transportTargetMetric.fingerprint;
    result.currentMetricFingerprint = metric.fingerprint;
    result.sourceAdjustmentStateFingerprint =
        result.adjustmentState.fingerprint;
    result.sourcePredictionFingerprint = result.prediction.fingerprint;
    result.sourcePressureWarmStartFingerprint = warmStart == nullptr
        ? 0 : warmStart->fingerprint;
    result.ownedStorageBytes = ownedStorageBytes(result);
    result.workingStorageBytes = workingStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening wall pressure-epoch aggregate limit exceeded");
    }
    result.fingerprint = resultFingerprint(result);
    validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
        result);
    return result;
}

} // namespace

SceneFluidRegionalOpeningMomentumWallPressureEpoch
acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        adjustmentSettings,
    const fluid::PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    return buildResult(
        exchange, transport, transportTargetMetric, nullptr,
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, baseMetric, metric, adjustmentSettings,
        pressureSettings, limits);
}

SceneFluidRegionalOpeningMomentumWallPressureEpoch
acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureWarmStart&
        warmStart,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        adjustmentSettings,
    const fluid::PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    return buildResult(
        exchange, transport, transportTargetMetric, &warmStart,
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, baseMetric, metric, adjustmentSettings,
        pressureSettings, limits);
}

void validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result) {
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
        result.adjustmentState);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
        result.prediction);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
        result.pressureEpoch);
    if (result.version
            != sceneFluidRegionalOpeningMomentumWallPressureEpochVersion
        || result.fingerprint == 0
        || result.sourceWallExchangeFingerprint == 0
        || result.sourceTransportFingerprint == 0
        || result.sourceTransportMetricFingerprint == 0
        || result.currentMetricFingerprint == 0
        || result.sourceAdjustmentStateFingerprint == 0
        || result.sourcePredictionFingerprint == 0
        || result.sourceWallExchangeFingerprint
            != result.adjustmentState.sourceAdjustmentFingerprint
        || result.sourceTransportFingerprint
            != result.adjustmentState.sourceTransportFingerprint
        || result.sourceTransportFingerprint
            != result.prediction.sourceTransportFingerprint
        || result.sourceTransportMetricFingerprint
            != result.adjustmentState.sourceMetricFingerprint
        || result.sourceTransportMetricFingerprint
            != result.prediction.sourceTransportMetricFingerprint
        || result.currentMetricFingerprint
            != result.prediction.currentMetricFingerprint
        || result.sourceAdjustmentStateFingerprint
            != result.adjustmentState.fingerprint
        || result.sourceAdjustmentStateFingerprint
            != result.prediction.sourceAdjustmentStateFingerprint
        || result.sourcePredictionFingerprint != result.prediction.fingerprint
        || result.sourcePredictionFingerprint
            != result.pressureEpoch.sourcePredictionFingerprint
        || result.sourcePressureWarmStartFingerprint
            != result.pressureEpoch.sourcePressureWarmStartFingerprint
        || result.prediction.predictedVelocityState.fingerprint
            != result.pressureEpoch
                   .sourcePredictedVelocityStateFingerprint
        || result.ownedStorageBytes != ownedStorageBytes(result)
        || result.workingStorageBytes != workingStorageBytes(result)
        || result.fingerprint != resultFingerprint(result)) {
        throw std::invalid_argument(
            "regional opening wall pressure-epoch integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallPressureEpoch(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result,
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(result);
    if (result != buildResult(
            exchange, transport, transportTargetMetric, nullptr,
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, baseMetric, metric,
            result.adjustmentState.settings, result.pressureEpoch.settings,
            limits)) {
        throw std::invalid_argument(
            "regional opening wall pressure-epoch sources are foreign");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallPressureEpoch(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result,
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureWarmStart&
        warmStart,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(result);
    if (result != buildResult(
            exchange, transport, transportTargetMetric, &warmStart,
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, baseMetric, metric,
            result.adjustmentState.settings, result.pressureEpoch.settings,
            limits)) {
        throw std::invalid_argument(
            "regional opening wall pressure-epoch sources are foreign");
    }
}

} // namespace simwing::fsi
