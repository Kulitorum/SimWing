#pragma once

#include "fluid/planar_region_fragment_opening_momentum_pressure_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_exchange.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallPressureEpochVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallPressureEpochLimits {
    fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits
        adjustmentLimits;
    fluid::PlanarPressureRegionFragmentOpeningMomentumPredictionLimits
        predictionLimits;
    fluid::PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits
        pressureEpochLimits;
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16384ULL * 1024ULL * 1024ULL;
};

// Complete immutable opt-in continuation from one accepted same-epoch wall
// exchange through identity-preserving collocated adjustment, consecutive
// opening-aware prediction, and private pressure acceptance. Every nested
// artifact is retained so provenance and numerical output can be validated
// independently. The existing momentum cycle and worker do not consume this
// receipt, and Structure traction remains unapplied.
struct SceneFluidRegionalOpeningMomentumWallPressureEpoch {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallPressureEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceWallExchangeFingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceTransportMetricFingerprint = 0;
    std::uint64_t currentMetricFingerprint = 0;
    std::uint64_t sourceAdjustmentStateFingerprint = 0;
    std::uint64_t sourcePredictionFingerprint = 0;
    std::uint64_t sourcePressureWarmStartFingerprint = 0;
    fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState
        adjustmentState;
    fluid::PlanarPressureRegionFragmentOpeningMomentumPrediction prediction;
    fluid::PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult
        pressureEpoch;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallPressureEpoch&) const =
        default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallPressureEpoch
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        adjustmentSettings = {},
    const fluid::PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings = {},
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits =
        {});

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallPressureEpoch
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        adjustmentSettings = {},
    const fluid::PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings = {},
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits =
        {});

void validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& result);

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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits =
        {});

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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidRegionalOpeningMomentumWallPressureEpochLimits& limits =
        {});

} // namespace simwing::fsi
