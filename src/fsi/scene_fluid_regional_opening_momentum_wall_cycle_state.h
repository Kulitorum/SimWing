#pragma once

#include "scene_fluid_regional_opening_momentum_wall_pressure_epoch.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallCycleStateVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallCycleStateLimits {
    fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits
        adjustmentLimits;
    fluid::PlanarPressureRegionFragmentOpeningAcceptedStateLimits
        acceptedStateLimits;
    std::size_t maximumWallTractions = 100'000'000;
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
};

// Compact in-memory accepted endpoint for the opt-in regional wall path.
// Adjusted collocated momentum belongs to transportMetric, accepted pressure
// belongs to acceptedMetric, and conservative wall traction belongs to the
// same source quadrature epoch. Prediction/warm/flux fingerprints retain the
// staggered bridge between those endpoints without retaining transient solver
// diagnostics or the complete wall exchange.
//
// This state is not serialized and is not consumed by a production worker.
struct SceneFluidRegionalOpeningMomentumWallCycleState {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallCycleStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceWallPressureEpochFingerprint = 0;
    std::uint64_t sourceWallExchangeFingerprint = 0;
    std::uint64_t transportMetricFingerprint = 0;
    std::uint64_t acceptedMetricFingerprint = 0;
    std::uint64_t predictionFingerprint = 0;
    std::uint64_t pressureWarmStartFingerprint = 0;
    std::uint64_t predictedOpeningFluxFingerprint = 0;
    fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState
        adjustedMomentum;
    fluid::PlanarPressureRegionFragmentOpeningAcceptedState acceptedPressure;
    SceneFluidAcceptedWallTractionSet wallTractions;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallCycleState&) const =
        default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallCycleState
captureSceneFluidRegionalOpeningMomentumWallCycleState(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& pressureEpoch,
    const SceneFluidRegionalOpeningMomentumWallExchange& wallExchange,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits = {});

void validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state);

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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits = {});

} // namespace simwing::fsi
