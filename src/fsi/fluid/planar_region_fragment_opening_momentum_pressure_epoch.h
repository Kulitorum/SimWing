#pragma once

#include "fluid/planar_region_fragment_opening_momentum_prediction.h"
#include "fluid/planar_region_fragment_opening_pressure_step.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumPressureEpochVersion = 1;

enum class
    PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage
    : std::uint8_t {
    None = 0,
    Resistance = 1,
    PressureProjection = 2,
    AggregateEnergy = 3,
};

struct PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits {
    PlanarPressureRegionFragmentOpeningVelocityStateLimits stateLimits;
    PlanarPressureRegionFragmentOpeningFluxLimits openingFluxLimits;
    PlanarPressureRegionFragmentOpeningPressureStepLimits pressureStepLimits;
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits acceptedStateLimits;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 8192ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningMomentumPressureEpochDiagnostics {
    bool accepted = false;
    bool usedTransportedPrediction = false;
    bool usedColdPressureStart = false;
    PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage
        failureStage =
            PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                None;
    PlanarPressureRegionFragmentOpeningPressureStepDiagnostics pressureStep;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumPressureEpochDiagnostics&)
        const = default;
};

// Atomic pressure transaction seeded from one transported consecutive
// predictor. Same-region relative velocities and aperture-relative samples
// are extracted from the predictor's immutable opening-aware state; retained
// pressure-layer links remain zero. The exact opening flux is rebuilt, a zero
// correction-pressure gauge is used as an explicit cold start, and resistance
// plus augmented projection run privately. Only a completely accepted step
// publishes an immutable accepted endpoint.
//
// This opt-in composition does not yet carry a pressure warm start, apply
// traction, perform viscosity/wall shear, rebase topology, or select the
// production worker.
struct PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumPressureEpochVersion;
    std::uint64_t sourcePredictionFingerprint = 0;
    std::uint64_t sourcePredictedVelocityStateFingerprint = 0;
    std::uint64_t sourcePredictedOpeningFluxFingerprint = 0;
    std::uint64_t currentPressureOperatorFingerprint = 0;
    std::uint64_t currentBasePressureOperatorFingerprint = 0;
    std::uint64_t currentOpeningFingerprint = 0;
    std::uint64_t currentFragmentFingerprint = 0;
    std::uint64_t currentTopologyFingerprint = 0;
    std::uint64_t currentVolumeRateFingerprint = 0;
    std::uint64_t currentResistanceDefinitionFingerprint = 0;
    PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
    PlanarPressureRegionFragmentOpeningMomentumPressureEpochDiagnostics
        diagnostics;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedState;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult&)
        const = default;
};

[[nodiscard]]
PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult
acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits = {});

void validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult&
        result);

void validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult& result,
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits = {});

} // namespace simwing::fsi::fluid
