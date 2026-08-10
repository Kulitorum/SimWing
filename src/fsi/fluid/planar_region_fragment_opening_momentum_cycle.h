#pragma once

#include "fluid/planar_region_fragment_opening_momentum_pressure_epoch.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumCycleVersion = 1;

enum class PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage
    : std::uint8_t {
    None = 0,
    MomentumTransport = 1,
    PressureResistance = 2,
    PressureProjection = 3,
    AggregateEnergy = 4,
};

struct PlanarPressureRegionFragmentOpeningMomentumCycleLimits {
    PlanarPressureRegionFragmentOpeningVelocityStateLimits stateLimits;
    PlanarPressureRegionFragmentOpeningMomentumTransportLimits transportLimits;
    PlanarPressureRegionFragmentOpeningMomentumPredictionLimits predictionLimits;
    PlanarPressureRegionFragmentOpeningPressureWarmStartLimits warmStartLimits;
    PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits pressureLimits;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 32ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningMomentumCycleDiagnostics {
    bool accepted = false;
    bool usedReentrantTransport = false;
    bool usedPressureWarmStart = false;
    PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage failureStage =
        PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::None;
    PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics transport;
    PlanarPressureRegionFragmentOpeningMomentumPredictionDiagnostics prediction;
    PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage
        pressureFailureStage =
            PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                None;
    std::size_t pressureIterationCount = 0;
    double pressureFinalResidualNormPascalsMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumCycleDiagnostics&)
        const = default;
};

// One atomic topology-stable transported opening cycle. The source transport
// owns collocated momentum at the previous geometry. The current accepted
// pressure endpoint supplies corrected target flow for the re-entrant ALE
// advance. That new momentum is predicted onto the next geometry; only the
// prior correction pressure is gauge-mapped there before pressure acceptance.
//
// The public endpoint is the staggered pair required by the next cycle:
// transported momentum at the current geometry and accepted pressure at the
// next geometry. A transport or pressure rejection retains scalar diagnostics
// and lineage but publishes neither endpoint.
struct PlanarPressureRegionFragmentOpeningMomentumCycleResult {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumCycleVersion;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceTransportMetricFingerprint = 0;
    std::uint64_t currentAcceptedStateFingerprint = 0;
    std::uint64_t currentPressureOperatorFingerprint = 0;
    std::uint64_t currentBasePressureOperatorFingerprint = 0;
    std::uint64_t currentOpeningFingerprint = 0;
    std::uint64_t currentFragmentFingerprint = 0;
    std::uint64_t currentTopologyFingerprint = 0;
    std::uint64_t currentVolumeRateFingerprint = 0;
    std::uint64_t currentMetricFingerprint = 0;
    std::uint64_t nextPressureOperatorFingerprint = 0;
    std::uint64_t nextBasePressureOperatorFingerprint = 0;
    std::uint64_t nextOpeningFingerprint = 0;
    std::uint64_t nextFragmentFingerprint = 0;
    std::uint64_t nextTopologyFingerprint = 0;
    std::uint64_t nextVolumeRateFingerprint = 0;
    std::uint64_t nextMetricFingerprint = 0;
    std::uint64_t currentAcceptedFlowStateFingerprint = 0;
    std::uint64_t transportAttemptFingerprint = 0;
    std::uint64_t predictionFingerprint = 0;
    std::uint64_t pressureWarmStartFingerprint = 0;
    std::uint64_t predictedOpeningFluxFingerprint = 0;
    std::uint64_t acceptedStateFingerprint = 0;
    PlanarPressureRegionFragmentOpeningMomentumTransportSettings
        transportSettings;
    PlanarPressureRegionFragmentOpeningPressureStepSettings pressureSettings;
    PlanarPressureRegionFragmentOpeningMomentumCycleDiagnostics diagnostics;
    PlanarPressureRegionFragmentOpeningMomentumTransport transport;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedState;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumCycleResult&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningMomentumCycleResult
advancePlanarPressureRegionFragmentOpeningMomentumCycle(
    const PlanarPressureRegionFragmentOpeningMomentumTransport&
        sourceTransport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        sourceTransportMetric,
    const PlanarPressureRegionFragmentOpeningAcceptedState&
        currentAcceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        nextPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        nextBasePressureOperator,
    const PlanarPressureRegionSweepLedger& nextSweep,
    const PlanarPressureRegionFragmentSet& nextFragments,
    const PlanarPressureRegionFragmentTopology& nextTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& nextVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        transportSettings = {},
    const PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings = {},
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result);

void validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const PlanarPressureRegionFragmentOpeningMomentumTransport&
        sourceTransport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        sourceTransportMetric,
    const PlanarPressureRegionFragmentOpeningAcceptedState&
        currentAcceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        nextPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        nextBasePressureOperator,
    const PlanarPressureRegionSweepLedger& nextSweep,
    const PlanarPressureRegionFragmentSet& nextFragments,
    const PlanarPressureRegionFragmentTopology& nextTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& nextVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits = {});

} // namespace simwing::fsi::fluid
