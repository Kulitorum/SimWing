#include "fluid/planar_region_fragment_opening_momentum_cycle.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening momentum-cycle storage overflows");
    }
    return first + second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening momentum-cycle limits are invalid");
    }
}

void addWorkingStorage(
    PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const std::size_t ownedBytes,
    const std::size_t workingBytes,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    result.workingStorageBytes = checkedAdd(
        result.workingStorageBytes,
        checkedAdd(ownedBytes, workingBytes));
    if (result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum-cycle working-storage limit exceeded");
    }
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result) {
    return checkedAdd(
        result.transport.ownedStorageBytes,
        result.acceptedState.ownedStorageBytes);
}

bool emptyTransport(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport) {
    return transport.fingerprint == 0
        && transport.sourceStateFingerprint == 0
        && transport.sourceTransportFingerprint == 0
        && transport.sourceMetricFingerprint == 0
        && transport.targetFlowStateFingerprint == 0
        && transport.targetMetricFingerprint == 0
        && transport.targetVolumeRateFingerprint == 0
        && transport.controls.empty()
        && transport.ownedStorageBytes == 0
        && transport.workingStorageBytes == 0;
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
        && state.pressureCorrectionPascals.empty()
        && state.ownedStorageBytes == 0;
}

PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage pressureFailure(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage
        stage) {
    switch (stage) {
    case PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
        Resistance:
        return PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
            PressureResistance;
    case PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
        PressureProjection:
        return PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
            PressureProjection;
    case PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
        AggregateEnergy:
        return PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
            AggregateEnergy;
    case PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
        None:
        return PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
            None;
    }
    throw std::logic_error(
        "opening momentum-cycle pressure failure stage is invalid");
}

template<typename Source>
PlanarPressureRegionFragmentOpeningMomentumCycleResult buildCycle(
    const Source& sourceArtifact,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        sourceMetric,
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        transportSettings,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    validateLimits(limits);

    PlanarPressureRegionFragmentOpeningMomentumCycleResult result;
    if constexpr (std::is_same_v<
                      Source,
                      PlanarPressureRegionFragmentOpeningVelocityState>) {
        result.sourceStateFingerprint = sourceArtifact.fingerprint;
        result.diagnostics.usedInitialStateTransport = true;
    } else {
        result.sourceTransportFingerprint = sourceArtifact.fingerprint;
        result.diagnostics.usedReentrantTransport = true;
    }
    result.sourceMetricFingerprint = sourceMetric.fingerprint;
    result.currentAcceptedStateFingerprint = currentAcceptedState.fingerprint;
    result.currentPressureOperatorFingerprint =
        currentPressureOperator.fingerprint;
    result.currentBasePressureOperatorFingerprint =
        currentBasePressureOperator.fingerprint;
    result.currentOpeningFingerprint = currentOpenings.fingerprint;
    result.currentFragmentFingerprint = currentFragments.fingerprint;
    result.currentTopologyFingerprint = currentTopology.fingerprint;
    result.currentVolumeRateFingerprint = currentVolumeRates.fingerprint;
    result.currentMetricFingerprint = currentMetric.fingerprint;
    result.nextPressureOperatorFingerprint = nextPressureOperator.fingerprint;
    result.nextBasePressureOperatorFingerprint =
        nextBasePressureOperator.fingerprint;
    result.nextOpeningFingerprint = nextOpenings.fingerprint;
    result.nextFragmentFingerprint = nextFragments.fingerprint;
    result.nextTopologyFingerprint = nextTopology.fingerprint;
    result.nextVolumeRateFingerprint = nextVolumeRates.fingerprint;
    result.nextMetricFingerprint = nextMetric.fingerprint;
    result.transportSettings = transportSettings;
    result.pressureSettings = pressureSettings;

    const auto currentFlow =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            currentAcceptedState, currentPressureOperator,
            currentBasePressureOperator, grid, currentSweep,
            currentFragments, currentTopology, currentVolumeRates,
            currentOpeningDefinitions, currentOpenings,
            currentResistanceDefinitions, currentBaseMetric, currentMetric,
            limits.stateLimits);
    result.currentAcceptedFlowStateFingerprint = currentFlow.fingerprint;
    addWorkingStorage(
        result, currentFlow.ownedStorageBytes,
        currentFlow.workingStorageBytes, limits);

    auto transportAttempt =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            sourceArtifact, sourceMetric, currentFlow,
            currentMetric, grid, currentSweep, currentFragments,
            currentTopology, currentVolumeRates, transportSettings,
            limits.transportLimits);
    result.transportAttemptFingerprint = transportAttempt.fingerprint;
    result.diagnostics.transport = transportAttempt.diagnostics;
    addWorkingStorage(
        result, transportAttempt.ownedStorageBytes,
        transportAttempt.workingStorageBytes, limits);
    if (!transportAttempt.diagnostics.accepted) {
        result.diagnostics.failureStage =
            PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                MomentumTransport;
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
            result);
        return result;
    }

    const auto prediction =
        predictPlanarPressureRegionFragmentOpeningMomentum(
            transportAttempt, currentMetric, grid, nextSweep,
            nextFragments, nextTopology, nextVolumeRates,
            nextOpeningDefinitions, nextOpenings, nextBaseMetric, nextMetric,
            limits.predictionLimits);
    result.predictionFingerprint = prediction.fingerprint;
    result.diagnostics.prediction = prediction.diagnostics;
    addWorkingStorage(
        result, prediction.ownedStorageBytes,
        prediction.workingStorageBytes, limits);

    const auto warmStart =
        buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
            currentAcceptedState, currentPressureOperator,
            currentBasePressureOperator, grid, currentSweep,
            currentFragments, currentTopology, currentVolumeRates,
            currentOpeningDefinitions, currentOpenings,
            currentResistanceDefinitions, nextPressureOperator,
            nextBasePressureOperator, nextSweep, nextFragments,
            nextTopology, nextVolumeRates, nextOpeningDefinitions,
            nextOpenings, limits.warmStartLimits);
    result.pressureWarmStartFingerprint = warmStart.fingerprint;
    result.diagnostics.usedPressureWarmStart = true;
    addWorkingStorage(
        result, warmStart.ownedStorageBytes,
        warmStart.workingStorageBytes, limits);

    auto pressureAttempt =
        acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
            prediction, warmStart, nextPressureOperator,
            nextBasePressureOperator, grid, nextSweep, nextFragments,
            nextTopology, nextVolumeRates, nextOpeningDefinitions,
            nextOpenings, nextResistanceDefinitions, nextBaseMetric,
            nextMetric, pressureSettings, limits.pressureLimits);
    result.predictedOpeningFluxFingerprint =
        pressureAttempt.sourcePredictedOpeningFluxFingerprint;
    result.diagnostics.pressureFailureStage =
        pressureAttempt.diagnostics.failureStage;
    result.diagnostics.pressureIterationCount =
        pressureAttempt.diagnostics.pressureStep.projection.pressureSolve
            .iterationCount;
    result.diagnostics.pressureFinalResidualNormPascalsMeters =
        pressureAttempt.diagnostics.pressureStep.projection.pressureSolve
            .finalResidualL2PascalsMeters;
    addWorkingStorage(
        result, pressureAttempt.ownedStorageBytes,
        pressureAttempt.workingStorageBytes, limits);
    if (!pressureAttempt.diagnostics.accepted) {
        result.diagnostics.failureStage = pressureFailure(
            pressureAttempt.diagnostics.failureStage);
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
            result);
        return result;
    }

    result.diagnostics.accepted = true;
    result.acceptedStateFingerprint =
        pressureAttempt.acceptedState.fingerprint;
    result.transport = std::move(transportAttempt);
    result.acceptedState = std::move(pressureAttempt.acceptedState);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening momentum-cycle owned-storage limit exceeded");
    }
    validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningMomentumCycleResult
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        transportSettings,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    return buildCycle(
        sourceTransport, sourceTransportMetric, currentAcceptedState,
        currentPressureOperator, currentBasePressureOperator, grid,
        currentSweep, currentFragments, currentTopology, currentVolumeRates,
        currentOpeningDefinitions, currentOpenings,
        currentResistanceDefinitions, currentBaseMetric, currentMetric,
        nextPressureOperator, nextBasePressureOperator, nextSweep,
        nextFragments, nextTopology, nextVolumeRates,
        nextOpeningDefinitions, nextOpenings, nextResistanceDefinitions,
        nextBaseMetric, nextMetric, transportSettings, pressureSettings,
        limits);
}

PlanarPressureRegionFragmentOpeningMomentumCycleResult
advancePlanarPressureRegionFragmentOpeningMomentumCycle(
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        transportSettings,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings&
        pressureSettings,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    return buildCycle(
        sourceState, sourceMetric, currentAcceptedState,
        currentPressureOperator, currentBasePressureOperator, grid,
        currentSweep, currentFragments, currentTopology, currentVolumeRates,
        currentOpeningDefinitions, currentOpenings,
        currentResistanceDefinitions, currentBaseMetric, currentMetric,
        nextPressureOperator, nextBasePressureOperator, nextSweep,
        nextFragments, nextTopology, nextVolumeRates,
        nextOpeningDefinitions, nextOpenings, nextResistanceDefinitions,
        nextBaseMetric, nextMetric, transportSettings, pressureSettings,
        limits);
}

void validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result) {
    const bool validSourceLineage =
        (result.sourceStateFingerprint != 0)
        != (result.sourceTransportFingerprint != 0);
    const bool validSourcePolicy =
        result.diagnostics.usedInitialStateTransport
            != result.diagnostics.usedReentrantTransport
        && result.diagnostics.usedInitialStateTransport
            == (result.sourceStateFingerprint != 0);
    const bool endpointsEmpty = emptyTransport(result.transport)
        && emptyAcceptedState(result.acceptedState)
        && result.acceptedStateFingerprint == 0
        && result.ownedStorageBytes == 0;
    const bool endpointsAccepted = result.diagnostics.accepted
        && result.transport.fingerprint == result.transportAttemptFingerprint
        && result.acceptedState.fingerprint == result.acceptedStateFingerprint
        && result.acceptedState.accepted
        && result.ownedStorageBytes == ownedStorageBytes(result);
    const bool stoppedAtTransport =
        result.diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                MomentumTransport
        && !result.diagnostics.transport.accepted
        && result.diagnostics.transport.failureStage
            != PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && result.predictionFingerprint == 0
        && result.pressureWarmStartFingerprint == 0
        && result.predictedOpeningFluxFingerprint == 0
        && result.diagnostics.prediction
            == PlanarPressureRegionFragmentOpeningMomentumPredictionDiagnostics{}
        && !result.diagnostics.usedPressureWarmStart
        && result.diagnostics.pressureFailureStage
            == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                None;
    const bool stoppedAtPressure =
        !result.diagnostics.accepted
        && result.diagnostics.failureStage
            != PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                None
        && result.diagnostics.failureStage
            != PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                MomentumTransport
        && result.diagnostics.transport.accepted
        && result.diagnostics.prediction.finite
        && result.predictionFingerprint != 0
        && result.pressureWarmStartFingerprint != 0
        && result.predictedOpeningFluxFingerprint != 0
        && result.diagnostics.usedPressureWarmStart
        && pressureFailure(result.diagnostics.pressureFailureStage)
            == result.diagnostics.failureStage;
    const bool completed = result.diagnostics.accepted
        && result.diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                None
        && result.diagnostics.transport.accepted
        && result.diagnostics.prediction.finite
        && result.predictionFingerprint != 0
        && result.pressureWarmStartFingerprint != 0
        && result.predictedOpeningFluxFingerprint != 0
        && result.diagnostics.usedPressureWarmStart
        && result.diagnostics.pressureFailureStage
            == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                None;
    if (result.version
            != planarPressureRegionFragmentOpeningMomentumCycleVersion
        || !validSourceLineage
        || !validSourcePolicy
        || result.sourceMetricFingerprint == 0
        || result.currentAcceptedStateFingerprint == 0
        || result.currentPressureOperatorFingerprint == 0
        || result.currentBasePressureOperatorFingerprint == 0
        || result.currentOpeningFingerprint == 0
        || result.currentFragmentFingerprint == 0
        || result.currentTopologyFingerprint == 0
        || result.currentVolumeRateFingerprint == 0
        || result.currentMetricFingerprint == 0
        || result.nextPressureOperatorFingerprint == 0
        || result.nextBasePressureOperatorFingerprint == 0
        || result.nextOpeningFingerprint == 0
        || result.nextFragmentFingerprint == 0
        || result.nextTopologyFingerprint == 0
        || result.nextVolumeRateFingerprint == 0
        || result.nextMetricFingerprint == 0
        || result.currentAcceptedFlowStateFingerprint == 0
        || result.transportAttemptFingerprint == 0
        || !std::isfinite(
            result.diagnostics.pressureFinalResidualNormPascalsMeters)
        || result.workingStorageBytes == 0
        || (result.diagnostics.accepted != completed)
        || (result.diagnostics.accepted
                ? !endpointsAccepted
                : !endpointsEmpty)
        || (!completed && !stoppedAtTransport && !stoppedAtPressure)) {
        throw std::invalid_argument(
            "opening momentum-cycle result integrity is invalid");
    }
    if (result.diagnostics.accepted) {
        validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
            result.transport);
        validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
            result.acceptedState);
        if (result.transport.sourceStateFingerprint
                != result.sourceStateFingerprint
            || result.transport.sourceTransportFingerprint
                != result.sourceTransportFingerprint
            || result.transport.sourceMetricFingerprint
                != result.sourceMetricFingerprint
            || result.transport.targetFlowStateFingerprint
                != result.currentAcceptedFlowStateFingerprint
            || result.transport.targetMetricFingerprint
                != result.currentMetricFingerprint
            || result.transport.targetVolumeRateFingerprint
                != result.currentVolumeRateFingerprint
            || result.transport.settings != result.transportSettings
            || result.acceptedState.sourcePressureOperatorFingerprint
                != result.nextPressureOperatorFingerprint
            || result.acceptedState.sourceBasePressureOperatorFingerprint
                != result.nextBasePressureOperatorFingerprint
            || result.acceptedState.sourceOpeningFingerprint
                != result.nextOpeningFingerprint
            || result.acceptedState.sourceFragmentFingerprint
                != result.nextFragmentFingerprint
            || result.acceptedState.sourceTopologyFingerprint
                != result.nextTopologyFingerprint
            || result.acceptedState.sourceVolumeRateFingerprint
                != result.nextVolumeRateFingerprint
            || result.acceptedState.sourceOpeningFluxFingerprint
                != result.predictedOpeningFluxFingerprint
            || result.acceptedState.settings != result.pressureSettings) {
            throw std::invalid_argument(
                "opening momentum-cycle endpoint lineage is invalid");
        }
    }
}

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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
        result);
    if (result != buildCycle(
            sourceTransport, sourceTransportMetric, currentAcceptedState,
            currentPressureOperator, currentBasePressureOperator, grid,
            currentSweep, currentFragments, currentTopology,
            currentVolumeRates, currentOpeningDefinitions, currentOpenings,
            currentResistanceDefinitions, currentBaseMetric, currentMetric,
            nextPressureOperator, nextBasePressureOperator, nextSweep,
            nextFragments, nextTopology, nextVolumeRates,
            nextOpeningDefinitions, nextOpenings, nextResistanceDefinitions,
            nextBaseMetric, nextMetric, result.transportSettings,
            result.pressureSettings, limits)) {
        throw std::invalid_argument(
            "opening momentum-cycle result is foreign to its source");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        nextOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& nextOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        nextResistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& nextBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& nextMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
        result);
    if (result != buildCycle(
            sourceState, sourceMetric, currentAcceptedState,
            currentPressureOperator, currentBasePressureOperator, grid,
            currentSweep, currentFragments, currentTopology,
            currentVolumeRates, currentOpeningDefinitions, currentOpenings,
            currentResistanceDefinitions, currentBaseMetric, currentMetric,
            nextPressureOperator, nextBasePressureOperator, nextSweep,
            nextFragments, nextTopology, nextVolumeRates,
            nextOpeningDefinitions, nextOpenings, nextResistanceDefinitions,
            nextBaseMetric, nextMetric, result.transportSettings,
            result.pressureSettings, limits)) {
        throw std::invalid_argument(
            "opening momentum-cycle result is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
