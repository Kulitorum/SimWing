#include "fluid/planar_region_fragment_opening_momentum_pressure_epoch.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace simwing::fsi::fluid {
namespace {

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening momentum pressure-epoch storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening momentum pressure-epoch storage overflows");
    }
    return first * second;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult&
        result) {
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

PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage
failureStage(
    const PlanarPressureRegionFragmentOpeningPressureStepDiagnostics& step) {
    if (step.accepted) {
        return PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
            None;
    }
    if (!step.resistance.accepted) {
        return PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
            Resistance;
    }
    if (!step.projection.accepted) {
        return PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
            PressureProjection;
    }
    return PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
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

PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult acceptEpoch(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningPressureWarmStart* warmStart,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening momentum pressure-epoch limits are invalid");
    }
    validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
        prediction);
    validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        metric, grid, sweep, fragments, topology, baseMetric,
        openingDefinitions, openings, limits.stateLimits.metricLimits);
    validatePlanarPressureRegionFragmentOpeningVelocityState(
        prediction.predictedVelocityState, metric, limits.stateLimits);
    if (prediction.currentMetricFingerprint != metric.fingerprint
        || prediction.currentVolumeRateFingerprint != volumeRates.fingerprint
        || prediction.densityKgPerCubicMeter
            != settings.projection.densityKgPerCubicMeter
        || prediction.timeStepSeconds
            != settings.projection.timeStepSeconds
        || prediction.predictedVelocityState.samples.size()
            != metric.dofs.size()) {
        throw std::invalid_argument(
            "opening momentum pressure-epoch prediction is foreign");
    }
    if (warmStart != nullptr) {
        validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
            *warmStart);
        const auto& warmLimits = limits.warmStartLimits;
        if (warmLimits.maximumPressureCorrections == 0
            || warmLimits.maximumComponents == 0
            || warmLimits.maximumOwnedBytes == 0
            || warmLimits.maximumWorkingBytes == 0) {
            throw std::invalid_argument(
                "opening momentum pressure-epoch warm-start limits are invalid");
        }
        if (warmStart->pressureCorrectionCount
                > warmLimits.maximumPressureCorrections
            || warmStart->connectedComponentCount
                > warmLimits.maximumComponents
            || warmStart->ownedStorageBytes
                > warmLimits.maximumOwnedBytes
            || warmStart->workingStorageBytes
                > warmLimits.maximumWorkingBytes) {
            throw std::length_error(
                "opening momentum pressure-epoch warm start exceeds limits");
        }
        if (warmStart->currentPressureOperatorFingerprint
                != pressureOperator.fingerprint
            || warmStart->currentBasePressureOperatorFingerprint
                != basePressureOperator.fingerprint
            || warmStart->currentOpeningFingerprint != openings.fingerprint
            || warmStart->currentFragmentFingerprint != fragments.fingerprint
            || warmStart->currentTopologyFingerprint != topology.fingerprint
            || warmStart->currentVolumeRateFingerprint
                != volumeRates.fingerprint
            || warmStart->pressureCorrectionPascals.size()
                != fragments.fragments.size()
            || warmStart->connectedComponentCount
                != pressureOperator.components.size()) {
            throw std::invalid_argument(
                "opening momentum pressure-epoch warm start is foreign");
        }
    }

    const std::size_t baseWorkingBytes = checkedAdd(
        checkedMultiply(topology.links.size(), sizeof(double)),
        checkedAdd(
            checkedMultiply(
                openings.patches.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocitySample)),
            checkedAdd(
                checkedMultiply(fragments.fragments.size(), sizeof(double)),
                checkedAdd(
                    checkedMultiply(topology.links.size(),
                                    sizeof(std::uint8_t)),
                    checkedMultiply(openings.patches.size(),
                                    sizeof(std::uint8_t))))));
    if (baseWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum pressure-epoch working-storage limit exceeded");
    }

    std::vector<double> topologyVelocity(topology.links.size(), 0.0);
    std::vector<std::uint8_t> assignedLinks(topology.links.size(), 0);
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        openingSamples(openings.patches.size());
    std::vector<std::uint8_t> assignedOpenings(openings.patches.size(), 0);
    for (const auto& sample : prediction.predictedVelocityState.samples) {
        const auto& dof = metric.dofs.at(sample.dofIndex);
        if (sample.stableId != dof.stableId
            || sample.kind != dof.kind
            || sample.sourceFaceLinkStableId
                != dof.sourceFaceLinkStableId) {
            throw std::logic_error(
                "opening momentum pressure-epoch sample binding is invalid");
        }
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            if (dof.sourceFaceLinkIndex >= topologyVelocity.size()
                || assignedLinks[dof.sourceFaceLinkIndex] != 0) {
                throw std::logic_error(
                    "opening momentum pressure-epoch grid coverage is invalid");
            }
            topologyVelocity[dof.sourceFaceLinkIndex] =
                sample.relativeNormalVelocityMetersPerSecond;
            assignedLinks[dof.sourceFaceLinkIndex] = 1;
        } else if (dof.kind
                   == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                       OpeningPatch) {
            if (dof.sourceOpeningPatchIndex >= openingSamples.size()
                || assignedOpenings[dof.sourceOpeningPatchIndex] != 0) {
                throw std::logic_error(
                    "opening momentum pressure-epoch aperture coverage is invalid");
            }
            openingSamples[dof.sourceOpeningPatchIndex] = {
                dof.sourceOpeningPatchStableId,
                sample.relativeNormalVelocityMetersPerSecond,
            };
            assignedOpenings[dof.sourceOpeningPatchIndex] = 1;
        }
    }
    for (const auto& link : topology.links) {
        const bool assigned = assignedLinks[link.linkIndex] != 0;
        if (assigned
                != (link.kind
                    == PlanarPressureRegionFragmentFaceKind::SameRegionGrid)
            || (link.kind
                    == PlanarPressureRegionFragmentFaceKind::PressureLayerWall
                && topologyVelocity[link.linkIndex] != 0.0)) {
            throw std::logic_error(
                "opening momentum pressure-epoch topology coverage is invalid");
        }
    }
    if (!std::ranges::all_of(
            assignedOpenings,
            [](const std::uint8_t assigned) { return assigned != 0; })) {
        throw std::logic_error(
            "opening momentum pressure-epoch aperture coverage is incomplete");
    }

    auto openingFlux =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            grid, sweep, fragments, topology, openingDefinitions,
            openings, openingSamples, limits.openingFluxLimits);
    const std::size_t workingBytes = checkedAdd(
        baseWorkingBytes, openingFlux.ownedStorageBytes);
    if (workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum pressure-epoch aggregate working storage exceeded");
    }
    std::vector<double> pressureCorrection = warmStart == nullptr
        ? std::vector<double>(fragments.fragments.size(), 0.0)
        : warmStart->pressureCorrectionPascals;

    PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult result;
    result.sourcePredictionFingerprint = prediction.fingerprint;
    result.sourcePredictedVelocityStateFingerprint =
        prediction.predictedVelocityState.fingerprint;
    result.sourcePredictedOpeningFluxFingerprint = openingFlux.fingerprint;
    result.sourcePressureWarmStartFingerprint = warmStart == nullptr
        ? 0 : warmStart->fingerprint;
    result.currentPressureOperatorFingerprint = pressureOperator.fingerprint;
    result.currentBasePressureOperatorFingerprint =
        basePressureOperator.fingerprint;
    result.currentOpeningFingerprint = openings.fingerprint;
    result.currentFragmentFingerprint = fragments.fingerprint;
    result.currentTopologyFingerprint = topology.fingerprint;
    result.currentVolumeRateFingerprint = volumeRates.fingerprint;
    result.settings = settings;
    result.workingStorageBytes = workingBytes;
    result.diagnostics.usedTransportedPrediction = true;
    result.diagnostics.usedColdPressureStart = warmStart == nullptr;
    result.diagnostics.usedWarmPressureStart = warmStart != nullptr;
    result.diagnostics.pressureStep =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, topologyVelocity, openingSamples,
            openingFlux, pressureCorrection, settings,
            limits.pressureStepLimits);
    result.currentResistanceDefinitionFingerprint =
        result.diagnostics.pressureStep.resistance
            .resistanceDefinitionFingerprint;
    result.diagnostics.accepted = result.diagnostics.pressureStep.accepted;
    result.diagnostics.failureStage = failureStage(
        result.diagnostics.pressureStep);
    if (result.diagnostics.accepted) {
        result.acceptedState =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                pressureOperator, basePressureOperator, grid, sweep,
                fragments, topology, volumeRates, openingDefinitions,
                openings, resistanceDefinitions,
                result.diagnostics.pressureStep, topologyVelocity,
                openingSamples, openingFlux, pressureCorrection, settings,
                limits.acceptedStateLimits);
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening momentum pressure-epoch owned-storage limit exceeded");
    }
    validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
        result);
    return result;
}

} // namespace

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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits) {
    return acceptEpoch(
        prediction, nullptr, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, baseMetric, metric, settings, limits);
}

PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult
acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits) {
    return acceptEpoch(
        prediction, &warmStart, pressureOperator, basePressureOperator, grid,
        sweep, fragments, topology, volumeRates, openingDefinitions,
        openings, resistanceDefinitions, baseMetric, metric, settings,
        limits);
}

void validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult&
        result) {
    const auto expectedFailureStage = failureStage(
        result.diagnostics.pressureStep);
    const bool acceptedPayload = result.acceptedState.fingerprint != 0
        && result.acceptedState.accepted;
    const bool emptyPayload = emptyAcceptedState(result.acceptedState);
    if (result.version
            != planarPressureRegionFragmentOpeningMomentumPressureEpochVersion
        || result.sourcePredictionFingerprint == 0
        || result.sourcePredictedVelocityStateFingerprint == 0
        || result.sourcePredictedOpeningFluxFingerprint == 0
        || result.currentPressureOperatorFingerprint == 0
        || result.currentBasePressureOperatorFingerprint == 0
        || result.currentOpeningFingerprint == 0
        || result.currentFragmentFingerprint == 0
        || result.currentTopologyFingerprint == 0
        || result.currentVolumeRateFingerprint == 0
        || result.currentResistanceDefinitionFingerprint == 0
        || !result.diagnostics.usedTransportedPrediction
        || result.diagnostics.usedColdPressureStart
            == result.diagnostics.usedWarmPressureStart
        || result.diagnostics.usedColdPressureStart
            != (result.sourcePressureWarmStartFingerprint == 0)
        || result.diagnostics.accepted
            != result.diagnostics.pressureStep.accepted
        || result.diagnostics.failureStage != expectedFailureStage
        || result.diagnostics.pressureStep.sourceOpeningFluxFingerprint
            != result.sourcePredictedOpeningFluxFingerprint
        || result.diagnostics.pressureStep.resistance
                .sourceOpeningFingerprint
            != result.currentOpeningFingerprint
        || result.diagnostics.pressureStep.resistance
                .sourceOpeningFluxFingerprint
            != result.sourcePredictedOpeningFluxFingerprint
        || result.diagnostics.pressureStep.resistance
                .resistanceDefinitionFingerprint
            != result.currentResistanceDefinitionFingerprint
        || (result.diagnostics.accepted
                ? (!acceptedPayload
                    || result.diagnostics.failureStage
                        != PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                            None)
                : (!emptyPayload
                    || result.diagnostics.failureStage
                        == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                            None))
        || result.ownedStorageBytes != ownedStorageBytes(result)
        || result.workingStorageBytes == 0) {
        throw std::invalid_argument(
            "opening momentum pressure-epoch result integrity is invalid");
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
                != result.sourcePredictedOpeningFluxFingerprint
            || result.acceptedState.resistanceDefinitionFingerprint
                != result.currentResistanceDefinitionFingerprint
            || result.acceptedState.settings != result.settings) {
            throw std::invalid_argument(
                "opening momentum pressure-epoch endpoint is invalid");
        }
    }
}

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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
        result);
    if (result != acceptEpoch(
            prediction, nullptr, pressureOperator, basePressureOperator, grid,
            sweep,
            fragments, topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, baseMetric, metric, settings, limits)) {
        throw std::invalid_argument(
            "opening momentum pressure-epoch result is foreign");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochResult& result,
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits&
        limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
        result);
    if (result != acceptEpoch(
            prediction, &warmStart, pressureOperator, basePressureOperator,
            grid, sweep, fragments, topology, volumeRates,
            openingDefinitions, openings, resistanceDefinitions, baseMetric,
            metric, settings, limits)) {
        throw std::invalid_argument(
            "opening momentum pressure-epoch result is foreign");
    }
}

} // namespace simwing::fsi::fluid
