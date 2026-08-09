#include "scene_fluid_mimetic_pressure_epoch.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace simwing::fsi {
namespace {

void initializeIdentity(
    SceneFluidMimeticPressureEpochResult& result,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources) {
    validateSceneFluidQuadratureDefinition(quadrature);
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidMimeticTraceSystem(fullSystem, controlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedSystem, fullSystem);
    validateSceneFluidMimeticPressureSources(sources, controlCells);
    if (controlCells.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || quadrature.fingerprint == 0
        || quadrature.surfaceDefinitionFingerprint
            != pressureVolumes.surfaceDefinitionFingerprint
        || quadrature.surfaceDefinitionFingerprint
            != controlCells.surfaceDefinitionFingerprint
        || quadrature.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || quadrature.acceptedStepCount != controlCells.acceptedStepCount
        || quadrature.simulationTimeSeconds
            != controlCells.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure epoch identity is invalid");
    }
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.fullTraceSystemFingerprint = fullSystem.fingerprint;
    result.condensedTraceSystemFingerprint = condensedSystem.fingerprint;
    result.pressureSourceFingerprint = sources.fingerprint;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
}

SceneFluidMimeticPressureEpochResult acceptWithWarmStart(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const std::vector<double>& warmStart,
    const std::uint64_t warmStartFingerprint,
    const std::uint64_t topologyTransitionFingerprint,
    const SceneFluidMimeticTraceSolveSettings& solveSettings,
    const SceneFluidMimeticPressureEpochLimits& limits) {
    SceneFluidMimeticPressureEpochResult result;
    initializeIdentity(
        result, quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, sources);
    result.warmStartFingerprint = warmStartFingerprint;
    result.topologyTransitionFingerprint = topologyTransitionFingerprint;
    result.diagnostics.usedConsecutiveWarmStart =
        warmStartFingerprint != 0;
    const auto pressure = solveSceneFluidMimeticPressureSystem(
        condensedSystem, fullSystem, sources, warmStart, solveSettings);
    result.diagnostics.pressureSolve = pressure.diagnostics;
    if (!pressure.diagnostics.accepted) {
        result.diagnostics.failureStage =
            SceneFluidMimeticPressureEpochFailureStage::PressureSolve;
        validateSceneFluidMimeticPressureEpochResultIntegrity(result);
        return result;
    }
    result.acceptedPressureState = captureSceneFluidMimeticPressureState(
        controlCells, fullSystem, condensedSystem, sources, pressure,
        limits.state);
    result.acceptedPressureSamples = sampleSceneFluidMimeticPressure(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, result.acceptedPressureState, limits.sampling);
    result.diagnostics.accepted = true;
    validateSceneFluidMimeticPressureEpochResult(
        result, quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, sources);
    return result;
}

} // namespace

SceneFluidMimeticPressureEpochResult acceptSceneFluidMimeticPressureEpoch(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticTraceSolveSettings& solveSettings,
    const SceneFluidMimeticPressureEpochLimits& limits) {
    if (condensedSystem.traces.size()
        > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        throw std::length_error(
            "scene fluid mimetic pressure epoch warm-start storage overflows");
    }
    const std::size_t warmBytes =
        condensedSystem.traces.size() * sizeof(double);
    if (warmBytes > limits.maximumBootstrapWarmStartBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure epoch warm-start limit exceeded");
    }
    const std::vector<double> zeroWarmStart(
        condensedSystem.traces.size(), 0.0);
    return acceptWithWarmStart(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, sources, zeroWarmStart, 0, 0, solveSettings,
        limits);
}

SceneFluidMimeticPressureEpochResult acceptSceneFluidMimeticPressureEpoch(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticTraceSolveSettings& solveSettings,
    const SceneFluidMimeticPressureEpochLimits& limits) {
    const auto warmStart = buildSceneFluidMimeticPressureWarmStart(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem, controlCells, fullSystem, condensedSystem,
        topologyTransition, limits.warmStart);
    return acceptWithWarmStart(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, sources, warmStart.reducedTracePascals,
        warmStart.fingerprint, topologyTransition.fingerprint,
        solveSettings, limits);
}

void validateSceneFluidMimeticPressureEpochResultIntegrity(
    const SceneFluidMimeticPressureEpochResult& result) {
    const bool acceptedPayload =
        result.acceptedPressureState.fingerprint != 0
        && result.acceptedPressureSamples.fingerprint != 0;
    const bool emptyPayload =
        result.acceptedPressureState.fingerprint == 0
        && result.acceptedPressureState.controls.empty()
        && result.acceptedPressureState.traces.empty()
        && result.acceptedPressureSamples.fingerprint == 0
        && result.acceptedPressureSamples.bindings.empty()
        && result.acceptedPressureSamples.pressures.empty();
    if (result.version != sceneFluidMimeticPressureEpochVersion
        || result.mimeticControlCellFingerprint == 0
        || result.fullTraceSystemFingerprint == 0
        || result.condensedTraceSystemFingerprint == 0
        || result.pressureSourceFingerprint == 0
        || result.quadratureFingerprint == 0
        || result.pressureControlVolumeFingerprint == 0
        || result.structureDefinitionFingerprint == 0
        || !std::isfinite(result.simulationTimeSeconds)
        || result.diagnostics.accepted
            != result.diagnostics.pressureSolve.accepted
        || result.diagnostics.pressureSolve.reducedTraceSolve
                .traceSystemFingerprint
            != result.condensedTraceSystemFingerprint
        || (result.diagnostics.accepted
            && (!result.diagnostics.pressureSolve
                    .reconstructedFullResidualConverged
                || !result.diagnostics.pressureSolve.reducedTraceSolve
                    .compatible
                || !result.diagnostics.pressureSolve.reducedTraceSolve
                    .converged
                || !result.diagnostics.pressureSolve.reducedTraceSolve
                    .finite))
        || result.diagnostics.usedConsecutiveWarmStart
            != (result.warmStartFingerprint != 0)
        || (result.warmStartFingerprint == 0)
            != (result.topologyTransitionFingerprint == 0)
        || (result.diagnostics.accepted
                ? (result.diagnostics.failureStage
                        != SceneFluidMimeticPressureEpochFailureStage::None
                    || !acceptedPayload)
                : (result.diagnostics.failureStage
                        != SceneFluidMimeticPressureEpochFailureStage::
                            PressureSolve
                    || !emptyPayload))) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure epoch result integrity is invalid");
    }
    if (result.diagnostics.accepted) {
        validateSceneFluidMimeticPressureStateIntegrity(
            result.acceptedPressureState);
        validateSceneFluidMimeticPressureSampleIntegrity(
            result.acceptedPressureSamples);
        if (result.acceptedPressureState.fingerprint
                != result.acceptedPressureSamples.pressureStateFingerprint
            || result.acceptedPressureState.pressureSourceFingerprint
                != result.pressureSourceFingerprint) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure epoch accepted payload is invalid");
        }
    }
}

void validateSceneFluidMimeticPressureEpochResult(
    const SceneFluidMimeticPressureEpochResult& result,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources) {
    SceneFluidMimeticPressureEpochResult identity;
    initializeIdentity(
        identity, quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, sources);
    validateSceneFluidMimeticPressureEpochResultIntegrity(result);
    if (result.mimeticControlCellFingerprint
            != identity.mimeticControlCellFingerprint
        || result.fullTraceSystemFingerprint
            != identity.fullTraceSystemFingerprint
        || result.condensedTraceSystemFingerprint
            != identity.condensedTraceSystemFingerprint
        || result.pressureSourceFingerprint
            != identity.pressureSourceFingerprint
        || result.quadratureFingerprint != identity.quadratureFingerprint
        || result.pressureControlVolumeFingerprint
            != identity.pressureControlVolumeFingerprint
        || result.structureDefinitionFingerprint
            != identity.structureDefinitionFingerprint
        || result.acceptedStepCount != identity.acceptedStepCount
        || result.simulationTimeSeconds != identity.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure epoch result identity is invalid");
    }
    if (result.diagnostics.accepted) {
        validateSceneFluidMimeticPressureState(
            result.acceptedPressureState, controlCells, fullSystem,
            condensedSystem);
        validateSceneFluidMimeticPressureSamples(
            result.acceptedPressureSamples, quadrature, pressureVolumes,
            controlCells, fullSystem, condensedSystem,
            result.acceptedPressureState);
    }
}

} // namespace simwing::fsi
