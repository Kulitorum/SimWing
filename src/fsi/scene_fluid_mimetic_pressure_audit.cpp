#include "scene_fluid_mimetic_pressure_audit.h"

#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

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

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    template<typename Enum>
    void enumeration(const Enum value) {
        integer(static_cast<std::underlying_type_t<Enum>>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedStorageSum(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint) {
    const auto& solve = endpoint.pressureEpoch.diagnostics.pressureSolve
        .reducedTraceSolve;
    if (solve.components.size()
        > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticTraceSolveComponentDiagnostics)) {
        throw std::length_error(
            "scene fluid mimetic pressure-audit storage overflows");
    }
    const std::size_t solveDiagnosticBytes = solve.components.size()
        * sizeof(SceneFluidMimeticTraceSolveComponentDiagnostics);
    const std::size_t values[]{
        endpoint.controlCells.ownedStorageBytes,
        endpoint.fullTraceSystem.ownedStorageBytes,
        endpoint.condensedTraceSystem.ownedStorageBytes,
        endpoint.predictedTraceFlows.ownedStorageBytes,
        endpoint.pressureSources.ownedStorageBytes,
        endpoint.pressureEpoch.acceptedPressureState.ownedStorageBytes,
        endpoint.pressureEpoch.acceptedPressureSamples.ownedStorageBytes,
        solveDiagnosticBytes,
    };
    std::size_t result = 0;
    for (const std::size_t value : values) {
        if (value > std::numeric_limits<std::size_t>::max() - result) {
            throw std::length_error(
                "scene fluid mimetic pressure-audit storage overflows");
        }
        result += value;
    }
    return result;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint) {
    Fingerprint fingerprint;
    fingerprint.integer(endpoint.version);
    fingerprint.integer(endpoint.scenePressureEpochFingerprint);
    fingerprint.integer(
        endpoint.pressureTopologyTransitionFingerprint);
    fingerprint.integer(endpoint.structureDefinitionFingerprint);
    fingerprint.integer(endpoint.acceptedStepCount);
    fingerprint.real(endpoint.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint8_t>(
        endpoint.usesRegionWallPrediction));
    fingerprint.integer(static_cast<std::uint8_t>(
        endpoint.usesConsecutiveWarmStart));
    fingerprint.integer(static_cast<std::uint64_t>(
        endpoint.ownedStorageBytes));
    fingerprint.integer(endpoint.controlCells.fingerprint);
    fingerprint.integer(endpoint.fullTraceSystem.fingerprint);
    fingerprint.integer(endpoint.condensedTraceSystem.fingerprint);
    fingerprint.integer(endpoint.predictedTraceFlows.fingerprint);
    fingerprint.integer(endpoint.pressureSources.fingerprint);
    fingerprint.integer(endpoint.pressureEpoch.version);
    fingerprint.integer(
        endpoint.pressureEpoch.mimeticControlCellFingerprint);
    fingerprint.integer(endpoint.pressureEpoch.fullTraceSystemFingerprint);
    fingerprint.integer(
        endpoint.pressureEpoch.condensedTraceSystemFingerprint);
    fingerprint.integer(endpoint.pressureEpoch.pressureSourceFingerprint);
    fingerprint.integer(endpoint.pressureEpoch.warmStartFingerprint);
    fingerprint.integer(
        endpoint.pressureEpoch.topologyTransitionFingerprint);
    fingerprint.integer(endpoint.pressureEpoch.quadratureFingerprint);
    fingerprint.integer(
        endpoint.pressureEpoch.pressureControlVolumeFingerprint);
    fingerprint.integer(
        endpoint.pressureEpoch.structureDefinitionFingerprint);
    fingerprint.integer(endpoint.pressureEpoch.acceptedStepCount);
    fingerprint.real(endpoint.pressureEpoch.simulationTimeSeconds);
    const auto& epochDiagnostics = endpoint.pressureEpoch.diagnostics;
    fingerprint.integer(static_cast<std::uint8_t>(
        epochDiagnostics.accepted));
    fingerprint.integer(static_cast<std::uint8_t>(
        epochDiagnostics.usedConsecutiveWarmStart));
    fingerprint.enumeration(epochDiagnostics.failureStage);
    const auto& pressureSolve = epochDiagnostics.pressureSolve;
    fingerprint.integer(static_cast<std::uint8_t>(
        pressureSolve.accepted));
    fingerprint.integer(static_cast<std::uint8_t>(
        pressureSolve.reconstructedFullResidualConverged));
    fingerprint.real(
        pressureSolve.reconstructedFullResidualTolerancePascalsMeters);
    fingerprint.real(
        pressureSolve.reconstructedFullResidualL2PascalsMeters);
    fingerprint.real(
        pressureSolve.reconstructedFullResidualMaximumPascalsMeters);
    fingerprint.real(pressureSolve.maximumCellConservationResidual);
    const auto& traceSolve = pressureSolve.reducedTraceSolve;
    fingerprint.integer(static_cast<std::uint8_t>(traceSolve.compatible));
    fingerprint.integer(static_cast<std::uint8_t>(traceSolve.converged));
    fingerprint.integer(static_cast<std::uint8_t>(traceSolve.finite));
    fingerprint.integer(traceSolve.traceSystemFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(traceSolve.traceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        traceSolve.componentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        traceSolve.iterationCount));
    fingerprint.real(
        traceSolve.maximumAbsoluteComponentCompatibilityPascalsMeters);
    fingerprint.real(traceSolve.initialResidualL2PascalsMeters);
    fingerprint.real(traceSolve.finalResidualL2PascalsMeters);
    fingerprint.real(traceSolve.finalResidualMaximumPascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        traceSolve.components.size()));
    for (const auto& component : traceSolve.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.traceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeTraceIndex));
        fingerprint.real(component.rightHandSideSumPascalsMeters);
        fingerprint.real(component.compatibilityCorrectionPascalsMeters);
        fingerprint.real(component.traceGaugeBeforePascals);
        fingerprint.real(component.traceGaugeAfterPascals);
    }
    fingerprint.integer(
        endpoint.pressureEpoch.acceptedPressureState.fingerprint);
    fingerprint.integer(
        endpoint.pressureEpoch.acceptedPressureSamples.fingerprint);
    return fingerprint.value();
}

void validateSettings(
    const SceneFluidMimeticPressureAuditSettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit physical settings are invalid");
    }
}

void validateCommonInputIdentity(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux) {
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidOpeningFluxIntegrity(openingFlux);
    if (pressureEpoch.fingerprint == 0
        || pressureEpoch.surfaceStateFingerprint != state.fingerprint
        || pressureEpoch.openingPatches.fingerprint
            != openingFlux.openingPatchFingerprint
        || pressureEpoch.cellCounts != grid.cellCounts()
        || pressureEpoch.lowerMeters != grid.lowerMeters()
        || pressureEpoch.upperMeters != grid.upperMeters()
        || openingFlux.cellCounts != grid.cellCounts()
        || openingFlux.lowerMeters != grid.lowerMeters()
        || openingFlux.upperMeters != grid.upperMeters()
        || pressureEpoch.surfaceDefinitionFingerprint
            != openingFlux.surfaceDefinitionFingerprint
        || pressureEpoch.surfaceStateFingerprint
            != openingFlux.surfaceStateFingerprint
        || pressureEpoch.acceptedStepCount != openingFlux.acceptedStepCount
        || pressureEpoch.simulationTimeSeconds
            != openingFlux.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit common input identity is invalid");
    }
}

void validateInputIdentity(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition) {
    validateCommonInputIdentity(
        surface, state, grid, pressureEpoch, openingFlux);
    validateSceneFluidPressureVolumeRateIntegrity(geometryVolumeRates);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        topologyTransition);
    if (pressureEpoch.pressureControlVolumes.fingerprint
            != geometryVolumeRates.currentPressureControlVolumeFingerprint
        || pressureEpoch.pressureControlVolumes.fingerprint
            != topologyTransition.currentPressureControlVolumeFingerprint
        || pressureEpoch.pressureFaceLinks.fingerprint
            != topologyTransition.currentPressureFaceLinkFingerprint
        || pressureEpoch.cellVolumes.fingerprint
            != geometryVolumeRates.currentCellVolumeFingerprint
        || pressureEpoch.acceptedStepCount
            != geometryVolumeRates.currentAcceptedStepCount
        || pressureEpoch.simulationTimeSeconds
            != geometryVolumeRates.currentSimulationTimeSeconds
        || geometryVolumeRates.pressureTopologyTransitionFingerprint
            != topologyTransition.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit input identity is invalid");
    }
}

SceneFluidMimeticPressureAuditEndpoint finishEndpoint(
    SceneFluidMimeticPressureAuditEndpoint result,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    if (!result.pressureEpoch.diagnostics.accepted) {
        throw std::runtime_error(
            "scene fluid mimetic pressure-audit solve was not accepted");
    }
    result.ownedStorageBytes = checkedStorageSum(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure-audit endpoint exceeds its byte limit");
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticPressureAuditEndpointIntegrity(result);
    return result;
}

SceneFluidMimeticPressureAuditEndpoint buildEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& graphPressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField* predictedVelocityMetersPerSecond,
    const SceneFluidRegionWallExchange* wallExchange,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditEndpoint* previousEndpoint,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    validateSettings(settings);
    validateInputIdentity(
        surface, state, grid, graphPressureEpoch, openingFlux,
        geometryVolumeRates, topologyTransition);
    if ((predictedVelocityMetersPerSecond == nullptr)
        == (wallExchange == nullptr)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit predictor is invalid");
    }
    if (previousEndpoint != nullptr) {
        validateSceneFluidMimeticPressureAuditEndpointIntegrity(
            *previousEndpoint);
        if (previousEndpoint->acceptedStepCount
                != topologyTransition.previousAcceptedStepCount
            || previousEndpoint->simulationTimeSeconds
                != topologyTransition.previousSimulationTimeSeconds
            || previousEndpoint->controlCells
                    .pressureControlVolumeFingerprint
                != topologyTransition
                    .previousPressureControlVolumeFingerprint
            || previousEndpoint->controlCells.pressureFaceLinkFingerprint
                != topologyTransition.previousPressureFaceLinkFingerprint) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-audit previous endpoint is foreign");
        }
    }

    SceneFluidMimeticPressureAuditEndpoint result;
    result.scenePressureEpochFingerprint = graphPressureEpoch.fingerprint;
    result.pressureTopologyTransitionFingerprint =
        topologyTransition.fingerprint;
    result.structureDefinitionFingerprint =
        graphPressureEpoch.structureDefinitionFingerprint;
    result.acceptedStepCount = graphPressureEpoch.acceptedStepCount;
    result.simulationTimeSeconds = graphPressureEpoch.simulationTimeSeconds;
    result.usesRegionWallPrediction = wallExchange != nullptr;
    result.usesConsecutiveWarmStart = previousEndpoint != nullptr;
    result.controlCells = buildSceneFluidMimeticControlCells(
        surface, state, grid, graphPressureEpoch.gridEpoch,
        graphPressureEpoch.openingCaps,
        graphPressureEpoch.openingQuadrature,
        graphPressureEpoch.openingPatches,
        graphPressureEpoch.pressureControlVolumes,
        graphPressureEpoch.pressureFaceLinks, settings.controlCells,
        limits.controlCells);
    result.fullTraceSystem = buildSceneFluidMimeticTraceSystem(
        result.controlCells, settings.traceSystem, limits.traceSystem);
    result.condensedTraceSystem =
        buildSceneFluidMimeticCondensedTraceSystem(
            result.fullTraceSystem, limits.condensedSystem);
    if (wallExchange != nullptr) {
        result.predictedTraceFlows = sampleSceneFluidMimeticTraceFlows(
            result.controlCells, result.fullTraceSystem,
            graphPressureEpoch.pressureFaceLinks, openingFlux,
            *wallExchange, limits.traceFlow);
    } else {
        result.predictedTraceFlows = sampleSceneFluidMimeticTraceFlows(
            result.controlCells, result.fullTraceSystem,
            graphPressureEpoch.pressureFaceLinks, openingFlux, grid,
            *predictedVelocityMetersPerSecond, limits.traceFlow);
    }
    SceneFluidMimeticPressureSourceSettings sourceSettings;
    sourceSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    sourceSettings.timeStepSeconds = geometryVolumeRates.durationSeconds;
    result.pressureSources = buildSceneFluidMimeticPressureSources(
        result.controlCells, result.fullTraceSystem,
        result.predictedTraceFlows, geometryVolumeRates, sourceSettings,
        limits.pressureSource);
    if (previousEndpoint != nullptr) {
        result.pressureEpoch = acceptSceneFluidMimeticPressureEpoch(
            graphPressureEpoch.gridEpoch.quadrature,
            graphPressureEpoch.pressureControlVolumes,
            result.controlCells, result.fullTraceSystem,
            result.condensedTraceSystem, result.pressureSources,
            previousEndpoint->pressureEpoch.acceptedPressureState,
            previousEndpoint->controlCells,
            previousEndpoint->fullTraceSystem,
            previousEndpoint->condensedTraceSystem, topologyTransition,
            settings.pressureSolve, limits.pressureEpoch);
    } else {
        result.pressureEpoch = acceptSceneFluidMimeticPressureEpoch(
            graphPressureEpoch.gridEpoch.quadrature,
            graphPressureEpoch.pressureControlVolumes,
            result.controlCells, result.fullTraceSystem,
            result.condensedTraceSystem, result.pressureSources,
            settings.pressureSolve, limits.pressureEpoch);
    }
    return finishEndpoint(std::move(result), limits);
}

} // namespace

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& gridEpoch,
    const SceneFluidOpeningCapSet& openingCaps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& pressureFaceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    validateSettings(settings);
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidOpeningFluxIntegrity(openingFlux);
    if (openingFlux.surfaceDefinitionFingerprint != surface.fingerprint
        || openingFlux.surfaceStateFingerprint != state.fingerprint
        || openingFlux.openingPatchFingerprint != openingPatches.fingerprint
        || openingFlux.cellCounts != grid.cellCounts()
        || openingFlux.lowerMeters != grid.lowerMeters()
        || openingFlux.upperMeters != grid.upperMeters()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit fixed input identity is invalid");
    }
    SceneFluidMimeticPressureAuditEndpoint result;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.controlCells = buildSceneFluidMimeticControlCells(
        surface, state, grid, gridEpoch, openingCaps,
        openingQuadrature, openingPatches, pressureVolumes,
        pressureFaceLinks, settings.controlCells, limits.controlCells);
    result.fullTraceSystem = buildSceneFluidMimeticTraceSystem(
        result.controlCells, settings.traceSystem, limits.traceSystem);
    result.condensedTraceSystem =
        buildSceneFluidMimeticCondensedTraceSystem(
            result.fullTraceSystem, limits.condensedSystem);
    result.predictedTraceFlows = sampleSceneFluidMimeticTraceFlows(
        result.controlCells, result.fullTraceSystem, pressureFaceLinks,
        openingFlux, grid, predictedVelocityMetersPerSecond,
        limits.traceFlow);
    SceneFluidMimeticPressureSourceSettings sourceSettings;
    sourceSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    sourceSettings.timeStepSeconds = settings.timeStepSeconds;
    result.pressureSources = buildSceneFluidMimeticPressureSources(
        result.controlCells, result.fullTraceSystem,
        result.predictedTraceFlows, sourceSettings,
        limits.pressureSource);
    result.pressureEpoch = acceptSceneFluidMimeticPressureEpoch(
        gridEpoch.quadrature, pressureVolumes, result.controlCells,
        result.fullTraceSystem, result.condensedTraceSystem,
        result.pressureSources, settings.pressureSolve,
        limits.pressureEpoch);
    return finishEndpoint(std::move(result), limits);
}

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    validateSettings(settings);
    validateCommonInputIdentity(
        surface, state, grid, pressureEpoch, openingFlux);
    SceneFluidMimeticPressureAuditEndpoint result;
    result.scenePressureEpochFingerprint = pressureEpoch.fingerprint;
    result.structureDefinitionFingerprint =
        pressureEpoch.structureDefinitionFingerprint;
    result.acceptedStepCount = pressureEpoch.acceptedStepCount;
    result.simulationTimeSeconds = pressureEpoch.simulationTimeSeconds;
    result.controlCells = buildSceneFluidMimeticControlCells(
        surface, state, grid, pressureEpoch.gridEpoch,
        pressureEpoch.openingCaps, pressureEpoch.openingQuadrature,
        pressureEpoch.openingPatches, pressureEpoch.pressureControlVolumes,
        pressureEpoch.pressureFaceLinks, settings.controlCells,
        limits.controlCells);
    result.fullTraceSystem = buildSceneFluidMimeticTraceSystem(
        result.controlCells, settings.traceSystem, limits.traceSystem);
    result.condensedTraceSystem =
        buildSceneFluidMimeticCondensedTraceSystem(
            result.fullTraceSystem, limits.condensedSystem);
    result.predictedTraceFlows = sampleSceneFluidMimeticTraceFlows(
        result.controlCells, result.fullTraceSystem,
        pressureEpoch.pressureFaceLinks, openingFlux, grid,
        predictedVelocityMetersPerSecond, limits.traceFlow);
    SceneFluidMimeticPressureSourceSettings sourceSettings;
    sourceSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    sourceSettings.timeStepSeconds = settings.timeStepSeconds;
    result.pressureSources = buildSceneFluidMimeticPressureSources(
        result.controlCells, result.fullTraceSystem,
        result.predictedTraceFlows, sourceSettings,
        limits.pressureSource);
    result.pressureEpoch = acceptSceneFluidMimeticPressureEpoch(
        pressureEpoch.gridEpoch.quadrature,
        pressureEpoch.pressureControlVolumes, result.controlCells,
        result.fullTraceSystem, result.condensedTraceSystem,
        result.pressureSources, settings.pressureSolve,
        limits.pressureEpoch);
    return finishEndpoint(std::move(result), limits);
}

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    return buildEndpoint(
        surface, state, grid, pressureEpoch, openingFlux,
        &predictedVelocityMetersPerSecond, nullptr, geometryVolumeRates,
        topologyTransition, nullptr, settings, limits);
}

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditEndpoint& previousEndpoint,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    return buildEndpoint(
        surface, state, grid, pressureEpoch, openingFlux,
        &predictedVelocityMetersPerSecond, nullptr, geometryVolumeRates,
        topologyTransition, &previousEndpoint, settings, limits);
}

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    return buildEndpoint(
        surface, state, grid, pressureEpoch, openingFlux, nullptr,
        &wallExchange, geometryVolumeRates, topologyTransition, nullptr,
        settings, limits);
}

SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditEndpoint& previousEndpoint,
    const SceneFluidMimeticPressureAuditSettings& settings,
    const SceneFluidMimeticPressureAuditLimits& limits) {
    return buildEndpoint(
        surface, state, grid, pressureEpoch, openingFlux, nullptr,
        &wallExchange, geometryVolumeRates, topologyTransition,
        &previousEndpoint, settings, limits);
}

void validateSceneFluidMimeticPressureAuditEndpointIntegrity(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint) {
    validateSceneFluidMimeticControlCellIntegrity(endpoint.controlCells);
    validateSceneFluidMimeticTraceSystem(
        endpoint.fullTraceSystem, endpoint.controlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        endpoint.condensedTraceSystem, endpoint.fullTraceSystem);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        endpoint.predictedTraceFlows);
    validateSceneFluidMimeticPressureSources(
        endpoint.pressureSources, endpoint.controlCells);
    validateSceneFluidMimeticPressureEpochResultIntegrity(
        endpoint.pressureEpoch);
    if (endpoint.version != sceneFluidMimeticPressureAuditVersion
        || endpoint.fingerprint == 0
        || endpoint.structureDefinitionFingerprint == 0
        || !std::isfinite(endpoint.simulationTimeSeconds)
        || !endpoint.pressureEpoch.diagnostics.accepted
        || endpoint.usesConsecutiveWarmStart
            != endpoint.pressureEpoch.diagnostics.usedConsecutiveWarmStart
        || endpoint.pressureEpoch.topologyTransitionFingerprint
            != (endpoint.usesConsecutiveWarmStart
                    ? endpoint.pressureTopologyTransitionFingerprint : 0)
        || endpoint.controlCells.structureDefinitionFingerprint
            != endpoint.structureDefinitionFingerprint
        || endpoint.controlCells.acceptedStepCount
            != endpoint.acceptedStepCount
        || endpoint.controlCells.simulationTimeSeconds
            != endpoint.simulationTimeSeconds
        || endpoint.fullTraceSystem.fingerprint
            != endpoint.pressureEpoch.fullTraceSystemFingerprint
        || endpoint.condensedTraceSystem.fingerprint
            != endpoint.pressureEpoch.condensedTraceSystemFingerprint
        || endpoint.pressureSources.fingerprint
            != endpoint.pressureEpoch.pressureSourceFingerprint
        || endpoint.predictedTraceFlows.fingerprint
            != endpoint.pressureSources.mimeticTraceFlowFingerprint
        || endpoint.usesRegionWallPrediction
            != (endpoint.predictedTraceFlows
                    .regionWallExchangeFingerprint != 0)
        || endpoint.ownedStorageBytes != checkedStorageSum(endpoint)
        || endpoint.fingerprint != productFingerprint(endpoint)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-audit endpoint integrity is invalid");
    }
}

} // namespace simwing::fsi
