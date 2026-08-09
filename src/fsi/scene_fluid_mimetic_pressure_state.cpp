#include "scene_fluid_mimetic_pressure_state.h"

#include <algorithm>
#include <bit>
#include <cmath>
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

std::size_t storageBytesForCounts(const std::size_t controlCount,
                                  const std::size_t traceCount) {
    if (controlCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticAcceptedControlPressure)
        || traceCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticAcceptedTracePressure)) {
        throw std::length_error(
            "scene fluid mimetic pressure-state storage overflows");
    }
    const std::size_t controlBytes = controlCount
        * sizeof(SceneFluidMimeticAcceptedControlPressure);
    const std::size_t traceBytes = traceCount
        * sizeof(SceneFluidMimeticAcceptedTracePressure);
    if (traceBytes > std::numeric_limits<std::size_t>::max() - controlBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure-state storage overflows");
    }
    return controlBytes + traceBytes;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticPressureState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.mimeticControlCellFingerprint);
    fingerprint.integer(state.fullTraceSystemFingerprint);
    fingerprint.integer(state.condensedTraceSystemFingerprint);
    fingerprint.integer(state.pressureSourceFingerprint);
    fingerprint.integer(state.structureDefinitionFingerprint);
    fingerprint.integer(state.acceptedStepCount);
    fingerprint.real(state.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(state.componentCount));
    fingerprint.real(state.maximumAbsoluteControlPressurePascals);
    fingerprint.real(state.maximumAbsoluteTracePressurePascals);
    fingerprint.integer(static_cast<std::uint64_t>(state.controls.size()));
    for (const auto& control : state.controls) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(control.cellIndex));
        fingerprint.integer(control.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.real(control.pressurePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.traces.size()));
    for (const auto& trace : state.traces) {
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.reducedTraceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.fullTraceIndex));
        fingerprint.integer(trace.stableId);
        fingerprint.enumeration(trace.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.componentIndex));
        fingerprint.integer(static_cast<std::uint8_t>(trace.isGauge));
        fingerprint.real(trace.pressurePascals);
    }
    return fingerprint.value();
}

bool evaluationMatches(
    const SceneFluidMimeticTraceEvaluation& first,
    const SceneFluidMimeticTraceEvaluation& second) {
    return first.cellScalars == second.cellScalars
        && first.halfFaceIntegratedOutwardFluxes
            == second.halfFaceIntegratedOutwardFluxes
        && first.traceIntegratedOutwardFluxSums
            == second.traceIntegratedOutwardFluxSums
        && first.maximumCellConservationResidual
            == second.maximumCellConservationResidual
        && first.maximumTraceFluxImbalance
            == second.maximumTraceFluxImbalance;
}

} // namespace

SceneFluidMimeticPressureState captureSceneFluidMimeticPressureState(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticPressureSolveResult& acceptedPressure,
    const SceneFluidMimeticPressureStateLimits& limits) {
    validateSceneFluidMimeticTraceSystem(fullSystem, controlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedSystem, fullSystem);
    validateSceneFluidMimeticPressureSources(sources, controlCells);
    if (!acceptedPressure.diagnostics.accepted
        || !acceptedPressure.diagnostics.reconstructedFullResidualConverged
        || !acceptedPressure.diagnostics.reducedTraceSolve.compatible
        || !acceptedPressure.diagnostics.reducedTraceSolve.converged
        || !acceptedPressure.diagnostics.reducedTraceSolve.finite
        || acceptedPressure.fullTraceSystemFingerprint
            != fullSystem.fingerprint
        || acceptedPressure.condensedTraceSystemFingerprint
            != condensedSystem.fingerprint
        || acceptedPressure.pressureSourceFingerprint != sources.fingerprint
        || acceptedPressure.reducedTracePascals.size()
            != condensedSystem.traces.size()
        || acceptedPressure.fullTracePascals.size()
            != fullSystem.traces.size()
        || acceptedPressure.evaluation.cellScalars.size()
            != controlCells.controlCells.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-state source is not accepted");
    }
    if (controlCells.controlCells.size() > limits.maximumControlCells
        || condensedSystem.traces.size() > limits.maximumReducedTraces) {
        throw std::length_error(
            "scene fluid mimetic pressure-state count limit exceeded");
    }
    const std::size_t expectedBytes = storageBytesForCounts(
        controlCells.controlCells.size(), condensedSystem.traces.size());
    if (expectedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure-state byte limit exceeded");
    }

    const auto integratedSources =
        sceneFluidMimeticIntegratedCellSources(sources);
    const auto fullRightHandSide =
        buildSceneFluidMimeticTraceRightHandSide(
            fullSystem, integratedSources);
    const auto reconstructed = reconstructSceneFluidMimeticFullTraces(
        condensedSystem, fullSystem, fullRightHandSide,
        acceptedPressure.reducedTracePascals);
    const auto reevaluated = evaluateSceneFluidMimeticTraceSystem(
        fullSystem, reconstructed, integratedSources);
    if (reconstructed != acceptedPressure.fullTracePascals
        || !evaluationMatches(reevaluated, acceptedPressure.evaluation)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-state reconstruction is invalid");
    }
    for (const std::size_t gauge :
         condensedSystem.componentGaugeTraceIndices) {
        if (acceptedPressure.reducedTracePascals[gauge] != 0.0) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-state gauge is not normalized");
        }
    }

    SceneFluidMimeticPressureState result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.fullTraceSystemFingerprint = fullSystem.fingerprint;
    result.condensedTraceSystemFingerprint = condensedSystem.fingerprint;
    result.pressureSourceFingerprint = sources.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.ownedStorageBytes = expectedBytes;
    result.componentCount = condensedSystem.componentCount;
    result.controls.reserve(controlCells.controlCells.size());
    for (const auto& cell : controlCells.controlCells) {
        const double pressure = acceptedPressure.evaluation.cellScalars[
            cell.controlCellIndex];
        if (!std::isfinite(pressure)) {
            throw std::overflow_error(
                "scene fluid mimetic pressure-state control is non-finite");
        }
        result.controls.push_back({
            cell.controlCellIndex,
            cell.controlVolumeIndex,
            cell.stableId,
            cell.cellIndex,
            cell.regionId,
            cell.componentIndex,
            pressure,
        });
        result.maximumAbsoluteControlPressurePascals = std::max(
            result.maximumAbsoluteControlPressurePascals,
            std::abs(pressure));
    }
    result.traces.reserve(condensedSystem.traces.size());
    for (const auto& trace : condensedSystem.traces) {
        const double pressure = acceptedPressure.reducedTracePascals[
            trace.traceIndex];
        if (!std::isfinite(pressure)) {
            throw std::overflow_error(
                "scene fluid mimetic pressure-state trace is non-finite");
        }
        result.traces.push_back({
            trace.traceIndex,
            trace.fullTraceIndex,
            trace.stableId,
            trace.kind,
            trace.componentIndex,
            trace.isGauge,
            pressure,
        });
        result.maximumAbsoluteTracePressurePascals = std::max(
            result.maximumAbsoluteTracePressurePascals,
            std::abs(pressure));
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticPressureState(
        result, controlCells, fullSystem, condensedSystem);
    return result;
}

void validateSceneFluidMimeticPressureStateIntegrity(
    const SceneFluidMimeticPressureState& state) {
    if (state.version != sceneFluidMimeticPressureStateVersion
        || state.fingerprint == 0
        || state.mimeticControlCellFingerprint == 0
        || state.fullTraceSystemFingerprint == 0
        || state.condensedTraceSystemFingerprint == 0
        || state.pressureSourceFingerprint == 0
        || state.structureDefinitionFingerprint == 0
        || !std::isfinite(state.simulationTimeSeconds)
        || state.componentCount == 0
        || state.controls.empty() || state.traces.empty()
        || state.ownedStorageBytes != storageBytesForCounts(
            state.controls.size(), state.traces.size())) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-state integrity is invalid");
    }
    double maximumControl = 0.0;
    for (const auto& control : state.controls) {
        if (control.controlCellIndex != &control - state.controls.data()
            || control.stableId == 0 || control.regionId == invalidStableId
            || control.componentIndex >= state.componentCount
            || !std::isfinite(control.pressurePascals)) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-state control is invalid");
        }
        maximumControl = std::max(
            maximumControl, std::abs(control.pressurePascals));
    }
    double maximumTrace = 0.0;
    for (const auto& trace : state.traces) {
        if (trace.reducedTraceIndex != &trace - state.traces.data()
            || trace.stableId == 0
            || trace.componentIndex >= state.componentCount
            || (trace.kind
                    != SceneFluidMimeticHalfFaceKind::CartesianTrace
                && trace.kind
                    != SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace)
            || !std::isfinite(trace.pressurePascals)
            || (trace.isGauge && trace.pressurePascals != 0.0)) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-state trace is invalid");
        }
        maximumTrace = std::max(
            maximumTrace, std::abs(trace.pressurePascals));
    }
    if (state.maximumAbsoluteControlPressurePascals != maximumControl
        || state.maximumAbsoluteTracePressurePascals != maximumTrace
        || productFingerprint(state) != state.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-state summary is invalid");
    }
}

void validateSceneFluidMimeticPressureState(
    const SceneFluidMimeticPressureState& state,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem) {
    validateSceneFluidMimeticTraceSystem(fullSystem, controlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedSystem, fullSystem);
    validateSceneFluidMimeticPressureStateIntegrity(state);
    if (state.mimeticControlCellFingerprint != controlCells.fingerprint
        || state.fullTraceSystemFingerprint != fullSystem.fingerprint
        || state.condensedTraceSystemFingerprint
            != condensedSystem.fingerprint
        || state.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || state.acceptedStepCount != controlCells.acceptedStepCount
        || state.simulationTimeSeconds != controlCells.simulationTimeSeconds
        || state.componentCount != condensedSystem.componentCount
        || state.controls.size() != controlCells.controlCells.size()
        || state.traces.size() != condensedSystem.traces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-state identity is invalid");
    }
    for (std::size_t index = 0; index < state.controls.size(); ++index) {
        const auto& value = state.controls[index];
        const auto& cell = controlCells.controlCells[index];
        if (value.controlCellIndex != cell.controlCellIndex
            || value.controlVolumeIndex != cell.controlVolumeIndex
            || value.stableId != cell.stableId
            || value.cellIndex != cell.cellIndex
            || value.regionId != cell.regionId
            || value.componentIndex != cell.componentIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-state control topology is invalid");
        }
    }
    for (std::size_t index = 0; index < state.traces.size(); ++index) {
        const auto& value = state.traces[index];
        const auto& trace = condensedSystem.traces[index];
        if (value.reducedTraceIndex != trace.traceIndex
            || value.fullTraceIndex != trace.fullTraceIndex
            || value.stableId != trace.stableId
            || value.kind != trace.kind
            || value.componentIndex != trace.componentIndex
            || value.isGauge != trace.isGauge) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-state trace topology is invalid");
        }
    }
}

} // namespace simwing::fsi
