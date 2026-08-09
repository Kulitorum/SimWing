#include "scene_fluid_mimetic_pressure_warm_start.h"

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
            "scene fluid mimetic pressure warm-start storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t count,
                            const std::size_t elementBytes) {
    if (count > std::numeric_limits<std::size_t>::max() / elementBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure warm-start storage overflows");
    }
    return count * elementBytes;
}

std::size_t ownedBytesForCounts(const std::size_t traceCount,
                                const std::size_t componentCount) {
    return checkedAdd(
        checkedMultiply(traceCount, sizeof(double)),
        checkedMultiply(componentCount, sizeof(double)));
}

std::size_t workingBytesForCounts(const std::size_t controlCount,
                                  const std::size_t traceCount,
                                  const std::size_t componentCount) {
    std::size_t result = ownedBytesForCounts(traceCount, componentCount);
    result = checkedAdd(
        result, checkedMultiply(controlCount, 2 * sizeof(double)));
    return checkedAdd(
        result, checkedMultiply(controlCount, sizeof(std::uint8_t)));
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticPressureWarmStart& warmStart) {
    Fingerprint fingerprint;
    fingerprint.integer(warmStart.version);
    fingerprint.integer(warmStart.sourcePressureStateFingerprint);
    fingerprint.integer(warmStart.sourceTopologyTransitionFingerprint);
    fingerprint.integer(warmStart.previousMimeticControlCellFingerprint);
    fingerprint.integer(warmStart.previousFullTraceSystemFingerprint);
    fingerprint.integer(
        warmStart.previousCondensedTraceSystemFingerprint);
    fingerprint.integer(warmStart.currentMimeticControlCellFingerprint);
    fingerprint.integer(warmStart.currentFullTraceSystemFingerprint);
    fingerprint.integer(warmStart.currentCondensedTraceSystemFingerprint);
    fingerprint.integer(warmStart.structureDefinitionFingerprint);
    fingerprint.integer(warmStart.previousAcceptedStepCount);
    fingerprint.integer(warmStart.currentAcceptedStepCount);
    fingerprint.real(warmStart.previousSimulationTimeSeconds);
    fingerprint.real(warmStart.currentSimulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.workingStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.componentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.previousReducedTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.currentReducedTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.retainedTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.appearedTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.disappearedTraceCount));
    fingerprint.real(warmStart.maximumAbsoluteGaugeShiftPascals);
    fingerprint.real(warmStart.maximumAbsolutePressurePascals);
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.componentGaugeShiftsPascals.size()));
    for (const double value : warmStart.componentGaugeShiftsPascals) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        warmStart.reducedTracePascals.size()));
    for (const double value : warmStart.reducedTracePascals) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

void validateInputs(
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition) {
    validateSceneFluidMimeticPressureState(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem);
    validateSceneFluidMimeticTraceSystem(
        currentFullSystem, currentControlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        currentCondensedSystem, currentFullSystem);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        topologyTransition);
    if (previousControlCells.pressureControlVolumeFingerprint
            != topologyTransition
                .previousPressureControlVolumeFingerprint
        || previousControlCells.pressureFaceLinkFingerprint
            != topologyTransition.previousPressureFaceLinkFingerprint
        || currentControlCells.pressureControlVolumeFingerprint
            != topologyTransition.currentPressureControlVolumeFingerprint
        || currentControlCells.pressureFaceLinkFingerprint
            != topologyTransition.currentPressureFaceLinkFingerprint
        || previousControlCells.surfaceDefinitionFingerprint
            != topologyTransition.surfaceDefinitionFingerprint
        || currentControlCells.surfaceDefinitionFingerprint
            != topologyTransition.surfaceDefinitionFingerprint
        || previousControlCells.structureDefinitionFingerprint
            != topologyTransition.structureDefinitionFingerprint
        || currentControlCells.structureDefinitionFingerprint
            != topologyTransition.structureDefinitionFingerprint
        || previousControlCells.acceptedStepCount
            != topologyTransition.previousAcceptedStepCount
        || currentControlCells.acceptedStepCount
            != topologyTransition.currentAcceptedStepCount
        || previousControlCells.simulationTimeSeconds
            != topologyTransition.previousSimulationTimeSeconds
        || currentControlCells.simulationTimeSeconds
            != topologyTransition.currentSimulationTimeSeconds
        || previousControlCells.cellCounts != topologyTransition.cellCounts
        || currentControlCells.cellCounts != topologyTransition.cellCounts
        || previousControlCells.lowerMeters != topologyTransition.lowerMeters
        || currentControlCells.lowerMeters != topologyTransition.lowerMeters
        || previousControlCells.upperMeters != topologyTransition.upperMeters
        || currentControlCells.upperMeters != topologyTransition.upperMeters
        || previousControlCells.controlCells.size()
            != topologyTransition.previousControlVolumeCount
        || currentControlCells.controlCells.size()
            != topologyTransition.currentControlVolumeCount) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure warm-start identity is invalid");
    }
}

SceneFluidMimeticPressureWarmStart buildWarmStart(
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureWarmStartLimits& limits) {
    validateInputs(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem, currentControlCells, currentFullSystem,
        currentCondensedSystem, topologyTransition);
    if (previousControlCells.controlCells.size()
            > limits.maximumControlCells
        || currentControlCells.controlCells.size()
            > limits.maximumControlCells
        || previousCondensedSystem.traces.size()
            > limits.maximumPreviousReducedTraces
        || currentCondensedSystem.traces.size()
            > limits.maximumCurrentReducedTraces
        || currentCondensedSystem.componentCount
            > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid mimetic pressure warm-start count limit exceeded");
    }
    const std::size_t ownedBytes = ownedBytesForCounts(
        currentCondensedSystem.traces.size(),
        currentCondensedSystem.componentCount);
    const std::size_t workingBytes = workingBytesForCounts(
        currentControlCells.controlCells.size(),
        currentCondensedSystem.traces.size(),
        currentCondensedSystem.componentCount);
    if (ownedBytes > limits.maximumOwnedBytes
        || workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure warm-start byte limit exceeded");
    }

    std::vector<double> rebasedControls(
        currentControlCells.controlCells.size(), 0.0);
    std::vector<double> appearanceDonorAreas(
        currentControlCells.controlCells.size(), 0.0);
    std::vector<std::uint8_t> initialized(
        currentControlCells.controlCells.size(), 0);
    for (const auto& retained : topologyTransition.retainedControls) {
        const auto& previous = previousState.controls[
            retained.previousControlVolumeIndex];
        const auto& current = currentControlCells.controlCells[
            retained.currentControlVolumeIndex];
        if (previous.stableId != retained.stableId
            || current.stableId != retained.stableId) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure warm-start retained control is invalid");
        }
        rebasedControls[retained.currentControlVolumeIndex] =
            previous.pressurePascals;
        initialized[retained.currentControlVolumeIndex] = 1;
    }
    for (const auto& donor : topologyTransition.appearanceDonors) {
        rebasedControls[donor.appearedCurrentControlVolumeIndex] +=
            donor.linkAreaSquareMeters
            * previousState.controls[
                donor.retainedPreviousControlVolumeIndex]
                  .pressurePascals;
        appearanceDonorAreas[
            donor.appearedCurrentControlVolumeIndex] +=
                donor.linkAreaSquareMeters;
        initialized[donor.appearedCurrentControlVolumeIndex] = 1;
    }
    for (std::size_t index = 0; index < rebasedControls.size(); ++index) {
        if (appearanceDonorAreas[index] > 0.0) {
            rebasedControls[index] /= appearanceDonorAreas[index];
        }
        if (initialized[index] == 0
            || !std::isfinite(rebasedControls[index])) {
            throw std::overflow_error(
                "scene fluid mimetic pressure warm-start control remap is invalid");
        }
    }

    SceneFluidMimeticPressureWarmStart result;
    result.sourcePressureStateFingerprint = previousState.fingerprint;
    result.sourceTopologyTransitionFingerprint =
        topologyTransition.fingerprint;
    result.previousMimeticControlCellFingerprint =
        previousControlCells.fingerprint;
    result.previousFullTraceSystemFingerprint =
        previousFullSystem.fingerprint;
    result.previousCondensedTraceSystemFingerprint =
        previousCondensedSystem.fingerprint;
    result.currentMimeticControlCellFingerprint =
        currentControlCells.fingerprint;
    result.currentFullTraceSystemFingerprint = currentFullSystem.fingerprint;
    result.currentCondensedTraceSystemFingerprint =
        currentCondensedSystem.fingerprint;
    result.structureDefinitionFingerprint =
        currentControlCells.structureDefinitionFingerprint;
    result.previousAcceptedStepCount =
        previousControlCells.acceptedStepCount;
    result.currentAcceptedStepCount = currentControlCells.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousControlCells.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentControlCells.simulationTimeSeconds;
    result.ownedStorageBytes = ownedBytes;
    result.workingStorageBytes = workingBytes;
    result.componentCount = currentCondensedSystem.componentCount;
    result.previousReducedTraceCount = previousState.traces.size();
    result.currentReducedTraceCount = currentCondensedSystem.traces.size();
    result.componentGaugeShiftsPascals.assign(result.componentCount, 0.0);
    result.reducedTracePascals.assign(
        result.currentReducedTraceCount, 0.0);

    std::size_t previousTraceIndex = 0;
    for (const auto& trace : currentCondensedSystem.traces) {
        while (previousTraceIndex < previousState.traces.size()
               && previousState.traces[previousTraceIndex].stableId
                   < trace.stableId) {
            ++previousTraceIndex;
        }
        double pressure = 0.0;
        if (previousTraceIndex < previousState.traces.size()
            && previousState.traces[previousTraceIndex].stableId
                == trace.stableId) {
            if (previousState.traces[previousTraceIndex].kind
                != trace.kind) {
                throw std::invalid_argument(
                    "scene fluid mimetic pressure warm-start retained trace kind is invalid");
            }
            pressure = previousState.traces[
                previousTraceIndex].pressurePascals;
            ++result.retainedTraceCount;
        } else {
            const auto& fullTrace =
                currentFullSystem.traces[trace.fullTraceIndex];
            if (fullTrace.incidenceCount == 0) {
                throw std::invalid_argument(
                    "scene fluid mimetic pressure warm-start appeared trace has no endpoint");
            }
            for (std::size_t offset = 0;
                 offset < fullTrace.incidenceCount; ++offset) {
                const auto& incidence = currentFullSystem.incidences[
                    fullTrace.firstIncidence + offset];
                pressure += rebasedControls[incidence.controlCellIndex];
            }
            pressure /= static_cast<double>(fullTrace.incidenceCount);
            ++result.appearedTraceCount;
        }
        if (!std::isfinite(pressure)) {
            throw std::overflow_error(
                "scene fluid mimetic pressure warm-start trace is non-finite");
        }
        result.reducedTracePascals[trace.traceIndex] = pressure;
    }
    result.disappearedTraceCount = previousState.traces.size()
        - result.retainedTraceCount;

    for (std::size_t component = 0;
         component < result.componentCount; ++component) {
        const std::size_t gauge =
            currentCondensedSystem.componentGaugeTraceIndices[component];
        const double shift = result.reducedTracePascals[gauge];
        result.componentGaugeShiftsPascals[component] = shift;
        result.maximumAbsoluteGaugeShiftPascals = std::max(
            result.maximumAbsoluteGaugeShiftPascals, std::abs(shift));
    }
    for (const auto& trace : currentCondensedSystem.traces) {
        result.reducedTracePascals[trace.traceIndex] -=
            result.componentGaugeShiftsPascals[trace.componentIndex];
        result.maximumAbsolutePressurePascals = std::max(
            result.maximumAbsolutePressurePascals,
            std::abs(result.reducedTracePascals[trace.traceIndex]));
    }
    for (const std::size_t gauge :
         currentCondensedSystem.componentGaugeTraceIndices) {
        result.reducedTracePascals[gauge] = 0.0;
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticPressureWarmStartIntegrity(result);
    return result;
}

} // namespace

SceneFluidMimeticPressureWarmStart
buildSceneFluidMimeticPressureWarmStart(
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureWarmStartLimits& limits) {
    return buildWarmStart(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem, currentControlCells, currentFullSystem,
        currentCondensedSystem, topologyTransition, limits);
}

void validateSceneFluidMimeticPressureWarmStartIntegrity(
    const SceneFluidMimeticPressureWarmStart& warmStart) {
    if (warmStart.version
            != sceneFluidMimeticPressureWarmStartVersion
        || warmStart.fingerprint == 0
        || warmStart.sourcePressureStateFingerprint == 0
        || warmStart.sourceTopologyTransitionFingerprint == 0
        || warmStart.previousMimeticControlCellFingerprint == 0
        || warmStart.previousFullTraceSystemFingerprint == 0
        || warmStart.previousCondensedTraceSystemFingerprint == 0
        || warmStart.currentMimeticControlCellFingerprint == 0
        || warmStart.currentFullTraceSystemFingerprint == 0
        || warmStart.currentCondensedTraceSystemFingerprint == 0
        || warmStart.structureDefinitionFingerprint == 0
        || !std::isfinite(warmStart.previousSimulationTimeSeconds)
        || !std::isfinite(warmStart.currentSimulationTimeSeconds)
        || warmStart.previousAcceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || warmStart.currentAcceptedStepCount
            != warmStart.previousAcceptedStepCount + 1
        || !(warmStart.currentSimulationTimeSeconds
            > warmStart.previousSimulationTimeSeconds)
        || warmStart.componentCount == 0
        || warmStart.previousReducedTraceCount == 0
        || warmStart.currentReducedTraceCount == 0
        || warmStart.componentGaugeShiftsPascals.size()
            != warmStart.componentCount
        || warmStart.reducedTracePascals.size()
            != warmStart.currentReducedTraceCount
        || warmStart.retainedTraceCount
                + warmStart.appearedTraceCount
            != warmStart.currentReducedTraceCount
        || warmStart.retainedTraceCount
                + warmStart.disappearedTraceCount
            != warmStart.previousReducedTraceCount
        || warmStart.ownedStorageBytes != ownedBytesForCounts(
            warmStart.currentReducedTraceCount, warmStart.componentCount)
        || warmStart.workingStorageBytes < warmStart.ownedStorageBytes) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure warm-start integrity is invalid");
    }
    double maximumShift = 0.0;
    for (const double shift : warmStart.componentGaugeShiftsPascals) {
        if (!std::isfinite(shift)) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure warm-start gauge shift is invalid");
        }
        maximumShift = std::max(maximumShift, std::abs(shift));
    }
    double maximumPressure = 0.0;
    for (const double pressure : warmStart.reducedTracePascals) {
        if (!std::isfinite(pressure)) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure warm-start pressure is invalid");
        }
        maximumPressure = std::max(maximumPressure, std::abs(pressure));
    }
    if (warmStart.maximumAbsoluteGaugeShiftPascals != maximumShift
        || warmStart.maximumAbsolutePressurePascals != maximumPressure
        || productFingerprint(warmStart) != warmStart.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure warm-start summary is invalid");
    }
}

void validateSceneFluidMimeticPressureWarmStart(
    const SceneFluidMimeticPressureWarmStart& warmStart,
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition) {
    validateSceneFluidMimeticPressureWarmStartIntegrity(warmStart);
    validateInputs(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem, currentControlCells, currentFullSystem,
        currentCondensedSystem, topologyTransition);
    const SceneFluidMimeticPressureWarmStartLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildWarmStart(
        previousState, previousControlCells, previousFullSystem,
        previousCondensedSystem, currentControlCells, currentFullSystem,
        currentCondensedSystem, topologyTransition, unlimited);
    if (warmStart != expected) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure warm-start payload is invalid");
    }
}

} // namespace simwing::fsi
