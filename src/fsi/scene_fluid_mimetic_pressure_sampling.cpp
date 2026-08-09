#include "scene_fluid_mimetic_pressure_sampling.h"

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

std::size_t storageBytesForCount(const std::size_t sampleCount) {
    if (sampleCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidPressureSampleBinding)
        || sampleCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidQuadraturePressure)) {
        throw std::length_error(
            "scene fluid mimetic pressure-sampling storage overflows");
    }
    const std::size_t bindingBytes = sampleCount
        * sizeof(SceneFluidPressureSampleBinding);
    const std::size_t pressureBytes = sampleCount
        * sizeof(SceneFluidQuadraturePressure);
    if (pressureBytes
        > std::numeric_limits<std::size_t>::max() - bindingBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure-sampling storage overflows");
    }
    return bindingBytes + pressureBytes;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticPressureSampleSet& samples) {
    Fingerprint fingerprint;
    fingerprint.integer(samples.version);
    fingerprint.integer(samples.quadratureFingerprint);
    fingerprint.integer(samples.pressureControlVolumeFingerprint);
    fingerprint.integer(samples.mimeticControlCellFingerprint);
    fingerprint.integer(samples.fullTraceSystemFingerprint);
    fingerprint.integer(samples.condensedTraceSystemFingerprint);
    fingerprint.integer(samples.pressureStateFingerprint);
    fingerprint.integer(samples.surfaceDefinitionFingerprint);
    fingerprint.integer(samples.structureDefinitionFingerprint);
    fingerprint.integer(samples.acceptedStepCount);
    fingerprint.real(samples.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.componentCount));
    fingerprint.real(samples.maximumAbsolutePressureDifferencePascals);
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.bindings.size()));
    for (const auto& binding : samples.bindings) {
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.sampleIndex));
        fingerprint.integer(binding.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.negativeSideControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.positiveSideControlVolumeIndex));
        fingerprint.integer(binding.negativeSideControlVolumeStableId);
        fingerprint.integer(binding.positiveSideControlVolumeStableId);
        fingerprint.integer(binding.negativeSideRegionId);
        fingerprint.integer(binding.positiveSideRegionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.componentIndex));
        fingerprint.real(binding.pressureDifferencePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.pressures.size()));
    for (const auto& pressure : samples.pressures) {
        fingerprint.integer(pressure.stableId);
        fingerprint.real(pressure.negativeSidePressurePascals);
        fingerprint.real(pressure.positiveSidePressurePascals);
    }
    return fingerprint.value();
}

void validateSources(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState) {
    validateSceneFluidQuadratureDefinition(quadrature);
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidMimeticPressureState(
        pressureState, controlCells, fullSystem, condensedSystem);
    if (controlCells.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || quadrature.surfaceDefinitionFingerprint
            != pressureVolumes.surfaceDefinitionFingerprint
        || quadrature.surfaceDefinitionFingerprint
            != controlCells.surfaceDefinitionFingerprint
        || quadrature.structureDefinitionFingerprint
            != pressureVolumes.structureDefinitionFingerprint
        || quadrature.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || quadrature.acceptedStepCount
            != pressureVolumes.acceptedStepCount
        || quadrature.acceptedStepCount
            != controlCells.acceptedStepCount
        || quadrature.simulationTimeSeconds
            != pressureVolumes.simulationTimeSeconds
        || quadrature.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || pressureVolumes.controlVolumes.size()
            != controlCells.controlCells.size()
        || pressureVolumes.controlVolumes.size()
            != pressureState.controls.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-sampling source identity is invalid");
    }
    for (std::size_t index = 0;
         index < pressureVolumes.controlVolumes.size(); ++index) {
        const auto& volume = pressureVolumes.controlVolumes[index];
        const auto& control = pressureState.controls[index];
        if (control.controlVolumeIndex != index
            || control.stableId != volume.stableId
            || control.cellIndex != volume.cellIndex
            || control.regionId != volume.regionId
            || control.componentIndex != volume.componentIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-sampling control topology is invalid");
        }
    }
}

const SceneFluidPressureControlVolume& findControl(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::size_t cellIndex,
    const StableId regionId) {
    if (cellIndex >= pressureVolumes.cells.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure sample is outside the pressure grid");
    }
    const auto& cell = pressureVolumes.cells[cellIndex];
    const SceneFluidPressureControlVolume* found = nullptr;
    for (std::size_t offset = 0; offset < cell.controlVolumeCount; ++offset) {
        const auto& control = pressureVolumes.controlVolumes[
            cell.firstControlVolume + offset];
        if (control.regionId == regionId) {
            if (found != nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic pressure sample has duplicate controls");
            }
            found = &control;
        }
    }
    if (found == nullptr) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure sample has no matching control");
    }
    return *found;
}

SceneFluidMimeticPressureSampleSet buildSamples(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState,
    const SceneFluidPressureSamplingLimits& limits) {
    if (quadrature.points.size() > limits.maximumSamples) {
        throw std::length_error(
            "scene fluid mimetic pressure sampling exceeds its sample limit");
    }
    const std::size_t expectedBytes = storageBytesForCount(
        quadrature.points.size());
    if (expectedBytes > limits.maximumSamplingBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure sampling exceeds its byte limit");
    }

    SceneFluidMimeticPressureSampleSet result;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.fullTraceSystemFingerprint = fullSystem.fingerprint;
    result.condensedTraceSystemFingerprint = condensedSystem.fingerprint;
    result.pressureStateFingerprint = pressureState.fingerprint;
    result.surfaceDefinitionFingerprint =
        quadrature.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        quadrature.structureDefinitionFingerprint;
    result.acceptedStepCount = quadrature.acceptedStepCount;
    result.simulationTimeSeconds = quadrature.simulationTimeSeconds;
    result.controlVolumeCount = pressureVolumes.controlVolumes.size();
    result.componentCount = condensedSystem.componentCount;
    result.bindings.reserve(quadrature.points.size());
    result.pressures.reserve(quadrature.points.size());
    for (std::size_t index = 0; index < quadrature.points.size(); ++index) {
        const auto& point = quadrature.points[index];
        const auto& negative = findControl(
            pressureVolumes, point.negativeSideCellIndex,
            point.negativeSideRegionId);
        const auto& positive = findControl(
            pressureVolumes, point.positiveSideCellIndex,
            point.positiveSideRegionId);
        if (negative.componentIndex != positive.componentIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic sheet sides use independent pressure gauges");
        }
        const double negativePressure = pressureState.controls[
            negative.controlVolumeIndex].pressurePascals;
        const double positivePressure = pressureState.controls[
            positive.controlVolumeIndex].pressurePascals;
        const double difference = negativePressure - positivePressure;
        if (!std::isfinite(difference)) {
            throw std::overflow_error(
                "scene fluid mimetic pressure sample is non-finite");
        }
        result.maximumAbsolutePressureDifferencePascals = std::max(
            result.maximumAbsolutePressureDifferencePascals,
            std::abs(difference));
        result.bindings.push_back({
            index,
            point.stableId,
            negative.controlVolumeIndex,
            positive.controlVolumeIndex,
            negative.stableId,
            positive.stableId,
            negative.regionId,
            positive.regionId,
            negative.componentIndex,
            difference,
        });
        result.pressures.push_back({
            point.stableId, negativePressure, positivePressure,
        });
    }
    result.ownedStorageBytes = expectedBytes;
    result.fingerprint = productFingerprint(result);
    return result;
}

} // namespace

SceneFluidMimeticPressureSampleSet sampleSceneFluidMimeticPressure(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState,
    const SceneFluidPressureSamplingLimits& limits) {
    validateSources(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, pressureState);
    auto result = buildSamples(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, pressureState, limits);
    validateSceneFluidMimeticPressureSamples(
        result, quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, pressureState);
    return result;
}

void validateSceneFluidMimeticPressureSampleIntegrity(
    const SceneFluidMimeticPressureSampleSet& samples) {
    if (samples.version != sceneFluidMimeticPressureSamplingVersion
        || samples.fingerprint == 0
        || samples.quadratureFingerprint == 0
        || samples.pressureControlVolumeFingerprint == 0
        || samples.mimeticControlCellFingerprint == 0
        || samples.fullTraceSystemFingerprint == 0
        || samples.condensedTraceSystemFingerprint == 0
        || samples.pressureStateFingerprint == 0
        || samples.surfaceDefinitionFingerprint == 0
        || samples.structureDefinitionFingerprint == 0
        || !std::isfinite(samples.simulationTimeSeconds)
        || samples.controlVolumeCount == 0
        || samples.componentCount == 0
        || samples.bindings.empty()
        || samples.bindings.size() != samples.pressures.size()
        || samples.ownedStorageBytes
            != storageBytesForCount(samples.bindings.size())) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-sampling integrity is invalid");
    }
    double maximumDifference = 0.0;
    for (std::size_t index = 0; index < samples.bindings.size(); ++index) {
        const auto& binding = samples.bindings[index];
        const auto& pressure = samples.pressures[index];
        if (binding.sampleIndex != index
            || binding.stableId == 0
            || binding.negativeSideControlVolumeIndex
                >= samples.controlVolumeCount
            || binding.positiveSideControlVolumeIndex
                >= samples.controlVolumeCount
            || binding.negativeSideControlVolumeStableId == 0
            || binding.positiveSideControlVolumeStableId == 0
            || binding.negativeSideRegionId == invalidStableId
            || binding.positiveSideRegionId == invalidStableId
            || binding.componentIndex >= samples.componentCount
            || !std::isfinite(binding.pressureDifferencePascals)
            || pressure.stableId != binding.stableId
            || !std::isfinite(pressure.negativeSidePressurePascals)
            || !std::isfinite(pressure.positiveSidePressurePascals)
            || binding.pressureDifferencePascals
                != pressure.negativeSidePressurePascals
                    - pressure.positiveSidePressurePascals) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure sample is invalid");
        }
        maximumDifference = std::max(
            maximumDifference,
            std::abs(binding.pressureDifferencePascals));
    }
    if (samples.maximumAbsolutePressureDifferencePascals
            != maximumDifference
        || productFingerprint(samples) != samples.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-sampling summary is invalid");
    }
}

void validateSceneFluidMimeticPressureSamples(
    const SceneFluidMimeticPressureSampleSet& samples,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState) {
    validateSources(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, pressureState);
    validateSceneFluidMimeticPressureSampleIntegrity(samples);
    const SceneFluidPressureSamplingLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildSamples(
        quadrature, pressureVolumes, controlCells, fullSystem,
        condensedSystem, pressureState, unlimited);
    if (samples != expected) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-sampling payload is invalid");
    }
}

ConservativeTransferResult evaluateSceneFluidMimeticPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticPressureSampleSet& samples,
    const ConservativeTransferSettings& settings) {
    validateSceneFluidMimeticPressureSampleIntegrity(samples);
    if (samples.quadratureFingerprint != quadrature.fingerprint
        || samples.surfaceDefinitionFingerprint != surface.fingerprint
        || samples.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || samples.acceptedStepCount != state.acceptedStepCount
        || samples.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-traction binding is invalid");
    }
    return evaluateSceneFluidPressureQuadrature(
        surface, state, transfer, quadrature, samples.pressures, settings);
}

} // namespace simwing::fsi
