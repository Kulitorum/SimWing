#include "scene_fluid_pressure_sampling.h"

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
            byte(static_cast<std::uint8_t>(value & 0xffU));
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
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

bool checkedMultiply(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

std::size_t storageBytesForCount(const std::size_t sampleCount) {
    std::size_t bindingBytes = 0;
    std::size_t pressureBytes = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            sampleCount, sizeof(SceneFluidPressureSampleBinding),
            bindingBytes)
        || !checkedMultiply(
            sampleCount, sizeof(SceneFluidQuadraturePressure), pressureBytes)
        || !checkedAdd(bindingBytes, pressureBytes, total)) {
        throw std::length_error(
            "scene fluid pressure-sampling storage size overflows");
    }
    return total;
}

std::size_t storageBytes(const SceneFluidPressureSampleSet& samples) {
    if (samples.bindings.size() != samples.pressures.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-sampling vector sizes differ");
    }
    return storageBytesForCount(samples.bindings.size());
}

std::uint64_t sampleFingerprint(
    const SceneFluidPressureSampleSet& samples) {
    Fingerprint fingerprint;
    fingerprint.integer(samples.version);
    fingerprint.integer(samples.quadratureFingerprint);
    fingerprint.integer(samples.pressureControlVolumeFingerprint);
    fingerprint.integer(samples.pressureProjectionFingerprint);
    fingerprint.integer(samples.surfaceDefinitionFingerprint);
    fingerprint.integer(samples.structureDefinitionFingerprint);
    fingerprint.integer(samples.acceptedStepCount);
    fingerprint.real(samples.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.ownedStorageBytes));
    fingerprint.real(samples.maximumAbsolutePressureDifferencePascals);
    fingerprint.integer(static_cast<std::uint64_t>(samples.bindings.size()));
    for (const auto& binding : samples.bindings) {
        fingerprint.integer(static_cast<std::uint64_t>(binding.sampleIndex));
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
    fingerprint.integer(static_cast<std::uint64_t>(samples.pressures.size()));
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
    const SceneFluidPressureProjection& projection) {
    validateSceneFluidQuadratureDefinition(quadrature);
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidPressureProjectionIntegrity(projection);
    if (!projection.diagnostics.accepted
        || projection.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || projection.pressurePascals.size()
            != pressureVolumes.controlVolumes.size()
        || projection.controlVolumes.size()
            != pressureVolumes.controlVolumes.size()
        || quadrature.surfaceDefinitionFingerprint
            != pressureVolumes.surfaceDefinitionFingerprint
        || quadrature.structureDefinitionFingerprint
            != pressureVolumes.structureDefinitionFingerprint
        || quadrature.acceptedStepCount != pressureVolumes.acceptedStepCount
        || quadrature.acceptedStepCount != projection.acceptedStepCount
        || quadrature.simulationTimeSeconds
            != pressureVolumes.simulationTimeSeconds
        || quadrature.simulationTimeSeconds
            != projection.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid projected-pressure sampling source identity is invalid");
    }
    for (std::size_t index = 0;
         index < pressureVolumes.controlVolumes.size(); ++index) {
        if (projection.controlVolumes[index].controlVolumeIndex != index
            || projection.controlVolumes[index].stableId
                != pressureVolumes.controlVolumes[index].stableId) {
            throw std::invalid_argument(
                "scene fluid projected-pressure control topology is invalid");
        }
    }
}

const SceneFluidPressureControlVolume& findControl(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::size_t cellIndex,
    const StableId regionId) {
    if (cellIndex >= pressureVolumes.cells.size()) {
        throw std::invalid_argument(
            "scene fluid pressure sample cell is outside the pressure grid");
    }
    const auto& cell = pressureVolumes.cells[cellIndex];
    const SceneFluidPressureControlVolume* found = nullptr;
    for (std::size_t offset = 0;
         offset < cell.controlVolumeCount; ++offset) {
        const auto& control = pressureVolumes.controlVolumes[
            cell.firstControlVolume + offset];
        if (control.regionId == regionId) {
            if (found != nullptr) {
                throw std::invalid_argument(
                    "scene fluid pressure sample has duplicate cell-region controls");
            }
            found = &control;
        }
    }
    if (found == nullptr) {
        throw std::invalid_argument(
            "scene fluid pressure sample has no matching cell-region control");
    }
    return *found;
}

SceneFluidPressureSampleSet buildSamples(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureProjection& projection,
    const SceneFluidPressureSamplingLimits& limits) {
    if (quadrature.points.size() > limits.maximumSamples) {
        throw std::length_error(
            "scene fluid pressure sampling exceeds its sample limit");
    }
    const std::size_t expectedStorageBytes = storageBytesForCount(
        quadrature.points.size());
    if (expectedStorageBytes > limits.maximumSamplingBytes) {
        throw std::length_error(
            "scene fluid pressure sampling exceeds its byte limit");
    }

    SceneFluidPressureSampleSet result;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.pressureProjectionFingerprint = projection.fingerprint;
    result.surfaceDefinitionFingerprint =
        quadrature.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        quadrature.structureDefinitionFingerprint;
    result.acceptedStepCount = quadrature.acceptedStepCount;
    result.simulationTimeSeconds = quadrature.simulationTimeSeconds;
    result.bindings.reserve(quadrature.points.size());
    result.pressures.reserve(quadrature.points.size());
    for (std::size_t index = 0;
         index < quadrature.points.size(); ++index) {
        const auto& point = quadrature.points[index];
        const auto& negative = findControl(
            pressureVolumes, point.negativeSideCellIndex,
            point.negativeSideRegionId);
        const auto& positive = findControl(
            pressureVolumes, point.positiveSideCellIndex,
            point.positiveSideRegionId);
        if (negative.componentIndex != positive.componentIndex) {
            throw std::invalid_argument(
                "scene fluid sheet sides use independently gauged pressure components");
        }
        const double negativePressure = projection.pressurePascals[
            negative.controlVolumeIndex];
        const double positivePressure = projection.pressurePascals[
            positive.controlVolumeIndex];
        const double difference = negativePressure - positivePressure;
        if (!std::isfinite(negativePressure)
            || !std::isfinite(positivePressure)
            || !std::isfinite(difference)) {
            throw std::overflow_error(
                "scene fluid projected pressure sample is not finite");
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
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes != expectedStorageBytes) {
        throw std::logic_error(
            "scene fluid pressure-sampling storage count changed");
    }
    result.fingerprint = sampleFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureSampleSet sampleSceneFluidProjectedPressure(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureProjection& projection,
    const SceneFluidPressureSamplingLimits& limits) {
    validateSources(quadrature, pressureVolumes, projection);
    auto result = buildSamples(
        quadrature, pressureVolumes, projection, limits);
    validateSceneFluidProjectedPressureSamples(
        result, quadrature, pressureVolumes, projection);
    return result;
}

void validateSceneFluidPressureSampleIntegrity(
    const SceneFluidPressureSampleSet& samples) {
    if (samples.version != sceneFluidPressureSamplingVersion
        || samples.fingerprint == 0
        || samples.quadratureFingerprint == 0
        || samples.pressureControlVolumeFingerprint == 0
        || samples.pressureProjectionFingerprint == 0
        || samples.surfaceDefinitionFingerprint == 0
        || samples.structureDefinitionFingerprint == 0
        || !std::isfinite(samples.simulationTimeSeconds)
        || !std::isfinite(samples.maximumAbsolutePressureDifferencePascals)
        || samples.bindings.empty()
        || samples.ownedStorageBytes != storageBytes(samples)
        || samples.fingerprint != sampleFingerprint(samples)) {
        throw std::invalid_argument(
            "scene fluid pressure-sampling integrity is invalid");
    }
}

void validateSceneFluidProjectedPressureSamples(
    const SceneFluidPressureSampleSet& samples,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureProjection& projection) {
    validateSources(quadrature, pressureVolumes, projection);
    validateSceneFluidPressureSampleIntegrity(samples);
    const SceneFluidPressureSamplingLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildSamples(
        quadrature, pressureVolumes, projection, unlimited);
    if (samples != expected) {
        throw std::invalid_argument(
            "scene fluid projected-pressure sample payload is invalid");
    }
}

ConservativeTransferResult evaluateSceneFluidProjectedPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureSampleSet& samples,
    const ConservativeTransferSettings& settings) {
    validateSceneFluidPressureSampleIntegrity(samples);
    if (samples.quadratureFingerprint != quadrature.fingerprint
        || samples.surfaceDefinitionFingerprint != surface.fingerprint
        || samples.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || samples.acceptedStepCount != state.acceptedStepCount
        || samples.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid projected-pressure traction binding is invalid");
    }
    return evaluateSceneFluidPressureQuadrature(
        surface, state, transfer, quadrature, samples.pressures, settings);
}

} // namespace simwing::fsi
