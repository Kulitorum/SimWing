#include "scene_fluid_mimetic_pressure_sampling.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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
    fingerprint.integer(static_cast<std::uint64_t>(
        samples.extrapolatedZeroVolumeSideCount));
    fingerprint.real(samples.maximumExtrapolationDistanceMeters);
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

const SceneFluidPressureControlVolume* findControl(
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
    return found;
}

struct MaterialSampleCentroid {
    bool present = false;
    fluid::Vector3 wrappedMeters;
    bool negativeSideOmitted = false;
    bool positiveSideOmitted = false;
};

std::vector<MaterialSampleCentroid> materialSampleCentroids(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticControlCellSet& controlCells) {
    std::vector<MaterialSampleCentroid> result(quadrature.points.size());
    for (const auto& omitted : controlCells.omittedMaterialSamples) {
        if (omitted.sourceIndex >= result.size()
            || omitted.sourceStableId
                != quadrature.points[omitted.sourceIndex].stableId) {
            throw std::invalid_argument(
                "scene fluid mimetic omitted-material sample is foreign");
        }
        result[omitted.sourceIndex] = {
            true, omitted.centroidMeters,
            omitted.negativeSideOmitted,
            omitted.positiveSideOmitted,
        };
    }
    return result;
}

double periodicSquaredDistance(
    const fluid::Vector3& point,
    const Vec3& controlCentroid,
    const fluid::Vector3& lower,
    const fluid::Vector3& upper) {
    double squared = 0.0;
    const std::array<double, 3> pointValues{
        point.x, point.y, point.z};
    const std::array<double, 3> controlValues{
        controlCentroid.x, controlCentroid.y, controlCentroid.z};
    const std::array<double, 3> lowerValues{
        lower.x, lower.y, lower.z};
    const std::array<double, 3> upperValues{
        upper.x, upper.y, upper.z};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double length = upperValues[axis] - lowerValues[axis];
        double delta = std::abs(pointValues[axis] - controlValues[axis]);
        delta = std::fmod(delta, length);
        delta = std::min(delta, length - delta);
        squared += delta * delta;
    }
    return squared;
}

const SceneFluidPressureControlVolume& resolveControl(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::size_t cellIndex,
    const StableId regionId,
    const SceneFluidQuadraturePoint& point,
    const MaterialSampleCentroid& sampleCentroid,
    const bool negativeSide,
    const char* const side,
    bool& extrapolated,
    double& extrapolationDistanceMeters) {
    if (const auto* exact = findControl(
            pressureVolumes, cellIndex, regionId)) {
        extrapolated = false;
        extrapolationDistanceMeters = 0.0;
        return *exact;
    }
    if (!sampleCentroid.present
        || (negativeSide
                ? !sampleCentroid.negativeSideOmitted
                : !sampleCentroid.positiveSideOmitted)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure sample has no matching omitted-side centroid");
    }
    const auto counts = pressureVolumes.cellCounts;
    const auto& cell = pressureVolumes.cells[cellIndex].cell;
    const std::array<std::size_t, 6> neighborCells{
        (cell.i + counts.x - 1) % counts.x
            + counts.x * (cell.j + counts.y * cell.k),
        (cell.i + 1) % counts.x
            + counts.x * (cell.j + counts.y * cell.k),
        cell.i + counts.x
            * ((cell.j + counts.y - 1) % counts.y + counts.y * cell.k),
        cell.i + counts.x
            * ((cell.j + 1) % counts.y + counts.y * cell.k),
        cell.i + counts.x
            * (cell.j + counts.y * ((cell.k + counts.z - 1) % counts.z)),
        cell.i + counts.x
            * (cell.j + counts.y * ((cell.k + 1) % counts.z)),
    };
    const SceneFluidPressureControlVolume* nearest = nullptr;
    double nearestSquaredDistance =
        std::numeric_limits<double>::infinity();
    std::array<std::size_t, 6> visited{};
    std::size_t visitedCount = 0;
    for (const std::size_t neighborCell : neighborCells) {
        if (neighborCell == cellIndex
            || std::find(
                   visited.begin(), visited.begin() + visitedCount,
                   neighborCell)
                != visited.begin() + visitedCount) {
            continue;
        }
        visited[visitedCount++] = neighborCell;
        const auto* candidate = findControl(
            pressureVolumes, neighborCell, regionId);
        if (candidate == nullptr) {
            continue;
        }
        const double squaredDistance = periodicSquaredDistance(
            sampleCentroid.wrappedMeters, candidate->centroidMeters,
            pressureVolumes.lowerMeters, pressureVolumes.upperMeters);
        if (squaredDistance < nearestSquaredDistance
            || (squaredDistance == nearestSquaredDistance
                && nearest != nullptr
                && candidate->stableId < nearest->stableId)) {
            nearest = candidate;
            nearestSquaredDistance = squaredDistance;
        }
    }
    if (nearest == nullptr || !std::isfinite(nearestSquaredDistance)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure sample has no matching "
            + std::string(side) + " control in its cell or immediate "
              "neighbors: sample=" + std::to_string(point.stableId)
            + ", triangle=" + std::to_string(point.triangleId)
            + ", cell=" + std::to_string(cellIndex)
            + ", region=" + std::to_string(regionId));
    }
    extrapolated = true;
    extrapolationDistanceMeters = std::sqrt(nearestSquaredDistance);
    return *nearest;
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
    const auto sampleCentroids = materialSampleCentroids(
        quadrature, controlCells);
    result.bindings.reserve(quadrature.points.size());
    result.pressures.reserve(quadrature.points.size());
    for (std::size_t index = 0; index < quadrature.points.size(); ++index) {
        const auto& point = quadrature.points[index];
        bool negativeExtrapolated = false;
        bool positiveExtrapolated = false;
        double negativeExtrapolationDistance = 0.0;
        double positiveExtrapolationDistance = 0.0;
        const auto& negative = resolveControl(
            pressureVolumes, point.negativeSideCellIndex,
            point.negativeSideRegionId, point, sampleCentroids[index],
            true, "negative-side", negativeExtrapolated,
            negativeExtrapolationDistance);
        const auto& positive = resolveControl(
            pressureVolumes, point.positiveSideCellIndex,
            point.positiveSideRegionId, point, sampleCentroids[index],
            false, "positive-side", positiveExtrapolated,
            positiveExtrapolationDistance);
        result.extrapolatedZeroVolumeSideCount +=
            static_cast<std::size_t>(negativeExtrapolated)
            + static_cast<std::size_t>(positiveExtrapolated);
        result.maximumExtrapolationDistanceMeters = std::max({
            result.maximumExtrapolationDistanceMeters,
            negativeExtrapolationDistance,
            positiveExtrapolationDistance,
        });
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
    if (result.extrapolatedZeroVolumeSideCount
        != controlCells.omittedZeroVolumeMaterialSideCount) {
        throw std::logic_error(
            "scene fluid mimetic pressure extrapolation count does not match omitted material sides");
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
        || samples.extrapolatedZeroVolumeSideCount
            > 2 * samples.bindings.size()
        || !std::isfinite(samples.maximumExtrapolationDistanceMeters)
        || samples.maximumExtrapolationDistanceMeters < 0.0
        || (samples.extrapolatedZeroVolumeSideCount == 0
            && samples.maximumExtrapolationDistanceMeters != 0.0)
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
