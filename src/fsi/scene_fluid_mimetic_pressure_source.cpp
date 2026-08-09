#include "scene_fluid_mimetic_pressure_source.h"

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

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double next = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) {
            correction_ += (sum_ - next) + value;
        } else {
            correction_ += (value - next) + sum_;
        }
        sum_ = next;
    }

    [[nodiscard]] double value() const noexcept {
        return sum_ + correction_;
    }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

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

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

std::size_t storageBytes(
    const SceneFluidMimeticPressureSourceSet& sources) {
    std::size_t controlBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            sources.controls.size(),
            sizeof(SceneFluidMimeticPressureControlSource), controlBytes)
        || !checkedMultiply(
            sources.componentCount, 2 * sizeof(double), componentBytes)
        || !checkedAdd(controlBytes, componentBytes, total)) {
        throw std::length_error(
            "scene fluid mimetic pressure-source storage overflows");
    }
    return total;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticPressureSourceSet& sources) {
    Fingerprint fingerprint;
    fingerprint.integer(sources.version);
    fingerprint.integer(sources.mimeticControlCellFingerprint);
    fingerprint.integer(sources.structureDefinitionFingerprint);
    fingerprint.integer(sources.acceptedStepCount);
    fingerprint.real(sources.simulationTimeSeconds);
    fingerprint.real(sources.settings.densityKgPerCubicMeter);
    fingerprint.real(sources.settings.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(sources.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(sources.componentCount));
    fingerprint.real(
        sources.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond);
    fingerprint.real(sources
        .maximumAbsolutePredictedNetOutwardVolumeRateCubicMetersPerSecond);
    fingerprint.real(sources
        .maximumAbsolutePredictedContinuityResidualCubicMetersPerSecond);
    fingerprint.real(sources
        .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond);
    fingerprint.real(
        sources.maximumAbsoluteComponentIntegratedSourcePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(sources.controls.size()));
    for (const auto& control : sources.controls) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.real(
            control.geometryVolumeChangeRateCubicMetersPerSecond);
        fingerprint.real(
            control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            control.predictedContinuityResidualCubicMetersPerSecond);
        fingerprint.real(control.integratedSourcePascalsMeters);
    }
    for (const auto& values : {
             &sources.componentContinuityResidualsCubicMetersPerSecond,
             &sources.componentIntegratedSourcesPascalsMeters}) {
        fingerprint.integer(static_cast<std::uint64_t>(values->size()));
        for (const double value : *values) fingerprint.real(value);
    }
    return fingerprint.value();
}

void validateSettings(
    const SceneFluidMimeticPressureSourceSettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-source settings are invalid");
    }
}

SceneFluidMimeticPressureSourceSet buildSources(
    const SceneFluidMimeticControlCellSet& controlCells,
    const std::span<const double> predictedRates,
    const std::span<const double> geometryRates,
    const SceneFluidMimeticPressureSourceSettings& settings,
    const SceneFluidMimeticPressureSourceLimits& limits) {
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    validateSettings(settings);
    if (controlCells.controlCells.size() > limits.maximumControlCells) {
        throw std::length_error(
            "scene fluid mimetic pressure-source control limit exceeded");
    }
    if (predictedRates.size() != controlCells.controlCells.size()
        || (!geometryRates.empty()
            && geometryRates.size() != controlCells.controlCells.size())
        || !std::ranges::all_of(predictedRates, [](const double value) {
               return std::isfinite(value);
           })
        || !std::ranges::all_of(geometryRates, [](const double value) {
               return std::isfinite(value);
           })) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-source fields are invalid");
    }

    SceneFluidMimeticPressureSourceSet result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.settings = settings;
    for (const auto& cell : controlCells.controlCells) {
        result.componentCount = std::max(
            result.componentCount, cell.componentIndex + 1);
    }
    if (result.componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid mimetic pressure-source component limit exceeded");
    }
    result.controls.reserve(controlCells.controlCells.size());
    std::vector<CompensatedSum> componentContinuity(
        result.componentCount);
    std::vector<CompensatedSum> componentSources(result.componentCount);
    const double scale = -settings.densityKgPerCubicMeter
        / settings.timeStepSeconds;
    for (const auto& cell : controlCells.controlCells) {
        SceneFluidMimeticPressureControlSource control;
        control.controlCellIndex = cell.controlCellIndex;
        control.controlVolumeIndex = cell.controlVolumeIndex;
        control.stableId = cell.stableId;
        control.componentIndex = cell.componentIndex;
        control.geometryVolumeChangeRateCubicMetersPerSecond =
            geometryRates.empty()
            ? 0.0 : geometryRates[cell.controlCellIndex];
        control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond =
            predictedRates[cell.controlCellIndex];
        control.predictedContinuityResidualCubicMetersPerSecond =
            control.geometryVolumeChangeRateCubicMetersPerSecond
            + control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        control.integratedSourcePascalsMeters = scale
            * control.predictedContinuityResidualCubicMetersPerSecond;
        if (!std::isfinite(control.integratedSourcePascalsMeters)) {
            throw std::overflow_error(
                "scene fluid mimetic pressure source overflowed");
        }
        result.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond =
            std::max(result
                         .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
                     std::abs(control
                         .geometryVolumeChangeRateCubicMetersPerSecond));
        result
            .maximumAbsolutePredictedNetOutwardVolumeRateCubicMetersPerSecond =
            std::max(result
                         .maximumAbsolutePredictedNetOutwardVolumeRateCubicMetersPerSecond,
                     std::abs(control
                         .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond));
        result
            .maximumAbsolutePredictedContinuityResidualCubicMetersPerSecond =
            std::max(result
                         .maximumAbsolutePredictedContinuityResidualCubicMetersPerSecond,
                     std::abs(control
                         .predictedContinuityResidualCubicMetersPerSecond));
        componentContinuity[control.componentIndex].add(
            control.predictedContinuityResidualCubicMetersPerSecond);
        componentSources[control.componentIndex].add(
            control.integratedSourcePascalsMeters);
        result.controls.push_back(control);
    }
    result.componentContinuityResidualsCubicMetersPerSecond.resize(
        result.componentCount);
    result.componentIntegratedSourcesPascalsMeters.resize(
        result.componentCount);
    for (std::size_t component = 0;
         component < result.componentCount; ++component) {
        result.componentContinuityResidualsCubicMetersPerSecond[component] =
            componentContinuity[component].value();
        result.componentIntegratedSourcesPascalsMeters[component] =
            componentSources[component].value();
        result
            .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
            std::max(result
                         .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond,
                     std::abs(result
                         .componentContinuityResidualsCubicMetersPerSecond[
                             component]));
        result.maximumAbsoluteComponentIntegratedSourcePascalsMeters =
            std::max(result
                         .maximumAbsoluteComponentIntegratedSourcePascalsMeters,
                     std::abs(result
                         .componentIntegratedSourcesPascalsMeters[component]));
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic pressure-source byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    return result;
}

} // namespace

SceneFluidMimeticPressureSourceSet buildSceneFluidMimeticPressureSources(
    const SceneFluidMimeticControlCellSet& controlCells,
    const std::span<const double> predictedRates,
    const SceneFluidMimeticPressureSourceSettings& settings,
    const SceneFluidMimeticPressureSourceLimits& limits) {
    const auto result = buildSources(
        controlCells, predictedRates, {}, settings, limits);
    validateSceneFluidMimeticPressureSources(result, controlCells);
    return result;
}

SceneFluidMimeticPressureSourceSet buildSceneFluidMimeticPressureSources(
    const SceneFluidMimeticControlCellSet& controlCells,
    const std::span<const double> predictedRates,
    const std::span<const double> geometryRates,
    const SceneFluidMimeticPressureSourceSettings& settings,
    const SceneFluidMimeticPressureSourceLimits& limits) {
    const auto result = buildSources(
        controlCells, predictedRates, geometryRates, settings, limits);
    validateSceneFluidMimeticPressureSources(result, controlCells);
    return result;
}

void validateSceneFluidMimeticPressureSourceIntegrity(
    const SceneFluidMimeticPressureSourceSet& sources) {
    validateSettings(sources.settings);
    if (sources.version != sceneFluidMimeticPressureSourceVersion
        || sources.fingerprint == 0
        || sources.mimeticControlCellFingerprint == 0
        || sources.structureDefinitionFingerprint == 0
        || !std::isfinite(sources.simulationTimeSeconds)
        || sources.controls.empty()
        || sources.componentCount == 0
        || sources.componentContinuityResidualsCubicMetersPerSecond.size()
            != sources.componentCount
        || sources.componentIntegratedSourcesPascalsMeters.size()
            != sources.componentCount
        || sources.ownedStorageBytes != storageBytes(sources)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-source integrity is invalid");
    }
    std::vector<CompensatedSum> componentContinuity(
        sources.componentCount);
    std::vector<CompensatedSum> componentSources(sources.componentCount);
    double maximumGeometry = 0.0;
    double maximumPredicted = 0.0;
    double maximumContinuity = 0.0;
    const double scale = -sources.settings.densityKgPerCubicMeter
        / sources.settings.timeStepSeconds;
    for (const auto& control : sources.controls) {
        if (control.controlCellIndex
                != &control - sources.controls.data()
            || control.stableId == 0
            || control.componentIndex >= sources.componentCount
            || !std::isfinite(
                control.geometryVolumeChangeRateCubicMetersPerSecond)
            || !std::isfinite(
                control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || control.predictedContinuityResidualCubicMetersPerSecond
                != control.geometryVolumeChangeRateCubicMetersPerSecond
                    + control
                        .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond
            || control.integratedSourcePascalsMeters
                != scale
                    * control
                        .predictedContinuityResidualCubicMetersPerSecond) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-source control is invalid");
        }
        maximumGeometry = std::max(
            maximumGeometry,
            std::abs(control
                .geometryVolumeChangeRateCubicMetersPerSecond));
        maximumPredicted = std::max(
            maximumPredicted,
            std::abs(control
                .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond));
        maximumContinuity = std::max(
            maximumContinuity,
            std::abs(control
                .predictedContinuityResidualCubicMetersPerSecond));
        componentContinuity[control.componentIndex].add(
            control.predictedContinuityResidualCubicMetersPerSecond);
        componentSources[control.componentIndex].add(
            control.integratedSourcePascalsMeters);
    }
    double maximumComponentContinuity = 0.0;
    double maximumComponentSource = 0.0;
    for (std::size_t component = 0;
         component < sources.componentCount; ++component) {
        const double continuity = componentContinuity[component].value();
        const double source = componentSources[component].value();
        if (sources
                    .componentContinuityResidualsCubicMetersPerSecond[
                        component]
                != continuity
            || sources.componentIntegratedSourcesPascalsMeters[component]
                != source) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-source component is invalid");
        }
        maximumComponentContinuity = std::max(
            maximumComponentContinuity, std::abs(continuity));
        maximumComponentSource = std::max(
            maximumComponentSource, std::abs(source));
    }
    if (sources.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond
            != maximumGeometry
        || sources
                .maximumAbsolutePredictedNetOutwardVolumeRateCubicMetersPerSecond
            != maximumPredicted
        || sources
                .maximumAbsolutePredictedContinuityResidualCubicMetersPerSecond
            != maximumContinuity
        || sources
                .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond
            != maximumComponentContinuity
        || sources.maximumAbsoluteComponentIntegratedSourcePascalsMeters
            != maximumComponentSource
        || productFingerprint(sources) != sources.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-source summary is invalid");
    }
}

void validateSceneFluidMimeticPressureSources(
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticControlCellSet& controlCells) {
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    validateSceneFluidMimeticPressureSourceIntegrity(sources);
    if (sources.mimeticControlCellFingerprint != controlCells.fingerprint
        || sources.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || sources.acceptedStepCount != controlCells.acceptedStepCount
        || sources.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || sources.controls.size() != controlCells.controlCells.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-source identity is invalid");
    }
    for (std::size_t index = 0;
         index < sources.controls.size(); ++index) {
        const auto& source = sources.controls[index];
        const auto& cell = controlCells.controlCells[index];
        if (source.controlCellIndex != cell.controlCellIndex
            || source.controlVolumeIndex != cell.controlVolumeIndex
            || source.stableId != cell.stableId
            || source.componentIndex != cell.componentIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic pressure-source topology is invalid");
        }
    }
}

std::vector<double> sceneFluidMimeticIntegratedCellSources(
    const SceneFluidMimeticPressureSourceSet& sources) {
    validateSceneFluidMimeticPressureSourceIntegrity(sources);
    std::vector<double> result(sources.controls.size(), 0.0);
    for (const auto& control : sources.controls) {
        result[control.controlCellIndex] =
            control.integratedSourcePascalsMeters;
    }
    return result;
}

} // namespace simwing::fsi
