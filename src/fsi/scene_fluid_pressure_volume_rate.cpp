#include "scene_fluid_pressure_volume_rate.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
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

std::size_t storageBytesForCounts(const std::size_t controlCount,
                                  const std::size_t componentCount) {
    std::size_t controlBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            controlCount, sizeof(SceneFluidPressureControlVolumeRate),
            controlBytes)
        || !checkedMultiply(
            componentCount,
            sizeof(SceneFluidPressureComponentVolumeRate), componentBytes)
        || !checkedAdd(controlBytes, componentBytes, total)) {
        throw std::length_error(
            "scene fluid pressure-volume-rate storage size overflows");
    }
    return total;
}

std::size_t storageBytes(const SceneFluidPressureVolumeRateSet& rates) {
    return storageBytesForCounts(
        rates.controlVolumes.size(), rates.components.size());
}

std::uint64_t volumeRateFingerprint(
    const SceneFluidPressureVolumeRateSet& rates) {
    Fingerprint fingerprint;
    fingerprint.integer(rates.version);
    for (const std::uint64_t value : {
             rates.surfaceDefinitionFingerprint,
             rates.structureDefinitionFingerprint,
             rates.previousSurfaceStateFingerprint,
             rates.currentSurfaceStateFingerprint,
             rates.previousCellVolumeFingerprint,
             rates.currentCellVolumeFingerprint,
             rates.currentPressureControlVolumeFingerprint,
             rates.previousAcceptedStepCount,
             rates.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    for (const double value : {
             rates.previousSimulationTimeSeconds,
             rates.currentSimulationTimeSeconds,
             rates.durationSeconds,
             rates.lowerMeters.x,
             rates.lowerMeters.y,
             rates.lowerMeters.z,
             rates.upperMeters.x,
             rates.upperMeters.y,
             rates.upperMeters.z,
             rates.maximumAbsoluteControlVolumeChangeCubicMeters,
             rates.maximumAbsoluteControlVolumeRateCubicMetersPerSecond,
             rates.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
             rates.globalVolumeChangeCubicMeters,
             rates.globalVolumeChangeRateCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(rates.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(rates.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(rates.cellCounts.z));
    fingerprint.integer(static_cast<std::uint64_t>(rates.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        rates.previousControlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        rates.retainedControlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        rates.appearedControlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        rates.controlVolumes.size()));
    for (const auto& control : rates.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(control.cellIndex));
        fingerprint.integer(control.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.integer(static_cast<std::uint8_t>(
            control.appearedThisEpoch));
        fingerprint.real(control.previousVolumeCubicMeters);
        fingerprint.real(control.currentVolumeCubicMeters);
        fingerprint.real(control.volumeChangeCubicMeters);
        fingerprint.real(control.volumeChangeRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(rates.components.size()));
    for (const auto& component : rates.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.controlVolumeCount));
        fingerprint.real(component.previousVolumeCubicMeters);
        fingerprint.real(component.currentVolumeCubicMeters);
        fingerprint.real(component.volumeChangeCubicMeters);
        fingerprint.real(component.volumeChangeRateCubicMetersPerSecond);
    }
    return fingerprint.value();
}

bool sameGrid(const SceneFluidCellVolumeSet& first,
              const SceneFluidCellVolumeSet& second) {
    return first.cellCounts == second.cellCounts
        && first.lowerMeters == second.lowerMeters
        && first.upperMeters == second.upperMeters;
}

bool sameGrid(const SceneFluidCellVolumeSet& volumes,
              const SceneFluidPressureControlVolumeSet& pressureVolumes) {
    return volumes.cellCounts == pressureVolumes.cellCounts
        && volumes.lowerMeters == pressureVolumes.lowerMeters
        && volumes.upperMeters == pressureVolumes.upperMeters;
}

void validateSources(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes) {
    validateSceneFluidCellVolumeIntegrity(previousVolumes);
    validateSceneFluidCellVolumeIntegrity(currentVolumes);
    validateSceneFluidPressureControlVolumeIntegrity(
        currentPressureVolumes);
    if (previousVolumes.surfaceDefinitionFingerprint == 0
        || previousVolumes.structureDefinitionFingerprint == 0
        || previousVolumes.surfaceStateFingerprint == 0
        || currentVolumes.surfaceStateFingerprint == 0
        || previousVolumes.surfaceDefinitionFingerprint
            != currentVolumes.surfaceDefinitionFingerprint
        || previousVolumes.surfaceDefinitionFingerprint
            != currentPressureVolumes.surfaceDefinitionFingerprint
        || previousVolumes.structureDefinitionFingerprint
            != currentVolumes.structureDefinitionFingerprint
        || previousVolumes.structureDefinitionFingerprint
            != currentPressureVolumes.structureDefinitionFingerprint
        || currentPressureVolumes.cellVolumeFingerprint
            != currentVolumes.fingerprint
        || previousVolumes.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentVolumes.acceptedStepCount
            != previousVolumes.acceptedStepCount + 1
        || currentPressureVolumes.acceptedStepCount
            != currentVolumes.acceptedStepCount
        || !std::isfinite(previousVolumes.simulationTimeSeconds)
        || !std::isfinite(currentVolumes.simulationTimeSeconds)
        || !(currentVolumes.simulationTimeSeconds
             > previousVolumes.simulationTimeSeconds)
        || currentPressureVolumes.simulationTimeSeconds
            != currentVolumes.simulationTimeSeconds
        || previousVolumes.settings != currentVolumes.settings
        || previousVolumes.outsideRegionId != currentVolumes.outsideRegionId
        || !sameGrid(previousVolumes, currentVolumes)
        || !sameGrid(currentVolumes, currentPressureVolumes)
        || previousVolumes.cells.size() != currentVolumes.cells.size()
        || currentVolumes.cellRegionVolumes.size()
            != currentPressureVolumes.controlVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-volume-rate source identity is invalid");
    }

    for (std::size_t cellIndex = 0;
         cellIndex < currentVolumes.cells.size(); ++cellIndex) {
        const auto& previousCell = previousVolumes.cells[cellIndex];
        const auto& currentCell = currentVolumes.cells[cellIndex];
        if (previousCell.cellIndex != cellIndex
            || currentCell.cellIndex != cellIndex
            || previousCell.cell != currentCell.cell) {
            throw std::invalid_argument(
                "scene fluid pressure-volume-rate cell identity changed between epochs");
        }
        for (std::size_t offset = 0;
             offset < currentCell.regionVolumeCount; ++offset) {
            const std::size_t sourceIndex =
                currentCell.firstRegionVolume + offset;
            const auto& current =
                currentVolumes.cellRegionVolumes[sourceIndex];
            const auto& pressure =
                currentPressureVolumes.controlVolumes[sourceIndex];
            if (pressure.controlVolumeIndex != sourceIndex
                || pressure.cellIndex != cellIndex
                || pressure.regionId != current.regionId
                || pressure.volumeCubicMeters != current.volumeCubicMeters
                || !std::isfinite(current.volumeCubicMeters)
                || !(current.volumeCubicMeters > 0.0)) {
                throw std::invalid_argument(
                    "scene fluid pressure-volume-rate sparse topology is inconsistent");
            }
        }
    }

    std::map<std::pair<std::size_t, StableId>, double> currentByOwner;
    for (const auto& pressure : currentPressureVolumes.controlVolumes) {
        if (!currentByOwner.emplace(
                std::pair{pressure.cellIndex, pressure.regionId},
                pressure.volumeCubicMeters).second) {
            throw std::invalid_argument(
                "scene fluid pressure-volume-rate current owner is duplicated");
        }
    }
    for (std::size_t cellIndex = 0;
         cellIndex < previousVolumes.cells.size(); ++cellIndex) {
        const auto& cell = previousVolumes.cells[cellIndex];
        for (std::size_t offset = 0;
             offset < cell.regionVolumeCount; ++offset) {
            const auto& previous = previousVolumes.cellRegionVolumes[
                cell.firstRegionVolume + offset];
            if (!std::isfinite(previous.volumeCubicMeters)
                || !(previous.volumeCubicMeters > 0.0)
                || !currentByOwner.contains(
                    {cellIndex, previous.regionId})) {
                throw std::invalid_argument(
                    "scene fluid pressure-volume-rate control disappearance is unsupported");
            }
        }
    }
}

SceneFluidPressureVolumeRateSet buildRates(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureVolumeRateLimits& limits) {
    const std::size_t controlCount =
        currentPressureVolumes.controlVolumes.size();
    const std::size_t componentCount =
        currentPressureVolumes.components.size();
    if (controlCount > limits.maximumControlVolumes
        || componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid pressure-volume rate exceeds its count limit");
    }
    const std::size_t expectedStorageBytes = storageBytesForCounts(
        controlCount, componentCount);
    if (expectedStorageBytes > limits.maximumVolumeRateBytes) {
        throw std::length_error(
            "scene fluid pressure-volume rate exceeds its byte limit");
    }

    SceneFluidPressureVolumeRateSet result;
    result.surfaceDefinitionFingerprint =
        currentVolumes.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        currentVolumes.structureDefinitionFingerprint;
    result.previousSurfaceStateFingerprint =
        previousVolumes.surfaceStateFingerprint;
    result.currentSurfaceStateFingerprint =
        currentVolumes.surfaceStateFingerprint;
    result.previousCellVolumeFingerprint = previousVolumes.fingerprint;
    result.currentCellVolumeFingerprint = currentVolumes.fingerprint;
    result.currentPressureControlVolumeFingerprint =
        currentPressureVolumes.fingerprint;
    result.previousAcceptedStepCount = previousVolumes.acceptedStepCount;
    result.currentAcceptedStepCount = currentVolumes.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousVolumes.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentVolumes.simulationTimeSeconds;
    result.durationSeconds = result.currentSimulationTimeSeconds
        - result.previousSimulationTimeSeconds;
    result.cellCounts = currentVolumes.cellCounts;
    result.lowerMeters = currentVolumes.lowerMeters;
    result.upperMeters = currentVolumes.upperMeters;
    result.previousControlVolumeCount =
        previousVolumes.cellRegionVolumes.size();
    std::map<std::pair<std::size_t, StableId>, double> previousByOwner;
    for (std::size_t cellIndex = 0;
         cellIndex < previousVolumes.cells.size(); ++cellIndex) {
        const auto& cell = previousVolumes.cells[cellIndex];
        for (std::size_t offset = 0;
             offset < cell.regionVolumeCount; ++offset) {
            const auto& previous = previousVolumes.cellRegionVolumes[
                cell.firstRegionVolume + offset];
            previousByOwner.emplace(
                std::pair{cellIndex, previous.regionId},
                previous.volumeCubicMeters);
        }
    }
    result.controlVolumes.reserve(controlCount);
    for (const auto& pressure : currentPressureVolumes.controlVolumes) {
        const auto previous = previousByOwner.find(
            {pressure.cellIndex, pressure.regionId});
        SceneFluidPressureControlVolumeRate rate;
        rate.controlVolumeIndex = pressure.controlVolumeIndex;
        rate.stableId = pressure.stableId;
        rate.cellIndex = pressure.cellIndex;
        rate.regionId = pressure.regionId;
        rate.componentIndex = pressure.componentIndex;
        rate.appearedThisEpoch = previous == previousByOwner.end();
        rate.previousVolumeCubicMeters = rate.appearedThisEpoch
            ? 0.0 : previous->second;
        if (rate.appearedThisEpoch) {
            ++result.appearedControlVolumeCount;
        } else {
            ++result.retainedControlVolumeCount;
        }
        rate.currentVolumeCubicMeters = pressure.volumeCubicMeters;
        rate.volumeChangeCubicMeters = rate.currentVolumeCubicMeters
            - rate.previousVolumeCubicMeters;
        rate.volumeChangeRateCubicMetersPerSecond =
            rate.volumeChangeCubicMeters / result.durationSeconds;
        if (!std::isfinite(rate.volumeChangeCubicMeters)
            || !std::isfinite(
                rate.volumeChangeRateCubicMetersPerSecond)) {
            throw std::overflow_error(
                "scene fluid pressure control-volume rate is not finite");
        }
        result.maximumAbsoluteControlVolumeChangeCubicMeters = std::max(
            result.maximumAbsoluteControlVolumeChangeCubicMeters,
            std::abs(rate.volumeChangeCubicMeters));
        result.maximumAbsoluteControlVolumeRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteControlVolumeRateCubicMetersPerSecond,
                std::abs(rate.volumeChangeRateCubicMetersPerSecond));
        result.globalVolumeChangeCubicMeters +=
            rate.volumeChangeCubicMeters;
        result.controlVolumes.push_back(rate);
    }
    result.globalVolumeChangeRateCubicMetersPerSecond =
        result.globalVolumeChangeCubicMeters / result.durationSeconds;

    result.components.reserve(componentCount);
    for (const auto& source : currentPressureVolumes.components) {
        SceneFluidPressureComponentVolumeRate component;
        component.componentIndex = source.componentIndex;
        component.controlVolumeCount = source.controlVolumeCount;
        for (std::size_t offset = 0;
             offset < source.controlVolumeCount; ++offset) {
            const auto& control = result.controlVolumes[
                currentPressureVolumes.componentControlVolumeIndices[
                    source.firstControlVolumeMember + offset]];
            component.previousVolumeCubicMeters +=
                control.previousVolumeCubicMeters;
            component.currentVolumeCubicMeters +=
                control.currentVolumeCubicMeters;
        }
        component.volumeChangeCubicMeters =
            component.currentVolumeCubicMeters
            - component.previousVolumeCubicMeters;
        component.volumeChangeRateCubicMetersPerSecond =
            component.volumeChangeCubicMeters / result.durationSeconds;
        if (!std::isfinite(component.volumeChangeCubicMeters)
            || !std::isfinite(
                component.volumeChangeRateCubicMetersPerSecond)) {
            throw std::overflow_error(
                "scene fluid pressure-component volume rate is not finite");
        }
        result.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
                std::abs(
                    component.volumeChangeRateCubicMetersPerSecond));
        result.components.push_back(component);
    }
    if (!std::isfinite(result.globalVolumeChangeCubicMeters)
        || !std::isfinite(
            result.globalVolumeChangeRateCubicMetersPerSecond)) {
        throw std::overflow_error(
            "scene fluid pressure-volume global rate is not finite");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes != expectedStorageBytes) {
        throw std::logic_error(
            "scene fluid pressure-volume-rate storage count changed");
    }
    result.fingerprint = volumeRateFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureVolumeRateSet buildSceneFluidPressureVolumeRates(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureVolumeRateLimits& limits) {
    validateSources(
        previousVolumes, currentVolumes, currentPressureVolumes);
    auto result = buildRates(
        previousVolumes, currentVolumes, currentPressureVolumes, limits);
    validateSceneFluidPressureVolumeRates(
        result, previousVolumes, currentVolumes, currentPressureVolumes);
    return result;
}

void validateSceneFluidPressureVolumeRateIntegrity(
    const SceneFluidPressureVolumeRateSet& rates) {
    if (rates.version != sceneFluidPressureVolumeRateVersion
        || rates.fingerprint == 0
        || !(rates.durationSeconds > 0.0)
        || !std::isfinite(rates.durationSeconds)
        || rates.controlVolumes.empty()
        || rates.components.empty()
        || rates.previousControlVolumeCount == 0
        || rates.retainedControlVolumeCount
            != rates.previousControlVolumeCount
        || rates.retainedControlVolumeCount
                + rates.appearedControlVolumeCount
            != rates.controlVolumes.size()
        || rates.ownedStorageBytes != storageBytes(rates)
        || rates.fingerprint != volumeRateFingerprint(rates)) {
        throw std::invalid_argument(
            "scene fluid pressure-volume-rate integrity is invalid");
    }
}

void validateSceneFluidPressureVolumeRates(
    const SceneFluidPressureVolumeRateSet& rates,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes) {
    validateSources(
        previousVolumes, currentVolumes, currentPressureVolumes);
    validateSceneFluidPressureVolumeRateIntegrity(rates);
    const SceneFluidPressureVolumeRateLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildRates(
        previousVolumes, currentVolumes, currentPressureVolumes, unlimited);
    if (rates != expected) {
        throw std::invalid_argument(
            "scene fluid pressure-volume-rate payload is invalid");
    }
}

} // namespace simwing::fsi
