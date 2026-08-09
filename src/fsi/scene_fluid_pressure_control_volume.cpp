#include "scene_fluid_pressure_control_volume.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
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

double tolerance(const SceneFluidCellVolumeSet& volumes,
                 const double reference) {
    return std::max(
        volumes.settings.absoluteVolumeToleranceCubicMeters,
        volumes.settings.relativeVolumeTolerance * std::abs(reference));
}

bool finite(const Vec3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

std::uint64_t controlVolumeStableId(const std::size_t cellIndex,
                                    const StableId regionId) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x7072657373766f6cULL});
    fingerprint.integer(static_cast<std::uint64_t>(cellIndex));
    fingerprint.integer(regionId);
    return fingerprint.value();
}

std::size_t storageBytesForCounts(const std::size_t cellCount,
                                  const std::size_t controlVolumeCount,
                                  const std::size_t regionCount,
                                  const std::size_t componentCount) {
    std::size_t cellBytes = 0;
    std::size_t controlVolumeBytes = 0;
    std::size_t regionBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t regionMemberBytes = 0;
    std::size_t componentMemberBytes = 0;
    std::size_t first = 0;
    std::size_t second = 0;
    std::size_t third = 0;
    std::size_t fourth = 0;
    std::size_t total = 0;
    if (!checkedMultiply(cellCount,
                         sizeof(SceneFluidPressureCell), cellBytes)
        || !checkedMultiply(
            controlVolumeCount,
            sizeof(SceneFluidPressureControlVolume), controlVolumeBytes)
        || !checkedMultiply(regionCount,
                            sizeof(SceneFluidPressureRegion), regionBytes)
        || !checkedMultiply(componentCount,
                            sizeof(SceneFluidPressureComponent), componentBytes)
        || !checkedMultiply(controlVolumeCount,
                            sizeof(std::size_t), regionMemberBytes)
        || !checkedMultiply(
            controlVolumeCount,
            sizeof(std::size_t), componentMemberBytes)
        || !checkedAdd(cellBytes, controlVolumeBytes, first)
        || !checkedAdd(regionBytes, componentBytes, second)
        || !checkedAdd(first, second, third)
        || !checkedAdd(regionMemberBytes, componentMemberBytes, fourth)
        || !checkedAdd(third, fourth, total)) {
        throw std::length_error(
            "scene fluid pressure-control-volume storage size overflows");
    }
    return total;
}

std::size_t storageBytes(
    const SceneFluidPressureControlVolumeSet& pressureVolumes) {
    if (pressureVolumes.regionControlVolumeIndices.size()
            != pressureVolumes.controlVolumes.size()
        || pressureVolumes.componentControlVolumeIndices.size()
            != pressureVolumes.controlVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume member counts are invalid");
    }
    return storageBytesForCounts(
        pressureVolumes.cells.size(), pressureVolumes.controlVolumes.size(),
        pressureVolumes.regions.size(), pressureVolumes.components.size());
}

std::uint64_t pressureVolumeFingerprint(
    const SceneFluidPressureControlVolumeSet& pressureVolumes) {
    Fingerprint fingerprint;
    fingerprint.integer(pressureVolumes.version);
    for (const std::uint64_t value : {
             pressureVolumes.surfaceDefinitionFingerprint,
             pressureVolumes.surfaceStateFingerprint,
             pressureVolumes.cellVolumeFingerprint,
             pressureVolumes.regionConnectivityFingerprint,
             pressureVolumes.structureDefinitionFingerprint,
             pressureVolumes.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(pressureVolumes.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.cellCounts.z));
    for (const double value : {
             pressureVolumes.lowerMeters.x,
             pressureVolumes.lowerMeters.y,
             pressureVolumes.lowerMeters.z,
             pressureVolumes.upperMeters.x,
             pressureVolumes.upperMeters.y,
             pressureVolumes.upperMeters.z,
             pressureVolumes.cellVolumeCubicMeters,
             pressureVolumes.domainVolumeCubicMeters,
             pressureVolumes.maximumCellVolumeResidualCubicMeters,
             pressureVolumes.maximumRegionVolumeResidualCubicMeters,
             pressureVolumes.domainVolumeResidualCubicMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.cells.size()));
    for (const auto& cell : pressureVolumes.cells) {
        fingerprint.integer(static_cast<std::uint64_t>(cell.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.firstControlVolume));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.controlVolumeCount));
        fingerprint.real(cell.assignedVolumeCubicMeters);
        fingerprint.real(cell.volumeResidualCubicMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.controlVolumes.size()));
    for (const auto& control : pressureVolumes.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(control.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(control.regionIndex));
        fingerprint.integer(control.regionId);
        fingerprint.enumeration(control.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.real(control.volumeCubicMeters);
        fingerprint.real(control.volumeFraction);
        fingerprint.real(control.centroidMeters.x);
        fingerprint.real(control.centroidMeters.y);
        fingerprint.real(control.centroidMeters.z);
        fingerprint.integer(static_cast<std::uint8_t>(
            control.belongsToGaugeRegion));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.regions.size()));
    for (const auto& region : pressureVolumes.regions) {
        fingerprint.integer(static_cast<std::uint64_t>(region.regionIndex));
        fingerprint.integer(region.regionId);
        fingerprint.enumeration(region.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            region.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            region.firstControlVolumeMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            region.controlVolumeCount));
        fingerprint.real(region.summedControlVolumeCubicMeters);
        fingerprint.real(region.sourceVolumeResidualCubicMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureVolumes.components.size()));
    for (const auto& component : pressureVolumes.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.gaugeRegionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.firstControlVolumeMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.controlVolumeCount));
        fingerprint.real(component.totalVolumeCubicMeters);
    }
    for (const auto& members : {
             &pressureVolumes.regionControlVolumeIndices,
             &pressureVolumes.componentControlVolumeIndices}) {
        fingerprint.integer(static_cast<std::uint64_t>(members->size()));
        for (const std::size_t member : *members) {
            fingerprint.integer(static_cast<std::uint64_t>(member));
        }
    }
    return fingerprint.value();
}

SceneFluidPressureControlVolumeSet buildPressureVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeLimits& limits) {
    if (volumes.cells.size() > limits.maximumCells
        || volumes.cellRegionVolumes.size() > limits.maximumControlVolumes
        || surface.regions.size() > limits.maximumRegions
        || surface.openings.size() > limits.maximumOpenings
        || connectivity.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid pressure-control-volume exceeds its count limits");
    }
    const std::size_t expectedStorageBytes = storageBytesForCounts(
        volumes.cells.size(), volumes.cellRegionVolumes.size(),
        connectivity.regions.size(), connectivity.components.size());
    if (expectedStorageBytes > limits.maximumControlVolumeBytes) {
        throw std::length_error(
            "scene fluid pressure-control-volume exceeds its byte limit");
    }
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidCellVolumeIntegrity(volumes);
    validateSceneFluidRegionConnectivity(connectivity, surface);
    if (volumes.surfaceDefinitionFingerprint != surface.fingerprint
        || connectivity.surfaceDefinitionFingerprint != surface.fingerprint
        || volumes.regionVolumes.size() != connectivity.regions.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume source identity is invalid");
    }
    SceneFluidPressureControlVolumeSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = volumes.surfaceStateFingerprint;
    result.cellVolumeFingerprint = volumes.fingerprint;
    result.regionConnectivityFingerprint = connectivity.fingerprint;
    result.structureDefinitionFingerprint =
        volumes.structureDefinitionFingerprint;
    result.acceptedStepCount = volumes.acceptedStepCount;
    result.simulationTimeSeconds = volumes.simulationTimeSeconds;
    result.cellCounts = volumes.cellCounts;
    result.lowerMeters = volumes.lowerMeters;
    result.upperMeters = volumes.upperMeters;
    result.cellVolumeCubicMeters = volumes.cellVolumeCubicMeters;
    result.domainVolumeCubicMeters =
        static_cast<double>(volumes.cells.size())
        * volumes.cellVolumeCubicMeters;
    result.cells.reserve(volumes.cells.size());
    result.controlVolumes.reserve(volumes.cellRegionVolumes.size());
    result.regions.resize(connectivity.regions.size());
    result.components.resize(connectivity.components.size());
    std::vector<std::vector<std::size_t>> controlsByRegion(
        connectivity.regions.size());
    std::vector<std::vector<std::size_t>> controlsByComponent(
        connectivity.components.size());
    std::unordered_set<std::uint64_t> stableIds;
    stableIds.reserve(volumes.cellRegionVolumes.size());
    std::size_t expectedFirstSource = 0;
    double summedDomainVolume = 0.0;
    for (const auto& sourceCell : volumes.cells) {
        if (sourceCell.cellIndex != result.cells.size()
            || sourceCell.firstRegionVolume != expectedFirstSource
            || sourceCell.firstRegionVolume
                    + sourceCell.regionVolumeCount
                > volumes.cellRegionVolumes.size()) {
            throw std::invalid_argument(
                "scene fluid pressure-control-volume cell source is invalid");
        }
        SceneFluidPressureCell cell;
        cell.cellIndex = sourceCell.cellIndex;
        cell.cell = sourceCell.cell;
        cell.firstControlVolume = result.controlVolumes.size();
        for (std::size_t offset = 0;
             offset < sourceCell.regionVolumeCount; ++offset) {
            const auto& source = volumes.cellRegionVolumes[
                sourceCell.firstRegionVolume + offset];
            const auto regionIndex = surface.mappings.regionIndex(
                source.regionId);
            if (!regionIndex || *regionIndex >= connectivity.regions.size()
                || !(source.volumeCubicMeters > 0.0)
                || !std::isfinite(source.volumeCubicMeters)
                || !std::isfinite(source.volumeFraction)
                || !finite(source.centroidMeters)) {
                throw std::invalid_argument(
                    "scene fluid pressure-control-volume region source is invalid");
            }
            const auto& connectedRegion = connectivity.regions[*regionIndex];
            const auto& component =
                connectivity.components[connectedRegion.componentIndex];
            SceneFluidPressureControlVolume control;
            control.controlVolumeIndex = result.controlVolumes.size();
            control.stableId = controlVolumeStableId(
                sourceCell.cellIndex, source.regionId);
            control.cellIndex = sourceCell.cellIndex;
            control.regionIndex = *regionIndex;
            control.regionId = source.regionId;
            control.kind = connectedRegion.kind;
            control.componentIndex = connectedRegion.componentIndex;
            control.volumeCubicMeters = source.volumeCubicMeters;
            control.volumeFraction = source.volumeFraction;
            control.centroidMeters = source.centroidMeters;
            control.belongsToGaugeRegion =
                source.regionId == component.gaugeRegionId;
            if (!stableIds.insert(control.stableId).second) {
                throw std::invalid_argument(
                    "scene fluid pressure-control-volume stable ID collides");
            }
            controlsByRegion[control.regionIndex].push_back(
                control.controlVolumeIndex);
            controlsByComponent[control.componentIndex].push_back(
                control.controlVolumeIndex);
            cell.assignedVolumeCubicMeters += control.volumeCubicMeters;
            result.controlVolumes.push_back(control);
        }
        expectedFirstSource += sourceCell.regionVolumeCount;
        cell.controlVolumeCount = result.controlVolumes.size()
            - cell.firstControlVolume;
        cell.volumeResidualCubicMeters = cell.assignedVolumeCubicMeters
            - result.cellVolumeCubicMeters;
        result.maximumCellVolumeResidualCubicMeters = std::max(
            result.maximumCellVolumeResidualCubicMeters,
            std::abs(cell.volumeResidualCubicMeters));
        if (cell.controlVolumeCount == 0
            || std::abs(cell.volumeResidualCubicMeters)
                > tolerance(volumes, result.cellVolumeCubicMeters)) {
            throw std::invalid_argument(
                "scene fluid pressure control volumes do not close a cell");
        }
        summedDomainVolume += cell.assignedVolumeCubicMeters;
        result.cells.push_back(cell);
    }
    if (expectedFirstSource != volumes.cellRegionVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume source range is incomplete");
    }

    result.regionControlVolumeIndices.reserve(result.controlVolumes.size());
    for (std::size_t regionIndex = 0;
         regionIndex < connectivity.regions.size(); ++regionIndex) {
        const auto& connected = connectivity.regions[regionIndex];
        const auto& source = volumes.regionVolumes[regionIndex];
        if (source.regionId != connected.regionId) {
            throw std::invalid_argument(
                "scene fluid pressure-control-volume region order is invalid");
        }
        SceneFluidPressureRegion region;
        region.regionIndex = regionIndex;
        region.regionId = connected.regionId;
        region.kind = connected.kind;
        region.componentIndex = connected.componentIndex;
        region.firstControlVolumeMember =
            result.regionControlVolumeIndices.size();
        region.controlVolumeCount = controlsByRegion[regionIndex].size();
        for (const std::size_t controlIndex : controlsByRegion[regionIndex]) {
            result.regionControlVolumeIndices.push_back(controlIndex);
            region.summedControlVolumeCubicMeters +=
                result.controlVolumes[controlIndex].volumeCubicMeters;
        }
        region.sourceVolumeResidualCubicMeters =
            region.summedControlVolumeCubicMeters
            - source.summedCellVolumeCubicMeters;
        result.maximumRegionVolumeResidualCubicMeters = std::max(
            result.maximumRegionVolumeResidualCubicMeters,
            std::abs(region.sourceVolumeResidualCubicMeters));
        if (region.controlVolumeCount == 0
            || std::abs(region.sourceVolumeResidualCubicMeters)
                > tolerance(volumes, source.summedCellVolumeCubicMeters)) {
            throw std::invalid_argument(
                "scene fluid pressure control volumes do not close a region");
        }
        result.regions[regionIndex] = region;
    }

    result.componentControlVolumeIndices.reserve(
        result.controlVolumes.size());
    for (std::size_t componentIndex = 0;
         componentIndex < connectivity.components.size(); ++componentIndex) {
        const auto& connected = connectivity.components[componentIndex];
        SceneFluidPressureComponent component;
        component.componentIndex = componentIndex;
        component.gaugeRegionId = connected.gaugeRegionId;
        component.firstControlVolumeMember =
            result.componentControlVolumeIndices.size();
        component.controlVolumeCount =
            controlsByComponent[componentIndex].size();
        component.gaugeControlVolumeIndex = result.controlVolumes.size();
        for (const std::size_t controlIndex :
             controlsByComponent[componentIndex]) {
            result.componentControlVolumeIndices.push_back(controlIndex);
            const auto& control = result.controlVolumes[controlIndex];
            component.totalVolumeCubicMeters += control.volumeCubicMeters;
            if (control.belongsToGaugeRegion
                && component.gaugeControlVolumeIndex
                    == result.controlVolumes.size()) {
                component.gaugeControlVolumeIndex = controlIndex;
            }
        }
        if (component.controlVolumeCount == 0
            || component.gaugeControlVolumeIndex
                == result.controlVolumes.size()
            || !(component.totalVolumeCubicMeters > 0.0)
            || !std::isfinite(component.totalVolumeCubicMeters)) {
            throw std::invalid_argument(
                "scene fluid pressure component has no valid gauge volume");
        }
        result.components[componentIndex] = component;
    }
    result.domainVolumeResidualCubicMeters = summedDomainVolume
        - result.domainVolumeCubicMeters;
    if (!std::isfinite(result.domainVolumeResidualCubicMeters)
        || std::abs(result.domainVolumeResidualCubicMeters)
            > tolerance(volumes, result.domainVolumeCubicMeters)) {
        throw std::invalid_argument(
            "scene fluid pressure control volumes do not close the domain");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes != expectedStorageBytes) {
        throw std::logic_error(
            "scene fluid pressure-control-volume storage count changed");
    }
    result.fingerprint = pressureVolumeFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureControlVolumeSet buildSceneFluidPressureControlVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeLimits& limits) {
    auto result = buildPressureVolumes(
        surface, volumes, connectivity, limits);
    validateSceneFluidPressureControlVolumes(
        result, surface, volumes, connectivity);
    return result;
}

void validateSceneFluidPressureControlVolumeIntegrity(
    const SceneFluidPressureControlVolumeSet& pressureVolumes) {
    if (pressureVolumes.version
            != sceneFluidPressureControlVolumeVersion
        || pressureVolumes.fingerprint == 0
        || pressureVolumes.ownedStorageBytes
            != storageBytes(pressureVolumes)
        || pressureVolumes.fingerprint
            != pressureVolumeFingerprint(pressureVolumes)) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume integrity is invalid");
    }
}

void validateSceneFluidPressureControlVolumes(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity) {
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    if (pressureVolumes.version != sceneFluidPressureControlVolumeVersion
        || pressureVolumes.fingerprint == 0
        || pressureVolumes.surfaceDefinitionFingerprint != surface.fingerprint
        || pressureVolumes.cellVolumeFingerprint != volumes.fingerprint
        || pressureVolumes.regionConnectivityFingerprint
            != connectivity.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume identity is invalid");
    }
    const SceneFluidPressureControlVolumeLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildPressureVolumes(
        surface, volumes, connectivity, unlimited);
    if (pressureVolumes != expected) {
        throw std::invalid_argument(
            "scene fluid pressure-control-volume payload is invalid");
    }
}

} // namespace simwing::fsi
