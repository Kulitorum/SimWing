#include "scene_fluid_region_connectivity.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
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

class DisjointSets final {
public:
    explicit DisjointSets(const std::size_t count)
        : parents_(count), ranks_(count) {
        std::iota(parents_.begin(), parents_.end(), std::size_t{0});
    }

    std::size_t find(const std::size_t value) {
        if (parents_[value] != value) {
            parents_[value] = find(parents_[value]);
        }
        return parents_[value];
    }

    void unite(const std::size_t first, const std::size_t second) {
        std::size_t firstRoot = find(first);
        std::size_t secondRoot = find(second);
        if (firstRoot == secondRoot) return;
        if (ranks_[firstRoot] < ranks_[secondRoot]) {
            std::swap(firstRoot, secondRoot);
        }
        parents_[secondRoot] = firstRoot;
        if (ranks_[firstRoot] == ranks_[secondRoot]) {
            ++ranks_[firstRoot];
        }
    }

private:
    std::vector<std::size_t> parents_;
    std::vector<std::uint8_t> ranks_;
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

std::size_t connectivityStorageBytes(
    const SceneFluidRegionConnectivity& connectivity) {
    std::size_t regionBytes = 0;
    std::size_t openingBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t regionIndexBytes = 0;
    std::size_t openingIndexBytes = 0;
    std::size_t first = 0;
    std::size_t second = 0;
    std::size_t third = 0;
    std::size_t total = 0;
    if (!checkedMultiply(connectivity.regions.size(),
                         sizeof(SceneFluidConnectedRegion), regionBytes)
        || !checkedMultiply(connectivity.openings.size(),
                            sizeof(SceneFluidConnectedOpening), openingBytes)
        || !checkedMultiply(connectivity.components.size(),
                            sizeof(SceneFluidRegionComponent), componentBytes)
        || !checkedMultiply(connectivity.componentRegionIndices.size(),
                            sizeof(std::size_t), regionIndexBytes)
        || !checkedMultiply(connectivity.componentOpeningIndices.size(),
                            sizeof(std::size_t), openingIndexBytes)
        || !checkedAdd(regionBytes, openingBytes, first)
        || !checkedAdd(componentBytes, regionIndexBytes, second)
        || !checkedAdd(first, second, third)
        || !checkedAdd(third, openingIndexBytes, total)) {
        throw std::length_error(
            "scene fluid region-connectivity storage size overflows");
    }
    return total;
}

std::size_t componentContinuityStorageBytes(const std::size_t count) {
    if (count != 0
        && count > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidRegionComponentContinuity)) {
        throw std::length_error(
            "scene fluid component-continuity storage size overflows");
    }
    return count * sizeof(SceneFluidRegionComponentContinuity);
}

std::uint64_t connectivityFingerprint(
    const SceneFluidRegionConnectivity& connectivity) {
    Fingerprint fingerprint;
    fingerprint.integer(connectivity.version);
    fingerprint.integer(connectivity.surfaceDefinitionFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        connectivity.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        connectivity.regions.size()));
    for (const auto& region : connectivity.regions) {
        fingerprint.integer(static_cast<std::uint64_t>(region.regionIndex));
        fingerprint.integer(region.regionId);
        fingerprint.enumeration(region.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            region.componentIndex));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        connectivity.openings.size()));
    for (const auto& opening : connectivity.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(opening.openingIndex));
        fingerprint.integer(opening.openingId);
        fingerprint.enumeration(opening.role);
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.negativeSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.positiveSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.componentIndex));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        connectivity.components.size()));
    for (const auto& component : connectivity.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.gaugeRegionId);
        fingerprint.integer(static_cast<std::uint8_t>(
            component.containsOutside));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.firstRegionMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.regionCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.firstOpeningMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.openingCount));
    }
    for (const auto& values : {
             std::pair{&connectivity.componentRegionIndices,
                       std::uint8_t{1}},
             std::pair{&connectivity.componentOpeningIndices,
                       std::uint8_t{2}}}) {
        fingerprint.integer(values.second);
        fingerprint.integer(static_cast<std::uint64_t>(
            values.first->size()));
        for (const std::size_t value : *values.first) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
    }
    return fingerprint.value();
}

SceneFluidRegionConnectivity buildConnectivity(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidRegionConnectivityLimits& limits) {
    validateSceneFluidSurfaceDefinition(surface);
    if (surface.regions.size() > limits.maximumRegions
        || surface.openings.size() > limits.maximumOpenings) {
        throw std::length_error(
            "scene fluid region-connectivity exceeds its input limits");
    }
    DisjointSets sets(surface.regions.size());
    for (const auto& opening : surface.openings) {
        sets.unite(
            opening.negativeSideRegionIndex,
            opening.positiveSideRegionIndex);
    }

    std::map<std::size_t, std::vector<std::size_t>> groupsByRoot;
    for (std::size_t region = 0; region < surface.regions.size(); ++region) {
        groupsByRoot[sets.find(region)].push_back(region);
    }
    struct Group {
        StableId gaugeRegionId = invalidStableId;
        std::vector<std::size_t> regions;
    };
    std::vector<Group> groups;
    groups.reserve(groupsByRoot.size());
    for (auto& [root, regions] : groupsByRoot) {
        static_cast<void>(root);
        std::ranges::sort(
            regions,
            [&](const std::size_t first, const std::size_t second) {
                return surface.regions[first].id < surface.regions[second].id;
            });
        groups.push_back({surface.regions[regions.front()].id,
                          std::move(regions)});
    }
    std::ranges::sort(
        groups, {}, &Group::gaugeRegionId);
    if (groups.size() > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid region-connectivity exceeds its component limit");
    }

    SceneFluidRegionConnectivity result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.regions.resize(surface.regions.size());
    result.openings.resize(surface.openings.size());
    result.components.reserve(groups.size());
    result.componentRegionIndices.reserve(surface.regions.size());
    result.componentOpeningIndices.reserve(surface.openings.size());
    std::vector<std::size_t> componentForRegion(
        surface.regions.size(), groups.size());
    for (std::size_t componentIndex = 0;
         componentIndex < groups.size(); ++componentIndex) {
        const auto& group = groups[componentIndex];
        SceneFluidRegionComponent component;
        component.componentIndex = componentIndex;
        component.gaugeRegionId = group.gaugeRegionId;
        component.firstRegionMember = result.componentRegionIndices.size();
        component.regionCount = group.regions.size();
        for (const std::size_t regionIndex : group.regions) {
            componentForRegion[regionIndex] = componentIndex;
            component.containsOutside = component.containsOutside
                || surface.regions[regionIndex].kind == RegionKind::Outside;
            result.componentRegionIndices.push_back(regionIndex);
        }
        result.components.push_back(component);
    }
    for (std::size_t regionIndex = 0;
         regionIndex < surface.regions.size(); ++regionIndex) {
        result.regions[regionIndex] = {
            regionIndex,
            surface.regions[regionIndex].id,
            surface.regions[regionIndex].kind,
            componentForRegion[regionIndex],
        };
    }
    std::vector<std::vector<std::size_t>> openingsByComponent(groups.size());
    for (std::size_t openingIndex = 0;
         openingIndex < surface.openings.size(); ++openingIndex) {
        const auto& opening = surface.openings[openingIndex];
        const std::size_t componentIndex =
            componentForRegion[opening.negativeSideRegionIndex];
        if (componentIndex
            != componentForRegion[opening.positiveSideRegionIndex]) {
            throw std::logic_error(
                "scene fluid opening escaped its connected component");
        }
        result.openings[openingIndex] = {
            openingIndex,
            opening.id,
            opening.role,
            opening.negativeSideRegionIndex,
            opening.positiveSideRegionIndex,
            componentIndex,
        };
        openingsByComponent[componentIndex].push_back(openingIndex);
    }
    for (std::size_t componentIndex = 0;
         componentIndex < result.components.size(); ++componentIndex) {
        auto& openingIndices = openingsByComponent[componentIndex];
        std::ranges::sort(
            openingIndices,
            [&](const std::size_t first, const std::size_t second) {
                return surface.openings[first].id
                    < surface.openings[second].id;
            });
        auto& component = result.components[componentIndex];
        component.firstOpeningMember =
            result.componentOpeningIndices.size();
        component.openingCount = openingIndices.size();
        result.componentOpeningIndices.insert(
            result.componentOpeningIndices.end(),
            openingIndices.begin(), openingIndices.end());
    }
    result.ownedStorageBytes = connectivityStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumConnectivityBytes) {
        throw std::length_error(
            "scene fluid region-connectivity exceeds its byte limit");
    }
    result.fingerprint = connectivityFingerprint(result);
    return result;
}

std::uint64_t componentContinuityFingerprint(
    const SceneFluidRegionComponentContinuitySet& components) {
    Fingerprint fingerprint;
    fingerprint.integer(components.version);
    for (const std::uint64_t value : {
             components.surfaceDefinitionFingerprint,
             components.regionConnectivityFingerprint,
             components.regionContinuityFingerprint,
             components.previousOpeningFluxFingerprint,
             components.currentOpeningFluxFingerprint,
             components.previousAcceptedStepCount,
             components.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(components.durationSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        components.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        components.failedComponentCount));
    fingerprint.real(
        components.maximumAbsoluteOpeningSourceResidualCubicMetersPerSecond);
    fingerprint.real(
        components.maximumAbsoluteContinuityResidualCubicMeters);
    fingerprint.real(components.globalContinuityResidualCubicMeters);
    fingerprint.integer(static_cast<std::uint8_t>(
        components.allComponentsCompatible));
    fingerprint.integer(static_cast<std::uint64_t>(components.components.size()));
    for (const auto& component : components.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.gaugeRegionId);
        fingerprint.integer(static_cast<std::uint8_t>(
            component.containsOutside));
        fingerprint.integer(static_cast<std::uint64_t>(component.regionCount));
        fingerprint.integer(static_cast<std::uint64_t>(component.openingCount));
        for (const double value : {
                 component.previousOpeningSourceResidualCubicMetersPerSecond,
                 component.currentOpeningSourceResidualCubicMetersPerSecond,
                 component.openingSourceToleranceCubicMetersPerSecond,
                 component.geometryVolumeChangeCubicMeters,
                 component.integratedOutwardRelativeVolumeCubicMeters,
                 component.continuityResidualCubicMeters,
                 component.toleranceCubicMeters}) {
            fingerprint.real(value);
        }
        fingerprint.integer(static_cast<std::uint8_t>(component.compatible));
    }
    return fingerprint.value();
}

SceneFluidRegionComponentContinuitySet buildComponentContinuity(
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionConnectivityLimits& limits) {
    validateSceneFluidRegionConnectivity(connectivity, surface);
    validateSceneFluidRegionContinuity(
        continuity, previousVolumes, currentVolumes,
        previousFlux, currentFlux);
    if (continuity.surfaceDefinitionFingerprint != surface.fingerprint
        || connectivity.surfaceDefinitionFingerprint != surface.fingerprint
        || connectivity.regions.size() != continuity.regions.size()) {
        throw std::invalid_argument(
            "scene fluid component-continuity identity is invalid");
    }
    if (connectivity.regions.size() > limits.maximumRegions
        || connectivity.openings.size() > limits.maximumOpenings
        || connectivity.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid component-continuity exceeds its count limits");
    }
    for (std::size_t regionIndex = 0;
         regionIndex < connectivity.regions.size(); ++regionIndex) {
        const auto& connected = connectivity.regions[regionIndex];
        const auto& audited = continuity.regions[regionIndex];
        if (connected.regionIndex != audited.regionIndex
            || connected.regionId != audited.regionId
            || connected.kind != audited.kind) {
            throw std::invalid_argument(
                "scene fluid component-continuity region identity is invalid");
        }
    }
    const std::size_t storageBytes = componentContinuityStorageBytes(
        connectivity.components.size());
    if (storageBytes > limits.maximumConnectivityBytes) {
        throw std::length_error(
            "scene fluid component-continuity exceeds its byte limit");
    }

    SceneFluidRegionComponentContinuitySet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.regionConnectivityFingerprint = connectivity.fingerprint;
    result.regionContinuityFingerprint = continuity.fingerprint;
    result.previousOpeningFluxFingerprint = previousFlux.fingerprint;
    result.currentOpeningFluxFingerprint = currentFlux.fingerprint;
    result.previousAcceptedStepCount = continuity.previousAcceptedStepCount;
    result.currentAcceptedStepCount = continuity.currentAcceptedStepCount;
    result.durationSeconds = continuity.durationSeconds;
    result.ownedStorageBytes = storageBytes;
    result.components.reserve(connectivity.components.size());
    for (const auto& connectedComponent : connectivity.components) {
        SceneFluidRegionComponentContinuity component;
        component.componentIndex = connectedComponent.componentIndex;
        component.gaugeRegionId = connectedComponent.gaugeRegionId;
        component.containsOutside = connectedComponent.containsOutside;
        component.regionCount = connectedComponent.regionCount;
        component.openingCount = connectedComponent.openingCount;
        for (std::size_t offset = 0;
             offset < connectedComponent.regionCount; ++offset) {
            const std::size_t regionIndex =
                connectivity.componentRegionIndices[
                    connectedComponent.firstRegionMember + offset];
            const auto& region = continuity.regions[regionIndex];
            component.previousOpeningSourceResidualCubicMetersPerSecond +=
                region.previousOutwardRelativeFlowRateCubicMetersPerSecond;
            component.currentOpeningSourceResidualCubicMetersPerSecond +=
                region.currentOutwardRelativeFlowRateCubicMetersPerSecond;
            component.geometryVolumeChangeCubicMeters +=
                region.geometryVolumeChangeCubicMeters;
            component.integratedOutwardRelativeVolumeCubicMeters +=
                region.integratedOutwardRelativeVolumeCubicMeters;
            component.continuityResidualCubicMeters +=
                region.continuityResidualCubicMeters;
            component.toleranceCubicMeters += region.toleranceCubicMeters;
        }
        component.openingSourceToleranceCubicMetersPerSecond =
            component.toleranceCubicMeters / result.durationSeconds;
        component.compatible =
            std::isfinite(component.continuityResidualCubicMeters)
            && std::isfinite(component.toleranceCubicMeters)
            && std::isfinite(
                component.openingSourceToleranceCubicMetersPerSecond)
            && std::abs(component.continuityResidualCubicMeters)
                <= component.toleranceCubicMeters
            && std::abs(component
                    .previousOpeningSourceResidualCubicMetersPerSecond)
                <= component.openingSourceToleranceCubicMetersPerSecond
            && std::abs(component
                    .currentOpeningSourceResidualCubicMetersPerSecond)
                <= component.openingSourceToleranceCubicMetersPerSecond;
        result.failedComponentCount += !component.compatible;
        result.maximumAbsoluteOpeningSourceResidualCubicMetersPerSecond =
            std::max({
                result.maximumAbsoluteOpeningSourceResidualCubicMetersPerSecond,
                std::abs(component
                    .previousOpeningSourceResidualCubicMetersPerSecond),
                std::abs(component
                    .currentOpeningSourceResidualCubicMetersPerSecond),
            });
        result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
            result.maximumAbsoluteContinuityResidualCubicMeters,
            std::abs(component.continuityResidualCubicMeters));
        result.globalContinuityResidualCubicMeters +=
            component.continuityResidualCubicMeters;
        result.components.push_back(component);
    }
    result.allComponentsCompatible = result.failedComponentCount == 0;
    if (!std::isfinite(
            result.maximumAbsoluteOpeningSourceResidualCubicMetersPerSecond)
        || !std::isfinite(
            result.maximumAbsoluteContinuityResidualCubicMeters)
        || !std::isfinite(result.globalContinuityResidualCubicMeters)) {
        throw std::invalid_argument(
            "scene fluid component-continuity ledger is non-finite");
    }
    result.fingerprint = componentContinuityFingerprint(result);
    return result;
}

} // namespace

SceneFluidRegionConnectivity buildSceneFluidRegionConnectivity(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidRegionConnectivityLimits& limits) {
    auto result = buildConnectivity(surface, limits);
    validateSceneFluidRegionConnectivity(result, surface);
    return result;
}

void validateSceneFluidRegionConnectivity(
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidSurfaceDefinition& surface) {
    validateSceneFluidSurfaceDefinition(surface);
    if (connectivity.version != sceneFluidRegionConnectivityVersion
        || connectivity.fingerprint == 0
        || connectivity.surfaceDefinitionFingerprint != surface.fingerprint) {
        throw std::invalid_argument(
            "scene fluid region-connectivity identity is invalid");
    }
    const SceneFluidRegionConnectivityLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildConnectivity(surface, unlimited);
    if (connectivity != expected
        || connectivity.ownedStorageBytes
            != connectivityStorageBytes(connectivity)
        || connectivity.fingerprint
            != connectivityFingerprint(connectivity)) {
        throw std::invalid_argument(
            "scene fluid region-connectivity payload is invalid");
    }
}

SceneFluidRegionComponentContinuitySet
auditSceneFluidRegionComponentContinuity(
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionConnectivityLimits& limits) {
    auto result = buildComponentContinuity(
        connectivity, continuity, surface, previousVolumes, currentVolumes,
        previousFlux, currentFlux, limits);
    validateSceneFluidRegionComponentContinuity(
        result, connectivity, continuity, surface,
        previousVolumes, currentVolumes, previousFlux, currentFlux);
    return result;
}

void validateSceneFluidRegionComponentContinuity(
    const SceneFluidRegionComponentContinuitySet& components,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux) {
    if (components.version != sceneFluidRegionComponentContinuityVersion
        || components.fingerprint == 0
        || components.regionConnectivityFingerprint
            != connectivity.fingerprint
        || components.regionContinuityFingerprint != continuity.fingerprint
        || components.surfaceDefinitionFingerprint != surface.fingerprint) {
        throw std::invalid_argument(
            "scene fluid component-continuity identity is invalid");
    }
    const SceneFluidRegionConnectivityLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildComponentContinuity(
        connectivity, continuity, surface, previousVolumes, currentVolumes,
        previousFlux, currentFlux, unlimited);
    if (components != expected
        || components.ownedStorageBytes
            != componentContinuityStorageBytes(components.components.size())
        || components.fingerprint
            != componentContinuityFingerprint(components)) {
        throw std::invalid_argument(
            "scene fluid component-continuity payload is invalid");
    }
}

} // namespace simwing::fsi
