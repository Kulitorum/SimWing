#include "fluid/planar_region_fragment_opening.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <compare>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        integer(static_cast<std::underlying_type_t<Enumeration>>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

struct WallKey {
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;

    auto operator<=>(const WallKey&) const = default;
};

class DisjointSet final {
public:
    explicit DisjointSet(const std::size_t count)
        : parent_(count), size_(count, 1) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    std::size_t find(const std::size_t value) {
        std::size_t root = value;
        while (parent_[root] != root) root = parent_[root];
        std::size_t current = value;
        while (parent_[current] != current) {
            const std::size_t next = parent_[current];
            parent_[current] = root;
            current = next;
        }
        return root;
    }

    void unite(const std::size_t first, const std::size_t second) {
        std::size_t firstRoot = find(first);
        std::size_t secondRoot = find(second);
        if (firstRoot == secondRoot) return;
        if (size_[firstRoot] < size_[secondRoot]
            || (size_[firstRoot] == size_[secondRoot]
                && firstRoot > secondRoot)) {
            std::swap(firstRoot, secondRoot);
        }
        parent_[secondRoot] = firstRoot;
        size_[firstRoot] += size_[secondRoot];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> size_;
};

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional fragment-opening storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional fragment-opening storage overflows");
    }
    return first + second;
}

template<typename... Counts>
std::size_t summedStorage(const Counts... counts) {
    std::size_t result = 0;
    ((result = checkedAdd(result, counts)), ...);
    return result;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    return summedStorage(
        checkedMultiply(
            openings.patches.size(),
            sizeof(PlanarPressureRegionFragmentOpeningPatch)),
        checkedMultiply(
            openings.partitions.size(),
            sizeof(PlanarPressureRegionFragmentOpeningWallPartition)),
        checkedMultiply(
            openings.openings.size(),
            sizeof(PlanarPressureRegionFragmentOpeningSummary)),
        checkedMultiply(
            openings.baseComponents.size(),
            sizeof(PlanarPressureRegionFragmentOpeningBaseComponentMapping)),
        checkedMultiply(
            openings.connectedComponents.size(),
            sizeof(PlanarPressureRegionFragmentOpeningConnectedComponent)));
}

std::size_t workingStorageBytes(const std::size_t patchCount,
                                const std::size_t pressureWallCount,
                                const std::size_t baseComponentCount) {
    return summedStorage(
        checkedMultiply(
            patchCount,
            sizeof(PlanarPressureRegionFragmentOpeningPatchDefinition)),
        checkedMultiply(patchCount, sizeof(std::uint64_t)),
        checkedMultiply(
            pressureWallCount,
            sizeof(std::pair<WallKey, std::size_t>)),
        checkedMultiply(patchCount, sizeof(std::size_t)),
        checkedMultiply(
            6 * sizeof(std::size_t), baseComponentCount));
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t openingFingerprint(
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    Fingerprint fingerprint;
    fingerprint.integer(openings.version);
    fingerprint.integer(openings.sourceFragmentFingerprint);
    fingerprint.integer(openings.sourceTopologyFingerprint);
    fingerprint.enumeration(openings.profileAxis);
    fingerprint.integer(static_cast<std::uint64_t>(openings.patches.size()));
    for (const auto& patch : openings.patches) {
        for (const std::size_t value : {
                 patch.patchIndex, patch.i, patch.j, patch.k,
                 patch.sourceFaceLinkIndex, patch.minusFragmentIndex,
                 patch.plusFragmentIndex, patch.minusBaseComponentIndex,
                 patch.plusBaseComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        for (const std::uint64_t value : {
                 patch.patchStableId, patch.openingStableId,
                 patch.surfaceStableId, patch.sourceFaceLinkStableId,
                 patch.minusFragmentStableId, patch.plusFragmentStableId,
                 patch.negativeSideRegionStableId,
                 patch.positiveSideRegionStableId}) {
            fingerprint.integer(value);
        }
        fingerprint.enumeration(patch.axis);
        fingerprint.real(patch.areaSquareMeters);
        fingerprint.real(patch.sourceWallAreaSquareMeters);
        fingerprint.real(patch.sourceWallAreaFraction);
        fingerprint.real(patch.centerDistanceMeters);
        fingerprintVector(fingerprint, patch.wrappedCentroidMeters);
        fingerprintVector(
            fingerprint, patch.unitNormalNegativeToPositive);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.partitions.size()));
    for (const auto& partition : openings.partitions) {
        for (const std::size_t value : {
                 partition.partitionIndex, partition.sourceFaceLinkIndex,
                 partition.i, partition.j, partition.k,
                 partition.openingPatchCount}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        for (const std::uint64_t value : {
                 partition.sourceFaceLinkStableId,
                 partition.surfaceStableId,
                 partition.negativeSideRegionStableId,
                 partition.positiveSideRegionStableId}) {
            fingerprint.integer(value);
        }
        fingerprint.enumeration(partition.axis);
        fingerprint.real(partition.wallAreaSquareMeters);
        fingerprint.real(partition.openingAreaSquareMeters);
        fingerprint.real(partition.solidAreaSquareMeters);
        fingerprint.real(partition.openingAreaFraction);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.openings.size()));
    for (const auto& opening : openings.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(opening.openingIndex));
        fingerprint.integer(opening.openingStableId);
        fingerprint.integer(opening.surfaceStableId);
        fingerprint.enumeration(opening.axis);
        fingerprint.integer(opening.negativeSideRegionStableId);
        fingerprint.integer(opening.positiveSideRegionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.minusBaseComponentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.plusBaseComponentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(opening.patchCount));
        fingerprint.real(opening.areaSquareMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.baseComponents.size()));
    for (const auto& mapping : openings.baseComponents) {
        fingerprint.integer(static_cast<std::uint64_t>(
            mapping.baseComponentIndex));
        fingerprint.integer(mapping.baseComponentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            mapping.connectedComponentIndex));
        fingerprint.integer(mapping.connectedComponentStableId);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.connectedComponents.size()));
    for (const auto& component : openings.connectedComponents) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.baseComponentCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.fragmentCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.openingPatchCount));
        fingerprint.real(component.volumeCubicMeters);
    }
    fingerprint.real(openings.totalOpeningAreaSquareMeters);
    fingerprint.real(openings.totalTouchedWallAreaSquareMeters);
    fingerprint.real(openings.totalSolidAreaOnTouchedWallsSquareMeters);
    fingerprint.real(openings.wallAreaPartitionResidualSquareMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        openings.workingStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningLimits& limits) {
    if (limits.maximumPatches == 0 || limits.maximumPartitions == 0
        || limits.maximumOpenings == 0
        || limits.maximumConnectedComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional fragment-opening limits are invalid");
    }
}

WallKey wallKey(
    const PlanarPressureRegionFragmentFaceLink& link) {
    return {link.surfaceStableId, link.axis, link.i, link.j, link.k};
}

WallKey wallKey(
    const PlanarPressureRegionFragmentOpeningPatchDefinition& definition) {
    return {definition.surfaceStableId, definition.axis,
            definition.i, definition.j, definition.k};
}

PlanarPressureRegionFragmentOpeningSet buildOpeningSet(
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningLimits& limits) {
    validateLimits(limits);
    const std::size_t patchCount = definitions.size();
    const std::size_t baseComponentCount = topology.components.size();
    if (patchCount > limits.maximumPatches) {
        throw std::length_error(
            "planar regional fragment-opening count limit exceeded");
    }
    std::size_t pressureWallCount = 0;
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            ++pressureWallCount;
        }
    }
    const std::size_t expectedWorkingBytes = workingStorageBytes(
        patchCount, pressureWallCount, baseComponentCount);
    const std::size_t maximumExpectedOwnedBytes = summedStorage(
        checkedMultiply(
            patchCount, sizeof(PlanarPressureRegionFragmentOpeningPatch)),
        checkedMultiply(
            patchCount,
            sizeof(PlanarPressureRegionFragmentOpeningWallPartition)),
        checkedMultiply(
            patchCount,
            sizeof(PlanarPressureRegionFragmentOpeningSummary)),
        checkedMultiply(
            baseComponentCount,
            sizeof(PlanarPressureRegionFragmentOpeningBaseComponentMapping)),
        checkedMultiply(
            baseComponentCount,
            sizeof(PlanarPressureRegionFragmentOpeningConnectedComponent)));
    if (expectedWorkingBytes > limits.maximumWorkingBytes
        || maximumExpectedOwnedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional fragment-opening byte limit exceeded");
    }

    std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition> canonical(
        definitions.begin(), definitions.end());
    std::ranges::sort(canonical, [](const auto& first, const auto& second) {
        return std::pair{first.openingStableId, first.patchStableId}
            < std::pair{second.openingStableId, second.patchStableId};
    });
    std::vector<std::uint64_t> patchIds;
    patchIds.reserve(patchCount);
    for (const auto& definition : canonical) {
        patchIds.push_back(definition.patchStableId);
    }
    std::ranges::sort(patchIds);
    if (std::ranges::adjacent_find(patchIds) != patchIds.end()) {
        throw std::invalid_argument(
            "planar regional fragment-opening patch IDs are duplicated");
    }

    std::vector<std::pair<WallKey, std::size_t>> wallLookup;
    wallLookup.reserve(pressureWallCount);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            wallLookup.emplace_back(wallKey(link), link.linkIndex);
        }
    }
    std::ranges::sort(
        wallLookup, {}, &std::pair<WallKey, std::size_t>::first);
    if (std::ranges::adjacent_find(
            wallLookup,
            [](const auto& first, const auto& second) {
                return first.first == second.first;
            }) != wallLookup.end()) {
        throw std::logic_error(
            "planar regional fragment pressure-wall keys are ambiguous");
    }

    PlanarPressureRegionFragmentOpeningSet result;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.profileAxis = topology.profileAxis;
    result.patches.reserve(patchCount);
    result.partitions.reserve(patchCount);
    result.openings.reserve(patchCount);
    result.connectedComponents.reserve(baseComponentCount);
    DisjointSet connected(baseComponentCount);
    for (const auto& definition : canonical) {
        if (definition.patchStableId == 0
            || definition.openingStableId == 0
            || definition.surfaceStableId == 0
            || definition.axis != topology.profileAxis
            || definition.negativeSideRegionStableId == 0
            || definition.positiveSideRegionStableId == 0
            || definition.negativeSideRegionStableId
                == definition.positiveSideRegionStableId
            || !std::isfinite(definition.areaSquareMeters)
            || !(definition.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "planar regional fragment-opening definition is invalid");
        }
        const auto found = std::ranges::lower_bound(
            wallLookup, wallKey(definition), {},
            &std::pair<WallKey, std::size_t>::first);
        if (found == wallLookup.end()
            || found->first != wallKey(definition)) {
            throw std::invalid_argument(
                "planar regional fragment-opening patch has no wall tile");
        }
        const auto& link = topology.links[found->second];
        if (link.minusRegionStableId
                != definition.negativeSideRegionStableId
            || link.plusRegionStableId
                != definition.positiveSideRegionStableId
            || link.minusComponentIndex == link.plusComponentIndex
            || definition.areaSquareMeters > link.areaSquareMeters) {
            throw std::invalid_argument(
                "planar regional fragment-opening patch is incompatible with its wall");
        }
        result.patches.push_back({
            result.patches.size(),
            definition.patchStableId,
            definition.openingStableId,
            definition.surfaceStableId,
            definition.axis,
            definition.i,
            definition.j,
            definition.k,
            link.linkIndex,
            link.stableId,
            link.minusFragmentIndex,
            link.plusFragmentIndex,
            link.minusFragmentStableId,
            link.plusFragmentStableId,
            link.minusComponentIndex,
            link.plusComponentIndex,
            link.minusRegionStableId,
            link.plusRegionStableId,
            definition.areaSquareMeters,
            link.areaSquareMeters,
            definition.areaSquareMeters / link.areaSquareMeters,
            link.centerDistanceMeters,
            link.wrappedCentroidMeters,
            link.unitNormalMinusToPlus,
        });
        result.totalOpeningAreaSquareMeters += definition.areaSquareMeters;
        connected.unite(
            link.minusComponentIndex, link.plusComponentIndex);
    }

    std::vector<std::size_t> partitionOrder(patchCount);
    std::iota(partitionOrder.begin(), partitionOrder.end(), 0);
    std::ranges::sort(partitionOrder, [&](const auto first, const auto second) {
        return std::pair{
                   result.patches[first].sourceFaceLinkIndex,
                   result.patches[first].patchStableId}
            < std::pair{
                   result.patches[second].sourceFaceLinkIndex,
                   result.patches[second].patchStableId};
    });
    for (std::size_t offset = 0; offset < partitionOrder.size();) {
        const auto& first = result.patches[partitionOrder[offset]];
        const auto& link = topology.links[first.sourceFaceLinkIndex];
        PlanarPressureRegionFragmentOpeningWallPartition partition;
        partition.partitionIndex = result.partitions.size();
        partition.sourceFaceLinkIndex = link.linkIndex;
        partition.sourceFaceLinkStableId = link.stableId;
        partition.surfaceStableId = link.surfaceStableId;
        partition.axis = link.axis;
        partition.i = link.i;
        partition.j = link.j;
        partition.k = link.k;
        partition.negativeSideRegionStableId = link.minusRegionStableId;
        partition.positiveSideRegionStableId = link.plusRegionStableId;
        partition.wallAreaSquareMeters = link.areaSquareMeters;
        while (offset < partitionOrder.size()
               && result.patches[partitionOrder[offset]].sourceFaceLinkIndex
                   == link.linkIndex) {
            ++partition.openingPatchCount;
            partition.openingAreaSquareMeters +=
                result.patches[partitionOrder[offset]].areaSquareMeters;
            ++offset;
        }
        if (partition.openingAreaSquareMeters
            > partition.wallAreaSquareMeters) {
            throw std::invalid_argument(
                "planar regional fragment-opening patches exceed wall area");
        }
        partition.solidAreaSquareMeters =
            partition.wallAreaSquareMeters
            - partition.openingAreaSquareMeters;
        partition.openingAreaFraction =
            partition.openingAreaSquareMeters
            / partition.wallAreaSquareMeters;
        result.totalTouchedWallAreaSquareMeters +=
            partition.wallAreaSquareMeters;
        result.totalSolidAreaOnTouchedWallsSquareMeters +=
            partition.solidAreaSquareMeters;
        result.partitions.push_back(partition);
    }
    if (result.partitions.size() > limits.maximumPartitions) {
        throw std::length_error(
            "planar regional fragment-opening partition limit exceeded");
    }

    for (std::size_t offset = 0; offset < result.patches.size();) {
        const auto& first = result.patches[offset];
        PlanarPressureRegionFragmentOpeningSummary opening;
        opening.openingIndex = result.openings.size();
        opening.openingStableId = first.openingStableId;
        opening.surfaceStableId = first.surfaceStableId;
        opening.axis = first.axis;
        opening.negativeSideRegionStableId =
            first.negativeSideRegionStableId;
        opening.positiveSideRegionStableId =
            first.positiveSideRegionStableId;
        opening.minusBaseComponentIndex = first.minusBaseComponentIndex;
        opening.plusBaseComponentIndex = first.plusBaseComponentIndex;
        while (offset < result.patches.size()
               && result.patches[offset].openingStableId
                   == opening.openingStableId) {
            const auto& patch = result.patches[offset];
            if (patch.surfaceStableId != opening.surfaceStableId
                || patch.axis != opening.axis
                || patch.negativeSideRegionStableId
                    != opening.negativeSideRegionStableId
                || patch.positiveSideRegionStableId
                    != opening.positiveSideRegionStableId
                || patch.minusBaseComponentIndex
                    != opening.minusBaseComponentIndex
                || patch.plusBaseComponentIndex
                    != opening.plusBaseComponentIndex) {
                throw std::invalid_argument(
                    "planar regional fragment-opening patches disagree within one opening");
            }
            ++opening.patchCount;
            opening.areaSquareMeters += patch.areaSquareMeters;
            ++offset;
        }
        result.openings.push_back(opening);
    }
    if (result.openings.size() > limits.maximumOpenings) {
        throw std::length_error(
            "planar regional fragment-opening opening limit exceeded");
    }

    std::vector<std::size_t> roots(baseComponentCount);
    for (std::size_t index = 0; index < baseComponentCount; ++index) {
        roots[index] = connected.find(index);
    }
    const std::size_t unassigned = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> rootToConnected(
        baseComponentCount, unassigned);
    std::vector<std::size_t> baseToConnected(baseComponentCount, unassigned);
    for (std::size_t baseIndex = 0;
         baseIndex < baseComponentCount; ++baseIndex) {
        const std::size_t root = roots[baseIndex];
        if (rootToConnected[root] == unassigned) {
            rootToConnected[root] = result.connectedComponents.size();
            PlanarPressureRegionFragmentOpeningConnectedComponent component;
            component.componentIndex = result.connectedComponents.size();
            result.connectedComponents.push_back(component);
        }
        const std::size_t connectedIndex = rootToConnected[root];
        baseToConnected[baseIndex] = connectedIndex;
        auto& component = result.connectedComponents[connectedIndex];
        ++component.baseComponentCount;
        component.fragmentCount += topology.components[baseIndex].fragmentCount;
        component.volumeCubicMeters +=
            topology.components[baseIndex].volumeCubicMeters;
    }
    result.baseComponents.resize(baseComponentCount);
    std::vector<Fingerprint> componentFingerprints(
        result.connectedComponents.size());
    for (std::size_t connectedIndex = 0;
         connectedIndex < result.connectedComponents.size();
         ++connectedIndex) {
        componentFingerprints[connectedIndex].integer(
            planarPressureRegionFragmentOpeningVersion);
        componentFingerprints[connectedIndex].integer(
            static_cast<std::uint64_t>(
                result.connectedComponents[connectedIndex]
                    .baseComponentCount));
    }
    for (std::size_t baseIndex = 0;
         baseIndex < baseComponentCount; ++baseIndex) {
        componentFingerprints[baseToConnected[baseIndex]].integer(
            topology.components[baseIndex].stableId);
    }
    for (std::size_t connectedIndex = 0;
         connectedIndex < result.connectedComponents.size();
         ++connectedIndex) {
        result.connectedComponents[connectedIndex].stableId =
            componentFingerprints[connectedIndex].value();
    }
    for (std::size_t baseIndex = 0;
         baseIndex < baseComponentCount; ++baseIndex) {
        const auto& base = topology.components[baseIndex];
        const std::size_t connectedIndex = baseToConnected[baseIndex];
        result.baseComponents[baseIndex] = {
            baseIndex,
            base.stableId,
            connectedIndex,
            result.connectedComponents[connectedIndex].stableId,
        };
    }
    for (const auto& patch : result.patches) {
        const std::size_t connectedIndex =
            result.baseComponents[patch.minusBaseComponentIndex]
                .connectedComponentIndex;
        if (connectedIndex
            != result.baseComponents[patch.plusBaseComponentIndex]
                   .connectedComponentIndex) {
            throw std::logic_error(
                "planar regional fragment-opening component union failed");
        }
        ++result.connectedComponents[connectedIndex].openingPatchCount;
    }
    if (result.connectedComponents.size()
        > limits.maximumConnectedComponents) {
        throw std::length_error(
            "planar regional fragment-opening connected-component limit exceeded");
    }

    result.wallAreaPartitionResidualSquareMeters =
        result.totalTouchedWallAreaSquareMeters
        - result.totalOpeningAreaSquareMeters
        - result.totalSolidAreaOnTouchedWallsSquareMeters;
    const double areaScale = std::max({
        result.totalTouchedWallAreaSquareMeters,
        result.totalOpeningAreaSquareMeters,
        result.totalSolidAreaOnTouchedWallsSquareMeters, 1.0});
    if (!std::isfinite(result.wallAreaPartitionResidualSquareMeters)
        || std::abs(result.wallAreaPartitionResidualSquareMeters)
            > 4096.0 * std::numeric_limits<double>::epsilon() * areaScale) {
        throw std::overflow_error(
            "planar regional fragment-opening wall partition does not close");
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    result.workingStorageBytes = expectedWorkingBytes;
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional fragment-opening owned byte limit exceeded");
    }
    result.fingerprint = openingFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningSet
buildPlanarPressureRegionFragmentOpenings(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningLimits& limits) {
    validatePlanarPressureRegionFragmentTopology(
        topology, grid, sweep, fragments, limits.topologyLimits);
    auto result = buildOpeningSet(
        fragments, topology, definitions, limits);
    validatePlanarPressureRegionFragmentOpenings(
        result, grid, sweep, fragments, topology, definitions, limits);
    return result;
}

void validatePlanarPressureRegionFragmentOpenings(
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningLimits& limits) {
    validatePlanarPressureRegionFragmentTopology(
        topology, grid, sweep, fragments, limits.topologyLimits);
    validateLimits(limits);
    if (openings.version
            != planarPressureRegionFragmentOpeningVersion
        || openings.fingerprint == 0
        || openings.sourceFragmentFingerprint != fragments.fingerprint
        || openings.sourceTopologyFingerprint != topology.fingerprint
        || openings.profileAxis != topology.profileAxis
        || !std::isfinite(openings.totalOpeningAreaSquareMeters)
        || !std::isfinite(openings.totalTouchedWallAreaSquareMeters)
        || !std::isfinite(
            openings.totalSolidAreaOnTouchedWallsSquareMeters)
        || !std::isfinite(openings.wallAreaPartitionResidualSquareMeters)
        || openings.ownedStorageBytes != ownedStorageBytes(openings)
        || openings.fingerprint != openingFingerprint(openings)
        || openings != buildOpeningSet(
                           fragments, topology, definitions, limits)) {
        throw std::invalid_argument(
            "planar regional fragment-opening set is invalid");
    }
}

} // namespace simwing::fsi::fluid
