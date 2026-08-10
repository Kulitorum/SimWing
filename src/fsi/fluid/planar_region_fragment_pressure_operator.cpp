#include "fluid/planar_region_fragment_pressure_operator.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t rowIdentityDomain = 0x5357'5246'5052'4f57ULL;
constexpr std::uint64_t entryIdentityDomain = 0x5357'5246'5045'4e54ULL;

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

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

void validateLimits(
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits) {
    if (limits.maximumRows == 0 || limits.maximumEntries == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar regional fragment pressure-operator limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second,
                           const char* message) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(message);
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second,
                       const char* message) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(message);
    }
    return first + second;
}

std::size_t storageBytes(const std::size_t rowCount,
                         const std::size_t entryCount,
                         const std::size_t componentCount,
                         const std::size_t memberCount) {
    constexpr const char* message =
        "planar regional fragment pressure-operator storage overflows";
    const std::size_t rows = checkedProduct(
        rowCount,
        sizeof(PlanarPressureRegionFragmentPressureOperatorRow), message);
    const std::size_t entries = checkedProduct(
        entryCount,
        sizeof(PlanarPressureRegionFragmentPressureOperatorEntry), message);
    const std::size_t components = checkedProduct(
        componentCount,
        sizeof(PlanarPressureRegionFragmentPressureOperatorComponent),
        message);
    const std::size_t members = checkedProduct(
        memberCount, sizeof(std::size_t), message);
    return checkedSum(
        checkedSum(rows, entries, message),
        checkedSum(components, members, message), message);
}

std::size_t storageBytes(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator) {
    return storageBytes(
        pressureOperator.rows.size(), pressureOperator.entries.size(),
        pressureOperator.components.size(),
        pressureOperator.componentFragmentIndices.size());
}

std::uint64_t rowStableId(const std::uint64_t fragmentStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(rowIdentityDomain);
    fingerprint.integer(fragmentStableId);
    return fingerprint.value();
}

std::uint64_t entryStableId(
    const std::uint64_t sourceFaceLinkStableId,
    const std::uint64_t rowFragmentStableId,
    const std::uint64_t columnFragmentStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(entryIdentityDomain);
    fingerprint.integer(sourceFaceLinkStableId);
    fingerprint.integer(rowFragmentStableId);
    fingerprint.integer(columnFragmentStableId);
    return fingerprint.value();
}

std::uint64_t operatorFingerprint(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator) {
    Fingerprint fingerprint;
    fingerprint.integer(pressureOperator.version);
    fingerprint.integer(pressureOperator.sourceFragmentFingerprint);
    fingerprint.integer(pressureOperator.sourceTopologyFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.rows.size()));
    for (const auto& row : pressureOperator.rows) {
        fingerprint.integer(static_cast<std::uint64_t>(row.rowIndex));
        fingerprint.integer(row.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(row.fragmentIndex));
        fingerprint.integer(row.fragmentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(row.componentIndex));
        fingerprint.integer(static_cast<std::uint8_t>(row.isGauge));
        fingerprint.integer(static_cast<std::uint64_t>(row.firstEntry));
        fingerprint.integer(static_cast<std::uint64_t>(row.entryCount));
        fingerprint.real(row.diagonalGeometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.entries.size()));
    for (const auto& entry : pressureOperator.entries) {
        fingerprint.integer(static_cast<std::uint64_t>(entry.entryIndex));
        fingerprint.integer(entry.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            entry.sourceFaceLinkIndex));
        fingerprint.integer(entry.sourceFaceLinkStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            entry.columnFragmentIndex));
        fingerprint.integer(entry.columnFragmentStableId);
        fingerprint.real(entry.geometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.components.size()));
    for (const auto& component : pressureOperator.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeFragmentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.firstFragmentMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.fragmentCount));
        fingerprint.real(component.totalVolumeCubicMeters);
        fingerprint.real(component.totalGeometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.componentFragmentIndices.size()));
    for (const std::size_t index
         : pressureOperator.componentFragmentIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(index));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.includedSameRegionGridLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.excludedPressureLayerWallLinkCount));
    fingerprint.real(pressureOperator.totalGeometryWeightMeters);
    fingerprint.real(
        pressureOperator.totalDiagonalGeometryWeightMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.ownedStorageBytes));
    return fingerprint.value();
}

void validatePayload(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator) {
    const std::size_t expectedEntryCount = checkedProduct(
        pressureOperator.includedSameRegionGridLinkCount, 2,
        "planar regional pressure-operator entry count overflows");
    if (pressureOperator.version
            != planarPressureRegionFragmentPressureOperatorVersion
        || pressureOperator.fingerprint == 0
        || pressureOperator.sourceFragmentFingerprint == 0
        || pressureOperator.sourceTopologyFingerprint == 0
        || pressureOperator.rows.empty()
        || pressureOperator.components.empty()
        || pressureOperator.entries.size() != expectedEntryCount
        || !std::isfinite(pressureOperator.totalGeometryWeightMeters)
        || !(pressureOperator.totalGeometryWeightMeters > 0.0)
        || !std::isfinite(
            pressureOperator.totalDiagonalGeometryWeightMeters)
        || !(pressureOperator.totalDiagonalGeometryWeightMeters > 0.0)
        || pressureOperator.ownedStorageBytes
            != storageBytes(pressureOperator)
        || pressureOperator.fingerprint
            != operatorFingerprint(pressureOperator)) {
        throw std::invalid_argument(
            "planar regional fragment pressure-operator identity is invalid");
    }

    std::size_t expectedFirstEntry = 0;
    double diagonalTotal = 0.0;
    std::set<std::uint64_t> rowStableIds;
    std::set<std::uint64_t> entryStableIds;
    for (std::size_t rowIndex = 0;
         rowIndex < pressureOperator.rows.size(); ++rowIndex) {
        const auto& row = pressureOperator.rows[rowIndex];
        const std::size_t rowEnd = checkedSum(
            row.firstEntry, row.entryCount,
            "planar regional pressure-operator row range overflows");
        if (row.rowIndex != rowIndex || row.fragmentIndex != rowIndex
            || row.stableId == 0 || row.fragmentStableId == 0
            || !rowStableIds.insert(row.stableId).second
            || row.componentIndex >= pressureOperator.components.size()
            || row.firstEntry != expectedFirstEntry
            || rowEnd > pressureOperator.entries.size()
            || row.entryCount == 0
            || !std::isfinite(row.diagonalGeometryWeightMeters)
            || !(row.diagonalGeometryWeightMeters > 0.0)) {
            throw std::invalid_argument(
                "planar regional fragment pressure-operator row is invalid");
        }
        double diagonal = 0.0;
        std::size_t previousColumn = 0;
        std::size_t previousLink = 0;
        bool hasPrevious = false;
        for (std::size_t entryIndex = row.firstEntry;
             entryIndex < rowEnd; ++entryIndex) {
            const auto& entry = pressureOperator.entries[entryIndex];
            if (entry.entryIndex != entryIndex || entry.stableId == 0
                || !entryStableIds.insert(entry.stableId).second
                || entry.sourceFaceLinkStableId == 0
                || entry.columnFragmentIndex
                    >= pressureOperator.rows.size()
                || entry.columnFragmentIndex == rowIndex
                || entry.columnFragmentStableId
                    != pressureOperator.rows[
                        entry.columnFragmentIndex].fragmentStableId
                || !std::isfinite(entry.geometryWeightMeters)
                || !(entry.geometryWeightMeters > 0.0)
                || (hasPrevious
                    && (entry.columnFragmentIndex < previousColumn
                        || (entry.columnFragmentIndex == previousColumn
                            && entry.sourceFaceLinkIndex < previousLink)))) {
                throw std::invalid_argument(
                    "planar regional fragment pressure-operator entry is invalid");
            }
            diagonal += entry.geometryWeightMeters;
            previousColumn = entry.columnFragmentIndex;
            previousLink = entry.sourceFaceLinkIndex;
            hasPrevious = true;
        }
        if (diagonal != row.diagonalGeometryWeightMeters) {
            throw std::invalid_argument(
                "planar regional pressure-operator diagonal is invalid");
        }
        diagonalTotal += diagonal;
        expectedFirstEntry = rowEnd;
    }
    if (expectedFirstEntry != pressureOperator.entries.size()
        || diagonalTotal
            != pressureOperator.totalDiagonalGeometryWeightMeters) {
        throw std::invalid_argument(
            "planar regional pressure-operator row ledger is invalid");
    }

    std::size_t expectedFirstMember = 0;
    std::vector<bool> seenRows(pressureOperator.rows.size(), false);
    double componentWeightTotal = 0.0;
    for (std::size_t componentIndex = 0;
         componentIndex < pressureOperator.components.size();
         ++componentIndex) {
        const auto& component = pressureOperator.components[componentIndex];
        const std::size_t componentEnd = checkedSum(
            component.firstFragmentMember, component.fragmentCount,
            "planar regional pressure-operator component range overflows");
        if (component.componentIndex != componentIndex
            || component.stableId == 0 || component.regionStableId == 0
            || component.gaugeFragmentIndex
                >= pressureOperator.rows.size()
            || component.firstFragmentMember != expectedFirstMember
            || component.fragmentCount == 0
            || componentEnd
                > pressureOperator.componentFragmentIndices.size()
            || !std::isfinite(component.totalVolumeCubicMeters)
            || !(component.totalVolumeCubicMeters > 0.0)
            || !std::isfinite(component.totalGeometryWeightMeters)
            || !(component.totalGeometryWeightMeters > 0.0)
            || !pressureOperator.rows[
                    component.gaugeFragmentIndex].isGauge
            || pressureOperator.rows[
                    component.gaugeFragmentIndex].componentIndex
                != componentIndex) {
            throw std::invalid_argument(
                "planar regional fragment pressure-operator component is invalid");
        }
        std::uint64_t previousStableId = 0;
        for (std::size_t member = component.firstFragmentMember;
             member < componentEnd; ++member) {
            const std::size_t rowIndex =
                pressureOperator.componentFragmentIndices[member];
            if (rowIndex >= pressureOperator.rows.size()
                || seenRows[rowIndex]
                || pressureOperator.rows[rowIndex].componentIndex
                    != componentIndex
                || (member != component.firstFragmentMember
                    && pressureOperator.rows[rowIndex].fragmentStableId
                        <= previousStableId)) {
                throw std::invalid_argument(
                    "planar regional pressure-operator component member is invalid");
            }
            seenRows[rowIndex] = true;
            previousStableId =
                pressureOperator.rows[rowIndex].fragmentStableId;
        }
        componentWeightTotal += component.totalGeometryWeightMeters;
        expectedFirstMember = componentEnd;
    }
    const double componentWeightScale = std::max({
        std::abs(componentWeightTotal),
        std::abs(pressureOperator.totalGeometryWeightMeters), 1.0});
    const double componentWeightTolerance = 2048.0
        * std::numeric_limits<double>::epsilon() * componentWeightScale;
    if (expectedFirstMember
            != pressureOperator.componentFragmentIndices.size()
        || expectedFirstMember != pressureOperator.rows.size()
        || std::ranges::find(seenRows, false) != seenRows.end()
        || std::abs(componentWeightTotal
                    - pressureOperator.totalGeometryWeightMeters)
            > componentWeightTolerance) {
        throw std::invalid_argument(
            "planar regional pressure-operator component ledger is invalid");
    }
}

PlanarPressureRegionFragmentPressureOperator buildOperator(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentTopology(
        topology, grid, sweep, fragments, limits.topologyLimits);
    const std::size_t rowCount = fragments.fragments.size();
    const std::size_t entryCount = checkedProduct(
        topology.sameRegionGridLinkCount, 2,
        "planar regional fragment pressure-operator entry count overflows");
    const std::size_t componentCount = topology.components.size();
    if (rowCount > limits.maximumRows) {
        throw std::length_error(
            "planar regional fragment pressure-operator row limit exceeded");
    }
    if (entryCount > limits.maximumEntries) {
        throw std::length_error(
            "planar regional fragment pressure-operator entry limit exceeded");
    }
    if (componentCount > limits.maximumComponents) {
        throw std::length_error(
            "planar regional fragment pressure-operator component limit exceeded");
    }
    const std::size_t ownedBytes = storageBytes(
        rowCount, entryCount, componentCount, rowCount);
    if (ownedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional fragment pressure-operator byte limit exceeded");
    }

    PlanarPressureRegionFragmentPressureOperator result;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.includedSameRegionGridLinkCount =
        topology.sameRegionGridLinkCount;
    result.excludedPressureLayerWallLinkCount =
        topology.pressureLayerWallLinkCount;
    result.ownedStorageBytes = ownedBytes;
    result.rows.resize(rowCount);
    result.entries.resize(entryCount);
    result.components.reserve(componentCount);
    result.componentFragmentIndices.reserve(rowCount);

    std::vector<std::size_t> degrees(rowCount, 0);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            if (link.sameRegionGeometryWeightMeters != 0.0) {
                throw std::invalid_argument(
                    "planar regional pressure wall has conductance");
            }
            continue;
        }
        if (link.minusFragmentIndex >= rowCount
            || link.plusFragmentIndex >= rowCount
            || link.minusFragmentIndex == link.plusFragmentIndex
            || link.minusComponentIndex != link.plusComponentIndex
            || !std::isfinite(link.sameRegionGeometryWeightMeters)
            || !(link.sameRegionGeometryWeightMeters > 0.0)) {
            throw std::invalid_argument(
                "planar regional grid link cannot enter the pressure operator");
        }
        ++degrees[link.minusFragmentIndex];
        ++degrees[link.plusFragmentIndex];
        result.totalGeometryWeightMeters +=
            link.sameRegionGeometryWeightMeters;
    }

    std::size_t nextEntry = 0;
    for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const auto& fragment = fragments.fragments[rowIndex];
        const auto& topologyFragment = topology.fragments[rowIndex];
        auto& row = result.rows[rowIndex];
        row.rowIndex = rowIndex;
        row.stableId = rowStableId(fragment.stableId);
        row.fragmentIndex = rowIndex;
        row.fragmentStableId = fragment.stableId;
        row.componentIndex = topologyFragment.componentIndex;
        row.firstEntry = nextEntry;
        row.entryCount = degrees[rowIndex];
        nextEntry = checkedSum(
            nextEntry, row.entryCount,
            "planar regional pressure-operator row range overflows");
    }
    if (nextEntry != entryCount) {
        throw std::logic_error(
            "planar regional pressure-operator entry count is inconsistent");
    }

    std::vector<std::size_t> cursors(rowCount, 0);
    for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        cursors[rowIndex] = result.rows[rowIndex].firstEntry;
    }
    std::set<std::uint64_t> entryStableIds;
    const auto appendEntry = [&](
        const PlanarPressureRegionFragmentFaceLink& link,
        const std::size_t sourceIndex,
        const std::size_t destinationIndex) {
        auto& entry = result.entries[cursors[sourceIndex]++];
        entry.stableId = entryStableId(
            link.stableId, fragments.fragments[sourceIndex].stableId,
            fragments.fragments[destinationIndex].stableId);
        entry.sourceFaceLinkIndex = link.linkIndex;
        entry.sourceFaceLinkStableId = link.stableId;
        entry.columnFragmentIndex = destinationIndex;
        entry.columnFragmentStableId =
            fragments.fragments[destinationIndex].stableId;
        entry.geometryWeightMeters =
            link.sameRegionGeometryWeightMeters;
        if (!entryStableIds.insert(entry.stableId).second) {
            throw std::invalid_argument(
                "planar regional pressure-operator entry stable ID collides");
        }
    };
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        appendEntry(
            link, link.minusFragmentIndex, link.plusFragmentIndex);
        appendEntry(
            link, link.plusFragmentIndex, link.minusFragmentIndex);
    }
    for (auto& row : result.rows) {
        auto first = result.entries.begin()
            + static_cast<std::ptrdiff_t>(row.firstEntry);
        auto last = first + static_cast<std::ptrdiff_t>(row.entryCount);
        std::sort(
            first, last,
            [](const auto& firstEntry, const auto& secondEntry) {
                if (firstEntry.columnFragmentIndex
                    != secondEntry.columnFragmentIndex) {
                    return firstEntry.columnFragmentIndex
                        < secondEntry.columnFragmentIndex;
                }
                return firstEntry.sourceFaceLinkIndex
                    < secondEntry.sourceFaceLinkIndex;
            });
        for (auto entry = first; entry != last; ++entry) {
            entry->entryIndex = static_cast<std::size_t>(
                entry - result.entries.begin());
            row.diagonalGeometryWeightMeters +=
                entry->geometryWeightMeters;
        }
        result.totalDiagonalGeometryWeightMeters +=
            row.diagonalGeometryWeightMeters;
    }

    for (const auto& sourceComponent : topology.components) {
        std::vector<std::size_t> members;
        for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
            if (topology.fragments[rowIndex].componentIndex
                == sourceComponent.componentIndex) {
                members.push_back(rowIndex);
            }
        }
        std::ranges::sort(
            members,
            [&](const std::size_t first, const std::size_t second) {
                return fragments.fragments[first].stableId
                    < fragments.fragments[second].stableId;
            });
        if (members.empty()) {
            throw std::invalid_argument(
                "planar regional pressure-operator component is empty");
        }
        PlanarPressureRegionFragmentPressureOperatorComponent component;
        component.componentIndex = sourceComponent.componentIndex;
        component.stableId = sourceComponent.stableId;
        component.regionStableId = sourceComponent.regionStableId;
        component.gaugeFragmentIndex = members.front();
        component.firstFragmentMember =
            result.componentFragmentIndices.size();
        component.fragmentCount = members.size();
        component.totalVolumeCubicMeters =
            sourceComponent.volumeCubicMeters;
        for (const std::size_t member : members) {
            result.componentFragmentIndices.push_back(member);
            component.totalGeometryWeightMeters +=
                0.5 * result.rows[member]
                    .diagonalGeometryWeightMeters;
        }
        result.rows[component.gaugeFragmentIndex].isGauge = true;
        result.components.push_back(component);
    }
    result.fingerprint = operatorFingerprint(result);
    validatePayload(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentPressureOperator
buildPlanarPressureRegionFragmentPressureOperator(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits) {
    auto result = buildOperator(
        grid, sweep, fragments, topology, limits);
    validatePlanarPressureRegionFragmentPressureOperator(
        result, grid, sweep, fragments, topology, limits);
    return result;
}

void validatePlanarPressureRegionFragmentPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits) {
    validateLimits(limits);
    if (pressureOperator.rows.size() > limits.maximumRows
        || pressureOperator.entries.size() > limits.maximumEntries
        || pressureOperator.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional pressure-operator validation limit exceeded");
    }
    validatePayload(pressureOperator);
    const auto expected = buildOperator(
        grid, sweep, fragments, topology, limits);
    if (expected != pressureOperator) {
        throw std::invalid_argument(
            "planar regional fragment pressure operator is invalid");
    }
}

std::vector<double> applyPlanarPressureRegionFragmentPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const std::span<const double> pressurePascals) {
    validatePayload(pressureOperator);
    if (pressurePascals.size() != pressureOperator.rows.size()) {
        throw std::invalid_argument(
            "planar regional pressure vector has the wrong size");
    }
    for (const double value : pressurePascals) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "planar regional pressure vector is not finite");
        }
    }
    std::vector<double> result(pressureOperator.rows.size(), 0.0);
    for (const auto& row : pressureOperator.rows) {
        for (std::size_t offset = 0; offset < row.entryCount; ++offset) {
            const auto& entry = pressureOperator.entries[
                row.firstEntry + offset];
            result[row.rowIndex] += entry.geometryWeightMeters
                * (pressurePascals[row.rowIndex]
                   - pressurePascals[entry.columnFragmentIndex]);
        }
        if (!std::isfinite(result[row.rowIndex])) {
            throw std::overflow_error(
                "planar regional pressure-operator action is not finite");
        }
    }
    return result;
}

} // namespace simwing::fsi::fluid
