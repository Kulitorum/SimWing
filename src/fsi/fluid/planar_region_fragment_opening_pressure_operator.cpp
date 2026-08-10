#include "fluid/planar_region_fragment_opening_pressure_operator.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t rowDomain = 0x5357'4f50'5052'4f57ULL;
constexpr std::uint64_t entryDomain = 0x5357'4f50'5045'4e54ULL;

class Fingerprint final {
public:
    template<typename Unsigned> void integer(Unsigned value) {
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
    template<typename Enum> void enumeration(const Enum value) {
        integer(static_cast<std::underlying_type_t<Enum>>(value));
    }
    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }
private:
    std::uint64_t value_ = fnvOffsetBasis;
};

struct DirectedEntry {
    std::size_t row = 0;
    std::size_t column = 0;
    PlanarPressureRegionFragmentOpeningPressureEntryKind kind =
        PlanarPressureRegionFragmentOpeningPressureEntryKind::SameRegionGrid;
    std::size_t sourceIndex = 0;
    std::uint64_t sourceStableId = 0;
    double weight = 0.0;
};

std::size_t checkedMultiply(const std::size_t a, const std::size_t b) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a)
        throw std::length_error("opening pressure-operator storage overflows");
    return a * b;
}
std::size_t checkedAdd(const std::size_t a, const std::size_t b) {
    if (b > std::numeric_limits<std::size_t>::max() - a)
        throw std::length_error("opening pressure-operator storage overflows");
    return a + b;
}
template<typename... Values> std::size_t storageSum(const Values... values) {
    std::size_t result = 0;
    ((result = checkedAdd(result, values)), ...);
    return result;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits) {
    if (limits.maximumRows == 0 || limits.maximumEntries == 0
        || limits.maximumComponents == 0 || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0)
        throw std::invalid_argument("opening pressure-operator limits are invalid");
}

std::size_t ownedBytes(
    const PlanarPressureRegionFragmentOpeningPressureOperator& op) {
    return storageSum(
        checkedMultiply(op.rows.size(), sizeof(op.rows.front())),
        checkedMultiply(op.entries.size(), sizeof(op.entries.front())),
        checkedMultiply(op.components.size(), sizeof(op.components.front())),
        checkedMultiply(op.componentFragmentIndices.size(), sizeof(std::size_t)));
}

std::uint64_t rowId(const std::uint64_t fragmentId) {
    Fingerprint f; f.integer(rowDomain); f.integer(fragmentId); return f.value();
}
std::uint64_t entryId(
    const PlanarPressureRegionFragmentOpeningPressureEntryKind kind,
    const std::uint64_t sourceId, const std::uint64_t rowFragmentId,
    const std::uint64_t columnFragmentId) {
    Fingerprint f; f.integer(entryDomain); f.enumeration(kind);
    f.integer(sourceId); f.integer(rowFragmentId); f.integer(columnFragmentId);
    return f.value();
}

std::uint64_t operatorFingerprint(
    const PlanarPressureRegionFragmentOpeningPressureOperator& op) {
    Fingerprint f;
    f.integer(op.version); f.integer(op.sourceBaseOperatorFingerprint);
    f.integer(op.sourceOpeningFingerprint); f.integer(op.sourceFragmentFingerprint);
    f.integer(op.sourceTopologyFingerprint);
    f.integer(static_cast<std::uint64_t>(op.rows.size()));
    for (const auto& row : op.rows) {
        for (const std::size_t v : {row.rowIndex, row.fragmentIndex,
             row.connectedComponentIndex, row.firstEntry, row.entryCount})
            f.integer(static_cast<std::uint64_t>(v));
        f.integer(row.stableId); f.integer(row.fragmentStableId);
        f.integer(static_cast<std::uint8_t>(row.isGauge));
        f.real(row.diagonalGeometryWeightMeters);
    }
    f.integer(static_cast<std::uint64_t>(op.entries.size()));
    for (const auto& entry : op.entries) {
        f.integer(static_cast<std::uint64_t>(entry.entryIndex));
        f.integer(entry.stableId); f.enumeration(entry.kind);
        f.integer(static_cast<std::uint64_t>(entry.sourceIndex));
        f.integer(entry.sourceStableId);
        f.integer(static_cast<std::uint64_t>(entry.columnFragmentIndex));
        f.integer(entry.columnFragmentStableId);
        f.real(entry.geometryWeightMeters);
    }
    f.integer(static_cast<std::uint64_t>(op.components.size()));
    for (const auto& c : op.components) {
        for (const std::size_t v : {c.componentIndex, c.baseComponentCount,
             c.gaugeFragmentIndex, c.firstFragmentMember, c.fragmentCount})
            f.integer(static_cast<std::uint64_t>(v));
        f.integer(c.stableId); f.real(c.totalVolumeCubicMeters);
        f.real(c.totalGeometryWeightMeters);
    }
    f.integer(static_cast<std::uint64_t>(op.componentFragmentIndices.size()));
    for (const auto v : op.componentFragmentIndices)
        f.integer(static_cast<std::uint64_t>(v));
    f.integer(static_cast<std::uint64_t>(op.includedSameRegionGridLinkCount));
    f.integer(static_cast<std::uint64_t>(op.includedOpeningPatchCount));
    for (const double v : {op.totalPressureLayerWallAreaSquareMeters,
         op.totalOpeningAreaSquareMeters, op.totalSolidPressureLayerWallAreaSquareMeters,
         op.wallAreaPartitionResidualSquareMeters, op.sameRegionGeometryWeightMeters,
         op.openingGeometryWeightMeters, op.totalGeometryWeightMeters,
         op.totalDiagonalGeometryWeightMeters}) f.real(v);
    f.integer(static_cast<std::uint64_t>(op.ownedStorageBytes));
    f.integer(static_cast<std::uint64_t>(op.workingStorageBytes));
    return f.value();
}

PlanarPressureRegionFragmentOpeningPressureOperator buildOperator(
    const PlanarPressureRegionFragmentPressureOperator& base,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits) {
    validateLimits(limits);
    const std::size_t rowCount = base.rows.size();
    const std::size_t entryCount = checkedAdd(
        base.entries.size(), checkedMultiply(openings.patches.size(), 2));
    const std::size_t componentCount = openings.connectedComponents.size();
    if (rowCount > limits.maximumRows || entryCount > limits.maximumEntries
        || componentCount > limits.maximumComponents)
        throw std::length_error("opening pressure-operator count limit exceeded");
    const std::size_t expectedOwned = storageSum(
        checkedMultiply(rowCount, sizeof(PlanarPressureRegionFragmentOpeningPressureOperatorRow)),
        checkedMultiply(entryCount, sizeof(PlanarPressureRegionFragmentOpeningPressureOperatorEntry)),
        checkedMultiply(componentCount, sizeof(PlanarPressureRegionFragmentOpeningPressureOperatorComponent)),
        checkedMultiply(rowCount, sizeof(std::size_t)));
    const std::size_t working = checkedMultiply(entryCount, sizeof(DirectedEntry));
    if (expectedOwned > limits.maximumOwnedBytes || working > limits.maximumWorkingBytes)
        throw std::length_error("opening pressure-operator byte limit exceeded");

    std::vector<DirectedEntry> directed;
    directed.reserve(entryCount);
    for (const auto& entry : base.entries) {
        const std::size_t row = base.rows[entry.columnFragmentIndex].fragmentIndex;
        // Recover this directed entry's owning row from the base row ranges.
        const auto owner = std::ranges::upper_bound(
            base.rows, entry.entryIndex, {},
            [](const auto& r) { return r.firstEntry; });
        const std::size_t ownerIndex = static_cast<std::size_t>(
            std::distance(base.rows.begin(), owner)) - 1;
        directed.push_back({ownerIndex, row,
            PlanarPressureRegionFragmentOpeningPressureEntryKind::SameRegionGrid,
            entry.sourceFaceLinkIndex, entry.sourceFaceLinkStableId,
            entry.geometryWeightMeters});
    }
    for (const auto& patch : openings.patches) {
        const double weight = patch.areaSquareMeters / patch.centerDistanceMeters;
        if (!std::isfinite(weight) || !(weight > 0.0))
            throw std::overflow_error("opening pressure-operator patch weight is invalid");
        directed.push_back({patch.minusFragmentIndex, patch.plusFragmentIndex,
            PlanarPressureRegionFragmentOpeningPressureEntryKind::OpeningPatch,
            patch.patchIndex, patch.patchStableId, weight});
        directed.push_back({patch.plusFragmentIndex, patch.minusFragmentIndex,
            PlanarPressureRegionFragmentOpeningPressureEntryKind::OpeningPatch,
            patch.patchIndex, patch.patchStableId, weight});
    }
    std::ranges::sort(directed, [](const auto& a, const auto& b) {
        return std::tuple{a.row, a.column, a.kind, a.sourceStableId}
            < std::tuple{b.row, b.column, b.kind, b.sourceStableId};
    });

    PlanarPressureRegionFragmentOpeningPressureOperator result;
    result.sourceBaseOperatorFingerprint = base.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.rows.resize(rowCount); result.entries.reserve(entryCount);
    for (std::size_t i = 0; i < rowCount; ++i) {
        const auto& source = topology.fragments[i];
        result.rows[i] = {i, rowId(source.fragmentStableId), i,
            source.fragmentStableId,
            openings.baseComponents[source.componentIndex].connectedComponentIndex,
            false, result.entries.size(), 0, 0.0};
        while (result.entries.size() < directed.size()
               && directed[result.entries.size()].row == i) {
            const auto& d = directed[result.entries.size()];
            const auto index = result.entries.size();
            result.entries.push_back({index,
                entryId(d.kind, d.sourceStableId, source.fragmentStableId,
                    topology.fragments[d.column].fragmentStableId),
                d.kind, d.sourceIndex, d.sourceStableId, d.column,
                topology.fragments[d.column].fragmentStableId, d.weight});
            ++result.rows[i].entryCount;
            result.rows[i].diagonalGeometryWeightMeters += d.weight;
        }
    }
    result.components.resize(componentCount);
    for (std::size_t c = 0; c < componentCount; ++c) {
        const auto& source = openings.connectedComponents[c];
        auto& out = result.components[c];
        out.componentIndex = c; out.stableId = source.stableId;
        out.baseComponentCount = source.baseComponentCount;
        out.firstFragmentMember = result.componentFragmentIndices.size();
        out.totalVolumeCubicMeters = source.volumeCubicMeters;
        out.gaugeFragmentIndex = std::numeric_limits<std::size_t>::max();
        for (std::size_t i = 0; i < rowCount; ++i) {
            if (result.rows[i].connectedComponentIndex != c) continue;
            result.componentFragmentIndices.push_back(i); ++out.fragmentCount;
            if (out.gaugeFragmentIndex == std::numeric_limits<std::size_t>::max()
                || result.rows[i].fragmentStableId
                    < result.rows[out.gaugeFragmentIndex].fragmentStableId)
                out.gaugeFragmentIndex = i;
            out.totalGeometryWeightMeters +=
                0.5 * result.rows[i].diagonalGeometryWeightMeters;
        }
        result.rows[out.gaugeFragmentIndex].isGauge = true;
    }
    result.includedSameRegionGridLinkCount = base.includedSameRegionGridLinkCount;
    result.includedOpeningPatchCount = openings.patches.size();
    for (const auto& link : topology.links)
        if (link.kind == PlanarPressureRegionFragmentFaceKind::PressureLayerWall)
            result.totalPressureLayerWallAreaSquareMeters += link.areaSquareMeters;
    result.totalOpeningAreaSquareMeters = openings.totalOpeningAreaSquareMeters;
    result.totalSolidPressureLayerWallAreaSquareMeters =
        result.totalPressureLayerWallAreaSquareMeters - result.totalOpeningAreaSquareMeters;
    result.wallAreaPartitionResidualSquareMeters =
        result.totalPressureLayerWallAreaSquareMeters
        - result.totalOpeningAreaSquareMeters
        - result.totalSolidPressureLayerWallAreaSquareMeters;
    result.sameRegionGeometryWeightMeters = base.totalGeometryWeightMeters;
    for (const auto& patch : openings.patches)
        result.openingGeometryWeightMeters +=
            patch.areaSquareMeters / patch.centerDistanceMeters;
    result.totalGeometryWeightMeters = result.sameRegionGeometryWeightMeters
        + result.openingGeometryWeightMeters;
    result.totalDiagonalGeometryWeightMeters = 2.0 * result.totalGeometryWeightMeters;
    result.ownedStorageBytes = ownedBytes(result); result.workingStorageBytes = working;
    if (result.ownedStorageBytes != expectedOwned)
        throw std::logic_error("opening pressure-operator storage is inconsistent");
    result.fingerprint = operatorFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureOperator
buildPlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& base,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits) {
    validatePlanarPressureRegionFragmentPressureOperator(
        base, grid, sweep, fragments, topology, limits.baseOperatorLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, definitions, limits.openingLimits);
    auto result = buildOperator(base, fragments, topology, openings, limits);
    validatePlanarPressureRegionFragmentOpeningPressureOperator(
        result, base, grid, sweep, fragments, topology, definitions, openings, limits);
    return result;
}

void validatePlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentOpeningPressureOperator& op,
    const PlanarPressureRegionFragmentPressureOperator& base,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits) {
    validatePlanarPressureRegionFragmentPressureOperator(
        base, grid, sweep, fragments, topology, limits.baseOperatorLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, definitions, limits.openingLimits);
    validateLimits(limits);
    if (op.version != planarPressureRegionFragmentOpeningPressureOperatorVersion
        || op.fingerprint == 0 || op.sourceBaseOperatorFingerprint != base.fingerprint
        || op.sourceOpeningFingerprint != openings.fingerprint
        || op.sourceFragmentFingerprint != fragments.fingerprint
        || op.sourceTopologyFingerprint != topology.fingerprint
        || op.ownedStorageBytes != ownedBytes(op)
        || op.fingerprint != operatorFingerprint(op)
        || op != buildOperator(base, fragments, topology, openings, limits))
        throw std::invalid_argument("opening pressure-operator is invalid");
}

std::vector<double> applyPlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentOpeningPressureOperator& op,
    const std::span<const double> pressure) {
    if (pressure.size() != op.rows.size()
        || !std::ranges::all_of(pressure, [](const double v) { return std::isfinite(v); }))
        throw std::invalid_argument("opening pressure-operator input is invalid");
    std::vector<double> result(op.rows.size(), 0.0);
    for (const auto& row : op.rows) {
        result[row.rowIndex] = row.diagonalGeometryWeightMeters * pressure[row.rowIndex];
        for (std::size_t e = row.firstEntry; e < row.firstEntry + row.entryCount; ++e)
            result[row.rowIndex] -= op.entries[e].geometryWeightMeters
                * pressure[op.entries[e].columnFragmentIndex];
        if (!std::isfinite(result[row.rowIndex]))
            throw std::overflow_error("opening pressure-operator result is non-finite");
    }
    return result;
}

} // namespace simwing::fsi::fluid
