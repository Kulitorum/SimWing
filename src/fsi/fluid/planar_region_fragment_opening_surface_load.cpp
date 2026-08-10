#include "fluid/planar_region_fragment_opening_surface_load.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    if (limits.maximumTiles == 0 || limits.maximumSurfaces == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening-aware regional surface-load limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening-aware regional surface-load storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening-aware regional surface-load storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger) {
    return checkedSum(
        checkedProduct(
            ledger.tiles.size(),
            sizeof(PlanarPressureRegionFragmentOpeningSurfaceLoadTile)),
        checkedProduct(
            ledger.surfaces.size(),
            sizeof(PlanarPressureRegionFragmentOpeningSurfaceLoadSummary)));
}

std::size_t workingStorageBytes(const std::size_t tileCount) {
    return checkedSum(
        checkedProduct(
            tileCount,
            sizeof(std::pair<std::size_t, std::size_t>)),
        checkedProduct(tileCount, sizeof(std::size_t)));
}

Vector3 scaledVector(const Vector3& value, const double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vector3 vectorSum(const Vector3& first, const Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

Vector3 vectorDifference(const Vector3& first, const Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

Vector3 crossProduct(const Vector3& first, const Vector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

void addVector(Vector3& target, const Vector3& value) {
    target.x += value.x;
    target.y += value.y;
    target.z += value.z;
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double closureTolerance(const double scale) {
    return std::max(
        1.0e-12,
        128.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, std::abs(scale)));
}

struct VectorPartition {
    Vector3 opening;
    Vector3 solid;
};

VectorPartition partitionVector(const Vector3& source,
                                const double openingArea,
                                const double solidArea,
                                const double wallArea) {
    if (openingArea == 0.0) return {{}, source};
    if (solidArea == 0.0) return {source, {}};
    const Vector3 opening = scaledVector(source, openingArea / wallArea);
    return {opening, vectorDifference(source, opening)};
}

struct ScalarPartition {
    double opening = 0.0;
    double solid = 0.0;
};

ScalarPartition partitionScalar(const double source,
                                const double openingArea,
                                const double solidArea,
                                const double wallArea) {
    if (openingArea == 0.0) return {0.0, source};
    if (solidArea == 0.0) return {source, 0.0};
    const double opening = source * (openingArea / wallArea);
    return {opening, source - opening};
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void fingerprintTile(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadTile& tile) {
    for (const std::size_t value : {
             tile.tileIndex,
             tile.sourceSurfaceLoadTileIndex,
             tile.sourcePressureWallIndex,
             tile.sourceFaceLinkIndex}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(tile.sourceFaceLinkStableId);
    fingerprint.integer(tile.surfaceStableId);
    fingerprint.enumeration(tile.axis);
    fingerprint.integer(tile.minusRegionStableId);
    fingerprint.integer(tile.plusRegionStableId);
    fingerprint.boolean(tile.touchedByOpening);
    fingerprint.integer(static_cast<std::uint64_t>(
        tile.sourceOpeningPartitionIndex));
    fingerprint.integer(static_cast<std::uint64_t>(
        tile.openingPatchCount));
    fingerprint.real(tile.wallAreaSquareMeters);
    fingerprint.real(tile.openingAreaSquareMeters);
    fingerprint.real(tile.solidAreaSquareMeters);
    fingerprint.real(tile.openingAreaFraction);
    fingerprint.boolean(tile.hasExactSubtileCentroids);
    fingerprintVector(fingerprint, tile.wrappedCentroidMeters);
    fingerprintVector(
        fingerprint, tile.openingAreaWeightedCentroidMeters);
    fingerprintVector(
        fingerprint, tile.solidAreaWeightedCentroidMeters);
    fingerprintVector(fingerprint, tile.unitNormalMinusToPlus);
    fingerprintVector(
        fingerprint, tile.authoredPressureTractionOnSheetPascals);
    fingerprintVector(
        fingerprint, tile.correctionPressureTractionOnSheetPascals);
    fingerprintVector(
        fingerprint, tile.totalPressureTractionOnSheetPascals);
    fingerprintVector(
        fingerprint,
        tile.openingRemovedAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        tile.openingRemovedCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        tile.openingRemovedTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, tile.solidAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, tile.solidCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, tile.solidTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        tile.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint, tile.solidTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint,
        tile.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
    fingerprintVector(
        fingerprint, tile.solidTotalPressureMomentOnSheetNewtonMeters);
    fingerprint.real(tile.openingRemovedTotalPressureWorkToSheetJoules);
    fingerprint.real(tile.solidTotalPressureWorkToSheetJoules);
}

void fingerprintSurface(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadSummary& surface) {
    fingerprint.integer(static_cast<std::uint64_t>(surface.surfaceIndex));
    fingerprint.integer(static_cast<std::uint64_t>(
        surface.sourceSurfaceLoadSurfaceIndex));
    fingerprint.integer(surface.surfaceStableId);
    fingerprint.enumeration(surface.axis);
    fingerprint.integer(surface.minusRegionStableId);
    fingerprint.integer(surface.plusRegionStableId);
    fingerprint.integer(static_cast<std::uint64_t>(surface.tileCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        surface.openingTouchedTileCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        surface.fullyOpenTileCount));
    fingerprint.real(surface.wallAreaSquareMeters);
    fingerprint.real(surface.openingAreaSquareMeters);
    fingerprint.real(surface.solidAreaSquareMeters);
    fingerprintVector(
        fingerprint, surface.openingAreaWeightedCentroidMeters);
    fingerprintVector(
        fingerprint, surface.solidAreaWeightedCentroidMeters);
    fingerprintVector(
        fingerprint,
        surface.openingRemovedAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        surface.openingRemovedCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        surface.openingRemovedTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, surface.solidAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, surface.solidCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, surface.solidTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        surface.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint,
        surface.solidTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint,
        surface.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
    fingerprintVector(
        fingerprint,
        surface.solidTotalPressureMomentOnSheetNewtonMeters);
    fingerprint.real(surface.openingRemovedTotalPressureWorkToSheetJoules);
    fingerprint.real(surface.solidTotalPressureWorkToSheetJoules);
}

std::uint64_t ledgerFingerprint(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger) {
    Fingerprint fingerprint;
    fingerprint.integer(ledger.version);
    fingerprint.integer(ledger.sourceSurfaceLoadFingerprint);
    fingerprint.integer(ledger.sourcePressureStateFingerprint);
    fingerprint.integer(ledger.sourceOpeningFingerprint);
    fingerprint.integer(ledger.sourceTopologyFingerprint);
    fingerprint.boolean(ledger.staticGeometry);
    fingerprint.boolean(ledger.usesMovingVolumeRates);
    fingerprint.real(ledger.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(ledger.tiles.size()));
    for (const auto& tile : ledger.tiles) {
        fingerprintTile(fingerprint, tile);
    }
    fingerprint.integer(static_cast<std::uint64_t>(ledger.surfaces.size()));
    for (const auto& surface : ledger.surfaces) {
        fingerprintSurface(fingerprint, surface);
    }
    fingerprint.real(ledger.totalWallAreaSquareMeters);
    fingerprint.real(ledger.totalOpeningAreaSquareMeters);
    fingerprint.real(ledger.totalSolidAreaSquareMeters);
    fingerprint.real(ledger.wallAreaPartitionResidualSquareMeters);
    fingerprintVector(
        fingerprint, ledger.sourceAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.sourceCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.sourceTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        ledger.openingRemovedAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        ledger.openingRemovedCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        ledger.openingRemovedTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.solidAuthoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.solidCorrectionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.solidTotalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint,
        ledger.sourceTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint,
        ledger.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint,
        ledger.solidTotalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint, ledger.sourceTotalPressureMomentOnSheetNewtonMeters);
    fingerprintVector(
        fingerprint,
        ledger.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
    fingerprintVector(
        fingerprint,
        ledger.solidTotalPressureMomentOnSheetNewtonMeters);
    fingerprint.real(ledger.sourceTotalPressureWorkToSheetJoules);
    fingerprint.real(ledger.openingRemovedTotalPressureWorkToSheetJoules);
    fingerprint.real(ledger.solidTotalPressureWorkToSheetJoules);
    fingerprint.real(
        ledger.maximumAbsoluteAreaPartitionResidualSquareMeters);
    fingerprint.real(ledger.maximumAbsoluteForcePartitionResidualNewtons);
    fingerprint.real(ledger.maximumAbsoluteForceSplitResidualNewtons);
    fingerprint.real(
        ledger.maximumAbsoluteImpulsePartitionResidualNewtonSeconds);
    fingerprint.real(
        ledger.maximumAbsoluteMomentPartitionResidualNewtonMeters);
    fingerprint.real(
        ledger.maximumAbsoluteSurfaceAggregationResidualNewtons);
    fingerprint.real(ledger.workPartitionResidualJoules);
    fingerprint.boolean(ledger.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        ledger.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        ledger.workingStorageBytes));
    return fingerprint.value();
}

void addTileToSurface(
    PlanarPressureRegionFragmentOpeningSurfaceLoadSummary& surface,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadTile& tile) {
    ++surface.tileCount;
    if (tile.touchedByOpening) ++surface.openingTouchedTileCount;
    if (tile.solidAreaSquareMeters == 0.0) ++surface.fullyOpenTileCount;
    surface.wallAreaSquareMeters += tile.wallAreaSquareMeters;
    surface.openingAreaSquareMeters += tile.openingAreaSquareMeters;
    surface.solidAreaSquareMeters += tile.solidAreaSquareMeters;
    addVector(
        surface.openingAreaWeightedCentroidMeters,
        scaledVector(
            tile.openingAreaWeightedCentroidMeters,
            tile.openingAreaSquareMeters));
    addVector(
        surface.solidAreaWeightedCentroidMeters,
        scaledVector(
            tile.solidAreaWeightedCentroidMeters,
            tile.solidAreaSquareMeters));
    addVector(
        surface.openingRemovedAuthoredPressureForceOnSheetNewtons,
        tile.openingRemovedAuthoredPressureForceOnSheetNewtons);
    addVector(
        surface.openingRemovedCorrectionPressureForceOnSheetNewtons,
        tile.openingRemovedCorrectionPressureForceOnSheetNewtons);
    addVector(
        surface.openingRemovedTotalPressureForceOnSheetNewtons,
        tile.openingRemovedTotalPressureForceOnSheetNewtons);
    addVector(
        surface.solidAuthoredPressureForceOnSheetNewtons,
        tile.solidAuthoredPressureForceOnSheetNewtons);
    addVector(
        surface.solidCorrectionPressureForceOnSheetNewtons,
        tile.solidCorrectionPressureForceOnSheetNewtons);
    addVector(
        surface.solidTotalPressureForceOnSheetNewtons,
        tile.solidTotalPressureForceOnSheetNewtons);
    addVector(
        surface.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds,
        tile.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
    addVector(
        surface.solidTotalPressureImpulseOnSheetNewtonSeconds,
        tile.solidTotalPressureImpulseOnSheetNewtonSeconds);
    addVector(
        surface.openingRemovedTotalPressureMomentOnSheetNewtonMeters,
        tile.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
    addVector(
        surface.solidTotalPressureMomentOnSheetNewtonMeters,
        tile.solidTotalPressureMomentOnSheetNewtonMeters);
    surface.openingRemovedTotalPressureWorkToSheetJoules +=
        tile.openingRemovedTotalPressureWorkToSheetJoules;
    surface.solidTotalPressureWorkToSheetJoules +=
        tile.solidTotalPressureWorkToSheetJoules;
}

PlanarPressureRegionFragmentOpeningSurfaceLoadLedger buildLedger(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    validateLimits(limits);
    if (surfaceLoads.tiles.size() > limits.maximumTiles
        || surfaceLoads.surfaces.size() > limits.maximumSurfaces) {
        throw std::length_error(
            "opening-aware regional surface-load entity limit exceeded");
    }
    if (surfaceLoads.sourceTopologyFingerprint
            != openings.sourceTopologyFingerprint
        || surfaceLoads.sourceTopologyFingerprint == 0) {
        throw std::invalid_argument(
            "opening-aware regional surface-load sources disagree");
    }

    PlanarPressureRegionFragmentOpeningSurfaceLoadLedger result;
    result.sourceSurfaceLoadFingerprint = surfaceLoads.fingerprint;
    result.sourcePressureStateFingerprint =
        surfaceLoads.sourcePressureStateFingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceTopologyFingerprint =
        surfaceLoads.sourceTopologyFingerprint;
    result.staticGeometry = surfaceLoads.staticGeometry;
    result.usesMovingVolumeRates = surfaceLoads.usesMovingVolumeRates;
    result.timeStepSeconds = surfaceLoads.timeStepSeconds;
    result.tiles.reserve(surfaceLoads.tiles.size());
    result.surfaces.resize(surfaceLoads.surfaces.size());
    result.ownedStorageBytes = checkedSum(
        checkedProduct(
            surfaceLoads.tiles.size(),
            sizeof(PlanarPressureRegionFragmentOpeningSurfaceLoadTile)),
        checkedProduct(
            surfaceLoads.surfaces.size(),
            sizeof(PlanarPressureRegionFragmentOpeningSurfaceLoadSummary)));
    result.workingStorageBytes = workingStorageBytes(
        surfaceLoads.tiles.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening-aware regional surface-load storage limit exceeded");
    }

    for (std::size_t index = 0;
         index < surfaceLoads.surfaces.size(); ++index) {
        const auto& source = surfaceLoads.surfaces[index];
        auto& surface = result.surfaces[index];
        surface.surfaceIndex = index;
        surface.sourceSurfaceLoadSurfaceIndex = source.surfaceIndex;
        surface.surfaceStableId = source.surfaceStableId;
        surface.axis = source.axis;
        surface.minusRegionStableId = source.minusRegionStableId;
        surface.plusRegionStableId = source.plusRegionStableId;
    }

    using TileLookup = std::pair<std::size_t, std::size_t>;
    std::vector<TileLookup> tileLookup;
    tileLookup.reserve(surfaceLoads.tiles.size());
    for (std::size_t index = 0; index < surfaceLoads.tiles.size(); ++index) {
        tileLookup.emplace_back(
            surfaceLoads.tiles[index].sourceFaceLinkIndex, index);
    }
    std::ranges::sort(tileLookup, {}, &TileLookup::first);
    if (std::ranges::adjacent_find(
            tileLookup,
            [](const auto& first, const auto& second) {
                return first.first == second.first;
            }) != tileLookup.end()) {
        throw std::invalid_argument(
            "opening-aware regional surface-load source links repeat");
    }
    const std::size_t noPartition = openings.partitions.size();
    std::vector<std::size_t> partitionByTile(
        surfaceLoads.tiles.size(), noPartition);
    for (std::size_t index = 0;
         index < openings.partitions.size(); ++index) {
        const auto& partition = openings.partitions[index];
        const auto found = std::ranges::lower_bound(
            tileLookup, partition.sourceFaceLinkIndex, {},
            &TileLookup::first);
        if (found == tileLookup.end()
            || found->first != partition.sourceFaceLinkIndex
            || partitionByTile[found->second] != noPartition) {
            throw std::invalid_argument(
                "opening-aware regional surface-load partition has no unique tile");
        }
        const auto& tile = surfaceLoads.tiles[found->second];
        if (partition.sourceFaceLinkStableId
                != tile.sourceFaceLinkStableId
            || partition.surfaceStableId != tile.surfaceStableId
            || partition.axis != tile.axis
            || partition.negativeSideRegionStableId
                != tile.minusRegionStableId
            || partition.positiveSideRegionStableId
                != tile.plusRegionStableId
            || partition.wallAreaSquareMeters != tile.areaSquareMeters) {
            throw std::invalid_argument(
                "opening-aware regional surface-load partition disagrees with its tile");
        }
        partitionByTile[found->second] = index;
    }

    double maximumForceScale = 0.0;
    double maximumMomentScale = 0.0;
    double maximumWorkScale = 0.0;
    for (std::size_t index = 0; index < surfaceLoads.tiles.size(); ++index) {
        const auto& source = surfaceLoads.tiles[index];
        const bool touched = partitionByTile[index] != noPartition;
        const auto* partition = touched
            ? &openings.partitions[partitionByTile[index]]
            : nullptr;
        const double openingArea = touched
            ? partition->openingAreaSquareMeters : 0.0;
        const double solidArea = touched
            ? partition->solidAreaSquareMeters : source.areaSquareMeters;
        const double openingFraction = touched
            ? partition->openingAreaFraction : 0.0;
        if (!std::isfinite(openingArea) || !std::isfinite(solidArea)
            || !std::isfinite(openingFraction)
            || openingArea < 0.0 || solidArea < 0.0
            || openingFraction < 0.0 || openingFraction > 1.0) {
            throw std::invalid_argument(
                "opening-aware regional surface-load area is invalid");
        }

        const auto authored = partitionVector(
            source.authoredPressureForceOnSheetNewtons,
            openingArea, solidArea, source.areaSquareMeters);
        const auto correction = partitionVector(
            source.correctionPressureForceOnSheetNewtons,
            openingArea, solidArea, source.areaSquareMeters);
        const auto total = partitionVector(
            source.totalPressureForceOnSheetNewtons,
            openingArea, solidArea, source.areaSquareMeters);
        const auto impulse = partitionVector(
            source.totalPressureImpulseOnSheetNewtonSeconds,
            openingArea, solidArea, source.areaSquareMeters);
        const auto work = partitionScalar(
            source.totalPressureWorkToSheetJoules,
            openingArea, solidArea, source.areaSquareMeters);
        const bool exactSubtileCentroids = touched
            && partition->hasExactSubtileCentroids;
        const Vector3 openingCentroid = touched
            ? partition->openingAreaWeightedCentroidMeters
            : source.wrappedCentroidMeters;
        const Vector3 solidCentroid = touched
            ? partition->solidAreaWeightedCentroidMeters
            : source.wrappedCentroidMeters;
        const Vector3 sourceMoment = crossProduct(
            source.wrappedCentroidMeters,
            source.totalPressureForceOnSheetNewtons);
        const Vector3 openingMoment = crossProduct(
            openingCentroid, total.opening);
        const Vector3 solidMoment = crossProduct(
            solidCentroid, total.solid);
        if (!finiteVector(authored.opening) || !finiteVector(authored.solid)
            || !finiteVector(correction.opening)
            || !finiteVector(correction.solid)
            || !finiteVector(total.opening) || !finiteVector(total.solid)
            || !finiteVector(impulse.opening)
            || !finiteVector(impulse.solid)
            || !finiteVector(sourceMoment) || !finiteVector(openingMoment)
            || !finiteVector(solidMoment)
            || !std::isfinite(work.opening)
            || !std::isfinite(work.solid)) {
            throw std::overflow_error(
                "opening-aware regional surface-load partition is non-finite");
        }

        result.tiles.push_back({
            result.tiles.size(),
            source.tileIndex,
            source.sourcePressureWallIndex,
            source.sourceFaceLinkIndex,
            source.sourceFaceLinkStableId,
            source.surfaceStableId,
            source.axis,
            source.minusRegionStableId,
            source.plusRegionStableId,
            touched,
            touched ? partition->partitionIndex : 0,
            touched ? partition->openingPatchCount : 0,
            source.areaSquareMeters,
            openingArea,
            solidArea,
            openingFraction,
            exactSubtileCentroids,
            source.wrappedCentroidMeters,
            openingCentroid,
            solidCentroid,
            source.unitNormalMinusToPlus,
            source.authoredPressureTractionOnSheetPascals,
            source.correctionPressureTractionOnSheetPascals,
            source.totalPressureTractionOnSheetPascals,
            authored.opening,
            correction.opening,
            total.opening,
            authored.solid,
            correction.solid,
            total.solid,
            impulse.opening,
            impulse.solid,
            openingMoment,
            solidMoment,
            work.opening,
            work.solid,
        });
        const auto& tile = result.tiles.back();

        result.totalWallAreaSquareMeters += tile.wallAreaSquareMeters;
        result.totalOpeningAreaSquareMeters += tile.openingAreaSquareMeters;
        result.totalSolidAreaSquareMeters += tile.solidAreaSquareMeters;
        addVector(
            result.sourceTotalPressureMomentOnSheetNewtonMeters,
            sourceMoment);
        addVector(
            result.openingRemovedAuthoredPressureForceOnSheetNewtons,
            tile.openingRemovedAuthoredPressureForceOnSheetNewtons);
        addVector(
            result.openingRemovedCorrectionPressureForceOnSheetNewtons,
            tile.openingRemovedCorrectionPressureForceOnSheetNewtons);
        addVector(
            result.openingRemovedTotalPressureForceOnSheetNewtons,
            tile.openingRemovedTotalPressureForceOnSheetNewtons);
        addVector(
            result.solidAuthoredPressureForceOnSheetNewtons,
            tile.solidAuthoredPressureForceOnSheetNewtons);
        addVector(
            result.solidCorrectionPressureForceOnSheetNewtons,
            tile.solidCorrectionPressureForceOnSheetNewtons);
        addVector(
            result.solidTotalPressureForceOnSheetNewtons,
            tile.solidTotalPressureForceOnSheetNewtons);
        addVector(
            result.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds,
            tile.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
        addVector(
            result.solidTotalPressureImpulseOnSheetNewtonSeconds,
            tile.solidTotalPressureImpulseOnSheetNewtonSeconds);
        addVector(
            result.openingRemovedTotalPressureMomentOnSheetNewtonMeters,
            tile.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
        addVector(
            result.solidTotalPressureMomentOnSheetNewtonMeters,
            tile.solidTotalPressureMomentOnSheetNewtonMeters);
        result.openingRemovedTotalPressureWorkToSheetJoules +=
            tile.openingRemovedTotalPressureWorkToSheetJoules;
        result.solidTotalPressureWorkToSheetJoules +=
            tile.solidTotalPressureWorkToSheetJoules;

        const auto surfacePosition = std::ranges::lower_bound(
            surfaceLoads.surfaces, source.surfaceStableId, {},
            &PlanarPressureRegionFragmentSurfaceLoadSummary::surfaceStableId);
        if (surfacePosition == surfaceLoads.surfaces.end()
            || surfacePosition->surfaceStableId != source.surfaceStableId) {
            throw std::invalid_argument(
                "opening-aware regional surface-load tile has no surface");
        }
        const std::size_t surfaceIndex = static_cast<std::size_t>(
            surfacePosition - surfaceLoads.surfaces.begin());
        addTileToSurface(result.surfaces[surfaceIndex], tile);

        result.maximumAbsoluteAreaPartitionResidualSquareMeters = std::max(
            result.maximumAbsoluteAreaPartitionResidualSquareMeters,
            std::abs(source.areaSquareMeters - openingArea - solidArea));
        for (const auto& force : {
                 std::pair{source.authoredPressureForceOnSheetNewtons,
                           vectorSum(authored.opening, authored.solid)},
                 std::pair{source.correctionPressureForceOnSheetNewtons,
                           vectorSum(correction.opening, correction.solid)},
                 std::pair{source.totalPressureForceOnSheetNewtons,
                           vectorSum(total.opening, total.solid)},
                 std::pair{total.opening,
                           vectorSum(authored.opening, correction.opening)},
                 std::pair{total.solid,
                           vectorSum(authored.solid, correction.solid)},
                 std::pair{total.solid,
                           scaledVector(
                               source.totalPressureTractionOnSheetPascals,
                               solidArea)}}) {
            result.maximumAbsoluteForcePartitionResidualNewtons = std::max(
                result.maximumAbsoluteForcePartitionResidualNewtons,
                maximumAbsoluteComponent(
                    vectorDifference(force.first, force.second)));
        }
        result.maximumAbsoluteForceSplitResidualNewtons = std::max({
            result.maximumAbsoluteForceSplitResidualNewtons,
            maximumAbsoluteComponent(vectorDifference(
                total.opening,
                vectorSum(authored.opening, correction.opening))),
            maximumAbsoluteComponent(vectorDifference(
                total.solid,
                vectorSum(authored.solid, correction.solid))),
        });
        result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds =
            std::max({
                result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds,
                maximumAbsoluteComponent(vectorDifference(
                    source.totalPressureImpulseOnSheetNewtonSeconds,
                    vectorSum(impulse.opening, impulse.solid))),
                maximumAbsoluteComponent(vectorDifference(
                    impulse.opening,
                    scaledVector(total.opening, result.timeStepSeconds))),
                maximumAbsoluteComponent(vectorDifference(
                    impulse.solid,
                    scaledVector(total.solid, result.timeStepSeconds))),
            });
        result.maximumAbsoluteMomentPartitionResidualNewtonMeters = std::max(
            result.maximumAbsoluteMomentPartitionResidualNewtonMeters,
            maximumAbsoluteComponent(vectorDifference(
                sourceMoment, vectorSum(openingMoment, solidMoment))));
        result.workPartitionResidualJoules +=
            source.totalPressureWorkToSheetJoules
            - work.opening - work.solid;
        maximumForceScale = std::max({
            maximumForceScale,
            maximumAbsoluteComponent(
                source.authoredPressureForceOnSheetNewtons),
            maximumAbsoluteComponent(
                source.correctionPressureForceOnSheetNewtons),
            maximumAbsoluteComponent(
                source.totalPressureForceOnSheetNewtons),
        });
        maximumMomentScale = std::max(
            maximumMomentScale,
            maximumAbsoluteComponent(sourceMoment));
        maximumWorkScale = std::max(
            maximumWorkScale,
            std::abs(source.totalPressureWorkToSheetJoules));
    }

    result.sourceAuthoredPressureForceOnSheetNewtons =
        surfaceLoads.authoredPressureForceOnSheetNewtons;
    result.sourceCorrectionPressureForceOnSheetNewtons =
        surfaceLoads.correctionPressureForceOnSheetNewtons;
    result.sourceTotalPressureForceOnSheetNewtons =
        surfaceLoads.totalPressureForceOnSheetNewtons;
    result.sourceTotalPressureImpulseOnSheetNewtonSeconds =
        surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds;
    result.sourceTotalPressureWorkToSheetJoules =
        surfaceLoads.totalPressureWorkToSheetJoules;
    result.wallAreaPartitionResidualSquareMeters =
        result.totalWallAreaSquareMeters
        - result.totalOpeningAreaSquareMeters
        - result.totalSolidAreaSquareMeters;

    Vector3 surfaceOpeningAuthored;
    Vector3 surfaceOpeningCorrection;
    Vector3 surfaceOpeningTotal;
    Vector3 surfaceSolidAuthored;
    Vector3 surfaceSolidCorrection;
    Vector3 surfaceSolidTotal;
    Vector3 surfaceOpeningImpulse;
    Vector3 surfaceSolidImpulse;
    Vector3 surfaceOpeningMoment;
    Vector3 surfaceSolidMoment;
    double surfaceOpeningWork = 0.0;
    double surfaceSolidWork = 0.0;
    double surfaceWallArea = 0.0;
    double surfaceOpeningArea = 0.0;
    double surfaceSolidArea = 0.0;
    for (auto& surface : result.surfaces) {
        if (surface.tileCount == 0) {
            throw std::invalid_argument(
                "opening-aware regional surface-load surface is empty");
        }
        if (surface.openingAreaSquareMeters > 0.0) {
            surface.openingAreaWeightedCentroidMeters = scaledVector(
                surface.openingAreaWeightedCentroidMeters,
                1.0 / surface.openingAreaSquareMeters);
        }
        if (surface.solidAreaSquareMeters > 0.0) {
            surface.solidAreaWeightedCentroidMeters = scaledVector(
                surface.solidAreaWeightedCentroidMeters,
                1.0 / surface.solidAreaSquareMeters);
        }
        surfaceWallArea += surface.wallAreaSquareMeters;
        surfaceOpeningArea += surface.openingAreaSquareMeters;
        surfaceSolidArea += surface.solidAreaSquareMeters;
        addVector(surfaceOpeningAuthored,
                  surface.openingRemovedAuthoredPressureForceOnSheetNewtons);
        addVector(surfaceOpeningCorrection,
                  surface.openingRemovedCorrectionPressureForceOnSheetNewtons);
        addVector(surfaceOpeningTotal,
                  surface.openingRemovedTotalPressureForceOnSheetNewtons);
        addVector(surfaceSolidAuthored,
                  surface.solidAuthoredPressureForceOnSheetNewtons);
        addVector(surfaceSolidCorrection,
                  surface.solidCorrectionPressureForceOnSheetNewtons);
        addVector(surfaceSolidTotal,
                  surface.solidTotalPressureForceOnSheetNewtons);
        addVector(surfaceOpeningImpulse,
                  surface.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds);
        addVector(surfaceSolidImpulse,
                  surface.solidTotalPressureImpulseOnSheetNewtonSeconds);
        addVector(surfaceOpeningMoment,
                  surface.openingRemovedTotalPressureMomentOnSheetNewtonMeters);
        addVector(surfaceSolidMoment,
                  surface.solidTotalPressureMomentOnSheetNewtonMeters);
        surfaceOpeningWork +=
            surface.openingRemovedTotalPressureWorkToSheetJoules;
        surfaceSolidWork += surface.solidTotalPressureWorkToSheetJoules;
    }

    const auto updateSurfaceResidual = [&](const Vector3& first,
                                           const Vector3& second) {
        result.maximumAbsoluteSurfaceAggregationResidualNewtons = std::max(
            result.maximumAbsoluteSurfaceAggregationResidualNewtons,
            maximumAbsoluteComponent(vectorDifference(first, second)));
    };
    updateSurfaceResidual(
        surfaceOpeningAuthored,
        result.openingRemovedAuthoredPressureForceOnSheetNewtons);
    updateSurfaceResidual(
        surfaceOpeningCorrection,
        result.openingRemovedCorrectionPressureForceOnSheetNewtons);
    updateSurfaceResidual(
        surfaceOpeningTotal,
        result.openingRemovedTotalPressureForceOnSheetNewtons);
    updateSurfaceResidual(
        surfaceSolidAuthored,
        result.solidAuthoredPressureForceOnSheetNewtons);
    updateSurfaceResidual(
        surfaceSolidCorrection,
        result.solidCorrectionPressureForceOnSheetNewtons);
    updateSurfaceResidual(
        surfaceSolidTotal,
        result.solidTotalPressureForceOnSheetNewtons);
    result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds = std::max({
        result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds,
        maximumAbsoluteComponent(vectorDifference(
            surfaceOpeningImpulse,
            result.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds)),
        maximumAbsoluteComponent(vectorDifference(
            surfaceSolidImpulse,
            result.solidTotalPressureImpulseOnSheetNewtonSeconds)),
    });

    const double areaTolerance = closureTolerance(
        result.totalWallAreaSquareMeters);
    const double forceTolerance = closureTolerance(maximumForceScale);
    const double impulseTolerance =
        forceTolerance * result.timeStepSeconds;
    const double momentTolerance = closureTolerance(maximumMomentScale);
    const double workTolerance = closureTolerance(std::max(
        maximumWorkScale,
        std::abs(result.sourceTotalPressureWorkToSheetJoules)));
    const Vector3 authoredPartitionResidual = vectorDifference(
        result.sourceAuthoredPressureForceOnSheetNewtons,
        vectorSum(
            result.openingRemovedAuthoredPressureForceOnSheetNewtons,
            result.solidAuthoredPressureForceOnSheetNewtons));
    const Vector3 correctionPartitionResidual = vectorDifference(
        result.sourceCorrectionPressureForceOnSheetNewtons,
        vectorSum(
            result.openingRemovedCorrectionPressureForceOnSheetNewtons,
            result.solidCorrectionPressureForceOnSheetNewtons));
    const Vector3 totalPartitionResidual = vectorDifference(
        result.sourceTotalPressureForceOnSheetNewtons,
        vectorSum(
            result.openingRemovedTotalPressureForceOnSheetNewtons,
            result.solidTotalPressureForceOnSheetNewtons));
    const Vector3 impulsePartitionResidual = vectorDifference(
        result.sourceTotalPressureImpulseOnSheetNewtonSeconds,
        vectorSum(
            result.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds,
            result.solidTotalPressureImpulseOnSheetNewtonSeconds));
    const Vector3 momentPartitionResidual = vectorDifference(
        result.sourceTotalPressureMomentOnSheetNewtonMeters,
        vectorSum(
            result.openingRemovedTotalPressureMomentOnSheetNewtonMeters,
            result.solidTotalPressureMomentOnSheetNewtonMeters));
    result.maximumAbsoluteForcePartitionResidualNewtons = std::max({
        result.maximumAbsoluteForcePartitionResidualNewtons,
        maximumAbsoluteComponent(authoredPartitionResidual),
        maximumAbsoluteComponent(correctionPartitionResidual),
        maximumAbsoluteComponent(totalPartitionResidual),
    });
    result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds = std::max(
        result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds,
        maximumAbsoluteComponent(impulsePartitionResidual));
    result.maximumAbsoluteMomentPartitionResidualNewtonMeters = std::max(
        result.maximumAbsoluteMomentPartitionResidualNewtonMeters,
        maximumAbsoluteComponent(momentPartitionResidual));

    if (!finiteVector(result.sourceAuthoredPressureForceOnSheetNewtons)
        || !finiteVector(result.sourceCorrectionPressureForceOnSheetNewtons)
        || !finiteVector(result.sourceTotalPressureForceOnSheetNewtons)
        || !finiteVector(
            result.openingRemovedAuthoredPressureForceOnSheetNewtons)
        || !finiteVector(
            result.openingRemovedCorrectionPressureForceOnSheetNewtons)
        || !finiteVector(
            result.openingRemovedTotalPressureForceOnSheetNewtons)
        || !finiteVector(result.solidAuthoredPressureForceOnSheetNewtons)
        || !finiteVector(result.solidCorrectionPressureForceOnSheetNewtons)
        || !finiteVector(result.solidTotalPressureForceOnSheetNewtons)
        || !finiteVector(
            result.sourceTotalPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(
            result.openingRemovedTotalPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(
            result.solidTotalPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(result.sourceTotalPressureMomentOnSheetNewtonMeters)
        || !finiteVector(
            result.openingRemovedTotalPressureMomentOnSheetNewtonMeters)
        || !finiteVector(
            result.solidTotalPressureMomentOnSheetNewtonMeters)
        || !std::isfinite(result.sourceTotalPressureWorkToSheetJoules)
        || !std::isfinite(
            result.openingRemovedTotalPressureWorkToSheetJoules)
        || !std::isfinite(result.solidTotalPressureWorkToSheetJoules)
        || !std::isfinite(result.wallAreaPartitionResidualSquareMeters)
        || std::abs(result.wallAreaPartitionResidualSquareMeters)
            > areaTolerance
        || std::abs(result.totalWallAreaSquareMeters
                    - surfaceLoads.totalAreaSquareMeters)
            > areaTolerance
        || std::abs(result.totalOpeningAreaSquareMeters
                    - openings.totalOpeningAreaSquareMeters)
            > areaTolerance
        || std::abs(surfaceWallArea - result.totalWallAreaSquareMeters)
            > areaTolerance
        || std::abs(surfaceOpeningArea - result.totalOpeningAreaSquareMeters)
            > areaTolerance
        || std::abs(surfaceSolidArea - result.totalSolidAreaSquareMeters)
            > areaTolerance
        || result.maximumAbsoluteAreaPartitionResidualSquareMeters
            > areaTolerance
        || result.maximumAbsoluteForcePartitionResidualNewtons
            > forceTolerance
        || result.maximumAbsoluteForceSplitResidualNewtons > forceTolerance
        || result.maximumAbsoluteImpulsePartitionResidualNewtonSeconds
            > impulseTolerance
        || result.maximumAbsoluteMomentPartitionResidualNewtonMeters
            > momentTolerance
        || result.maximumAbsoluteSurfaceAggregationResidualNewtons
            > forceTolerance
        || maximumAbsoluteComponent(vectorDifference(
               surfaceOpeningMoment,
               result.openingRemovedTotalPressureMomentOnSheetNewtonMeters))
            > momentTolerance
        || maximumAbsoluteComponent(vectorDifference(
               surfaceSolidMoment,
               result.solidTotalPressureMomentOnSheetNewtonMeters))
            > momentTolerance
        || std::abs(surfaceOpeningWork
                    - result.openingRemovedTotalPressureWorkToSheetJoules)
            > workTolerance
        || std::abs(surfaceSolidWork
                    - result.solidTotalPressureWorkToSheetJoules)
            > workTolerance
        || std::abs(result.workPartitionResidualJoules) > workTolerance
        || std::abs(result.sourceTotalPressureWorkToSheetJoules
                    - result.openingRemovedTotalPressureWorkToSheetJoules
                    - result.solidTotalPressureWorkToSheetJoules)
            > workTolerance) {
        throw std::invalid_argument(
            "opening-aware regional surface-load closure failed");
    }

    result.accepted = true;
    result.fingerprint = ledgerFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningSurfaceLoadLedger
capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState, limits.surfaceLoadLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, openingDefinitions,
        limits.openingLimits);
    const auto result = buildLedger(surfaceLoads, openings, limits);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        result, surfaceLoads, pressureState, grid, sweep, fragments,
        topology, openingDefinitions, openings, limits);
    return result;
}

PlanarPressureRegionFragmentOpeningSurfaceLoadLedger
capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState, limits.surfaceLoadLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, openingDefinitions,
        limits.openingLimits);
    const auto result = buildLedger(surfaceLoads, openings, limits);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        result, surfaceLoads, pressureState, grid, sweep, fragments,
        topology, openingDefinitions, openings, limits);
    return result;
}

void validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger) {
    if (ledger.version
            != planarPressureRegionFragmentOpeningSurfaceLoadVersion
        || ledger.fingerprint == 0 || !ledger.accepted
        || ledger.sourceSurfaceLoadFingerprint == 0
        || ledger.sourcePressureStateFingerprint == 0
        || ledger.sourceOpeningFingerprint == 0
        || ledger.sourceTopologyFingerprint == 0
        || !std::isfinite(ledger.timeStepSeconds)
        || !(ledger.timeStepSeconds > 0.0)
        || !std::isfinite(ledger.totalWallAreaSquareMeters)
        || !std::isfinite(ledger.totalOpeningAreaSquareMeters)
        || !std::isfinite(ledger.totalSolidAreaSquareMeters)
        || !std::isfinite(ledger.wallAreaPartitionResidualSquareMeters)
        || !finiteVector(ledger.sourceTotalPressureForceOnSheetNewtons)
        || !finiteVector(ledger.solidTotalPressureForceOnSheetNewtons)
        || !finiteVector(
            ledger.openingRemovedTotalPressureForceOnSheetNewtons)
        || !std::isfinite(ledger.sourceTotalPressureWorkToSheetJoules)
        || !std::isfinite(ledger.solidTotalPressureWorkToSheetJoules)
        || !std::isfinite(
            ledger.openingRemovedTotalPressureWorkToSheetJoules)
        || ledger.ownedStorageBytes != ownedStorageBytes(ledger)
        || ledger.workingStorageBytes
            != workingStorageBytes(ledger.tiles.size())
        || ledger.fingerprint != ledgerFingerprint(ledger)) {
        throw std::invalid_argument(
            "opening-aware regional surface-load ledger integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState, limits.surfaceLoadLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, openingDefinitions,
        limits.openingLimits);
    if (ledger.tiles.size() > limits.maximumTiles
        || ledger.surfaces.size() > limits.maximumSurfaces
        || ledger.ownedStorageBytes > limits.maximumOwnedBytes
        || ledger.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening-aware regional surface-load validation limit exceeded");
    }
    validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
        ledger);
    if (ledger != buildLedger(surfaceLoads, openings, limits)) {
        throw std::invalid_argument(
            "opening-aware regional surface-load ledger is corrupted");
    }
}

void validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState, limits.surfaceLoadLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, openingDefinitions,
        limits.openingLimits);
    if (ledger.tiles.size() > limits.maximumTiles
        || ledger.surfaces.size() > limits.maximumSurfaces
        || ledger.ownedStorageBytes > limits.maximumOwnedBytes
        || ledger.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening-aware regional surface-load validation limit exceeded");
    }
    validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
        ledger);
    if (ledger != buildLedger(surfaceLoads, openings, limits)) {
        throw std::invalid_argument(
            "opening-aware regional surface-load ledger is corrupted");
    }
}

} // namespace simwing::fsi::fluid
