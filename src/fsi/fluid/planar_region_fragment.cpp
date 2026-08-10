#include "fluid/planar_region_fragment.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t fragmentIdentityDomain = 0x5357'5246'5241'4731ULL;

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

    void signedInteger(const std::int64_t value) {
        integer(std::bit_cast<std::uint64_t>(value));
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

struct AxisGeometry {
    std::size_t cellCount = 0;
    double lowerMeters = 0.0;
    double upperMeters = 0.0;
    double spacingMeters = 0.0;
    double domainLengthMeters = 0.0;
    double transverseCellAreaSquareMeters = 0.0;
    std::size_t transverseTileCount = 0;
};

void validateLimits(const PlanarPressureRegionFragmentLimits& limits) {
    if (limits.maximumIntervals == 0 || limits.maximumRegions == 0
        || limits.maximumCells == 0 || limits.maximumFragments == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region fragment limits are invalid");
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

std::size_t ownedStorageBytes(const std::size_t fragmentCount,
                              const std::size_t regionCount,
                              const std::size_t cellCount) {
    constexpr const char* message =
        "planar pressure region fragment storage size overflows";
    const std::size_t fragments = checkedProduct(
        fragmentCount, sizeof(PlanarPressureRegionFragment), message);
    const std::size_t regions = checkedProduct(
        regionCount, sizeof(PlanarPressureRegionFragmentRegionSummary),
        message);
    const std::size_t cells = checkedProduct(
        cellCount, sizeof(PlanarPressureRegionFragmentCellSummary),
        message);
    return checkedSum(checkedSum(fragments, regions, message),
                      cells, message);
}

AxisGeometry axisGeometry(const PeriodicCartesianGrid& grid,
                          const GridFaceAxis axis) {
    const GridCellCounts counts = grid.cellCounts();
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    AxisGeometry result;
    switch (axis) {
    case GridFaceAxis::X:
        result = {
            counts.x, lower.x, upper.x, spacing.x, upper.x - lower.x,
            spacing.y * spacing.z,
            checkedProduct(
                counts.y, counts.z,
                "planar pressure region fragment tile count overflows"),
        };
        break;
    case GridFaceAxis::Y:
        result = {
            counts.y, lower.y, upper.y, spacing.y, upper.y - lower.y,
            spacing.x * spacing.z,
            checkedProduct(
                counts.x, counts.z,
                "planar pressure region fragment tile count overflows"),
        };
        break;
    case GridFaceAxis::Z:
        result = {
            counts.z, lower.z, upper.z, spacing.z, upper.z - lower.z,
            spacing.x * spacing.y,
            checkedProduct(
                counts.x, counts.y,
                "planar pressure region fragment tile count overflows"),
        };
        break;
    default:
        throw std::invalid_argument(
            "planar pressure region fragment axis is invalid");
    }
    if (result.cellCount == 0
        || !std::isfinite(result.transverseCellAreaSquareMeters)
        || !(result.transverseCellAreaSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region fragment axis geometry is invalid");
    }
    return result;
}

std::int64_t floorToInt64(const double value) {
    if (!std::isfinite(value)
        || value < static_cast<double>(
            std::numeric_limits<std::int64_t>::min())
        || value > static_cast<double>(
            std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error(
            "planar pressure region fragment cell ordinal overflows");
    }
    return static_cast<std::int64_t>(std::floor(value));
}

std::int64_t floorDiv(const std::int64_t numerator,
                      const std::int64_t positiveDenominator) {
    std::int64_t quotient = numerator / positiveDenominator;
    const std::int64_t remainder = numerator % positiveDenominator;
    if (remainder < 0) --quotient;
    return quotient;
}

std::size_t positiveModulo(const std::int64_t value,
                           const std::int64_t positiveModulus) {
    std::int64_t remainder = value % positiveModulus;
    if (remainder < 0) remainder += positiveModulus;
    return static_cast<std::size_t>(remainder);
}

PlanarPressureRegionFragmentBoundary layerBoundary(
    const std::uint64_t surfaceStableId,
    const double coordinateMeters) {
    return {
        PlanarPressureRegionFragmentBoundaryKind::PressureLayer,
        surfaceStableId,
        0,
        0,
        coordinateMeters,
    };
}

PlanarPressureRegionFragmentBoundary gridBoundary(
    const std::int64_t boundaryOrdinal,
    const std::size_t axisCellCount,
    const double coordinateMeters) {
    const auto count = static_cast<std::int64_t>(axisCellCount);
    return {
        PlanarPressureRegionFragmentBoundaryKind::GridFace,
        0,
        positiveModulo(boundaryOrdinal, count),
        floorDiv(boundaryOrdinal, count),
        coordinateMeters,
    };
}

std::uint64_t fragmentStableId(
    const GridFaceAxis axis,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k,
    const std::uint64_t regionStableId,
    const PlanarPressureRegionFragmentBoundary& lower,
    const PlanarPressureRegionFragmentBoundary& upper) {
    Fingerprint fingerprint;
    fingerprint.integer(fragmentIdentityDomain);
    fingerprint.enumeration(axis);
    fingerprint.integer(static_cast<std::uint64_t>(i));
    fingerprint.integer(static_cast<std::uint64_t>(j));
    fingerprint.integer(static_cast<std::uint64_t>(k));
    fingerprint.integer(regionStableId);
    for (const auto* boundary : {&lower, &upper}) {
        fingerprint.enumeration(boundary->kind);
        fingerprint.integer(boundary->surfaceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            boundary->faceCoordinate));
        fingerprint.signedInteger(boundary->periodicImage);
    }
    return fingerprint.value();
}

void assignAxisCoordinate(Vector3& value,
                          const GridFaceAxis axis,
                          const double coordinate) {
    switch (axis) {
    case GridFaceAxis::X: value.x = coordinate; break;
    case GridFaceAxis::Y: value.y = coordinate; break;
    case GridFaceAxis::Z: value.z = coordinate; break;
    }
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({std::abs(value.x), std::abs(value.y),
                     std::abs(value.z)});
}

std::uint64_t setFingerprint(
    const PlanarPressureRegionFragmentSet& set) {
    Fingerprint fingerprint;
    fingerprint.integer(set.version);
    fingerprint.integer(set.sourceSweepVersion);
    fingerprint.enumeration(set.axis);
    fingerprint.integer(static_cast<std::uint64_t>(set.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(set.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(set.cellCounts.z));
    for (const double value : {
             set.lowerMeters.x, set.lowerMeters.y, set.lowerMeters.z,
             set.upperMeters.x, set.upperMeters.y, set.upperMeters.z,
             set.spacingMeters.x, set.spacingMeters.y, set.spacingMeters.z,
             set.profileWindowLowerCoordinateMeters,
             set.profileWindowUpperCoordinateMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.fragments.size()));
    for (const auto& fragment : set.fragments) {
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(fragment.i));
        fingerprint.integer(static_cast<std::uint64_t>(fragment.j));
        fingerprint.integer(static_cast<std::uint64_t>(fragment.k));
        fingerprint.signedInteger(fragment.axisCellPeriodicImage);
        for (const auto* boundary
             : {&fragment.lowerBoundary, &fragment.upperBoundary}) {
            fingerprint.enumeration(boundary->kind);
            fingerprint.integer(boundary->surfaceStableId);
            fingerprint.integer(static_cast<std::uint64_t>(
                boundary->faceCoordinate));
            fingerprint.signedInteger(boundary->periodicImage);
            fingerprint.real(boundary->unwrappedCoordinateMeters);
        }
        for (const double value : {
                 fragment.unwrappedLowerCoordinateMeters,
                 fragment.unwrappedUpperCoordinateMeters,
                 fragment.unwrappedAxisCentroidMeters,
                 fragment.transverseAreaSquareMeters,
                 fragment.volumeCubicMeters,
                 fragment.wrappedCentroidMeters.x,
                 fragment.wrappedCentroidMeters.y,
                 fragment.wrappedCentroidMeters.z,
                 fragment.pressurePascals}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.regions.size()));
    for (const auto& region : set.regions) {
        fingerprint.integer(region.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(region.fragmentCount));
        fingerprint.real(region.volumeCubicMeters);
        fingerprint.real(region.sourceProfileVolumeCubicMeters);
        fingerprint.real(region.volumeClosureResidualCubicMeters);
        fingerprint.real(region.pressurePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.cells.size()));
    for (const auto& cell : set.cells) {
        fingerprint.integer(static_cast<std::uint64_t>(cell.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(cell.fragmentCount));
        for (const double value : {
                 cell.fragmentVolumeCubicMeters,
                 cell.volumeClosureResidualCubicMeters,
                 cell.fragmentFirstMomentCubicMetersSquared.x,
                 cell.fragmentFirstMomentCubicMetersSquared.y,
                 cell.fragmentFirstMomentCubicMetersSquared.z,
                 cell.firstMomentClosureResidualCubicMetersSquared.x,
                 cell.firstMomentClosureResidualCubicMetersSquared.y,
                 cell.firstMomentClosureResidualCubicMetersSquared.z}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        set.maximumFragmentsPerCell));
    fingerprint.real(set.fragmentVolumeCubicMeters);
    fingerprint.real(set.domainVolumeClosureResidualCubicMeters);
    fingerprint.real(
        set.maximumAbsoluteRegionVolumeClosureResidualCubicMeters);
    fingerprint.real(
        set.maximumAbsoluteCellVolumeClosureResidualCubicMeters);
    fingerprint.real(
        set.maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared);
    fingerprint.integer(static_cast<std::uint64_t>(set.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentSet buildSet(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionSweepLedger(
        sweep,
        {limits.maximumIntervals, limits.maximumRegions,
         std::numeric_limits<std::size_t>::max()});
    if (grid.cellCount() > limits.maximumCells) {
        throw std::length_error(
            "planar pressure region fragment cell limit exceeded");
    }
    const auto& profile = sweep.currentProfile;
    const AxisGeometry geometry = axisGeometry(grid, profile.axis);
    const double windowLength =
        profile.windowUpperCoordinateMeters
        - profile.windowLowerCoordinateMeters;
    const double domainScale = std::max(
        std::abs(windowLength), std::abs(geometry.domainLengthMeters));
    const double coordinateTolerance = 64.0
        * std::numeric_limits<double>::epsilon() * domainScale;
    if (sweep.axis != profile.axis
        || std::abs(windowLength - geometry.domainLengthMeters)
            > coordinateTolerance
        || geometry.cellCount
            > static_cast<std::size_t>(
                std::numeric_limits<std::int64_t>::max())) {
        throw std::invalid_argument(
            "planar pressure region fragment profile/grid identity is invalid");
    }

    PlanarPressureRegionFragmentSet result;
    result.sourceSweepVersion = sweep.version;
    result.axis = profile.axis;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.spacingMeters = grid.cellSpacingMeters();
    result.profileWindowLowerCoordinateMeters =
        profile.windowLowerCoordinateMeters;
    result.profileWindowUpperCoordinateMeters =
        profile.windowUpperCoordinateMeters;
    result.cells.resize(grid.cellCount());
    for (std::size_t k = 0; k < result.cellCounts.z; ++k) {
        for (std::size_t j = 0; j < result.cellCounts.y; ++j) {
            for (std::size_t i = 0; i < result.cellCounts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                result.cells[index].cellIndex = index;
                result.cells[index].i = i;
                result.cells[index].j = j;
                result.cells[index].k = k;
            }
        }
    }

    std::set<std::uint64_t> fragmentIds;
    for (const auto& interval : profile.intervals) {
        std::vector<double> points;
        points.push_back(interval.lowerCoordinateMeters);
        std::int64_t boundaryOrdinal = floorToInt64(
            (interval.lowerCoordinateMeters - geometry.lowerMeters)
            / geometry.spacingMeters) + 1;
        for (;;) {
            const double boundaryCoordinate = geometry.lowerMeters
                + static_cast<double>(boundaryOrdinal)
                    * geometry.spacingMeters;
            if (!(boundaryCoordinate < interval.upperCoordinateMeters)) {
                break;
            }
            if (!(boundaryCoordinate > points.back())) {
                throw std::invalid_argument(
                    "planar pressure region fragment grid split is not ordered");
            }
            points.push_back(boundaryCoordinate);
            if (boundaryOrdinal
                == std::numeric_limits<std::int64_t>::max()) {
                throw std::overflow_error(
                    "planar pressure region fragment boundary ordinal overflows");
            }
            ++boundaryOrdinal;
        }
        points.push_back(interval.upperCoordinateMeters);
        for (std::size_t segment = 0; segment + 1 < points.size();
             ++segment) {
            const double lowerCoordinate = points[segment];
            const double upperCoordinate = points[segment + 1];
            const double midpoint = 0.5 * lowerCoordinate
                + 0.5 * upperCoordinate;
            const std::int64_t cellOrdinal = floorToInt64(
                (midpoint - geometry.lowerMeters)
                / geometry.spacingMeters);
            const std::int64_t axisCount = static_cast<std::int64_t>(
                geometry.cellCount);
            const std::size_t wrappedAxisCell = positiveModulo(
                cellOrdinal, axisCount);
            const std::int64_t cellImage = floorDiv(
                cellOrdinal, axisCount);
            const auto lowerBoundary = segment == 0
                ? layerBoundary(
                    interval.lowerSurfaceStableId, lowerCoordinate)
                : gridBoundary(
                    cellOrdinal, geometry.cellCount, lowerCoordinate);
            const auto upperBoundary = segment + 2 == points.size()
                ? layerBoundary(
                    interval.upperSurfaceStableId, upperCoordinate)
                : gridBoundary(
                    cellOrdinal + 1, geometry.cellCount, upperCoordinate);
            const double axialWidth = upperCoordinate - lowerCoordinate;
            const double fragmentVolume = axialWidth
                * geometry.transverseCellAreaSquareMeters;
            if (!std::isfinite(fragmentVolume)
                || !(fragmentVolume > 0.0)) {
                throw std::invalid_argument(
                    "planar pressure region fragment volume is invalid");
            }
            const std::size_t prospective = checkedSum(
                result.fragments.size(), geometry.transverseTileCount,
                "planar pressure region fragment count overflows");
            if (prospective > limits.maximumFragments) {
                throw std::length_error(
                    "planar pressure region fragment count limit exceeded");
            }

            const auto append = [&](const std::size_t i,
                                    const std::size_t j,
                                    const std::size_t k) {
                PlanarPressureRegionFragment fragment;
                fragment.regionStableId = interval.regionStableId;
                fragment.i = i;
                fragment.j = j;
                fragment.k = k;
                fragment.axisCellPeriodicImage = cellImage;
                fragment.lowerBoundary = lowerBoundary;
                fragment.upperBoundary = upperBoundary;
                fragment.unwrappedLowerCoordinateMeters = lowerCoordinate;
                fragment.unwrappedUpperCoordinateMeters = upperCoordinate;
                fragment.unwrappedAxisCentroidMeters = midpoint;
                fragment.transverseAreaSquareMeters =
                    geometry.transverseCellAreaSquareMeters;
                fragment.volumeCubicMeters = fragmentVolume;
                fragment.wrappedCentroidMeters =
                    grid.cellCenterMeters(i, j, k);
                assignAxisCoordinate(
                    fragment.wrappedCentroidMeters, result.axis,
                    midpoint - static_cast<double>(cellImage)
                        * geometry.domainLengthMeters);
                fragment.pressurePascals = interval.pressurePascals;
                fragment.stableId = fragmentStableId(
                    result.axis, i, j, k, fragment.regionStableId,
                    lowerBoundary, upperBoundary);
                if (!fragmentIds.insert(fragment.stableId).second) {
                    throw std::invalid_argument(
                        "planar pressure region fragment stable ID collides");
                }
                result.fragments.push_back(fragment);
            };
            switch (result.axis) {
            case GridFaceAxis::X:
                for (std::size_t k = 0; k < result.cellCounts.z; ++k) {
                    for (std::size_t j = 0; j < result.cellCounts.y; ++j) {
                        append(wrappedAxisCell, j, k);
                    }
                }
                break;
            case GridFaceAxis::Y:
                for (std::size_t k = 0; k < result.cellCounts.z; ++k) {
                    for (std::size_t i = 0; i < result.cellCounts.x; ++i) {
                        append(i, wrappedAxisCell, k);
                    }
                }
                break;
            case GridFaceAxis::Z:
                for (std::size_t j = 0; j < result.cellCounts.y; ++j) {
                    for (std::size_t i = 0; i < result.cellCounts.x; ++i) {
                        append(i, j, wrappedAxisCell);
                    }
                }
                break;
            }
        }
    }

    std::map<std::uint64_t, PlanarPressureRegionFragmentRegionSummary>
        regions;
    for (const auto& fragment : result.fragments) {
        auto [entry, inserted] = regions.try_emplace(
            fragment.regionStableId);
        auto& region = entry->second;
        if (inserted) {
            region.regionStableId = fragment.regionStableId;
            region.pressurePascals = fragment.pressurePascals;
        } else if (region.pressurePascals != fragment.pressurePascals) {
            throw std::invalid_argument(
                "planar pressure region fragment pressure is inconsistent");
        }
        ++region.fragmentCount;
        region.volumeCubicMeters += fragment.volumeCubicMeters;
        result.fragmentVolumeCubicMeters += fragment.volumeCubicMeters;

        auto& cell = result.cells[grid.cellIndex(
            fragment.i, fragment.j, fragment.k)];
        ++cell.fragmentCount;
        cell.fragmentVolumeCubicMeters += fragment.volumeCubicMeters;
        cell.fragmentFirstMomentCubicMetersSquared.x +=
            fragment.volumeCubicMeters * fragment.wrappedCentroidMeters.x;
        cell.fragmentFirstMomentCubicMetersSquared.y +=
            fragment.volumeCubicMeters * fragment.wrappedCentroidMeters.y;
        cell.fragmentFirstMomentCubicMetersSquared.z +=
            fragment.volumeCubicMeters * fragment.wrappedCentroidMeters.z;
    }
    result.regions.reserve(regions.size());
    for (const auto& source : profile.regions) {
        const auto found = regions.find(source.regionStableId);
        if (found == regions.end()) {
            throw std::invalid_argument(
                "planar pressure region fragment region is missing");
        }
        auto region = found->second;
        region.sourceProfileVolumeCubicMeters = source.volumeCubicMeters;
        region.volumeClosureResidualCubicMeters =
            region.volumeCubicMeters
            - region.sourceProfileVolumeCubicMeters;
        if (region.pressurePascals != source.pressurePascals) {
            throw std::invalid_argument(
                "planar pressure region fragment source pressure changed");
        }
        result.maximumAbsoluteRegionVolumeClosureResidualCubicMeters =
            std::max(
                result.maximumAbsoluteRegionVolumeClosureResidualCubicMeters,
                std::abs(region.volumeClosureResidualCubicMeters));
        result.regions.push_back(region);
    }
    if (result.regions.size() != regions.size()) {
        throw std::invalid_argument(
            "planar pressure region fragment region identity changed");
    }

    const double cellVolume = grid.cellVolumeCubicMeters();
    for (auto& cell : result.cells) {
        const Vector3 center = grid.cellCenterMeters(cell.i, cell.j, cell.k);
        cell.volumeClosureResidualCubicMeters =
            cell.fragmentVolumeCubicMeters - cellVolume;
        cell.firstMomentClosureResidualCubicMetersSquared = {
            cell.fragmentFirstMomentCubicMetersSquared.x
                - cellVolume * center.x,
            cell.fragmentFirstMomentCubicMetersSquared.y
                - cellVolume * center.y,
            cell.fragmentFirstMomentCubicMetersSquared.z
                - cellVolume * center.z,
        };
        result.maximumFragmentsPerCell = std::max(
            result.maximumFragmentsPerCell, cell.fragmentCount);
        result.maximumAbsoluteCellVolumeClosureResidualCubicMeters =
            std::max(
                result.maximumAbsoluteCellVolumeClosureResidualCubicMeters,
                std::abs(cell.volumeClosureResidualCubicMeters));
        result
            .maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared =
            std::max(
                result
                    .maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared,
                maximumAbsoluteComponent(
                    cell.firstMomentClosureResidualCubicMetersSquared));
    }
    result.domainVolumeClosureResidualCubicMeters =
        result.fragmentVolumeCubicMeters
        - profile.geometricDomainVolumeCubicMeters;
    const double volumeTolerance = 1024.0
        * std::numeric_limits<double>::epsilon()
        * profile.geometricDomainVolumeCubicMeters;
    const double coordinateMagnitude = std::max({
        std::abs(result.lowerMeters.x), std::abs(result.lowerMeters.y),
        std::abs(result.lowerMeters.z), std::abs(result.upperMeters.x),
        std::abs(result.upperMeters.y), std::abs(result.upperMeters.z),
        1.0,
    });
    const double momentTolerance = volumeTolerance * coordinateMagnitude;
    if (std::abs(result.domainVolumeClosureResidualCubicMeters)
            > volumeTolerance
        || result.maximumAbsoluteRegionVolumeClosureResidualCubicMeters
            > volumeTolerance
        || result.maximumAbsoluteCellVolumeClosureResidualCubicMeters
            > volumeTolerance
        || result
            .maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared
            > momentTolerance) {
        throw std::invalid_argument(
            "planar pressure region fragment geometry does not close");
    }
    result.ownedStorageBytes = ownedStorageBytes(
        result.fragments.size(), result.regions.size(), result.cells.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region fragment byte limit exceeded");
    }
    result.fingerprint = setFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentSet buildPlanarPressureRegionFragments(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentLimits& limits) {
    auto result = buildSet(grid, sweep, limits);
    validatePlanarPressureRegionFragments(result, grid, sweep, limits);
    return result;
}

void validatePlanarPressureRegionFragments(
    const PlanarPressureRegionFragmentSet& fragments,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentLimits& limits) {
    validateLimits(limits);
    if (fragments.fragments.size() > limits.maximumFragments
        || fragments.regions.size() > limits.maximumRegions
        || fragments.cells.size() > limits.maximumCells) {
        throw std::length_error(
            "planar pressure region fragment validation limit exceeded");
    }
    const auto expected = buildSet(grid, sweep, limits);
    if (expected != fragments) {
        throw std::invalid_argument(
            "planar pressure region fragment set is invalid");
    }
}

} // namespace simwing::fsi::fluid
