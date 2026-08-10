#include "fluid/planar_region_fragment_topology.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <compare>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t linkIdentityDomain = 0x5357'5246'4c49'4e4bULL;
constexpr std::uint64_t componentIdentityDomain =
    0x5357'5246'434f'4d50ULL;

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

struct BoundaryIdentity {
    PlanarPressureRegionFragmentBoundaryKind kind =
        PlanarPressureRegionFragmentBoundaryKind::GridFace;
    std::uint64_t surfaceStableId = 0;
    std::size_t faceCoordinate = 0;
    std::int64_t periodicImage = 0;

    auto operator<=>(const BoundaryIdentity&) const = default;
};

struct HalfFaceKey {
    PlanarPressureRegionFragmentFaceKind kind =
        PlanarPressureRegionFragmentFaceKind::SameRegionGrid;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::int64_t facePeriodicImage = 0;
    std::uint64_t surfaceStableId = 0;
    std::uint64_t regionStableId = 0;
    std::int64_t profileCellPeriodicImage = 0;
    BoundaryIdentity profileLowerBoundary;
    BoundaryIdentity profileUpperBoundary;

    auto operator<=>(const HalfFaceKey&) const = default;
};

struct HalfFace {
    std::size_t fragmentIndex = 0;
    bool isMinusSide = false;
    double areaSquareMeters = 0.0;
    Vector3 wrappedCentroidMeters;
    double halfCenterDistanceMeters = 0.0;
};

class DisjointSet final {
public:
    explicit DisjointSet(const std::size_t count)
        : parents_(count), ranks_(count, 0) {
        std::iota(parents_.begin(), parents_.end(), 0);
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

struct ComponentWork {
    std::size_t root = 0;
    std::uint64_t regionStableId = 0;
    std::vector<std::uint64_t> fragmentStableIds;
    std::size_t sameRegionGridLinkCount = 0;
    double volumeCubicMeters = 0.0;
    std::uint64_t stableId = 0;
};

constexpr GridFaceAxis axes[] = {
    GridFaceAxis::X,
    GridFaceAxis::Y,
    GridFaceAxis::Z,
};

void validateLimits(
    const PlanarPressureRegionFragmentTopologyLimits& limits) {
    if (limits.maximumLinks == 0 || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region fragment topology limits are invalid");
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

std::size_t storageBytes(const std::size_t linkCount,
                         const std::size_t fragmentCount,
                         const std::size_t componentCount) {
    constexpr const char* message =
        "planar pressure region fragment topology storage overflows";
    const std::size_t links = checkedProduct(
        linkCount, sizeof(PlanarPressureRegionFragmentFaceLink), message);
    const std::size_t fragments = checkedProduct(
        fragmentCount,
        sizeof(PlanarPressureRegionFragmentTopologySummary), message);
    const std::size_t components = checkedProduct(
        componentCount,
        sizeof(PlanarPressureRegionFragmentComponent), message);
    return checkedSum(checkedSum(links, fragments, message),
                      components, message);
}

std::size_t axisCellCount(const GridCellCounts counts,
                          const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return counts.x;
    case GridFaceAxis::Y: return counts.y;
    case GridFaceAxis::Z: return counts.z;
    }
    throw std::invalid_argument(
        "planar pressure region fragment topology axis is invalid");
}

std::size_t fragmentCellCoordinate(
    const PlanarPressureRegionFragment& fragment,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return fragment.i;
    case GridFaceAxis::Y: return fragment.j;
    case GridFaceAxis::Z: return fragment.k;
    }
    throw std::invalid_argument(
        "planar pressure region fragment topology axis is invalid");
}

double vectorCoordinate(const Vector3 value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "planar pressure region fragment topology axis is invalid");
}

void assignVectorCoordinate(Vector3& value,
                            const GridFaceAxis axis,
                            const double coordinate) {
    switch (axis) {
    case GridFaceAxis::X: value.x = coordinate; break;
    case GridFaceAxis::Y: value.y = coordinate; break;
    case GridFaceAxis::Z: value.z = coordinate; break;
    }
}

void assignLatticeCoordinate(HalfFaceKey& key,
                             const GridFaceAxis axis,
                             const std::size_t coordinate) {
    switch (axis) {
    case GridFaceAxis::X: key.i = coordinate; break;
    case GridFaceAxis::Y: key.j = coordinate; break;
    case GridFaceAxis::Z: key.k = coordinate; break;
    }
}

void assignLatticeCoordinate(
    PlanarPressureRegionFragmentFaceLink& link,
    const GridFaceAxis axis,
    const std::size_t coordinate) {
    switch (axis) {
    case GridFaceAxis::X: link.i = coordinate; break;
    case GridFaceAxis::Y: link.j = coordinate; break;
    case GridFaceAxis::Z: link.k = coordinate; break;
    }
}

std::size_t latticeCoordinate(const HalfFaceKey& key,
                              const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return key.i;
    case GridFaceAxis::Y: return key.j;
    case GridFaceAxis::Z: return key.k;
    }
    throw std::invalid_argument(
        "planar pressure region fragment topology axis is invalid");
}

BoundaryIdentity boundaryIdentity(
    const PlanarPressureRegionFragmentBoundary& boundary) {
    return {
        boundary.kind,
        boundary.surfaceStableId,
        boundary.faceCoordinate,
        boundary.periodicImage,
    };
}

double periodicCoordinate(const double coordinate,
                          const double lower,
                          const double upper) {
    const double length = upper - lower;
    const double image = std::floor((coordinate - lower) / length);
    double wrapped = coordinate - image * length;
    if (wrapped < lower) wrapped += length;
    if (!(wrapped < upper)) wrapped -= length;
    return wrapped;
}

std::int64_t periodicImage(const double coordinate,
                           const double lower,
                           const double upper) {
    const double image = std::floor((coordinate - lower) / (upper - lower));
    if (!std::isfinite(image)
        || image < static_cast<double>(
            std::numeric_limits<std::int64_t>::min())
        || image > static_cast<double>(
            std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error(
            "planar pressure region fragment topology image overflows");
    }
    return static_cast<std::int64_t>(image);
}

Vector3 unitAxis(const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return {1.0, 0.0, 0.0};
    case GridFaceAxis::Y: return {0.0, 1.0, 0.0};
    case GridFaceAxis::Z: return {0.0, 0.0, 1.0};
    }
    throw std::invalid_argument(
        "planar pressure region fragment topology axis is invalid");
}

double fragmentAxisLength(
    const PlanarPressureRegionFragment& fragment,
    const GridFaceAxis axis,
    const GridFaceAxis profileAxis,
    const Vector3 spacing) {
    return axis == profileAxis
        ? fragment.unwrappedUpperCoordinateMeters
            - fragment.unwrappedLowerCoordinateMeters
        : vectorCoordinate(spacing, axis);
}

double fragmentFaceArea(
    const PlanarPressureRegionFragment& fragment,
    const GridFaceAxis faceAxis,
    const GridFaceAxis profileAxis,
    const Vector3 spacing) {
    double result = 1.0;
    for (const GridFaceAxis axis : axes) {
        if (axis != faceAxis) {
            result *= fragmentAxisLength(
                fragment, axis, profileAxis, spacing);
        }
    }
    return result;
}

double fragmentBoundaryArea(
    const PlanarPressureRegionFragment& fragment,
    const GridFaceAxis profileAxis,
    const Vector3 spacing) {
    const double x = fragmentAxisLength(
        fragment, GridFaceAxis::X, profileAxis, spacing);
    const double y = fragmentAxisLength(
        fragment, GridFaceAxis::Y, profileAxis, spacing);
    const double z = fragmentAxisLength(
        fragment, GridFaceAxis::Z, profileAxis, spacing);
    return 2.0 * (x * y + x * z + y * z);
}

Vector3 halfFaceCentroid(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionFragment& fragment,
    const GridFaceAxis faceAxis,
    const GridFaceAxis profileAxis,
    const std::size_t faceCoordinate,
    const PlanarPressureRegionFragmentBoundary* profileBoundary) {
    Vector3 result = grid.cellCenterMeters(
        fragment.i, fragment.j, fragment.k);
    if (faceAxis == profileAxis) {
        const Vector3 lower = grid.lowerMeters();
        const Vector3 upper = grid.upperMeters();
        assignVectorCoordinate(
            result, profileAxis,
            periodicCoordinate(
                profileBoundary->unwrappedCoordinateMeters,
                vectorCoordinate(lower, profileAxis),
                vectorCoordinate(upper, profileAxis)));
    } else {
        const Vector3 lower = grid.lowerMeters();
        const Vector3 spacing = grid.cellSpacingMeters();
        assignVectorCoordinate(
            result, faceAxis,
            vectorCoordinate(lower, faceAxis)
                + static_cast<double>(faceCoordinate)
                    * vectorCoordinate(spacing, faceAxis));
        assignVectorCoordinate(
            result, profileAxis,
            vectorCoordinate(fragment.wrappedCentroidMeters, profileAxis));
    }
    return result;
}

double maximumAbsoluteDifference(const Vector3 first,
                                 const Vector3 second) {
    return std::max({
        std::abs(first.x - second.x),
        std::abs(first.y - second.y),
        std::abs(first.z - second.z),
    });
}

std::uint64_t linkStableId(
    const PlanarPressureRegionFragmentFaceLink& link) {
    Fingerprint fingerprint;
    fingerprint.integer(linkIdentityDomain);
    fingerprint.enumeration(link.kind);
    fingerprint.enumeration(link.axis);
    fingerprint.integer(static_cast<std::uint64_t>(link.i));
    fingerprint.integer(static_cast<std::uint64_t>(link.j));
    fingerprint.integer(static_cast<std::uint64_t>(link.k));
    fingerprint.signedInteger(link.facePeriodicImage);
    fingerprint.integer(link.surfaceStableId);
    fingerprint.integer(link.minusFragmentStableId);
    fingerprint.integer(link.plusFragmentStableId);
    return fingerprint.value();
}

std::uint64_t componentStableId(
    const std::uint64_t regionStableId,
    std::vector<std::uint64_t> fragmentStableIds) {
    std::ranges::sort(fragmentStableIds);
    Fingerprint fingerprint;
    fingerprint.integer(componentIdentityDomain);
    fingerprint.integer(regionStableId);
    fingerprint.integer(static_cast<std::uint64_t>(
        fragmentStableIds.size()));
    for (const std::uint64_t stableId : fragmentStableIds) {
        fingerprint.integer(stableId);
    }
    return fingerprint.value();
}

std::uint64_t topologyFingerprint(
    const PlanarPressureRegionFragmentTopology& topology) {
    Fingerprint fingerprint;
    fingerprint.integer(topology.version);
    fingerprint.integer(topology.sourceFragmentFingerprint);
    fingerprint.enumeration(topology.profileAxis);
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.z));
    for (const double value : {
             topology.lowerMeters.x, topology.lowerMeters.y,
             topology.lowerMeters.z, topology.upperMeters.x,
             topology.upperMeters.y, topology.upperMeters.z,
             topology.spacingMeters.x, topology.spacingMeters.y,
             topology.spacingMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(topology.links.size()));
    for (const auto& link : topology.links) {
        fingerprint.integer(static_cast<std::uint64_t>(link.linkIndex));
        fingerprint.integer(link.stableId);
        fingerprint.enumeration(link.kind);
        fingerprint.enumeration(link.axis);
        fingerprint.integer(static_cast<std::uint64_t>(link.i));
        fingerprint.integer(static_cast<std::uint64_t>(link.j));
        fingerprint.integer(static_cast<std::uint64_t>(link.k));
        fingerprint.signedInteger(link.facePeriodicImage);
        fingerprint.integer(link.surfaceStableId);
        for (const std::size_t value : {
                 link.minusFragmentIndex, link.plusFragmentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        for (const std::uint64_t value : {
                 link.minusFragmentStableId, link.plusFragmentStableId,
                 link.minusRegionStableId, link.plusRegionStableId}) {
            fingerprint.integer(value);
        }
        fingerprint.integer(static_cast<std::uint64_t>(
            link.minusComponentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            link.plusComponentIndex));
        for (const double value : {
                 link.areaSquareMeters, link.wrappedCentroidMeters.x,
                 link.wrappedCentroidMeters.y,
                 link.wrappedCentroidMeters.z,
                 link.centerDistanceMeters,
                 link.sameRegionGeometryWeightMeters,
                 link.pressureJumpPascals,
                 link.unitNormalMinusToPlus.x,
                 link.unitNormalMinusToPlus.y,
                 link.unitNormalMinusToPlus.z}) {
            fingerprint.real(value);
        }
        fingerprint.integer(static_cast<std::uint8_t>(
            link.crossesPeriodicBoundary));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        topology.fragments.size()));
    for (const auto& fragment : topology.fragments) {
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.fragmentIndex));
        fingerprint.integer(fragment.fragmentStableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.incidentFaceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.sameRegionGridFaceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.pressureLayerWallFaceCount));
        fingerprint.real(fragment.incidentFaceAreaSquareMeters);
        fingerprint.real(fragment.expectedBoundaryAreaSquareMeters);
        fingerprint.real(
            fragment.boundaryAreaClosureResidualSquareMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        topology.components.size()));
    for (const auto& component : topology.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.fragmentCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.sameRegionGridLinkCount));
        fingerprint.real(component.volumeCubicMeters);
    }
    for (const std::size_t value : {
             topology.sameRegionGridLinkCount,
             topology.pressureLayerWallLinkCount,
             topology.periodicGridLinkCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    for (const double value : {
             topology.totalUniqueFaceAreaSquareMeters,
             topology.totalIncidentFaceAreaSquareMeters,
             topology.totalExpectedFragmentBoundaryAreaSquareMeters,
             topology
                 .maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        topology.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentTopology buildTopology(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& source,
    const PlanarPressureRegionFragmentTopologyLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragments(
        source, grid, sweep, limits.fragmentLimits);
    const std::size_t expectedLinks = checkedProduct(
        source.fragments.size(), 3,
        "planar pressure region fragment topology link count overflows");
    if (expectedLinks > limits.maximumLinks) {
        throw std::length_error(
            "planar pressure region fragment topology link limit exceeded");
    }

    PlanarPressureRegionFragmentTopology result;
    result.sourceFragmentFingerprint = source.fingerprint;
    result.profileAxis = source.axis;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.spacingMeters = grid.cellSpacingMeters();
    result.fragments.resize(source.fragments.size());
    for (std::size_t index = 0; index < source.fragments.size(); ++index) {
        const auto& fragment = source.fragments[index];
        auto& summary = result.fragments[index];
        summary.fragmentIndex = index;
        summary.fragmentStableId = fragment.stableId;
        summary.regionStableId = fragment.regionStableId;
        summary.expectedBoundaryAreaSquareMeters = fragmentBoundaryArea(
            fragment, source.axis, result.spacingMeters);
        if (!std::isfinite(summary.expectedBoundaryAreaSquareMeters)
            || !(summary.expectedBoundaryAreaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "planar pressure region fragment boundary area is invalid");
        }
    }

    std::map<HalfFaceKey, std::vector<HalfFace>> halfFaces;
    for (std::size_t fragmentIndex = 0;
         fragmentIndex < source.fragments.size(); ++fragmentIndex) {
        const auto& fragment = source.fragments[fragmentIndex];
        for (const GridFaceAxis faceAxis : axes) {
            for (const bool lowerSide : {true, false}) {
                HalfFaceKey key;
                key.axis = faceAxis;
                key.i = fragment.i;
                key.j = fragment.j;
                key.k = fragment.k;
                key.regionStableId = fragment.regionStableId;
                const PlanarPressureRegionFragmentBoundary* boundary =
                    nullptr;
                std::size_t faceCoordinate = 0;
                if (faceAxis == source.axis) {
                    boundary = lowerSide
                        ? &fragment.lowerBoundary
                        : &fragment.upperBoundary;
                    if (boundary->kind
                        == PlanarPressureRegionFragmentBoundaryKind::
                            PressureLayer) {
                        key.kind =
                            PlanarPressureRegionFragmentFaceKind::
                                PressureLayerWall;
                        key.surfaceStableId = boundary->surfaceStableId;
                        key.regionStableId = 0;
                        assignLatticeCoordinate(key, source.axis, 0);
                    } else {
                        key.kind =
                            PlanarPressureRegionFragmentFaceKind::
                                SameRegionGrid;
                        faceCoordinate = boundary->faceCoordinate;
                        assignLatticeCoordinate(
                            key, faceAxis, faceCoordinate);
                        key.facePeriodicImage = boundary->periodicImage;
                    }
                } else {
                    key.kind =
                        PlanarPressureRegionFragmentFaceKind::
                            SameRegionGrid;
                    const std::size_t cellCoordinate =
                        fragmentCellCoordinate(fragment, faceAxis);
                    const std::size_t count = axisCellCount(
                        result.cellCounts, faceAxis);
                    faceCoordinate = lowerSide
                        ? cellCoordinate
                        : (cellCoordinate + 1) % count;
                    assignLatticeCoordinate(
                        key, faceAxis, faceCoordinate);
                    key.profileCellPeriodicImage =
                        fragment.axisCellPeriodicImage;
                    key.profileLowerBoundary = boundaryIdentity(
                        fragment.lowerBoundary);
                    key.profileUpperBoundary = boundaryIdentity(
                        fragment.upperBoundary);
                }

                HalfFace half;
                half.fragmentIndex = fragmentIndex;
                half.isMinusSide = !lowerSide;
                half.areaSquareMeters = fragmentFaceArea(
                    fragment, faceAxis, source.axis,
                    result.spacingMeters);
                half.wrappedCentroidMeters = halfFaceCentroid(
                    grid, fragment, faceAxis, source.axis,
                    faceCoordinate, boundary);
                half.halfCenterDistanceMeters = faceAxis == source.axis
                    ? (lowerSide
                        ? fragment.unwrappedAxisCentroidMeters
                            - fragment.unwrappedLowerCoordinateMeters
                        : fragment.unwrappedUpperCoordinateMeters
                            - fragment.unwrappedAxisCentroidMeters)
                    : 0.5 * vectorCoordinate(
                        result.spacingMeters, faceAxis);
                if (!std::isfinite(half.areaSquareMeters)
                    || !(half.areaSquareMeters > 0.0)
                    || !std::isfinite(half.halfCenterDistanceMeters)
                    || !(half.halfCenterDistanceMeters > 0.0)) {
                    throw std::invalid_argument(
                        "planar pressure region fragment half-face is invalid");
                }
                halfFaces[key].push_back(half);
            }
        }
    }

    DisjointSet components(source.fragments.size());
    std::set<std::uint64_t> linkStableIds;
    result.links.reserve(expectedLinks);
    for (const auto& [key, halves] : halfFaces) {
        if (halves.size() != 2
            || halves[0].isMinusSide == halves[1].isMinusSide) {
            throw std::invalid_argument(
                "planar pressure region fragment half-faces do not pair");
        }
        const HalfFace& minusHalf = halves[0].isMinusSide
            ? halves[0] : halves[1];
        const HalfFace& plusHalf = halves[0].isMinusSide
            ? halves[1] : halves[0];
        if (minusHalf.fragmentIndex == plusHalf.fragmentIndex) {
            throw std::invalid_argument(
                "planar pressure region fragment face self-connects");
        }
        const auto& minus = source.fragments[minusHalf.fragmentIndex];
        const auto& plus = source.fragments[plusHalf.fragmentIndex];
        const double areaScale = std::max({
            minusHalf.areaSquareMeters, plusHalf.areaSquareMeters, 1.0});
        const double coordinateScale = std::max({
            std::abs(minusHalf.wrappedCentroidMeters.x),
            std::abs(minusHalf.wrappedCentroidMeters.y),
            std::abs(minusHalf.wrappedCentroidMeters.z),
            std::abs(plusHalf.wrappedCentroidMeters.x),
            std::abs(plusHalf.wrappedCentroidMeters.y),
            std::abs(plusHalf.wrappedCentroidMeters.z), 1.0});
        const double tolerance = 1024.0
            * std::numeric_limits<double>::epsilon();
        if (std::abs(minusHalf.areaSquareMeters
                     - plusHalf.areaSquareMeters)
                > tolerance * areaScale
            || maximumAbsoluteDifference(
                    minusHalf.wrappedCentroidMeters,
                    plusHalf.wrappedCentroidMeters)
                > tolerance * coordinateScale) {
            throw std::invalid_argument(
                "planar pressure region fragment paired geometry disagrees");
        }

        PlanarPressureRegionFragmentFaceLink link;
        link.linkIndex = result.links.size();
        link.kind = key.kind;
        link.axis = key.axis;
        link.i = key.i;
        link.j = key.j;
        link.k = key.k;
        link.facePeriodicImage = key.facePeriodicImage;
        link.surfaceStableId = key.surfaceStableId;
        link.minusFragmentIndex = minusHalf.fragmentIndex;
        link.plusFragmentIndex = plusHalf.fragmentIndex;
        link.minusFragmentStableId = minus.stableId;
        link.plusFragmentStableId = plus.stableId;
        link.minusRegionStableId = minus.regionStableId;
        link.plusRegionStableId = plus.regionStableId;
        link.areaSquareMeters = minusHalf.areaSquareMeters;
        link.wrappedCentroidMeters = plusHalf.wrappedCentroidMeters;
        link.centerDistanceMeters =
            minusHalf.halfCenterDistanceMeters
            + plusHalf.halfCenterDistanceMeters;
        link.pressureJumpPascals =
            plus.pressurePascals - minus.pressurePascals;
        link.unitNormalMinusToPlus = unitAxis(link.axis);
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            if (link.surfaceStableId != 0
                || link.minusRegionStableId != link.plusRegionStableId
                || minus.pressurePascals != plus.pressurePascals) {
                throw std::invalid_argument(
                    "planar pressure region grid face crosses a region");
            }
            link.sameRegionGeometryWeightMeters =
                link.areaSquareMeters / link.centerDistanceMeters;
            link.crossesPeriodicBoundary =
                latticeCoordinate(key, link.axis) == 0;
            ++result.sameRegionGridLinkCount;
            if (link.crossesPeriodicBoundary) {
                ++result.periodicGridLinkCount;
            }
            components.unite(
                link.minusFragmentIndex, link.plusFragmentIndex);
        } else {
            if (link.surfaceStableId == 0) {
                throw std::invalid_argument(
                    "planar pressure region wall surface is invalid");
            }
            const auto& plusBoundary = plus.lowerBoundary;
            if (plusBoundary.kind
                    != PlanarPressureRegionFragmentBoundaryKind::PressureLayer
                || plusBoundary.surfaceStableId != link.surfaceStableId) {
                throw std::invalid_argument(
                    "planar pressure region wall orientation is invalid");
            }
            assignLatticeCoordinate(
                link, source.axis,
                fragmentCellCoordinate(plus, source.axis));
            link.facePeriodicImage = periodicImage(
                plusBoundary.unwrappedCoordinateMeters,
                vectorCoordinate(result.lowerMeters, source.axis),
                vectorCoordinate(result.upperMeters, source.axis));
            ++result.pressureLayerWallLinkCount;
        }
        if (!std::isfinite(link.centerDistanceMeters)
            || !(link.centerDistanceMeters > 0.0)
            || !std::isfinite(link.sameRegionGeometryWeightMeters)
            || !std::isfinite(link.pressureJumpPascals)) {
            throw std::invalid_argument(
                "planar pressure region face-link metric is invalid");
        }
        link.stableId = linkStableId(link);
        if (!linkStableIds.insert(link.stableId).second) {
            throw std::invalid_argument(
                "planar pressure region face-link stable ID collides");
        }
        for (const std::size_t index : {
                 link.minusFragmentIndex, link.plusFragmentIndex}) {
            auto& summary = result.fragments[index];
            ++summary.incidentFaceCount;
            summary.incidentFaceAreaSquareMeters += link.areaSquareMeters;
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
                ++summary.sameRegionGridFaceCount;
            } else {
                ++summary.pressureLayerWallFaceCount;
            }
        }
        result.totalUniqueFaceAreaSquareMeters += link.areaSquareMeters;
        result.links.push_back(link);
    }
    if (result.links.size() != expectedLinks
        || result.sameRegionGridLinkCount
                + result.pressureLayerWallLinkCount
            != result.links.size()) {
        throw std::invalid_argument(
            "planar pressure region fragment topology count does not close");
    }

    std::map<std::size_t, ComponentWork> componentWork;
    for (std::size_t index = 0; index < source.fragments.size(); ++index) {
        const std::size_t root = components.find(index);
        const auto& fragment = source.fragments[index];
        auto [entry, inserted] = componentWork.try_emplace(root);
        auto& component = entry->second;
        if (inserted) {
            component.root = root;
            component.regionStableId = fragment.regionStableId;
        } else if (component.regionStableId != fragment.regionStableId) {
            throw std::invalid_argument(
                "planar pressure region topology component crosses a region");
        }
        component.fragmentStableIds.push_back(fragment.stableId);
        component.volumeCubicMeters += fragment.volumeCubicMeters;
    }
    for (const auto& link : result.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            ++componentWork.at(
                components.find(link.minusFragmentIndex))
                    .sameRegionGridLinkCount;
        }
    }
    std::vector<ComponentWork> orderedComponents;
    orderedComponents.reserve(componentWork.size());
    for (auto& [root, component] : componentWork) {
        static_cast<void>(root);
        component.stableId = componentStableId(
            component.regionStableId, component.fragmentStableIds);
        orderedComponents.push_back(std::move(component));
    }
    std::ranges::sort(
        orderedComponents,
        [](const ComponentWork& first, const ComponentWork& second) {
            if (first.regionStableId != second.regionStableId) {
                return first.regionStableId < second.regionStableId;
            }
            return first.stableId < second.stableId;
        });
    if (orderedComponents.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar pressure region fragment topology component limit exceeded");
    }
    std::map<std::size_t, std::size_t> componentIndexByRoot;
    result.components.reserve(orderedComponents.size());
    for (std::size_t index = 0; index < orderedComponents.size(); ++index) {
        const auto& sourceComponent = orderedComponents[index];
        componentIndexByRoot.emplace(sourceComponent.root, index);
        result.components.push_back({
            index,
            sourceComponent.stableId,
            sourceComponent.regionStableId,
            sourceComponent.fragmentStableIds.size(),
            sourceComponent.sameRegionGridLinkCount,
            sourceComponent.volumeCubicMeters,
        });
    }
    for (std::size_t index = 0; index < result.fragments.size(); ++index) {
        result.fragments[index].componentIndex = componentIndexByRoot.at(
            components.find(index));
    }
    for (auto& link : result.links) {
        link.minusComponentIndex = componentIndexByRoot.at(
            components.find(link.minusFragmentIndex));
        link.plusComponentIndex = componentIndexByRoot.at(
            components.find(link.plusFragmentIndex));
        if (link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid
            && link.minusComponentIndex != link.plusComponentIndex) {
            throw std::invalid_argument(
                "planar pressure region grid link leaves its component");
        }
    }

    const double areaScale = std::max(
        result.totalUniqueFaceAreaSquareMeters, 1.0);
    const double areaTolerance = 4096.0
        * std::numeric_limits<double>::epsilon() * areaScale;
    for (auto& fragment : result.fragments) {
        fragment.boundaryAreaClosureResidualSquareMeters =
            fragment.incidentFaceAreaSquareMeters
            - fragment.expectedBoundaryAreaSquareMeters;
        result.totalIncidentFaceAreaSquareMeters +=
            fragment.incidentFaceAreaSquareMeters;
        result.totalExpectedFragmentBoundaryAreaSquareMeters +=
            fragment.expectedBoundaryAreaSquareMeters;
        result.maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters =
            std::max(
                result
                    .maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters,
                std::abs(
                    fragment.boundaryAreaClosureResidualSquareMeters));
        if (fragment.incidentFaceCount != 6
            || fragment.sameRegionGridFaceCount
                    + fragment.pressureLayerWallFaceCount
                != 6) {
            throw std::invalid_argument(
                "planar pressure region fragment face incidence does not close");
        }
    }
    if (std::abs(result.totalIncidentFaceAreaSquareMeters
                 - 2.0 * result.totalUniqueFaceAreaSquareMeters)
            > areaTolerance
        || std::abs(result.totalExpectedFragmentBoundaryAreaSquareMeters
                    - result.totalIncidentFaceAreaSquareMeters)
            > areaTolerance
        || result
                .maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters
            > areaTolerance) {
        throw std::invalid_argument(
            "planar pressure region fragment face area does not close");
    }
    result.ownedStorageBytes = storageBytes(
        result.links.size(), result.fragments.size(),
        result.components.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region fragment topology byte limit exceeded");
    }
    result.fingerprint = topologyFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentTopology
buildPlanarPressureRegionFragmentTopology(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopologyLimits& limits) {
    auto result = buildTopology(grid, sweep, fragments, limits);
    validatePlanarPressureRegionFragmentTopology(
        result, grid, sweep, fragments, limits);
    return result;
}

void validatePlanarPressureRegionFragmentTopology(
    const PlanarPressureRegionFragmentTopology& topology,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopologyLimits& limits) {
    validateLimits(limits);
    if (topology.links.size() > limits.maximumLinks
        || topology.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar pressure region fragment topology validation limit exceeded");
    }
    const auto expected = buildTopology(grid, sweep, fragments, limits);
    if (expected != topology) {
        throw std::invalid_argument(
            "planar pressure region fragment topology is invalid");
    }
}

} // namespace simwing::fsi::fluid
