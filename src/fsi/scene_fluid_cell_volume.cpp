#include "scene_fluid_cell_volume.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
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

struct Contribution {
    std::size_t cellIndex = 0;
    std::size_t regionIndex = 0;
    std::size_t ordinal = 0;
    double volumeCubicMeters = 0.0;
    Vec3 firstMomentMeters4;
};

struct PolyhedronMeasure {
    double volumeCubicMeters = 0.0;
    Vec3 firstMomentMeters4;
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

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vec3 add(const Vec3& first, const Vec3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

Vec3 scale(const Vec3& value, const double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

bool finite(const Vec3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

double length(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

double tolerance(const SceneFluidCellVolumeSettings& settings,
                 const double referenceVolume) {
    return std::max(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance * std::abs(referenceVolume));
}

void validateSettings(const SceneFluidCellVolumeSettings& settings) {
    if (!std::isfinite(settings.absoluteVolumeToleranceCubicMeters)
        || settings.absoluteVolumeToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeVolumeTolerance)
        || settings.relativeVolumeTolerance < 0.0
        || (settings.absoluteVolumeToleranceCubicMeters == 0.0
            && settings.relativeVolumeTolerance == 0.0)) {
        throw std::invalid_argument(
            "scene fluid cell-volume tolerances are invalid");
    }
}

std::size_t outsideRegionIndex(
    const SceneFluidSurfaceDefinition& surface) {
    std::size_t result = surface.regions.size();
    for (std::size_t index = 0; index < surface.regions.size(); ++index) {
        if (surface.regions[index].kind != RegionKind::Outside) {
            continue;
        }
        if (result != surface.regions.size()) {
            throw std::invalid_argument(
                "scene fluid cell volumes require exactly one Outside region");
        }
        result = index;
    }
    if (result == surface.regions.size()) {
        throw std::invalid_argument(
            "scene fluid cell volumes require one Outside region");
    }
    return result;
}

void appendContribution(std::vector<Contribution>& contributions,
                        const SceneFluidCellVolumeLimits& limits,
                        const std::size_t cellIndex,
                        const std::size_t region,
                        const double value,
                        const Vec3& firstMoment) {
    if (!std::isfinite(value) || !finite(firstMoment)) {
        throw std::invalid_argument(
            "scene fluid cell-volume contribution is non-finite");
    }
    if (contributions.size() == limits.maximumContributionEvents) {
        throw std::length_error(
            "scene fluid cell-volume contributions exceed their limit");
    }
    contributions.push_back({cellIndex, region, contributions.size(), value,
                             firstMoment});
}

using Polygon = std::vector<Vec3>;

struct ConvexPolyhedron {
    std::vector<Polygon> faces;
};

double coordinate(const Vec3& point, const std::size_t axis) {
    if (axis == 0) return point.x;
    if (axis == 1) return point.y;
    return point.z;
}

void setCoordinate(Vec3& point,
                   const std::size_t axis,
                   const double value) {
    if (axis == 0) point.x = value;
    else if (axis == 1) point.y = value;
    else point.z = value;
}

bool pointsNear(const Vec3& first, const Vec3& second) {
    const double scale = std::max({
        1.0,
        std::abs(first.x), std::abs(first.y), std::abs(first.z),
        std::abs(second.x), std::abs(second.y), std::abs(second.z),
    });
    const double bound = 256.0
        * std::numeric_limits<double>::epsilon() * scale;
    return length(subtract(first, second)) <= bound;
}

Vec3 planeIntersection(const Vec3& first,
                       const Vec3& second,
                       const std::size_t axis,
                       const double plane) {
    const double firstCoordinate = coordinate(first, axis);
    const double denominator = coordinate(second, axis) - firstCoordinate;
    if (!std::isfinite(denominator) || denominator == 0.0) {
        throw std::invalid_argument(
            "scene fluid tetrahedron clipping edge is parallel to its crossing plane");
    }
    const double parameter = (plane - firstCoordinate) / denominator;
    Vec3 result{
        first.x + parameter * (second.x - first.x),
        first.y + parameter * (second.y - first.y),
        first.z + parameter * (second.z - first.z),
    };
    setCoordinate(result, axis, plane);
    if (!std::isfinite(result.x) || !std::isfinite(result.y)
        || !std::isfinite(result.z)) {
        throw std::invalid_argument(
            "scene fluid tetrahedron clipping intersection is non-finite");
    }
    return result;
}

void appendDistinct(Polygon& polygon, const Vec3& point) {
    if (polygon.empty() || !pointsNear(polygon.back(), point)) {
        polygon.push_back(point);
    }
}

Vec3 planeNormal(const std::size_t axis, const bool keepGreater) {
    Vec3 result;
    setCoordinate(result, axis, keepGreater ? -1.0 : 1.0);
    return result;
}

Vec3 faceBasisU(const std::size_t axis) {
    if (axis == 0) return {0.0, 1.0, 0.0};
    if (axis == 1) return {0.0, 0.0, 1.0};
    return {1.0, 0.0, 0.0};
}

void clipPolyhedron(ConvexPolyhedron& polyhedron,
                    const std::size_t axis,
                    const double plane,
                    const bool keepGreater) {
    std::vector<Polygon> clippedFaces;
    clippedFaces.reserve(polyhedron.faces.size() + 1);
    Polygon capPoints;
    const auto inside = [&](const Vec3& point) {
        return keepGreater
            ? coordinate(point, axis) >= plane
            : coordinate(point, axis) <= plane;
    };
    for (const auto& face : polyhedron.faces) {
        Polygon clipped;
        clipped.reserve(face.size() + 1);
        for (std::size_t index = 0; index < face.size(); ++index) {
            const Vec3& first = face[index];
            const Vec3& second = face[(index + 1) % face.size()];
            const bool firstInside = inside(first);
            const bool secondInside = inside(second);
            if (firstInside) appendDistinct(clipped, first);
            if (firstInside != secondInside) {
                const Vec3 intersection = planeIntersection(
                    first, second, axis, plane);
                appendDistinct(clipped, intersection);
                bool duplicate = false;
                for (const auto& existing : capPoints) {
                    duplicate = duplicate || pointsNear(existing, intersection);
                }
                if (!duplicate) capPoints.push_back(intersection);
            }
        }
        if (clipped.size() > 1
            && pointsNear(clipped.front(), clipped.back())) {
            clipped.pop_back();
        }
        if (clipped.size() >= 3) clippedFaces.push_back(std::move(clipped));
    }
    if (capPoints.size() >= 3) {
        Vec3 centre;
        for (const auto& point : capPoints) {
            centre.x += point.x;
            centre.y += point.y;
            centre.z += point.z;
        }
        const double inverseCount = 1.0
            / static_cast<double>(capPoints.size());
        centre.x *= inverseCount;
        centre.y *= inverseCount;
        centre.z *= inverseCount;
        const Vec3 normal = planeNormal(axis, keepGreater);
        const Vec3 basisU = faceBasisU(axis);
        const Vec3 basisV = cross(normal, basisU);
        std::ranges::sort(
            capPoints,
            [&](const Vec3& first, const Vec3& second) {
                const Vec3 firstRelative = subtract(first, centre);
                const Vec3 secondRelative = subtract(second, centre);
                const double firstAngle = std::atan2(
                    dot(firstRelative, basisV), dot(firstRelative, basisU));
                const double secondAngle = std::atan2(
                    dot(secondRelative, basisV), dot(secondRelative, basisU));
                if (firstAngle != secondAngle) return firstAngle < secondAngle;
                return std::tie(first.x, first.y, first.z)
                    < std::tie(second.x, second.y, second.z);
            });
        clippedFaces.push_back(std::move(capPoints));
    }
    polyhedron.faces = std::move(clippedFaces);
}

ConvexPolyhedron tetrahedron(const Vec3& origin,
                             const Vec3& first,
                             Vec3 second,
                             Vec3 third,
                             const double signedSixVolume) {
    if (signedSixVolume < 0.0) std::swap(second, third);
    const std::array<Vec3, 4> vertices{{origin, first, second, third}};
    ConvexPolyhedron result;
    result.faces = {
        {vertices[0], vertices[2], vertices[1]},
        {vertices[0], vertices[1], vertices[3]},
        {vertices[0], vertices[3], vertices[2]},
        {vertices[1], vertices[2], vertices[3]},
    };
    return result;
}

PolyhedronMeasure polyhedronMeasure(const ConvexPolyhedron& polyhedron) {
    Vec3 centre;
    std::size_t pointCount = 0;
    for (const auto& face : polyhedron.faces) {
        for (const auto& point : face) {
            centre.x += point.x;
            centre.y += point.y;
            centre.z += point.z;
            ++pointCount;
        }
    }
    if (pointCount == 0) return {};
    const double inverseCount = 1.0 / static_cast<double>(pointCount);
    centre.x *= inverseCount;
    centre.y *= inverseCount;
    centre.z *= inverseCount;
    PolyhedronMeasure result;
    double summedSixVolume = 0.0;
    for (const auto& face : polyhedron.faces) {
        for (std::size_t vertex = 1; vertex + 1 < face.size(); ++vertex) {
            const double sixVolume = dot(
                subtract(face[0], centre),
                cross(subtract(face[vertex], centre),
                      subtract(face[vertex + 1], centre)));
            const double volume = sixVolume / 6.0;
            const Vec3 tetrahedronCentroid = scale(
                add(add(centre, face[0]),
                    add(face[vertex], face[vertex + 1])),
                0.25);
            summedSixVolume += sixVolume;
            result.firstMomentMeters4 = add(
                result.firstMomentMeters4,
                scale(tetrahedronCentroid, volume));
        }
    }
    result.volumeCubicMeters = summedSixVolume / 6.0;
    if (!std::isfinite(result.volumeCubicMeters)
        || !finite(result.firstMomentMeters4)) {
        throw std::invalid_argument(
            "scene fluid clipped tetrahedron measure is non-finite");
    }
    return result;
}

std::pair<std::size_t, std::size_t> cellRange(
    const double minimum,
    const double maximum,
    const double gridLower,
    const double spacing,
    const std::size_t count) {
    const auto index = [&](const double coordinateValue) {
        const double scaled = (coordinateValue - gridLower) / spacing;
        if (!std::isfinite(scaled)) {
            throw std::invalid_argument(
                "scene fluid tetrahedron cell range is non-finite");
        }
        if (scaled <= 0.0) return std::size_t{0};
        if (scaled >= static_cast<double>(count)) return count - 1;
        return static_cast<std::size_t>(std::floor(scaled));
    };
    return {index(minimum), index(maximum)};
}

std::vector<double> wholeSurfaceVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const std::size_t outsideIndex,
    const fluid::PeriodicCartesianGrid& grid) {
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    const double domainVolume = (upper.x - lower.x)
        * (upper.y - lower.y) * (upper.z - lower.z);
    std::vector<double> result(surface.regions.size(), 0.0);
    result[outsideIndex] = domainVolume;
    const auto addTriangle = [&](
        const std::array<std::size_t, 3>& vertices,
        const std::size_t negativeRegion,
        const std::size_t positiveRegion) {
        if (negativeRegion == positiveRegion) return;
        const auto& first =
            state.vertices[vertices[0]].positionMeters;
        const auto& second =
            state.vertices[vertices[1]].positionMeters;
        const auto& third =
            state.vertices[vertices[2]].positionMeters;
        const Vec3 areaVector = cross(
            subtract(second, first), subtract(third, first));
        const Vec3 centroid{
            (first.x + second.x + third.x) / 3.0,
            (first.y + second.y + third.y) / 3.0,
            (first.z + second.z + third.z) / 3.0,
        };
        const double contribution = dot(centroid, areaVector) / 6.0;
        result[negativeRegion] += contribution;
        result[positiveRegion] -= contribution;
    };
    for (const auto& triangle : surface.triangles) {
        addTriangle(
            triangle.vertexIndices,
            triangle.negativeSideRegionIndex,
            triangle.positiveSideRegionIndex);
    }
    for (const auto& triangle : caps.triangles) {
        const auto& cap = caps.caps[triangle.openingIndex];
        addTriangle(
            triangle.vertexIndices,
            cap.negativeSideRegionIndex,
            cap.positiveSideRegionIndex);
    }
    return result;
}

std::uint64_t volumeFingerprint(const SceneFluidCellVolumeSet& volumes) {
    Fingerprint fingerprint;
    fingerprint.integer(volumes.version);
    fingerprint.integer(volumes.surfaceDefinitionFingerprint);
    fingerprint.integer(volumes.surfaceStateFingerprint);
    fingerprint.integer(volumes.gridEpochFingerprint);
    fingerprint.integer(volumes.openingCapFingerprint);
    fingerprint.integer(volumes.structureDefinitionFingerprint);
    fingerprint.integer(volumes.acceptedStepCount);
    fingerprint.real(volumes.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.z));
    for (const double value : {
             volumes.lowerMeters.x, volumes.lowerMeters.y,
             volumes.lowerMeters.z, volumes.upperMeters.x,
             volumes.upperMeters.y, volumes.upperMeters.z,
             volumes.settings.absoluteVolumeToleranceCubicMeters,
             volumes.settings.relativeVolumeTolerance,
             volumes.settings.openingCaps.planarityToleranceMeters,
             volumes.settings.openingCaps.minimumTriangleAreaSquareMeters,
             volumes.settings.openingCaps.convexityTolerance,
             volumes.cellVolumeCubicMeters,
             volumes.openingCapAreaSquareMeters,
             volumes.maximumTetrahedronVolumeResidualCubicMeters,
             volumes.maximumCellVolumeResidualCubicMeters,
             volumes.maximumCellFirstMomentResidualMeters4,
             volumes.maximumRegionVolumeResidualCubicMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(volumes.outsideRegionId);
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.openingCapCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.tetrahedronCellClipCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.nonzeroTetrahedronCellClipCount));
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cells.size()));
    for (const auto& cell : volumes.cells) {
        fingerprint.integer(static_cast<std::uint64_t>(cell.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.firstRegionVolume));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.regionVolumeCount));
        fingerprint.real(cell.assignedVolumeCubicMeters);
        fingerprint.real(cell.volumeResidualCubicMeters);
        for (const double value : {
                 cell.assignedFirstMomentMeters4.x,
                 cell.assignedFirstMomentMeters4.y,
                 cell.assignedFirstMomentMeters4.z,
                 cell.firstMomentResidualMeters4.x,
                 cell.firstMomentResidualMeters4.y,
                 cell.firstMomentResidualMeters4.z}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.cellRegionVolumes.size()));
    for (const auto& value : volumes.cellRegionVolumes) {
        fingerprint.integer(value.regionId);
        fingerprint.real(value.volumeCubicMeters);
        fingerprint.real(value.volumeFraction);
        for (const double component : {
                 value.firstMomentMeters4.x,
                 value.firstMomentMeters4.y,
                 value.firstMomentMeters4.z,
                 value.centroidMeters.x,
                 value.centroidMeters.y,
                 value.centroidMeters.z}) {
            fingerprint.real(component);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.regionVolumes.size()));
    for (const auto& value : volumes.regionVolumes) {
        fingerprint.integer(value.regionId);
        fingerprint.real(value.summedCellVolumeCubicMeters);
        fingerprint.real(value.wholeSurfaceVolumeCubicMeters);
        fingerprint.real(value.volumeResidualCubicMeters);
    }
    return fingerprint.value();
}

SceneFluidCellVolumeSet buildVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidCellVolumeSettings& settings,
    const SceneFluidCellVolumeLimits& limits) {
    validateSettings(settings);
    validateSceneFluidOpeningCaps(caps, surface, state);
    const std::size_t outsideIndex = outsideRegionIndex(surface);
    if (grid.cellCount() > limits.maximumCells
        || grid.cellCount() > limits.maximumCellRegionVolumes) {
        throw std::length_error(
            "scene fluid cell volumes exceed their cell limits");
    }

    SceneFluidCellVolumeSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.gridEpochFingerprint = epoch.fingerprint;
    result.openingCapFingerprint = caps.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    result.outsideRegionId = surface.regions[outsideIndex].id;
    result.cellVolumeCubicMeters = grid.cellVolumeCubicMeters();
    result.openingCapCount = caps.caps.size();
    result.openingCapAreaSquareMeters = caps.totalAreaSquareMeters;

    const auto counts = grid.cellCounts();
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    std::vector<Contribution> contributions;
    contributions.reserve(std::min(
        limits.maximumContributionEvents, grid.cellCount()));
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        const std::size_t i = cell % counts.x;
        const std::size_t j = (cell / counts.x) % counts.y;
        const std::size_t k = (cell / counts.x) / counts.y;
        const Vec3 centre{
            lower.x + (static_cast<double>(i) + 0.5) * spacing.x,
            lower.y + (static_cast<double>(j) + 0.5) * spacing.y,
            lower.z + (static_cast<double>(k) + 0.5) * spacing.z,
        };
        appendContribution(
            contributions, limits, cell, outsideIndex,
            result.cellVolumeCubicMeters,
            scale(centre, result.cellVolumeCubicMeters));
    }

    const Vec3 origin{lower.x, lower.y, lower.z};
    const double clipVolumeTolerance = tolerance(
        settings, result.cellVolumeCubicMeters);
    // A closed oriented surface is the signed sum of tetrahedra from any fixed
    // reference point. Clipping each tetrahedron distributes that same chain
    // cell-by-cell without requiring face-local contours to close in one tile.
    const auto distributeTriangle = [&](
        const std::array<std::size_t, 3>& vertices,
        const std::size_t negativeRegion,
        const std::size_t positiveRegion) {
        if (negativeRegion == positiveRegion) return;
        const auto& first =
            state.vertices[vertices[0]].positionMeters;
        const auto& second =
            state.vertices[vertices[1]].positionMeters;
        const auto& third =
            state.vertices[vertices[2]].positionMeters;
        const double signedSixVolume = dot(
            subtract(first, origin),
            cross(subtract(second, origin), subtract(third, origin)));
        if (!std::isfinite(signedSixVolume)) {
            throw std::invalid_argument(
                "scene fluid signed tetrahedron volume is non-finite");
        }
        if (signedSixVolume == 0.0) return;
        const double minimumX = std::min({origin.x, first.x, second.x, third.x});
        const double maximumX = std::max({origin.x, first.x, second.x, third.x});
        const double minimumY = std::min({origin.y, first.y, second.y, third.y});
        const double maximumY = std::max({origin.y, first.y, second.y, third.y});
        const double minimumZ = std::min({origin.z, first.z, second.z, third.z});
        const double maximumZ = std::max({origin.z, first.z, second.z, third.z});
        const auto [firstI, lastI] = cellRange(
            minimumX, maximumX, lower.x, spacing.x, counts.x);
        const auto [firstJ, lastJ] = cellRange(
            minimumY, maximumY, lower.y, spacing.y, counts.y);
        const auto [firstK, lastK] = cellRange(
            minimumZ, maximumZ, lower.z, spacing.z, counts.z);
        const ConvexPolyhedron source = tetrahedron(
            origin, first, second, third, signedSixVolume);
        double summedClippedVolume = 0.0;
        for (std::size_t k = firstK; k <= lastK; ++k) {
            for (std::size_t j = firstJ; j <= lastJ; ++j) {
                for (std::size_t i = firstI; i <= lastI; ++i) {
                    if (result.tetrahedronCellClipCount
                        == limits.maximumTetrahedronCellClips) {
                        throw std::length_error(
                            "scene fluid tetrahedron/cell clips exceed their limit");
                    }
                    ++result.tetrahedronCellClipCount;
                    ConvexPolyhedron clipped = source;
                    const double cellLowerX = lower.x
                        + static_cast<double>(i) * spacing.x;
                    const double cellLowerY = lower.y
                        + static_cast<double>(j) * spacing.y;
                    const double cellLowerZ = lower.z
                        + static_cast<double>(k) * spacing.z;
                    clipPolyhedron(clipped, 0, cellLowerX, true);
                    clipPolyhedron(clipped, 0, cellLowerX + spacing.x, false);
                    clipPolyhedron(clipped, 1, cellLowerY, true);
                    clipPolyhedron(clipped, 1, cellLowerY + spacing.y, false);
                    clipPolyhedron(clipped, 2, cellLowerZ, true);
                    clipPolyhedron(clipped, 2, cellLowerZ + spacing.z, false);
                    const PolyhedronMeasure clippedMeasure =
                        polyhedronMeasure(clipped);
                    double clippedVolume =
                        clippedMeasure.volumeCubicMeters;
                    if (clippedVolume < -clipVolumeTolerance
                        || clippedVolume
                            > result.cellVolumeCubicMeters
                                + clipVolumeTolerance) {
                        throw std::invalid_argument(
                            "scene fluid clipped tetrahedron volume is outside its cell");
                    }
                    if (clippedVolume <= 0.0) continue;
                    clippedVolume = std::min(
                        clippedVolume, result.cellVolumeCubicMeters);
                    summedClippedVolume += clippedVolume;
                    ++result.nonzeroTetrahedronCellClipCount;
                    const double signedVolume = std::copysign(
                        clippedVolume, signedSixVolume);
                    const double signedMeasureScale = std::copysign(
                        clippedVolume
                            / clippedMeasure.volumeCubicMeters,
                        signedSixVolume);
                    const Vec3 signedFirstMoment = scale(
                        clippedMeasure.firstMomentMeters4,
                        signedMeasureScale);
                    const std::size_t cell = grid.cellIndex(i, j, k);
                    appendContribution(
                        contributions, limits, cell,
                        negativeRegion, signedVolume,
                        signedFirstMoment);
                    appendContribution(
                        contributions, limits, cell,
                        positiveRegion, -signedVolume,
                        scale(signedFirstMoment, -1.0));
                }
            }
        }
        const double tetrahedronVolume = std::abs(signedSixVolume) / 6.0;
        const double tetrahedronResidual =
            summedClippedVolume - tetrahedronVolume;
        result.maximumTetrahedronVolumeResidualCubicMeters = std::max(
            result.maximumTetrahedronVolumeResidualCubicMeters,
            std::abs(tetrahedronResidual));
        if (!std::isfinite(tetrahedronResidual)
            || std::abs(tetrahedronResidual)
                > tolerance(settings, tetrahedronVolume)) {
            throw std::invalid_argument(
                "scene fluid clipped cells do not close their signed tetrahedron");
        }
    };
    for (const auto& triangle : surface.triangles) {
        distributeTriangle(
            triangle.vertexIndices,
            triangle.negativeSideRegionIndex,
            triangle.positiveSideRegionIndex);
    }
    for (const auto& triangle : caps.triangles) {
        const auto& cap = caps.caps[triangle.openingIndex];
        distributeTriangle(
            triangle.vertexIndices,
            cap.negativeSideRegionIndex,
            cap.positiveSideRegionIndex);
    }

    std::ranges::sort(
        contributions,
        [](const Contribution& first, const Contribution& second) {
            return std::tie(
                       first.cellIndex, first.regionIndex, first.ordinal)
                < std::tie(
                    second.cellIndex, second.regionIndex, second.ordinal);
        });
    const double cellTolerance = tolerance(
        settings, result.cellVolumeCubicMeters);
    std::vector<double> summedRegions(surface.regions.size(), 0.0);
    result.cells.reserve(grid.cellCount());
    std::size_t contributionIndex = 0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t cellIndex = grid.cellIndex(i, j, k);
                SceneFluidCellVolume cell;
                cell.cellIndex = cellIndex;
                cell.cell = {i, j, k};
                cell.firstRegionVolume = result.cellRegionVolumes.size();
                while (contributionIndex < contributions.size()
                       && contributions[contributionIndex].cellIndex
                           == cellIndex) {
                    const std::size_t region =
                        contributions[contributionIndex].regionIndex;
                    double volume = 0.0;
                    Vec3 firstMoment;
                    while (contributionIndex < contributions.size()
                           && contributions[contributionIndex].cellIndex
                               == cellIndex
                           && contributions[contributionIndex].regionIndex
                               == region) {
                        volume += contributions[contributionIndex]
                                      .volumeCubicMeters;
                        firstMoment = add(
                            firstMoment,
                            contributions[contributionIndex]
                                .firstMomentMeters4);
                        ++contributionIndex;
                    }
                    if (!std::isfinite(volume) || volume < -cellTolerance) {
                        throw std::invalid_argument(
                            "scene fluid cell region has a negative or non-finite volume");
                    }
                    if (std::abs(volume) <= cellTolerance) {
                        volume = 0.0;
                    }
                    if (volume == 0.0) {
                        continue;
                    }
                    if (!finite(firstMoment)) {
                        throw std::invalid_argument(
                            "scene fluid cell-region first moment is non-finite");
                    }
                    Vec3 centroid = scale(firstMoment, 1.0 / volume);
                    const Vec3 cellLower{
                        lower.x + static_cast<double>(i) * spacing.x,
                        lower.y + static_cast<double>(j) * spacing.y,
                        lower.z + static_cast<double>(k) * spacing.z,
                    };
                    const Vec3 cellUpper{
                        cellLower.x + spacing.x,
                        cellLower.y + spacing.y,
                        cellLower.z + spacing.z,
                    };
                    const double coordinateScale = std::max({
                        1.0,
                        std::abs(cellLower.x), std::abs(cellLower.y),
                        std::abs(cellLower.z), std::abs(cellUpper.x),
                        std::abs(cellUpper.y), std::abs(cellUpper.z),
                    });
                    const double coordinateTolerance = 4096.0
                        * std::numeric_limits<double>::epsilon()
                        * coordinateScale;
                    if (!finite(centroid)
                        || centroid.x < cellLower.x - coordinateTolerance
                        || centroid.x > cellUpper.x + coordinateTolerance
                        || centroid.y < cellLower.y - coordinateTolerance
                        || centroid.y > cellUpper.y + coordinateTolerance
                        || centroid.z < cellLower.z - coordinateTolerance
                        || centroid.z > cellUpper.z + coordinateTolerance) {
                        throw std::invalid_argument(
                            "scene fluid cell-region centroid is outside its cell");
                    }
                    centroid.x = std::clamp(
                        centroid.x, cellLower.x, cellUpper.x);
                    centroid.y = std::clamp(
                        centroid.y, cellLower.y, cellUpper.y);
                    centroid.z = std::clamp(
                        centroid.z, cellLower.z, cellUpper.z);
                    if (result.cellRegionVolumes.size()
                        == limits.maximumCellRegionVolumes) {
                        throw std::length_error(
                            "scene fluid cell-region volumes exceed their limit");
                    }
                    result.cellRegionVolumes.push_back({
                        surface.regions[region].id,
                        volume,
                        volume / result.cellVolumeCubicMeters,
                        firstMoment,
                        centroid,
                    });
                    cell.assignedVolumeCubicMeters += volume;
                    cell.assignedFirstMomentMeters4 = add(
                        cell.assignedFirstMomentMeters4, firstMoment);
                    summedRegions[region] += volume;
                }
                cell.regionVolumeCount = result.cellRegionVolumes.size()
                    - cell.firstRegionVolume;
                cell.volumeResidualCubicMeters =
                    cell.assignedVolumeCubicMeters
                    - result.cellVolumeCubicMeters;
                const Vec3 cellCentroid{
                    lower.x + (static_cast<double>(i) + 0.5) * spacing.x,
                    lower.y + (static_cast<double>(j) + 0.5) * spacing.y,
                    lower.z + (static_cast<double>(k) + 0.5) * spacing.z,
                };
                cell.firstMomentResidualMeters4 = subtract(
                    cell.assignedFirstMomentMeters4,
                    scale(cellCentroid, result.cellVolumeCubicMeters));
                result.maximumCellVolumeResidualCubicMeters = std::max(
                    result.maximumCellVolumeResidualCubicMeters,
                    std::abs(cell.volumeResidualCubicMeters));
                result.maximumCellFirstMomentResidualMeters4 = std::max(
                    result.maximumCellFirstMomentResidualMeters4,
                    length(cell.firstMomentResidualMeters4));
                const double cellCoordinateScale = std::max({
                    1.0, std::abs(cellCentroid.x),
                    std::abs(cellCentroid.y),
                    std::abs(cellCentroid.z),
                });
                const double cellMomentTolerance = 16.0
                    * cellTolerance * cellCoordinateScale;
                if (!std::isfinite(cell.volumeResidualCubicMeters)
                    || std::abs(cell.volumeResidualCubicMeters)
                        > cellTolerance
                    || !finite(cell.firstMomentResidualMeters4)
                    || length(cell.firstMomentResidualMeters4)
                        > cellMomentTolerance) {
                    throw std::invalid_argument(
                        "scene fluid cell-region measures do not close the cell");
                }
                result.cells.push_back(cell);
            }
        }
    }
    if (contributionIndex != contributions.size()) {
        throw std::logic_error(
            "scene fluid cell-volume contribution escaped the grid");
    }

    const auto whole = wholeSurfaceVolumes(
        surface, state, caps, outsideIndex, grid);
    const auto upper = grid.upperMeters();
    const double domainVolume = (upper.x - lower.x)
        * (upper.y - lower.y) * (upper.z - lower.z);
    const double regionTolerance = tolerance(settings, domainVolume);
    result.regionVolumes.reserve(surface.regions.size());
    for (std::size_t region = 0; region < surface.regions.size(); ++region) {
        double expected = whole[region];
        if (!std::isfinite(expected) || expected < -regionTolerance) {
            throw std::invalid_argument(
                "scene fluid whole-surface region volume is invalid");
        }
        if (std::abs(expected) <= regionTolerance) {
            expected = 0.0;
        }
        if (region != outsideIndex && !(expected > 0.0)) {
            throw std::invalid_argument(
                "scene fluid non-Outside region has no positive closed volume");
        }
        const double residual = summedRegions[region] - expected;
        result.maximumRegionVolumeResidualCubicMeters = std::max(
            result.maximumRegionVolumeResidualCubicMeters,
            std::abs(residual));
        if (!std::isfinite(residual)
            || std::abs(residual) > regionTolerance) {
            throw std::invalid_argument(
                "scene fluid cell volumes do not close their whole-surface region volume");
        }
        result.regionVolumes.push_back({
            surface.regions[region].id,
            summedRegions[region],
            expected,
            residual,
        });
    }

    std::size_t cellBytes = 0;
    std::size_t cellRegionBytes = 0;
    std::size_t regionBytes = 0;
    std::size_t first = 0;
    std::size_t total = 0;
    if (!checkedMultiply(result.cells.size(),
                         sizeof(SceneFluidCellVolume), cellBytes)
        || !checkedMultiply(result.cellRegionVolumes.size(),
                            sizeof(SceneFluidCellRegionVolume),
                            cellRegionBytes)
        || !checkedMultiply(result.regionVolumes.size(),
                            sizeof(SceneFluidRegionVolume), regionBytes)
        || !checkedAdd(cellBytes, cellRegionBytes, first)
        || !checkedAdd(first, regionBytes, total)
        || total > limits.maximumVolumeBytes) {
        throw std::length_error(
            "scene fluid cell-volume result exceeds its byte limit");
    }
    result.fingerprint = volumeFingerprint(result);
    return result;
}

} // namespace

SceneFluidCellVolumeSet buildSceneFluidCellVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidCellVolumeSettings& settings,
    const SceneFluidCellVolumeLimits& limits) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    const auto caps = buildSceneFluidOpeningCaps(
        surface, state, settings.openingCaps, limits.openingCaps);
    auto result = buildVolumes(
        surface, state, grid, epoch, caps, settings, limits);
    validateSceneFluidCellVolumes(
        result, surface, state, grid, transfer, epoch);
    return result;
}

void validateSceneFluidCellVolumeIntegrity(
    const SceneFluidCellVolumeSet& volumes) {
    if (volumes.version != sceneFluidCellVolumeVersion
        || volumes.fingerprint == 0
        || volumes.fingerprint != volumeFingerprint(volumes)) {
        throw std::invalid_argument(
            "scene fluid cell-volume integrity is invalid");
    }
}

void validateSceneFluidCellVolumes(
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidCellVolumeIntegrity(volumes);
    if (volumes.version != sceneFluidCellVolumeVersion
        || volumes.surfaceDefinitionFingerprint != surface.fingerprint
        || volumes.surfaceStateFingerprint != state.fingerprint
        || volumes.gridEpochFingerprint != epoch.fingerprint
        || volumes.openingCapFingerprint == 0
        || volumes.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || volumes.acceptedStepCount != state.acceptedStepCount
        || volumes.simulationTimeSeconds != state.simulationTimeSeconds
        || volumes.cellCounts != grid.cellCounts()
        || volumes.lowerMeters != grid.lowerMeters()
        || volumes.upperMeters != grid.upperMeters()) {
        throw std::invalid_argument(
            "scene fluid cell-volume identity is invalid");
    }
    const SceneFluidOpeningCapLimits unlimitedCapLimits{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto caps = buildSceneFluidOpeningCaps(
        surface, state, volumes.settings.openingCaps, unlimitedCapLimits);
    SceneFluidCellVolumeLimits unlimited;
    unlimited.openingCaps = unlimitedCapLimits;
    unlimited.maximumCells = std::numeric_limits<std::size_t>::max();
    unlimited.maximumContributionEvents =
        std::numeric_limits<std::size_t>::max();
    unlimited.maximumTetrahedronCellClips =
        std::numeric_limits<std::size_t>::max();
    unlimited.maximumCellRegionVolumes =
        std::numeric_limits<std::size_t>::max();
    unlimited.maximumVolumeBytes = std::numeric_limits<std::size_t>::max();
    const auto expected = buildVolumes(
        surface, state, grid, epoch, caps, volumes.settings, unlimited);
    if (volumes != expected) {
        throw std::invalid_argument(
            "scene fluid cell-volume payload is invalid");
    }
}

} // namespace simwing::fsi
