#include "fluid/scene_surface_clipping.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
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

    template<typename Enum>
    void enumeration(const Enum value) {
        using Underlying = std::underlying_type_t<Enum>;
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

Vec3 add(const Vec3& first, const Vec3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vec3 scale(const Vec3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double norm(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

double coordinate(const Vec3& value, const std::size_t axis) {
    if (axis == 0) {
        return value.x;
    }
    if (axis == 1) {
        return value.y;
    }
    return value.z;
}

void setCoordinate(Vec3& value,
                   const std::size_t axis,
                   const double coordinateValue) {
    if (axis == 0) {
        value.x = coordinateValue;
    } else if (axis == 1) {
        value.y = coordinateValue;
    } else {
        value.z = coordinateValue;
    }
}

bool samePosition(const SceneFluidClippedVertex& first,
                  const SceneFluidClippedVertex& second) {
    return first.positionMeters.x == second.positionMeters.x
        && first.positionMeters.y == second.positionMeters.y
        && first.positionMeters.z == second.positionMeters.z;
}

void removeConsecutiveDuplicates(
    std::vector<SceneFluidClippedVertex>& polygon) {
    if (polygon.empty()) {
        return;
    }
    std::vector<SceneFluidClippedVertex> unique;
    unique.reserve(polygon.size());
    for (const auto& vertex : polygon) {
        if (unique.empty() || !samePosition(unique.back(), vertex)) {
            unique.push_back(vertex);
        }
    }
    if (unique.size() > 1 && samePosition(unique.front(), unique.back())) {
        unique.pop_back();
    }
    polygon = std::move(unique);
}

SceneFluidClippedVertex interpolate(
    const SceneFluidClippedVertex& first,
    const SceneFluidClippedVertex& second,
    const double firstDistance,
    const double secondDistance,
    const std::size_t axis,
    const double planeCoordinate) {
    const double denominator = firstDistance - secondDistance;
    if (!std::isfinite(denominator) || denominator == 0.0) {
        throw std::runtime_error(
            "scene fluid cell clipping encountered a parallel crossing");
    }
    const double parameter = std::clamp(
        firstDistance / denominator, 0.0, 1.0);
    SceneFluidClippedVertex result;
    result.positionMeters = add(
        scale(first.positionMeters, 1.0 - parameter),
        scale(second.positionMeters, parameter));
    // A segment already lying on an earlier clip plane must stay there when
    // a later-axis plane splits it. Re-evaluating c*((1-t)+t) can move the
    // shared coordinate by an ulp and erase exact face ownership downstream.
    for (std::size_t preservedAxis = 0;
         preservedAxis < 3; ++preservedAxis) {
        const double firstCoordinate = coordinate(
            first.positionMeters, preservedAxis);
        if (firstCoordinate
            == coordinate(second.positionMeters, preservedAxis)) {
            setCoordinate(result.positionMeters, preservedAxis,
                          firstCoordinate);
        }
    }
    setCoordinate(result.positionMeters, axis, planeCoordinate);
    for (std::size_t corner = 0; corner < 3; ++corner) {
        result.barycentricCoordinates[corner] =
            (1.0 - parameter) * first.barycentricCoordinates[corner]
            + parameter * second.barycentricCoordinates[corner];
    }
    return result;
}

void clipAgainstPlane(
    std::vector<SceneFluidClippedVertex>& polygon,
    const std::size_t axis,
    const double planeCoordinate,
    const bool keepGreater) {
    if (polygon.empty()) {
        return;
    }
    const auto distance = [&](const SceneFluidClippedVertex& vertex) {
        const double difference =
            coordinate(vertex.positionMeters, axis) - planeCoordinate;
        return keepGreater ? difference : -difference;
    };
    std::vector<SceneFluidClippedVertex> clipped;
    clipped.reserve(polygon.size() + 1);
    SceneFluidClippedVertex first = polygon.back();
    double firstDistance = distance(first);
    bool firstInside = firstDistance >= 0.0;
    for (const auto& second : polygon) {
        const double secondDistance = distance(second);
        const bool secondInside = secondDistance >= 0.0;
        if (secondInside != firstInside) {
            clipped.push_back(interpolate(
                first, second, firstDistance, secondDistance,
                axis, planeCoordinate));
        }
        if (secondInside) {
            clipped.push_back(second);
        }
        first = second;
        firstDistance = secondDistance;
        firstInside = secondInside;
    }
    removeConsecutiveDuplicates(clipped);
    polygon = std::move(clipped);
}

std::uint8_t coincidentBoundaryPlanes(
    const std::vector<SceneFluidClippedVertex>& vertices,
    const Vec3& lower,
    const Vec3& upper) {
    std::uint8_t result = CellBoundaryNone;
    const auto allAt = [&](const std::size_t axis, const double value) {
        return std::all_of(
            vertices.begin(), vertices.end(),
            [&](const SceneFluidClippedVertex& vertex) {
                return coordinate(vertex.positionMeters, axis) == value;
            });
    };
    if (allAt(0, lower.x)) {
        result |= CellBoundaryXMinus;
    }
    if (allAt(0, upper.x)) {
        result |= CellBoundaryXPlus;
    }
    if (allAt(1, lower.y)) {
        result |= CellBoundaryYMinus;
    }
    if (allAt(1, upper.y)) {
        result |= CellBoundaryYPlus;
    }
    if (allAt(2, lower.z)) {
        result |= CellBoundaryZMinus;
    }
    if (allAt(2, upper.z)) {
        result |= CellBoundaryZPlus;
    }
    return result;
}

void measurePatch(SceneFluidTriangleBoxClip& result) {
    const auto& vertices = result.vertices;
    if (vertices.empty()) {
        throw std::runtime_error(
            "scene fluid exact intersection clipped to an empty patch");
    }
    if (vertices.size() == 1) {
        result.dimension = SceneFluidPatchDimension::Point;
        result.centroidMeters = vertices.front().positionMeters;
        result.centroidBarycentricCoordinates =
            vertices.front().barycentricCoordinates;
        return;
    }

    Vec3 weightedCentroid;
    std::array<double, 3> weightedBarycentric{};
    double area = 0.0;
    for (std::size_t index = 1; index + 1 < vertices.size(); ++index) {
        const auto& first = vertices[0];
        const auto& second = vertices[index];
        const auto& third = vertices[index + 1];
        const double triangleArea = 0.5 * norm(cross(
            subtract(second.positionMeters, first.positionMeters),
            subtract(third.positionMeters, first.positionMeters)));
        if (!std::isfinite(triangleArea)) {
            throw std::overflow_error(
                "scene fluid clipped patch area is not finite");
        }
        if (!(triangleArea > 0.0)) {
            continue;
        }
        const Vec3 triangleCentroid = scale(
            add(add(first.positionMeters, second.positionMeters),
                third.positionMeters),
            1.0 / 3.0);
        weightedCentroid = add(
            weightedCentroid, scale(triangleCentroid, triangleArea));
        for (std::size_t corner = 0; corner < 3; ++corner) {
            weightedBarycentric[corner] += triangleArea
                * (first.barycentricCoordinates[corner]
                   + second.barycentricCoordinates[corner]
                   + third.barycentricCoordinates[corner])
                / 3.0;
        }
        area += triangleArea;
    }
    if (area > 0.0) {
        result.dimension = SceneFluidPatchDimension::Area;
        result.areaSquareMeters = area;
        result.centroidMeters = scale(weightedCentroid, 1.0 / area);
        for (std::size_t corner = 0; corner < 3; ++corner) {
            result.centroidBarycentricCoordinates[corner] =
                weightedBarycentric[corner] / area;
        }
        return;
    }

    result.dimension = SceneFluidPatchDimension::Segment;
    for (const auto& vertex : vertices) {
        result.centroidMeters = add(
            result.centroidMeters, vertex.positionMeters);
        for (std::size_t corner = 0; corner < 3; ++corner) {
            result.centroidBarycentricCoordinates[corner] +=
                vertex.barycentricCoordinates[corner];
        }
    }
    const double inverseCount = 1.0 / static_cast<double>(vertices.size());
    result.centroidMeters = scale(
        result.centroidMeters, inverseCount);
    for (double& coordinateValue :
         result.centroidBarycentricCoordinates) {
        coordinateValue *= inverseCount;
    }
}

struct ClippedPatch {
    SceneFluidCellPatch patch;
    std::vector<SceneFluidClippedVertex> vertices;
};

ClippedPatch clippedPatch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidCellIntersection& intersection) {
    const auto& triangle = surface.triangles[intersection.triangleIndex];
    ClippedPatch result;
    result.patch.cellIndex = intersection.cellIndex;
    result.patch.cell = intersection.cell;
    result.patch.triangleIndex = intersection.triangleIndex;
    result.patch.triangleId = intersection.triangleId;
    std::array<Vec3, 3> positions;
    for (std::size_t corner = 0; corner < 3; ++corner) {
        positions[corner] =
            state.vertices[triangle.vertexIndices[corner]].positionMeters;
    }

    const Vector3 gridLower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vec3 lower{
        gridLower.x + static_cast<double>(intersection.cell.i) * spacing.x,
        gridLower.y + static_cast<double>(intersection.cell.j) * spacing.y,
        gridLower.z + static_cast<double>(intersection.cell.k) * spacing.z,
    };
    const Vec3 upper{
        lower.x + spacing.x,
        lower.y + spacing.y,
        lower.z + spacing.z,
    };
    auto clipped = clipSceneFluidTriangleToAxisAlignedBox(
        positions, lower, upper);
    if (!clipped) {
        throw std::runtime_error(
            "scene fluid exact intersection clipped to an empty patch");
    }
    result.patch.dimension = clipped->dimension;
    result.patch.coincidentBoundaryPlanes =
        clipped->coincidentBoundaryPlanes;
    result.patch.areaSquareMeters = clipped->areaSquareMeters;
    result.patch.centroidMeters = clipped->centroidMeters;
    result.patch.centroidBarycentricCoordinates =
        clipped->centroidBarycentricCoordinates;
    result.vertices = std::move(clipped->vertices);
    result.patch.vertexCount = result.vertices.size();
    return result;
}

std::uint64_t patchFingerprint(const SceneFluidGridPatchSet& patches) {
    Fingerprint fingerprint;
    fingerprint.integer(patches.version);
    fingerprint.integer(patches.surfaceDefinitionFingerprint);
    fingerprint.integer(patches.surfaceStateFingerprint);
    fingerprint.integer(patches.intersectionSetFingerprint);
    fingerprint.integer(patches.structureDefinitionFingerprint);
    fingerprint.integer(patches.acceptedStepCount);
    fingerprint.real(patches.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(patches.patches.size()));
    for (const auto& patch : patches.patches) {
        fingerprint.integer(static_cast<std::uint64_t>(patch.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(patch.triangleIndex));
        fingerprint.integer(patch.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(patch.firstVertex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.vertexCount));
        fingerprint.enumeration(patch.dimension);
        fingerprint.integer(patch.coincidentBoundaryPlanes);
        fingerprint.real(patch.areaSquareMeters);
        fingerprint.real(patch.centroidMeters.x);
        fingerprint.real(patch.centroidMeters.y);
        fingerprint.real(patch.centroidMeters.z);
        for (const double value : patch.centroidBarycentricCoordinates) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(patches.vertices.size()));
    for (const auto& vertex : patches.vertices) {
        fingerprint.real(vertex.positionMeters.x);
        fingerprint.real(vertex.positionMeters.y);
        fingerprint.real(vertex.positionMeters.z);
        for (const double value : vertex.barycentricCoordinates) {
            fingerprint.real(value);
        }
    }
    return fingerprint.value();
}

} // namespace

std::optional<SceneFluidTriangleBoxClip>
clipSceneFluidTriangleToAxisAlignedBox(
    const std::array<Vec3, 3>& triangle,
    const Vec3& lowerMeters,
    const Vec3& upperMeters) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double lower = coordinate(lowerMeters, axis);
        const double upper = coordinate(upperMeters, axis);
        if (!std::isfinite(lower) || !std::isfinite(upper)
            || !(lower < upper)) {
            throw std::invalid_argument(
                "scene fluid triangle clip box is invalid");
        }
        for (const Vec3& vertex : triangle) {
            if (!std::isfinite(coordinate(vertex, axis))) {
                throw std::invalid_argument(
                    "scene fluid triangle clip input is not finite");
            }
        }
    }

    SceneFluidTriangleBoxClip result;
    result.vertices.reserve(9);
    for (std::size_t corner = 0; corner < 3; ++corner) {
        SceneFluidClippedVertex vertex;
        vertex.positionMeters = triangle[corner];
        vertex.barycentricCoordinates[corner] = 1.0;
        result.vertices.push_back(vertex);
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
        clipAgainstPlane(
            result.vertices, axis, coordinate(lowerMeters, axis), true);
        clipAgainstPlane(
            result.vertices, axis, coordinate(upperMeters, axis), false);
    }
    removeConsecutiveDuplicates(result.vertices);
    if (result.vertices.empty()) {
        return std::nullopt;
    }
    measurePatch(result);
    result.coincidentBoundaryPlanes = coincidentBoundaryPlanes(
        result.vertices, lowerMeters, upperMeters);
    return result;
}

std::span<const SceneFluidCellPatch>
SceneFluidGridPatchSet::patchesForCell(
    const std::size_t cellIndex) const noexcept {
    const auto first = std::lower_bound(
        patches.begin(), patches.end(), cellIndex,
        [](const SceneFluidCellPatch& patch, const std::size_t expected) {
            return patch.cellIndex < expected;
        });
    const auto last = std::upper_bound(
        first, patches.end(), cellIndex,
        [](const std::size_t expected, const SceneFluidCellPatch& patch) {
            return expected < patch.cellIndex;
        });
    return {first, last};
}

std::span<const SceneFluidClippedVertex>
SceneFluidGridPatchSet::verticesForPatch(
    const SceneFluidCellPatch& patch) const {
    if (patch.firstVertex > vertices.size()
        || patch.vertexCount > vertices.size() - patch.firstVertex) {
        throw std::out_of_range(
            "scene fluid patch vertex range is out of bounds");
    }
    return std::span<const SceneFluidClippedVertex>(
        vertices).subspan(patch.firstVertex, patch.vertexCount);
}

SceneFluidGridPatchSet clipSceneFluidSurfaceToCells(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchLimits& limits) {
    validateSceneFluidGridIntersections(
        intersections, surface, state, grid, candidates);
    if (intersections.settings.separationToleranceMeters != 0.0) {
        throw std::invalid_argument(
            "scene fluid cell clipping requires exact zero-tolerance intersections");
    }
    if (intersections.intersections.size() > limits.maximumPatches) {
        throw std::length_error(
            "scene fluid clipped patches exceed their count limit");
    }

    std::size_t vertexCount = 0;
    for (const auto& intersection : intersections.intersections) {
        const ClippedPatch clipped = clippedPatch(
            surface, state, grid, intersection);
        if (clipped.vertices.size()
                > std::numeric_limits<std::size_t>::max() - vertexCount) {
            throw std::length_error(
                "scene fluid clipped vertex count overflows");
        }
        vertexCount += clipped.vertices.size();
        if (vertexCount > limits.maximumVertices) {
            throw std::length_error(
                "scene fluid clipped vertices exceed their count limit");
        }
    }
    std::size_t patchBytes = 0;
    std::size_t vertexBytes = 0;
    if (!checkedMultiply(intersections.intersections.size(),
                         sizeof(SceneFluidCellPatch), patchBytes)
        || !checkedMultiply(vertexCount,
                            sizeof(SceneFluidClippedVertex), vertexBytes)
        || vertexBytes > std::numeric_limits<std::size_t>::max() - patchBytes
        || patchBytes + vertexBytes > limits.maximumPatchBytes) {
        throw std::length_error(
            "scene fluid clipped patch storage exceeds its byte limit");
    }

    SceneFluidGridPatchSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.intersectionSetFingerprint = intersections.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.patches.reserve(intersections.intersections.size());
    result.vertices.reserve(vertexCount);
    for (const auto& intersection : intersections.intersections) {
        ClippedPatch clipped = clippedPatch(
            surface, state, grid, intersection);
        clipped.patch.firstVertex = result.vertices.size();
        result.patches.push_back(clipped.patch);
        result.vertices.insert(
            result.vertices.end(),
            clipped.vertices.begin(), clipped.vertices.end());
    }
    result.fingerprint = patchFingerprint(result);
    validateSceneFluidGridPatches(
        result, surface, state, grid, candidates, intersections);
    return result;
}

void validateSceneFluidGridPatches(
    const SceneFluidGridPatchSet& patches,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections) {
    validateSceneFluidGridIntersections(
        intersections, surface, state, grid, candidates);
    if (intersections.settings.separationToleranceMeters != 0.0
        || patches.version != sceneFluidGridPatchVersion
        || patches.fingerprint == 0
        || patches.surfaceDefinitionFingerprint != surface.fingerprint
        || patches.surfaceStateFingerprint != state.fingerprint
        || patches.intersectionSetFingerprint != intersections.fingerprint
        || patches.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || patches.acceptedStepCount != state.acceptedStepCount
        || patches.simulationTimeSeconds != state.simulationTimeSeconds
        || patches.patches.size() != intersections.intersections.size()) {
        throw std::invalid_argument(
            "scene fluid clipped patch identity is invalid");
    }

    std::size_t expectedFirstVertex = 0;
    for (std::size_t index = 0; index < patches.patches.size(); ++index) {
        const auto& intersection = intersections.intersections[index];
        const auto& patch = patches.patches[index];
        const ClippedPatch expected = clippedPatch(
            surface, state, grid, intersection);
        SceneFluidCellPatch expectedPatch = expected.patch;
        expectedPatch.firstVertex = expectedFirstVertex;
        if (patch != expectedPatch
            || patch.firstVertex > patches.vertices.size()
            || patch.vertexCount
                > patches.vertices.size() - patch.firstVertex) {
            throw std::invalid_argument(
                "scene fluid clipped patch does not match its exact intersection");
        }
        const auto actualVertices = patches.verticesForPatch(patch);
        if (!std::equal(actualVertices.begin(), actualVertices.end(),
                        expected.vertices.begin(), expected.vertices.end())) {
            throw std::invalid_argument(
                "scene fluid clipped patch vertices are invalid");
        }
        expectedFirstVertex += patch.vertexCount;
    }
    if (expectedFirstVertex != patches.vertices.size()
        || patches.fingerprint != patchFingerprint(patches)) {
        throw std::invalid_argument(
            "scene fluid clipped patch storage or fingerprint is invalid");
    }
}

} // namespace simwing::fsi::fluid
