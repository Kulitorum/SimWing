#include "fluid/scene_surface_intersection.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

struct Vector {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool operator==(const Vector&) const = default;
};

Vector subtract(const Vector& first, const Vector& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vector cross(const Vector& first, const Vector& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double dot(const Vector& first, const Vector& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

double norm(const Vector& value) {
    return std::sqrt(dot(value, value));
}

Vector normalized(const Vector& value) {
    const double scale = std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
    if (!std::isfinite(scale) || !(scale > 0.0)) {
        return {};
    }
    const Vector scaled{
        value.x / scale,
        value.y / scale,
        value.z / scale,
    };
    const double length = norm(scaled);
    if (!std::isfinite(length) || !(length > 0.0)) {
        return {};
    }
    return {scaled.x / length,
            scaled.y / length,
            scaled.z / length};
}

Vector converted(const Vec3& value) {
    return {value.x, value.y, value.z};
}

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

void validateSettings(
    const SceneFluidGridIntersectionSettings& settings,
    const SceneFluidGridCandidateSet& candidates) {
    if (!std::isfinite(settings.separationToleranceMeters)
        || settings.separationToleranceMeters < 0.0
        || settings.separationToleranceMeters
            > candidates.settings.boundingPaddingMeters) {
        throw std::invalid_argument(
            "scene fluid grid-intersection settings exceed broad-phase support");
    }
}

bool overlapsOnAxis(
    const std::array<Vector, 3>& centeredTriangle,
    const Vector& halfExtent,
    const Vector& rawAxis,
    const double tolerance) {
    const Vector axis = normalized(rawAxis);
    if (axis == Vector{}) {
        return true;
    }
    double minimum = dot(centeredTriangle[0], axis);
    double maximum = minimum;
    for (std::size_t corner = 1; corner < 3; ++corner) {
        const double projection = dot(centeredTriangle[corner], axis);
        minimum = std::min(minimum, projection);
        maximum = std::max(maximum, projection);
    }
    const double radius = halfExtent.x * std::abs(axis.x)
        + halfExtent.y * std::abs(axis.y)
        + halfExtent.z * std::abs(axis.z);
    return minimum <= radius + tolerance
        && maximum >= -radius - tolerance;
}

bool triangleIntersectsCell(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidCellCandidate& candidate,
    const double tolerance) {
    const auto& triangle = surface.triangles[candidate.triangleIndex];
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vector3 lower = grid.lowerMeters();
    const Vector center{
        lower.x + (static_cast<double>(candidate.cell.i) + 0.5) * spacing.x,
        lower.y + (static_cast<double>(candidate.cell.j) + 0.5) * spacing.y,
        lower.z + (static_cast<double>(candidate.cell.k) + 0.5) * spacing.z,
    };
    const Vector halfExtent{
        0.5 * spacing.x,
        0.5 * spacing.y,
        0.5 * spacing.z,
    };
    std::array<Vector, 3> vertices{};
    for (std::size_t corner = 0; corner < 3; ++corner) {
        vertices[corner] = subtract(
            converted(state.vertices[
                triangle.vertexIndices[corner]].positionMeters),
            center);
    }
    const std::array<Vector, 3> edges{
        subtract(vertices[1], vertices[0]),
        subtract(vertices[2], vertices[1]),
        subtract(vertices[0], vertices[2]),
    };
    const Vector normal = cross(
        normalized(edges[0]),
        normalized(subtract(vertices[2], vertices[0])));
    if (!(norm(normal) > 0.0)) {
        throw std::invalid_argument(
            "scene fluid grid intersection received a degenerate triangle");
    }

    const std::array<Vector, 3> boxAxes{
        Vector{1.0, 0.0, 0.0},
        Vector{0.0, 1.0, 0.0},
        Vector{0.0, 0.0, 1.0},
    };
    for (const Vector& axis : boxAxes) {
        if (!overlapsOnAxis(vertices, halfExtent, axis, tolerance)) {
            return false;
        }
    }
    if (!overlapsOnAxis(vertices, halfExtent, normal, tolerance)) {
        return false;
    }
    for (const Vector& edge : edges) {
        for (const Vector& axis : boxAxes) {
            if (!overlapsOnAxis(
                    vertices, halfExtent, cross(edge, axis), tolerance)) {
                return false;
            }
        }
    }
    return true;
}

std::uint64_t intersectionFingerprint(
    const SceneFluidGridIntersectionSet& intersections) {
    Fingerprint fingerprint;
    fingerprint.integer(intersections.version);
    fingerprint.integer(intersections.surfaceDefinitionFingerprint);
    fingerprint.integer(intersections.surfaceStateFingerprint);
    fingerprint.integer(intersections.candidateSetFingerprint);
    fingerprint.integer(intersections.structureDefinitionFingerprint);
    fingerprint.integer(intersections.acceptedStepCount);
    fingerprint.real(intersections.simulationTimeSeconds);
    fingerprint.real(intersections.settings.separationToleranceMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        intersections.triangleIntersectionCounts.size()));
    for (const std::size_t count : intersections.triangleIntersectionCounts) {
        fingerprint.integer(static_cast<std::uint64_t>(count));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        intersections.intersections.size()));
    for (const auto& intersection : intersections.intersections) {
        fingerprint.integer(static_cast<std::uint64_t>(
            intersection.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(intersection.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(intersection.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(intersection.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            intersection.triangleIndex));
        fingerprint.integer(intersection.triangleId);
    }
    return fingerprint.value();
}

bool samePair(const SceneFluidCellIntersection& intersection,
              const SceneFluidCellCandidate& candidate) {
    return intersection.cellIndex == candidate.cellIndex
        && intersection.cell == candidate.cell
        && intersection.triangleIndex == candidate.triangleIndex
        && intersection.triangleId == candidate.triangleId;
}

} // namespace

std::span<const SceneFluidCellIntersection>
SceneFluidGridIntersectionSet::intersectionsForCell(
    const std::size_t cellIndex) const noexcept {
    const auto first = std::lower_bound(
        intersections.begin(), intersections.end(), cellIndex,
        [](const SceneFluidCellIntersection& intersection,
           const std::size_t expected) {
            return intersection.cellIndex < expected;
        });
    const auto last = std::upper_bound(
        first, intersections.end(), cellIndex,
        [](const std::size_t expected,
           const SceneFluidCellIntersection& intersection) {
            return expected < intersection.cellIndex;
        });
    return {first, last};
}

SceneFluidGridIntersectionSet intersectSceneFluidSurfaceWithGrid(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSettings& settings,
    const SceneFluidGridIntersectionLimits& limits) {
    validateSceneFluidGridCandidates(candidates, surface, state, grid);
    validateSettings(settings, candidates);

    std::size_t intersectionCount = 0;
    for (const auto& candidate : candidates.candidates) {
        if (triangleIntersectsCell(
                surface, state, grid, candidate,
                settings.separationToleranceMeters)) {
            if (intersectionCount == limits.maximumIntersections) {
                throw std::length_error(
                    "scene fluid grid intersections exceed their count limit");
            }
            ++intersectionCount;
        }
    }
    std::size_t intersectionBytes = 0;
    std::size_t countBytes = 0;
    if (!checkedMultiply(
            intersectionCount,
            sizeof(SceneFluidCellIntersection),
            intersectionBytes)
        || !checkedMultiply(
            surface.triangles.size(), sizeof(std::size_t), countBytes)
        || countBytes > std::numeric_limits<std::size_t>::max()
            - intersectionBytes
        || countBytes + intersectionBytes
            > limits.maximumIntersectionBytes) {
        throw std::length_error(
            "scene fluid grid intersections exceed their byte limit");
    }

    SceneFluidGridIntersectionSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.candidateSetFingerprint = candidates.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.settings = settings;
    result.triangleIntersectionCounts.assign(
        surface.triangles.size(), 0);
    result.intersections.reserve(intersectionCount);
    for (const auto& candidate : candidates.candidates) {
        if (!triangleIntersectsCell(
                surface, state, grid, candidate,
                settings.separationToleranceMeters)) {
            continue;
        }
        result.intersections.push_back({
            candidate.cellIndex,
            candidate.cell,
            candidate.triangleIndex,
            candidate.triangleId,
        });
        ++result.triangleIntersectionCounts[candidate.triangleIndex];
    }
    result.fingerprint = intersectionFingerprint(result);
    validateSceneFluidGridIntersections(
        result, surface, state, grid, candidates);
    return result;
}

void validateSceneFluidGridIntersections(
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates) {
    validateSceneFluidGridCandidates(candidates, surface, state, grid);
    validateSettings(intersections.settings, candidates);
    if (intersections.version != sceneFluidGridIntersectionVersion
        || intersections.fingerprint == 0
        || intersections.surfaceDefinitionFingerprint != surface.fingerprint
        || intersections.surfaceStateFingerprint != state.fingerprint
        || intersections.candidateSetFingerprint != candidates.fingerprint
        || intersections.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || intersections.acceptedStepCount != state.acceptedStepCount
        || intersections.simulationTimeSeconds != state.simulationTimeSeconds
        || intersections.triangleIntersectionCounts.size()
            != surface.triangles.size()) {
        throw std::invalid_argument(
            "scene fluid grid-intersection identity is invalid");
    }

    std::vector<std::size_t> observedCounts(surface.triangles.size(), 0);
    std::size_t intersectionIndex = 0;
    for (const auto& candidate : candidates.candidates) {
        const bool expected = triangleIntersectsCell(
            surface, state, grid, candidate,
            intersections.settings.separationToleranceMeters);
        const bool present = intersectionIndex < intersections.intersections.size()
            && samePair(
                intersections.intersections[intersectionIndex], candidate);
        if (expected != present) {
            throw std::invalid_argument(
                "scene fluid grid intersections do not match the narrow phase");
        }
        if (present) {
            ++observedCounts[candidate.triangleIndex];
            ++intersectionIndex;
        }
    }
    if (intersectionIndex != intersections.intersections.size()
        || observedCounts != intersections.triangleIntersectionCounts) {
        throw std::invalid_argument(
            "scene fluid grid intersections are incomplete or out of order");
    }
    if (intersections.fingerprint != intersectionFingerprint(intersections)) {
        throw std::invalid_argument(
            "scene fluid grid-intersection fingerprint does not match its payload");
    }
}

} // namespace simwing::fsi::fluid
