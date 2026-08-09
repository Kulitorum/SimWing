#include "scene_fluid_opening_cap.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <set>
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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
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

struct EdgeIncidence {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t negativeRegion = 0;
    std::size_t positiveRegion = 0;
};

using Edge = std::pair<std::size_t, std::size_t>;

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

Vec3 scaled(const Vec3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

Vec3 add(const Vec3& first, const Vec3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

Edge edgeKey(const std::size_t first, const std::size_t second) {
    return {std::min(first, second), std::max(first, second)};
}

void validateSettings(const SceneFluidOpeningCapSettings& settings) {
    if (!std::isfinite(settings.planarityToleranceMeters)
        || settings.planarityToleranceMeters < 0.0
        || !std::isfinite(settings.minimumTriangleAreaSquareMeters)
        || !(settings.minimumTriangleAreaSquareMeters > 0.0)
        || !std::isfinite(settings.convexityTolerance)
        || settings.convexityTolerance < 0.0
        || settings.convexityTolerance >= 1.0) {
        throw std::invalid_argument(
            "scene fluid opening-cap settings are invalid");
    }
}

std::uint64_t capFingerprint(const SceneFluidOpeningCapSet& caps) {
    Fingerprint fingerprint;
    fingerprint.integer(caps.version);
    fingerprint.integer(caps.surfaceDefinitionFingerprint);
    fingerprint.integer(caps.surfaceStateFingerprint);
    fingerprint.integer(caps.structureDefinitionFingerprint);
    fingerprint.integer(caps.acceptedStepCount);
    fingerprint.real(caps.simulationTimeSeconds);
    fingerprint.real(caps.settings.planarityToleranceMeters);
    fingerprint.real(caps.settings.minimumTriangleAreaSquareMeters);
    fingerprint.real(caps.settings.convexityTolerance);
    fingerprint.integer(static_cast<std::uint64_t>(
        caps.separatingBoundaryEdgeCount));
    fingerprint.real(caps.totalAreaSquareMeters);
    fingerprint.integer(static_cast<std::uint64_t>(caps.caps.size()));
    for (const auto& cap : caps.caps) {
        fingerprint.integer(static_cast<std::uint64_t>(cap.openingIndex));
        fingerprint.integer(cap.openingId);
        fingerprint.integer(static_cast<std::uint64_t>(
            cap.negativeSideRegionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            cap.positiveSideRegionIndex));
        fingerprint.enumeration(cap.role);
        fingerprint.integer(static_cast<std::uint64_t>(cap.firstTriangle));
        fingerprint.integer(static_cast<std::uint64_t>(cap.triangleCount));
        for (const double value : {
                 cap.unitNormalNegativeToPositive.x,
                 cap.unitNormalNegativeToPositive.y,
                 cap.unitNormalNegativeToPositive.z,
                 cap.centroidMeters.x, cap.centroidMeters.y,
                 cap.centroidMeters.z, cap.areaSquareMeters}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(caps.triangles.size()));
    for (const auto& triangle : caps.triangles) {
        fingerprint.integer(static_cast<std::uint64_t>(
            triangle.openingIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            triangle.triangleOrdinal));
        for (const std::size_t vertex : triangle.vertexIndices) {
            fingerprint.integer(static_cast<std::uint64_t>(vertex));
        }
        for (const double value : {
                 triangle.unitNormalNegativeToPositive.x,
                 triangle.unitNormalNegativeToPositive.y,
                 triangle.unitNormalNegativeToPositive.z,
                 triangle.centroidMeters.x, triangle.centroidMeters.y,
                 triangle.centroidMeters.z, triangle.areaSquareMeters}) {
            fingerprint.real(value);
        }
    }
    return fingerprint.value();
}

SceneFluidOpeningCapSet buildCaps(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSettings& settings,
    const SceneFluidOpeningCapLimits& limits) {
    validateSettings(settings);
    if (surface.openings.size() > limits.maximumCaps) {
        throw std::length_error(
            "scene fluid opening caps exceed their cap limit");
    }

    SceneFluidOpeningCapSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.settings = settings;

    std::map<Edge, std::vector<EdgeIncidence>> surfaceEdges;
    for (const auto& triangle : surface.triangles) {
        if (triangle.negativeSideRegionIndex
            == triangle.positiveSideRegionIndex) {
            continue;
        }
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const std::size_t from = triangle.vertexIndices[edge];
            const std::size_t to = triangle.vertexIndices[(edge + 1) % 3];
            auto& incidences = surfaceEdges[edgeKey(from, to)];
            incidences.push_back({
                from, to,
                triangle.negativeSideRegionIndex,
                triangle.positiveSideRegionIndex,
            });
            if (surfaceEdges.size() > limits.maximumBoundaryEdges) {
                throw std::length_error(
                    "scene fluid opening-cap surface edges exceed their limit");
            }
        }
    }

    std::map<Edge, EdgeIncidence> boundaryEdges;
    for (const auto& [edge, incidences] : surfaceEdges) {
        if (incidences.size() == 1) {
            boundaryEdges.emplace(edge, incidences.front());
            continue;
        }
        if (incidences.size() != 2
            || incidences[0].from != incidences[1].to
            || incidences[0].to != incidences[1].from
            || incidences[0].negativeRegion != incidences[1].negativeRegion
            || incidences[0].positiveRegion != incidences[1].positiveRegion) {
            throw std::invalid_argument(
                "scene fluid opening caps require consistently wound two-sided manifolds away from openings");
        }
    }
    result.separatingBoundaryEdgeCount = boundaryEdges.size();

    std::set<Edge> claimedBoundaryEdges;
    result.caps.reserve(surface.openings.size());
    for (std::size_t openingIndex = 0;
         openingIndex < surface.openings.size(); ++openingIndex) {
        const auto& opening = surface.openings[openingIndex];
        if (opening.orderedVertexIndices.size() < 3) {
            throw std::invalid_argument(
                "scene fluid opening cap requires at least three vertices");
        }
        bool capOrderMatchesOpening = true;
        bool orientationSet = false;
        for (std::size_t edge = 0;
             edge < opening.orderedVertexIndices.size(); ++edge) {
            const std::size_t from = opening.orderedVertexIndices[edge];
            const std::size_t to = opening.orderedVertexIndices[
                (edge + 1) % opening.orderedVertexIndices.size()];
            const Edge key = edgeKey(from, to);
            const auto found = boundaryEdges.find(key);
            if (found == boundaryEdges.end()
                || found->second.negativeRegion
                    != opening.negativeSideRegionIndex
                || found->second.positiveRegion
                    != opening.positiveSideRegionIndex) {
                throw std::invalid_argument(
                    "scene fluid opening loop does not match one separating surface boundary");
            }
            if (!claimedBoundaryEdges.insert(key).second) {
                throw std::invalid_argument(
                    "scene fluid opening boundary edge is claimed more than once");
            }
            const bool loopMatchesSurface =
                found->second.from == from && found->second.to == to;
            const bool thisCapMatchesOpening = !loopMatchesSurface;
            if (!orientationSet) {
                capOrderMatchesOpening = thisCapMatchesOpening;
                orientationSet = true;
            } else if (capOrderMatchesOpening != thisCapMatchesOpening) {
                throw std::invalid_argument(
                    "scene fluid opening loop orientation conflicts with its surface boundary");
            }
        }

        std::vector<std::size_t> capOrder = opening.orderedVertexIndices;
        if (!capOrderMatchesOpening) {
            std::ranges::reverse(capOrder);
        }
        Vec3 newell;
        Vec3 minimum = state.vertices[capOrder.front()].positionMeters;
        Vec3 maximum = minimum;
        for (std::size_t index = 0; index < capOrder.size(); ++index) {
            const Vec3& current =
                state.vertices[capOrder[index]].positionMeters;
            const Vec3& next = state.vertices[
                capOrder[(index + 1) % capOrder.size()]].positionMeters;
            newell.x += (current.y - next.y) * (current.z + next.z);
            newell.y += (current.z - next.z) * (current.x + next.x);
            newell.z += (current.x - next.x) * (current.y + next.y);
            minimum.x = std::min(minimum.x, current.x);
            minimum.y = std::min(minimum.y, current.y);
            minimum.z = std::min(minimum.z, current.z);
            maximum.x = std::max(maximum.x, current.x);
            maximum.y = std::max(maximum.y, current.y);
            maximum.z = std::max(maximum.z, current.z);
        }
        const double normalLength = length(newell);
        const double scale = length(subtract(maximum, minimum));
        if (!std::isfinite(normalLength) || !(normalLength > 0.0)
            || !std::isfinite(scale) || !(scale > 0.0)) {
            throw std::invalid_argument(
                "scene fluid opening cap has a degenerate boundary");
        }
        const Vec3 normal = scaled(newell, 1.0 / normalLength);
        const Vec3& planePoint =
            state.vertices[capOrder.front()].positionMeters;
        for (const std::size_t vertex : capOrder) {
            const double distance = dot(
                subtract(state.vertices[vertex].positionMeters, planePoint),
                normal);
            if (!std::isfinite(distance)
                || std::abs(distance) > settings.planarityToleranceMeters) {
                throw std::invalid_argument(
                    "scene fluid opening cap is not planar within tolerance");
            }
        }
        for (std::size_t index = 0; index < capOrder.size(); ++index) {
            const Vec3& previous = state.vertices[capOrder[
                (index + capOrder.size() - 1) % capOrder.size()]].positionMeters;
            const Vec3& current =
                state.vertices[capOrder[index]].positionMeters;
            const Vec3& next = state.vertices[
                capOrder[(index + 1) % capOrder.size()]].positionMeters;
            const double turn = dot(
                cross(subtract(current, previous), subtract(next, current)),
                normal);
            if (!std::isfinite(turn)
                || !(turn > settings.convexityTolerance * scale * scale)) {
                throw std::invalid_argument(
                    "scene fluid opening cap requires a strictly convex boundary");
            }
        }

        SceneFluidOpeningCap cap;
        cap.openingIndex = openingIndex;
        cap.openingId = opening.id;
        cap.negativeSideRegionIndex = opening.negativeSideRegionIndex;
        cap.positiveSideRegionIndex = opening.positiveSideRegionIndex;
        cap.role = opening.role;
        cap.firstTriangle = result.triangles.size();
        cap.unitNormalNegativeToPositive = normal;
        Vec3 weightedCentroid;
        for (std::size_t ordinal = 0; ordinal + 2 < capOrder.size();
             ++ordinal) {
            if (result.triangles.size() == limits.maximumCapTriangles) {
                throw std::length_error(
                    "scene fluid opening-cap triangles exceed their limit");
            }
            const std::array<std::size_t, 3> vertices{
                capOrder[0], capOrder[ordinal + 1], capOrder[ordinal + 2]};
            const Vec3& first = state.vertices[vertices[0]].positionMeters;
            const Vec3& second = state.vertices[vertices[1]].positionMeters;
            const Vec3& third = state.vertices[vertices[2]].positionMeters;
            const Vec3 areaVector = cross(
                subtract(second, first), subtract(third, first));
            const double area = 0.5 * length(areaVector);
            if (!std::isfinite(area)
                || !(area > settings.minimumTriangleAreaSquareMeters)
                || !(dot(areaVector, normal) > 0.0)) {
                throw std::invalid_argument(
                    "scene fluid opening-cap triangle is invalid");
            }
            const Vec3 centroid{
                (first.x + second.x + third.x) / 3.0,
                (first.y + second.y + third.y) / 3.0,
                (first.z + second.z + third.z) / 3.0,
            };
            result.triangles.push_back({
                openingIndex, ordinal, vertices, normal, centroid, area});
            cap.areaSquareMeters += area;
            weightedCentroid = add(
                weightedCentroid, scaled(centroid, area));
        }
        cap.triangleCount = result.triangles.size() - cap.firstTriangle;
        cap.centroidMeters = scaled(
            weightedCentroid, 1.0 / cap.areaSquareMeters);
        result.totalAreaSquareMeters += cap.areaSquareMeters;
        result.caps.push_back(cap);
    }
    if (claimedBoundaryEdges.size() != boundaryEdges.size()) {
        throw std::invalid_argument(
            "scene fluid separating surface has an unauthored boundary edge");
    }

    std::size_t capBytes = 0;
    std::size_t triangleBytes = 0;
    std::size_t totalBytes = 0;
    if (!checkedMultiply(result.caps.size(), sizeof(SceneFluidOpeningCap),
                         capBytes)
        || !checkedMultiply(
            result.triangles.size(), sizeof(SceneFluidOpeningCapTriangle),
            triangleBytes)
        || !checkedAdd(capBytes, triangleBytes, totalBytes)
        || totalBytes > limits.maximumCapBytes) {
        throw std::length_error(
            "scene fluid opening-cap result exceeds its byte limit");
    }
    result.fingerprint = capFingerprint(result);
    return result;
}

} // namespace

SceneFluidOpeningCapSet buildSceneFluidOpeningCaps(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSettings& settings,
    const SceneFluidOpeningCapLimits& limits) {
    validateSceneFluidSurfaceState(surface, state);
    auto result = buildCaps(surface, state, settings, limits);
    validateSceneFluidOpeningCaps(result, surface, state);
    return result;
}

void validateSceneFluidOpeningCaps(
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state) {
    validateSceneFluidSurfaceState(surface, state);
    if (caps.version != sceneFluidOpeningCapVersion
        || caps.fingerprint == 0
        || caps.surfaceDefinitionFingerprint != surface.fingerprint
        || caps.surfaceStateFingerprint != state.fingerprint
        || caps.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || caps.acceptedStepCount != state.acceptedStepCount
        || caps.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid opening-cap identity is invalid");
    }
    const SceneFluidOpeningCapLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildCaps(surface, state, caps.settings, unlimited);
    if (caps != expected || caps.fingerprint != capFingerprint(caps)) {
        throw std::invalid_argument(
            "scene fluid opening-cap payload is invalid");
    }
}

} // namespace simwing::fsi
