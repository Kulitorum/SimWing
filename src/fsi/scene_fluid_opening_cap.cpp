#include "scene_fluid_opening_cap.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <sstream>
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

bool validClosedRegionCycle(
    const Edge edge,
    const std::span<const EdgeIncidence> incidences) {
    if (incidences.size() < 2) {
        return false;
    }
    std::map<std::size_t, std::size_t> nextRegion;
    std::set<std::size_t> incomingRegions;
    for (const EdgeIncidence& incidence : incidences) {
        const bool followsCanonicalEdge =
            incidence.from == edge.first && incidence.to == edge.second;
        const bool opposesCanonicalEdge =
            incidence.from == edge.second && incidence.to == edge.first;
        if (!followsCanonicalEdge && !opposesCanonicalEdge) {
            return false;
        }
        const std::size_t fromRegion = followsCanonicalEdge
            ? incidence.negativeRegion : incidence.positiveRegion;
        const std::size_t toRegion = followsCanonicalEdge
            ? incidence.positiveRegion : incidence.negativeRegion;
        if (fromRegion == toRegion
            || !nextRegion.emplace(fromRegion, toRegion).second
            || !incomingRegions.insert(toRegion).second) {
            return false;
        }
    }
    if (nextRegion.size() != incomingRegions.size()) {
        return false;
    }
    const std::size_t start = nextRegion.begin()->first;
    std::set<std::size_t> visited;
    std::size_t current = start;
    for (std::size_t step = 0; step < nextRegion.size(); ++step) {
        if (!visited.insert(current).second) {
            return false;
        }
        const auto next = nextRegion.find(current);
        if (next == nextRegion.end()) {
            return false;
        }
        current = next->second;
    }
    return current == start && visited.size() == nextRegion.size();
}

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

template<typename PositionAt>
bool collapsedLoopGeometry(
    const std::vector<std::size_t>& loop,
    PositionAt&& positionAt,
    const double relativeAreaTolerance) {
    if (loop.size() < 3) {
        return false;
    }
    Vec3 minimum = positionAt(loop.front());
    Vec3 maximum = minimum;
    for (const std::size_t vertex : loop) {
        const Vec3 position = positionAt(vertex);
        minimum.x = std::min(minimum.x, position.x);
        minimum.y = std::min(minimum.y, position.y);
        minimum.z = std::min(minimum.z, position.z);
        maximum.x = std::max(maximum.x, position.x);
        maximum.y = std::max(maximum.y, position.y);
        maximum.z = std::max(maximum.z, position.z);
    }
    const double scale = length(subtract(maximum, minimum));
    if (!std::isfinite(scale) || !(scale > 0.0)) {
        return false;
    }

    const Vec3 origin = positionAt(loop.front());
    Vec3 normalizedAreaVector;
    for (std::size_t index = 1; index + 1 < loop.size(); ++index) {
        const Vec3 firstDelta = subtract(
            positionAt(loop[index]), origin);
        const Vec3 secondDelta = subtract(
            positionAt(loop[index + 1]), origin);
        const Vec3 first{firstDelta.x / scale,
                         firstDelta.y / scale,
                         firstDelta.z / scale};
        const Vec3 second{secondDelta.x / scale,
                          secondDelta.y / scale,
                          secondDelta.z / scale};
        const Vec3 contribution = cross(first, second);
        normalizedAreaVector.x += contribution.x;
        normalizedAreaVector.y += contribution.y;
        normalizedAreaVector.z += contribution.z;
    }
    const double relativeArea = 0.5 * length(normalizedAreaVector);
    return std::isfinite(relativeArea)
        && relativeArea <= relativeAreaTolerance;
}

bool validCollapsedBoundaryComponents(
    const std::map<Edge, std::vector<EdgeIncidence>>& unresolvedEdges,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const double relativeAreaTolerance) {
    using RegionPair = std::pair<std::size_t, std::size_t>;
    std::map<RegionPair, std::vector<EdgeIncidence>> regionBoundaries;
    for (const auto& [edge, incidences] : unresolvedEdges) {
        static_cast<void>(edge);
        if (incidences.size() != 1
            || incidences.front().negativeRegion
                == incidences.front().positiveRegion) {
            return false;
        }
        const EdgeIncidence& incidence = incidences.front();
        regionBoundaries[{incidence.negativeRegion,
                          incidence.positiveRegion}]
            .push_back(incidence);
    }

    for (const auto& [regions, incidences] : regionBoundaries) {
        static_cast<void>(regions);
        std::map<std::size_t, std::size_t> nextVertex;
        std::set<std::size_t> incomingVertices;
        for (const EdgeIncidence& incidence : incidences) {
            if (incidence.from == incidence.to
                || !nextVertex.emplace(
                        incidence.from, incidence.to).second
                || !incomingVertices.insert(incidence.to).second) {
                return false;
            }
        }
        if (nextVertex.size() != incomingVertices.size()) {
            return false;
        }

        std::set<std::size_t> unvisited;
        for (const auto& [from, to] : nextVertex) {
            static_cast<void>(to);
            unvisited.insert(from);
        }
        while (!unvisited.empty()) {
            const std::size_t start = *unvisited.begin();
            std::vector<std::size_t> loop;
            std::size_t current = start;
            do {
                if (!unvisited.erase(current)) {
                    return false;
                }
                loop.push_back(current);
                const auto next = nextVertex.find(current);
                if (next == nextVertex.end()) {
                    return false;
                }
                current = next->second;
            } while (current != start);

            if (!collapsedLoopGeometry(
                    loop,
                    [&surface](const std::size_t vertex) {
                        return surface.vertices[vertex]
                            .referencePositionMeters;
                    },
                    relativeAreaTolerance)
                || !collapsedLoopGeometry(
                    loop,
                    [&state](const std::size_t vertex) {
                        return state.vertices[vertex].positionMeters;
                    },
                    relativeAreaTolerance)) {
                return false;
            }
        }
    }
    return true;
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

struct LoopGeometry {
    Vec3 unitNormal;
    double scaleMeters = 0.0;
};

LoopGeometry loopGeometry(const std::span<const Vec3> positions) {
    Vec3 newell;
    Vec3 minimum = positions.front();
    Vec3 maximum = minimum;
    for (std::size_t index = 0; index < positions.size(); ++index) {
        const Vec3& current = positions[index];
        const Vec3& next = positions[(index + 1) % positions.size()];
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
    return {normal, scale};
}

void validatePlanarity(
    const std::span<const Vec3> positions,
    const LoopGeometry geometry,
    const SceneFluidOpeningCapSettings& settings) {
    const Vec3& planePoint = positions.front();
    for (const auto& position : positions) {
        const double distance = dot(
            subtract(position, planePoint), geometry.unitNormal);
        if (!std::isfinite(distance)
            || std::abs(distance) > settings.planarityToleranceMeters) {
            throw std::invalid_argument(
                "scene fluid opening cap is not planar within tolerance");
        }
    }
}

double signedTurn(const Vec3& first,
                  const Vec3& second,
                  const Vec3& third,
                  const Vec3& normal) {
    return dot(
        cross(subtract(second, first), subtract(third, second)),
        normal);
}

bool segmentIntersection(
    const Vec3& first,
    const Vec3& second,
    const Vec3& third,
    const Vec3& fourth,
    const Vec3& normal,
    const double tolerance) {
    const double firstThird = signedTurn(first, second, third, normal);
    const double firstFourth = signedTurn(first, second, fourth, normal);
    const double thirdFirst = signedTurn(third, fourth, first, normal);
    const double thirdSecond = signedTurn(third, fourth, second, normal);
    const auto withinSegment = [normal, tolerance](
        const Vec3& begin, const Vec3& end, const Vec3& point) {
        const auto projected = [normal](const Vec3& value) {
            return subtract(value, scaled(normal, dot(value, normal)));
        };
        return dot(projected(subtract(point, begin)),
                   projected(subtract(point, end))) <= tolerance;
    };
    if ((std::abs(firstThird) <= tolerance
            && withinSegment(first, second, third))
        || (std::abs(firstFourth) <= tolerance
            && withinSegment(first, second, fourth))
        || (std::abs(thirdFirst) <= tolerance
            && withinSegment(third, fourth, first))
        || (std::abs(thirdSecond) <= tolerance
            && withinSegment(third, fourth, second))) {
        return true;
    }
    const auto strictlyStraddles = [tolerance](
        const double firstSide, const double secondSide) {
        return (firstSide > tolerance && secondSide < -tolerance)
            || (firstSide < -tolerance && secondSide > tolerance);
    };
    return strictlyStraddles(firstThird, firstFourth)
        && strictlyStraddles(thirdFirst, thirdSecond);
}

void validateSimpleLoop(
    const std::span<const Vec3> positions,
    const LoopGeometry geometry,
    const SceneFluidOpeningCapSettings& settings,
    std::size_t& intersectionTestCount,
    const std::size_t maximumIntersectionTests,
    const bool requireNondegenerateTurns) {
    const double tolerance = settings.convexityTolerance
        * geometry.scaleMeters * geometry.scaleMeters;
    for (std::size_t index = 0; index < positions.size(); ++index) {
        const Vec3& previous = positions[
            (index + positions.size() - 1) % positions.size()];
        const Vec3& current = positions[index];
        const Vec3& next = positions[(index + 1) % positions.size()];
        const Vec3 incoming = subtract(current, previous);
        const Vec3 outgoing = subtract(next, current);
        const Vec3 projectedIncoming = subtract(
            incoming,
            scaled(geometry.unitNormal,
                   dot(incoming, geometry.unitNormal)));
        const Vec3 projectedOutgoing = subtract(
            outgoing,
            scaled(geometry.unitNormal,
                   dot(outgoing, geometry.unitNormal)));
        const double turn = signedTurn(
            previous, current, next, geometry.unitNormal);
        if (!std::isfinite(turn)
            || dot(projectedIncoming, projectedIncoming) <= tolerance
            || dot(projectedOutgoing, projectedOutgoing) <= tolerance
            || (std::abs(turn) <= tolerance
                && (requireNondegenerateTurns
                    || dot(projectedIncoming, projectedOutgoing) <= 0.0))) {
            throw std::invalid_argument(
                "scene fluid opening cap has a degenerate boundary turn");
        }
    }
    for (std::size_t first = 0; first < positions.size(); ++first) {
        const std::size_t firstNext = (first + 1) % positions.size();
        for (std::size_t second = first + 1;
             second < positions.size(); ++second) {
            const std::size_t secondNext =
                (second + 1) % positions.size();
            if (first == second || firstNext == second
                || secondNext == first) {
                continue;
            }
            if (intersectionTestCount == maximumIntersectionTests) {
                throw std::length_error(
                    "scene fluid opening-cap boundary intersection tests exceed their limit");
            }
            ++intersectionTestCount;
            if (segmentIntersection(
                    positions[first], positions[firstNext],
                    positions[second], positions[secondNext],
                    geometry.unitNormal, tolerance)) {
                throw std::invalid_argument(
                    "scene fluid opening cap boundary is not simple");
            }
        }
    }
}

bool pointInsideTriangle(
    const Vec3& point,
    const Vec3& first,
    const Vec3& second,
    const Vec3& third,
    const Vec3& normal,
    const double tolerance) {
    return signedTurn(first, second, point, normal) >= -tolerance
        && signedTurn(second, third, point, normal) >= -tolerance
        && signedTurn(third, first, point, normal) >= -tolerance;
}

std::vector<std::array<std::size_t, 3>> triangulateReferenceLoop(
    const std::span<const std::size_t> vertexIndices,
    const std::span<const Vec3> positions,
    const LoopGeometry geometry,
    const SceneFluidOpeningCapSettings& settings,
    std::size_t& pointTestCount,
    const std::size_t maximumPointTests) {
    const double tolerance = settings.convexityTolerance
        * geometry.scaleMeters * geometry.scaleMeters;
    bool strictlyConvex = true;
    for (std::size_t index = 0; index < positions.size(); ++index) {
        strictlyConvex = strictlyConvex
            && signedTurn(
                   positions[(index + positions.size() - 1)
                       % positions.size()],
                   positions[index],
                   positions[(index + 1) % positions.size()],
                   geometry.unitNormal) > tolerance;
    }
    std::vector<std::array<std::size_t, 3>> result;
    result.reserve(vertexIndices.size() - 2);
    if (strictlyConvex) {
        for (std::size_t ordinal = 0; ordinal + 2 < vertexIndices.size();
             ++ordinal) {
            result.push_back({
                vertexIndices[0], vertexIndices[ordinal + 1],
                vertexIndices[ordinal + 2]});
        }
        return result;
    }

    std::vector<std::size_t> active(positions.size());
    for (std::size_t index = 0; index < active.size(); ++index) {
        active[index] = index;
    }
    while (active.size() > 3) {
        bool clipped = false;
        for (std::size_t activeIndex = 0;
             activeIndex < active.size(); ++activeIndex) {
            const std::size_t previous = active[
                (activeIndex + active.size() - 1) % active.size()];
            const std::size_t current = active[activeIndex];
            const std::size_t next = active[
                (activeIndex + 1) % active.size()];
            if (!(signedTurn(
                    positions[previous], positions[current], positions[next],
                    geometry.unitNormal) > tolerance)) {
                continue;
            }
            bool containsVertex = false;
            for (const std::size_t candidate : active) {
                if (candidate == previous || candidate == current
                    || candidate == next) {
                    continue;
                }
                if (pointTestCount == maximumPointTests) {
                    throw std::length_error(
                        "scene fluid opening-cap triangulation point tests exceed their limit");
                }
                ++pointTestCount;
                if (pointInsideTriangle(
                        positions[candidate], positions[previous],
                        positions[current], positions[next],
                        geometry.unitNormal, tolerance)) {
                    containsVertex = true;
                    break;
                }
            }
            if (containsVertex) {
                continue;
            }
            result.push_back({
                vertexIndices[previous], vertexIndices[current],
                vertexIndices[next]});
            active.erase(active.begin()
                         + static_cast<std::ptrdiff_t>(activeIndex));
            clipped = true;
            break;
        }
        if (!clipped) {
            throw std::invalid_argument(
                "scene fluid opening cap triangulation is ambiguous");
        }
    }
    if (!(signedTurn(
            positions[active[0]], positions[active[1]], positions[active[2]],
            geometry.unitNormal) > tolerance)) {
        throw std::invalid_argument(
            "scene fluid opening cap triangulation is invalid");
    }
    result.push_back({
        vertexIndices[active[0]], vertexIndices[active[1]],
        vertexIndices[active[2]]});
    return result;
}

void validateSettings(const SceneFluidOpeningCapSettings& settings) {
    if (!std::isfinite(settings.planarityToleranceMeters)
        || settings.planarityToleranceMeters < 0.0
        || !std::isfinite(settings.minimumTriangleAreaSquareMeters)
        || !(settings.minimumTriangleAreaSquareMeters > 0.0)
        || !std::isfinite(settings.convexityTolerance)
        || settings.convexityTolerance < 0.0
        || settings.convexityTolerance >= 1.0
        || !std::isfinite(
            settings.collapsedBoundaryRelativeAreaTolerance)
        || settings.collapsedBoundaryRelativeAreaTolerance < 0.0
        || settings.collapsedBoundaryRelativeAreaTolerance >= 1.0) {
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
    fingerprint.real(
        caps.settings.collapsedBoundaryRelativeAreaTolerance);
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

    std::size_t incompleteSurfaceEdgeCount = 0;
    for (const auto& [edge, incidences] : surfaceEdges) {
        if (!validClosedRegionCycle(edge, incidences)) {
            ++incompleteSurfaceEdgeCount;
        }
    }
    result.separatingBoundaryEdgeCount = incompleteSurfaceEdgeCount;

    auto resolvedEdges = surfaceEdges;
    std::size_t boundaryIntersectionTestCount = 0;
    std::size_t triangulationPointTestCount = 0;
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
            const auto found = surfaceEdges.find(key);
            if (found == surfaceEdges.end()) {
                throw std::invalid_argument(
                    "scene fluid opening loop edge is absent from the material surface");
            }
            for (const EdgeIncidence& incidence : found->second) {
                if (incidence.negativeRegion
                        != opening.negativeSideRegionIndex
                    || incidence.positiveRegion
                        != opening.positiveSideRegionIndex) {
                    continue;
                }
                const bool loopMatchesSurface =
                    incidence.from == from && incidence.to == to;
                const bool thisCapMatchesOpening = !loopMatchesSurface;
                if (!orientationSet) {
                    capOrderMatchesOpening = thisCapMatchesOpening;
                    orientationSet = true;
                } else if (capOrderMatchesOpening
                           != thisCapMatchesOpening) {
                    throw std::invalid_argument(
                        "scene fluid opening loop orientation conflicts with its matching material boundary");
                }
            }
        }
        if (!orientationSet) {
            throw std::invalid_argument(
                "scene fluid opening has no material edge with the same region pair to fix cap winding");
        }

        std::vector<std::size_t> capOrder = opening.orderedVertexIndices;
        if (!capOrderMatchesOpening) {
            std::ranges::reverse(capOrder);
        }
        std::vector<Vec3> currentPositions;
        std::vector<Vec3> referencePositions;
        currentPositions.reserve(capOrder.size());
        referencePositions.reserve(capOrder.size());
        for (const std::size_t vertex : capOrder) {
            currentPositions.push_back(
                state.vertices[vertex].positionMeters);
            referencePositions.push_back(
                surface.vertices[vertex].referencePositionMeters);
        }
        const bool authoredTriangulation =
            !opening.capTriangleVertexIndices.empty();
        const auto currentGeometry = loopGeometry(currentPositions);
        const auto referenceGeometry = loopGeometry(referencePositions);
        if (!authoredTriangulation) {
            validatePlanarity(
                currentPositions, currentGeometry, settings);
            validatePlanarity(
                referencePositions, referenceGeometry, settings);
        }
        validateSimpleLoop(
            currentPositions, currentGeometry, settings,
            boundaryIntersectionTestCount,
            limits.maximumBoundaryIntersectionTests,
            !authoredTriangulation);
        validateSimpleLoop(
            referencePositions, referenceGeometry, settings,
            boundaryIntersectionTestCount,
            limits.maximumBoundaryIntersectionTests,
            !authoredTriangulation);
        const std::size_t requiredTriangles = capOrder.size() - 2;
        if (result.triangles.size() > limits.maximumCapTriangles
            || requiredTriangles
                > limits.maximumCapTriangles - result.triangles.size()) {
            throw std::length_error(
                "scene fluid opening-cap triangles exceed their limit");
        }
        std::vector<std::array<std::size_t, 3>> triangleVertices;
        if (authoredTriangulation) {
            triangleVertices = opening.capTriangleVertexIndices;
            if (!capOrderMatchesOpening) {
                for (auto& triangle : triangleVertices) {
                    std::swap(triangle[1], triangle[2]);
                }
            }
        } else {
            triangleVertices = triangulateReferenceLoop(
                capOrder, referencePositions, referenceGeometry, settings,
                triangulationPointTestCount,
                limits.maximumTriangulationPointTests);
        }
        const Vec3 normal = currentGeometry.unitNormal;

        SceneFluidOpeningCap cap;
        cap.openingIndex = openingIndex;
        cap.openingId = opening.id;
        cap.negativeSideRegionIndex = opening.negativeSideRegionIndex;
        cap.positiveSideRegionIndex = opening.positiveSideRegionIndex;
        cap.role = opening.role;
        cap.firstTriangle = result.triangles.size();
        cap.unitNormalNegativeToPositive = normal;
        Vec3 weightedCentroid;
        for (std::size_t ordinal = 0;
             ordinal < triangleVertices.size(); ++ordinal) {
            if (result.triangles.size() == limits.maximumCapTriangles) {
                throw std::length_error(
                    "scene fluid opening-cap triangles exceed their limit");
            }
            const auto vertices = triangleVertices[ordinal];
            const Vec3& first = state.vertices[vertices[0]].positionMeters;
            const Vec3& second = state.vertices[vertices[1]].positionMeters;
            const Vec3& third = state.vertices[vertices[2]].positionMeters;
            const Vec3 areaVector = cross(
                subtract(second, first), subtract(third, first));
            const double area = 0.5 * length(areaVector);
            bool referenceValid = true;
            if (authoredTriangulation) {
                const Vec3& referenceFirst =
                    surface.vertices[vertices[0]].referencePositionMeters;
                const Vec3& referenceSecond =
                    surface.vertices[vertices[1]].referencePositionMeters;
                const Vec3& referenceThird =
                    surface.vertices[vertices[2]].referencePositionMeters;
                const Vec3 referenceAreaVector = cross(
                    subtract(referenceSecond, referenceFirst),
                    subtract(referenceThird, referenceFirst));
                const double referenceArea =
                    0.5 * length(referenceAreaVector);
                const double referenceTolerance =
                    settings.convexityTolerance
                    * referenceGeometry.scaleMeters
                    * referenceGeometry.scaleMeters;
                referenceValid = std::isfinite(referenceArea)
                    && referenceArea
                        > settings.minimumTriangleAreaSquareMeters
                    && dot(referenceAreaVector,
                           referenceGeometry.unitNormal)
                        > referenceTolerance;
            }
            const double currentOrientationTolerance =
                authoredTriangulation
                ? settings.convexityTolerance
                    * currentGeometry.scaleMeters
                    * currentGeometry.scaleMeters
                : 0.0;
            if (!std::isfinite(area)
                || !(area > settings.minimumTriangleAreaSquareMeters)
                || !(dot(areaVector, normal)
                     > currentOrientationTolerance)
                || !referenceValid) {
                throw std::invalid_argument(
                    "scene fluid opening-cap triangle is invalid");
            }
            const Vec3 triangleNormal = authoredTriangulation
                ? scaled(areaVector, 0.5 / area)
                : normal;
            const Vec3 centroid{
                (first.x + second.x + third.x) / 3.0,
                (first.y + second.y + third.y) / 3.0,
                (first.z + second.z + third.z) / 3.0,
            };
            result.triangles.push_back({
                openingIndex, ordinal, vertices, triangleNormal,
                centroid, area});
            for (std::size_t edge = 0; edge < 3; ++edge) {
                const std::size_t from = vertices[edge];
                const std::size_t to = vertices[(edge + 1) % 3];
                resolvedEdges[edgeKey(from, to)].push_back({
                    from, to,
                    opening.negativeSideRegionIndex,
                    opening.positiveSideRegionIndex,
                });
                if (resolvedEdges.size() > limits.maximumBoundaryEdges) {
                    throw std::length_error(
                        "scene fluid opening-cap resolved edges exceed their limit");
                }
            }
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
    std::map<Edge, std::vector<EdgeIncidence>> unresolvedEdges;
    for (const auto& [edge, incidences] : resolvedEdges) {
        if (!validClosedRegionCycle(edge, incidences)) {
            unresolvedEdges.emplace(edge, incidences);
        }
    }
    if (!unresolvedEdges.empty()
        && !validCollapsedBoundaryComponents(
            unresolvedEdges, surface, state,
            settings.collapsedBoundaryRelativeAreaTolerance)) {
        const auto& [edge, incidences] = *unresolvedEdges.begin();
        std::ostringstream message;
        message << "scene fluid material and cap incidences do not form one closed region cycle around edge ("
                << edge.first << ", " << edge.second << "), stable IDs ("
                << surface.vertices[edge.first].id << ", "
                << surface.vertices[edge.second].id << "), positions [("
                << surface.vertices[edge.first].referencePositionMeters.x
                << ", "
                << surface.vertices[edge.first].referencePositionMeters.y
                << ", "
                << surface.vertices[edge.first].referencePositionMeters.z
                << "), ("
                << surface.vertices[edge.second].referencePositionMeters.x
                << ", "
                << surface.vertices[edge.second].referencePositionMeters.y
                << ", "
                << surface.vertices[edge.second].referencePositionMeters.z
                << ")]:";
        constexpr std::size_t maximumReportedIncidences = 8;
        for (std::size_t index = 0;
             index < std::min(incidences.size(),
                              maximumReportedIncidences);
             ++index) {
            const EdgeIncidence& incidence = incidences[index];
            message << " [" << incidence.from << "->"
                    << incidence.to << ", regions "
                    << incidence.negativeRegion << "->"
                    << incidence.positiveRegion << ']';
        }
        if (incidences.size() > maximumReportedIncidences) {
            message << " ...";
        }
        throw std::invalid_argument(message.str());
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
