#include "scene_fluid_capped_face_partition.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t partitionIdentityDomain = 0x6361706661636570ULL;

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
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

using FaceKey = std::tuple<fluid::GridFaceAxis, std::size_t,
                           std::size_t, std::size_t>;

FaceKey faceKey(const fluid::GridFaceAxis axis,
                const std::size_t i,
                const std::size_t j,
                const std::size_t k) {
    return {axis, i, j, k};
}

fluid::GridFaceAxis gridFaceAxis(
    const SceneFluidOpeningPatchFaceAxis axis) {
    switch (axis) {
    case SceneFluidOpeningPatchFaceAxis::X:
        return fluid::GridFaceAxis::X;
    case SceneFluidOpeningPatchFaceAxis::Y:
        return fluid::GridFaceAxis::Y;
    case SceneFluidOpeningPatchFaceAxis::Z:
        return fluid::GridFaceAxis::Z;
    }
    throw std::invalid_argument(
        "capped face partition opening patch has invalid axis");
}

struct Point2 {
    double u = 0.0;
    double v = 0.0;
};

Point2 facePoint(const fluid::GridFaceAxis axis, const Vec3& point) {
    if (axis == fluid::GridFaceAxis::X) return {point.y, point.z};
    if (axis == fluid::GridFaceAxis::Y) return {point.z, point.x};
    return {point.x, point.y};
}

struct FaceBounds {
    double minimumU = 0.0;
    double maximumU = 0.0;
    double minimumV = 0.0;
    double maximumV = 0.0;
};

FaceBounds faceBounds(const fluid::PeriodicCartesianGrid& grid,
                      const fluid::GridFaceAxis axis,
                      const std::size_t i,
                      const std::size_t j,
                      const std::size_t k) {
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    if (axis == fluid::GridFaceAxis::X) {
        const double minimumU = lower.y + static_cast<double>(j) * spacing.y;
        const double minimumV = lower.z + static_cast<double>(k) * spacing.z;
        return {minimumU, minimumU + spacing.y,
                minimumV, minimumV + spacing.z};
    }
    if (axis == fluid::GridFaceAxis::Y) {
        const double minimumU = lower.z + static_cast<double>(k) * spacing.z;
        const double minimumV = lower.x + static_cast<double>(i) * spacing.x;
        return {minimumU, minimumU + spacing.z,
                minimumV, minimumV + spacing.x};
    }
    const double minimumU = lower.x + static_cast<double>(i) * spacing.x;
    const double minimumV = lower.y + static_cast<double>(j) * spacing.y;
    return {minimumU, minimumU + spacing.x,
            minimumV, minimumV + spacing.y};
}

double cross2(const Point2& first,
              const Point2& second,
              const Point2& third) {
    return (second.u - first.u) * (third.v - first.v)
        - (second.v - first.v) * (third.u - first.u);
}

bool near(const double first,
          const double second,
          const double tolerance) {
    return std::isfinite(first) && std::isfinite(second)
        && std::abs(first - second) <= tolerance;
}

bool samePoint(const Point2& first,
               const Point2& second,
               const double tolerance) {
    return std::hypot(first.u - second.u, first.v - second.v)
        <= tolerance;
}

bool onSegment(const Point2& first,
               const Point2& second,
               const Point2& point,
               const double tolerance) {
    const double length = std::hypot(
        second.u - first.u, second.v - first.v);
    if (std::abs(cross2(first, second, point))
        > tolerance * std::max(tolerance, length)) {
        return false;
    }
    return point.u >= std::min(first.u, second.u) - tolerance
        && point.u <= std::max(first.u, second.u) + tolerance
        && point.v >= std::min(first.v, second.v) - tolerance
        && point.v <= std::max(first.v, second.v) + tolerance;
}

int orientation(const Point2& first,
                const Point2& second,
                const Point2& point,
                const double tolerance) {
    const double value = cross2(first, second, point);
    const double length = std::hypot(
        second.u - first.u, second.v - first.v);
    const double bound = tolerance * std::max(tolerance, length);
    return value > bound ? 1 : (value < -bound ? -1 : 0);
}

bool segmentsIntersect(const Point2& a,
                       const Point2& b,
                       const Point2& c,
                       const Point2& d,
                       const double tolerance) {
    const int ac = orientation(a, b, c, tolerance);
    const int ad = orientation(a, b, d, tolerance);
    const int ca = orientation(c, d, a, tolerance);
    const int cb = orientation(c, d, b, tolerance);
    if (ac != ad && ca != cb
        && ac != 0 && ad != 0 && ca != 0 && cb != 0) {
        return true;
    }
    return (ac == 0 && onSegment(a, b, c, tolerance))
        || (ad == 0 && onSegment(a, b, d, tolerance))
        || (ca == 0 && onSegment(c, d, a, tolerance))
        || (cb == 0 && onSegment(c, d, b, tolerance));
}

bool pointInside(const Point2& point,
                 const std::vector<Point2>& polygon,
                 const double tolerance) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1;
         i < polygon.size(); j = i++) {
        if (onSegment(polygon[j], polygon[i], point, tolerance)) {
            throw std::invalid_argument(
                "capped face partition components touch");
        }
        const bool straddles = (polygon[i].v > point.v)
            != (polygon[j].v > point.v);
        if (straddles) {
            const double crossingU = polygon[i].u
                + (point.v - polygon[i].v)
                    * (polygon[j].u - polygon[i].u)
                    / (polygon[j].v - polygon[i].v);
            if (crossingU > point.u) inside = !inside;
        }
    }
    return inside;
}

bool sharesEdge(const std::vector<std::size_t>& first,
                const std::vector<std::size_t>& second) {
    std::size_t a = 0;
    std::size_t b = 0;
    while (a < first.size() && b < second.size()) {
        if (first[a] == second[b]) return true;
        if (first[a] < second[b]) ++a;
        else ++b;
    }
    return false;
}

double boundaryParameter(const Point2& point,
                         const FaceBounds& bounds,
                         const double tolerance) {
    const double width = bounds.maximumU - bounds.minimumU;
    const double height = bounds.maximumV - bounds.minimumV;
    if (near(point.v, bounds.minimumV, tolerance)) {
        return point.u - bounds.minimumU;
    }
    if (near(point.u, bounds.maximumU, tolerance)) {
        return width + point.v - bounds.minimumV;
    }
    if (near(point.v, bounds.maximumV, tolerance)) {
        return width + height + bounds.maximumU - point.u;
    }
    if (near(point.u, bounds.minimumU, tolerance)) {
        return 2.0 * width + height + bounds.maximumV - point.v;
    }
    throw std::invalid_argument(
        "capped face partition boundary point has no owner");
}

Point2 normalizedPoint(const Point2& point,
                       const FaceBounds& bounds,
                       const double tolerance) {
    if (!std::isfinite(point.u) || !std::isfinite(point.v)
        || point.u < bounds.minimumU - tolerance
        || point.u > bounds.maximumU + tolerance
        || point.v < bounds.minimumV - tolerance
        || point.v > bounds.maximumV + tolerance) {
        throw std::invalid_argument(
            "capped face partition point lies outside its face");
    }
    Point2 result{std::clamp(point.u, bounds.minimumU, bounds.maximumU),
                  std::clamp(point.v, bounds.minimumV, bounds.maximumV)};
    if (near(result.u, bounds.minimumU, tolerance)) result.u = bounds.minimumU;
    if (near(result.u, bounds.maximumU, tolerance)) result.u = bounds.maximumU;
    if (near(result.v, bounds.minimumV, tolerance)) result.v = bounds.minimumV;
    if (near(result.v, bounds.maximumV, tolerance)) result.v = bounds.maximumV;
    return result;
}

bool isBoundaryPoint(const Point2& point,
                     const FaceBounds& bounds) {
    return point.u == bounds.minimumU || point.u == bounds.maximumU
        || point.v == bounds.minimumV || point.v == bounds.maximumV;
}

struct InputSegment {
    Point2 first;
    Point2 second;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
};

struct ArrangementNode {
    Point2 point;
    bool onBoundary = false;
};

struct ArrangementEdge {
    std::size_t firstNode = 0;
    std::size_t secondNode = 0;
    bool source = false;
    std::size_t directedFromNode = 0;
    std::size_t directedToNode = 0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
};

struct ArrangementHalfEdge {
    std::size_t edgeIndex = 0;
    std::size_t fromNode = 0;
    std::size_t toNode = 0;
    double angleRadians = 0.0;
};

struct ArrangementCycle {
    double signedAreaSquareMeters = 0.0;
    StableId regionId = invalidStableId;
    std::vector<std::size_t> edgeIndices;
    std::vector<Point2> polygon;
};

void consumePairTest(std::size_t& count,
                     const SceneFluidCappedFacePartitionLimits& limits) {
    if (count == limits.maximumSegmentPairTests) {
        throw std::length_error(
            "capped face partition exceeds its pair-test limit");
    }
    ++count;
}

std::map<StableId, double> arrangementAreas(
    const FaceBounds& bounds,
    const std::vector<InputSegment>& segments,
    const SceneFluidCappedFacePartitionSettings& settings,
    const SceneFluidCappedFacePartitionLimits& limits,
    std::size_t& pairTestCount) {
    const double width = bounds.maximumU - bounds.minimumU;
    const double height = bounds.maximumV - bounds.minimumV;
    const double perimeter = 2.0 * (width + height);
    const double faceArea = width * height;
    if (!(width > 0.0) || !(height > 0.0)
        || !std::isfinite(faceArea) || segments.empty()) {
        throw std::invalid_argument(
            "capped face partition has invalid bounds or no segments");
    }

    std::vector<ArrangementNode> nodes;
    std::vector<ArrangementEdge> edges;
    std::set<std::pair<std::size_t, std::size_t>> edgeKeys;
    auto nodeFor = [&](const Point2& source) {
        const Point2 point = normalizedPoint(
            source, bounds, settings.geometryToleranceMeters);
        std::size_t match = nodes.size();
        for (std::size_t index = 0; index < nodes.size(); ++index) {
            consumePairTest(pairTestCount, limits);
            if (samePoint(nodes[index].point, point,
                          settings.geometryToleranceMeters)) {
                if (match != nodes.size()) {
                    throw std::invalid_argument(
                        "capped face partition endpoint matches multiple nodes");
                }
                match = index;
            }
        }
        if (match != nodes.size()) return match;
        nodes.push_back({point, isBoundaryPoint(point, bounds)});
        return nodes.size() - 1;
    };

    for (const auto& segment : segments) {
        if (segment.negativeSideRegionId == invalidStableId
            || segment.positiveSideRegionId == invalidStableId
            || segment.negativeSideRegionId
                == segment.positiveSideRegionId) {
            throw std::invalid_argument(
                "capped face partition source regions are invalid");
        }
        const std::size_t from = nodeFor(segment.first);
        const std::size_t to = nodeFor(segment.second);
        if (from == to) {
            throw std::invalid_argument(
                "capped face partition source edge is degenerate");
        }
        const auto key = std::minmax(from, to);
        if (!edgeKeys.emplace(key.first, key.second).second) {
            throw std::invalid_argument(
                "capped face partition repeats a source edge");
        }
        ArrangementEdge edge;
        edge.firstNode = key.first;
        edge.secondNode = key.second;
        edge.source = true;
        edge.directedFromNode = from;
        edge.directedToNode = to;
        edge.negativeSideRegionId = segment.negativeSideRegionId;
        edge.positiveSideRegionId = segment.positiveSideRegionId;
        edges.push_back(edge);
    }

    const std::array<Point2, 4> corners{{
        {bounds.minimumU, bounds.minimumV},
        {bounds.maximumU, bounds.minimumV},
        {bounds.maximumU, bounds.maximumV},
        {bounds.minimumU, bounds.maximumV},
    }};
    for (const Point2& corner : corners) {
        static_cast<void>(nodeFor(corner));
    }

    std::vector<std::pair<double, std::size_t>> boundaryNodes;
    for (std::size_t nodeIndex = 0;
         nodeIndex < nodes.size(); ++nodeIndex) {
        if (!nodes[nodeIndex].onBoundary) continue;
        boundaryNodes.emplace_back(boundaryParameter(
            nodes[nodeIndex].point, bounds,
            settings.geometryToleranceMeters), nodeIndex);
    }
    std::ranges::sort(boundaryNodes);
    if (boundaryNodes.size() < 4) {
        throw std::invalid_argument(
            "capped face partition boundary is incomplete");
    }
    for (std::size_t index = 1; index < boundaryNodes.size(); ++index) {
        if (boundaryNodes[index].first - boundaryNodes[index - 1].first
            <= settings.geometryToleranceMeters) {
            throw std::invalid_argument(
                "capped face partition boundary nodes overlap");
        }
    }
    if (boundaryNodes.front().first + perimeter
            - boundaryNodes.back().first
        <= settings.geometryToleranceMeters) {
        throw std::invalid_argument(
            "capped face partition boundary nodes overlap");
    }
    for (std::size_t index = 0; index < boundaryNodes.size(); ++index) {
        const std::size_t first = boundaryNodes[index].second;
        const std::size_t second = boundaryNodes[
            (index + 1) % boundaryNodes.size()].second;
        const auto key = std::minmax(first, second);
        if (!edgeKeys.emplace(key.first, key.second).second) {
            throw std::invalid_argument(
                "capped face partition overlaps its rectangular boundary");
        }
        ArrangementEdge edge;
        edge.firstNode = key.first;
        edge.secondNode = key.second;
        edges.push_back(edge);
    }

    std::vector<std::size_t> sourceDegree(nodes.size(), 0);
    for (const auto& edge : edges) {
        if (!edge.source) continue;
        ++sourceDegree[edge.firstNode];
        ++sourceDegree[edge.secondNode];
    }
    for (std::size_t nodeIndex = 0;
         nodeIndex < nodes.size(); ++nodeIndex) {
        if ((sourceDegree[nodeIndex] == 1 && !nodes[nodeIndex].onBoundary)
            || (nodes[nodeIndex].onBoundary
                && sourceDegree[nodeIndex] == 0
                && std::ranges::none_of(
                    corners, [&](const Point2& corner) {
                        return nodes[nodeIndex].point.u == corner.u
                            && nodes[nodeIndex].point.v == corner.v;
                    }))) {
            throw std::invalid_argument(
                "capped face partition has an unpaired source endpoint");
        }
    }

    for (std::size_t first = 0; first < edges.size(); ++first) {
        for (std::size_t second = first + 1;
             second < edges.size(); ++second) {
            const auto& a = edges[first];
            const auto& b = edges[second];
            if (a.firstNode == b.firstNode
                || a.firstNode == b.secondNode
                || a.secondNode == b.firstNode
                || a.secondNode == b.secondNode) {
                continue;
            }
            consumePairTest(pairTestCount, limits);
            if (segmentsIntersect(
                    nodes[a.firstNode].point, nodes[a.secondNode].point,
                    nodes[b.firstNode].point, nodes[b.secondNode].point,
                    settings.geometryToleranceMeters)) {
                throw std::invalid_argument(
                    "capped face partition contains an unstitched crossing");
            }
        }
    }

    std::vector<ArrangementHalfEdge> halfEdges;
    std::vector<std::vector<std::size_t>> outgoing(nodes.size());
    halfEdges.reserve(2 * edges.size());
    for (std::size_t edgeIndex = 0;
         edgeIndex < edges.size(); ++edgeIndex) {
        const auto& edge = edges[edgeIndex];
        for (const auto [from, to] : {
                 std::pair{edge.firstNode, edge.secondNode},
                 std::pair{edge.secondNode, edge.firstNode}}) {
            const Point2 delta{
                nodes[to].point.u - nodes[from].point.u,
                nodes[to].point.v - nodes[from].point.v,
            };
            if (!(std::hypot(delta.u, delta.v)
                    > settings.geometryToleranceMeters)) {
                throw std::invalid_argument(
                    "capped face partition has a short edge");
            }
            const std::size_t halfEdgeIndex = halfEdges.size();
            halfEdges.push_back({
                edgeIndex, from, to, std::atan2(delta.v, delta.u)});
            outgoing[from].push_back(halfEdgeIndex);
        }
    }
    const double angularTolerance = settings.geometryToleranceMeters
        / std::min(width, height);
    for (auto& nodeOutgoing : outgoing) {
        std::ranges::sort(nodeOutgoing, {},
                          [&](const std::size_t halfEdgeIndex) {
                              return halfEdges[halfEdgeIndex].angleRadians;
                          });
        for (std::size_t index = 0;
             index < nodeOutgoing.size(); ++index) {
            const double first = halfEdges[nodeOutgoing[index]].angleRadians;
            double second = halfEdges[nodeOutgoing[
                (index + 1) % nodeOutgoing.size()]].angleRadians;
            if (index + 1 == nodeOutgoing.size()) {
                second += 2.0 * std::acos(-1.0);
            }
            if (second - first <= angularTolerance) {
                throw std::invalid_argument(
                    "capped face partition has overlapping rays");
            }
        }
    }

    std::vector<bool> visited(halfEdges.size(), false);
    std::vector<ArrangementCycle> cycles;
    for (std::size_t start = 0;
         start < halfEdges.size(); ++start) {
        if (visited[start]) continue;
        ArrangementCycle cycle;
        std::size_t current = start;
        while (true) {
            if (cycle.edgeIndices.size() > halfEdges.size()
                || visited[current]) {
                throw std::invalid_argument(
                    "capped face partition cycle is invalid");
            }
            visited[current] = true;
            const auto& halfEdge = halfEdges[current];
            cycle.edgeIndices.push_back(halfEdge.edgeIndex);
            cycle.polygon.push_back(nodes[halfEdge.fromNode].point);
            const auto& edge = edges[halfEdge.edgeIndex];
            if (edge.source) {
                const StableId candidate =
                    halfEdge.fromNode == edge.directedFromNode
                        && halfEdge.toNode == edge.directedToNode
                    ? edge.positiveSideRegionId
                    : edge.negativeSideRegionId;
                if (cycle.regionId == invalidStableId) {
                    cycle.regionId = candidate;
                } else if (cycle.regionId != candidate) {
                    throw std::invalid_argument(
                        "capped face partition authored winding conflicts");
                }
            }
            const std::size_t reverse = current ^ 1U;
            const auto& destinationOutgoing = outgoing[halfEdge.toNode];
            const auto found = std::ranges::find(
                destinationOutgoing, reverse);
            if (found == destinationOutgoing.end()) {
                throw std::logic_error(
                    "capped face partition reverse edge is missing");
            }
            const std::size_t position = static_cast<std::size_t>(
                found - destinationOutgoing.begin());
            current = destinationOutgoing[
                position == 0
                    ? destinationOutgoing.size() - 1
                    : position - 1];
            if (current == start) break;
        }
        double twiceSignedArea = 0.0;
        for (std::size_t index = 0;
             index < cycle.polygon.size(); ++index) {
            const auto& first = cycle.polygon[index];
            const auto& second = cycle.polygon[
                (index + 1) % cycle.polygon.size()];
            twiceSignedArea +=
                first.u * second.v - second.u * first.v;
        }
        cycle.signedAreaSquareMeters = 0.5 * twiceSignedArea;
        if (!std::isfinite(cycle.signedAreaSquareMeters)
            || cycle.signedAreaSquareMeters == 0.0) {
            throw std::invalid_argument(
                "capped face partition cycle area is invalid");
        }
        std::ranges::sort(cycle.edgeIndices);
        cycles.push_back(std::move(cycle));
    }

    std::map<StableId, double> areas;
    std::vector<std::size_t> unlabeledPositiveCycles;
    std::size_t unlabeledNegativeCycleCount = 0;
    for (std::size_t index = 0; index < cycles.size(); ++index) {
        const auto& cycle = cycles[index];
        if (cycle.regionId != invalidStableId) {
            areas[cycle.regionId] += cycle.signedAreaSquareMeters;
        } else if (cycle.signedAreaSquareMeters > 0.0) {
            unlabeledPositiveCycles.push_back(index);
        } else {
            ++unlabeledNegativeCycleCount;
        }
    }
    if (unlabeledNegativeCycleCount != 1
        || unlabeledPositiveCycles.size() > 1) {
        throw std::invalid_argument(
            "capped face partition has ambiguous unlabeled cycles");
    }
    if (!unlabeledPositiveCycles.empty()) {
        const auto& base = cycles[unlabeledPositiveCycles.front()];
        if (std::abs(base.signedAreaSquareMeters - faceArea)
            > settings.areaClosureToleranceSquareMeters) {
            throw std::invalid_argument(
                "capped face partition unlabeled base is not the full face");
        }
        StableId rootRegionId = invalidStableId;
        for (const auto& candidate : cycles) {
            if (candidate.signedAreaSquareMeters >= 0.0
                || candidate.regionId == invalidStableId) {
                continue;
            }
            bool nested = false;
            for (const auto& container : cycles) {
                if (container.signedAreaSquareMeters <= 0.0
                    || container.regionId == invalidStableId) {
                    continue;
                }
                if (sharesEdge(candidate.edgeIndices,
                               container.edgeIndices)) {
                    continue;
                }
                consumePairTest(pairTestCount, limits);
                if (pointInside(candidate.polygon.front(),
                                container.polygon,
                                settings.geometryToleranceMeters)) {
                    nested = true;
                    break;
                }
            }
            if (!nested) {
                if (rootRegionId == invalidStableId) {
                    rootRegionId = candidate.regionId;
                } else if (rootRegionId != candidate.regionId) {
                    throw std::invalid_argument(
                        "capped face partition has conflicting root regions");
                }
            }
        }
        if (rootRegionId == invalidStableId) {
            throw std::invalid_argument(
                "capped face partition has no root region");
        }
        areas[rootRegionId] += faceArea;
    }

    double assignedArea = 0.0;
    for (auto iterator = areas.begin(); iterator != areas.end();) {
        if (!(iterator->second > 0.0)
            || !std::isfinite(iterator->second)) {
            if (std::isfinite(iterator->second)
                && std::abs(iterator->second)
                    <= settings.areaClosureToleranceSquareMeters) {
                iterator = areas.erase(iterator);
                continue;
            }
            throw std::invalid_argument(
                "capped face partition region area is invalid");
        }
        assignedArea += iterator->second;
        ++iterator;
    }
    if (areas.size() < 2 || !std::isfinite(assignedArea)
        || std::abs(assignedArea - faceArea)
            > settings.areaClosureToleranceSquareMeters) {
        throw std::invalid_argument(
            "capped face partition does not close exact region area");
    }
    return areas;
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double dot(const Point2& first, const Point2& second) {
    return first.u * second.u + first.v * second.v;
}

Vec3 faceNormal(const fluid::GridFaceAxis axis) {
    if (axis == fluid::GridFaceAxis::X) return {1.0, 0.0, 0.0};
    if (axis == fluid::GridFaceAxis::Y) return {0.0, 1.0, 0.0};
    return {0.0, 0.0, 1.0};
}

InputSegment openingSegment(
    const SceneFluidOpeningFaceCrossing& crossing) {
    Point2 first = facePoint(crossing.axis, crossing.first.positionMeters);
    Point2 second = facePoint(crossing.axis, crossing.second.positionMeters);
    const Vec3 tangent3 = cross(
        crossing.negativeToPositiveDirectionInFace,
        faceNormal(crossing.axis));
    const Point2 preferred = facePoint(crossing.axis, tangent3);
    const Point2 actual{second.u - first.u, second.v - first.v};
    const double alignment = dot(actual, preferred);
    if (!std::isfinite(alignment) || alignment == 0.0) {
        throw std::invalid_argument(
            "capped face partition opening crossing has no direction");
    }
    if (alignment < 0.0) std::swap(first, second);
    return {first, second,
            crossing.negativeSideRegionId,
            crossing.positiveSideRegionId};
}

std::uint64_t partitionStableId(
    const SceneFluidCappedFacePartition& partition,
    const std::vector<std::size_t>& materialReferences,
    const std::vector<std::size_t>& openingReferences,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings) {
    Fingerprint fingerprint;
    fingerprint.integer(partitionIdentityDomain);
    fingerprint.enumeration(partition.axis);
    fingerprint.integer(static_cast<std::uint64_t>(partition.i));
    fingerprint.integer(static_cast<std::uint64_t>(partition.j));
    fingerprint.integer(static_cast<std::uint64_t>(partition.k));
    fingerprint.integer(static_cast<std::uint64_t>(
        partition.materialChainReferenceCount));
    for (std::size_t offset = 0;
         offset < partition.materialChainReferenceCount; ++offset) {
        const std::size_t reference = materialReferences[
            partition.firstMaterialChainReference + offset];
        fingerprint.integer(epoch.faceChains.chains[reference].stableId);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        partition.openingCrossingReferenceCount));
    for (std::size_t offset = 0;
         offset < partition.openingCrossingReferenceCount; ++offset) {
        const std::size_t reference = openingReferences[
            partition.firstOpeningCrossingReference + offset];
        fingerprint.integer(openingCrossings.crossings[reference].stableId);
    }
    return fingerprint.value();
}

void validateSettings(
    const SceneFluidCappedFacePartitionSettings& settings) {
    if (!std::isfinite(settings.geometryToleranceMeters)
        || !(settings.geometryToleranceMeters > 0.0)
        || !std::isfinite(settings.areaClosureToleranceSquareMeters)
        || !(settings.areaClosureToleranceSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "capped face partition settings are invalid");
    }
}

template<typename Value>
void addCountBytes(const std::size_t count,
                   std::size_t& total) {
    std::size_t bytes = 0;
    std::size_t next = 0;
    if (!checkedMultiply(count, sizeof(Value), bytes)
        || !checkedAdd(total, bytes, next)) {
        throw std::length_error(
            "capped face partition storage size overflows");
    }
    total = next;
}

std::size_t storageBytesForCounts(
    const std::size_t faceCount,
    const std::size_t partitionCount,
    const std::size_t materialReferenceCount,
    const std::size_t openingReferenceCount,
    const std::size_t regionAreaCount) {
    std::size_t result = 0;
    addCountBytes<SceneFluidCappedFace>(faceCount, result);
    addCountBytes<SceneFluidCappedFacePartition>(partitionCount, result);
    addCountBytes<std::size_t>(materialReferenceCount, result);
    addCountBytes<std::size_t>(openingReferenceCount, result);
    addCountBytes<fluid::SceneFluidFaceRegionArea>(
        regionAreaCount, result);
    return result;
}

std::size_t storageBytes(
    const SceneFluidCappedFacePartitionSet& partitions) {
    return storageBytesForCounts(
        partitions.faces.size(), partitions.partitions.size(),
        partitions.materialChainReferences.size(),
        partitions.openingCrossingReferences.size(),
        partitions.regionAreas.size());
}

std::uint64_t partitionsFingerprint(
    const SceneFluidCappedFacePartitionSet& partitions) {
    Fingerprint fingerprint;
    fingerprint.integer(partitions.version);
    fingerprint.integer(partitions.surfaceDefinitionFingerprint);
    fingerprint.integer(partitions.surfaceStateFingerprint);
    fingerprint.integer(partitions.gridEpochFingerprint);
    fingerprint.integer(partitions.openingPatchFingerprint);
    fingerprint.integer(partitions.openingFaceCrossingFingerprint);
    fingerprint.integer(partitions.structureDefinitionFingerprint);
    fingerprint.integer(partitions.acceptedStepCount);
    fingerprint.real(partitions.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(partitions.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(partitions.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(partitions.cellCounts.z));
    for (const double value : {
             partitions.lowerMeters.x, partitions.lowerMeters.y,
             partitions.lowerMeters.z, partitions.upperMeters.x,
             partitions.upperMeters.y, partitions.upperMeters.z,
             partitions.settings.geometryToleranceMeters,
             partitions.settings.areaClosureToleranceSquareMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.touchedFaceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.unresolvedTouchedFaceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.segmentPairTestCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.faces.size()));
    for (const auto& face : partitions.faces) {
        fingerprint.enumeration(face.axis);
        fingerprint.integer(static_cast<std::uint64_t>(face.i));
        fingerprint.integer(static_cast<std::uint64_t>(face.j));
        fingerprint.integer(static_cast<std::uint64_t>(face.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            face.partitionIndex));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        partitions.partitions.size()));
    for (const auto& partition : partitions.partitions) {
        fingerprint.integer(partition.stableId);
        fingerprint.enumeration(partition.axis);
        fingerprint.integer(static_cast<std::uint64_t>(partition.i));
        fingerprint.integer(static_cast<std::uint64_t>(partition.j));
        fingerprint.integer(static_cast<std::uint64_t>(partition.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.activeFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.firstMaterialChainReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.materialChainReferenceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.firstOpeningCrossingReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.openingCrossingReferenceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.firstRegionArea));
        fingerprint.integer(static_cast<std::uint64_t>(
            partition.regionAreaCount));
        fingerprint.real(partition.faceAreaSquareMeters);
        fingerprint.real(partition.assignedAreaSquareMeters);
        fingerprint.real(partition.areaResidualSquareMeters);
    }
    for (const std::size_t reference :
         partitions.materialChainReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    for (const std::size_t reference :
         partitions.openingCrossingReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    for (const auto& area : partitions.regionAreas) {
        fingerprint.integer(area.regionId);
        fingerprint.real(area.areaSquareMeters);
    }
    return fingerprint.value();
}

bool sameGrid(const SceneFluidCappedFacePartitionSet& partitions,
              const fluid::PeriodicCartesianGrid& grid) {
    return partitions.cellCounts == grid.cellCounts()
        && partitions.lowerMeters == grid.lowerMeters()
        && partitions.upperMeters == grid.upperMeters();
}

SceneFluidCappedFacePartitionSet buildPartitions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings,
    const SceneFluidCappedFacePartitionSettings& settings,
    const SceneFluidCappedFacePartitionLimits& limits) {
    validateSettings(settings);
    SceneFluidCappedFacePartitionSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.gridEpochFingerprint = epoch.fingerprint;
    result.openingPatchFingerprint = openingPatches.fingerprint;
    result.openingFaceCrossingFingerprint = openingCrossings.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;

    std::map<FaceKey, std::vector<std::size_t>> crossingReferencesByFace;
    for (std::size_t crossingIndex = 0;
         crossingIndex < openingCrossings.crossings.size(); ++crossingIndex) {
        const auto& crossing = openingCrossings.crossings[crossingIndex];
        crossingReferencesByFace[faceKey(
            crossing.axis, crossing.i, crossing.j, crossing.k)]
            .push_back(crossingIndex);
    }
    result.touchedFaceCount = crossingReferencesByFace.size();
    if (result.touchedFaceCount > limits.maximumTouchedFaces) {
        throw std::length_error(
            "capped face partitions exceed their touched-face limit");
    }
    if (storageBytesForCounts(result.touchedFaceCount, 0, 0, 0, 0)
        > limits.maximumPartitionBytes) {
        throw std::length_error(
            "capped face partitions exceed their byte limit");
    }

    std::map<FaceKey, std::size_t> activeFaces;
    for (std::size_t activeIndex = 0;
         activeIndex < epoch.faceTopology.activeFaces.size(); ++activeIndex) {
        const auto& face = epoch.faceTopology.activeFaces[activeIndex];
        activeFaces.emplace(faceKey(face.axis, face.i, face.j, face.k),
                            activeIndex);
    }
    std::set<FaceKey> faceOwnedOpeningArea;
    for (const auto& patch : openingPatches.patches) {
        if (patch.ownerKind != SceneFluidOpeningPatchOwnerKind::Face) continue;
        faceOwnedOpeningArea.insert(faceKey(
            gridFaceAxis(patch.faceAxis),
            patch.faceI, patch.faceJ, patch.faceK));
    }

    std::size_t sourceSegmentCount = 0;
    std::unordered_set<std::uint64_t> partitionStableIds;
    for (const auto& [key, openingReferences] : crossingReferencesByFace) {
        const auto axis = std::get<0>(key);
        const std::size_t i = std::get<1>(key);
        const std::size_t j = std::get<2>(key);
        const std::size_t k = std::get<3>(key);
        SceneFluidCappedFace cappedFace;
        cappedFace.axis = axis;
        cappedFace.i = i;
        cappedFace.j = j;
        cappedFace.k = k;
        const auto active = activeFaces.find(key);
        const std::size_t activeFaceIndex = active == activeFaces.end()
            ? invalidSceneFluidActiveFaceIndex : active->second;
        if (faceOwnedOpeningArea.contains(key)
            || (activeFaceIndex != invalidSceneFluidActiveFaceIndex
                && epoch.faceTopology.activeFaces[activeFaceIndex]
                        .coplanarPatchReferenceCount != 0)) {
            ++result.unresolvedTouchedFaceCount;
            result.faces.push_back(cappedFace);
            continue;
        }

        std::vector<std::size_t> materialReferences;
        std::vector<InputSegment> segments;
        for (std::size_t chainIndex = 0;
             chainIndex < epoch.faceChains.chains.size(); ++chainIndex) {
            const auto& chain = epoch.faceChains.chains[chainIndex];
            if (chain.activeFaceIndex != activeFaceIndex
                || chain.negativeSideRegionId
                    == chain.positiveSideRegionId) {
                continue;
            }
            materialReferences.push_back(chainIndex);
            const bool closed = chain.kind
                == fluid::SceneFluidFaceChainKind::Closed;
            const std::size_t edgeCount = closed
                ? chain.nodeReferenceCount
                : chain.nodeReferenceCount - 1;
            for (std::size_t offset = 0; offset < edgeCount; ++offset) {
                const std::size_t firstNode = epoch.faceChains.nodeReferences[
                    chain.firstNodeReference + offset];
                const std::size_t secondNode = epoch.faceChains.nodeReferences[
                    chain.firstNodeReference
                        + ((offset + 1) % chain.nodeReferenceCount)];
                segments.push_back({
                    facePoint(axis, epoch.faceGraph.nodes[firstNode]
                                        .positionMeters),
                    facePoint(axis, epoch.faceGraph.nodes[secondNode]
                                        .positionMeters),
                    chain.negativeSideRegionId,
                    chain.positiveSideRegionId,
                });
            }
        }
        for (const std::size_t crossingIndex : openingReferences) {
            segments.push_back(openingSegment(
                openingCrossings.crossings[crossingIndex]));
        }
        std::size_t newSourceCount = 0;
        if (!checkedAdd(sourceSegmentCount, segments.size(), newSourceCount)
            || newSourceCount > limits.maximumReferences) {
            throw std::length_error(
                "capped face partition source work exceeds its limit");
        }
        sourceSegmentCount = newSourceCount;

        std::map<StableId, double> areas;
        try {
            areas = arrangementAreas(
                faceBounds(grid, axis, i, j, k), segments,
                settings, limits, result.segmentPairTestCount);
        } catch (const std::invalid_argument&) {
            ++result.unresolvedTouchedFaceCount;
            result.faces.push_back(cappedFace);
            continue;
        }
        if (result.partitions.size() == limits.maximumPartitions) {
            throw std::length_error(
                "capped face partitions exceed their count limit");
        }
        std::size_t newReferenceCount = 0;
        if (!checkedAdd(result.materialChainReferences.size(),
                        materialReferences.size(), newReferenceCount)
            || !checkedAdd(newReferenceCount,
                           result.openingCrossingReferences.size(),
                           newReferenceCount)
            || !checkedAdd(newReferenceCount, openingReferences.size(),
                           newReferenceCount)
            || !checkedAdd(newReferenceCount, result.regionAreas.size(),
                           newReferenceCount)
            || !checkedAdd(newReferenceCount, areas.size(),
                           newReferenceCount)
            || newReferenceCount > limits.maximumReferences) {
            throw std::length_error(
                "capped face partition references exceed their limit");
        }
        std::size_t prospectivePartitionCount = 0;
        std::size_t prospectiveMaterialReferenceCount = 0;
        std::size_t prospectiveOpeningReferenceCount = 0;
        std::size_t prospectiveRegionAreaCount = 0;
        if (!checkedAdd(result.partitions.size(), std::size_t{1},
                        prospectivePartitionCount)
            || !checkedAdd(result.materialChainReferences.size(),
                           materialReferences.size(),
                           prospectiveMaterialReferenceCount)
            || !checkedAdd(result.openingCrossingReferences.size(),
                           openingReferences.size(),
                           prospectiveOpeningReferenceCount)
            || !checkedAdd(result.regionAreas.size(), areas.size(),
                           prospectiveRegionAreaCount)) {
            throw std::length_error(
                "capped face partition retained counts overflow");
        }
        if (storageBytesForCounts(
                result.touchedFaceCount, prospectivePartitionCount,
                prospectiveMaterialReferenceCount,
                prospectiveOpeningReferenceCount,
                prospectiveRegionAreaCount)
            > limits.maximumPartitionBytes) {
            throw std::length_error(
                "capped face partitions exceed their byte limit");
        }

        SceneFluidCappedFacePartition partition;
        partition.axis = axis;
        partition.i = i;
        partition.j = j;
        partition.k = k;
        partition.activeFaceIndex = activeFaceIndex;
        partition.firstMaterialChainReference =
            result.materialChainReferences.size();
        partition.materialChainReferenceCount = materialReferences.size();
        result.materialChainReferences.insert(
            result.materialChainReferences.end(),
            materialReferences.begin(), materialReferences.end());
        partition.firstOpeningCrossingReference =
            result.openingCrossingReferences.size();
        partition.openingCrossingReferenceCount = openingReferences.size();
        result.openingCrossingReferences.insert(
            result.openingCrossingReferences.end(),
            openingReferences.begin(), openingReferences.end());
        partition.firstRegionArea = result.regionAreas.size();
        partition.regionAreaCount = areas.size();
        const auto bounds = faceBounds(grid, axis, i, j, k);
        partition.faceAreaSquareMeters =
            (bounds.maximumU - bounds.minimumU)
            * (bounds.maximumV - bounds.minimumV);
        for (const auto [regionId, area] : areas) {
            result.regionAreas.push_back({regionId, area});
            partition.assignedAreaSquareMeters += area;
        }
        partition.areaResidualSquareMeters =
            partition.assignedAreaSquareMeters
            - partition.faceAreaSquareMeters;
        partition.stableId = partitionStableId(
            partition, result.materialChainReferences,
            result.openingCrossingReferences, epoch, openingCrossings);
        if (!partitionStableIds.insert(partition.stableId).second) {
            throw std::invalid_argument(
                "capped face partition stable ID collides");
        }
        cappedFace.partitionIndex = result.partitions.size();
        result.partitions.push_back(partition);
        result.faces.push_back(cappedFace);
    }

    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumPartitionBytes) {
        throw std::length_error(
            "capped face partitions exceed their byte limit");
    }
    result.fingerprint = partitionsFingerprint(result);
    return result;
}

} // namespace

SceneFluidCappedFacePartitionSet buildSceneFluidCappedFacePartitions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings,
    const SceneFluidCappedFacePartitionSettings& settings,
    const SceneFluidCappedFacePartitionLimits& limits) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidOpeningFaceCrossings(
        openingCrossings, surface, state, caps, quadrature,
        openingPatches, grid);
    auto result = buildPartitions(
        surface, state, grid, epoch, openingPatches,
        openingCrossings, settings, limits);
    validateSceneFluidCappedFacePartitions(
        result, surface, state, grid, transfer, epoch, caps, quadrature,
        openingPatches, openingCrossings);
    return result;
}

void validateSceneFluidCappedFacePartitions(
    const SceneFluidCappedFacePartitionSet& partitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFaceCrossingSet& openingCrossings) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidOpeningFaceCrossings(
        openingCrossings, surface, state, caps, quadrature,
        openingPatches, grid);
    validateSettings(partitions.settings);
    if (partitions.version != sceneFluidCappedFacePartitionVersion
        || partitions.fingerprint == 0
        || partitions.surfaceDefinitionFingerprint != surface.fingerprint
        || partitions.surfaceStateFingerprint != state.fingerprint
        || partitions.gridEpochFingerprint != epoch.fingerprint
        || partitions.openingPatchFingerprint != openingPatches.fingerprint
        || partitions.openingFaceCrossingFingerprint
            != openingCrossings.fingerprint
        || partitions.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || partitions.acceptedStepCount != state.acceptedStepCount
        || partitions.simulationTimeSeconds != state.simulationTimeSeconds
        || !sameGrid(partitions, grid)) {
        throw std::invalid_argument(
            "capped face partition identity is invalid");
    }
    const SceneFluidCappedFacePartitionLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildPartitions(
        surface, state, grid, epoch, openingPatches,
        openingCrossings, partitions.settings, unlimited);
    if (partitions != expected
        || partitions.ownedStorageBytes != storageBytes(partitions)
        || partitions.fingerprint != partitionsFingerprint(partitions)) {
        throw std::invalid_argument(
            "capped face partition payload is invalid");
    }
}

} // namespace simwing::fsi
