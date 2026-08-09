#include "fluid/scene_surface_face_partition.h"

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

struct Point2 { double u = 0.0; double v = 0.0; };

Point2 facePoint(const GridFaceAxis axis, const Vec3& point) {
    if (axis == GridFaceAxis::X) return {point.y, point.z};
    if (axis == GridFaceAxis::Y) return {point.z, point.x};
    return {point.x, point.y};
}

double cross2(const Point2& a, const Point2& b, const Point2& c) {
    return (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
}

bool onSegment(const Point2& a, const Point2& b, const Point2& p,
               const double tolerance) {
    const double length = std::hypot(b.u - a.u, b.v - a.v);
    if (std::abs(cross2(a, b, p))
        > tolerance * std::max(tolerance, length)) return false;
    return p.u >= std::min(a.u, b.u) - tolerance
        && p.u <= std::max(a.u, b.u) + tolerance
        && p.v >= std::min(a.v, b.v) - tolerance
        && p.v <= std::max(a.v, b.v) + tolerance;
}

int orientation(const Point2& a, const Point2& b, const Point2& c,
                const double tolerance) {
    const double value = cross2(a, b, c);
    const double length = std::hypot(b.u - a.u, b.v - a.v);
    const double bound = tolerance * std::max(tolerance, length);
    return value > bound ? 1 : (value < -bound ? -1 : 0);
}

bool segmentsIntersect(const Point2& a, const Point2& b,
                       const Point2& c, const Point2& d,
                       const double tolerance) {
    const int ac = orientation(a, b, c, tolerance);
    const int ad = orientation(a, b, d, tolerance);
    const int ca = orientation(c, d, a, tolerance);
    const int cb = orientation(c, d, b, tolerance);
    if (ac != ad && ca != cb && ac != 0 && ad != 0 && ca != 0 && cb != 0)
        return true;
    return (ac == 0 && onSegment(a, b, c, tolerance))
        || (ad == 0 && onSegment(a, b, d, tolerance))
        || (ca == 0 && onSegment(c, d, a, tolerance))
        || (cb == 0 && onSegment(c, d, b, tolerance));
}

bool pointInside(const Point2& point, const std::vector<Point2>& polygon,
                 const double tolerance) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1;
         i < polygon.size(); j = i++) {
        if (onSegment(polygon[j], polygon[i], point, tolerance)) {
            throw std::invalid_argument("scene fluid face loops touch");
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

double faceArea(const PeriodicCartesianGrid& grid, const GridFaceAxis axis) {
    const Vector3 h = grid.cellSpacingMeters();
    if (axis == GridFaceAxis::X) return h.y * h.z;
    if (axis == GridFaceAxis::Y) return h.x * h.z;
    return h.x * h.y;
}

struct FaceBounds {
    double minimumU = 0.0;
    double maximumU = 0.0;
    double minimumV = 0.0;
    double maximumV = 0.0;
};

struct AreaMoment2 {
    double areaSquareMeters = 0.0;
    double firstMomentUMeters3 = 0.0;
    double firstMomentVMeters3 = 0.0;
};

AreaMoment2 polygonAreaMoment(const std::vector<Point2>& polygon) {
    AreaMoment2 result;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1) % polygon.size()];
        const double cross = first.u * second.v - second.u * first.v;
        result.areaSquareMeters += 0.5 * cross;
        result.firstMomentUMeters3 +=
            (first.u + second.u) * cross / 6.0;
        result.firstMomentVMeters3 +=
            (first.v + second.v) * cross / 6.0;
    }
    return result;
}

AreaMoment2 rectangleAreaMoment(const FaceBounds& bounds) {
    const double area = (bounds.maximumU - bounds.minimumU)
        * (bounds.maximumV - bounds.minimumV);
    return {
        area,
        area * 0.5 * (bounds.minimumU + bounds.maximumU),
        area * 0.5 * (bounds.minimumV + bounds.maximumV),
    };
}

double facePlaneCoordinate(const PeriodicCartesianGrid& grid,
                           const SceneFluidActiveFace& face) {
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    if (face.axis == GridFaceAxis::X) {
        return lower.x + static_cast<double>(face.i) * spacing.x;
    }
    if (face.axis == GridFaceAxis::Y) {
        return lower.y + static_cast<double>(face.j) * spacing.y;
    }
    return lower.z + static_cast<double>(face.k) * spacing.z;
}

SceneFluidFaceRegionArea makeRegionArea(
    const StableId regionId,
    const GridFaceAxis axis,
    const double planeCoordinateMeters,
    const AreaMoment2& moment) {
    SceneFluidFaceRegionArea result;
    result.regionId = regionId;
    result.areaSquareMeters = moment.areaSquareMeters;
    if (axis == GridFaceAxis::X) {
        result.firstMomentMeters3 = {
            planeCoordinateMeters * moment.areaSquareMeters,
            moment.firstMomentUMeters3,
            moment.firstMomentVMeters3,
        };
    } else if (axis == GridFaceAxis::Y) {
        result.firstMomentMeters3 = {
            moment.firstMomentVMeters3,
            planeCoordinateMeters * moment.areaSquareMeters,
            moment.firstMomentUMeters3,
        };
    } else {
        result.firstMomentMeters3 = {
            moment.firstMomentUMeters3,
            moment.firstMomentVMeters3,
            planeCoordinateMeters * moment.areaSquareMeters,
        };
    }
    if (moment.areaSquareMeters > 0.0) {
        result.centroidMeters = {
            result.firstMomentMeters3.x / moment.areaSquareMeters,
            result.firstMomentMeters3.y / moment.areaSquareMeters,
            result.firstMomentMeters3.z / moment.areaSquareMeters,
        };
    }
    return result;
}

FaceBounds faceBounds(const PeriodicCartesianGrid& grid,
                      const SceneFluidActiveFace& face) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    if (face.axis == GridFaceAxis::X) {
        const double minimumU = lower.y
            + static_cast<double>(face.j) * spacing.y;
        const double minimumV = lower.z
            + static_cast<double>(face.k) * spacing.z;
        return {minimumU, minimumU + spacing.y,
                minimumV, minimumV + spacing.z};
    }
    if (face.axis == GridFaceAxis::Y) {
        const double minimumU = lower.z
            + static_cast<double>(face.k) * spacing.z;
        const double minimumV = lower.x
            + static_cast<double>(face.i) * spacing.x;
        return {minimumU, minimumU + spacing.z,
                minimumV, minimumV + spacing.x};
    }
    const double minimumU = lower.x
        + static_cast<double>(face.i) * spacing.x;
    const double minimumV = lower.y
        + static_cast<double>(face.j) * spacing.y;
    return {minimumU, minimumU + spacing.x,
            minimumV, minimumV + spacing.y};
}

bool near(const double first, const double second,
          const double tolerance) {
    return std::abs(first - second) <= tolerance;
}

double boundaryParameter(const Point2& point,
                         const FaceBounds& bounds,
                         const double tolerance) {
    const double width = bounds.maximumU - bounds.minimumU;
    const double height = bounds.maximumV - bounds.minimumV;
    if (point.u < bounds.minimumU - tolerance
        || point.u > bounds.maximumU + tolerance
        || point.v < bounds.minimumV - tolerance
        || point.v > bounds.maximumV + tolerance) {
        throw std::invalid_argument(
            "scene fluid boundary-chain point leaves its face");
    }
    const bool minimumU = near(point.u, bounds.minimumU, tolerance);
    const bool maximumU = near(point.u, bounds.maximumU, tolerance);
    const bool minimumV = near(point.v, bounds.minimumV, tolerance);
    const bool maximumV = near(point.v, bounds.maximumV, tolerance);
    if (minimumU && minimumV) return 0.0;
    if (maximumU && minimumV) return width;
    if (maximumU && maximumV) return width + height;
    if (minimumU && maximumV) return 2.0 * width + height;
    if (minimumV) return point.u - bounds.minimumU;
    if (maximumU) return width + point.v - bounds.minimumV;
    if (maximumV) return width + height + bounds.maximumU - point.u;
    if (minimumU)
        return 2.0 * width + height + bounds.maximumV - point.v;
    throw std::invalid_argument(
        "scene fluid boundary-chain endpoint is not on its face boundary");
}

AreaMoment2 boundaryChainPositiveAreaMoment(
    const SceneFluidActiveFace& face,
    const SceneFluidFaceChain& chain,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceGraph& graph,
    const PeriodicCartesianGrid& grid,
    const SceneFluidFacePartitionSettings& settings,
    const SceneFluidFacePartitionLimits& limits,
    std::size_t& segmentPairTestCount) {
    if (chain.kind != SceneFluidFaceChainKind::Open
        || chain.nodeReferenceCount < 2
        || chain.segmentReferenceCount + 1 != chain.nodeReferenceCount
        || chain.endpointFaceBoundaryMasks[0] == FaceBoundaryNone
        || chain.endpointFaceBoundaryMasks[1] == FaceBoundaryNone
        || chain.endpointOnAuthoredOpening[0]
        || chain.endpointOnAuthoredOpening[1]) {
        throw std::invalid_argument(
            "scene fluid boundary face chain is invalid");
    }
    const FaceBounds bounds = faceBounds(grid, face);
    const double width = bounds.maximumU - bounds.minimumU;
    const double height = bounds.maximumV - bounds.minimumV;
    const double perimeter = 2.0 * (width + height);
    std::vector<Point2> polygon;
    polygon.reserve(chain.nodeReferenceCount + 4);
    for (std::size_t offset = 0;
         offset < chain.nodeReferenceCount; ++offset) {
        const auto& node = graph.nodes[chains.nodeReferences[
            chain.firstNodeReference + offset]];
        if (node.activeFaceIndex != chain.activeFaceIndex
            || !std::isfinite(node.positionMeters.x)
            || !std::isfinite(node.positionMeters.y)
            || !std::isfinite(node.positionMeters.z)
            || (offset != 0 && offset + 1 != chain.nodeReferenceCount
                && (node.faceBoundaryMask != FaceBoundaryNone
                    || node.authoredOpeningBoundary))) {
            throw std::invalid_argument(
                "scene fluid boundary face chain has invalid nodes");
        }
        const Point2 point = facePoint(face.axis, node.positionMeters);
        if (point.u < bounds.minimumU - settings.geometryToleranceMeters
            || point.u > bounds.maximumU + settings.geometryToleranceMeters
            || point.v < bounds.minimumV - settings.geometryToleranceMeters
            || point.v > bounds.maximumV + settings.geometryToleranceMeters) {
            throw std::invalid_argument(
                "scene fluid boundary face-chain node leaves its face");
        }
        polygon.push_back(point);
    }
    const double startParameter = boundaryParameter(
        polygon.front(), bounds, settings.geometryToleranceMeters);
    const double endParameter = boundaryParameter(
        polygon.back(), bounds, settings.geometryToleranceMeters);
    double targetParameter = startParameter;
    if (targetParameter <= endParameter
        + settings.geometryToleranceMeters) {
        targetParameter += perimeter;
    }
    const std::array<std::pair<double, Point2>, 4> corners{{
        {0.0, {bounds.minimumU, bounds.minimumV}},
        {width, {bounds.maximumU, bounds.minimumV}},
        {width + height, {bounds.maximumU, bounds.maximumV}},
        {2.0 * width + height, {bounds.minimumU, bounds.maximumV}},
    }};
    std::vector<std::pair<double, Point2>> traversedCorners;
    traversedCorners.reserve(corners.size());
    for (const auto& [baseParameter, point] : corners) {
        double parameter = baseParameter;
        while (parameter <= endParameter
               + settings.geometryToleranceMeters) {
            parameter += perimeter;
        }
        if (parameter < targetParameter
            - settings.geometryToleranceMeters) {
            traversedCorners.emplace_back(parameter, point);
        }
    }
    std::ranges::sort(traversedCorners, {}, &std::pair<double, Point2>::first);
    for (const auto& [parameter, point] : traversedCorners) {
        static_cast<void>(parameter);
        polygon.push_back(point);
    }
    for (std::size_t first = 0; first < polygon.size(); ++first) {
        const std::size_t firstNext = (first + 1) % polygon.size();
        for (std::size_t second = first + 1;
             second < polygon.size(); ++second) {
            const std::size_t secondNext = (second + 1) % polygon.size();
            if (firstNext == second || secondNext == first) continue;
            if (segmentPairTestCount == limits.maximumSegmentPairTests) {
                throw std::length_error(
                    "scene fluid face partition exceeds pair-test limit");
            }
            ++segmentPairTestCount;
            if (segmentsIntersect(
                    polygon[first], polygon[firstNext], polygon[second],
                    polygon[secondNext], settings.geometryToleranceMeters)) {
                throw std::invalid_argument(
                    "scene fluid boundary face chain self-intersects");
            }
        }
    }
    const auto positiveMoment = polygonAreaMoment(polygon);
    const double fullArea = width * height;
    if (!std::isfinite(positiveMoment.areaSquareMeters)
        || !std::isfinite(positiveMoment.firstMomentUMeters3)
        || !std::isfinite(positiveMoment.firstMomentVMeters3)
        || !(positiveMoment.areaSquareMeters > 0.0)
        || !(positiveMoment.areaSquareMeters < fullArea)) {
        throw std::invalid_argument(
            "scene fluid boundary face-chain area is invalid");
    }
    return positiveMoment;
}

struct ArrangementNode {
    Point2 point;
    std::size_t sourceGraphNode = std::numeric_limits<std::size_t>::max();
    bool onBoundary = false;
};

struct ArrangementEdge {
    std::size_t firstNode = 0;
    std::size_t secondNode = 0;
    bool material = false;
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

bool supportsRegionSeparatingArrangement(
    const std::size_t activeFaceIndex,
    const std::vector<std::size_t>& chainIndices,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceGraphRange& graphRange) {
    std::vector<std::size_t> degree(graphRange.nodeCount, 0);
    for (const std::size_t chainIndex : chainIndices) {
        if (chainIndex >= chains.chains.size()) return false;
        const auto& chain = chains.chains[chainIndex];
        if (chain.activeFaceIndex != activeFaceIndex
            || chain.kind != SceneFluidFaceChainKind::Open
            || chain.nodeReferenceCount < 2
            || chain.segmentReferenceCount + 1
                != chain.nodeReferenceCount
            || chain.negativeSideRegionId == invalidStableId
            || chain.positiveSideRegionId == invalidStableId
            || chain.negativeSideRegionId == chain.positiveSideRegionId
            || chain.endpointOnAuthoredOpening[0]
            || chain.endpointOnAuthoredOpening[1]) {
            return false;
        }
        for (std::size_t offset = 0;
             offset < chain.nodeReferenceCount; ++offset) {
            const std::size_t nodeIndex = chains.nodeReferences[
                chain.firstNodeReference + offset];
            if (nodeIndex < graphRange.firstNode
                || nodeIndex >= graphRange.firstNode + graphRange.nodeCount
                || graph.nodes[nodeIndex].authoredOpeningBoundary) {
                return false;
            }
            if (offset != 0) ++degree[nodeIndex - graphRange.firstNode];
            if (offset + 1 != chain.nodeReferenceCount) {
                ++degree[nodeIndex - graphRange.firstNode];
            }
        }
    }
    for (std::size_t offset = 0; offset < degree.size(); ++offset) {
        if (degree[offset] == 1
            && graph.nodes[graphRange.firstNode + offset].faceBoundaryMask
                == FaceBoundaryNone) {
            return false;
        }
    }
    return true;
}

std::map<StableId, AreaMoment2> boundaryChainArrangementAreaMoments(
    const std::size_t activeFaceIndex,
    const SceneFluidActiveFace& face,
    const std::vector<std::size_t>& chainIndices,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceGraphRange& graphRange,
    const PeriodicCartesianGrid& grid,
    const SceneFluidFacePartitionSettings& settings,
    const SceneFluidFacePartitionLimits& limits,
    std::size_t& segmentPairTestCount) {
    const FaceBounds bounds = faceBounds(grid, face);
    const double width = bounds.maximumU - bounds.minimumU;
    const double height = bounds.maximumV - bounds.minimumV;
    const double perimeter = 2.0 * (width + height);
    std::vector<ArrangementNode> nodes;
    std::vector<ArrangementEdge> edges;
    std::map<std::size_t, std::size_t> localNodeByGraphNode;
    std::set<std::pair<std::size_t, std::size_t>> edgeKeys;

    auto localNode = [&](const std::size_t graphNodeIndex) {
        const auto found = localNodeByGraphNode.find(graphNodeIndex);
        if (found != localNodeByGraphNode.end()) return found->second;
        if (graphNodeIndex < graphRange.firstNode
            || graphNodeIndex
                >= graphRange.firstNode + graphRange.nodeCount) {
            throw std::invalid_argument(
                "scene fluid face arrangement has a foreign node");
        }
        const auto& source = graph.nodes[graphNodeIndex];
        ArrangementNode node;
        node.point = facePoint(face.axis, source.positionMeters);
        node.sourceGraphNode = graphNodeIndex;
        node.onBoundary = source.faceBoundaryMask != FaceBoundaryNone;
        if (!std::isfinite(node.point.u) || !std::isfinite(node.point.v)
            || node.point.u
                < bounds.minimumU - settings.geometryToleranceMeters
            || node.point.u
                > bounds.maximumU + settings.geometryToleranceMeters
            || node.point.v
                < bounds.minimumV - settings.geometryToleranceMeters
            || node.point.v
                > bounds.maximumV + settings.geometryToleranceMeters
            || source.authoredOpeningBoundary) {
            throw std::invalid_argument(
                "scene fluid face arrangement has invalid source nodes");
        }
        const std::size_t index = nodes.size();
        nodes.push_back(node);
        localNodeByGraphNode.emplace(graphNodeIndex, index);
        return index;
    };

    for (const std::size_t chainIndex : chainIndices) {
        if (chainIndex >= chains.chains.size()) {
            throw std::invalid_argument(
                "scene fluid face arrangement has an invalid chain");
        }
        const auto& chain = chains.chains[chainIndex];
        if (chain.activeFaceIndex != activeFaceIndex
            || chain.kind != SceneFluidFaceChainKind::Open
            || chain.nodeReferenceCount < 2
            || chain.segmentReferenceCount + 1
                != chain.nodeReferenceCount
            || chain.negativeSideRegionId == invalidStableId
            || chain.positiveSideRegionId == invalidStableId
            || chain.negativeSideRegionId
                == chain.positiveSideRegionId
            || chain.endpointOnAuthoredOpening[0]
            || chain.endpointOnAuthoredOpening[1]) {
            throw std::invalid_argument(
                "scene fluid face arrangement chain is invalid");
        }
        for (std::size_t offset = 0;
             offset + 1 < chain.nodeReferenceCount; ++offset) {
            const std::size_t from = localNode(chains.nodeReferences[
                chain.firstNodeReference + offset]);
            const std::size_t to = localNode(chains.nodeReferences[
                chain.firstNodeReference + offset + 1]);
            if (from == to) {
                throw std::invalid_argument(
                    "scene fluid face arrangement edge is degenerate");
            }
            const auto key = std::minmax(from, to);
            if (!edgeKeys.emplace(key.first, key.second).second) {
                throw std::invalid_argument(
                    "scene fluid face arrangement repeats an edge");
            }
            ArrangementEdge edge;
            edge.firstNode = key.first;
            edge.secondNode = key.second;
            edge.material = true;
            edge.directedFromNode = from;
            edge.directedToNode = to;
            edge.negativeSideRegionId = chain.negativeSideRegionId;
            edge.positiveSideRegionId = chain.positiveSideRegionId;
            edges.push_back(edge);
        }
    }
    const std::array<Point2, 4> cornerPoints{{
        {bounds.minimumU, bounds.minimumV},
        {bounds.maximumU, bounds.minimumV},
        {bounds.maximumU, bounds.maximumV},
        {bounds.minimumU, bounds.maximumV},
    }};
    for (const Point2& corner : cornerPoints) {
        std::size_t matchingNode = nodes.size();
        for (std::size_t nodeIndex = 0;
             nodeIndex < nodes.size(); ++nodeIndex) {
            if (near(nodes[nodeIndex].point.u, corner.u,
                     settings.geometryToleranceMeters)
                && near(nodes[nodeIndex].point.v, corner.v,
                        settings.geometryToleranceMeters)) {
                if (matchingNode != nodes.size()) {
                    throw std::invalid_argument(
                        "scene fluid face arrangement duplicates a corner");
                }
                matchingNode = nodeIndex;
            }
        }
        if (matchingNode == nodes.size()) {
            nodes.push_back({corner,
                std::numeric_limits<std::size_t>::max(), true});
        } else if (!nodes[matchingNode].onBoundary) {
            throw std::invalid_argument(
                "scene fluid face arrangement corner ownership is invalid");
        }
    }

    std::vector<std::pair<double, std::size_t>> boundaryNodes;
    for (std::size_t nodeIndex = 0;
         nodeIndex < nodes.size(); ++nodeIndex) {
        if (!nodes[nodeIndex].onBoundary) continue;
        boundaryNodes.emplace_back(
            boundaryParameter(nodes[nodeIndex].point, bounds,
                              settings.geometryToleranceMeters),
            nodeIndex);
    }
    std::ranges::sort(boundaryNodes);
    if (boundaryNodes.size() < 4) {
        throw std::invalid_argument(
            "scene fluid face arrangement boundary is incomplete");
    }
    for (std::size_t index = 1;
         index < boundaryNodes.size(); ++index) {
        if (boundaryNodes[index].first
                - boundaryNodes[index - 1].first
            <= settings.geometryToleranceMeters) {
            throw std::invalid_argument(
                "scene fluid face arrangement boundary nodes overlap");
        }
    }
    if (boundaryNodes.front().first + perimeter
            - boundaryNodes.back().first
        <= settings.geometryToleranceMeters) {
        throw std::invalid_argument(
            "scene fluid face arrangement boundary nodes overlap");
    }
    for (std::size_t index = 0;
         index < boundaryNodes.size(); ++index) {
        const std::size_t first = boundaryNodes[index].second;
        const std::size_t second = boundaryNodes[
            (index + 1) % boundaryNodes.size()].second;
        const auto key = std::minmax(first, second);
        if (!edgeKeys.emplace(key.first, key.second).second) {
            throw std::invalid_argument(
                "scene fluid face arrangement overlaps its boundary");
        }
        ArrangementEdge edge;
        edge.firstNode = key.first;
        edge.secondNode = key.second;
        edges.push_back(edge);
    }

    std::vector<std::size_t> materialDegree(nodes.size(), 0);
    for (const auto& edge : edges) {
        if (!edge.material) continue;
        ++materialDegree[edge.firstNode];
        ++materialDegree[edge.secondNode];
    }
    for (std::size_t nodeIndex = 0;
         nodeIndex < nodes.size(); ++nodeIndex) {
        const auto& node = nodes[nodeIndex];
        if (node.sourceGraphNode
            == std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        if ((materialDegree[nodeIndex] == 1 && !node.onBoundary)
            || (node.onBoundary && materialDegree[nodeIndex] == 0)) {
            throw std::invalid_argument(
                "scene fluid face arrangement has incomplete node ownership");
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
            if (segmentPairTestCount == limits.maximumSegmentPairTests) {
                throw std::length_error(
                    "scene fluid face partition exceeds pair-test limit");
            }
            ++segmentPairTestCount;
            if (segmentsIntersect(
                    nodes[a.firstNode].point, nodes[a.secondNode].point,
                    nodes[b.firstNode].point, nodes[b.secondNode].point,
                    settings.geometryToleranceMeters)) {
                throw std::invalid_argument(
                    "scene fluid face arrangement edges intersect");
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
                    "scene fluid face arrangement has a short edge");
            }
            const std::size_t halfEdgeIndex = halfEdges.size();
            halfEdges.push_back({
                edgeIndex, from, to, std::atan2(delta.v, delta.u)});
            outgoing[from].push_back(halfEdgeIndex);
        }
    }
    const double angularTolerance = settings.geometryToleranceMeters
        / std::min(width, height);
    for (std::size_t nodeIndex = 0;
         nodeIndex < nodes.size(); ++nodeIndex) {
        auto& nodeOutgoing = outgoing[nodeIndex];
        std::ranges::sort(nodeOutgoing, {},
                          [&](const std::size_t halfEdgeIndex) {
                              return halfEdges[halfEdgeIndex].angleRadians;
                          });
        for (std::size_t index = 0;
             index < nodeOutgoing.size(); ++index) {
            const double first = halfEdges[nodeOutgoing[index]].angleRadians;
            double second = halfEdges[nodeOutgoing[
                (index + 1) % nodeOutgoing.size()]].angleRadians;
            if (index + 1 == nodeOutgoing.size()) second += 2.0 * std::acos(-1.0);
            if (second - first <= angularTolerance) {
                throw std::invalid_argument(
                    "scene fluid face arrangement has overlapping rays");
            }
        }
    }

    std::vector<bool> visited(halfEdges.size(), false);
    std::map<StableId, AreaMoment2> areas;
    std::size_t exteriorCycleCount = 0;
    for (std::size_t start = 0;
         start < halfEdges.size(); ++start) {
        if (visited[start]) continue;
        std::vector<std::size_t> cycle;
        StableId regionId = invalidStableId;
        std::size_t current = start;
        while (true) {
            if (cycle.size() > halfEdges.size() || visited[current]) {
                throw std::invalid_argument(
                    "scene fluid face arrangement cycle is invalid");
            }
            visited[current] = true;
            cycle.push_back(current);
            const auto& halfEdge = halfEdges[current];
            const auto& edge = edges[halfEdge.edgeIndex];
            if (edge.material) {
                const StableId candidate =
                    halfEdge.fromNode == edge.directedFromNode
                        && halfEdge.toNode == edge.directedToNode
                    ? edge.positiveSideRegionId
                    : edge.negativeSideRegionId;
                if (regionId == invalidStableId) regionId = candidate;
                else if (regionId != candidate) {
                    throw std::invalid_argument(
                        "scene fluid face arrangement region winding conflicts");
                }
            }
            const std::size_t reverse = current ^ 1U;
            const auto& destinationOutgoing = outgoing[halfEdge.toNode];
            const auto found = std::ranges::find(
                destinationOutgoing, reverse);
            if (found == destinationOutgoing.end()) {
                throw std::logic_error(
                    "scene fluid face arrangement reverse edge is missing");
            }
            const std::size_t position = static_cast<std::size_t>(
                found - destinationOutgoing.begin());
            current = destinationOutgoing[
                position == 0
                    ? destinationOutgoing.size() - 1
                    : position - 1];
            if (current == start) break;
        }
        AreaMoment2 cycleMoment;
        for (const std::size_t halfEdgeIndex : cycle) {
            const auto& halfEdge = halfEdges[halfEdgeIndex];
            const Point2& first = nodes[halfEdge.fromNode].point;
            const Point2& second = nodes[halfEdge.toNode].point;
            const double cross =
                first.u * second.v - second.u * first.v;
            cycleMoment.areaSquareMeters += 0.5 * cross;
            cycleMoment.firstMomentUMeters3 +=
                (first.u + second.u) * cross / 6.0;
            cycleMoment.firstMomentVMeters3 +=
                (first.v + second.v) * cross / 6.0;
        }
        if (!std::isfinite(cycleMoment.areaSquareMeters)
            || !std::isfinite(cycleMoment.firstMomentUMeters3)
            || !std::isfinite(cycleMoment.firstMomentVMeters3)
            || cycleMoment.areaSquareMeters == 0.0) {
            throw std::invalid_argument(
                "scene fluid face arrangement cycle area is invalid");
        }
        if (cycleMoment.areaSquareMeters < 0.0) {
            ++exteriorCycleCount;
            continue;
        }
        if (regionId == invalidStableId) {
            throw std::invalid_argument(
                "scene fluid face arrangement has an unlabeled interior");
        }
        auto& regionMoment = areas[regionId];
        regionMoment.areaSquareMeters += cycleMoment.areaSquareMeters;
        regionMoment.firstMomentUMeters3 +=
            cycleMoment.firstMomentUMeters3;
        regionMoment.firstMomentVMeters3 +=
            cycleMoment.firstMomentVMeters3;
    }
    if (exteriorCycleCount != 1 || areas.size() < 2) {
        throw std::invalid_argument(
            "scene fluid face arrangement does not form one bounded partition");
    }
    double assignedArea = 0.0;
    for (const auto& [regionId, moment] : areas) {
        static_cast<void>(regionId);
        if (!(moment.areaSquareMeters > 0.0)
            || !std::isfinite(moment.areaSquareMeters)
            || !std::isfinite(moment.firstMomentUMeters3)
            || !std::isfinite(moment.firstMomentVMeters3)) {
            throw std::invalid_argument(
                "scene fluid face arrangement region area is invalid");
        }
        assignedArea += moment.areaSquareMeters;
    }
    const double fullArea = width * height;
    if (!std::isfinite(assignedArea)
        || std::abs(assignedArea - fullArea)
            > settings.areaClosureToleranceSquareMeters) {
        throw std::invalid_argument(
            "scene fluid face arrangement does not close area");
    }
    return areas;
}

std::uint64_t partitionFingerprint(
    const SceneFluidFacePartitionSet& value) {
    Fingerprint f;
    f.integer(value.version);
    f.integer(value.surfaceDefinitionFingerprint);
    f.integer(value.surfaceStateFingerprint);
    f.integer(value.faceLoopFingerprint);
    f.integer(value.structureDefinitionFingerprint);
    f.integer(value.acceptedStepCount);
    f.real(value.simulationTimeSeconds);
    f.real(value.settings.geometryToleranceMeters);
    f.real(value.settings.areaClosureToleranceSquareMeters);
    f.integer(static_cast<std::uint64_t>(value.unresolvedActiveFaceCount));
    f.integer(static_cast<std::uint64_t>(
        value.ignoredSameRegionChainCount));
    f.integer(static_cast<std::uint64_t>(value.segmentPairTestCount));
    f.integer(static_cast<std::uint64_t>(value.loopContainment.size()));
    for (const auto& item : value.loopContainment) {
        f.integer(static_cast<std::uint64_t>(item.loopIndex));
        f.integer(static_cast<std::uint64_t>(item.parentLoopIndex));
        f.integer(static_cast<std::uint64_t>(item.depth));
    }
    f.integer(static_cast<std::uint64_t>(value.partitions.size()));
    for (const auto& item : value.partitions) {
        f.integer(item.stableId);
        f.integer(static_cast<std::uint64_t>(item.activeFaceIndex));
        f.integer(static_cast<std::uint8_t>(item.kind));
        f.integer(item.rootExteriorRegionId);
        f.integer(static_cast<std::uint64_t>(item.firstLoopReference));
        f.integer(static_cast<std::uint64_t>(item.loopReferenceCount));
        f.integer(static_cast<std::uint64_t>(
            item.firstOpenChainReference));
        f.integer(static_cast<std::uint64_t>(
            item.openChainReferenceCount));
        f.integer(static_cast<std::uint64_t>(item.firstRegionArea));
        f.integer(static_cast<std::uint64_t>(item.regionAreaCount));
        f.real(item.faceAreaSquareMeters);
        f.real(item.assignedAreaSquareMeters);
        f.real(item.areaResidualSquareMeters);
    }
    f.integer(static_cast<std::uint64_t>(value.loopReferences.size()));
    for (auto item : value.loopReferences)
        f.integer(static_cast<std::uint64_t>(item));
    f.integer(static_cast<std::uint64_t>(value.openChainReferences.size()));
    for (auto item : value.openChainReferences)
        f.integer(static_cast<std::uint64_t>(item));
    f.integer(static_cast<std::uint64_t>(value.regionAreas.size()));
    for (const auto& item : value.regionAreas) {
        f.integer(item.regionId);
        f.real(item.areaSquareMeters);
        f.real(item.firstMomentMeters3.x);
        f.real(item.firstMomentMeters3.y);
        f.real(item.firstMomentMeters3.z);
        f.real(item.centroidMeters.x);
        f.real(item.centroidMeters.y);
        f.real(item.centroidMeters.z);
    }
    return f.value();
}

void validateSettings(const SceneFluidFacePartitionSettings& settings,
                      const PeriodicCartesianGrid& grid) {
    const Vector3 h = grid.cellSpacingMeters();
    const double minH = std::min({h.x, h.y, h.z});
    const double maxArea = std::max({h.x*h.y, h.x*h.z, h.y*h.z});
    if (!std::isfinite(settings.geometryToleranceMeters)
        || settings.geometryToleranceMeters < 0.0
        || settings.geometryToleranceMeters > 1.0e-6 * minH
        || !std::isfinite(settings.areaClosureToleranceSquareMeters)
        || settings.areaClosureToleranceSquareMeters < 0.0
        || settings.areaClosureToleranceSquareMeters > 1.0e-6 * maxArea)
        throw std::invalid_argument("scene fluid face-partition settings are invalid");
}

SceneFluidFacePartitionSet buildPartitions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSet& loops,
    const SceneFluidFacePartitionSettings& settings,
    const SceneFluidFacePartitionLimits& limits) {
    validateSettings(settings, grid);
    SceneFluidFacePartitionSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.faceLoopFingerprint = loops.fingerprint;
    result.structureDefinitionFingerprint = state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.settings = settings;
    result.loopContainment.resize(loops.loops.size());
    result.ignoredSameRegionChainCount = std::ranges::count_if(
        chains.chains, [](const SceneFluidFaceChain& chain) {
            return chain.negativeSideRegionId == chain.positiveSideRegionId;
        });

    std::vector<std::vector<std::size_t>> loopsByFace(topology.activeFaces.size());
    std::vector<std::vector<std::size_t>> sameRegionLoopsByFace(
        topology.activeFaces.size());
    std::vector<std::vector<Point2>> polygons(loops.loops.size());
    for (std::size_t index = 0; index < loops.loops.size(); ++index) {
        result.loopContainment[index].loopIndex = index;
        const auto& loop = loops.loops[index];
        const auto& chain = chains.chains[loop.chainIndex];
        if (chain.negativeSideRegionId != chain.positiveSideRegionId) {
            loopsByFace[loop.activeFaceIndex].push_back(index);
        } else {
            sameRegionLoopsByFace[loop.activeFaceIndex].push_back(index);
        }
        auto& polygon = polygons[index];
        polygon.reserve(chain.nodeReferenceCount);
        for (std::size_t offset = 0; offset < chain.nodeReferenceCount; ++offset) {
            const std::size_t node = chains.nodeReferences[
                chain.firstNodeReference + offset];
            polygon.push_back(facePoint(
                topology.activeFaces[loop.activeFaceIndex].axis,
                graph.nodes[node].positionMeters));
        }
    }

    std::vector<std::vector<std::size_t>> openChainsByFace(
        topology.activeFaces.size());
    std::vector<std::vector<std::size_t>> sameRegionOpenChainsByFace(
        topology.activeFaces.size());
    for (std::size_t chainIndex = 0;
         chainIndex < chains.chains.size(); ++chainIndex) {
        if (chains.chains[chainIndex].kind
                == SceneFluidFaceChainKind::Open
            && chains.chains[chainIndex].negativeSideRegionId
                != chains.chains[chainIndex].positiveSideRegionId) {
            openChainsByFace[chains.chains[chainIndex].activeFaceIndex]
                .push_back(chainIndex);
        } else if (chains.chains[chainIndex].kind
                       == SceneFluidFaceChainKind::Open) {
            sameRegionOpenChainsByFace[
                chains.chains[chainIndex].activeFaceIndex]
                .push_back(chainIndex);
        }
    }

    for (std::size_t faceIndex = 0;
         faceIndex < topology.activeFaces.size(); ++faceIndex) {
        const auto& face = topology.activeFaces[faceIndex];
        auto& faceLoops = loopsByFace[faceIndex];
        for (std::size_t first = 0; first < faceLoops.size(); ++first) {
            for (std::size_t second = first + 1;
                 second < faceLoops.size(); ++second) {
                const auto& a = polygons[faceLoops[first]];
                const auto& b = polygons[faceLoops[second]];
                for (std::size_t ai = 0; ai < a.size(); ++ai)
                    for (std::size_t bi = 0; bi < b.size(); ++bi) {
                        if (result.segmentPairTestCount
                            == limits.maximumSegmentPairTests)
                            throw std::length_error(
                                "scene fluid face partition exceeds pair-test limit");
                        ++result.segmentPairTestCount;
                        if (segmentsIntersect(a[ai], a[(ai+1)%a.size()],
                                              b[bi], b[(bi+1)%b.size()],
                                              settings.geometryToleranceMeters))
                            throw std::invalid_argument(
                                "scene fluid face loops intersect or touch");
                    }
            }
        }
        for (const std::size_t child : faceLoops) {
            std::size_t parent = noParentFaceLoop;
            for (const std::size_t candidate : faceLoops) {
                if (candidate == child
                    || loops.loops[candidate].areaSquareMeters
                        <= loops.loops[child].areaSquareMeters) continue;
                if (pointInside(polygons[child].front(), polygons[candidate],
                                settings.geometryToleranceMeters)
                    && (parent == noParentFaceLoop
                        || loops.loops[candidate].areaSquareMeters
                            < loops.loops[parent].areaSquareMeters))
                    parent = candidate;
            }
            result.loopContainment[child].parentLoopIndex = parent;
        }
        for (const std::size_t loopIndex : faceLoops) {
            auto& containment = result.loopContainment[loopIndex];
            std::size_t parent = containment.parentLoopIndex;
            while (parent != noParentFaceLoop) {
                ++containment.depth;
                if (containment.depth > faceLoops.size())
                    throw std::invalid_argument(
                        "scene fluid face-loop containment cycle");
                parent = result.loopContainment[parent].parentLoopIndex;
            }
            if (containment.parentLoopIndex != noParentFaceLoop
                && loops.loops[loopIndex].exteriorRegionId
                    != loops.loops[containment.parentLoopIndex].enclosedRegionId)
                throw std::invalid_argument(
                    "scene fluid nested face-loop regions are discontinuous");
        }

        bool boundaryTouch = false;
        for (const std::size_t loopIndex : faceLoops) {
            const auto& chain = chains.chains[loops.loops[loopIndex].chainIndex];
            for (std::size_t offset = 0; offset < chain.nodeReferenceCount; ++offset) {
                const auto& node = graph.nodes[chains.nodeReferences[
                    chain.firstNodeReference + offset]];
                boundaryTouch = boundaryTouch
                    || node.faceBoundaryMask != FaceBoundaryNone
                    || node.authoredOpeningBoundary;
            }
        }
        const auto& faceOpenChains = openChainsByFace[faceIndex];
        const bool soleBoundaryChain =
            faceLoops.empty() && faceOpenChains.size() == 1
            && face.coplanarPatchReferenceCount == 0
            && chains.chains[faceOpenChains.front()]
                    .endpointFaceBoundaryMasks[0]
                != FaceBoundaryNone
            && chains.chains[faceOpenChains.front()]
                    .endpointFaceBoundaryMasks[1]
                != FaceBoundaryNone
            && !chains.chains[faceOpenChains.front()]
                    .endpointOnAuthoredOpening[0]
            && !chains.chains[faceOpenChains.front()]
                    .endpointOnAuthoredOpening[1];
        if (soleBoundaryChain) {
            if (result.partitions.size() == limits.maximumPartitions) {
                throw std::length_error(
                    "scene fluid face partitions exceed their limit");
            }
            const std::size_t chainIndex = faceOpenChains.front();
            const auto& chain = chains.chains[chainIndex];
            SceneFluidFacePartition partition;
            partition.stableId = face.stableId;
            partition.activeFaceIndex = faceIndex;
            partition.kind =
                SceneFluidFacePartitionKind::BoundaryOpenChain;
            partition.firstLoopReference = result.loopReferences.size();
            partition.firstOpenChainReference =
                result.openChainReferences.size();
            partition.openChainReferenceCount = 1;
            result.openChainReferences.push_back(chainIndex);
            partition.faceAreaSquareMeters = faceArea(grid, face.axis);
            const auto positiveMoment = boundaryChainPositiveAreaMoment(
                face, chain, chains, graph, grid, settings, limits,
                result.segmentPairTestCount);
            const auto fullMoment = rectangleAreaMoment(
                faceBounds(grid, face));
            const AreaMoment2 negativeMoment{
                fullMoment.areaSquareMeters
                    - positiveMoment.areaSquareMeters,
                fullMoment.firstMomentUMeters3
                    - positiveMoment.firstMomentUMeters3,
                fullMoment.firstMomentVMeters3
                    - positiveMoment.firstMomentVMeters3,
            };
            std::map<StableId, AreaMoment2> areas;
            areas.emplace(chain.positiveSideRegionId, positiveMoment);
            areas.emplace(chain.negativeSideRegionId, negativeMoment);
            const double planeCoordinate =
                facePlaneCoordinate(grid, face);
            partition.firstRegionArea = result.regionAreas.size();
            for (const auto& [region, moment] : areas) {
                if (!(moment.areaSquareMeters > 0.0)
                    || !std::isfinite(moment.areaSquareMeters)
                    || !std::isfinite(moment.firstMomentUMeters3)
                    || !std::isfinite(moment.firstMomentVMeters3)) {
                    throw std::invalid_argument(
                        "scene fluid boundary face region area is invalid");
                }
                result.regionAreas.push_back(makeRegionArea(
                    region, face.axis, planeCoordinate, moment));
                partition.assignedAreaSquareMeters +=
                    moment.areaSquareMeters;
            }
            partition.regionAreaCount = result.regionAreas.size()
                - partition.firstRegionArea;
            partition.areaResidualSquareMeters =
                partition.assignedAreaSquareMeters
                - partition.faceAreaSquareMeters;
            if (!std::isfinite(partition.areaResidualSquareMeters)
                || std::abs(partition.areaResidualSquareMeters)
                    > settings.areaClosureToleranceSquareMeters) {
                throw std::invalid_argument(
                    "scene fluid boundary face partition does not close area");
            }
            result.partitions.push_back(partition);
            continue;
        }
        const auto& graphRange = graph.faceRanges[faceIndex];
        const bool boundaryArrangement =
            faceLoops.empty() && faceOpenChains.size() > 1
            && face.coplanarPatchReferenceCount == 0
            && supportsRegionSeparatingArrangement(
                faceIndex, faceOpenChains, chains, graph, graphRange);
        if (boundaryArrangement) {
            if (result.partitions.size() == limits.maximumPartitions) {
                throw std::length_error(
                    "scene fluid face partitions exceed their limit");
            }
            SceneFluidFacePartition partition;
            partition.stableId = face.stableId;
            partition.activeFaceIndex = faceIndex;
            partition.kind =
                SceneFluidFacePartitionKind::BoundaryChainArrangement;
            partition.firstLoopReference = result.loopReferences.size();
            partition.firstOpenChainReference =
                result.openChainReferences.size();
            partition.openChainReferenceCount = faceOpenChains.size();
            result.openChainReferences.insert(
                result.openChainReferences.end(),
                faceOpenChains.begin(), faceOpenChains.end());
            partition.faceAreaSquareMeters = faceArea(grid, face.axis);
            const auto areas = boundaryChainArrangementAreaMoments(
                faceIndex, face, faceOpenChains, chains, graph, graphRange,
                grid, settings, limits, result.segmentPairTestCount);
            const double planeCoordinate =
                facePlaneCoordinate(grid, face);
            partition.firstRegionArea = result.regionAreas.size();
            for (const auto& [region, moment] : areas) {
                result.regionAreas.push_back(makeRegionArea(
                    region, face.axis, planeCoordinate, moment));
                partition.assignedAreaSquareMeters +=
                    moment.areaSquareMeters;
            }
            partition.regionAreaCount = result.regionAreas.size()
                - partition.firstRegionArea;
            partition.areaResidualSquareMeters =
                partition.assignedAreaSquareMeters
                - partition.faceAreaSquareMeters;
            if (!std::isfinite(partition.areaResidualSquareMeters)
                || std::abs(partition.areaResidualSquareMeters)
                    > settings.areaClosureToleranceSquareMeters) {
                throw std::invalid_argument(
                    "scene fluid face arrangement does not close area");
            }
            result.partitions.push_back(partition);
            continue;
        }
        const auto& faceSameRegionLoops = sameRegionLoopsByFace[faceIndex];
        const auto& faceSameRegionOpenChains =
            sameRegionOpenChainsByFace[faceIndex];
        if (faceLoops.empty() && faceOpenChains.empty()
            && face.coplanarPatchReferenceCount == 0
            && (!faceSameRegionLoops.empty()
                || !faceSameRegionOpenChains.empty())) {
            StableId regionId = invalidStableId;
            auto acceptSameRegionChain = [&](const SceneFluidFaceChain& chain) {
                if (chain.negativeSideRegionId == invalidStableId
                    || chain.negativeSideRegionId
                        != chain.positiveSideRegionId) {
                    return false;
                }
                if (regionId == invalidStableId) {
                    regionId = chain.negativeSideRegionId;
                }
                return regionId == chain.negativeSideRegionId;
            };
            bool oneRegion = std::ranges::all_of(
                faceSameRegionLoops, [&](const std::size_t loopIndex) {
                    return acceptSameRegionChain(chains.chains[
                        loops.loops[loopIndex].chainIndex]);
                });
            oneRegion = oneRegion && std::ranges::all_of(
                faceSameRegionOpenChains,
                [&](const std::size_t chainIndex) {
                    return acceptSameRegionChain(chains.chains[chainIndex]);
                });
            if (oneRegion) {
                if (result.partitions.size() == limits.maximumPartitions) {
                    throw std::length_error(
                        "scene fluid face partitions exceed their limit");
                }
                SceneFluidFacePartition partition;
                partition.stableId = face.stableId;
                partition.activeFaceIndex = faceIndex;
                partition.kind =
                    SceneFluidFacePartitionKind::SameRegionSheets;
                partition.rootExteriorRegionId = regionId;
                partition.firstLoopReference = result.loopReferences.size();
                partition.loopReferenceCount = faceSameRegionLoops.size();
                result.loopReferences.insert(
                    result.loopReferences.end(), faceSameRegionLoops.begin(),
                    faceSameRegionLoops.end());
                partition.firstOpenChainReference =
                    result.openChainReferences.size();
                partition.openChainReferenceCount =
                    faceSameRegionOpenChains.size();
                result.openChainReferences.insert(
                    result.openChainReferences.end(),
                    faceSameRegionOpenChains.begin(),
                    faceSameRegionOpenChains.end());
                partition.firstRegionArea = result.regionAreas.size();
                partition.regionAreaCount = 1;
                partition.faceAreaSquareMeters = faceArea(grid, face.axis);
                partition.assignedAreaSquareMeters =
                    partition.faceAreaSquareMeters;
                result.regionAreas.push_back(makeRegionArea(
                    regionId, face.axis, facePlaneCoordinate(grid, face),
                    rectangleAreaMoment(faceBounds(grid, face))));
                result.partitions.push_back(partition);
                continue;
            }
        }
        if (faceLoops.empty() || !faceOpenChains.empty()
            || face.coplanarPatchReferenceCount != 0 || boundaryTouch) {
            ++result.unresolvedActiveFaceCount;
            continue;
        }
        StableId rootExterior = invalidStableId;
        for (const std::size_t loopIndex : faceLoops) {
            if (result.loopContainment[loopIndex].parentLoopIndex
                != noParentFaceLoop) continue;
            if (rootExterior == invalidStableId)
                rootExterior = loops.loops[loopIndex].exteriorRegionId;
            else if (rootExterior != loops.loops[loopIndex].exteriorRegionId)
                throw std::invalid_argument(
                    "scene fluid root face loops disagree on exterior region");
        }
        if (result.partitions.size() == limits.maximumPartitions)
            throw std::length_error("scene fluid face partitions exceed their limit");
        SceneFluidFacePartition partition;
        partition.stableId = face.stableId;
        partition.activeFaceIndex = faceIndex;
        partition.kind = SceneFluidFacePartitionKind::ClosedLoops;
        partition.rootExteriorRegionId = rootExterior;
        partition.firstLoopReference = result.loopReferences.size();
        partition.loopReferenceCount = faceLoops.size();
        partition.firstOpenChainReference =
            result.openChainReferences.size();
        result.loopReferences.insert(result.loopReferences.end(),
                                     faceLoops.begin(), faceLoops.end());
        std::map<StableId, AreaMoment2> areas;
        partition.faceAreaSquareMeters = faceArea(grid, face.axis);
        areas[rootExterior] = rectangleAreaMoment(faceBounds(grid, face));
        for (const std::size_t loopIndex : faceLoops) {
            const auto& loop = loops.loops[loopIndex];
            const auto point = facePoint(face.axis, loop.centroidMeters);
            const AreaMoment2 loopMoment{
                loop.areaSquareMeters,
                loop.areaSquareMeters * point.u,
                loop.areaSquareMeters * point.v,
            };
            auto& exterior = areas[loop.exteriorRegionId];
            exterior.areaSquareMeters -= loopMoment.areaSquareMeters;
            exterior.firstMomentUMeters3 -=
                loopMoment.firstMomentUMeters3;
            exterior.firstMomentVMeters3 -=
                loopMoment.firstMomentVMeters3;
            auto& enclosed = areas[loop.enclosedRegionId];
            enclosed.areaSquareMeters += loopMoment.areaSquareMeters;
            enclosed.firstMomentUMeters3 +=
                loopMoment.firstMomentUMeters3;
            enclosed.firstMomentVMeters3 +=
                loopMoment.firstMomentVMeters3;
        }
        const double planeCoordinate = facePlaneCoordinate(grid, face);
        partition.firstRegionArea = result.regionAreas.size();
        for (auto& [region, moment] : areas) {
            if (moment.areaSquareMeters
                < -settings.areaClosureToleranceSquareMeters) {
                throw std::invalid_argument("scene fluid face region area is negative");
            }
            if (moment.areaSquareMeters < 0.0) {
                moment = {};
            }
            result.regionAreas.push_back(makeRegionArea(
                region, face.axis, planeCoordinate, moment));
            partition.assignedAreaSquareMeters += moment.areaSquareMeters;
        }
        partition.regionAreaCount = result.regionAreas.size()
            - partition.firstRegionArea;
        partition.areaResidualSquareMeters =
            partition.assignedAreaSquareMeters - partition.faceAreaSquareMeters;
        if (!std::isfinite(partition.areaResidualSquareMeters)
            || std::abs(partition.areaResidualSquareMeters)
                > settings.areaClosureToleranceSquareMeters)
            throw std::invalid_argument("scene fluid face partition does not close area");
        result.partitions.push_back(partition);
    }

    std::size_t sourceReferenceCount = 0;
    std::size_t referenceCount = 0;
    std::size_t partitionBytes = 0;
    std::size_t containmentBytes = 0;
    std::size_t loopReferenceBytes = 0;
    std::size_t openChainReferenceBytes = 0;
    std::size_t regionBytes = 0;
    std::size_t fixedBytes = 0;
    std::size_t sourceReferenceBytes = 0;
    std::size_t referenceBytes = 0;
    std::size_t total = 0;
    if (!checkedAdd(result.loopReferences.size(),
                    result.openChainReferences.size(),
                    sourceReferenceCount)
        || !checkedAdd(sourceReferenceCount, result.regionAreas.size(),
                       referenceCount)
        || referenceCount > limits.maximumReferences
        || !checkedMultiply(result.partitions.size(),
                            sizeof(SceneFluidFacePartition), partitionBytes)
        || !checkedMultiply(result.loopContainment.size(),
                            sizeof(SceneFluidFaceLoopContainment),
                            containmentBytes)
        || !checkedMultiply(result.loopReferences.size(),
                            sizeof(std::size_t), loopReferenceBytes)
        || !checkedMultiply(result.openChainReferences.size(),
                            sizeof(std::size_t), openChainReferenceBytes)
        || !checkedMultiply(result.regionAreas.size(),
                            sizeof(SceneFluidFaceRegionArea), regionBytes)
        || !checkedAdd(partitionBytes, containmentBytes, fixedBytes)
        || !checkedAdd(loopReferenceBytes, openChainReferenceBytes,
                       sourceReferenceBytes)
        || !checkedAdd(sourceReferenceBytes, regionBytes, referenceBytes)
        || !checkedAdd(fixedBytes, referenceBytes, total)
        || total > limits.maximumPartitionBytes) {
        throw std::length_error(
            "scene fluid face partitions exceed storage limits");
    }
    result.fingerprint = partitionFingerprint(result);
    return result;
}

} // namespace

SceneFluidFacePartitionSet buildSceneFluidFacePartitions(
    const SceneFluidSurfaceDefinition& surface, const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid, const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches, const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings, const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph, const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSet& loops,
    const SceneFluidFacePartitionSettings& settings,
    const SceneFluidFacePartitionLimits& limits) {
    validateSceneFluidFaceLoops(loops, surface, state, grid, candidates,
        intersections, patches, ownership, crossings, topology, graph, chains);
    auto result = buildPartitions(surface, state, grid, topology, graph,
                                  chains, loops, settings, limits);
    validateSceneFluidFacePartitions(result, surface, state, grid, candidates,
        intersections, patches, ownership, crossings, topology, graph, chains, loops);
    return result;
}

void validateSceneFluidFacePartitions(
    const SceneFluidFacePartitionSet& partitions,
    const SceneFluidSurfaceDefinition& surface, const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid, const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches, const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings, const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph, const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSet& loops) {
    validateSceneFluidFaceLoops(loops, surface, state, grid, candidates,
        intersections, patches, ownership, crossings, topology, graph, chains);
    if (partitions.version != sceneFluidFacePartitionVersion
        || partitions.fingerprint == 0
        || partitions.surfaceDefinitionFingerprint != surface.fingerprint
        || partitions.surfaceStateFingerprint != state.fingerprint
        || partitions.faceLoopFingerprint != loops.fingerprint
        || partitions.structureDefinitionFingerprint != state.structureDefinitionFingerprint
        || partitions.acceptedStepCount != state.acceptedStepCount
        || partitions.simulationTimeSeconds != state.simulationTimeSeconds)
        throw std::invalid_argument("scene fluid face-partition identity is invalid");
    const auto expected = buildPartitions(surface, state, grid, topology, graph,
        chains, loops, partitions.settings,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (partitions != expected)
        throw std::invalid_argument(
            "scene fluid face partitions do not match their source loops");
}

} // namespace simwing::fsi::fluid
