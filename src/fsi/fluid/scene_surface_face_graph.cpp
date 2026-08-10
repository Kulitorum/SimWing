#include "fluid/scene_surface_face_graph.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
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

Vec3 interpolate(const Vec3& first,
    const Vec3& second,
    const double parameter) {
    return {
        std::lerp(first.x, second.x, parameter),
        std::lerp(first.y, second.y, parameter),
        std::lerp(first.z, second.z, parameter),
    };
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double distance(const Vec3& first, const Vec3& second) {
    return std::hypot(first.x - second.x,
                      first.y - second.y,
                      first.z - second.z);
}

std::size_t axisIndex(const GridFaceAxis axis) {
    return static_cast<std::size_t>(axis);
}

double facePlaneCoordinate(const PeriodicCartesianGrid& grid,
                           const SceneFluidActiveFace& face) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    if (face.axis == GridFaceAxis::X) {
        return lower.x + static_cast<double>(face.i) * spacing.x;
    }
    if (face.axis == GridFaceAxis::Y) {
        return lower.y + static_cast<double>(face.j) * spacing.y;
    }
    return lower.z + static_cast<double>(face.k) * spacing.z;
}

struct FaceRectangle {
    std::size_t uAxis = 0;
    std::size_t vAxis = 1;
    double uLower = 0.0;
    double uUpper = 0.0;
    double vLower = 0.0;
    double vUpper = 0.0;
};

FaceRectangle faceRectangle(const PeriodicCartesianGrid& grid,
                            const SceneFluidActiveFace& face) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    if (face.axis == GridFaceAxis::X) {
        const double uLower = lower.y
            + static_cast<double>(face.j) * spacing.y;
        const double vLower = lower.z
            + static_cast<double>(face.k) * spacing.z;
        return {1, 2, uLower, uLower + spacing.y,
                vLower, vLower + spacing.z};
    }
    if (face.axis == GridFaceAxis::Y) {
        const double uLower = lower.x
            + static_cast<double>(face.i) * spacing.x;
        const double vLower = lower.z
            + static_cast<double>(face.k) * spacing.z;
        return {0, 2, uLower, uLower + spacing.x,
                vLower, vLower + spacing.z};
    }
    const double uLower = lower.x
        + static_cast<double>(face.i) * spacing.x;
    const double vLower = lower.y
        + static_cast<double>(face.j) * spacing.y;
    return {0, 1, uLower, uLower + spacing.x,
            vLower, vLower + spacing.y};
}

bool near(const double first, const double second, const double tolerance) {
    return std::abs(first - second) <= tolerance;
}

std::uint8_t faceBoundaryMask(
    const Vec3& position,
    const FaceRectangle& rectangle,
    const double tolerance) {
    std::uint8_t result = FaceBoundaryNone;
    const double u = coordinate(position, rectangle.uAxis);
    const double v = coordinate(position, rectangle.vAxis);
    if (near(u, rectangle.uLower, tolerance)) {
        result |= FaceBoundaryUMinus;
    }
    if (near(u, rectangle.uUpper, tolerance)) {
        result |= FaceBoundaryUPlus;
    }
    if (near(v, rectangle.vLower, tolerance)) {
        result |= FaceBoundaryVMinus;
    }
    if (near(v, rectangle.vUpper, tolerance)) {
        result |= FaceBoundaryVPlus;
    }
    return result;
}

Vec3 canonicalGridEdgePoint(
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidActiveFace& face,
    const SceneFluidSurfaceTriangle& triangle,
    const FaceRectangle& rectangle,
    const std::uint8_t boundaryMask,
    const Vec3& supplied,
    const double tolerance) {
    Vec3 result = supplied;
    std::array<bool, 3> fixed{};
    const std::size_t faceAxis = axisIndex(face.axis);
    fixed[faceAxis] = true;
    setCoordinate(result, faceAxis, facePlaneCoordinate(grid, face));
    const auto fixBoundary = [&](const std::uint8_t bit,
                                 const std::size_t axis,
                                 const double value) {
        if ((boundaryMask & bit) == 0) {
            return;
        }
        if (fixed[axis]
            && coordinate(result, axis) != value) {
            throw std::invalid_argument(
                "scene fluid face-graph grid-edge endpoint has conflicting planes");
        }
        fixed[axis] = true;
        setCoordinate(result, axis, value);
    };
    fixBoundary(FaceBoundaryUMinus, rectangle.uAxis, rectangle.uLower);
    fixBoundary(FaceBoundaryUPlus, rectangle.uAxis, rectangle.uUpper);
    fixBoundary(FaceBoundaryVMinus, rectangle.vAxis, rectangle.vLower);
    fixBoundary(FaceBoundaryVPlus, rectangle.vAxis, rectangle.vUpper);

    const std::size_t fixedCount = std::ranges::count(fixed, true);
    if (fixedCount == 2) {
        const auto& first = state.vertices[
            triangle.vertexIndices[0]].positionMeters;
        const auto& second = state.vertices[
            triangle.vertexIndices[1]].positionMeters;
        const auto& third = state.vertices[
            triangle.vertexIndices[2]].positionMeters;
        const Vec3 normal = cross(
            subtract(second, first), subtract(third, first));
        const auto unknown = static_cast<std::size_t>(
            std::ranges::find(fixed, false) - fixed.begin());
        const double denominator = coordinate(normal, unknown);
        if (!std::isfinite(denominator) || denominator == 0.0) {
            throw std::invalid_argument(
                "scene fluid face-graph triangle is parallel to its grid edge");
        }
        double numerator = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (!fixed[axis]) {
                continue;
            }
            numerator = std::fma(
                coordinate(normal, axis),
                coordinate(result, axis) - coordinate(first, axis),
                numerator);
        }
        setCoordinate(
            result, unknown,
            coordinate(first, unknown) - numerator / denominator);
    } else if (fixedCount != 3) {
        throw std::invalid_argument(
            "scene fluid face-graph grid-edge endpoint has incomplete provenance");
    }
    if (!std::isfinite(result.x) || !std::isfinite(result.y)
        || !std::isfinite(result.z)
        || distance(result, supplied) > tolerance) {
        throw std::invalid_argument(
            "scene fluid face-graph canonical grid-edge endpoint is inconsistent");
    }
    return result;
}

using EdgeKey = std::pair<StableId, StableId>;

EdgeKey edgeKey(StableId first, StableId second) {
    if (second < first) {
        std::swap(first, second);
    }
    return {first, second};
}

struct OpeningTopology {
    std::set<StableId> vertices;
    std::set<EdgeKey> edges;
};

OpeningTopology openingTopology(
    const SceneFluidSurfaceDefinition& surface) {
    OpeningTopology result;
    for (const auto& opening : surface.openings) {
        for (std::size_t index = 0;
             index < opening.orderedVertexIndices.size(); ++index) {
            const std::size_t next =
                (index + 1) % opening.orderedVertexIndices.size();
            const StableId first =
                surface.vertices[opening.orderedVertexIndices[index]].id;
            const StableId second =
                surface.vertices[opening.orderedVertexIndices[next]].id;
            result.vertices.insert(first);
            result.edges.insert(edgeKey(first, second));
        }
    }
    return result;
}

// kind, then the stable provenance fields. The active face is supplied by the
// surrounding per-face map and deliberately is not repeated here.
using NodeKey = std::tuple<std::uint8_t,
                           StableId,
                           StableId,
                           StableId,
                           std::uint8_t>;

struct EndpointNode {
    NodeKey key;
    SceneFluidFaceGraphNode node;
};

std::uint64_t nodeStableId(const std::uint64_t activeFaceStableId,
                           const NodeKey& key) {
    Fingerprint fingerprint;
    fingerprint.integer(sceneFluidFaceGraphVersion);
    fingerprint.integer(activeFaceStableId);
    fingerprint.integer(std::get<0>(key));
    fingerprint.integer(std::get<1>(key));
    fingerprint.integer(std::get<2>(key));
    fingerprint.integer(std::get<3>(key));
    fingerprint.integer(std::get<4>(key));
    return fingerprint.value();
}

EndpointNode endpointNode(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidActiveFace& face,
    const std::size_t activeFaceIndex,
    const SceneFluidFaceCrossing& crossing,
    const SceneFluidClippedVertex& endpoint,
    const OpeningTopology& authoredOpenings,
    const SceneFluidFaceGraphSettings& settings) {
    const auto& triangle = surface.triangles[crossing.triangleIndex];
    const FaceRectangle rectangle = faceRectangle(grid, face);
    const std::uint8_t suppliedBoundary = faceBoundaryMask(
        endpoint.positionMeters, rectangle, settings.endpointToleranceMeters);
    std::array<bool, 3> zero{};
    std::size_t zeroCount = 0;
    for (std::size_t corner = 0; corner < 3; ++corner) {
        zero[corner] = std::abs(endpoint.barycentricCoordinates[corner])
            <= settings.barycentricZeroTolerance;
        zeroCount += zero[corner] ? 1 : 0;
    }

    EndpointNode result;
    result.node.activeFaceIndex = activeFaceIndex;
    if (zeroCount >= 2) {
        const auto greatest = std::max_element(
            endpoint.barycentricCoordinates.begin(),
            endpoint.barycentricCoordinates.end());
        const std::size_t corner = static_cast<std::size_t>(
            greatest - endpoint.barycentricCoordinates.begin());
        const std::size_t vertexIndex = triangle.vertexIndices[corner];
        const StableId vertexId = surface.vertices[vertexIndex].id;
        const Vec3 canonical = state.vertices[vertexIndex].positionMeters;
        if (distance(canonical, endpoint.positionMeters)
            > settings.endpointToleranceMeters) {
            throw std::invalid_argument(
                "scene fluid face-graph vertex endpoint is inconsistent");
        }
        result.node.kind = SceneFluidFaceNodeKind::SurfaceVertex;
        result.node.positionMeters = canonical;
        result.node.sourceVertexId = vertexId;
        result.node.faceBoundaryMask = faceBoundaryMask(
            canonical, rectangle, settings.endpointToleranceMeters);
        result.node.authoredOpeningBoundary =
            authoredOpenings.vertices.contains(vertexId);
        result.key = {static_cast<std::uint8_t>(result.node.kind),
                      vertexId, 0, 0, result.node.faceBoundaryMask};
    } else if (zeroCount == 1) {
        std::array<std::size_t, 2> corners{};
        std::size_t count = 0;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            if (!zero[corner]) {
                corners[count++] = corner;
            }
        }
        const std::size_t firstVertexIndex =
            triangle.vertexIndices[corners[0]];
        const std::size_t secondVertexIndex =
            triangle.vertexIndices[corners[1]];
        StableId firstId = surface.vertices[firstVertexIndex].id;
        StableId secondId = surface.vertices[secondVertexIndex].id;
        Vec3 first = state.vertices[firstVertexIndex].positionMeters;
        Vec3 second = state.vertices[secondVertexIndex].positionMeters;
        if (secondId < firstId) {
            std::swap(firstId, secondId);
            std::swap(first, second);
        }
        const std::size_t faceAxis = axisIndex(face.axis);
        const double plane = facePlaneCoordinate(grid, face);
        const double denominator =
            coordinate(second, faceAxis) - coordinate(first, faceAxis);
        if (!std::isfinite(denominator) || denominator == 0.0) {
            throw std::invalid_argument(
                "scene fluid face-graph source edge lies in the face");
        }
        const double parameter =
            (plane - coordinate(first, faceAxis)) / denominator;
        if (!std::isfinite(parameter)
            || parameter < -settings.barycentricZeroTolerance
            || parameter > 1.0 + settings.barycentricZeroTolerance) {
            throw std::invalid_argument(
                "scene fluid face-graph source edge misses the face");
        }
        Vec3 canonical = interpolate(
            first, second, std::clamp(parameter, 0.0, 1.0));
        setCoordinate(canonical, faceAxis, plane);
        if (distance(canonical, endpoint.positionMeters)
            > settings.endpointToleranceMeters) {
            throw std::invalid_argument(
                "scene fluid face-graph edge endpoint is inconsistent");
        }
        result.node.kind = SceneFluidFaceNodeKind::SurfaceEdge;
        result.node.positionMeters = canonical;
        result.node.sourceEdgeVertexIds = {firstId, secondId};
        result.node.faceBoundaryMask = faceBoundaryMask(
            canonical, rectangle, settings.endpointToleranceMeters);
        result.node.authoredOpeningBoundary =
            authoredOpenings.edges.contains({firstId, secondId});
        result.key = {static_cast<std::uint8_t>(result.node.kind),
                      firstId, secondId, 0, result.node.faceBoundaryMask};
    } else {
        if (suppliedBoundary == FaceBoundaryNone) {
            throw std::invalid_argument(
                "scene fluid face-graph endpoint has no surface or grid-edge provenance");
        }
        result.node.kind = SceneFluidFaceNodeKind::GridEdge;
        result.node.positionMeters = canonicalGridEdgePoint(
            state, grid, face, triangle, rectangle,
            suppliedBoundary, endpoint.positionMeters,
            settings.endpointToleranceMeters);
        result.node.sourceTriangleId = crossing.triangleId;
        result.node.faceBoundaryMask = suppliedBoundary;
        result.key = {static_cast<std::uint8_t>(result.node.kind),
                      crossing.triangleId, 0, 0, suppliedBoundary};
    }
    result.node.stableId = nodeStableId(face.stableId, result.key);
    return result;
}

struct PendingSegment {
    std::size_t crossingIndex = 0;
    NodeKey first;
    NodeKey second;
};

std::uint64_t graphFingerprint(const SceneFluidFaceGraph& graph) {
    Fingerprint fingerprint;
    fingerprint.integer(graph.version);
    fingerprint.integer(graph.surfaceDefinitionFingerprint);
    fingerprint.integer(graph.surfaceStateFingerprint);
    fingerprint.integer(graph.faceTopologyFingerprint);
    fingerprint.integer(graph.structureDefinitionFingerprint);
    fingerprint.integer(graph.acceptedStepCount);
    fingerprint.real(graph.simulationTimeSeconds);
    fingerprint.real(graph.settings.endpointToleranceMeters);
    fingerprint.real(graph.settings.barycentricZeroTolerance);
    fingerprint.integer(static_cast<std::uint64_t>(graph.degreeOneNodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(graph.degreeTwoNodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(graph.higherDegreeNodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        graph.authoredOpeningNodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(graph.faceRanges.size()));
    for (const auto& range : graph.faceRanges) {
        fingerprint.integer(static_cast<std::uint64_t>(range.activeFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(range.firstNode));
        fingerprint.integer(static_cast<std::uint64_t>(range.nodeCount));
        fingerprint.integer(static_cast<std::uint64_t>(range.firstSegment));
        fingerprint.integer(static_cast<std::uint64_t>(range.segmentCount));
    }
    fingerprint.integer(static_cast<std::uint64_t>(graph.nodes.size()));
    for (const auto& node : graph.nodes) {
        fingerprint.integer(node.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(node.activeFaceIndex));
        fingerprint.enumeration(node.kind);
        fingerprint.real(node.positionMeters.x);
        fingerprint.real(node.positionMeters.y);
        fingerprint.real(node.positionMeters.z);
        fingerprint.integer(node.sourceVertexId);
        fingerprint.integer(node.sourceEdgeVertexIds[0]);
        fingerprint.integer(node.sourceEdgeVertexIds[1]);
        fingerprint.integer(node.sourceTriangleId);
        fingerprint.integer(node.faceBoundaryMask);
        fingerprint.integer(static_cast<std::uint8_t>(
            node.authoredOpeningBoundary ? 1 : 0));
        fingerprint.integer(static_cast<std::uint64_t>(
            node.firstIncidentSegmentReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            node.incidentSegmentReferenceCount));
    }
    fingerprint.integer(static_cast<std::uint64_t>(graph.segments.size()));
    for (const auto& segment : graph.segments) {
        fingerprint.integer(segment.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(segment.activeFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(segment.crossingIndex));
        fingerprint.integer(static_cast<std::uint64_t>(segment.nodeIndices[0]));
        fingerprint.integer(static_cast<std::uint64_t>(segment.nodeIndices[1]));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        graph.incidentSegmentReferences.size()));
    for (const std::size_t reference : graph.incidentSegmentReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    return fingerprint.value();
}

void validateSettings(const SceneFluidFaceGraphSettings& settings,
                      const PeriodicCartesianGrid& grid) {
    const Vector3 spacing = grid.cellSpacingMeters();
    const double minimumSpacing = std::min({spacing.x, spacing.y, spacing.z});
    if (!std::isfinite(settings.endpointToleranceMeters)
        || settings.endpointToleranceMeters < 0.0
        || settings.endpointToleranceMeters > 1.0e-6 * minimumSpacing
        || !std::isfinite(settings.barycentricZeroTolerance)
        || settings.barycentricZeroTolerance < 0.0
        || settings.barycentricZeroTolerance > 1.0e-8) {
        throw std::invalid_argument(
            "scene fluid face-graph settings are invalid");
    }
}

SceneFluidFaceGraph buildGraph(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraphSettings& settings,
    const SceneFluidFaceGraphLimits& limits) {
    validateSettings(settings, grid);
    SceneFluidFaceGraph result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.faceTopologyFingerprint = topology.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.settings = settings;
    result.faceRanges.reserve(topology.activeFaces.size());
    result.segments.reserve(crossings.crossings.size());
    const auto authoredOpenings = openingTopology(surface);
    std::set<std::uint64_t> nodeStableIds;

    for (std::size_t activeFaceIndex = 0;
         activeFaceIndex < topology.activeFaces.size(); ++activeFaceIndex) {
        const auto& face = topology.activeFaces[activeFaceIndex];
        SceneFluidFaceGraphRange range;
        range.activeFaceIndex = activeFaceIndex;
        range.firstNode = result.nodes.size();
        range.firstSegment = result.segments.size();

        std::map<NodeKey, SceneFluidFaceGraphNode> nodes;
        std::vector<PendingSegment> pending;
        for (const std::size_t crossingIndex :
             topology.crossingsForFace(face)) {
            const auto& crossing = crossings.crossings[crossingIndex];
            EndpointNode first = endpointNode(
                surface, state, grid, face, activeFaceIndex, crossing,
                crossing.first, authoredOpenings, settings);
            EndpointNode second = endpointNode(
                surface, state, grid, face, activeFaceIndex, crossing,
                crossing.second, authoredOpenings, settings);
            if (first.key == second.key) {
                throw std::invalid_argument(
                    "scene fluid face-graph segment collapses to one node");
            }
            for (auto* endpoint : {&first, &second}) {
                const auto [found, inserted] = nodes.emplace(
                    endpoint->key, endpoint->node);
                if (!inserted && found->second != endpoint->node) {
                    throw std::invalid_argument(
                        "scene fluid face-graph node provenance is inconsistent");
                }
            }
            pending.push_back({crossingIndex, first.key, second.key});
        }
        if (result.nodes.size() > limits.maximumNodes
            || nodes.size() > limits.maximumNodes - result.nodes.size()
            || result.segments.size() > limits.maximumSegments
            || pending.size()
                > limits.maximumSegments - result.segments.size()) {
            throw std::length_error(
                "scene fluid face graph exceeds its node or segment limit");
        }

        std::map<NodeKey, std::size_t> nodeIndices;
        std::vector<std::vector<std::size_t>> incidents;
        incidents.reserve(nodes.size());
        for (const auto& [key, node] : nodes) {
            if (!nodeStableIds.insert(node.stableId).second) {
                throw std::invalid_argument(
                    "scene fluid face-graph stable-ID collision");
            }
            nodeIndices.emplace(key, result.nodes.size());
            result.nodes.push_back(node);
            incidents.emplace_back();
        }
        range.nodeCount = nodes.size();

        for (const auto& item : pending) {
            SceneFluidFaceGraphSegment segment;
            segment.stableId = crossings.crossings[item.crossingIndex].stableId;
            segment.activeFaceIndex = activeFaceIndex;
            segment.crossingIndex = item.crossingIndex;
            segment.nodeIndices = {
                nodeIndices.at(item.first),
                nodeIndices.at(item.second),
            };
            const std::size_t segmentIndex = result.segments.size();
            incidents[segment.nodeIndices[0] - range.firstNode].push_back(
                segmentIndex);
            incidents[segment.nodeIndices[1] - range.firstNode].push_back(
                segmentIndex);
            result.segments.push_back(segment);
        }
        range.segmentCount = pending.size();

        for (std::size_t localNode = 0;
             localNode < incidents.size(); ++localNode) {
            auto& node = result.nodes[range.firstNode + localNode];
            auto& nodeIncidents = incidents[localNode];
            std::sort(nodeIncidents.begin(), nodeIncidents.end(),
                      [&](const std::size_t first, const std::size_t second) {
                          return result.segments[first].stableId
                              < result.segments[second].stableId;
                      });
            if (result.incidentSegmentReferences.size()
                    > limits.maximumIncidentReferences
                || nodeIncidents.size()
                    > limits.maximumIncidentReferences
                        - result.incidentSegmentReferences.size()) {
                throw std::length_error(
                    "scene fluid face graph exceeds its incident-reference limit");
            }
            node.firstIncidentSegmentReference =
                result.incidentSegmentReferences.size();
            node.incidentSegmentReferenceCount = nodeIncidents.size();
            result.incidentSegmentReferences.insert(
                result.incidentSegmentReferences.end(),
                nodeIncidents.begin(), nodeIncidents.end());
            if (nodeIncidents.size() == 1) {
                ++result.degreeOneNodeCount;
            } else if (nodeIncidents.size() == 2) {
                ++result.degreeTwoNodeCount;
            } else if (nodeIncidents.size() > 2) {
                ++result.higherDegreeNodeCount;
            }
            if (node.authoredOpeningBoundary) {
                ++result.authoredOpeningNodeCount;
            }
        }
        result.faceRanges.push_back(range);
    }

    std::size_t rangeBytes = 0;
    std::size_t nodeBytes = 0;
    std::size_t segmentBytes = 0;
    std::size_t incidentBytes = 0;
    std::size_t firstTotal = 0;
    std::size_t secondTotal = 0;
    std::size_t totalBytes = 0;
    if (!checkedMultiply(result.faceRanges.size(),
                         sizeof(SceneFluidFaceGraphRange), rangeBytes)
        || !checkedMultiply(result.nodes.size(),
                            sizeof(SceneFluidFaceGraphNode), nodeBytes)
        || !checkedMultiply(result.segments.size(),
                            sizeof(SceneFluidFaceGraphSegment), segmentBytes)
        || !checkedMultiply(result.incidentSegmentReferences.size(),
                            sizeof(std::size_t), incidentBytes)
        || !checkedAdd(rangeBytes, nodeBytes, firstTotal)
        || !checkedAdd(segmentBytes, incidentBytes, secondTotal)
        || !checkedAdd(firstTotal, secondTotal, totalBytes)
        || totalBytes > limits.maximumGraphBytes) {
        throw std::length_error(
            "scene fluid face graph exceeds its storage limit");
    }
    result.fingerprint = graphFingerprint(result);
    return result;
}

} // namespace

SceneFluidFaceGraph buildSceneFluidFaceGraph(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraphSettings& settings,
    const SceneFluidFaceGraphLimits& limits) {
    validateSceneFluidFaceTopology(
        topology, surface, state, grid, candidates, intersections, patches,
        ownership, crossings);
    SceneFluidFaceGraph result = buildGraph(
        surface, state, grid, crossings, topology, settings, limits);
    validateSceneFluidFaceGraph(
        result, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology);
    return result;
}

void validateSceneFluidFaceGraph(
    const SceneFluidFaceGraph& graph,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology) {
    validateSceneFluidFaceTopology(
        topology, surface, state, grid, candidates, intersections, patches,
        ownership, crossings);
    if (graph.version != sceneFluidFaceGraphVersion
        || graph.fingerprint == 0
        || graph.surfaceDefinitionFingerprint != surface.fingerprint
        || graph.surfaceStateFingerprint != state.fingerprint
        || graph.faceTopologyFingerprint != topology.fingerprint
        || graph.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || graph.acceptedStepCount != state.acceptedStepCount
        || graph.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid face-graph identity is invalid");
    }
    const SceneFluidFaceGraph expected = buildGraph(
        surface,
        state,
        grid,
        crossings,
        topology,
        graph.settings,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (graph != expected) {
        throw std::invalid_argument(
            "scene fluid face graph does not match its source topology");
    }
}

} // namespace simwing::fsi::fluid
