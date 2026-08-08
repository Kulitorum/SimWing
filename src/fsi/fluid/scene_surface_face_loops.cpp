#include "fluid/scene_surface_face_loops.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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

struct Point2 {
    double u = 0.0;
    double v = 0.0;
};

Point2 facePoint(const GridFaceAxis axis, const Vec3& position) {
    if (axis == GridFaceAxis::X) {
        return {position.y, position.z};
    }
    if (axis == GridFaceAxis::Y) {
        // Z x X = +Y: this is the right-handed chart seen from +Y.
        return {position.z, position.x};
    }
    return {position.x, position.y};
}

Vec3 physicalPoint(const GridFaceAxis axis,
                   const double planeCoordinate,
                   const Point2& point) {
    if (axis == GridFaceAxis::X) {
        return {planeCoordinate, point.u, point.v};
    }
    if (axis == GridFaceAxis::Y) {
        return {point.v, planeCoordinate, point.u};
    }
    return {point.u, point.v, planeCoordinate};
}

double planeCoordinate(const PeriodicCartesianGrid& grid,
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

double cross2(const Point2& first,
              const Point2& second,
              const Point2& third) {
    return (second.u - first.u) * (third.v - first.v)
        - (second.v - first.v) * (third.u - first.u);
}

bool onSegment(const Point2& first,
               const Point2& second,
               const Point2& point,
               const double tolerance) {
    const double length = std::hypot(second.u - first.u,
                                     second.v - first.v);
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
                const Point2& third,
                const double tolerance) {
    const double value = cross2(first, second, third);
    const double length = std::hypot(second.u - first.u,
                                     second.v - first.v);
    const double bound = tolerance * std::max(tolerance, length);
    if (value > bound) {
        return 1;
    }
    if (value < -bound) {
        return -1;
    }
    return 0;
}

bool segmentsIntersect(const Point2& first,
                       const Point2& second,
                       const Point2& third,
                       const Point2& fourth,
                       const double tolerance) {
    const int firstThird = orientation(first, second, third, tolerance);
    const int firstFourth = orientation(first, second, fourth, tolerance);
    const int thirdFirst = orientation(third, fourth, first, tolerance);
    const int thirdSecond = orientation(third, fourth, second, tolerance);
    if (firstThird != firstFourth && thirdFirst != thirdSecond
        && firstThird != 0 && firstFourth != 0
        && thirdFirst != 0 && thirdSecond != 0) {
        return true;
    }
    return (firstThird == 0 && onSegment(first, second, third, tolerance))
        || (firstFourth == 0 && onSegment(first, second, fourth, tolerance))
        || (thirdFirst == 0 && onSegment(third, fourth, first, tolerance))
        || (thirdSecond == 0 && onSegment(third, fourth, second, tolerance));
}

void rejectSelfIntersection(const std::vector<Point2>& points,
                            const double tolerance) {
    for (std::size_t first = 0; first < points.size(); ++first) {
        const std::size_t firstNext = (first + 1) % points.size();
        for (std::size_t second = first + 1;
             second < points.size(); ++second) {
            const std::size_t secondNext = (second + 1) % points.size();
            if (first == second || firstNext == second
                || secondNext == first) {
                continue;
            }
            if (segmentsIntersect(points[first], points[firstNext],
                                  points[second], points[secondNext],
                                  tolerance)) {
                throw std::invalid_argument(
                    "scene fluid face loop self-intersects");
            }
        }
    }
}

std::uint64_t loopsFingerprint(const SceneFluidFaceLoopSet& loops) {
    Fingerprint fingerprint;
    fingerprint.integer(loops.version);
    fingerprint.integer(loops.surfaceDefinitionFingerprint);
    fingerprint.integer(loops.surfaceStateFingerprint);
    fingerprint.integer(loops.faceChainFingerprint);
    fingerprint.integer(loops.structureDefinitionFingerprint);
    fingerprint.integer(loops.acceptedStepCount);
    fingerprint.real(loops.simulationTimeSeconds);
    fingerprint.real(loops.settings.minimumAbsoluteAreaSquareMeters);
    fingerprint.real(loops.settings.intersectionToleranceMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        loops.unresolvedOpenChainCount));
    fingerprint.real(loops.summedLoopAreaSquareMeters);
    fingerprint.integer(static_cast<std::uint64_t>(loops.loops.size()));
    for (const auto& loop : loops.loops) {
        fingerprint.integer(loop.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(loop.chainIndex));
        fingerprint.integer(static_cast<std::uint64_t>(loop.activeFaceIndex));
        fingerprint.real(loop.signedAreaSquareMeters);
        fingerprint.real(loop.areaSquareMeters);
        fingerprint.real(loop.centroidMeters.x);
        fingerprint.real(loop.centroidMeters.y);
        fingerprint.real(loop.centroidMeters.z);
        fingerprint.integer(loop.enclosedRegionId);
        fingerprint.integer(loop.exteriorRegionId);
        fingerprint.integer(static_cast<std::uint8_t>(
            loop.positiveSideIsInterior ? 1 : 0));
    }
    return fingerprint.value();
}

void validateSettings(const SceneFluidFaceLoopSettings& settings,
                      const PeriodicCartesianGrid& grid) {
    const Vector3 spacing = grid.cellSpacingMeters();
    const double minimumSpacing = std::min({spacing.x, spacing.y, spacing.z});
    const double maximumFaceArea = std::max({
        spacing.x * spacing.y,
        spacing.x * spacing.z,
        spacing.y * spacing.z,
    });
    if (!std::isfinite(settings.minimumAbsoluteAreaSquareMeters)
        || settings.minimumAbsoluteAreaSquareMeters < 0.0
        || settings.minimumAbsoluteAreaSquareMeters
            > maximumFaceArea
        || !std::isfinite(settings.intersectionToleranceMeters)
        || settings.intersectionToleranceMeters < 0.0
        || settings.intersectionToleranceMeters > 1.0e-6 * minimumSpacing) {
        throw std::invalid_argument(
            "scene fluid face-loop settings are invalid");
    }
}

SceneFluidFaceLoopSet buildLoops(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSettings& settings,
    const SceneFluidFaceLoopLimits& limits) {
    validateSettings(settings, grid);
    SceneFluidFaceLoopSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.faceChainFingerprint = chains.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.settings = settings;
    result.unresolvedOpenChainCount = chains.openChainCount;
    result.loops.reserve(chains.closedChainCount);
    std::set<std::uint64_t> stableIds;

    for (std::size_t chainIndex = 0;
         chainIndex < chains.chains.size(); ++chainIndex) {
        const auto& chain = chains.chains[chainIndex];
        if (chain.kind == SceneFluidFaceChainKind::Open) {
            continue;
        }
        if (result.loops.size() == limits.maximumLoops) {
            throw std::length_error(
                "scene fluid face loops exceed their count limit");
        }
        if (chain.nodeReferenceCount < 3
            || chain.nodeReferenceCount != chain.segmentReferenceCount) {
            throw std::invalid_argument(
                "scene fluid closed face chain cannot form a polygon");
        }
        const auto& face = topology.activeFaces[chain.activeFaceIndex];
        std::vector<Point2> points;
        points.reserve(chain.nodeReferenceCount);
        for (std::size_t offset = 0;
             offset < chain.nodeReferenceCount; ++offset) {
            const std::size_t nodeIndex = chains.nodeReferences[
                chain.firstNodeReference + offset];
            points.push_back(facePoint(
                face.axis, graph.nodes[nodeIndex].positionMeters));
        }
        rejectSelfIntersection(points, settings.intersectionToleranceMeters);

        double twiceSignedArea = 0.0;
        double centroidUNumerator = 0.0;
        double centroidVNumerator = 0.0;
        for (std::size_t index = 0; index < points.size(); ++index) {
            const auto& first = points[index];
            const auto& second = points[(index + 1) % points.size()];
            const double term = first.u * second.v - second.u * first.v;
            twiceSignedArea += term;
            centroidUNumerator += (first.u + second.u) * term;
            centroidVNumerator += (first.v + second.v) * term;
        }
        const double signedArea = 0.5 * twiceSignedArea;
        const double area = std::abs(signedArea);
        if (!std::isfinite(signedArea)
            || !(area > settings.minimumAbsoluteAreaSquareMeters)) {
            throw std::invalid_argument(
                "scene fluid face loop area is degenerate");
        }
        const Point2 centroid{
            centroidUNumerator / (3.0 * twiceSignedArea),
            centroidVNumerator / (3.0 * twiceSignedArea),
        };
        if (!std::isfinite(centroid.u) || !std::isfinite(centroid.v)) {
            throw std::invalid_argument(
                "scene fluid face loop centroid is not finite");
        }

        SceneFluidFaceLoop loop;
        loop.stableId = chain.stableId;
        if (!stableIds.insert(loop.stableId).second) {
            throw std::invalid_argument(
                "scene fluid face-loop stable-ID collision");
        }
        loop.chainIndex = chainIndex;
        loop.activeFaceIndex = chain.activeFaceIndex;
        loop.signedAreaSquareMeters = signedArea;
        loop.areaSquareMeters = area;
        loop.centroidMeters = physicalPoint(
            face.axis, planeCoordinate(grid, face), centroid);
        loop.positiveSideIsInterior = signedArea > 0.0;
        loop.enclosedRegionId = loop.positiveSideIsInterior
            ? chain.positiveSideRegionId : chain.negativeSideRegionId;
        loop.exteriorRegionId = loop.positiveSideIsInterior
            ? chain.negativeSideRegionId : chain.positiveSideRegionId;
        result.summedLoopAreaSquareMeters += area;
        result.loops.push_back(loop);
    }
    if (!std::isfinite(result.summedLoopAreaSquareMeters)) {
        throw std::length_error(
            "scene fluid summed face-loop area is not finite");
    }
    std::size_t loopBytes = 0;
    if (!checkedMultiply(result.loops.size(),
                         sizeof(SceneFluidFaceLoop), loopBytes)
        || loopBytes > limits.maximumLoopBytes) {
        throw std::length_error(
            "scene fluid face loops exceed their storage limit");
    }
    result.fingerprint = loopsFingerprint(result);
    return result;
}

} // namespace

SceneFluidFaceLoopSet buildSceneFluidFaceLoops(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains,
    const SceneFluidFaceLoopSettings& settings,
    const SceneFluidFaceLoopLimits& limits) {
    validateSceneFluidFaceChains(
        chains, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology, graph);
    SceneFluidFaceLoopSet result = buildLoops(
        surface, state, grid, topology, graph, chains, settings, limits);
    validateSceneFluidFaceLoops(
        result, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology, graph, chains);
    return result;
}

void validateSceneFluidFaceLoops(
    const SceneFluidFaceLoopSet& loops,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainSet& chains) {
    validateSceneFluidFaceChains(
        chains, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology, graph);
    if (loops.version != sceneFluidFaceLoopVersion
        || loops.fingerprint == 0
        || loops.surfaceDefinitionFingerprint != surface.fingerprint
        || loops.surfaceStateFingerprint != state.fingerprint
        || loops.faceChainFingerprint != chains.fingerprint
        || loops.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || loops.acceptedStepCount != state.acceptedStepCount
        || loops.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid face-loop identity is invalid");
    }
    const SceneFluidFaceLoopSet expected = buildLoops(
        surface,
        state,
        grid,
        topology,
        graph,
        chains,
        loops.settings,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (loops != expected) {
        throw std::invalid_argument(
            "scene fluid face loops do not match their source chains");
    }
}

} // namespace simwing::fsi::fluid
