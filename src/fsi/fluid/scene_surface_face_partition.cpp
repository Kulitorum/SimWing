#include "fluid/scene_surface_face_partition.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
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
        f.integer(item.rootExteriorRegionId);
        f.integer(static_cast<std::uint64_t>(item.firstLoopReference));
        f.integer(static_cast<std::uint64_t>(item.loopReferenceCount));
        f.integer(static_cast<std::uint64_t>(item.firstRegionArea));
        f.integer(static_cast<std::uint64_t>(item.regionAreaCount));
        f.real(item.faceAreaSquareMeters);
        f.real(item.assignedAreaSquareMeters);
        f.real(item.areaResidualSquareMeters);
    }
    f.integer(static_cast<std::uint64_t>(value.loopReferences.size()));
    for (auto item : value.loopReferences)
        f.integer(static_cast<std::uint64_t>(item));
    f.integer(static_cast<std::uint64_t>(value.regionAreas.size()));
    for (const auto& item : value.regionAreas) {
        f.integer(item.regionId);
        f.real(item.areaSquareMeters);
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

    std::vector<std::vector<std::size_t>> loopsByFace(topology.activeFaces.size());
    std::vector<std::vector<Point2>> polygons(loops.loops.size());
    for (std::size_t index = 0; index < loops.loops.size(); ++index) {
        result.loopContainment[index].loopIndex = index;
        const auto& loop = loops.loops[index];
        loopsByFace[loop.activeFaceIndex].push_back(index);
        const auto& chain = chains.chains[loop.chainIndex];
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

    std::vector<std::size_t> openByFace(topology.activeFaces.size(), 0);
    for (const auto& chain : chains.chains)
        if (chain.kind == SceneFluidFaceChainKind::Open)
            ++openByFace[chain.activeFaceIndex];

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
        if (faceLoops.empty() || openByFace[faceIndex] != 0
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
        partition.rootExteriorRegionId = rootExterior;
        partition.firstLoopReference = result.loopReferences.size();
        partition.loopReferenceCount = faceLoops.size();
        result.loopReferences.insert(result.loopReferences.end(),
                                     faceLoops.begin(), faceLoops.end());
        std::map<StableId, double> areas;
        partition.faceAreaSquareMeters = faceArea(grid, face.axis);
        areas[rootExterior] = partition.faceAreaSquareMeters;
        for (const std::size_t loopIndex : faceLoops) {
            const auto& loop = loops.loops[loopIndex];
            areas[loop.exteriorRegionId] -= loop.areaSquareMeters;
            areas[loop.enclosedRegionId] += loop.areaSquareMeters;
        }
        partition.firstRegionArea = result.regionAreas.size();
        for (auto [region, area] : areas) {
            if (area < -settings.areaClosureToleranceSquareMeters)
                throw std::invalid_argument("scene fluid face region area is negative");
            if (area < 0.0) area = 0.0;
            result.regionAreas.push_back({region, area});
            partition.assignedAreaSquareMeters += area;
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

    std::size_t referenceCount = 0, partitionBytes = 0, containmentBytes = 0;
    std::size_t loopReferenceBytes = 0, regionBytes = 0, first = 0, second = 0, total = 0;
    if (!checkedAdd(result.loopReferences.size(), result.regionAreas.size(), referenceCount)
        || referenceCount > limits.maximumReferences
        || !checkedMultiply(result.partitions.size(), sizeof(SceneFluidFacePartition), partitionBytes)
        || !checkedMultiply(result.loopContainment.size(), sizeof(SceneFluidFaceLoopContainment), containmentBytes)
        || !checkedMultiply(result.loopReferences.size(), sizeof(std::size_t), loopReferenceBytes)
        || !checkedMultiply(result.regionAreas.size(), sizeof(SceneFluidFaceRegionArea), regionBytes)
        || !checkedAdd(partitionBytes, containmentBytes, first)
        || !checkedAdd(loopReferenceBytes, regionBytes, second)
        || !checkedAdd(first, second, total)
        || total > limits.maximumPartitionBytes)
        throw std::length_error("scene fluid face partitions exceed storage limits");
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
