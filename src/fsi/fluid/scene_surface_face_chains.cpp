#include "fluid/scene_surface_face_chains.h"

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
constexpr std::size_t missingReference =
    std::numeric_limits<std::size_t>::max();

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
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

Vec3 faceNormal(const GridFaceAxis axis) {
    if (axis == GridFaceAxis::X) {
        return {1.0, 0.0, 0.0};
    }
    if (axis == GridFaceAxis::Y) {
        return {0.0, 1.0, 0.0};
    }
    return {0.0, 0.0, 1.0};
}

struct DirectedSegment {
    std::size_t fromNode = 0;
    std::size_t toNode = 0;
};

DirectedSegment directedSegment(
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceGraphSegment& segment,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidActiveFace& face) {
    const auto& crossing = crossings.crossings[segment.crossingIndex];
    const Vec3 preferredTangent = cross(
        crossing.negativeToPositiveDirectionInFace,
        faceNormal(face.axis));
    const Vec3 actualTangent = subtract(
        graph.nodes[segment.nodeIndices[1]].positionMeters,
        graph.nodes[segment.nodeIndices[0]].positionMeters);
    const double alignment = dot(actualTangent, preferredTangent);
    if (!std::isfinite(alignment) || alignment == 0.0) {
        throw std::invalid_argument(
            "scene fluid face-chain segment has no winding direction");
    }
    if (alignment > 0.0) {
        return {segment.nodeIndices[0], segment.nodeIndices[1]};
    }
    return {segment.nodeIndices[1], segment.nodeIndices[0]};
}

std::uint64_t chainStableId(
    const SceneFluidFaceChain& chain,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const std::vector<std::size_t>& segmentReferences) {
    Fingerprint fingerprint;
    fingerprint.integer(sceneFluidFaceChainVersion);
    fingerprint.integer(
        topology.activeFaces[chain.activeFaceIndex].stableId);
    fingerprint.enumeration(chain.kind);
    fingerprint.integer(chain.negativeSideRegionId);
    fingerprint.integer(chain.positiveSideRegionId);
    fingerprint.integer(static_cast<std::uint64_t>(
        chain.segmentReferenceCount));
    for (std::size_t offset = 0;
         offset < chain.segmentReferenceCount; ++offset) {
        const std::size_t segmentIndex = segmentReferences[
            chain.firstSegmentReference + offset];
        fingerprint.integer(graph.segments[segmentIndex].stableId);
    }
    return fingerprint.value();
}

std::uint64_t chainsFingerprint(const SceneFluidFaceChainSet& chains) {
    Fingerprint fingerprint;
    fingerprint.integer(chains.version);
    fingerprint.integer(chains.surfaceDefinitionFingerprint);
    fingerprint.integer(chains.surfaceStateFingerprint);
    fingerprint.integer(chains.faceGraphFingerprint);
    fingerprint.integer(chains.structureDefinitionFingerprint);
    fingerprint.integer(chains.acceptedStepCount);
    fingerprint.real(chains.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(chains.openChainCount));
    fingerprint.integer(static_cast<std::uint64_t>(chains.closedChainCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        chains.openingEndpointCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        chains.gridBoundaryEndpointCount));
    fingerprint.integer(static_cast<std::uint64_t>(chains.chains.size()));
    for (const auto& chain : chains.chains) {
        fingerprint.integer(chain.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(chain.activeFaceIndex));
        fingerprint.enumeration(chain.kind);
        fingerprint.integer(chain.negativeSideRegionId);
        fingerprint.integer(chain.positiveSideRegionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            chain.firstNodeReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            chain.nodeReferenceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            chain.firstSegmentReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            chain.segmentReferenceCount));
        fingerprint.real(chain.lengthMeters);
        fingerprint.integer(chain.endpointFaceBoundaryMasks[0]);
        fingerprint.integer(chain.endpointFaceBoundaryMasks[1]);
        fingerprint.integer(static_cast<std::uint8_t>(
            chain.endpointOnAuthoredOpening[0] ? 1 : 0));
        fingerprint.integer(static_cast<std::uint8_t>(
            chain.endpointOnAuthoredOpening[1] ? 1 : 0));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        chains.nodeReferences.size()));
    for (const std::size_t reference : chains.nodeReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        chains.segmentReferences.size()));
    for (const std::size_t reference : chains.segmentReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    return fingerprint.value();
}

SceneFluidFaceChainSet buildChains(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph,
    const SceneFluidFaceChainLimits& limits) {
    SceneFluidFaceChainSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.faceGraphFingerprint = graph.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.segmentReferences.reserve(graph.segments.size());

    std::vector<bool> globallyVisited(graph.segments.size(), false);
    std::set<std::uint64_t> chainStableIds;
    for (const auto& range : graph.faceRanges) {
        const auto& face = topology.activeFaces[range.activeFaceIndex];
        std::vector<std::size_t> incoming(range.nodeCount, missingReference);
        std::vector<std::size_t> outgoing(range.nodeCount, missingReference);
        std::vector<DirectedSegment> directions(range.segmentCount);

        for (std::size_t offset = 0; offset < range.segmentCount; ++offset) {
            const std::size_t segmentIndex = range.firstSegment + offset;
            const auto& segment = graph.segments[segmentIndex];
            const DirectedSegment direction = directedSegment(
                graph, segment, crossings, face);
            directions[offset] = direction;
            const std::size_t from = direction.fromNode - range.firstNode;
            const std::size_t to = direction.toNode - range.firstNode;
            if (from >= range.nodeCount || to >= range.nodeCount
                || outgoing[from] != missingReference
                || incoming[to] != missingReference) {
                throw std::invalid_argument(
                    "scene fluid face graph has a branch or conflicting winding");
            }
            outgoing[from] = segmentIndex;
            incoming[to] = segmentIndex;
        }

        auto appendChain = [&](const std::size_t startNode,
                               const bool expectedClosed) {
            if (result.chains.size() == limits.maximumChains) {
                throw std::length_error(
                    "scene fluid face chains exceed their count limit");
            }
            SceneFluidFaceChain chain;
            chain.activeFaceIndex = range.activeFaceIndex;
            chain.firstNodeReference = result.nodeReferences.size();
            chain.firstSegmentReference = result.segmentReferences.size();
            std::size_t currentNode = startNode;
            result.nodeReferences.push_back(currentNode);

            while (true) {
                const std::size_t localNode = currentNode - range.firstNode;
                const std::size_t segmentIndex = outgoing[localNode];
                if (segmentIndex == missingReference) {
                    chain.kind = SceneFluidFaceChainKind::Open;
                    break;
                }
                if (globallyVisited[segmentIndex]) {
                    throw std::invalid_argument(
                        "scene fluid face chain revisits a segment");
                }
                globallyVisited[segmentIndex] = true;
                result.segmentReferences.push_back(segmentIndex);
                const auto& crossing = crossings.crossings[
                    graph.segments[segmentIndex].crossingIndex];
                if (chain.negativeSideRegionId == invalidStableId) {
                    chain.negativeSideRegionId =
                        crossing.negativeSideRegionId;
                    chain.positiveSideRegionId =
                        crossing.positiveSideRegionId;
                } else if (chain.negativeSideRegionId
                               != crossing.negativeSideRegionId
                           || chain.positiveSideRegionId
                               != crossing.positiveSideRegionId) {
                    throw std::invalid_argument(
                        "scene fluid face chain changes authored region pair");
                }
                chain.lengthMeters += crossing.lengthMeters;
                const DirectedSegment& direction = directions[
                    segmentIndex - range.firstSegment];
                currentNode = direction.toNode;
                if (currentNode == startNode) {
                    chain.kind = SceneFluidFaceChainKind::Closed;
                    break;
                }
                result.nodeReferences.push_back(currentNode);
            }

            chain.nodeReferenceCount = result.nodeReferences.size()
                - chain.firstNodeReference;
            chain.segmentReferenceCount = result.segmentReferences.size()
                - chain.firstSegmentReference;
            if (chain.segmentReferenceCount == 0
                || !std::isfinite(chain.lengthMeters)
                || !(chain.lengthMeters > 0.0)
                || (expectedClosed
                    != (chain.kind == SceneFluidFaceChainKind::Closed))) {
                throw std::invalid_argument(
                    "scene fluid face chain is incomplete or misclassified");
            }
            if (result.nodeReferences.size() > limits.maximumNodeReferences
                || result.segmentReferences.size()
                    > limits.maximumSegmentReferences) {
                throw std::length_error(
                    "scene fluid face-chain references exceed their limits");
            }
            if (chain.kind == SceneFluidFaceChainKind::Open) {
                ++result.openChainCount;
                const std::size_t firstNode = result.nodeReferences[
                    chain.firstNodeReference];
                const std::size_t lastNode = result.nodeReferences[
                    chain.firstNodeReference + chain.nodeReferenceCount - 1];
                const auto& first = graph.nodes[firstNode];
                const auto& last = graph.nodes[lastNode];
                chain.endpointFaceBoundaryMasks = {
                    first.faceBoundaryMask, last.faceBoundaryMask};
                chain.endpointOnAuthoredOpening = {
                    first.authoredOpeningBoundary,
                    last.authoredOpeningBoundary};
                for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
                    result.gridBoundaryEndpointCount +=
                        chain.endpointFaceBoundaryMasks[endpoint]
                            != FaceBoundaryNone ? 1 : 0;
                    result.openingEndpointCount +=
                        chain.endpointOnAuthoredOpening[endpoint] ? 1 : 0;
                }
            } else {
                ++result.closedChainCount;
            }
            chain.stableId = chainStableId(
                chain, topology, graph, result.segmentReferences);
            if (!chainStableIds.insert(chain.stableId).second) {
                throw std::invalid_argument(
                    "scene fluid face-chain stable-ID collision");
            }
            result.chains.push_back(chain);
        };

        std::vector<std::size_t> starts;
        for (std::size_t localNode = 0;
             localNode < range.nodeCount; ++localNode) {
            if (incoming[localNode] == missingReference
                && outgoing[localNode] != missingReference) {
                starts.push_back(range.firstNode + localNode);
            }
        }
        std::sort(starts.begin(), starts.end(),
                  [&](const std::size_t first, const std::size_t second) {
                      return graph.nodes[first].stableId
                          < graph.nodes[second].stableId;
                  });
        for (const std::size_t start : starts) {
            const std::size_t firstSegment = outgoing[start - range.firstNode];
            if (!globallyVisited[firstSegment]) {
                appendChain(start, false);
            }
        }

        while (true) {
            std::size_t cycleStart = missingReference;
            for (std::size_t localNode = 0;
                 localNode < range.nodeCount; ++localNode) {
                const std::size_t segmentIndex = outgoing[localNode];
                if (segmentIndex != missingReference
                    && !globallyVisited[segmentIndex]
                    && (cycleStart == missingReference
                        || graph.nodes[range.firstNode + localNode].stableId
                            < graph.nodes[cycleStart].stableId)) {
                    cycleStart = range.firstNode + localNode;
                }
            }
            if (cycleStart == missingReference) {
                break;
            }
            appendChain(cycleStart, true);
        }
    }

    if (std::find(globallyVisited.begin(), globallyVisited.end(), false)
        != globallyVisited.end()) {
        throw std::invalid_argument(
            "scene fluid face chains do not cover every graph segment");
    }
    std::size_t chainBytes = 0;
    std::size_t nodeReferenceBytes = 0;
    std::size_t segmentReferenceBytes = 0;
    std::size_t referenceBytes = 0;
    std::size_t totalBytes = 0;
    if (!checkedMultiply(result.chains.size(),
                         sizeof(SceneFluidFaceChain), chainBytes)
        || !checkedMultiply(result.nodeReferences.size(),
                            sizeof(std::size_t), nodeReferenceBytes)
        || !checkedMultiply(result.segmentReferences.size(),
                            sizeof(std::size_t), segmentReferenceBytes)
        || !checkedAdd(nodeReferenceBytes, segmentReferenceBytes,
                       referenceBytes)
        || !checkedAdd(chainBytes, referenceBytes, totalBytes)
        || totalBytes > limits.maximumChainBytes) {
        throw std::length_error(
            "scene fluid face chains exceed their storage limit");
    }
    result.fingerprint = chainsFingerprint(result);
    return result;
}

} // namespace

SceneFluidFaceChainSet buildSceneFluidFaceChains(
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
    const SceneFluidFaceChainLimits& limits) {
    validateSceneFluidFaceGraph(
        graph, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology);
    SceneFluidFaceChainSet result = buildChains(
        surface, state, crossings, topology, graph, limits);
    validateSceneFluidFaceChains(
        result, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology, graph);
    return result;
}

void validateSceneFluidFaceChains(
    const SceneFluidFaceChainSet& chains,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopology& topology,
    const SceneFluidFaceGraph& graph) {
    validateSceneFluidFaceGraph(
        graph, surface, state, grid, candidates, intersections, patches,
        ownership, crossings, topology);
    if (chains.version != sceneFluidFaceChainVersion
        || chains.fingerprint == 0
        || chains.surfaceDefinitionFingerprint != surface.fingerprint
        || chains.surfaceStateFingerprint != state.fingerprint
        || chains.faceGraphFingerprint != graph.fingerprint
        || chains.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || chains.acceptedStepCount != state.acceptedStepCount
        || chains.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid face-chain identity is invalid");
    }
    const SceneFluidFaceChainSet expected = buildChains(
        surface,
        state,
        crossings,
        topology,
        graph,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (chains != expected) {
        throw std::invalid_argument(
            "scene fluid face chains do not match their source graph");
    }
}

} // namespace simwing::fsi::fluid
