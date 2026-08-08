#include "scene_fluid_grid_epoch.h"

#include <limits>
#include <stdexcept>
#include <type_traits>

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

template<typename Value>
bool addVectorBytes(const std::vector<Value>& values,
                    std::size_t& total) {
    if (values.size()
        > std::numeric_limits<std::size_t>::max() / sizeof(Value)) {
        return false;
    }
    std::size_t next = 0;
    if (!checkedAdd(total, values.size() * sizeof(Value), next)) {
        return false;
    }
    total = next;
    return true;
}

std::size_t ownedStorageBytes(const SceneFluidGridEpoch& epoch) {
    std::size_t total = 0;
    const bool valid =
        addVectorBytes(epoch.candidates.triangleBounds, total)
        && addVectorBytes(epoch.candidates.candidates, total)
        && addVectorBytes(
            epoch.intersections.triangleIntersectionCounts, total)
        && addVectorBytes(epoch.intersections.intersections, total)
        && addVectorBytes(epoch.patches.patches, total)
        && addVectorBytes(epoch.patches.vertices, total)
        && addVectorBytes(epoch.ownership.cellPatches, total)
        && addVectorBytes(epoch.ownership.facePatches, total)
        && addVectorBytes(epoch.crossings.crossings, total)
        && addVectorBytes(epoch.faceTopology.activeFaces, total)
        && addVectorBytes(epoch.faceTopology.crossingReferences, total)
        && addVectorBytes(epoch.faceTopology.coplanarPatchReferences, total)
        && addVectorBytes(epoch.faceGraph.faceRanges, total)
        && addVectorBytes(epoch.faceGraph.nodes, total)
        && addVectorBytes(epoch.faceGraph.segments, total)
        && addVectorBytes(
            epoch.faceGraph.incidentSegmentReferences, total)
        && addVectorBytes(epoch.faceChains.chains, total)
        && addVectorBytes(epoch.faceChains.nodeReferences, total)
        && addVectorBytes(epoch.faceChains.segmentReferences, total)
        && addVectorBytes(epoch.faceLoops.loops, total)
        && addVectorBytes(epoch.facePartitions.loopContainment, total)
        && addVectorBytes(epoch.facePartitions.partitions, total)
        && addVectorBytes(epoch.facePartitions.loopReferences, total)
        && addVectorBytes(epoch.facePartitions.regionAreas, total)
        && addVectorBytes(epoch.quadrature.points, total);
    if (!valid) {
        throw std::length_error(
            "scene fluid grid-epoch storage size overflows");
    }
    return total;
}

std::uint64_t epochFingerprint(const SceneFluidGridEpoch& epoch) {
    Fingerprint fingerprint;
    fingerprint.integer(epoch.version);
    fingerprint.integer(static_cast<std::uint64_t>(
        epoch.ownedStorageBytes));
    fingerprint.integer(epoch.candidates.fingerprint);
    fingerprint.integer(epoch.intersections.fingerprint);
    fingerprint.integer(epoch.patches.fingerprint);
    fingerprint.integer(epoch.ownership.fingerprint);
    fingerprint.integer(epoch.crossings.fingerprint);
    fingerprint.integer(epoch.faceTopology.fingerprint);
    fingerprint.integer(epoch.faceGraph.fingerprint);
    fingerprint.integer(epoch.faceChains.fingerprint);
    fingerprint.integer(epoch.faceLoops.fingerprint);
    fingerprint.integer(epoch.facePartitions.fingerprint);
    fingerprint.integer(epoch.quadrature.fingerprint);
    return fingerprint.value();
}

} // namespace

SceneFluidGridEpoch buildSceneFluidGridEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpochSettings& settings,
    const SceneFluidGridEpochLimits& limits) {
    SceneFluidGridEpoch result;
    result.candidates = fluid::buildSceneFluidGridCandidates(
        surface, state, grid, settings.candidates, limits.candidates);
    result.intersections = fluid::intersectSceneFluidSurfaceWithGrid(
        surface, state, grid, result.candidates,
        settings.intersections, limits.intersections);
    result.patches = fluid::clipSceneFluidSurfaceToCells(
        surface, state, grid, result.candidates, result.intersections,
        limits.patches);
    result.ownership = fluid::ownSceneFluidSurfacePatches(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, limits.ownership);
    result.crossings = fluid::buildSceneFluidFaceCrossings(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, limits.crossings);
    result.faceTopology = fluid::buildSceneFluidFaceTopology(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, result.crossings,
        limits.faceTopology);
    result.faceGraph = fluid::buildSceneFluidFaceGraph(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, result.crossings,
        result.faceTopology, settings.faceGraph, limits.faceGraph);
    result.faceChains = fluid::buildSceneFluidFaceChains(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, result.crossings,
        result.faceTopology, result.faceGraph, limits.faceChains);
    result.faceLoops = fluid::buildSceneFluidFaceLoops(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, result.crossings,
        result.faceTopology, result.faceGraph, result.faceChains,
        settings.faceLoops, limits.faceLoops);
    result.facePartitions = fluid::buildSceneFluidFacePartitions(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, result.crossings,
        result.faceTopology, result.faceGraph, result.faceChains,
        result.faceLoops, settings.facePartitions, limits.facePartitions);
    result.quadrature = buildSceneFluidQuadrature(
        surface, state, grid, result.candidates, result.intersections,
        result.patches, result.ownership, transfer);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumEpochBytes) {
        throw std::length_error(
            "scene fluid grid epoch exceeds its aggregate byte limit");
    }
    result.fingerprint = epochFingerprint(result);
    validateSceneFluidGridEpoch(result, surface, state, grid, transfer);
    return result;
}

void validateSceneFluidGridEpoch(
    const SceneFluidGridEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer) {
    if (epoch.version != sceneFluidGridEpochVersion
        || epoch.fingerprint == 0) {
        throw std::invalid_argument(
            "scene fluid grid-epoch identity is invalid");
    }
    fluid::validateSceneFluidGridCandidates(
        epoch.candidates, surface, state, grid);
    fluid::validateSceneFluidGridIntersections(
        epoch.intersections, surface, state, grid, epoch.candidates);
    fluid::validateSceneFluidGridPatches(
        epoch.patches, surface, state, grid, epoch.candidates,
        epoch.intersections);
    fluid::validateSceneFluidPatchOwnership(
        epoch.ownership, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches);
    fluid::validateSceneFluidFaceCrossings(
        epoch.crossings, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership);
    fluid::validateSceneFluidFaceTopology(
        epoch.faceTopology, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership,
        epoch.crossings);
    fluid::validateSceneFluidFaceGraph(
        epoch.faceGraph, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership,
        epoch.crossings, epoch.faceTopology);
    fluid::validateSceneFluidFaceChains(
        epoch.faceChains, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership,
        epoch.crossings, epoch.faceTopology, epoch.faceGraph);
    fluid::validateSceneFluidFaceLoops(
        epoch.faceLoops, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership,
        epoch.crossings, epoch.faceTopology, epoch.faceGraph,
        epoch.faceChains);
    fluid::validateSceneFluidFacePartitions(
        epoch.facePartitions, surface, state, grid, epoch.candidates,
        epoch.intersections, epoch.patches, epoch.ownership,
        epoch.crossings, epoch.faceTopology, epoch.faceGraph,
        epoch.faceChains, epoch.faceLoops);
    validateSceneFluidQuadratureDefinition(epoch.quadrature);
    if (epoch.quadrature.surfaceDefinitionFingerprint
            != surface.fingerprint
        || epoch.quadrature.surfaceStateFingerprint != state.fingerprint
        || epoch.quadrature.patchOwnershipFingerprint
            != epoch.ownership.fingerprint
        || epoch.quadrature.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || epoch.quadrature.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || epoch.quadrature.acceptedStepCount != state.acceptedStepCount
        || epoch.quadrature.simulationTimeSeconds
            != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid grid-epoch quadrature binding is invalid");
    }
    const auto expectedQuadrature = buildSceneFluidQuadrature(
        surface, state, grid, epoch.candidates, epoch.intersections,
        epoch.patches, epoch.ownership, transfer);
    if (epoch.quadrature != expectedQuadrature
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "scene fluid grid-epoch payload is invalid");
    }
}

} // namespace simwing::fsi
