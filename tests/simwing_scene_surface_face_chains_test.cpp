#include "fluid/scene_surface_face_chains.h"

#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <utility>

namespace {

using namespace simwing::fsi;
using namespace simwing::fsi::fluid;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

template<typename Callback>
void expectLimited(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected, message);
}

FabricMaterial material() {
    return {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
            0.041, 0.02, 0.0125, 2.5e-12};
}

Scene openQuadScene(const bool withOpening = false) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-chains-open";
    scene.metadata.exporterVersion = "scene-surface-face-chains-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.2, 1.2}},
        {11, {2.8, 1.2, 1.2}},
        {12, {2.8, 1.8, 1.8}},
        {13, {1.2, 1.8, 1.8}},
    };
    scene.fabricMaterials = {material()};
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {1.6, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 12, 13},
         {{{0.0, 0.0}, {1.6, 0.6}, {0.0, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    if (withOpening) {
        scene.openings = {
            {700, {10, 11, 12}, 1, 2, OpeningRole::Intake},
        };
    }
    return scene;
}

Scene closedTetrahedronScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-chains-closed";
    scene.metadata.exporterVersion = "scene-surface-face-chains-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.5, 1.5}},
        {11, {2.8, 1.2, 1.2}},
        {12, {2.8, 1.8, 1.2}},
        {13, {2.8, 1.5, 1.8}},
    };
    scene.fabricMaterials = {material()};
    const std::array<Vec2, 3> chart{{{0.0, 0.0},
                                      {1.0, 0.0},
                                      {0.0, 1.0}}};
    scene.triangles = {
        {500, {10, 12, 11}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
        {501, {10, 11, 13}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
        {502, {10, 13, 12}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
        {503, {11, 12, 13}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
    };
    return scene;
}

Scene threeRegionJunctionScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-surface-face-chains-three-region";
    scene.metadata.exporterVersion = "scene-surface-face-chains-test/2";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell-a"},
        {3, RegionKind::Cell, "cell-b"},
    };
    scene.vertices = {
        {10, {1.5, 1.5, 1.5}},
        {11, {2.5, 1.5, 1.5}},
        {12, {2.5, 1.8, 1.5}},
        {13, {2.5, 1.35, 1.7598076211353315}},
        {14, {2.5, 1.35, 1.2401923788646685}},
    };
    scene.fabricMaterials = {material()};
    const std::array<Vec2, 3> chart{{
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.3}}};
    scene.triangles = {
        {500, {10, 11, 12}, chart,
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 11, 13}, chart,
         2, 3, 100, 901, SurfaceRole::Rib},
        {502, {10, 11, 14}, chart,
         3, 1, 100, 902, SurfaceRole::Skin},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return PeriodicCartesianGrid({4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

struct Pipeline {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceState state;
    SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections;
    SceneFluidGridPatchSet patches;
    SceneFluidPatchOwnership ownership;
    SceneFluidFaceCrossingSet crossings;
    SceneFluidFaceTopology topology;
    SceneFluidFaceGraph graph;

    explicit Pipeline(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          state(captureSceneFluidSurfaceState(
              surface.definition, structureAssembly.mappings, structure)),
          candidates(buildSceneFluidGridCandidates(
              surface.definition, state, grid())),
          intersections(intersectSceneFluidSurfaceWithGrid(
              surface.definition, state, grid(), candidates)),
          patches(clipSceneFluidSurfaceToCells(
              surface.definition, state, grid(), candidates, intersections)),
          ownership(ownSceneFluidSurfacePatches(
              surface.definition, state, grid(), candidates, intersections,
              patches)),
          crossings(buildSceneFluidFaceCrossings(
              surface.definition, state, grid(), candidates, intersections,
              patches, ownership)),
          topology(buildSceneFluidFaceTopology(
              surface.definition, state, grid(), candidates, intersections,
              patches, ownership, crossings)),
          graph(buildSceneFluidFaceGraph(
              surface.definition, state, grid(), candidates, intersections,
              patches, ownership, crossings, topology)) {}
};

SceneFluidFaceChainSet chains(
    const Pipeline& pipeline,
    const SceneFluidFaceChainLimits& limits = {}) {
    return buildSceneFluidFaceChains(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology,
        pipeline.graph,
        limits);
}

void testDirectedOpenChain() {
    Pipeline pipeline(openQuadScene(true));
    const auto first = chains(pipeline);
    const auto second = chains(pipeline);
    check(first == second && first.fingerprint != 0
              && first.openChainCount == 1
              && first.closedChainCount == 0
              && first.chains.size() == 1
              && first.nodeReferences.size() == 3
              && first.segmentReferences.size() == 2,
          "scene face chains: stitched quad becomes one deterministic open chain");
    const auto& chain = first.chains.front();
    check(chain.stableId != 0
              && chain.kind == SceneFluidFaceChainKind::Open
              && chain.negativeSideRegionId == 1
              && chain.positiveSideRegionId == 2
              && chain.nodeReferenceCount == 3
              && chain.segmentReferenceCount == 2,
          "scene face chains: winding and authored region pair direct the chain");
    check(first.openingEndpointCount == 1
              && first.gridBoundaryEndpointCount == 0
              && (chain.endpointOnAuthoredOpening[0]
                  != chain.endpointOnAuthoredOpening[1]),
          "scene face chains: authored opening provenance reaches an open endpoint");
    checkNear(chain.lengthMeters, std::sqrt(0.72), 3.0e-15,
              "scene face chains: open chain length sums each segment once");
    validateSceneFluidFaceChains(
        first,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology,
        pipeline.graph);
}

void testClosedLoop() {
    Pipeline pipeline(closedTetrahedronScene());
    const auto result = chains(pipeline);
    check(result.openChainCount == 0
              && result.closedChainCount == 1
              && result.chains.size() == 1
              && result.nodeReferences.size() == 3
              && result.segmentReferences.size() == 3,
          "scene face chains: closed tetrahedron slice becomes one loop");
    const auto& loop = result.chains.front();
    check(loop.kind == SceneFluidFaceChainKind::Closed
              && loop.negativeSideRegionId == 2
              && loop.positiveSideRegionId == 1
              && loop.nodeReferenceCount == loop.segmentReferenceCount
              && loop.endpointFaceBoundaryMasks
                  == std::array<std::uint8_t, 2>{0, 0}
              && loop.endpointOnAuthoredOpening
                  == std::array<bool, 2>{false, false},
          "scene face chains: closed loops have no synthetic endpoints");
    checkNear(loop.lengthMeters,
              0.3 + 2.0 * std::sqrt(0.1125), 4.0e-15,
              "scene face chains: closed loop perimeter is analytic");
}

void testThreeRegionJunction() {
    Pipeline pipeline(threeRegionJunctionScene());
    check(pipeline.graph.higherDegreeNodeCount == 1
              && pipeline.graph.segments.size() == 3,
          "scene face chains: geometric graph retains one degree-three junction");
    const auto result = chains(pipeline);
    std::set<std::pair<StableId, StableId>> regionPairs;
    bool everyChainEndsAtJunction = true;
    std::size_t junctionNode = pipeline.graph.nodes.size();
    for (std::size_t node = 0; node < pipeline.graph.nodes.size(); ++node) {
        if (pipeline.graph.nodes[node].incidentSegmentReferenceCount == 3) {
            junctionNode = node;
        }
    }
    for (const auto& chain : result.chains) {
        regionPairs.emplace(
            chain.negativeSideRegionId, chain.positiveSideRegionId);
        const std::size_t first = result.nodeReferences[
            chain.firstNodeReference];
        const std::size_t last = result.nodeReferences[
            chain.firstNodeReference + chain.nodeReferenceCount - 1];
        everyChainEndsAtJunction = everyChainEndsAtJunction
            && chain.kind == SceneFluidFaceChainKind::Open
            && chain.segmentReferenceCount == 1
            && chain.nodeReferenceCount == 2
            && (first == junctionNode || last == junctionNode);
    }
    check(junctionNode < pipeline.graph.nodes.size()
              && result.openChainCount == 3
              && result.closedChainCount == 0
              && result.segmentReferences.size() == 3
              && regionPairs == std::set<std::pair<StableId, StableId>>{
                  {1, 2}, {2, 3}, {3, 1}}
              && everyChainEndsAtJunction,
          "scene face chains: one physical junction terminates three region-pair chains");
    validateSceneFluidFaceChains(
        result,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology,
        pipeline.graph);

    auto branchScene = threeRegionJunctionScene();
    for (auto& triangle : branchScene.triangles) {
        triangle.negativeSideRegionId = 1;
        triangle.positiveSideRegionId = 2;
    }
    Pipeline branch(std::move(branchScene));
    expectInvalid(
        [&] { static_cast<void>(chains(branch)); },
        "scene face chains: a degree-three branch within one region pair rejects");
}

void testLimitsAndTransactionalValidation() {
    Pipeline pipeline(openQuadScene());
    SceneFluidFaceChainLimits limits;
    limits.maximumChains = 0;
    expectLimited(
        [&] { static_cast<void>(chains(pipeline, limits)); },
        "scene face chains: chain count is bounded");
    limits = {};
    limits.maximumNodeReferences = 2;
    expectLimited(
        [&] { static_cast<void>(chains(pipeline, limits)); },
        "scene face chains: ordered node references are bounded");
    limits = {};
    limits.maximumSegmentReferences = 1;
    expectLimited(
        [&] { static_cast<void>(chains(pipeline, limits)); },
        "scene face chains: ordered segment references are bounded");
    limits = {};
    limits.maximumChainBytes = sizeof(SceneFluidFaceChain)
        + 5 * sizeof(std::size_t) - 1;
    expectLimited(
        [&] { static_cast<void>(chains(pipeline, limits)); },
        "scene face chains: chain storage is byte-bounded");

    const auto accepted = chains(pipeline);
    auto corrupt = accepted;
    corrupt.chains.front().negativeSideRegionId = 2;
    expectInvalid(
        [&] { validateSceneFluidFaceChains(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            pipeline.ownership,
            pipeline.crossings,
            pipeline.topology,
            pipeline.graph); },
        "scene face chains: corrupt authored-side identity is rejected");
    validateSceneFluidFaceChains(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology,
        pipeline.graph);
}

} // namespace

int main() {
    testDirectedOpenChain();
    testClosedLoop();
    testThreeRegionJunction();
    testLimitsAndTransactionalValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-face-chain test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-face-chain tests passed\n");
    return 0;
}
