#include "fluid/scene_surface_face_graph.h"

#include <cstdio>
#include <stdexcept>

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

Scene quadScene(const bool withOpening = false) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-graph";
    scene.metadata.exporterVersion = "scene-surface-face-graph-test/1";
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
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
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

    explicit Pipeline(const bool withOpening = false)
        : scene(quadScene(withOpening)),
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
              patches, ownership, crossings)) {}
};

SceneFluidFaceGraph graph(
    const Pipeline& pipeline,
    const SceneFluidFaceGraphSettings& settings = {},
    const SceneFluidFaceGraphLimits& limits = {}) {
    return buildSceneFluidFaceGraph(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology,
        settings,
        limits);
}

void testTopologicalStitching() {
    Pipeline pipeline;
    const auto first = graph(pipeline);
    const auto second = graph(pipeline);
    check(first == second && first.fingerprint != 0
              && first.faceRanges.size() == 1
              && first.nodes.size() == 3
              && first.segments.size() == 2
              && first.incidentSegmentReferences.size() == 4
              && first.degreeOneNodeCount == 2
              && first.degreeTwoNodeCount == 1
              && first.higherDegreeNodeCount == 0,
          "scene face graph: adjacent triangle segments stitch deterministically");
    check(first.faceRanges.front().nodeCount == 3
              && first.faceRanges.front().segmentCount == 2,
          "scene face graph: face-local ranges own the stitched graph");

    std::size_t sharedNode = first.nodes.size();
    for (std::size_t index = 0; index < first.nodes.size(); ++index) {
        const auto& node = first.nodes[index];
        check(node.stableId != 0
                  && node.kind == SceneFluidFaceNodeKind::SurfaceEdge
                  && node.sourceVertexId == invalidStableId
                  && node.sourceTriangleId == invalidStableId,
              "scene face graph: interior endpoints use stable surface-edge provenance");
        if (node.sourceEdgeVertexIds
            == std::array<StableId, 2>{10, 12}) {
            sharedNode = index;
            check(node.incidentSegmentReferenceCount == 2,
                  "scene face graph: shared authored edge has degree two");
        }
    }
    check(sharedNode < first.nodes.size(),
          "scene face graph: shared authored edge becomes one node");
    if (sharedNode < first.nodes.size()) {
        check(first.segments[0].nodeIndices[0] == sharedNode
                  || first.segments[0].nodeIndices[1] == sharedNode,
              "scene face graph: first triangle reaches the shared node");
        check(first.segments[1].nodeIndices[0] == sharedNode
                  || first.segments[1].nodeIndices[1] == sharedNode,
              "scene face graph: second triangle reaches the shared node");
    }
    validateSceneFluidFaceGraph(
        first,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology);
}

void testOpeningProvenance() {
    Pipeline pipeline(true);
    const auto result = graph(pipeline);
    check(result.authoredOpeningNodeCount == 2,
          "scene face graph: opening vertices and edges remain identifiable");
    bool foundBoundary = false;
    for (const auto& node : result.nodes) {
        foundBoundary = foundBoundary || node.authoredOpeningBoundary;
    }
    check(foundBoundary,
          "scene face graph: opening provenance is attached to graph nodes");
}

void testLimitsSettingsAndTransactionalValidation() {
    Pipeline pipeline;
    SceneFluidFaceGraphLimits limits;
    limits.maximumNodes = 2;
    expectLimited(
        [&] { static_cast<void>(graph(pipeline, {}, limits)); },
        "scene face graph: node count is bounded");
    limits = {};
    limits.maximumSegments = 1;
    expectLimited(
        [&] { static_cast<void>(graph(pipeline, {}, limits)); },
        "scene face graph: segment count is bounded");
    limits = {};
    limits.maximumIncidentReferences = 3;
    expectLimited(
        [&] { static_cast<void>(graph(pipeline, {}, limits)); },
        "scene face graph: incident-reference count is bounded");
    limits = {};
    limits.maximumGraphBytes = sizeof(SceneFluidFaceGraphRange)
        + 3 * sizeof(SceneFluidFaceGraphNode)
        + 2 * sizeof(SceneFluidFaceGraphSegment)
        + 4 * sizeof(std::size_t) - 1;
    expectLimited(
        [&] { static_cast<void>(graph(pipeline, {}, limits)); },
        "scene face graph: graph storage is byte-bounded");
    SceneFluidFaceGraphSettings settings;
    settings.endpointToleranceMeters = 1.0e-3;
    expectInvalid(
        [&] { static_cast<void>(graph(pipeline, settings)); },
        "scene face graph: over-broad endpoint welding is rejected");

    const auto accepted = graph(pipeline);
    auto corrupt = accepted;
    ++corrupt.segments.front().crossingIndex;
    expectInvalid(
        [&] { validateSceneFluidFaceGraph(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            pipeline.ownership,
            pipeline.crossings,
            pipeline.topology); },
        "scene face graph: corrupt connectivity is rejected transactionally");
    validateSceneFluidFaceGraph(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        pipeline.topology);
}

} // namespace

int main() {
    testTopologicalStitching();
    testOpeningProvenance();
    testLimitsSettingsAndTransactionalValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-face-graph test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-face-graph tests passed\n");
    return 0;
}
