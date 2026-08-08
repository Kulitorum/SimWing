#include "fluid/scene_surface_face_loops.h"

#include <cmath>
#include <cstdio>
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

Scene tetrahedronScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-loops";
    scene.metadata.exporterVersion = "scene-surface-face-loops-test/1";
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

Scene openTriangleScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-loops-open";
    scene.metadata.exporterVersion = "scene-surface-face-loops-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.2, 1.2}},
        {11, {2.8, 1.2, 1.2}},
        {12, {1.2, 1.8, 1.8}},
    };
    scene.fabricMaterials = {material()};
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {0.0, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
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
    SceneFluidFaceChainSet chains;

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
              patches, ownership, crossings, topology)),
          chains(buildSceneFluidFaceChains(
              surface.definition, state, grid(), candidates, intersections,
              patches, ownership, crossings, topology, graph)) {}
};

SceneFluidFaceLoopSet loops(
    const Pipeline& pipeline,
    const SceneFluidFaceLoopSettings& settings = {},
    const SceneFluidFaceLoopLimits& limits = {}) {
    return buildSceneFluidFaceLoops(
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
        pipeline.chains,
        settings,
        limits);
}

void testOrientedClosedLoopGeometry() {
    Pipeline pipeline(tetrahedronScene());
    const auto first = loops(pipeline);
    const auto second = loops(pipeline);
    check(first == second && first.fingerprint != 0
              && first.unresolvedOpenChainCount == 0
              && first.loops.size() == 1,
          "scene face loops: closed chain geometry is deterministic");
    const auto& loop = first.loops.front();
    check(loop.stableId == pipeline.chains.chains.front().stableId
              && loop.chainIndex == 0
              && loop.enclosedRegionId == 2
              && loop.exteriorRegionId == 1
              && !loop.positiveSideIsInterior
              && loop.signedAreaSquareMeters < 0.0,
          "scene face loops: signed winding identifies the enclosed authored region");
    checkNear(loop.areaSquareMeters, 0.045, 2.0e-15,
              "scene face loops: section area is analytic");
    checkNear(first.summedLoopAreaSquareMeters, 0.045, 2.0e-15,
              "scene face loops: loop area is counted once");
    checkNear(loop.centroidMeters.x, 2.0, 0.0,
              "scene face loops: centroid stays on its MAC face");
    checkNear(loop.centroidMeters.y, 1.5, 2.0e-15,
              "scene face loops: centroid Y is analytic");
    checkNear(loop.centroidMeters.z, 1.45, 2.0e-15,
              "scene face loops: centroid Z is analytic");
    validateSceneFluidFaceLoops(
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
        pipeline.graph,
        pipeline.chains);
}

void testOpenChainsRemainUnresolved() {
    Pipeline pipeline(openTriangleScene());
    const auto result = loops(pipeline);
    check(result.loops.empty()
              && result.unresolvedOpenChainCount == 1
              && result.summedLoopAreaSquareMeters == 0.0,
          "scene face loops: open chains do not invent sealed area");
}

void testSettingsLimitsAndTransactionalValidation() {
    Pipeline pipeline(tetrahedronScene());
    SceneFluidFaceLoopSettings settings;
    settings.minimumAbsoluteAreaSquareMeters = 0.05;
    expectInvalid(
        [&] { static_cast<void>(loops(pipeline, settings)); },
        "scene face loops: configured degenerate-area floor is enforced");
    settings = {};
    settings.intersectionToleranceMeters = 1.0e-3;
    expectInvalid(
        [&] { static_cast<void>(loops(pipeline, settings)); },
        "scene face loops: over-broad intersection tolerance is rejected");

    SceneFluidFaceLoopLimits limits;
    limits.maximumLoops = 0;
    expectLimited(
        [&] { static_cast<void>(loops(pipeline, {}, limits)); },
        "scene face loops: loop count is bounded");
    limits = {};
    limits.maximumLoopBytes = sizeof(SceneFluidFaceLoop) - 1;
    expectLimited(
        [&] { static_cast<void>(loops(pipeline, {}, limits)); },
        "scene face loops: loop storage is byte-bounded");

    const auto accepted = loops(pipeline);
    auto corrupt = accepted;
    corrupt.loops.front().enclosedRegionId = 1;
    expectInvalid(
        [&] { validateSceneFluidFaceLoops(
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
            pipeline.graph,
            pipeline.chains); },
        "scene face loops: corrupt region classification is rejected");
    validateSceneFluidFaceLoops(
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
        pipeline.graph,
        pipeline.chains);
}

} // namespace

int main() {
    testOrientedClosedLoopGeometry();
    testOpenChainsRemainUnresolved();
    testSettingsLimitsAndTransactionalValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-face-loop test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-face-loop tests passed\n");
    return 0;
}
