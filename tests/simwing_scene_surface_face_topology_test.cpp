#include "fluid/scene_surface_face_topology.h"

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

Scene multiSurfaceScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-face-topology";
    scene.metadata.exporterVersion = "scene-surface-face-topology-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.10, 1.10}},
        {11, {2.8, 1.10, 1.10}},
        {12, {1.2, 1.35, 1.35}},
        {20, {1.2, 1.60, 1.10}},
        {21, {2.8, 1.60, 1.10}},
        {22, {1.2, 1.85, 1.35}},
        {30, {1.2, 1.2, 2.0}},
        {31, {1.8, 1.2, 2.0}},
        {32, {1.2, 1.8, 2.0}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {0.0, 0.25}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {20, 21, 22},
         {{{0.0, 0.0}, {1.6, 0.0}, {0.0, 0.25}}},
         1, 2, 100, 901, SurfaceRole::Skin},
        {502, {30, 31, 32},
         {{{0.0, 0.0}, {0.6, 0.0}, {0.0, 0.6}}},
         1, 2, 100, 902, SurfaceRole::Rib},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return PeriodicCartesianGrid({4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

struct Pipeline {
    Scene scene = multiSurfaceScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structureAssembly = assembleSceneStructure(scene);
    Structure structure{structureAssembly.definition};
    SceneFluidSurfaceState state = captureSceneFluidSurfaceState(
        surface.definition, structureAssembly.mappings, structure);
    SceneFluidGridCandidateSet candidates = buildSceneFluidGridCandidates(
        surface.definition, state, grid());
    SceneFluidGridIntersectionSet intersections =
        intersectSceneFluidSurfaceWithGrid(
            surface.definition, state, grid(), candidates);
    SceneFluidGridPatchSet patches = clipSceneFluidSurfaceToCells(
        surface.definition, state, grid(), candidates, intersections);
    SceneFluidPatchOwnership ownership = ownSceneFluidSurfacePatches(
        surface.definition, state, grid(), candidates, intersections, patches);
    SceneFluidFaceCrossingSet crossings = buildSceneFluidFaceCrossings(
        surface.definition, state, grid(), candidates, intersections, patches,
        ownership);
};

SceneFluidFaceTopology topology(
    const Pipeline& pipeline,
    const SceneFluidFaceTopologyLimits& limits = {}) {
    return buildSceneFluidFaceTopology(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings,
        limits);
}

void testCanonicalFaceIndex() {
    Pipeline pipeline;
    const auto first = topology(pipeline);
    const auto second = topology(pipeline);
    check(first == second && first.fingerprint != 0
              && first.activeFaces.size() == 2
              && first.crossingReferences.size() == 2
              && first.coplanarPatchReferences.size() == 1,
          "scene face topology: mixed face content is indexed deterministically");

    const auto* transverse = first.activeFace(GridFaceAxis::X, 2, 1, 1);
    check(transverse != nullptr && transverse->stableId != 0
              && transverse->crossingReferenceCount == 2
              && transverse->coplanarPatchReferenceCount == 0,
          "scene face topology: multiple sheet crossings share one active face");
    if (transverse != nullptr) {
        const auto references = first.crossingsForFace(*transverse);
        check(references.size() == 2
                  && pipeline.crossings.crossings[references[0]].triangleId
                      != pipeline.crossings.crossings[references[1]].triangleId,
              "scene face topology: crossing references retain separate triangles");
        checkNear(transverse->summedCrossingLengthMeters,
                  2.0 * std::sqrt(0.03125), 3.0e-15,
                  "scene face topology: crossing length is summed without union claims");
    }

    const auto* coplanar = first.activeFace(GridFaceAxis::Z, 1, 1, 2);
    check(coplanar != nullptr && coplanar->stableId != 0
              && coplanar->crossingReferenceCount == 0
              && coplanar->coplanarPatchReferenceCount == 1,
          "scene face topology: coplanar area remains a distinct face reference");
    if (coplanar != nullptr) {
        const auto references = first.coplanarPatchesForFace(*coplanar);
        check(references.size() == 1
                  && pipeline.ownership.facePatches[references[0]].triangleId
                      == 502,
              "scene face topology: coplanar reference resolves to its owner");
        checkNear(coplanar->summedCoplanarAreaSquareMeters,
                  0.18, 2.0e-15,
                  "scene face topology: coplanar area is preserved exactly once");
    }
    check(first.activeFace(GridFaceAxis::Y, 1, 1, 1) == nullptr,
          "scene face topology: inactive faces do not acquire synthetic records");
    validateSceneFluidFaceTopology(
        first,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings);
}

void testLimitsAndTransactionalValidation() {
    Pipeline pipeline;
    SceneFluidFaceTopologyLimits limits;
    limits.maximumActiveFaces = 1;
    expectLimited(
        [&] { static_cast<void>(topology(pipeline, limits)); },
        "scene face topology: active-face count is bounded");
    limits = {};
    limits.maximumReferences = 2;
    expectLimited(
        [&] { static_cast<void>(topology(pipeline, limits)); },
        "scene face topology: flattened reference count is bounded");
    limits = {};
    limits.maximumTopologyBytes =
        2 * sizeof(SceneFluidActiveFace)
        + 3 * sizeof(std::size_t) - 1;
    expectLimited(
        [&] { static_cast<void>(topology(pipeline, limits)); },
        "scene face topology: topology storage is byte-bounded");

    const auto accepted = topology(pipeline);
    auto corrupt = accepted;
    ++corrupt.crossingReferences.front();
    expectInvalid(
        [&] { validateSceneFluidFaceTopology(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            pipeline.ownership,
            pipeline.crossings); },
        "scene face topology: corrupt references are rejected transactionally");
    validateSceneFluidFaceTopology(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.crossings);
}

} // namespace

int main() {
    testCanonicalFaceIndex();
    testLimitsAndTransactionalValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-face-topology test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-face-topology tests passed\n");
    return 0;
}
