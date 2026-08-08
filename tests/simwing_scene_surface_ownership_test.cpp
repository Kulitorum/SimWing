#include "fluid/scene_surface_ownership.h"

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

Scene triangleScene(const double farCoordinate = 2.7,
                    const double zMeters = 1.5) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-ownership";
    scene.metadata.exporterVersion = "scene-surface-ownership-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.1, 1.1, zMeters}},
        {11, {farCoordinate, 1.1, zMeters}},
        {12, {1.1, farCoordinate, zMeters}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0},
           {farCoordinate - 1.1, 0.0},
           {0.0, farCoordinate - 1.1}}},
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
              surface.definition,
              state,
              grid(),
              candidates,
              intersections)) {}
};

SceneFluidPatchOwnership ownership(const Pipeline& pipeline) {
    return ownSceneFluidSurfacePatches(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches);
}

void testOrdinaryCellOwnership() {
    Pipeline pipeline(triangleScene());
    const auto first = ownership(pipeline);
    const auto second = ownership(pipeline);
    check(first == second && first.fingerprint != 0
              && first.patchSetFingerprint == pipeline.patches.fingerprint
              && first.cellPatches.size() == 3
              && first.facePatches.empty()
              && first.pointContactPatchCount == 0
              && first.segmentContactPatchCount == 0,
          "scene patch ownership: ordinary clipped areas remain deterministic cell owners");
    checkNear(first.ownedAreaSquareMeters, 1.28, 3.0e-15,
              "scene patch ownership: ordinary cell ownership partitions total area");
    for (std::size_t index = 0; index < first.cellPatches.size(); ++index) {
        const auto& owned = first.cellPatches[index];
        check(owned.sourcePatchIndex == index
                  && owned.triangleId == 500
                  && owned.negativeSideRegionId == 1
                  && owned.positiveSideRegionId == 2,
              "scene patch ownership: cell patch preserves source and authored sides");
    }
    validateSceneFluidPatchOwnership(
        first,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches);
}

void testUniqueFaceOwnershipAndWinding() {
    Pipeline boundary(triangleScene(2.7, 2.0));
    const auto owned = ownership(boundary);
    check(owned.cellPatches.empty() && owned.facePatches.size() == 3
              && owned.pointContactPatchCount == 0,
          "scene patch ownership: duplicated grid-plane areas become three unique faces");
    checkNear(owned.ownedAreaSquareMeters, 1.28, 3.0e-15,
              "scene patch ownership: paired face ownership removes area duplication");
    for (const auto& face : owned.facePatches) {
        check(face.axis == GridFaceAxis::Z && face.k == 2
                  && face.lowerCellSourcePatchIndex
                      < face.upperCellSourcePatchIndex
                  && face.triangleId == 500
                  && face.negativeSideRegionId == 1
                  && face.positiveSideRegionId == 2
                  && face.triangleNormalAxisSign == 1,
              "scene patch ownership: face location, pair, sides, and +Z winding survive");
    }

    Scene reversedScene = triangleScene(2.7, 2.0);
    std::swap(reversedScene.triangles.front().vertexIds[1],
              reversedScene.triangles.front().vertexIds[2]);
    Pipeline reversed(std::move(reversedScene));
    const auto reversedOwned = ownership(reversed);
    check(reversedOwned.facePatches.size() == 3
              && reversedOwned.facePatches.front().triangleNormalAxisSign == -1
              && reversedOwned.facePatches.front().negativeSideRegionId == 1
              && reversedOwned.facePatches.front().positiveSideRegionId == 2,
          "scene patch ownership: reversed winding is explicit without swapping authored sides");
}

void testContactAndPeriodicBoundaryRejection() {
    Pipeline touching(triangleScene(2.9));
    const auto touchingOwned = ownership(touching);
    check(touchingOwned.cellPatches.size() == 3
              && touchingOwned.pointContactPatchCount == 1
              && touchingOwned.segmentContactPatchCount == 0,
          "scene patch ownership: zero-area contact remains diagnostic, never an owner");
    checkNear(touchingOwned.ownedAreaSquareMeters, 1.62, 4.0e-15,
              "scene patch ownership: point contact cannot duplicate physical area");

    Pipeline domainBoundary(triangleScene(2.7, 0.0));
    expectInvalid(
        [&] { static_cast<void>(ownership(domainBoundary)); },
        "scene patch ownership: periodic-domain boundary needs explicit image ownership");
}

void testTransactionalLimitsAndValidation() {
    Pipeline pipeline(triangleScene());
    SceneFluidPatchOwnershipLimits limits;
    limits.maximumOwnedAreaPatches = 2;
    expectLimited(
        [&] { static_cast<void>(ownSceneFluidSurfacePatches(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            limits)); },
        "scene patch ownership: unique area-owner count is bounded");
    limits = {};
    limits.maximumOwnershipBytes =
        2 * sizeof(SceneFluidOwnedCellPatch);
    expectLimited(
        [&] { static_cast<void>(ownSceneFluidSurfacePatches(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            limits)); },
        "scene patch ownership: owner storage is byte-bounded");

    const auto accepted = ownership(pipeline);
    auto corrupt = accepted;
    corrupt.cellPatches.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidPatchOwnership(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches); },
        "scene patch ownership: corruption is rejected transactionally");
    validateSceneFluidPatchOwnership(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches);
}

} // namespace

int main() {
    testOrdinaryCellOwnership();
    testUniqueFaceOwnershipAndWinding();
    testContactAndPeriodicBoundaryRejection();
    testTransactionalLimitsAndValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-ownership test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-ownership tests passed\n");
    return 0;
}
