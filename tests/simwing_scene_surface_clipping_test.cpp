#include "fluid/scene_surface_clipping.h"

#include <array>
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
    scene.metadata.designChecksum = "sha256:scene-surface-clipping";
    scene.metadata.exporterVersion = "scene-surface-clipping-test/1";
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
              surface.definition, state, grid(), candidates)) {}
};

void testReusableTriangleBoxPrimitive() {
    const std::array<Vec3, 3> triangle{{
        {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {0.0, 1.0, 0.5},
    }};
    const auto clipped = clipSceneFluidTriangleToAxisAlignedBox(
        triangle, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    check(clipped
              && clipped->dimension == SceneFluidPatchDimension::Area
              && clipped->vertices.size() == 3
              && clipped->coincidentBoundaryPlanes == CellBoundaryNone,
          "triangle-box primitive preserves an interior area triangle");
    checkNear(clipped ? clipped->areaSquareMeters : 0.0, 0.5, 0.0,
              "triangle-box primitive preserves analytic area");
    checkNear(clipped ? clipped->centroidMeters.x : 0.0,
              1.0 / 3.0, 0.0,
              "triangle-box primitive preserves analytic centroid");

    const auto disjoint = clipSceneFluidTriangleToAxisAlignedBox(
        triangle, {2.0, 2.0, 2.0}, {3.0, 3.0, 3.0});
    check(!disjoint,
          "triangle-box primitive reports a disjoint box without geometry");
    expectInvalid(
        [&] { static_cast<void>(clipSceneFluidTriangleToAxisAlignedBox(
            triangle, {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0})); },
        "triangle-box primitive rejects a degenerate box");
}

void testExactClippedPatches() {
    Pipeline pipeline(triangleScene());
    const auto first = clipSceneFluidSurfaceToCells(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections);
    const auto second = clipSceneFluidSurfaceToCells(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections);
    check(first == second && first.fingerprint != 0
              && first.surfaceStateFingerprint == pipeline.state.fingerprint
              && first.intersectionSetFingerprint
                  == pipeline.intersections.fingerprint,
          "scene grid clipping: repeated exact patches are deterministic and bound");
    check(first.patches.size() == 3 && first.vertices.size() == 11
              && first.patches[0].dimension == SceneFluidPatchDimension::Area
              && first.patches[0].vertexCount == 5
              && first.patches[1].vertexCount == 3
              && first.patches[2].vertexCount == 3,
          "scene grid clipping: each exact intersection owns one flattened polygon");
    checkNear(first.patches[0].areaSquareMeters, 0.79, 2.0e-15,
              "scene grid clipping: central clipped pentagon area is analytic");
    checkNear(first.patches[1].areaSquareMeters, 0.245, 2.0e-15,
              "scene grid clipping: first boundary triangle area is analytic");
    checkNear(first.patches[2].areaSquareMeters, 0.245, 2.0e-15,
              "scene grid clipping: second boundary triangle area is analytic");
    double totalArea = 0.0;
    for (const auto& patch : first.patches) {
        totalArea += patch.areaSquareMeters;
        const double barycentricSum =
            patch.centroidBarycentricCoordinates[0]
            + patch.centroidBarycentricCoordinates[1]
            + patch.centroidBarycentricCoordinates[2];
        checkNear(barycentricSum, 1.0, 2.0e-15,
                  "scene grid clipping: centroid barycentric coordinates close");
    }
    checkNear(totalArea, 1.28, 3.0e-15,
              "scene grid clipping: noncoincident patches partition triangle area");
    check(first.patchesForCell(grid().cellIndex(1, 1, 1)).size() == 1
              && first.verticesForPatch(first.patches[0]).size() == 5,
          "scene grid clipping: cell and flattened-vertex views are bounded");
    validateSceneFluidGridPatches(
        first,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections);
}

void testContactDimensionAndBoundaryAmbiguity() {
    Pipeline touching(triangleScene(2.9));
    const auto touchingPatches = clipSceneFluidSurfaceToCells(
        touching.surface.definition,
        touching.state,
        grid(),
        touching.candidates,
        touching.intersections);
    const auto contact = touchingPatches.patchesForCell(
        grid().cellIndex(2, 2, 1));
    check(contact.size() == 1
              && contact.front().dimension == SceneFluidPatchDimension::Point
              && contact.front().vertexCount == 1
              && contact.front().areaSquareMeters == 0.0
              && (contact.front().coincidentBoundaryPlanes
                  & (CellBoundaryXMinus | CellBoundaryYMinus))
                  == (CellBoundaryXMinus | CellBoundaryYMinus),
          "scene grid clipping: exact point contact remains explicit and zero-area");

    Pipeline boundary(triangleScene(2.7, 2.0));
    const auto boundaryPatches = clipSceneFluidSurfaceToCells(
        boundary.surface.definition,
        boundary.state,
        grid(),
        boundary.candidates,
        boundary.intersections);
    double duplicatedArea = 0.0;
    std::size_t lowerOwners = 0;
    std::size_t upperOwners = 0;
    for (const auto& patch : boundaryPatches.patches) {
        duplicatedArea += patch.areaSquareMeters;
        lowerOwners += (patch.coincidentBoundaryPlanes
                        & CellBoundaryZMinus) != 0;
        upperOwners += (patch.coincidentBoundaryPlanes
                        & CellBoundaryZPlus) != 0;
    }
    checkNear(duplicatedArea, 2.56, 6.0e-15,
              "scene grid clipping: shared-plane ambiguity remains visibly duplicated");
    check(boundaryPatches.patches.size() == 6
              && lowerOwners == 3 && upperOwners == 3,
          "scene grid clipping: both adjacent boundary owners are explicitly flagged");
}

void testToleranceAndTransactionalLimits() {
    Pipeline pipeline(triangleScene());
    SceneFluidGridCandidateSettings broadSettings;
    broadSettings.boundingPaddingMeters = 0.2;
    const auto broad = buildSceneFluidGridCandidates(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        broadSettings);
    SceneFluidGridIntersectionSettings narrowSettings;
    narrowSettings.separationToleranceMeters = 0.15;
    const auto tolerant = intersectSceneFluidSurfaceWithGrid(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        broad,
        narrowSettings);
    expectInvalid(
        [&] { static_cast<void>(clipSceneFluidSurfaceToCells(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            broad,
            tolerant)); },
        "scene grid clipping: proximity-only intersections cannot become physical patches");

    SceneFluidGridPatchLimits limits;
    limits.maximumPatches = 2;
    expectLimited(
        [&] { static_cast<void>(clipSceneFluidSurfaceToCells(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            limits)); },
        "scene grid clipping: patch count is bounded before publication");
    limits = {};
    limits.maximumVertices = 10;
    expectLimited(
        [&] { static_cast<void>(clipSceneFluidSurfaceToCells(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            limits)); },
        "scene grid clipping: flattened vertex count is bounded");
    limits = {};
    limits.maximumPatchBytes =
        3 * sizeof(SceneFluidCellPatch)
        + 10 * sizeof(SceneFluidClippedVertex);
    expectLimited(
        [&] { static_cast<void>(clipSceneFluidSurfaceToCells(
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections,
            limits)); },
        "scene grid clipping: patch storage is byte-bounded");

    const auto accepted = clipSceneFluidSurfaceToCells(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections);
    auto corrupt = accepted;
    corrupt.vertices.front().positionMeters.x += 0.01;
    expectInvalid(
        [&] { validateSceneFluidGridPatches(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            grid(),
            pipeline.candidates,
            pipeline.intersections); },
        "scene grid clipping: vertex corruption is rejected transactionally");
    validateSceneFluidGridPatches(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections);
}

} // namespace

int main() {
    testReusableTriangleBoxPrimitive();
    testExactClippedPatches();
    testContactDimensionAndBoundaryAmbiguity();
    testToleranceAndTransactionalLimits();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-clipping test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-clipping tests passed\n");
    return 0;
}
