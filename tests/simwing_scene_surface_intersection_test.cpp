#include "fluid/scene_surface_intersection.h"

#include <algorithm>
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

Scene triangleScene(const double farCoordinate = 2.7,
                    const double zMeters = 1.5) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-intersection";
    scene.metadata.exporterVersion = "scene-surface-intersection-test/1";
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

struct Fixture {
    Scene scene = triangleScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structure = assembleSceneStructure(scene);
};

void testExactNarrowPhase() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structure.ok(),
          "scene grid intersection: authoritative triangle assembles");
    Structure structure(fixture.structure.definition);
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto candidates = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid());
    const auto first = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition, state, grid(), candidates);
    const auto second = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition, state, grid(), candidates);

    check(candidates.candidates.size() == 4
              && first == second
              && first.fingerprint != 0
              && first.surfaceStateFingerprint == state.fingerprint
              && first.candidateSetFingerprint == candidates.fingerprint,
          "scene grid intersection: exact results are deterministic and epoch-bound");
    check(first.intersections.size() == 3
              && first.triangleIntersectionCounts
                  == std::vector<std::size_t>{3}
              && first.intersections[0].cellIndex
                  == grid().cellIndex(1, 1, 1)
              && first.intersections[1].cellIndex
                  == grid().cellIndex(2, 1, 1)
              && first.intersections[2].cellIndex
                  == grid().cellIndex(1, 2, 1),
          "scene grid intersection: separating axes prune the AABB false positive");
    check(first.intersectionsForCell(grid().cellIndex(1, 1, 1)).size() == 1
              && first.intersectionsForCell(
                  grid().cellIndex(2, 2, 1)).empty(),
          "scene grid intersection: cell-major exact lookup is explicit");
    validateSceneFluidGridIntersections(
        first, fixture.surface.definition, state, grid(), candidates);

    Scene reversedScene = triangleScene();
    std::swap(reversedScene.triangles.front().vertexIds[1],
              reversedScene.triangles.front().vertexIds[2]);
    const auto reversedSurface = assembleSceneFluidSurface(reversedScene);
    const auto reversedStructure = assembleSceneStructure(reversedScene);
    check(reversedSurface.ok() && reversedStructure.ok(),
          "scene grid intersection: reversed-winding fixture assembles");
    if (!reversedSurface.ok() || !reversedStructure.ok()) {
        return;
    }
    Structure reversedOwner(reversedStructure.definition);
    const auto reversedState = captureSceneFluidSurfaceState(
        reversedSurface.definition,
        reversedStructure.mappings,
        reversedOwner);
    const auto reversedCandidates = buildSceneFluidGridCandidates(
        reversedSurface.definition, reversedState, grid());
    const auto reversed = intersectSceneFluidSurfaceWithGrid(
        reversedSurface.definition,
        reversedState,
        grid(),
        reversedCandidates);
    check(reversed.intersections == first.intersections,
          "scene grid intersection: triangle winding cannot change geometric contact");
}

void testTouchingAndGridPlaneCoverage() {
    Scene touchingScene = triangleScene(2.9);
    const auto touchingSurface = assembleSceneFluidSurface(touchingScene);
    const auto touchingStructure = assembleSceneStructure(touchingScene);
    Structure touchingOwner(touchingStructure.definition);
    const auto touchingState = captureSceneFluidSurfaceState(
        touchingSurface.definition,
        touchingStructure.mappings,
        touchingOwner);
    const auto touchingCandidates = buildSceneFluidGridCandidates(
        touchingSurface.definition, touchingState, grid());
    const auto touching = intersectSceneFluidSurfaceWithGrid(
        touchingSurface.definition,
        touchingState,
        grid(),
        touchingCandidates);
    check(touching.intersections.size() == 4
              && !touching.intersectionsForCell(
                  grid().cellIndex(2, 2, 1)).empty(),
          "scene grid intersection: exact point contact is retained conservatively");

    Scene boundaryScene = triangleScene(2.7, 2.0);
    const auto boundarySurface = assembleSceneFluidSurface(boundaryScene);
    const auto boundaryStructure = assembleSceneStructure(boundaryScene);
    Structure boundaryOwner(boundaryStructure.definition);
    const auto boundaryState = captureSceneFluidSurfaceState(
        boundarySurface.definition,
        boundaryStructure.mappings,
        boundaryOwner);
    const auto boundaryCandidates = buildSceneFluidGridCandidates(
        boundarySurface.definition, boundaryState, grid());
    const auto boundary = intersectSceneFluidSurfaceWithGrid(
        boundarySurface.definition,
        boundaryState,
        grid(),
        boundaryCandidates);
    check(boundaryCandidates.candidates.size() == 8
              && boundary.intersections.size() == 6
              && boundary.triangleIntersectionCounts[0] == 6,
          "scene grid intersection: an internal grid-plane surface keeps both exact sides");
}

void testToleranceAndAcceptedEpochBinding() {
    Fixture fixture;
    Structure structure(fixture.structure.definition);
    const auto initialState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto initialCandidates = buildSceneFluidGridCandidates(
        fixture.surface.definition, initialState, grid());
    const auto initialExact = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition,
        initialState,
        grid(),
        initialCandidates);
    SceneFluidGridCandidateSettings broadSettings;
    broadSettings.boundingPaddingMeters = 0.2;
    const auto broad = buildSceneFluidGridCandidates(
        fixture.surface.definition,
        initialState,
        grid(),
        broadSettings);
    SceneFluidGridIntersectionSettings narrowSettings;
    narrowSettings.separationToleranceMeters = 0.15;
    const auto tolerant = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition,
        initialState,
        grid(),
        broad,
        narrowSettings);
    check(!tolerant.intersectionsForCell(
              grid().cellIndex(2, 2, 1)).empty(),
          "scene grid intersection: explicit supported tolerance retains a near contact");
    narrowSettings.separationToleranceMeters = 0.21;
    expectInvalid(
        [&] { static_cast<void>(intersectSceneFluidSurfaceWithGrid(
            fixture.surface.definition,
            initialState,
            grid(),
            broad,
            narrowSettings)); },
        "scene grid intersection: tolerance cannot exceed broad-phase padding");

    for (std::size_t node = 0;
         node < fixture.structure.definition.nodes.size(); ++node) {
        structure.addExternalForce(
            node,
            {fixture.structure.definition.nodes[node].massKg, 0.0, 0.0});
    }
    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = 0.01;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    static_cast<void>(structure.step(stepSettings));
    const auto movedState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto movedCandidates = buildSceneFluidGridCandidates(
        fixture.surface.definition, movedState, grid());
    const auto moved = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition,
        movedState,
        grid(),
        movedCandidates);
    check(moved.intersections == initialExact.intersections
              && moved.surfaceStateFingerprint
                  != initialExact.surfaceStateFingerprint
              && moved.fingerprint != initialExact.fingerprint,
          "scene grid intersection: accepted Structure motion creates a new exact epoch");
    expectInvalid(
        [&] { validateSceneFluidGridIntersections(
            moved,
            fixture.surface.definition,
            initialState,
            grid(),
            broad); },
        "scene grid intersection: exact results cannot cross surface epochs");
}

void testTransactionalValidationAndLimits() {
    Fixture fixture;
    Structure structure(fixture.structure.definition);
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto candidates = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid());
    const auto accepted = intersectSceneFluidSurfaceWithGrid(
        fixture.surface.definition, state, grid(), candidates);

    SceneFluidGridIntersectionLimits limits;
    limits.maximumIntersections = 2;
    expectLimited(
        [&] { static_cast<void>(intersectSceneFluidSurfaceWithGrid(
            fixture.surface.definition,
            state,
            grid(),
            candidates,
            {},
            limits)); },
        "scene grid intersection: exact-pair count is bounded before publication");
    limits = {};
    limits.maximumIntersectionBytes =
        sizeof(std::size_t)
        + 2 * sizeof(SceneFluidCellIntersection);
    expectLimited(
        [&] { static_cast<void>(intersectSceneFluidSurfaceWithGrid(
            fixture.surface.definition,
            state,
            grid(),
            candidates,
            {},
            limits)); },
        "scene grid intersection: exact-pair storage is byte-bounded");

    auto corrupt = accepted;
    corrupt.intersections.pop_back();
    expectInvalid(
        [&] { validateSceneFluidGridIntersections(
            corrupt,
            fixture.surface.definition,
            state,
            grid(),
            candidates); },
        "scene grid intersection: missing exact pairs are rejected transactionally");
    validateSceneFluidGridIntersections(
        accepted, fixture.surface.definition, state, grid(), candidates);
}

} // namespace

int main() {
    testExactNarrowPhase();
    testTouchingAndGridPlaneCoverage();
    testToleranceAndAcceptedEpochBinding();
    testTransactionalValidationAndLimits();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-intersection test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-intersection tests passed\n");
    return 0;
}
