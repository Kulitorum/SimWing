#include "fluid/scene_surface_geometry.h"

#include <cmath>
#include <cstdio>
#include <limits>
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

Scene surfaceScene(const double zMeters = 1.5) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-geometry";
    scene.metadata.exporterVersion = "scene-surface-geometry-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.2, zMeters}},
        {11, {2.8, 1.2, zMeters}},
        {12, {2.8, 2.8, zMeters}},
        {13, {1.2, 2.8, zMeters}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12}, {{{0.0, 0.0}, {1.6, 0.0}, {1.6, 1.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 12, 13}, {{{0.0, 0.0}, {1.6, 1.6}, {0.0, 1.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    scene.openings = {
        {600, {10, 11, 12, 13}, 1, 2, OpeningRole::Intake},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return PeriodicCartesianGrid({4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

struct Fixture {
    Scene scene = surfaceScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structure = assembleSceneStructure(scene);
};

void testCanonicalBroadPhase() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structure.ok(),
          "scene grid geometry: authoritative fixture assembles");
    Structure structure(fixture.structure.definition);
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto first = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid());
    const auto second = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid());

    check(first == second && first.fingerprint != 0
              && first.version == sceneFluidGridCandidateVersion
              && first.surfaceDefinitionFingerprint
                  == fixture.surface.definition.fingerprint
              && first.surfaceStateFingerprint == state.fingerprint,
          "scene grid geometry: repeated broad phases are bit-deterministic and bound");
    check(first.triangleBounds.size() == 2
              && first.candidates.size() == 8
              && first.triangleBounds[0]
                  == SceneFluidTriangleCandidateBounds{
                      0, 500, {1, 1, 1}, {2, 2, 1}, 4}
              && first.triangleBounds[1]
                  == SceneFluidTriangleCandidateBounds{
                      1, 501, {1, 1, 1}, {2, 2, 1}, 4},
          "scene grid geometry: padded AABB ranges retain exact triangle ownership");
    const auto cellCandidates = first.candidatesForCell(
        grid().cellIndex(1, 1, 1));
    check(cellCandidates.size() == 2
              && cellCandidates[0].triangleId == 500
              && cellCandidates[1].triangleId == 501
              && first.candidatesForCell(grid().cellIndex(0, 0, 0)).empty(),
          "scene grid geometry: cell-major lookup is canonical and bounded");
    validateSceneFluidGridCandidates(
        first, fixture.surface.definition, state, grid());
}

void testBoundaryAndPaddingCoverage() {
    Fixture fixture;
    Structure structure(fixture.structure.definition);
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    SceneFluidGridCandidateSettings padded;
    padded.boundingPaddingMeters = 0.25;
    const auto expanded = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid(), padded);
    check(expanded.candidates.size() == 32
              && expanded.triangleBounds[0].firstCell
                  == GridCellCoordinate{0, 0, 1}
              && expanded.triangleBounds[0].lastCell
                  == GridCellCoordinate{3, 3, 1},
          "scene grid geometry: explicit padding conservatively expands candidates");

    Scene boundaryScene = surfaceScene(2.0);
    const auto boundarySurface = assembleSceneFluidSurface(boundaryScene);
    const auto boundaryStructure = assembleSceneStructure(boundaryScene);
    Structure boundaryOwner(boundaryStructure.definition);
    const auto boundaryState = captureSceneFluidSurfaceState(
        boundarySurface.definition,
        boundaryStructure.mappings,
        boundaryOwner);
    const auto boundary = buildSceneFluidGridCandidates(
        boundarySurface.definition, boundaryState, grid());
    check(boundary.candidates.size() == 16
              && boundary.triangleBounds[0].firstCell.k == 1
              && boundary.triangleBounds[0].lastCell.k == 2,
          "scene grid geometry: a zero-thickness surface on a grid plane reaches both cells");
}

void testAcceptedEpochIdentity() {
    Fixture fixture;
    Structure structure(fixture.structure.definition);
    const auto initialState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto initial = buildSceneFluidGridCandidates(
        fixture.surface.definition, initialState, grid());

    for (std::size_t node = 0;
         node < fixture.structure.definition.nodes.size(); ++node) {
        structure.addExternalForce(
            node,
            {fixture.structure.definition.nodes[node].massKg, 0.0, 0.0});
    }
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.01;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    const auto diagnostics = structure.step(settings);
    const auto movedState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto moved = buildSceneFluidGridCandidates(
        fixture.surface.definition, movedState, grid());
    check(diagnostics.finite
              && movedState.acceptedStepCount == 1
              && movedState.fingerprint != initialState.fingerprint
              && moved.surfaceStateFingerprint == movedState.fingerprint
              && moved.fingerprint != initial.fingerprint
              && moved.candidates == initial.candidates,
          "scene grid geometry: accepted motion changes epoch identity even within one cell set");
    expectInvalid(
        [&] { validateSceneFluidGridCandidates(
            initial, fixture.surface.definition, movedState, grid()); },
        "scene grid geometry: candidates cannot cross accepted Structure epochs");
}

void testTransactionalValidationAndLimits() {
    Fixture fixture;
    Structure structure(fixture.structure.definition);
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structure.mappings, structure);
    const auto accepted = buildSceneFluidGridCandidates(
        fixture.surface.definition, state, grid());

    auto corruptState = state;
    corruptState.vertices.front().positionMeters.x += 0.01;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidGridCandidates(
            fixture.surface.definition, corruptState, grid())); },
        "scene grid geometry: corrupted state fingerprint is rejected");

    const PeriodicCartesianGrid tooSmall(
        {4, 4, 4}, {}, {2.5, 2.5, 2.5});
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidGridCandidates(
            fixture.surface.definition, state, tooSmall)); },
        "scene grid geometry: out-of-domain geometry requires an explicit topology decision");

    SceneFluidGridCandidateSettings invalidSettings;
    invalidSettings.boundingPaddingMeters = -1.0;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidGridCandidates(
            fixture.surface.definition, state, grid(), invalidSettings)); },
        "scene grid geometry: negative padding is rejected");

    SceneFluidGridCandidateLimits limits;
    limits.maximumCandidates = 7;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidGridCandidates(
            fixture.surface.definition, state, grid(), {}, limits)); },
        "scene grid geometry: candidate count is bounded before publication");
    limits = {};
    limits.maximumCandidateBytes =
        2 * sizeof(SceneFluidTriangleCandidateBounds)
        + 7 * sizeof(SceneFluidCellCandidate);
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidGridCandidates(
            fixture.surface.definition, state, grid(), {}, limits)); },
        "scene grid geometry: candidate storage is byte-bounded");

    auto corruptCandidates = accepted;
    ++corruptCandidates.candidates.front().triangleId;
    expectInvalid(
        [&] { validateSceneFluidGridCandidates(
            corruptCandidates, fixture.surface.definition, state, grid()); },
        "scene grid geometry: candidate corruption is rejected transactionally");
    validateSceneFluidGridCandidates(
        accepted, fixture.surface.definition, state, grid());
}

} // namespace

int main() {
    testCanonicalBroadPhase();
    testBoundaryAndPaddingCoverage();
    testAcceptedEpochIdentity();
    testTransactionalValidationAndLimits();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-geometry test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-geometry tests passed\n");
    return 0;
}
