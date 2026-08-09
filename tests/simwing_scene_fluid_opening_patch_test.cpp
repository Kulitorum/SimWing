#include "scene_fluid_opening_patch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace {

using namespace simwing::fsi;

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
                     "FAIL: %s (actual %.17g expected %.17g)\n",
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

template<std::size_t VertexCount>
std::array<Vec2, 3> intrinsicChart(
    const std::array<Vec3, VertexCount>& positions,
    const std::array<std::size_t, 3>& vertices) {
    const Vec3& first = positions[vertices[0]];
    const Vec3& second = positions[vertices[1]];
    const Vec3& third = positions[vertices[2]];
    const Vec3 edge{second.x - first.x,
                    second.y - first.y,
                    second.z - first.z};
    const Vec3 diagonal{third.x - first.x,
                        third.y - first.y,
                        third.z - first.z};
    const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
    const double projected = (diagonal.x * edge.x
                              + diagonal.y * edge.y
                              + diagonal.z * edge.z) / edgeLength;
    const double diagonalSquared = diagonal.x * diagonal.x
        + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
    return {{{0.0, 0.0},
             {edgeLength, 0.0},
             {projected, std::sqrt(std::max(
                 0.0, diagonalSquared - projected * projected))}}};
}

Scene openSquarePyramidScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-opening-patch";
    scene.metadata.exporterVersion = "scene-fluid-opening-patch-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 5> positions{{
        {1.0, 1.0, 1.0},
        {2.0, 0.5, 0.5},
        {2.0, 1.5, 0.5},
        {2.0, 1.5, 1.5},
        {2.0, 0.5, 1.5},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        scene.triangles.push_back({
            500 + face,
            {10 + faces[face][0], 10 + faces[face][1],
             10 + faces[face][2]},
            intrinsicChart(positions, faces[face]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {0.0, 0.0, 0.0}, {4.0, 4.0, 4.0}};
}

struct Fixture {
    Scene scene = openSquarePyramidScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structureAssembly = assembleSceneStructure(scene);
    Structure structure{structureAssembly.definition};

    SceneFluidSurfaceState state() const {
        return captureSceneFluidSurfaceState(
            surface.definition, structureAssembly.mappings, structure);
    }
};

struct OpeningEpoch {
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;
};

OpeningEpoch buildEpoch(const Fixture& fixture,
                        const SceneFluidSurfaceState& state,
                        const fluid::PeriodicCartesianGrid& targetGrid = grid()) {
    OpeningEpoch result;
    result.caps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    result.quadrature = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, state, result.caps);
    result.patches = buildSceneFluidOpeningGridPatches(
        fixture.surface.definition, state, result.caps,
        result.quadrature, targetGrid);
    return result;
}

void testFaceAndCellOwnershipAcrossAcceptedMotion() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "opening-patch square-pyramid fixture assembles");
    const auto initialState = fixture.state();
    const auto initial = buildEpoch(fixture, initialState);
    const auto repeated = buildEpoch(fixture, initialState);
    check(initial.patches == repeated.patches
              && initial.patches.version
                  == sceneFluidOpeningGridPatchVersion
              && initial.patches.fingerprint != 0
              && initial.patches.pointRanges.size() == 2
              && initial.patches.patches.size() == 6,
          "grid-aligned square mouth has deterministic unique patch ownership");
    check(initial.patches.candidateCellCount == 16,
          "grid-aligned fan triangles retain bounded adjacent candidates");
    check(std::ranges::all_of(
              initial.patches.patches,
              [](const auto& patch) {
                  return patch.ownerKind
                          == SceneFluidOpeningPatchOwnerKind::Face
                      && patch.faceAxis
                          == SceneFluidOpeningPatchFaceAxis::X
                      && patch.faceI == 2
                      && patch.cell.i == 1
                      && patch.negativeSideRegionId == 2
                      && patch.positiveSideRegionId == 1;
              }),
          "coincident mouth tiles are owned once by canonical X faces");
    checkNear(initial.patches.totalAreaSquareMeters, 1.0, 1.0e-14,
              "grid face partitions preserve exact mouth area");
    checkNear(
        initial.patches.totalSurfaceSweepRateCubicMetersPerSecond,
        0.0, 0.0,
        "stationary grid face partitions have zero cap sweep");
    for (std::size_t point = 0;
         point < initial.patches.pointRanges.size(); ++point) {
        const auto& range = initial.patches.pointRanges[point];
        check(range.sourcePointStableId
                  == initial.quadrature.points[point].stableId
                  && range.patchCount == 3
                  && initial.patches.patchesForPoint(range).size() == 3,
              "each source fan triangle owns three grid-face patches");
        checkNear(range.areaSquareMeters, 0.5, 1.0e-14,
                  "each source fan triangle closes its area partition");
    }
    for (const auto& patch : initial.patches.patches) {
        check(patch.stableId != 0
                  && initial.patches.verticesForPatch(patch).size()
                      == patch.vertexCount,
              "owned opening patch exposes bounded exact polygon geometry");
    }
    validateSceneFluidOpeningGridPatches(
        initial.patches, fixture.surface.definition, initialState,
        initial.caps, initial.quadrature, grid());

    for (std::size_t node = 0;
         node < fixture.structure.definition().nodes.size(); ++node) {
        fixture.structure.addExternalForce(
            node,
            {fixture.structure.definition().nodes[node].massKg,
             0.0, 0.0});
    }
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.1;
    settings.substeps = 1;
    settings.constraintIterations = 4;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(settings);
    const auto movedState = fixture.state();
    const auto moved = buildEpoch(fixture, movedState);
    check(diagnostics.finite
              && moved.patches.fingerprint != initial.patches.fingerprint
              && moved.patches.candidateCellCount == 8
              && moved.patches.patches.size() == 6,
          "accepted off-face motion creates one new bounded patch epoch");
    check(std::ranges::all_of(
              moved.patches.patches,
              [](const auto& patch) {
                  return patch.ownerKind
                          == SceneFluidOpeningPatchOwnerKind::Cell
                      && patch.cell.i == 2;
              }),
          "off-face mouth partitions are uniquely cell-owned");
    checkNear(moved.patches.totalAreaSquareMeters, 1.0, 2.0e-14,
              "off-face patch partition preserves mouth area");
    checkNear(moved.patches.totalSurfaceSweepRateCubicMetersPerSecond,
              0.1, 5.0e-13,
              "off-face patch partition preserves analytic cap sweep");
}

void testBoundaryCorruptionAndLimits() {
    Fixture fixture;
    const auto state = fixture.state();
    const auto accepted = buildEpoch(fixture, state);

    auto corrupt = accepted.patches;
    corrupt.patches.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidOpeningGridPatches(
            corrupt, fixture.surface.definition, state,
            accepted.caps, accepted.quadrature, grid()); },
        "opening-patch validation rejects payload corruption");

    const fluid::PeriodicCartesianGrid foreignGrid{
        {8, 4, 4}, {0.0, 0.0, 0.0}, {4.0, 4.0, 4.0}};
    expectInvalid(
        [&] { validateSceneFluidOpeningGridPatches(
            accepted.patches, fixture.surface.definition, state,
            accepted.caps, accepted.quadrature, foreignGrid); },
        "opening patches reject a foreign grid");

    const fluid::PeriodicCartesianGrid boundaryGrid{
        {2, 4, 4}, {0.0, 0.0, 0.0}, {2.0, 4.0, 4.0}};
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, boundaryGrid)); },
        "opening patches reject unresolved periodic-boundary ownership");

    SceneFluidOpeningGridPatchSettings invalidSettings;
    invalidSettings.relativeAreaTolerance = -1.0;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, grid(), invalidSettings)); },
        "opening patches reject invalid tolerances");

    SceneFluidOpeningGridPatchLimits limits;
    limits.maximumCandidateCells = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, grid(), {}, limits)); },
        "opening patches bound candidate work");
    limits = {};
    limits.maximumPatches = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, grid(), {}, limits)); },
        "opening patches bound patch count");
    limits = {};
    limits.maximumVertices = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, grid(), {}, limits)); },
        "opening patches bound polygon vertices");
    limits = {};
    limits.maximumPatchBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningGridPatches(
            fixture.surface.definition, state, accepted.caps,
            accepted.quadrature, grid(), {}, limits)); },
        "opening patches bound aggregate storage");
}

} // namespace

int main() {
    try {
        testFaceAndCellOwnershipAcrossAcceptedMotion();
        testBoundaryCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d opening-patch check(s) failed\n", failures);
        return 1;
    }
    std::puts("all scene fluid opening-patch checks passed");
    return 0;
}
