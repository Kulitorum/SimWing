#include "scene_fluid_opening_quadrature.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>

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

double materialSurfaceSweepRate(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state) {
    double result = 0.0;
    for (const auto& triangle : surface.triangles) {
        if (triangle.negativeSideRegionIndex
            == triangle.positiveSideRegionIndex) {
            continue;
        }
        const auto& first = state.vertices[triangle.vertexIndices[0]];
        const auto& second = state.vertices[triangle.vertexIndices[1]];
        const auto& third = state.vertices[triangle.vertexIndices[2]];
        const Vec3 firstEdge{
            second.positionMeters.x - first.positionMeters.x,
            second.positionMeters.y - first.positionMeters.y,
            second.positionMeters.z - first.positionMeters.z,
        };
        const Vec3 secondEdge{
            third.positionMeters.x - first.positionMeters.x,
            third.positionMeters.y - first.positionMeters.y,
            third.positionMeters.z - first.positionMeters.z,
        };
        const Vec3 areaVector{
            0.5 * (firstEdge.y * secondEdge.z
                   - firstEdge.z * secondEdge.y),
            0.5 * (firstEdge.z * secondEdge.x
                   - firstEdge.x * secondEdge.z),
            0.5 * (firstEdge.x * secondEdge.y
                   - firstEdge.y * secondEdge.x),
        };
        const Vec3 centroidVelocity{
            (first.velocityMetersPerSecond.x
             + second.velocityMetersPerSecond.x
             + third.velocityMetersPerSecond.x) / 3.0,
            (first.velocityMetersPerSecond.y
             + second.velocityMetersPerSecond.y
             + third.velocityMetersPerSecond.y) / 3.0,
            (first.velocityMetersPerSecond.z
             + second.velocityMetersPerSecond.z
             + third.velocityMetersPerSecond.z) / 3.0,
        };
        result += areaVector.x * centroidVelocity.x
            + areaVector.y * centroidVelocity.y
            + areaVector.z * centroidVelocity.z;
    }
    return result;
}

Scene openSquarePyramidScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-opening-quadrature";
    scene.metadata.exporterVersion =
        "scene-fluid-opening-quadrature-test/1";
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

void testAcceptedOpeningKinematics() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "opening quadrature square-pyramid fixture assembles");
    const auto initialState = fixture.state();
    const auto initialCaps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, initialState);
    const auto initial = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, initialState, initialCaps);
    const auto repeated = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, initialState, initialCaps);
    check(initial == repeated
              && initial.version == sceneFluidOpeningQuadratureVersion
              && initial.fingerprint != 0
              && initial.openings.size() == 1
              && initial.points.size() == 2,
          "opening quadrature is immutable and deterministic");
    const auto& opening = initial.openings.front();
    check(opening.openingId == 700
              && opening.negativeSideRegionId == 2
              && opening.positiveSideRegionId == 1
              && opening.role == OpeningRole::Intake
              && opening.firstPoint == 0
              && opening.pointCount == 2
              && opening.unitNormalNegativeToPositive.x > 0.99,
          "opening quadrature preserves stable identity and region direction");
    checkNear(opening.areaSquareMeters, 1.0, 1.0e-14,
              "opening quadrature preserves exact square mouth area");
    checkNear(initial.totalAreaSquareMeters, 1.0, 1.0e-14,
              "opening quadrature total area closes");
    checkNear(opening.surfaceSweepRateCubicMetersPerSecond, 0.0, 0.0,
              "stationary opening has zero surface-sweep rate");
    check(initial.points[0].stableId != 0
              && initial.points[1].stableId != 0
              && initial.points[0].stableId != initial.points[1].stableId
              && initial.points[0].triangleOrdinal == 0
              && initial.points[1].triangleOrdinal == 1,
          "fan triangles receive unique stable quadrature identities");
    checkNear(initial.points[0].areaSquareMeters, 0.5, 1.0e-14,
              "first fan triangle owns half the square mouth");
    checkNear(initial.points[1].areaSquareMeters, 0.5, 1.0e-14,
              "second fan triangle owns half the square mouth");
    validateSceneFluidOpeningQuadrature(
        initial, fixture.surface.definition, initialState, initialCaps);

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
    const auto movedCaps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, movedState);
    const auto moved = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, movedState, movedCaps);
    check(diagnostics.finite
              && moved.fingerprint != initial.fingerprint
              && moved.openingCapFingerprint != initial.openingCapFingerprint,
          "accepted Structure motion creates a distinct opening quadrature");
    check(moved.points[0].stableId == initial.points[0].stableId
              && moved.points[1].stableId == initial.points[1].stableId,
          "accepted motion preserves topology-stable point identities");
    checkNear(moved.openings.front().centroidMeters.x,
              opening.centroidMeters.x + 0.01, 2.0e-13,
              "opening quadrature follows accepted cap position");
    checkNear(moved.points[0].velocityMetersPerSecond.x,
              0.1, 2.0e-13,
              "opening point samples accepted linear vertex velocity");
    checkNear(moved.points[0].surfaceSweepRateCubicMetersPerSecond,
              0.05, 2.0e-13,
              "first cap triangle integrates analytic normal sweep rate");
    checkNear(moved.openings.front()
                  .surfaceSweepRateCubicMetersPerSecond,
              0.1, 4.0e-13,
              "opening aggregate integrates analytic normal sweep rate");
    checkNear(moved.totalSurfaceSweepRateCubicMetersPerSecond,
              0.1, 4.0e-13,
              "global opening sweep ledger closes");
    checkNear(materialSurfaceSweepRate(
                  fixture.surface.definition, movedState)
                  + moved.totalSurfaceSweepRateCubicMetersPerSecond,
              0.0, 5.0e-14,
              "material and virtual-cap sweeps close under rigid translation");
}

void testBindingCorruptionAndLimits() {
    Fixture fixture;
    const auto state = fixture.state();
    const auto caps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    const auto accepted = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, state, caps);

    auto corrupt = accepted;
    corrupt.points.front().velocityMetersPerSecond.x += 0.01;
    expectInvalid(
        [&] { validateSceneFluidOpeningQuadrature(
            corrupt, fixture.surface.definition, state, caps); },
        "opening quadrature rejects kinematic payload corruption");

    for (std::size_t node = 0;
         node < fixture.structure.definition().nodes.size(); ++node) {
        fixture.structure.addExternalForce(
            node,
            {0.0, fixture.structure.definition().nodes[node].massKg, 0.0});
    }
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.1;
    settings.substeps = 1;
    settings.constraintIterations = 4;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    static_cast<void>(fixture.structure.step(settings));
    const auto foreignState = fixture.state();
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningQuadrature(
            fixture.surface.definition, foreignState, caps)); },
        "opening quadrature rejects a cap from another accepted state");

    SceneFluidOpeningQuadratureLimits limits;
    limits.maximumOpenings = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningQuadrature(
            fixture.surface.definition, state, caps, limits)); },
        "opening quadrature bounds opening count");
    limits = {};
    limits.maximumPoints = 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningQuadrature(
            fixture.surface.definition, state, caps, limits)); },
        "opening quadrature bounds point count");
    limits = {};
    limits.maximumQuadratureBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningQuadrature(
            fixture.surface.definition, state, caps, limits)); },
        "opening quadrature bounds aggregate storage");
}

} // namespace

int main() {
    try {
        testAcceptedOpeningKinematics();
        testBindingCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d opening-quadrature check(s) failed\n", failures);
        return 1;
    }
    std::puts("all scene fluid opening-quadrature checks passed");
    return 0;
}
