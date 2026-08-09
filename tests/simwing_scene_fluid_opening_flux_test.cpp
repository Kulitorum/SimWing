#include "scene_fluid_opening_flux.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
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
    scene.metadata.designChecksum = "sha256:scene-fluid-opening-flux";
    scene.metadata.exporterVersion = "scene-fluid-opening-flux-test/1";
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

struct OpeningGeometry {
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;
};

OpeningGeometry buildGeometry(const Fixture& fixture,
                              const SceneFluidSurfaceState& state) {
    OpeningGeometry result;
    result.caps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    result.quadrature = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, state, result.caps);
    result.patches = buildSceneFluidOpeningGridPatches(
        fixture.surface.definition, state, result.caps,
        result.quadrature, grid());
    return result;
}

SceneFluidOpeningFluxSet evaluate(
    const Fixture& fixture,
    const SceneFluidSurfaceState& state,
    const OpeningGeometry& geometry,
    const fluid::MacVelocityField& velocity,
    const SceneFluidOpeningFluxLimits& limits = {}) {
    return evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, state, geometry.caps,
        geometry.quadrature, geometry.patches, grid(), velocity, limits);
}

fluid::MacVelocityField tiledFaceVelocity(const double sign = 1.0) {
    fluid::MacVelocityField velocity(grid());
    const std::array<std::array<double, 2>, 2> values{{
        {{1.0, 3.0}}, {{2.0, 4.0}},
    }};
    for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t j = 0; j < 2; ++j) {
            velocity.xFaces()[grid().cellIndex(2, j, k)] =
                sign * values[j][k];
        }
    }
    return velocity;
}

fluid::MacVelocityField linearXVelocity() {
    fluid::MacVelocityField velocity(grid());
    const auto counts = grid().cellCounts();
    const auto lower = grid().lowerMeters();
    const auto spacing = grid().cellSpacingMeters();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                velocity.xFaces()[grid().cellIndex(i, j, k)] =
                    lower.x + static_cast<double>(i) * spacing.x;
            }
        }
    }
    return velocity;
}

void testResolvedFaceFluxAndOrientation() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "opening-flux square-pyramid fixture assembles");
    const auto state = fixture.state();
    const auto geometry = buildGeometry(fixture, state);
    const auto velocity = tiledFaceVelocity();
    const auto first = evaluate(fixture, state, geometry, velocity);
    const auto repeated = evaluate(fixture, state, geometry, velocity);
    check(first == repeated
              && first.version == sceneFluidOpeningFluxVersion
              && first.fingerprint != 0
              && first.velocityFingerprint != 0
              && first.openings.size() == 1
              && first.samples.size() == 6
              && first.velocityEvaluationCount == 6,
          "resolved opening flux is immutable, deterministic, and face sampled");
    const auto& opening = first.openings.front();
    check(opening.openingId == 700
              && opening.negativeSideRegionId == 2
              && opening.positiveSideRegionId == 1
              && opening.sampleCount == 6,
          "opening flux preserves authored negative-to-positive identity");
    checkNear(opening.areaSquareMeters, 1.0, 1.0e-14,
              "resolved opening flux preserves mouth area");
    checkNear(opening.fluidVolumeFlowRateCubicMetersPerSecond,
              2.5, 2.0e-14,
              "partial MAC tiles integrate their exact resolved face flow");
    checkNear(opening.surfaceSweepRateCubicMetersPerSecond,
              0.0, 0.0,
              "stationary mouth has no geometric sweep");
    checkNear(opening.relativeVolumeFlowRateCubicMetersPerSecond,
              2.5, 2.0e-14,
              "stationary mouth relative flow equals fluid flow");
    validateSceneFluidOpeningFlux(
        first, fixture.surface.definition, state, geometry.caps,
        geometry.quadrature, geometry.patches, grid(), velocity);

    const auto reversedVelocity = tiledFaceVelocity(-1.0);
    const auto reversed = evaluate(
        fixture, state, geometry, reversedVelocity);
    checkNear(reversed.totalRelativeVolumeFlowRateCubicMetersPerSecond,
              -2.5, 2.0e-14,
              "negative MAC velocity reverses authored opening-flow sign");
}

void testMovingOffFaceRelativeFlux() {
    Fixture fixture;
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
    const auto state = fixture.state();
    const auto geometry = buildGeometry(fixture, state);
    fluid::MacVelocityField coMovingVelocity(grid());
    std::ranges::fill(coMovingVelocity.xFaces(), 0.1);
    const auto coMoving = evaluate(
        fixture, state, geometry, coMovingVelocity);
    checkNear(coMoving.totalFluidVolumeFlowRateCubicMetersPerSecond,
              0.1, 5.0e-14,
              "co-moving uniform air retains its analytic opening flow");
    checkNear(coMoving.totalRelativeVolumeFlowRateCubicMetersPerSecond,
              0.0, 5.0e-14,
              "co-moving air and mouth have Galilean-invariant zero relative flow");

    const auto velocity = linearXVelocity();
    const auto flux = evaluate(fixture, state, geometry, velocity);
    check(diagnostics.finite
              && std::ranges::all_of(
                  geometry.patches.patches,
                  [](const auto& patch) {
                      return patch.ownerKind
                          == SceneFluidOpeningPatchOwnerKind::Cell;
                  })
              && flux.velocityEvaluationCount > flux.samples.size(),
          "moving off-face opening uses cell-polygon velocity quadrature");
    checkNear(flux.totalAreaSquareMeters, 1.0, 2.0e-14,
              "moving flux preserves accepted opening area");
    checkNear(flux.totalFluidVolumeFlowRateCubicMetersPerSecond,
              2.01, 8.0e-13,
              "staggered interpolation integrates an analytic linear field");
    checkNear(flux.totalSurfaceSweepRateCubicMetersPerSecond,
              0.1, 5.0e-13,
              "moving flux retains accepted cap sweep");
    checkNear(flux.totalRelativeVolumeFlowRateCubicMetersPerSecond,
              1.91, 1.2e-12,
              "relative opening flow subtracts accepted boundary motion");
}

void testCorruptionBindingAndLimits() {
    Fixture fixture;
    const auto state = fixture.state();
    const auto geometry = buildGeometry(fixture, state);
    const auto velocity = tiledFaceVelocity();
    const auto accepted = evaluate(fixture, state, geometry, velocity);

    auto corrupt = accepted;
    corrupt.samples.front().relativeVolumeFlowRateCubicMetersPerSecond += 0.01;
    expectInvalid(
        [&] { validateSceneFluidOpeningFlux(
            corrupt, fixture.surface.definition, state, geometry.caps,
            geometry.quadrature, geometry.patches, grid(), velocity); },
        "opening flux rejects ledger corruption");

    auto changedVelocity = velocity;
    changedVelocity.yFaces().front() += 0.25;
    expectInvalid(
        [&] { validateSceneFluidOpeningFlux(
            accepted, fixture.surface.definition, state, geometry.caps,
            geometry.quadrature, geometry.patches, grid(),
            changedVelocity); },
        "opening flux binds the complete MAC field");

    auto nonFiniteVelocity = velocity;
    nonFiniteVelocity.yFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(evaluate(
            fixture, state, geometry, nonFiniteVelocity)); },
        "opening flux rejects a non-finite MAC field");

    SceneFluidOpeningFluxLimits limits;
    limits.maximumOpenings = 0;
    expectLimited(
        [&] { static_cast<void>(evaluate(
            fixture, state, geometry, velocity, limits)); },
        "opening flux bounds opening count");
    limits = {};
    limits.maximumPatchSamples = 0;
    expectLimited(
        [&] { static_cast<void>(evaluate(
            fixture, state, geometry, velocity, limits)); },
        "opening flux bounds patch sample count");
    limits = {};
    limits.maximumVelocityEvaluations = 0;
    expectLimited(
        [&] { static_cast<void>(evaluate(
            fixture, state, geometry, velocity, limits)); },
        "opening flux bounds field-evaluation work");
    limits = {};
    limits.maximumFluxBytes = 0;
    expectLimited(
        [&] { static_cast<void>(evaluate(
            fixture, state, geometry, velocity, limits)); },
        "opening flux bounds aggregate storage");
}

} // namespace

int main() {
    try {
        testResolvedFaceFluxAndOrientation();
        testMovingOffFaceRelativeFlux();
        testCorruptionBindingAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d opening-flux check(s) failed\n", failures);
        return 1;
    }
    std::puts("all scene fluid opening-flux checks passed");
    return 0;
}
