#include "surface_aerodynamics.h"

#include <cmath>
#include <cstdio>

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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

Scene tetrahedronScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:surface-aerodynamics-tetrahedron";
    scene.metadata.exporterVersion = "surface-aerodynamics-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {0.0, 0.0, 0.0}},
        {11, {1.0, 0.0, 0.0}},
        {12, {0.0, 1.0, 0.0}},
        {13, {0.0, 0.0, 1.0}},
    };
    scene.fabricMaterials = {
        {100, "test-fabric", 800.0, 600.0, 200.0, 0.01,
         0.04, 0.01, 0.0, 0.0},
    };
    const auto chart = std::array<Vec2, 3>{
        Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}};
    // Every face points from the enclosed cell to the outside.
    scene.triangles = {
        {500, {11, 12, 13}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
        {501, {10, 13, 12}, chart, 2, 1, 100, 901, SurfaceRole::Skin},
        {502, {10, 11, 13}, chart, 2, 1, 100, 902, SurfaceRole::Skin},
        {503, {10, 12, 11}, chart, 2, 1, 100, 903, SurfaceRole::Skin},
    };
    return scene;
}

Scene intakeScene() {
    Scene scene = tetrahedronScene();
    scene.triangles.erase(scene.triangles.begin());
    scene.openings = {
        {600, {11, 12, 13}, 2, 1, OpeningRole::Intake,
         {{{11, 12, 13}}}},
    };
    return scene;
}

SceneFluidSurfaceState referenceState(
    const SceneFluidSurfaceDefinition& surface) {
    SceneFluidSurfaceState state;
    state.definitionFingerprint = surface.fingerprint;
    state.structureDefinitionFingerprint = 1;
    for (const auto& vertex : surface.vertices) {
        state.vertices.push_back(
            {vertex.id, vertex.referencePositionMeters, {}});
    }
    // A captured adapter state normally owns this fingerprint. The model only
    // accepts such canonical states, so generate one through real adapters in
    // the test rather than duplicating the private hash.
    const Scene scene = tetrahedronScene();
    const auto structureAssembly = assembleSceneStructure(scene);
    Structure structure(structureAssembly.definition);
    return captureSceneFluidSurfaceState(
        surface, structureAssembly.mappings, structure);
}

void testClosedCellAndPolarForce() {
    const Scene scene = tetrahedronScene();
    const auto assembled = assembleSceneFluidSurface(scene);
    check(assembled.ok(), "surface aerodynamics: fixture assembles");
    const auto state = referenceState(assembled.definition);

    SurfaceAerodynamicsSettings settings;
    settings.timeStepSeconds = 0.01;
    settings.windRampSeconds = 0.0;
    settings.targetWindMetersPerSecond = {0.0, 10.0, 0.0};
    settings.trimIncidenceRadians = 0.10;
    SurfaceAerodynamicsModel model(assembled.definition, settings);
    const SurfaceAerodynamicsState initial = model.initialState(state);
    check(initial.cells.size() == 1,
          "surface aerodynamics: one authored cell owns one gas state");
    checkNear(initial.cells[0].referenceVolumeCubicMeters,
              1.0 / 6.0, 1.0e-14,
              "surface aerodynamics: closed tetrahedron volume is exact");
    checkNear(initial.cells[0].gaugePressurePascals, 0.0, 1.0e-12,
              "surface aerodynamics: startup begins without a pressure impulse");

    const SurfaceAerodynamicsCandidate candidate = model.advance(initial, state);
    const SurfaceAerodynamicsCandidate repeated = model.advance(initial, state);
    check(candidate == repeated,
          "surface aerodynamics: identical immutable inputs replay exactly");
    check(initial.acceptedStepCount == 0,
          "surface aerodynamics: evaluating a candidate does not mutate input state");
    check(candidate.nextState.acceptedStepCount == 1
              && candidate.triangleTractions.size() == 4
              && candidate.trianglePressureJumpPascals.size() == 4,
          "surface aerodynamics: candidate advances state and covers every triangle");
    check(candidate.nextState.cells[0].gaugePressurePascals > 0.0,
          "surface aerodynamics: bootstrap pressure enters through the ramped candidate");
    check(candidate.diagnostics.finite
              && candidate.diagnostics.dynamicPressurePascals > 0.0
              && candidate.diagnostics.liftNewtons > 0.0
              && candidate.diagnostics.dragNewtons > 0.0,
          "surface aerodynamics: bounded polar produces finite lift and drag");

    StructureVector3 integrated{};
    const auto& vertices = state.vertices;
    for (std::size_t index = 0; index < assembled.definition.triangles.size();
         ++index) {
        const auto& triangle = assembled.definition.triangles[index];
        const Vec3 a = vertices[triangle.vertexIndices[0]].positionMeters;
        const Vec3 b = vertices[triangle.vertexIndices[1]].positionMeters;
        const Vec3 c = vertices[triangle.vertexIndices[2]].positionMeters;
        const double abx = b.x - a.x;
        const double aby = b.y - a.y;
        const double abz = b.z - a.z;
        const double acx = c.x - a.x;
        const double acy = c.y - a.y;
        const double acz = c.z - a.z;
        const double area = 0.5 * std::sqrt(
            std::pow(aby * acz - abz * acy, 2.0)
            + std::pow(abz * acx - abx * acz, 2.0)
            + std::pow(abx * acy - aby * acx, 2.0));
        integrated.x += candidate.triangleTractions[index].tractionPascals.x * area;
        integrated.y += candidate.triangleTractions[index].tractionPascals.y * area;
        integrated.z += candidate.triangleTractions[index].tractionPascals.z * area;
    }
    checkNear(integrated.x,
              candidate.diagnostics.aerodynamicForceNewtons.x, 1.0e-11,
              "surface aerodynamics: closed pressure has zero net X force");
    checkNear(integrated.y,
              candidate.diagnostics.aerodynamicForceNewtons.y, 1.0e-11,
              "surface aerodynamics: integrated drag matches the polar target");
    checkNear(integrated.z,
              candidate.diagnostics.aerodynamicForceNewtons.z, 1.0e-11,
              "surface aerodynamics: integrated lift matches the polar target");
}

void testAuthoredIntakeAddsMass() {
    const Scene scene = intakeScene();
    const auto surface = assembleSceneFluidSurface(scene);
    const auto structureAssembly = assembleSceneStructure(scene);
    check(surface.ok() && structureAssembly.ok(),
          "surface aerodynamics: open-cell fixture assembles");
    Structure structure(structureAssembly.definition);
    const auto moving = captureSceneFluidSurfaceState(
        surface.definition, structureAssembly.mappings, structure);

    SurfaceAerodynamicsSettings settings;
    settings.windRampSeconds = 0.0;
    settings.targetWindMetersPerSecond = {-5.773502691896258,
                                          -5.773502691896258,
                                          -5.773502691896258};
    settings.initialCellPressureDynamicFraction = 0.0;
    SurfaceAerodynamicsModel model(surface.definition, settings);
    const SurfaceAerodynamicsState initial = model.initialState(moving);
    const SurfaceAerodynamicsCandidate candidate = model.advance(initial, moving);
    check(candidate.nextState.cells[0].airMassKilograms
              > initial.cells[0].airMassKilograms
              && candidate.diagnostics.totalOpeningMassFlowKilogramsPerSecond
                  > 0.0,
          "surface aerodynamics: inward ram flow adds mass through an authored intake");
}

} // namespace

int main() {
    testClosedCellAndPolarForce();
    testAuthoredIntakeAddsMass();
    if (failures != 0) {
        std::fprintf(stderr, "%d surface-aerodynamics test(s) failed\n", failures);
        return 1;
    }
    std::puts("surface-aerodynamics tests passed");
    return 0;
}
