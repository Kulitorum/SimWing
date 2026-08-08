#include "scene_fluid_surface_transfer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

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

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

Scene surfaceScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-surface-transfer";
    scene.metadata.exporterVersion = "scene-fluid-surface-transfer-test/1";
    scene.regions = {
        {2, RegionKind::Cell, "cell"},
        {1, RegionKind::Outside, "outside"},
    };
    scene.vertices = {
        {12, {1.0, 1.0, 0.0}},
        {10, {0.0, 0.0, 0.0}},
        {13, {0.0, 1.0, 0.0}},
        {11, {1.0, 0.0, 0.0}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {501, {10, 12, 13}, {{{0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {500, {10, 11, 12}, {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    scene.openings = {
        {600, {10, 11, 12, 13}, 1, 2, OpeningRole::Intake},
    };
    return scene;
}

void reverseRelevantCollections(Scene& scene) {
    std::reverse(scene.regions.begin(), scene.regions.end());
    std::reverse(scene.vertices.begin(), scene.vertices.end());
    std::reverse(scene.fabricMaterials.begin(), scene.fabricMaterials.end());
    std::reverse(scene.triangles.begin(), scene.triangles.end());
    std::reverse(scene.openings.begin(), scene.openings.end());
}

void testAuthoritativeTopologyAndUniformTransfer() {
    const Scene scene = surfaceScene();
    const auto surfaceAssembly = assembleSceneFluidSurface(scene);
    const auto structureAssembly = assembleSceneStructure(scene);
    check(surfaceAssembly.ok() && structureAssembly.ok(),
          "scene surface transfer: both source adapters accept the fixture");
    Structure structure(structureAssembly.definition);
    SceneFluidSurfaceTransfer transfer(
        surfaceAssembly.definition,
        structureAssembly.mappings,
        structure);

    check(transfer.surfaceDefinitionFingerprint()
              == surfaceAssembly.definition.fingerprint
              && transfer.couplingSurfaceFingerprint() != 0
              && transfer.targetDefinitionFingerprint()
                  == structure.definitionFingerprint(),
          "scene surface transfer: all three immutable identities are exposed");
    check(transfer.nodes().size() == 4
              && transfer.nodes()[0]
                  == CouplingSurfaceNodeDefinition{10, 0}
              && transfer.nodes()[3]
                  == CouplingSurfaceNodeDefinition{13, 3},
          "scene surface transfer: canonical vertices map one-to-one to Structure");
    check(transfer.triangles().size() == 2
              && transfer.triangles()[0]
                  == CouplingSurfaceTriangleDefinition{500, {10, 11, 12}}
              && transfer.triangles()[1]
                  == CouplingSurfaceTriangleDefinition{501, {10, 12, 13}},
          "scene surface transfer: authored oriented triangles reach coupling topology");

    const auto state = captureSceneFluidSurfaceState(
        surfaceAssembly.definition,
        structureAssembly.mappings,
        structure);
    const std::vector<CouplingTriangleTraction> tractions = {
        {500, {0.0, 0.0, 10.0}},
        {501, {0.0, 0.0, 10.0}},
    };
    const auto result = transfer.evaluate(state, tractions);
    const auto& diagnostics = result.diagnostics();
    checkNear(diagnostics.surfaceAreaSquareMeters, 1.0, 1.0e-15,
              "scene surface transfer: current triangle area is exact");
    checkNear(diagnostics.integratedSurfaceForceNewtons.z, 10.0, 1.0e-14,
              "scene surface transfer: uniform pressure traction integrates exactly");
    check(diagnostics.forceResidualNormNewtons < 1.0e-14
              && diagnostics.momentResidualNormNewtonMeters < 1.0e-14
              && std::abs(diagnostics.powerResidualWatts) < 1.0e-14
              && diagnostics.finite,
          "scene surface transfer: conservative force, moment, and power ledgers close");
    check(result.nodeLoads().size() == 4,
          "scene surface transfer: every authoritative vertex receives one load");
    checkNear(result.nodeLoads()[0].forceNewtons.z, 10.0 / 3.0, 1.0e-14,
              "scene surface transfer: shared first vertex receives both shares");
    checkNear(result.nodeLoads()[1].forceNewtons.z, 5.0 / 3.0, 1.0e-14,
              "scene surface transfer: boundary vertex receives one share");
    checkNear(result.nodeLoads()[2].forceNewtons.z, 10.0 / 3.0, 1.0e-14,
              "scene surface transfer: shared diagonal vertex receives both shares");
    checkNear(result.nodeLoads()[3].forceNewtons.z, 5.0 / 3.0, 1.0e-14,
              "scene surface transfer: final boundary vertex receives one share");

    transfer.addLoadsTo(structure, result);
    const auto loaded = structure.diagnostics();
    checkNear(loaded.pendingExternalForceNewtons.z, 10.0, 1.0e-14,
              "scene surface transfer: validated result enters pending XPBD loads");

    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = 0.01;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    const auto stepDiagnostics = structure.step(stepSettings);
    const auto movedState = captureSceneFluidSurfaceState(
        surfaceAssembly.definition,
        structureAssembly.mappings,
        structure);
    const auto movedResult = transfer.evaluate(movedState, tractions);
    check(stepDiagnostics.finite
              && movedState.acceptedStepCount == 1
              && movedResult.diagnostics().integratedSurfacePowerWatts > 0.0
              && std::abs(movedResult.diagnostics().powerResidualWatts)
                  < 1.0e-13,
          "scene surface transfer: accepted XPBD motion closes the next power transfer");
}

void testQuadratureAndDeterminism() {
    Scene firstScene = surfaceScene();
    Scene secondScene = firstScene;
    reverseRelevantCollections(secondScene);
    const auto firstSurface = assembleSceneFluidSurface(firstScene);
    const auto secondSurface = assembleSceneFluidSurface(secondScene);
    const auto firstStructure = assembleSceneStructure(firstScene);
    const auto secondStructure = assembleSceneStructure(secondScene);
    Structure first(firstStructure.definition);
    Structure second(secondStructure.definition);
    SceneFluidSurfaceTransfer firstTransfer(
        firstSurface.definition, firstStructure.mappings, first);
    SceneFluidSurfaceTransfer secondTransfer(
        secondSurface.definition, secondStructure.mappings, second);
    check(firstTransfer.surfaceDefinitionFingerprint()
              == secondTransfer.surfaceDefinitionFingerprint()
              && firstTransfer.couplingSurfaceFingerprint()
                  == secondTransfer.couplingSurfaceFingerprint(),
          "scene surface transfer: source collection order cannot change topology identity");

    const auto state = captureSceneFluidSurfaceState(
        firstSurface.definition, firstStructure.mappings, first);
    const std::vector<CouplingTriangleTractionQuadrature> quadrature = {
        {700, 500, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
         0.5, {0.0, 0.0, 10.0}},
        {701, 501, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
         0.5, {0.0, 0.0, 10.0}},
    };
    const auto result = firstTransfer.evaluateQuadrature(state, quadrature);
    checkNear(result.diagnostics().integratedSurfaceForceNewtons.z,
              10.0, 1.0e-14,
              "scene surface transfer: face-resolved quadrature delegates exactly");
    check(result.diagnostics().forceResidualNormNewtons < 1.0e-14
              && result.diagnostics().momentResidualNormNewtonMeters < 1.0e-14,
          "scene surface transfer: delegated quadrature remains conservative");
}

void testStateAndGeometryBinding() {
    const Scene scene = surfaceScene();
    const auto surfaceAssembly = assembleSceneFluidSurface(scene);
    const auto structureAssembly = assembleSceneStructure(scene);
    Structure structure(structureAssembly.definition);
    SceneFluidSurfaceTransfer transfer(
        surfaceAssembly.definition,
        structureAssembly.mappings,
        structure);
    const std::vector<CouplingTriangleTraction> tractions = {
        {500, {0.0, 0.0, 1.0}},
        {501, {0.0, 0.0, 1.0}},
    };
    const auto state = captureSceneFluidSurfaceState(
        surfaceAssembly.definition,
        structureAssembly.mappings,
        structure);

    auto invalid = state;
    ++invalid.version;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(invalid, tractions)); },
        "scene surface transfer: state version mismatch is rejected");
    invalid = state;
    ++invalid.definitionFingerprint;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(invalid, tractions)); },
        "scene surface transfer: foreign surface identity is rejected");
    invalid = state;
    ++invalid.structureDefinitionFingerprint;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(invalid, tractions)); },
        "scene surface transfer: foreign Structure identity is rejected");
    invalid = state;
    std::reverse(invalid.vertices.begin(), invalid.vertices.end());
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(invalid, tractions)); },
        "scene surface transfer: reordered state vertices are rejected");
    invalid = state;
    invalid.vertices.front().velocityMetersPerSecond.x =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(invalid, tractions)); },
        "scene surface transfer: non-finite state is rejected");

    Scene foreignScene = scene;
    foreignScene.vertices.front().positionMeters.z = 0.1;
    const auto foreignSurface = assembleSceneFluidSurface(foreignScene);
    check(foreignSurface.ok(),
          "scene surface transfer: foreign geometry fixture remains valid");
    expectRejected(
        [&] {
            SceneFluidSurfaceTransfer foreignTransfer(
                foreignSurface.definition,
                structureAssembly.mappings,
                structure);
        },
        "scene surface transfer: constructor rejects foreign reference geometry");
}

} // namespace

int main() {
    testAuthoritativeTopologyAndUniformTransfer();
    testQuadratureAndDeterminism();
    testStateAndGeometryBinding();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid-surface transfer test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene fluid-surface transfer tests passed\n");
    return 0;
}
