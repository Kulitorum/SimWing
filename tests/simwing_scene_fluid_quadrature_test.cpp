#include "scene_fluid_quadrature.h"

#include <cmath>
#include <cstdio>
#include <limits>
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

Scene triangleScene(const double zMeters = 1.5) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-quadrature";
    scene.metadata.exporterVersion = "scene-fluid-quadrature-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.1, 1.1, zMeters}},
        {11, {2.7, 1.1, zMeters}},
        {12, {1.1, 2.7, zMeters}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {0.0, 1.6}}},
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
    SceneFluidSurfaceTransfer transfer;
    SceneFluidSurfaceState state;
    SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections;
    SceneFluidGridPatchSet patches;
    SceneFluidPatchOwnership ownership;
    SceneFluidQuadratureDefinition quadrature;

    explicit Pipeline(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(surface.definition, structureAssembly.mappings, structure),
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
              intersections)),
          ownership(ownSceneFluidSurfacePatches(
              surface.definition,
              state,
              grid(),
              candidates,
              intersections,
              patches)),
          quadrature(buildSceneFluidQuadrature(
              surface.definition,
              state,
              grid(),
              candidates,
              intersections,
              patches,
              ownership,
              transfer)) {}
};

std::vector<SceneFluidQuadratureTraction> uniformTractions(
    const SceneFluidQuadratureDefinition& definition,
    const StructureVector3 traction = {0.0, 0.0, 10.0}) {
    std::vector<SceneFluidQuadratureTraction> result;
    result.reserve(definition.points.size());
    for (const auto& point : definition.points) {
        result.push_back({point.stableId, traction});
    }
    return result;
}

void testCellQuadratureTransfer() {
    Pipeline pipeline(triangleScene());
    const auto rebuilt = buildSceneFluidQuadrature(
        pipeline.surface.definition,
        pipeline.state,
        grid(),
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        pipeline.transfer);
    check(pipeline.quadrature == rebuilt
              && pipeline.quadrature.fingerprint != 0
              && pipeline.quadrature.points.size() == 3,
          "scene fluid quadrature: cell-owned points are deterministic and bound");
    double area = 0.0;
    for (const auto& point : pipeline.quadrature.points) {
        area += point.areaSquareMeters;
        check(point.triangleId == 500
                  && point.negativeSideRegionId == 1
                  && point.positiveSideRegionId == 2
                  && point.materialId == 100
                  && point.sheetId == 900
                  && point.ownerKind == SceneFluidQuadratureOwnerKind::Cell
                  && point.negativeSideCellIndex
                      == point.positiveSideCellIndex,
              "scene fluid quadrature: authored physical ownership reaches every point");
    }
    checkNear(area, 1.28, 3.0e-15,
              "scene fluid quadrature: point areas partition the triangle");

    const auto tractions = uniformTractions(pipeline.quadrature);
    const auto result = evaluateSceneFluidQuadrature(
        pipeline.transfer,
        pipeline.state,
        pipeline.quadrature,
        tractions);
    checkNear(result.diagnostics().surfaceAreaSquareMeters,
              1.28, 3.0e-15,
              "scene fluid quadrature: conservative transfer sees unique area");
    checkNear(result.diagnostics().integratedSurfaceForceNewtons.z,
              12.8, 3.0e-14,
              "scene fluid quadrature: uniform traction integrates exactly");
    check(result.diagnostics().forceResidualNormNewtons < 2.0e-14
              && result.diagnostics().momentResidualNormNewtonMeters < 2.0e-14
              && std::abs(result.diagnostics().powerResidualWatts) < 1.0e-14,
          "scene fluid quadrature: conservative force, moment, and power close");
    pipeline.transfer.addLoadsTo(pipeline.structure, result);
    checkNear(pipeline.structure.diagnostics().pendingExternalForceNewtons.z,
              12.8, 3.0e-14,
              "scene fluid quadrature: unique loads enter the real Structure");
}

void testFaceQuadratureDoesNotDoubleCount() {
    Pipeline boundary(triangleScene(2.0));
    check(boundary.ownership.facePatches.size() == 3
              && boundary.quadrature.points.size() == 3,
          "scene fluid quadrature: three paired grid faces create three points");
    for (const auto& point : boundary.quadrature.points) {
        check(point.ownerKind == SceneFluidQuadratureOwnerKind::Face
                  && point.positiveSideCellIndex
                      == point.negativeSideCellIndex + 16,
              "scene fluid quadrature: shared-plane ownership retains exact side cells");
    }
    const auto result = evaluateSceneFluidQuadrature(
        boundary.transfer,
        boundary.state,
        boundary.quadrature,
        uniformTractions(boundary.quadrature));
    checkNear(result.diagnostics().surfaceAreaSquareMeters,
              1.28, 3.0e-15,
              "scene fluid quadrature: paired cell patches are integrated once");
    checkNear(result.diagnostics().integratedSurfaceForceNewtons.z,
              12.8, 3.0e-14,
              "scene fluid quadrature: shared-plane force is not doubled");
}

void testTransactionalBindingAndTractionValidation() {
    Pipeline pipeline(triangleScene());
    const auto accepted = uniformTractions(pipeline.quadrature);
    auto corruptDefinition = pipeline.quadrature;
    ++corruptDefinition.fingerprint;
    expectInvalid(
        [&] { static_cast<void>(evaluateSceneFluidQuadrature(
            pipeline.transfer,
            pipeline.state,
            corruptDefinition,
            accepted)); },
        "scene fluid quadrature: corrupt definition identity is rejected");

    auto corruptTractions = accepted;
    std::swap(corruptTractions[0], corruptTractions[1]);
    expectInvalid(
        [&] { static_cast<void>(evaluateSceneFluidQuadrature(
            pipeline.transfer,
            pipeline.state,
            pipeline.quadrature,
            corruptTractions)); },
        "scene fluid quadrature: reordered tractions are rejected");
    corruptTractions = accepted;
    corruptTractions.front().tractionPascals.x =
        std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(evaluateSceneFluidQuadrature(
            pipeline.transfer,
            pipeline.state,
            pipeline.quadrature,
            corruptTractions)); },
        "scene fluid quadrature: non-finite traction is rejected");

    auto foreignState = pipeline.state;
    ++foreignState.acceptedStepCount;
    expectInvalid(
        [&] { static_cast<void>(evaluateSceneFluidQuadrature(
            pipeline.transfer,
            foreignState,
            pipeline.quadrature,
            accepted)); },
        "scene fluid quadrature: foreign accepted epoch is rejected");
}

void testAcceptedKinematicsSampling() {
    Pipeline pipeline(triangleScene());
    const auto samples = sampleSceneFluidQuadratureKinematics(
        pipeline.surface.definition,
        pipeline.state,
        pipeline.quadrature);
    check(samples.size() == pipeline.quadrature.points.size(),
          "scene fluid quadrature: every point receives accepted kinematics");
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& sample = samples[index];
        check(sample.stableId == pipeline.quadrature.points[index].stableId,
              "scene fluid quadrature: kinematics retain stable point order");
        checkNear(sample.positionMeters.z, 1.5, 5.0e-16,
                  "scene fluid quadrature: barycentric position stays on the surface");
        checkNear(sample.velocityMetersPerSecond.x, 0.0, 0.0,
                  "scene fluid quadrature: accepted zero X velocity is exact");
        checkNear(sample.velocityMetersPerSecond.y, 0.0, 0.0,
                  "scene fluid quadrature: accepted zero Y velocity is exact");
        checkNear(sample.velocityMetersPerSecond.z, 0.0, 0.0,
                  "scene fluid quadrature: accepted zero Z velocity is exact");
    }
    auto foreign = pipeline.state;
    ++foreign.acceptedStepCount;
    expectInvalid(
        [&] { static_cast<void>(sampleSceneFluidQuadratureKinematics(
            pipeline.surface.definition, foreign, pipeline.quadrature)); },
        "scene fluid quadrature: foreign kinematic epoch is rejected");
}

} // namespace

int main() {
    testCellQuadratureTransfer();
    testFaceQuadratureDoesNotDoubleCount();
    testTransactionalBindingAndTractionValidation();
    testAcceptedKinematicsSampling();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid-quadrature test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene fluid-quadrature tests passed\n");
    return 0;
}
