#include "scene_fluid_grid_epoch.h"

#include <algorithm>
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

Scene movingTriangleScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-grid-epoch";
    scene.metadata.exporterVersion = "scene-fluid-grid-epoch-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.1, 1.1, 1.9}},
        {11, {2.7, 1.1, 1.9}},
        {12, {2.7, 2.7, 1.9}},
        {13, {1.1, 2.7, 1.9}},
    };
    scene.fabricMaterials = {
        {100, "ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {1.6, 1.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 12, 13},
         {{{0.0, 0.0}, {1.6, 1.6}, {0.0, 1.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

struct Fixture {
    Scene scene = movingTriangleScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structureAssembly = assembleSceneStructure(scene);
    Structure structure{structureAssembly.definition};
    SceneFluidSurfaceTransfer transfer{
        surface.definition, structureAssembly.mappings, structure};
};

std::vector<SceneFluidQuadratureTraction> uniformTractions(
    const SceneFluidQuadratureDefinition& quadrature) {
    std::vector<SceneFluidQuadratureTraction> result;
    result.reserve(quadrature.points.size());
    for (const auto& point : quadrature.points) {
        result.push_back({point.stableId, {10.0, 0.0, 0.0}});
    }
    return result;
}

void testAcceptedGridEpochAndConservativeQuadrature() {
    Fixture fixture;
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "grid epoch fixture assembles from authoritative scene identities");
    const auto initialState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto first = buildSceneFluidGridEpoch(
        fixture.surface.definition, initialState, grid(), fixture.transfer);
    const auto repeated = buildSceneFluidGridEpoch(
        fixture.surface.definition, initialState, grid(), fixture.transfer);
    check(first == repeated
              && first.version == sceneFluidGridEpochVersion
              && first.fingerprint != 0
              && first.ownedStorageBytes != 0
              && first.quadrature.points.size() > 3
              && first.faceLoops.unresolvedOpenChainCount != 0
              && first.facePartitions.unresolvedActiveFaceCount != 0,
          "one immutable epoch composes every deterministic geometry stage and keeps open partitions explicit");

    const auto initialTransfer = evaluateSceneFluidQuadrature(
        fixture.transfer, initialState, first.quadrature,
        uniformTractions(first.quadrature));
    checkNear(
        initialTransfer.diagnostics().surfaceAreaSquareMeters,
        2.56, 6.0e-15,
        "grid epoch quadrature owns the complete initial sheet area once");
    checkNear(
        initialTransfer.diagnostics().integratedSurfaceForceNewtons.x,
        25.6, 6.0e-14,
        "grid epoch quadrature integrates initial force exactly");
    check(initialTransfer.diagnostics().forceResidualNormNewtons < 2.0e-14
              && initialTransfer.diagnostics()
                     .momentResidualNormNewtonMeters < 2.0e-14,
          "initial epoch closes conservative force and moment transfer");

    const std::array<double, 4> accelerations{1.4, 1.0, 0.6, 1.0};
    for (std::size_t node = 0;
         node < fixture.structure.definition().nodes.size(); ++node) {
        fixture.structure.addExternalForce(
            node,
            {0.0, 0.0,
             accelerations[node]
                 * fixture.structure.definition().nodes[node].massKg});
    }
    StructureStepSettings step;
    step.timeStepSeconds = 0.5;
    step.substeps = 1;
    step.constraintIterations = 4;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto structuralDiagnostics = fixture.structure.step(step);
    const auto movedState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto moved = buildSceneFluidGridEpoch(
        fixture.surface.definition, movedState, grid(), fixture.transfer);

    double minimumZ = movedState.vertices.front().positionMeters.z;
    double maximumZ = minimumZ;
    for (const auto& vertex : movedState.vertices) {
        minimumZ = std::min(minimumZ, vertex.positionMeters.z);
        maximumZ = std::max(maximumZ, vertex.positionMeters.z);
    }
    check(structuralDiagnostics.finite
              && movedState.acceptedStepCount == 1
              && minimumZ > 2.0
              && maximumZ - minimumZ > 1.0e-4
              && moved.fingerprint != first.fingerprint
              && moved.candidates.fingerprint
                     != first.candidates.fingerprint
              && moved.ownership.fingerprint != first.ownership.fingerprint
              && moved.faceTopology.fingerprint
                     != first.faceTopology.fingerprint
              && moved.quadrature.fingerprint
                     != first.quadrature.fingerprint,
          "accepted non-rigid fabric motion across a grid plane remaps the complete bound epoch");
    for (const auto& point : moved.quadrature.points) {
        check((point.triangleId == 500 || point.triangleId == 501)
                  && point.negativeSideRegionId == 1
                  && point.positiveSideRegionId == 2
                  && point.materialId == 100
                  && point.sheetId == 900,
              "moving remap preserves authored triangle, side, material, and sheet identity");
    }

    const auto movedTransfer = evaluateSceneFluidQuadrature(
        fixture.transfer, movedState, moved.quadrature,
        uniformTractions(moved.quadrature));
    checkNear(
        movedTransfer.diagnostics().integratedSurfaceForceNewtons.x,
        10.0 * movedTransfer.diagnostics().surfaceAreaSquareMeters,
        6.0e-14,
        "grid crossing integrates the deformed unique area exactly");
    check(movedTransfer.diagnostics().surfaceAreaSquareMeters > 2.5
              && movedTransfer.diagnostics().surfaceAreaSquareMeters < 2.7,
          "remapped quadrature retains the complete deformed sheet area");
    check(movedTransfer.diagnostics().forceResidualNormNewtons < 2.0e-14
              && movedTransfer.diagnostics()
                     .momentResidualNormNewtonMeters < 2.0e-14,
          "remapped quadrature preserves force and moment closure");
    validateSceneFluidGridEpoch(
        moved, fixture.surface.definition, movedState, grid(),
        fixture.transfer);
}

void testTransactionalEpochRejectionAndBounds() {
    Fixture fixture;
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto accepted = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer);

    auto corrupt = accepted;
    ++corrupt.crossings.crossings.front().triangleId;
    expectInvalid(
        [&] { validateSceneFluidGridEpoch(
            corrupt, fixture.surface.definition, state, grid(),
            fixture.transfer); },
        "grid epoch rejects nested stage corruption before publication");

    auto foreignState = state;
    ++foreignState.acceptedStepCount;
    expectInvalid(
        [&] { validateSceneFluidGridEpoch(
            accepted, fixture.surface.definition, foreignState, grid(),
            fixture.transfer); },
        "grid epoch cannot mix geometry from different accepted steps");

    SceneFluidGridEpochLimits limits;
    limits.candidates.maximumCandidates = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidGridEpoch(
            fixture.surface.definition, state, grid(), fixture.transfer,
            {}, limits)); },
        "grid epoch preserves early per-stage allocation bounds");
    limits = {};
    limits.maximumEpochBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidGridEpoch(
            fixture.surface.definition, state, grid(), fixture.transfer,
            {}, limits)); },
        "grid epoch bounds aggregate owned payload before return");
    validateSceneFluidGridEpoch(
        accepted, fixture.surface.definition, state, grid(),
        fixture.transfer);
}

} // namespace

int main() {
    try {
        testAcceptedGridEpochAndConservativeQuadrature();
        testTransactionalEpochRejectionAndBounds();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid grid-epoch check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid grid-epoch checks passed");
    return 0;
}
