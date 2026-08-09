#include "scene_fluid_pressure_coupling.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <ranges>
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

double norm(const StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
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

Scene tetraScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-coupling";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-coupling-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 4> positions{{
        {1.2, 1.5, 1.45},
        {2.0, 1.2, 1.15},
        {2.0, 1.8, 1.15},
        {2.0, 1.5, 1.75},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 3> faces{{
        {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}},
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
        {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

SceneStructureAssembly fixedMouthAssembly(const Scene& scene) {
    auto assembly = assembleSceneStructure(scene);
    for (std::size_t node = 0;
         node < assembly.mappings.nodeVertexIds.size(); ++node) {
        if (assembly.mappings.nodeVertexIds[node] != 10) {
            assembly.definition.nodes[node].fixed = true;
        }
    }
    return assembly;
}

SceneFluidPressureCouplingSettings couplingSettings() {
    SceneFluidPressureCouplingSettings settings;
    settings.structure.timeStepSeconds = 0.05;
    settings.structure.substeps = 1;
    settings.structure.constraintIterations = 0;
    settings.structure.gravityMetersPerSecondSquared = {-0.8, 0.0, 0.0};
    settings.structure.velocityDampingPerSecond = 0.0;
    settings.pressureProjection.timeStepSeconds =
        settings.structure.timeStepSeconds;
    settings.pressureProjection
        .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond = 2.0e-11;
    settings.pressureProjection.relativeCorrectedVolumeRateTolerance =
        1.0e-11;
    settings.pressureProjection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.pressureProjection.pressureSolve.relativeResidualTolerance =
        1.0e-13;
    settings.pressureProjection.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    settings.relaxation.initialRelaxation = 0.5;
    settings.relaxation.minimumRelaxation = 0.02;
    settings.relaxation.maximumRelaxation = 1.0;
    settings.convergence.minimumIterations = 2;
    settings.convergence.maximumIterations = 20;
    settings.convergence.absoluteDisplacementToleranceMetres = 1.0e-9;
    settings.convergence.relativeDisplacementTolerance = 1.0e-7;
    settings.convergence.absoluteVelocityToleranceMetersPerSecond = 1.0e-8;
    settings.convergence.relativeVelocityTolerance = 1.0e-7;
    settings.convergence.absoluteTractionToleranceNewtons = 1.0e-8;
    settings.convergence.relativeTractionTolerance = 1.0e-7;
    return settings;
}

struct Fixture {
    Scene scene = tetraScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly assembly = fixedMouthAssembly(scene);
    Structure structure{assembly.definition};
};

std::size_t apexNode(const Fixture& fixture) {
    const auto found = std::ranges::find(
        fixture.assembly.mappings.nodeVertexIds, StableId{10});
    if (found == fixture.assembly.mappings.nodeVertexIds.end()) {
        throw std::runtime_error("pressure coupling apex is missing");
    }
    return static_cast<std::size_t>(
        found - fixture.assembly.mappings.nodeVertexIds.begin());
}

bool samePublicCheckpoint(const StructureCheckpoint& first,
                          const StructureCheckpoint& second) {
    return first.version == second.version
        && first.definitionFingerprint == second.definitionFingerprint
        && first.acceptedStepCount == second.acceptedStepCount
        && first.simulationTimeSeconds == second.simulationTimeSeconds
        && first.nodes == second.nodes
        && first.pendingExternalForcesNewtons
            == second.pendingExternalForcesNewtons
        && first.lastAppliedExternalForceNewtons
            == second.lastAppliedExternalForceNewtons;
}

void testStrongPressureFeedbackConvergesDeterministically() {
    Fixture first;
    Fixture second;
    const auto settings = couplingSettings();
    const fluid::PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    SceneFluidPressureCoupling firstCoupling(
        first.surface.definition, first.assembly.mappings,
        first.structure, grid, settings);
    SceneFluidPressureCoupling secondCoupling(
        second.surface.definition, second.assembly.mappings,
        second.structure, grid, settings);
    fluid::MacVelocityField velocity(grid);

    Structure control(first.assembly.definition);
    const auto controlDiagnostics = control.step(settings.structure);
    const double controlApexX =
        control.nodeStates()[apexNode(first)].positionMeters.x;

    const auto firstResult = firstCoupling.advance(first.structure, velocity);
    const auto secondResult = secondCoupling.advance(
        second.structure, velocity);
    const auto firstState = first.structure.checkpoint();
    const auto secondState = second.structure.checkpoint();
    const auto* projection = firstCoupling.acceptedPressureProjection();
    const double coupledApexX =
        first.structure.nodeStates()[apexNode(first)].positionMeters.x;
    const auto& transfer = firstCoupling.acceptedPressureTransfer();

    check(firstResult == secondResult
              && samePublicCheckpoint(firstState, secondState)
              && firstCoupling.acceptedPressureEpoch()
                  == secondCoupling.acceptedPressureEpoch()
              && firstCoupling.acceptedPressureTransfer()
                  == secondCoupling.acceptedPressureTransfer(),
          "strong scene pressure feedback is bit-deterministic");
    check(firstResult.accepted
              && firstResult.iteration.status
                  == StrongCouplingIterationStatus::Converged
              && firstResult.solverRunCount >= 2
              && firstResult.solverRunCount
                  <= settings.convergence.maximumIterations
              && first.structure.acceptedStepCount() == 1
              && first.structure.simulationTimeSeconds()
                  == settings.structure.timeStepSeconds,
          "strong scene pressure feedback commits one converged Structure step");
    check(projection != nullptr
              && projection->diagnostics.accepted
              && projection->diagnostics.usesMovingVolumeRates
              && projection->diagnostics
                  .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond
                  > 1.0e-5
              && projection->diagnostics
                  .correctedContinuityResidualMaximumCubicMetersPerSecond
                  < 2.0e-11,
          "accepted feedback retains its moving-volume pressure projection");
    check(transfer.diagnostics().finite
              && norm(transfer.diagnostics().integratedSurfaceForceNewtons)
                  > 1.0e-5
              && transfer.diagnostics().forceResidualNormNewtons < 1.0e-12
              && transfer.diagnostics().momentResidualNormNewtonMeters
                  < 1.0e-12
              && firstResult.interfaceForceClosureNewtons
                  <= settings.convergence
                         .absoluteTractionToleranceNewtons,
          "accepted feedback closes pressure load, force, moment, and interface iteration");
    check(controlDiagnostics.finite
              && coupledApexX > controlApexX
              && coupledApexX - controlApexX > 1.0e-8,
          "scene-derived pressure opposes the gravity-driven cell expansion");

    const std::uint64_t firstEpochFingerprint =
        firstCoupling.acceptedPressureEpoch().fingerprint;
    const auto firstContinuation = firstCoupling.advance(
        first.structure, velocity);
    const auto secondContinuation = secondCoupling.advance(
        second.structure, velocity);
    check(firstContinuation == secondContinuation
              && firstContinuation.accepted
              && firstContinuation.previousPressureEpochFingerprint
                  == firstEpochFingerprint
              && firstContinuation.currentPressureEpochFingerprint
                  == firstCoupling.acceptedPressureEpoch().fingerprint
              && firstContinuation.currentPressureEpochFingerprint
                  != firstEpochFingerprint
              && first.structure.acceptedStepCount() == 2
              && first.structure.checkpoint().nodes
                  == second.structure.checkpoint().nodes,
          "accepted pressure and warm state continue into the exact next coupled macro-step");
}

void testExhaustionAndProjectionFailureRollback() {
    Fixture exhaustedFixture;
    auto exhaustedSettings = couplingSettings();
    exhaustedSettings.convergence.minimumIterations = 1;
    exhaustedSettings.convergence.maximumIterations = 1;
    exhaustedSettings.convergence.absoluteDisplacementToleranceMetres =
        1.0e-30;
    exhaustedSettings.convergence.relativeDisplacementTolerance = 1.0e-30;
    exhaustedSettings.convergence
        .absoluteVelocityToleranceMetersPerSecond = 1.0e-30;
    exhaustedSettings.convergence.relativeVelocityTolerance = 1.0e-30;
    exhaustedSettings.convergence.absoluteTractionToleranceNewtons =
        1.0e-30;
    exhaustedSettings.convergence.relativeTractionTolerance = 1.0e-30;
    const fluid::PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    SceneFluidPressureCoupling exhaustedCoupling(
        exhaustedFixture.surface.definition,
        exhaustedFixture.assembly.mappings,
        exhaustedFixture.structure, grid, exhaustedSettings);
    exhaustedFixture.structure.addExternalForce(
        apexNode(exhaustedFixture), {0.01, 0.02, 0.03});
    const auto exhaustedBefore = exhaustedFixture.structure.checkpoint();
    const auto exhaustedEpoch =
        exhaustedCoupling.acceptedPressureEpoch().fingerprint;
    fluid::MacVelocityField velocity(grid);
    const auto exhausted = exhaustedCoupling.advance(
        exhaustedFixture.structure, velocity);
    check(!exhausted.accepted
              && exhausted.iteration.status
                  == StrongCouplingIterationStatus::Exhausted
              && exhausted.solverRunCount == 1
              && samePublicCheckpoint(
                  exhaustedFixture.structure.checkpoint(), exhaustedBefore)
              && exhaustedCoupling.acceptedPressureEpoch().fingerprint
                  == exhaustedEpoch
              && exhaustedCoupling.acceptedPressureProjection() == nullptr,
          "iteration exhaustion restores Structure and preserves the accepted pressure owner");

    Fixture failedFixture;
    auto failedSettings = couplingSettings();
    failedSettings.pressureProjection.pressureSolve.maximumIterations = 0;
    failedSettings.pressureProjection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-30;
    failedSettings.pressureProjection.pressureSolve
        .relativeResidualTolerance = 0.0;
    SceneFluidPressureCoupling failedCoupling(
        failedFixture.surface.definition, failedFixture.assembly.mappings,
        failedFixture.structure, grid, failedSettings);
    const auto failedBefore = failedFixture.structure.checkpoint();
    bool rejected = false;
    try {
        static_cast<void>(failedCoupling.advance(
            failedFixture.structure, velocity));
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check(rejected
              && samePublicCheckpoint(
                  failedFixture.structure.checkpoint(), failedBefore)
              && failedCoupling.acceptedPressureProjection() == nullptr
              && failedCoupling.acceptedSurfaceState().acceptedStepCount == 0,
          "projection failure restores the exact Structure and pressure baseline");

    Fixture topologyFixture;
    auto topologySettings = couplingSettings();
    topologySettings.structure.gravityMetersPerSecondSquared = {
        -800.0, 0.0, 0.0};
    SceneFluidPressureCoupling topologyCoupling(
        topologyFixture.surface.definition,
        topologyFixture.assembly.mappings,
        topologyFixture.structure, grid, topologySettings);
    const auto topologyBefore = topologyFixture.structure.checkpoint();
    bool topologyRejected = false;
    try {
        static_cast<void>(topologyCoupling.advance(
            topologyFixture.structure, velocity));
    } catch (const std::exception&) {
        topologyRejected = true;
    }
    check(topologyRejected
              && samePublicCheckpoint(
                  topologyFixture.structure.checkpoint(), topologyBefore)
              && topologyCoupling.acceptedPressureProjection() == nullptr
              && topologyCoupling.acceptedPressureEpoch().acceptedStepCount
                  == 0,
          "out-of-subset moving topology rejects without publishing a partial coupled state");
}

void testCouplingInterfaceLimits() {
    Fixture fixture;
    const fluid::PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    SceneFluidPressureCouplingLimits limits;
    limits.maximumCouplingNodes = 3;
    bool rejected = false;
    try {
        SceneFluidPressureCoupling coupling(
            fixture.surface.definition, fixture.assembly.mappings,
            fixture.structure, grid, couplingSettings(), limits);
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected && fixture.structure.acceptedStepCount() == 0,
          "scene pressure coupling bounds its nonlinear interface storage");
}

} // namespace

int main() {
    try {
        testStrongPressureFeedbackConvergesDeterministically();
        testExhaustionAndProjectionFailureRollback();
        testCouplingInterfaceLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene pressure coupling check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene pressure coupling checks passed");
    return 0;
}
