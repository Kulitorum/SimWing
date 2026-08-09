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
    const auto firstMac =
        firstCoupling.acceptedPressureCorrectedMacVelocity();
    const auto secondMac =
        secondCoupling.acceptedPressureCorrectedMacVelocity();

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
    check(firstMac == secondMac
              && firstMac.diagnostics.finite
              && firstMac.diagnostics.acceptedStepCount == 1
              && firstMac.diagnostics.simulationTimeSeconds
                  == settings.structure.timeStepSeconds
              && firstMac.diagnostics.pressureProjectionFingerprint
                  == projection->fingerprint
              && firstMac.diagnostics.openingPatchFingerprint
                  == firstCoupling.acceptedPressureEpoch()
                         .openingPatches.fingerprint
              && firstMac.diagnostics.faceCount
                  == 3 * grid.cellCount()
              && firstMac.diagnostics.linkCount
                  == projection->links.size()
              && firstMac.diagnostics.openingLinkCount > 0
              && firstMac.diagnostics.multiLinkFaceCount > 0
              && firstMac.diagnostics
                     .maximumAbsoluteVelocityMetersPerSecond > 1.0e-6
              && firstMac.diagnostics
                     .maximumVolumeFlowClosureCubicMetersPerSecond
                  < 1.0e-15,
          "accepted link flow collapses deterministically onto a conservative bulk MAC continuation");
    check(controlDiagnostics.finite
              && coupledApexX > controlApexX
              && coupledApexX - controlApexX > 1.0e-8,
          "scene-derived pressure opposes the gravity-driven cell expansion");

    const std::uint64_t firstEpochFingerprint =
        firstCoupling.acceptedPressureEpoch().fingerprint;
    const auto firstContinuation = firstCoupling.advance(
        first.structure, firstMac.velocityMetersPerSecond);
    const auto secondContinuation = secondCoupling.advance(
        second.structure, secondMac.velocityMetersPerSecond);
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

    Fixture appearanceFixture;
    auto appearanceSettings = couplingSettings();
    appearanceSettings.structure.gravityMetersPerSecondSquared = {
        -100.0, 0.0, 0.0};
    appearanceSettings.convergence.minimumIterations = 1;
    appearanceSettings.convergence.maximumIterations = 1;
    appearanceSettings.convergence.absoluteDisplacementToleranceMetres =
        1.0e6;
    appearanceSettings.convergence.relativeDisplacementTolerance = 1.0e6;
    appearanceSettings.convergence
        .absoluteVelocityToleranceMetersPerSecond = 1.0e6;
    appearanceSettings.convergence.relativeVelocityTolerance = 1.0e6;
    appearanceSettings.convergence.absoluteTractionToleranceNewtons = 1.0e6;
    appearanceSettings.convergence.relativeTractionTolerance = 1.0e6;
    SceneFluidPressureCoupling appearanceCoupling(
        appearanceFixture.surface.definition,
        appearanceFixture.assembly.mappings,
        appearanceFixture.structure, grid, appearanceSettings);
    const auto appearance = appearanceCoupling.advance(
        appearanceFixture.structure, velocity);
    check(appearance.accepted
              && appearanceCoupling.acceptedPressureEpoch()
                     .pressureControlVolumes.controlVolumes.size() == 66
              && appearanceCoupling.acceptedPressureProjection() != nullptr
              && appearanceCoupling.acceptedPressureProjection()
                     ->diagnostics.usesMovingVolumeRates,
          "one-row pressure topology appearance advances through coupling and warm-start rebase");

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

void testTransportedTopologyAppearance() {
    Fixture fixture;
    auto settings = couplingSettings();
    settings.structure.gravityMetersPerSecondSquared = {};
    settings.convergence.minimumIterations = 1;
    settings.convergence.maximumIterations = 1;
    settings.convergence.absoluteDisplacementToleranceMetres = 1.0e6;
    settings.convergence.relativeDisplacementTolerance = 1.0e6;
    settings.convergence.absoluteVelocityToleranceMetersPerSecond = 1.0e6;
    settings.convergence.relativeVelocityTolerance = 1.0e6;
    settings.convergence.absoluteTractionToleranceNewtons = 1.0e6;
    settings.convergence.relativeTractionTolerance = 1.0e6;
    const fluid::PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    SceneFluidPressureCoupling coupling(
        fixture.surface.definition, fixture.assembly.mappings,
        fixture.structure, grid, settings);
    fluid::MacVelocityField velocity(grid);
    const auto first = coupling.advance(fixture.structure, velocity);
    check(first.accepted
              && coupling.acceptedPressureEpoch()
                     .pressureControlVolumes.controlVolumes.size() == 65,
          "transported topology fixture first accepts before its cell crossing");
    if (!first.accepted || coupling.acceptedPressureProjection() == nullptr) {
        return;
    }
    const auto& previousEpoch = coupling.acceptedPressureEpoch();
    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid, previousEpoch.pressureControlVolumes,
        previousEpoch.pressureFaceLinks, previousEpoch.openingPatches,
        *coupling.acceptedPressureProjection(), velocity);
    SceneFluidRegionTransportSettings transportSettings;
    transportSettings.timeStepSeconds =
        settings.structure.timeStepSeconds;
    const auto transport = advanceSceneFluidRegionMomentum(
        momentum, previousEpoch.pressureFaceLinks,
        *coupling.acceptedPressureProjection(), transportSettings);
    check(transport.diagnostics.accepted,
          "transported topology fixture advances its accepted region momentum");
    const std::size_t apex = apexNode(fixture);
    const double apexMass = fixture.structure.definition().nodes[apex].massKg;
    fixture.structure.addExternalForce(
        apex, {-100.0 * apexMass, 0.0, 0.0});
    const auto second = coupling.advance(
        fixture.structure, velocity, transport);
    check(second.accepted
              && second.usesRegionRebase
              && second.regionRebase.appearedControlVolumeCount == 1
              && second.usesRegionWall
              && coupling.acceptedPressureEpoch()
                     .pressureControlVolumes.controlVolumes.size() == 66
              && coupling.acceptedPressureProjection() != nullptr,
          "transported coupling accepts a one-row appearance through region rebase and wall exchange");
    if (!second.accepted || coupling.acceptedPressureProjection() == nullptr) {
        return;
    }
    const auto& crossingEpoch = coupling.acceptedPressureEpoch();
    const auto crossingMomentum = reconstructSceneFluidRegionMomentumState(
        grid, crossingEpoch.pressureControlVolumes,
        crossingEpoch.pressureFaceLinks, crossingEpoch.openingPatches,
        *coupling.acceptedPressureProjection(), velocity);
    const auto crossingTransport = advanceSceneFluidRegionMomentum(
        crossingMomentum, crossingEpoch.pressureFaceLinks,
        *coupling.acceptedPressureProjection(), transportSettings);
    check(crossingTransport.diagnostics.accepted,
          "transported retirement fixture advances its crossing momentum");
    fixture.structure.addExternalForce(
        apex, {200.0 * apexMass, 0.0, 0.0});
    const auto third = coupling.advance(
        fixture.structure, velocity, crossingTransport);
    check(third.accepted
              && third.usesRegionRebase
              && third.regionRebase.appearedControlVolumeCount == 0
              && third.regionRebase.disappearedControlVolumeCount == 1
              && third.regionRebase
                     .sourceMomentumMappingResidualNormKilogramMetersPerSecond
                  < 1.0e-12
              && third.usesRegionWall
              && coupling.acceptedPressureEpoch()
                     .pressureControlVolumes.controlVolumes.size() == 65,
          "transported coupling conservatively retires a disappearing row through wall exchange and pressure projection");
}

void testCheckpointReplayAndTransactionalRejection() {
    Fixture fixture;
    const auto settings = couplingSettings();
    const fluid::PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    SceneFluidPressureCoupling coupling(
        fixture.surface.definition, fixture.assembly.mappings,
        fixture.structure, grid, settings);
    fluid::MacVelocityField velocity(grid);

    const auto initial = coupling.checkpoint(fixture.structure);
    const auto expectedFirst = coupling.advance(fixture.structure, velocity);
    const auto expectedFirstStructure = fixture.structure.checkpoint();
    coupling.restore(fixture.structure, initial);
    check(fixture.structure.acceptedStepCount() == 0
              && coupling.acceptedPressureProjection() == nullptr,
          "initial pressure-coupling checkpoint restores the zero-pressure baseline");
    const auto replayFirst = coupling.advance(fixture.structure, velocity);
    check(replayFirst == expectedFirst
              && samePublicCheckpoint(
                  fixture.structure.checkpoint(), expectedFirstStructure),
          "initial pressure-coupling checkpoint reproduces the exact first macro-step");

    const auto accepted = coupling.checkpoint(fixture.structure);
    const auto acceptedMac =
        coupling.acceptedPressureCorrectedMacVelocity();
    const auto expectedNext = coupling.advance(
        fixture.structure, acceptedMac.velocityMetersPerSecond);
    const auto expectedNextStructure = fixture.structure.checkpoint();
    coupling.restore(fixture.structure, accepted);
    const auto replayMac =
        coupling.acceptedPressureCorrectedMacVelocity();
    const auto replayNext = coupling.advance(
        fixture.structure, replayMac.velocityMetersPerSecond);
    check(replayNext == expectedNext
              && replayMac == acceptedMac
              && samePublicCheckpoint(
                  fixture.structure.checkpoint(), expectedNextStructure),
          "accepted pressure-coupling checkpoint reproduces the exact next macro-step");

    Fixture equivalentFixture;
    SceneFluidPressureCoupling equivalent(
        equivalentFixture.surface.definition,
        equivalentFixture.assembly.mappings,
        equivalentFixture.structure, grid, settings);
    equivalent.restore(equivalentFixture.structure, accepted);
    const auto equivalentMac =
        equivalent.acceptedPressureCorrectedMacVelocity();
    const auto equivalentNext = equivalent.advance(
        equivalentFixture.structure,
        equivalentMac.velocityMetersPerSecond);
    check(equivalentNext == expectedNext
              && equivalentMac == acceptedMac
              && samePublicCheckpoint(
                  equivalentFixture.structure.checkpoint(),
                  expectedNextStructure),
          "an equivalent owner restores accepted pressure state and replays exactly");

    const auto preservedStructure = fixture.structure.checkpoint();
    const std::uint64_t preservedEpoch =
        coupling.acceptedPressureEpoch().fingerprint;
    const std::uint64_t preservedProjection =
        coupling.acceptedPressureProjection()->fingerprint;
    const auto rejectedPreservesState = [&] {
        return samePublicCheckpoint(
                   fixture.structure.checkpoint(), preservedStructure)
            && coupling.acceptedPressureEpoch().fingerprint
                == preservedEpoch
            && coupling.acceptedPressureProjection() != nullptr
            && coupling.acceptedPressureProjection()->fingerprint
                == preservedProjection;
    };

    auto corrupt = accepted;
    ++corrupt.version;
    bool rejected = false;
    try {
        coupling.restore(fixture.structure, corrupt);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && rejectedPreservesState(),
          "checkpoint version rejection preserves current coupled state");

    corrupt = accepted;
    corrupt.settings.structure.gravityMetersPerSecondSquared.x -= 1.0;
    rejected = false;
    try {
        coupling.restore(fixture.structure, corrupt);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && rejectedPreservesState(),
          "foreign solver settings reject without changing coupled state");

    corrupt = accepted;
    corrupt.pressureProjection->pressurePascals.front() += 1.0;
    rejected = false;
    try {
        coupling.restore(fixture.structure, corrupt);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && rejectedPreservesState(),
          "corrupt accepted pressure payload rejects transactionally");

    corrupt = accepted;
    corrupt.pressureProjection.reset();
    rejected = false;
    try {
        coupling.restore(fixture.structure, corrupt);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && rejectedPreservesState(),
          "noninitial checkpoint cannot silently discard accepted pressure");
}

} // namespace

int main() {
    try {
        testStrongPressureFeedbackConvergesDeterministically();
        testExhaustionAndProjectionFailureRollback();
        testCouplingInterfaceLimits();
        testTransportedTopologyAppearance();
        testCheckpointReplayAndTransactionalRejection();
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
