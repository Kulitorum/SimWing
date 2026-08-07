#include "structure.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>

namespace {

using simwing::fsi::Structure;
using simwing::fsi::StructureConstraintDefinition;
using simwing::fsi::StructureConstraintKind;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureNodeDefinition;
using simwing::fsi::StructureStepSettings;
using simwing::fsi::StructureVector3;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(double actual,
               double expected,
               double tolerance,
               const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message,
                     actual,
                     expected);
        ++failures;
    }
}

StructureStepSettings unconstrainedSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.25;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {0.0, 0.0, -2.0};
    settings.velocityDampingPerSecond = 0.0;
    return settings;
}

void testGravityAndNodalForceResponse() {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 0.0}, 2.0, false});
    Structure structure(std::move(definition));
    structure.addExternalForce(0, {4.0, 0.0, 0.0});

    const auto diagnostics = structure.step(unconstrainedSettings());
    const auto state = structure.nodeStates();
    checkNear(state[0].velocityMetersPerSecond.x, 0.5, 1.0e-15,
              "force: F/m produces the expected x velocity");
    checkNear(state[0].velocityMetersPerSecond.z, -0.5, 1.0e-15,
              "gravity: acceleration produces the expected z velocity");
    checkNear(state[0].positionMeters.x, 0.125, 1.0e-15,
              "force: semi-implicit Euler produces the expected x position");
    checkNear(state[0].positionMeters.z, -0.125, 1.0e-15,
              "gravity: semi-implicit Euler produces the expected z position");
    checkNear(diagnostics.linearMomentumKgMetersPerSecond.x, 1.0, 1.0e-15,
              "diagnostics: linear momentum is defined from node state");
    checkNear(diagnostics.lastAppliedExternalForceNewtons.x, 4.0, 0.0,
              "diagnostics: accepted force is reported");
    checkNear(diagnostics.pendingExternalForceNewtons.x, 0.0, 0.0,
              "force: accepted step consumes pending loads");
    check(diagnostics.finite, "force: diagnostics remain finite");
}

Structure constrainedPair(StructureConstraintKind kind,
                          double initialLength) {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 0.0}, 0.0, true});
    definition.nodes.push_back({{initialLength, 0.0, 0.0}, 1.0, false});
    definition.constraints.push_back(
        {kind, 0, 1, 1.0, 0.0});
    return Structure(std::move(definition));
}

StructureStepSettings constraintSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 60.0;
    settings.substeps = 1;
    settings.constraintIterations = 1;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    return settings;
}

void testDistanceCableAndSuspensionTie() {
    Structure distance = constrainedPair(
        StructureConstraintKind::Distance, 2.0);
    auto diagnostics = distance.step(constraintSettings());
    checkNear(distance.nodeStates()[1].positionMeters.x, 1.0, 1.0e-12,
              "distance: hard bilateral constraint restores its length");
    checkNear(diagnostics.maximumDistanceErrorMeters, 0.0, 1.0e-12,
              "distance: residual diagnostic matches the restored pair");

    Structure slackCable = constrainedPair(
        StructureConstraintKind::Cable, 0.5);
    static_cast<void>(slackCable.step(constraintSettings()));
    checkNear(slackCable.nodeStates()[1].positionMeters.x, 0.5, 0.0,
              "cable: unilateral constraint does not push a slack line");

    Structure tautCable = constrainedPair(
        StructureConstraintKind::Cable, 2.0);
    diagnostics = tautCable.step(constraintSettings());
    checkNear(tautCable.nodeStates()[1].positionMeters.x, 1.0, 1.0e-12,
              "cable: unilateral constraint limits extension");
    checkNear(diagnostics.maximumCableExtensionMeters, 0.0, 1.0e-12,
              "cable: extension diagnostic is unilateral");

    Structure suspensionTie = constrainedPair(
        StructureConstraintKind::SuspensionTie, 2.0);
    static_cast<void>(suspensionTie.step(constraintSettings()));
    checkNear(suspensionTie.nodeStates()[1].positionMeters.x, 1.0, 1.0e-12,
              "suspension: generic tie uses the trusted bilateral primitive");
}

void testCheckpointRestoreAndDeterminism() {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 0.0}, 1.0, false});
    Structure structure(definition);

    StructureStepSettings settings = unconstrainedSettings();
    settings.gravityMetersPerSecondSquared = {};
    structure.addExternalForce(0, {2.0, -1.0, 0.5});
    static_cast<void>(structure.step(settings));
    structure.addExternalForce(0, {-0.25, 3.0, 0.0});
    const auto saved = structure.checkpoint();

    const auto firstDiagnostics = structure.step(settings);
    const auto firstState = structure.nodeStates();
    structure.restore(saved);
    const auto restoredDiagnostics = structure.step(settings);
    const auto restoredState = structure.nodeStates();

    check(firstState == restoredState,
          "checkpoint: restore gives bit-identical continuation state");
    check(firstDiagnostics == restoredDiagnostics,
          "checkpoint: restore gives deterministic diagnostics");
    check(structure.acceptedStepCount() == saved.acceptedStepCount + 1,
          "checkpoint: accepted step count rolls back and advances once");
    checkNear(structure.simulationTimeSeconds(),
              saved.simulationTimeSeconds + settings.timeStepSeconds,
              0.0,
              "checkpoint: simulation time rolls back and advances once");

    const auto committed = structure.checkpoint();
    auto tampered = committed;
    tampered.nodes[0].positionMeters.x += 1.0;
    bool tamperRejected = false;
    try {
        structure.restore(tampered);
    } catch (const std::invalid_argument&) {
        tamperRejected = true;
    }
    check(tamperRejected && structure.nodeStates() == committed.nodes,
          "checkpoint: public node tampering cannot diverge from the opaque complete state");

    Structure other(constrainedPair(
        StructureConstraintKind::Distance, 1.0).definition());
    bool rejected = false;
    try {
        other.restore(saved);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "checkpoint: definition fingerprint rejects a foreign structure");
}

void testFailedStepRollsBackLoadsAndState() {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 0.0}, 1.0, false});
    Structure structure(std::move(definition));
    structure.addExternalForce(0, {1.0, 2.0, 3.0});
    const auto before = structure.checkpoint();
    StructureStepSettings invalid = unconstrainedSettings();
    invalid.timeStepSeconds = 0.0;
    bool rejected = false;
    try {
        static_cast<void>(structure.step(invalid));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, "rollback: invalid core step is reported");
    const auto after = structure.checkpoint();
    check(before.acceptedStepCount == after.acceptedStepCount
              && before.simulationTimeSeconds == after.simulationTimeSeconds
              && before.nodes == after.nodes
              && before.pendingExternalForcesNewtons
                  == after.pendingExternalForcesNewtons,
          "rollback: failed step preserves committed state and pending loads");
}

void testMembraneBoundaryAndValidation() {
    StructureDefinition definition;
    definition.nodes = {
        StructureNodeDefinition{{0.0, 0.0, 0.0}, 0.0, true},
        StructureNodeDefinition{{1.0, 0.0, 0.0}, 1.0, false},
        StructureNodeDefinition{{0.0, 1.0, 0.0}, 1.0, false},
    };
    definition.triangles.push_back({{0, 1, 2}});
    simwing::fsi::StructureMembraneDefinition membrane;
    membrane.triangle = 0;
    membrane.materialCoordinates = {
        simwing::fsi::StructureVector2{0.0, 0.0},
        simwing::fsi::StructureVector2{1.0, 0.0},
        simwing::fsi::StructureVector2{0.0, 1.0},
    };
    membrane.material.warpStiffnessNewtonsPerMeter = 800.0;
    membrane.material.weftStiffnessNewtonsPerMeter = 500.0;
    membrane.material.couplingStiffnessNewtonsPerMeter = 100.0;
    membrane.material.shearStiffnessNewtonsPerMeter = 180.0;
    definition.membranes.push_back(membrane);

    Structure membranePatch(std::move(definition));
    StructureStepSettings settings = constraintSettings();
    settings.constraintIterations = 4;
    const auto diagnostics = membranePatch.step(settings);
    check(diagnostics.membraneCount == 1 && diagnostics.finite,
          "membrane: orthotropic primitive is wrapped and finite");
    checkNear(diagnostics.maximumAbsoluteMembraneStrain, 0.0, 1.0e-15,
              "membrane: rest chart reports zero strain");

    StructureDefinition invalid;
    invalid.nodes.push_back({{0.0, 0.0, 0.0}, 0.0, false});
    bool rejected = false;
    try {
        Structure shouldFail(std::move(invalid));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, "validation: dynamic zero-mass node is rejected");
}

void testContactCheckpointReplay() {
    StructureDefinition definition;
    definition.nodes = {
        {{-1.0, -1.0, 0.0}, 1.0, false},
        {{1.0, -1.0, 0.0}, 1.0, false},
        {{0.0, 1.0, 0.0}, 1.0, false},
        {{-1.0, -1.0, 0.01}, 1.0, false},
        {{0.0, 1.0, 0.01}, 1.0, false},
        {{1.0, -1.0, 0.01}, 1.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{3, 4, 5}}};
    definition.fabricSelfContact = {
        0.01, 0.0, 0.0, 0.0};
    Structure structure(std::move(definition));
    StructureStepSettings settings = constraintSettings();
    settings.constraintIterations = 4;
    const auto contactDiagnostics = structure.step(settings);
    const auto saved = structure.checkpoint();
    const auto firstDiagnostics = structure.step(settings);
    const auto firstState = structure.nodeStates();
    structure.restore(saved);
    const auto replayDiagnostics = structure.step(settings);
    const auto replayState = structure.nodeStates();

    check(contactDiagnostics.contactPairCount == 1,
          "contact: explicit fabric self-pair is registered");
    check(contactDiagnostics.activeContactCount > 0,
          "contact: overlapping disconnected sheets are active");
    check(firstState == replayState,
          "contact checkpoint: warm-start continuation is bit-identical");
    check(firstDiagnostics == replayDiagnostics,
          "contact checkpoint: diagnostics replay bit-identically");
}

void testSuspensionPilotCheckpointReplay() {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 1.0}, 0.0, true});
    simwing::fsi::StructureSuspensionDefinition suspension;
    suspension.pilotStableId = 100;
    suspension.pilotMassKg = 2.0;
    suspension.pilotInitialCenterOfMassWorldMeters = {};
    suspension.pilotPrincipalInertiaKgSquareMeters = {0.2, 0.3, 0.4};
    suspension.attachments.push_back({200, 0});
    suspension.harnessPoints.push_back({300, {}});
    suspension.segments.push_back(
        {400,
         {simwing::fsi::StructureSuspensionEndpointKind::SurfaceAttachment,
          200},
         {simwing::fsi::StructureSuspensionEndpointKind::PilotHarness, 300},
         0.75, 10000.0, 0.0, 1});
    definition.suspension = suspension;
    Structure structure(std::move(definition));
    StructureStepSettings settings = constraintSettings();
    settings.constraintIterations = suspension.solverIterations;
    static_cast<void>(structure.step(settings));
    const auto saved = structure.checkpoint();
    const auto firstDiagnostics = structure.step(settings);
    const auto firstNodes = structure.nodeStates();
    const auto firstSuspension = structure.suspensionState();
    structure.restore(saved);
    const auto replayDiagnostics = structure.step(settings);
    const auto replayNodes = structure.nodeStates();
    const auto replaySuspension = structure.suspensionState();

    check(firstDiagnostics.suspensionSegmentCount == 1,
          "suspension: composite adapter registers the line graph");
    check(firstDiagnostics.totalDynamicMassKg == suspension.pilotMassKg,
          "suspension: rigid pilot contributes to total mass diagnostics");
    check(firstSuspension.has_value()
              && firstSuspension->harnessPositionsMeters.size() == 1
              && firstSuspension->segments.size() == 1,
          "suspension: immutable public payload/line state is available");
    check(firstNodes == replayNodes && firstSuspension == replaySuspension,
          "suspension checkpoint: payload and line replay is bit-identical");
    check(firstDiagnostics == replayDiagnostics,
          "suspension checkpoint: composite diagnostics replay exactly");
}

} // namespace

int main() {
    testGravityAndNodalForceResponse();
    testDistanceCableAndSuspensionTie();
    testCheckpointRestoreAndDeterminism();
    testFailedStepRollsBackLoadsAndState();
    testMembraneBoundaryAndValidation();
    testContactCheckpointReplay();
    testSuspensionPilotCheckpointReplay();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing structure check(s) failed\n", failures);
        return 1;
    }
    std::puts("all SimWing structure checks passed");
    return 0;
}
