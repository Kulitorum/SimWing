#include "structure.h"
#include "structure_checkpoint_persistence.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

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

std::vector<std::uint8_t> serializedCheckpoint(
    const Structure& owner,
    const simwing::fsi::StructureCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::StructureCheckpointPersistenceError error;
    check(simwing::fsi::serializeStructureCheckpoint(
              owner, checkpoint, bytes, &error),
          "persistent Structure checkpoint serializes");
    return bytes;
}

std::uint64_t readU64(const std::vector<std::uint8_t>& bytes,
                      const std::size_t offset) {
    std::uint64_t result = 0;
    for (std::size_t index = 0; index < 8; ++index) {
        result |= static_cast<std::uint64_t>(bytes[offset + index])
            << (8U * index);
    }
    return result;
}

void writeU64(std::vector<std::uint8_t>& bytes,
              const std::size_t offset,
              std::uint64_t value) {
    for (std::size_t index = 0; index < 8; ++index) {
        bytes[offset + index] = static_cast<std::uint8_t>(value & 0xffU);
        value >>= 8U;
    }
}

std::uint64_t wireChecksum(
    const std::span<const std::uint8_t> bytes) noexcept {
    std::uint64_t result = 14695981039346656037ULL;
    for (const std::uint8_t value : bytes) {
        result ^= value;
        result *= 1099511628211ULL;
    }
    return result;
}

void refreshChecksum(std::vector<std::uint8_t>& bytes) {
    const std::size_t size = static_cast<std::size_t>(readU64(bytes, 8));
    writeU64(bytes, 16,
             wireChecksum(std::span<const std::uint8_t>{bytes}
                              .subspan(24, size)));
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

StructureDefinition suspensionDefinition() {
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
    return definition;
}

void testSuspensionPilotCheckpointReplay() {
    StructureDefinition definition = suspensionDefinition();
    const auto suspension = *definition.suspension;
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

void testPersistentCheckpointComposition() {
    const StructureDefinition definition = suspensionDefinition();
    Structure source(definition);
    StructureStepSettings settings = constraintSettings();
    settings.constraintIterations = definition.suspension->solverIterations;
    static_cast<void>(source.step(settings));
    source.addExternalForce(0, {0.25, -0.5, 0.75});
    const auto saved = source.checkpoint();
    const auto bytes = serializedCheckpoint(source, saved);
    check(serializedCheckpoint(source, saved) == bytes,
          "persistent Structure encoding is byte deterministic");

    Structure rebuilt(definition);
    simwing::fsi::StructureCheckpoint decoded;
    simwing::fsi::StructureCheckpointPersistenceError error;
    check(simwing::fsi::deserializeStructureCheckpoint(
              bytes, rebuilt, decoded, &error),
          "persistent Structure checkpoint decodes on rebuilt topology");
    check(serializedCheckpoint(rebuilt, decoded) == bytes,
          "decoded Structure checkpoint re-encodes byte identically");
    rebuilt.restore(decoded);
    const auto expectedDiagnostics = source.step(settings);
    const auto replayDiagnostics = rebuilt.step(settings);
    check(expectedDiagnostics == replayDiagnostics
              && source.nodeStates() == rebuilt.nodeStates()
              && source.suspensionState() == rebuilt.suspensionState()
              && serializedCheckpoint(source, source.checkpoint())
                  == serializedCheckpoint(rebuilt, rebuilt.checkpoint()),
          "decoded composite body/suspension continuation is bit-identical");

    simwing::fsi::StructureCheckpoint output = rebuilt.checkpoint();
    const auto preserved = serializedCheckpoint(rebuilt, output);
    const auto reject = [&](std::vector<std::uint8_t> corrupt,
                            const simwing::fsi::
                                StructureCheckpointPersistenceErrorCode expected,
                            const char* message,
                            const simwing::fsi::StructureCheckpointPersistenceLimits&
                                limits = {}) {
        check(!simwing::fsi::deserializeStructureCheckpoint(
                  corrupt, rebuilt, output, &error, limits)
                  && error.code == expected
                  && serializedCheckpoint(rebuilt, output) == preserved,
              message);
    };

    auto corrupt = bytes;
    corrupt[0] ^= 0xffU;
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::InvalidMagic,
           "persistent Structure rejects bad magic transactionally");
    corrupt = bytes;
    corrupt[4] = 2;
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::
               UnsupportedVersion,
           "persistent Structure rejects unsupported protocol versions");
    corrupt = bytes;
    corrupt[6] = 1;
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::InvalidData,
           "persistent Structure rejects reserved envelope bits");
    corrupt = bytes;
    corrupt.pop_back();
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::Truncated,
           "persistent Structure rejects truncation transactionally");
    corrupt = bytes;
    corrupt.push_back(0);
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::TrailingData,
           "persistent Structure rejects trailing data transactionally");
    corrupt = bytes;
    corrupt.back() ^= 1U;
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::
               ChecksumMismatch,
           "persistent Structure detects payload corruption");

    // The public-node copy starts after the 40-byte state header and its
    // count. Changing it without changing the nested body must be rejected by
    // the rebuilt-owner validation pass even with a refreshed wire checksum.
    corrupt = bytes;
    corrupt[72] ^= 1U;
    refreshChecksum(corrupt);
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::InvalidData,
           "persistent Structure rejects public/body state disagreement");

    corrupt = bytes;
    writeU64(corrupt, 64, 5'000'001);
    refreshChecksum(corrupt);
    reject(corrupt,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::LimitExceeded,
           "persistent Structure bounds node counts before allocation");

    simwing::fsi::StructureCheckpointPersistenceLimits limits;
    limits.maximumEncodedBytes = bytes.size() - 1;
    reject(bytes,
           simwing::fsi::StructureCheckpointPersistenceErrorCode::LimitExceeded,
           "persistent Structure enforces the byte limit", limits);

    Structure foreign(constrainedPair(
        StructureConstraintKind::Distance, 1.0).definition());
    check(!simwing::fsi::deserializeStructureCheckpoint(
              bytes, foreign, output, &error)
              && error.code
                  == simwing::fsi::StructureCheckpointPersistenceErrorCode::
                      TopologyMismatch
              && serializedCheckpoint(rebuilt, output) == preserved,
          "persistent Structure rejects a foreign definition transactionally");

    std::vector<std::uint8_t> callerBytes{1, 2, 3};
    limits = {};
    limits.maximumNodes = 0;
    check(!simwing::fsi::serializeStructureCheckpoint(
              source, saved, callerBytes, &error, limits)
              && callerBytes == std::vector<std::uint8_t>({1, 2, 3}),
          "failed Structure serialization preserves caller output");
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
    testPersistentCheckpointComposition();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing structure check(s) failed\n", failures);
        return 1;
    }
    std::puts("all SimWing structure checks passed");
    return 0;
}
