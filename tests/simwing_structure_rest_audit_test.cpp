#include "structure_rest_audit.h"

#include <cstdio>
#include <stdexcept>
#include <utility>

namespace {

using simwing::fsi::Structure;
using simwing::fsi::StructureConstraintKind;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureStepSettings;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

Structure makePair(const double initialLength) {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 0.0}, 0.0, true});
    definition.nodes.push_back({{initialLength, 0.0, 0.0}, 1.0, false});
    definition.constraints.push_back({
        StructureConstraintKind::Distance, 0, 1, 1.0, 0.0,
    });
    return Structure(std::move(definition));
}

Structure makeExtendedSuspension() {
    StructureDefinition definition;
    definition.nodes.push_back({{0.0, 0.0, 1.0}, 0.0, true});
    simwing::fsi::StructureSuspensionDefinition suspension;
    suspension.pilotStableId = 100;
    suspension.pilotMassKg = 2.0;
    suspension.pilotPrincipalInertiaKgSquareMeters = {0.2, 0.3, 0.4};
    suspension.solverIterations = 16;
    suspension.attachments.push_back({200, 0});
    suspension.harnessPoints.push_back({300, {}});
    suspension.segments.push_back({
        400,
        {simwing::fsi::StructureSuspensionEndpointKind::SurfaceAttachment,
         200},
        {simwing::fsi::StructureSuspensionEndpointKind::PilotHarness, 300},
        0.75, 10000.0, 0.0, 1,
    });
    definition.suspension = std::move(suspension);
    return Structure(std::move(definition));
}

StructureStepSettings auditStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 60.0;
    settings.substeps = 8;
    settings.constraintIterations = 16;
    settings.gravityMetersPerSecondSquared = {};
    return settings;
}

void testStationaryAndPrestrainedStates() {
    auto stationary = makePair(1.0);
    const auto stationaryBefore = stationary.checkpoint();
    const auto first = simwing::fsi::auditStructureRestState(
        stationary, auditStepSettings());
    const auto second = simwing::fsi::auditStructureRestState(
        stationary, auditStepSettings());
    const auto stationaryAfter = stationary.checkpoint();
    check(first.stationary && second.stationary
              && first.maximumNodeDisplacementMeters == 0.0
              && first.maximumNodeSpeedMetersPerSecond == 0.0
              && first.maximumSuspensionResidualMeters == 0.0
              && first.maximumNodeDisplacementMeters
                  == second.maximumNodeDisplacementMeters,
          "rest audit deterministically accepts an exact rest state");
    check(stationaryBefore.acceptedStepCount
                  == stationaryAfter.acceptedStepCount
              && stationaryBefore.simulationTimeSeconds
                  == stationaryAfter.simulationTimeSeconds
              && stationaryBefore.nodes == stationaryAfter.nodes
              && stationaryBefore.pendingExternalForcesNewtons
                  == stationaryAfter.pendingExternalForcesNewtons
              && stationaryBefore.lastAppliedExternalForceNewtons
                  == stationaryAfter.lastAppliedExternalForceNewtons,
          "rest audit restores the exact public Structure checkpoint");

    auto prestrained = makePair(2.0);
    const auto prestrainedBefore = prestrained.checkpoint();
    const auto displaced = simwing::fsi::auditStructureRestState(
        prestrained, auditStepSettings());
    const auto prestrainedAfter = prestrained.checkpoint();
    check(!displaced.stationary
              && displaced.maximumDisplacementNode == 1
              && displaced.maximumNodeDisplacementMeters > 0.9,
          "rest audit diagnoses a nonstationary prestrained state");
    check(prestrainedBefore.acceptedStepCount
                  == prestrainedAfter.acceptedStepCount
              && prestrainedBefore.simulationTimeSeconds
                  == prestrainedAfter.simulationTimeSeconds
              && prestrainedBefore.nodes == prestrainedAfter.nodes
              && prestrainedBefore.pendingExternalForcesNewtons
                  == prestrainedAfter.pendingExternalForcesNewtons,
          "nonstationary rest audit remains transactional");
}

void testInvalidProbeInputsAreRejected() {
    auto loaded = makePair(1.0);
    loaded.addExternalForce(1, {1.0, 0.0, 0.0});
    const auto before = loaded.checkpoint();
    bool loadedRejected = false;
    try {
        static_cast<void>(simwing::fsi::auditStructureRestState(
            loaded, auditStepSettings()));
    } catch (const std::invalid_argument&) {
        loadedRejected = true;
    }
    check(loadedRejected
              && loaded.checkpoint().pendingExternalForcesNewtons
                  == before.pendingExternalForcesNewtons,
          "rest audit rejects and preserves pending loads");

    auto invalidSettings = auditStepSettings();
    invalidSettings.gravityMetersPerSecondSquared.z = -9.80665;
    bool gravityRejected = false;
    try {
        static_cast<void>(simwing::fsi::auditStructureRestState(
            loaded, invalidSettings));
    } catch (const std::invalid_argument&) {
        gravityRejected = true;
    }
    check(gravityRejected,
          "rest audit rejects a nonzero gravity probe");
}

void testInitialSuspensionExtensionIsIdentified() {
    auto structure = makeExtendedSuspension();
    const auto before = structure.checkpoint();
    const auto audit = simwing::fsi::auditStructureRestState(
        structure, auditStepSettings());
    check(!audit.stationary
              && audit.maximumInitialSuspensionExtensionSegmentStableId
                  == 400
              && audit.maximumInitialSuspensionExtensionMeters == 0.25
              && structure.checkpoint().nodes == before.nodes,
          "rest audit identifies the initially extended suspension segment transactionally");
}

} // namespace

int main() {
    try {
        testStationaryAndPrestrainedStates();
        testInvalidProbeInputsAreRejected();
        testInitialSuspensionExtensionIsIdentified();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "FAIL: %s\n", exception.what());
        ++failures;
    }
    if (failures == 0) {
        std::printf("simwing Structure rest audit tests passed\n");
    }
    return failures == 0 ? 0 : 1;
}
