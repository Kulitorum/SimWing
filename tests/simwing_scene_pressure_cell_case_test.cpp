#include "scene_pressure_cell_case.h"
#include "scene_pressure_cell_checkpoint_persistence.h"
#include "scene_pressure_cell_operator_phase_audit.h"
#include "scene_pressure_cell_operator_phase_refinement_audit.h"
#include "scene_pressure_cell_operator_refinement_audit.h"
#include "scene_pressure_cell_mimetic_conductance_phase_refinement_audit.h"
#include "scene_fluid_mimetic_region_conductance_audit.h"
#include "scene_fluid_pressure_operator_response_audit.h"
#include "scene_fluid_pressure_owner_transition.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <ranges>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

fsi::SceneFluidPressureOperatorResponseAuditSettings
pressureCellResponseAuditSettings() {
    fsi::SceneFluidPressureOperatorResponseAuditSettings settings;
    settings.graphSolve.absoluteResidualTolerancePascalsMeters = 1.0e-10;
    settings.graphSolve.relativeResidualTolerance = 1.0e-11;
    settings.graphSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    settings.shadowSolve.absoluteResidualTolerancePascalsMeters = 1.0e-10;
    settings.shadowSolve.relativeResidualTolerance = 1.0e-11;
    settings.shadowSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    return settings;
}

const std::vector<fsi::fluid::Vector3>& pressureCellGridPhases() {
    static const std::vector<fsi::fluid::Vector3> phases{
        {0.0, 0.0, 0.0},
        {-0.5, 0.0, 0.0},
        {0.0, -0.5, 0.0},
        {0.0, 0.0, -0.5},
        {-0.5, -0.5, 0.0},
        {-0.5, 0.0, -0.5},
        {0.0, -0.5, -0.5},
        {-0.5, -0.5, -0.5},
    };
    return phases;
}

double norm(const fsi::StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "scene pressure cell frame serializes");
    return bytes;
}

std::vector<std::uint8_t> serializedCheckpoint(
    const fsi::ScenePressureCellCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    fsi::ScenePressureCellCheckpointPersistenceError error;
    check(fsi::serializeScenePressureCellCheckpoint(
              checkpoint, bytes, &error),
          "scene pressure cell checkpoint serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

const viewer::VectorField* vectorField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.vectorFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.vectorFields.end() ? nullptr : &*found;
}

void testVisibleStrongPressureCellAndReplay() {
    fsi::ScenePressureCellCase first;
    fsi::ScenePressureCellCase second;
    viewer::DiagnosticFrame frame;
    double peakDisplacement = 0.0;
    double peakPressure = 0.0;
    double peakPressureForce = 0.0;
    double peakWallForce = 0.0;
    double peakWallLoss = 0.0;
    double peakWallMomentumResidual = 0.0;
    double peakMacSpeed = 0.0;
    double peakMacSubfaceDeviation = 0.0;
    double peakBulkFlowChange = 0.0;
    double peakBulkViscousLoss = 0.0;
    double peakRegionEnergyLoss = 0.0;
    double peakRegionGclVolumeChange = 0.0;
    double peakRegionMomentumResidual = 0.0;
    double peakFlowPumpForce = 0.0;
    std::uint64_t peakIterations = 0;
    for (std::size_t step = 0; step < 120; ++step) {
        frame = first.advance();
        const auto repeated = second.advance();
        check(serialized(frame) == serialized(repeated),
              "scene pressure cell accepted frames are deterministic");
        const auto& diagnostics = first.diagnostics();
        peakDisplacement = std::max(
            peakDisplacement, diagnostics.maximumDisplacementMeters);
        peakPressure = std::max(
            peakPressure, diagnostics.maximumAbsolutePressurePascals);
        peakPressureForce = std::max(
            peakPressureForce, norm(diagnostics.pressureForceNewtons));
        peakWallForce = std::max(
            peakWallForce, norm(diagnostics.wallForceNewtons));
        if (diagnostics.coupling.usesRegionWall) {
            peakWallLoss = std::max(
                peakWallLoss,
                diagnostics.coupling.regionWall.viscousDissipationJoules);
            peakWallMomentumResidual = std::max(
                peakWallMomentumResidual,
                diagnostics.coupling.regionWall
                    .momentumResidualNormKilogramMetersPerSecond);
        }
        peakMacSpeed = std::max(
            peakMacSpeed, diagnostics.macVelocity
                .maximumAbsoluteVelocityMetersPerSecond);
        peakMacSubfaceDeviation = std::max(
            peakMacSubfaceDeviation, diagnostics.macVelocity
                .maximumSubfaceVelocityDeviationMetersPerSecond);
        peakBulkFlowChange = std::max(
            peakBulkFlowChange, diagnostics.bulkFlow
                .maximumVelocityChangeMetersPerSecond);
        peakBulkViscousLoss = std::max(
            peakBulkViscousLoss,
            diagnostics.bulkFlow.firstHalfViscousEnergyLossJoules
                + diagnostics.bulkFlow.secondHalfViscousEnergyLossJoules);
        if (diagnostics.usesRegionTransport) {
            peakRegionEnergyLoss = std::max(
                peakRegionEnergyLoss,
                diagnostics.regionTransport.advectiveEnergyLossJoules
                    + diagnostics.regionTransport.viscousEnergyLossJoules);
            peakRegionGclVolumeChange = std::max(
                peakRegionGclVolumeChange,
                diagnostics.regionTransport
                    .maximumAbsoluteGeometryVolumeChangeCubicMeters);
            peakRegionMomentumResidual = std::max(
                peakRegionMomentumResidual,
                diagnostics.regionTransport
                    .momentumResidualNormKilogramMetersPerSecond);
        }
        peakFlowPumpForce = std::max(
            peakFlowPumpForce, std::abs(diagnostics.flowPumpForceNewtons));
        peakIterations = std::max(
            peakIterations, diagnostics.coupling.solverRunCount);
        check(diagnostics.finite
                  && diagnostics.coupling.accepted
                  && diagnostics.coupling.iteration.status
                      == fsi::StrongCouplingIterationStatus::Converged
                  && diagnostics.coupling.pressureProjection.accepted
                  && diagnostics.coupling.pressureProjection
                         .correctedContinuityResidualMaximumCubicMetersPerSecond
                      < 2.0e-11
                  && diagnostics.coupling.pressureTransfer
                         .forceResidualNormNewtons < 1.0e-12
                  && diagnostics.coupling.pressureTransfer
                         .momentResidualNormNewtonMeters < 1.0e-12
                  && diagnostics.coupling.totalFluidTransfer
                         .forceResidualNormNewtons < 1.0e-12
                  && diagnostics.coupling.totalFluidTransfer
                         .momentResidualNormNewtonMeters < 1.0e-12,
              "scene pressure cell closes every pumped bulk-flow and strong pressure/wall-feedback step");
        check(diagnostics.usesRegionTransport == (step != 0)
                  && (!diagnostics.usesRegionTransport
                      || (diagnostics.regionTransport.accepted
                          && diagnostics.regionTransport
                                 .usesMovingVolumeRates
                          && diagnostics.regionTransport
                                 .usesBulkVelocityIncrement)),
              "scene pressure cell advances accepted moving-volume region momentum after bootstrap");
        check(diagnostics.coupling.usesRegionWall == (step != 0)
                  && (!diagnostics.coupling.usesRegionWall
                      || (diagnostics.coupling.regionWall.accepted
                          && diagnostics.coupling.regionWall.finite
                          && diagnostics.coupling.regionWall
                                 .viscousDissipationJoules >= 0.0
                          && diagnostics.coupling.regionWall
                                 .momentumResidualNormKilogramMetersPerSecond
                              < 1.0e-10)),
              "scene pressure cell exchanges conservative tangential material-wall momentum after bootstrap");
        check(diagnostics.bulkFlow.accepted
                  && diagnostics.bulkFlow.finite
                  && diagnostics.bulkFlow.finalDivergenceL2PerSecond
                      < 1.0e-11,
              "scene pressure cell viscous bulk predictor remains projected and finite");
    }

    check(peakDisplacement > 1.0e-5
              && peakPressure > 1.0e-4
              && peakPressureForce > 1.0e-6
              && peakWallForce > 1.0e-9
              && peakWallLoss > 0.0
              && peakWallMomentumResidual < 1.0e-10
              && peakMacSpeed > 1.0e-6
              && peakMacSubfaceDeviation > 1.0e-8
              && peakBulkFlowChange > 1.0e-8
              && peakBulkViscousLoss > 0.0
              && peakRegionEnergyLoss >= 0.0
              && peakRegionGclVolumeChange > 0.0
              && peakRegionMomentumResidual < 1.0e-10
              && peakFlowPumpForce > 1.0e-8
              && peakIterations >= 2,
          "visible cell develops sustained wind-driven motion, sparse pressure, conservative load, and an evolving bulk MAC predictor");
    check(frame.step == 120
              && frame.vertices.size() == 4
              && frame.triangles.size() == 3
              && frame.sceneChecksum == fsi::scenePressureCellCaseChecksum
              && frame.solverCommit == fsi::scenePressureCellCaseSolverId,
          "scene pressure cell publishes its accepted scene-v2 shell");

    const auto* displacement = scalarField(
        frame, "pressure_cell.displacement");
    const auto* pressureJump = scalarField(
        frame, "pressure_cell.pressure_jump");
    const auto* maximumPressure = scalarField(
        frame, "pressure_cell.maximum_pressure");
    const auto* iterations = scalarField(
        frame, "pressure_cell.coupling_iterations");
    const auto* macSpeed = scalarField(
        frame, "pressure_cell.maximum_mac_speed");
    const auto* macDeviation = scalarField(
        frame, "pressure_cell.mac_subface_deviation");
    const auto* macEmbeddedOpenings = scalarField(
        frame, "pressure_cell.mac_embedded_openings");
    const auto* bulkChange = scalarField(
        frame, "pressure_cell.bulk_flow_change");
    const auto* bulkDivergence = scalarField(
        frame, "pressure_cell.bulk_divergence");
    const auto* bulkViscousLoss = scalarField(
        frame, "pressure_cell.bulk_viscous_loss");
    const auto* regionEnergyLoss = scalarField(
        frame, "pressure_cell.region_transport_energy_loss");
    const auto* regionMomentumResidual = scalarField(
        frame, "pressure_cell.region_transport_momentum_residual");
    const auto* regionGclVolumeChange = scalarField(
        frame, "pressure_cell.region_gcl_volume_change");
    const auto* wallLoss = scalarField(
        frame, "pressure_cell.wall_viscous_loss");
    const auto* wallMomentumResidual = scalarField(
        frame, "pressure_cell.wall_momentum_residual");
    const auto* nodalForce = vectorField(
        frame, "pressure_cell.nodal_pressure_force");
    const auto* totalForce = vectorField(
        frame, "pressure_cell.total_pressure_force");
    const auto* nodalWallForce = vectorField(
        frame, "pressure_cell.nodal_wall_force");
    const auto* totalWallForce = vectorField(
        frame, "pressure_cell.total_wall_force");
    const auto* nodalTotalFluidForce = vectorField(
        frame, "pressure_cell.nodal_total_fluid_force");
    const auto* totalFluidForce = vectorField(
        frame, "pressure_cell.total_fluid_force");
    const auto* flowPump = vectorField(
        frame, "pressure_cell.flow_pump_force");
    check(displacement != nullptr
              && displacement->association
                  == viewer::FieldAssociation::Vertex
              && displacement->values.size() == frame.vertices.size()
              && pressureJump != nullptr
              && pressureJump->association
                  == viewer::FieldAssociation::Triangle
              && pressureJump->values.size() == frame.triangles.size()
              && maximumPressure != nullptr
              && maximumPressure->values.size() == 1
              && iterations != nullptr && iterations->values.size() == 1
              && macSpeed != nullptr && macSpeed->values.size() == 1
              && macDeviation != nullptr
              && macDeviation->values.size() == 1
              && macEmbeddedOpenings != nullptr
              && macEmbeddedOpenings->values.size() == 1
              && bulkChange != nullptr && bulkChange->values.size() == 1
              && bulkDivergence != nullptr
              && bulkDivergence->values.size() == 1
              && bulkViscousLoss != nullptr
              && bulkViscousLoss->values.size() == 1
              && regionEnergyLoss != nullptr
              && regionEnergyLoss->values.size() == 1
              && regionMomentumResidual != nullptr
              && regionMomentumResidual->values.size() == 1
              && regionGclVolumeChange != nullptr
              && regionGclVolumeChange->values.size() == 1
              && wallLoss != nullptr && wallLoss->values.size() == 1
              && wallMomentumResidual != nullptr
              && wallMomentumResidual->values.size() == 1
              && nodalForce != nullptr
              && nodalForce->association
                  == viewer::FieldAssociation::Vertex
              && nodalForce->values.size() == frame.vertices.size()
              && totalForce != nullptr && totalForce->values.size() == 1
              && nodalWallForce != nullptr
              && nodalWallForce->values.size() == frame.vertices.size()
              && totalWallForce != nullptr
              && totalWallForce->values.size() == 1
              && nodalTotalFluidForce != nullptr
              && nodalTotalFluidForce->values.size()
                  == frame.vertices.size()
              && totalFluidForce != nullptr
              && totalFluidForce->values.size() == 1
              && flowPump != nullptr && flowPump->values.size() == 1,
          "scene pressure cell frame exposes deformation, pressure, pump forcing, viscous bulk flow, MAC continuation, and iteration count");

    const auto checkpoint = first.checkpoint();
    check(checkpoint.regionMomentum.has_value()
              && first.acceptedRegionMomentum() != nullptr
              && checkpoint.regionMomentum
                  == *first.acceptedRegionMomentum()
              && checkpoint.coupling.pressureProjection
                     ->regionLinkFlowPredictionFingerprint != 0
              && checkpoint.coupling.pressureProjection
                     ->regionWallExchangeFingerprint != 0
              && checkpoint.coupling.wallTractions.has_value()
              && checkpoint.coupling.wallTractions->wallExchangeFingerprint
                  == checkpoint.coupling.pressureProjection
                         ->regionWallExchangeFingerprint,
          "scene pressure cell checkpoint owns its accepted transported-region and material-wall continuation state");
    const auto expected = first.advance();
    const auto expectedDiagnostics = first.diagnostics();
    first.restore(checkpoint);
    const auto replay = first.advance();
    check(serialized(replay) == serialized(expected)
              && first.diagnostics() == expectedDiagnostics,
          "scene pressure cell checkpoint reproduces the exact next frame");
}

void testPersistentCheckpointAndRejection() {
    fsi::ScenePressureCellCase initial;
    const auto initialBytes = serializedCheckpoint(initial.checkpoint());
    fsi::ScenePressureCellCheckpoint initialDecoded;
    fsi::ScenePressureCellCheckpointPersistenceError error;
    check(fsi::deserializeScenePressureCellCheckpoint(
              initialBytes, initialDecoded, &error),
          "initial scene pressure cell checkpoint decodes without pressure state");
    fsi::ScenePressureCellCase initialReplay;
    initialReplay.restore(initialDecoded);
    check(serialized(initialReplay.advance()) == serialized(initial.advance()),
          "persisted initial scene pressure cell reproduces the first frame");

    fsi::ScenePressureCellCase source;
    for (std::size_t step = 0; step < 35; ++step) {
        static_cast<void>(source.advance());
    }
    const auto saved = source.checkpoint();
    const auto bytes = serializedCheckpoint(saved);
    check(!bytes.empty() && serializedCheckpoint(saved) == bytes,
          "scene pressure cell checkpoint encoding is deterministic");

    fsi::ScenePressureCellCheckpoint decoded;
    check(fsi::deserializeScenePressureCellCheckpoint(
              bytes, decoded, &error)
              && serializedCheckpoint(decoded) == bytes
              && decoded.regionMomentum == saved.regionMomentum
              && decoded.coupling.wallTractions
                  == saved.coupling.wallTractions,
          "scene pressure cell checkpoint has a canonical bounded round trip");
    fsi::ScenePressureCellCase restored;
    restored.restore(decoded);
    check(restored.predictedVelocity() == source.predictedVelocity(),
          "persisted scene pressure cell reconstructs its accepted MAC predictor");
    const auto expected = source.advance();
    const auto replay = restored.advance();
    check(serialized(replay) == serialized(expected)
              && restored.diagnostics() == source.diagnostics(),
          "persisted scene pressure cell reproduces the exact next frame");

    auto corruptMomentum = saved;
    corruptMomentum.regionMomentum->controlVolumes.front()
        .velocityMetersPerSecond.x += 0.01;
    std::vector<std::uint8_t> unchanged{1, 2, 3};
    const auto original = unchanged;
    check(!fsi::serializeScenePressureCellCheckpoint(
              corruptMomentum, unchanged, &error)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidData
              && unchanged == original,
          "scene pressure cell rejects corrupt region momentum before publishing bytes");
    const auto restoredBeforeRejection = restored.checkpoint();
    bool restoreRejected = false;
    try {
        restored.restore(corruptMomentum);
    } catch (const std::exception&) {
        restoreRejected = true;
    }
    check(restoreRejected
              && serializedCheckpoint(restored.checkpoint())
                  == serializedCheckpoint(restoredBeforeRejection),
          "scene pressure cell rejects corrupt region momentum without mutating its accepted owners");

    auto corruptWall = saved;
    corruptWall.coupling.wallTractions->tractions.front()
        .tractionPascals.x += 0.01;
    unchanged = {1, 2, 3};
    check(!fsi::serializeScenePressureCellCheckpoint(
              corruptWall, unchanged, &error)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidData
              && unchanged == original,
          "scene pressure cell rejects corrupt accepted wall traction before publishing bytes");

    const auto preserved = serializedCheckpoint(decoded);
    const auto rejects = [&](std::vector<std::uint8_t> damaged,
                             const fsi::ScenePressureCellCheckpointPersistenceErrorCode
                                 expectedCode,
                             const char* message,
                             const fsi::ScenePressureCellCheckpointPersistenceLimits&
                                 limits = {}) {
        fsi::ScenePressureCellCheckpointPersistenceError rejection;
        check(!fsi::deserializeScenePressureCellCheckpoint(
                  damaged, decoded, &rejection, limits)
                  && rejection.code == expectedCode
                  && serializedCheckpoint(decoded) == preserved,
              message);
    };

    auto damaged = bytes;
    damaged.front() ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidMagic,
            "scene pressure cell rejects foreign magic transactionally");
    damaged = bytes;
    damaged[9] ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::UnsupportedVersion,
            "scene pressure cell rejects an unsupported wire version transactionally");
    damaged = bytes;
    damaged.back() ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::ChecksumMismatch,
            "scene pressure cell rejects payload corruption transactionally");
    damaged = bytes;
    damaged.pop_back();
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
            "scene pressure cell rejects truncation transactionally");
    damaged = bytes;
    damaged.push_back(0);
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::TrailingData,
            "scene pressure cell rejects trailing data transactionally");

    auto smallLimits =
        fsi::ScenePressureCellCheckpointPersistenceLimits{};
    smallLimits.maximumEncodedBytes = bytes.size() - 1;
    rejects(bytes,
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell rejects an oversized checkpoint transactionally",
            smallLimits);

    auto momentumLimits =
        fsi::ScenePressureCellCheckpointPersistenceLimits{};
    momentumLimits.maximumMomentumStorageBytes =
        saved.regionMomentum->ownedStorageBytes - 1;
    rejects(bytes,
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell bounds persisted region momentum before publication",
            momentumLimits);

    auto wallLimits =
        fsi::ScenePressureCellCheckpointPersistenceLimits{};
    wallLimits.maximumWallTractionStorageBytes =
        saved.coupling.wallTractions->ownedStorageBytes - 1;
    rejects(bytes,
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell bounds persisted wall traction before publication",
            wallLimits);

    auto foreign = saved;
    foreign.version += 1;
    unchanged = {1, 2, 3};
    check(!fsi::serializeScenePressureCellCheckpoint(
              foreign, unchanged, &error)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidData
              && unchanged == original,
          "scene pressure cell rejects foreign state before publishing bytes");
}

void testOptInMimeticPressureAuditIsShadowOnly() {
    fsi::ScenePressureCellCase production;
    fsi::ScenePressureCellCase audited(true);
    std::uint64_t previousAuditFingerprint = 0;
    std::uint64_t previousComparisonFingerprint = 0;
    for (std::size_t step = 0; step < 4; ++step) {
        const auto productionFrame = production.advance();
        const auto auditedFrame = audited.advance();
        const auto* endpoint = audited.acceptedMimeticPressureAudit();
        const auto* comparison =
            audited.acceptedMimeticPressureComparison();
        const auto* ownerTransition =
            audited.acceptedMimeticPressureOwnerTransition();
        const auto productionCheckpoint = production.checkpoint();
        const auto auditedCheckpoint = audited.checkpoint();
        auto auditedGraphOnly = auditedCheckpoint;
        auditedGraphOnly.coupling
            .mimeticPressureAuditSettingsFingerprint = 0;
        auditedGraphOnly.coupling.mimeticPressureState.reset();
        check(serialized(auditedFrame) == serialized(productionFrame)
                  && auditedCheckpoint.coupling
                         .mimeticPressureAuditSettingsFingerprint != 0
                  && auditedCheckpoint.coupling
                         .mimeticPressureState.has_value()
                  && serializedCheckpoint(auditedGraphOnly)
                      == serializedCheckpoint(productionCheckpoint),
              "opt-in mimetic audit leaves frames and the graph-pressure portion of its composite checkpoint byte-identical");
        check(endpoint != nullptr
                  && audited.diagnostics().coupling
                      .usesMimeticPressureAudit
                  && audited.diagnostics().coupling
                         .mimeticPressureAuditFingerprint
                      == endpoint->fingerprint
                  && endpoint->pressureEpoch.diagnostics.accepted
                  && comparison != nullptr
                  && ownerTransition != nullptr
                  && audited.diagnostics().coupling
                         .mimeticPressureComparisonFingerprint
                      == comparison->fingerprint
                  && audited.diagnostics().coupling
                         .mimeticPressureOwnerTransition
                      == *ownerTransition
                  && ownerTransition->comparisonFingerprint
                      == comparison->fingerprint
                  && ownerTransition->selectedOwner
                      == fsi::SceneFluidPressureOwner::ReferenceGraph
                  && comparison->diagnostics.finite
                  && comparison->includesSourceComparison
                  && comparison->sourceDiagnostics.finite
                  && comparison->sourceDiagnostics.geometryVolumeRate.exact
                  && comparison->sourceDiagnostics
                         .predictedNetOutwardVolumeRate.relativeDeltaL2
                      < 1.0e-14
                  && comparison->sourceDiagnostics.continuityResidual
                         .relativeDeltaL2
                      < 1.0e-14
                  && comparison->sourceDiagnostics.integratedSource
                         .relativeDeltaL2
                      < 1.0e-14
                  && comparison->diagnostics.bestFitShadowPressureScale > 2.0
                  && comparison->diagnostics.bestFitShadowPressureScale < 3.0
                  && comparison->diagnostics
                         .pressureDifferenceCosineSimilarity
                      > 1.0 - 1.0e-14
                  && comparison->diagnostics
                         .relativeBestFitPressureShapeResidualL2
                      < 1.0e-14
                  && comparison->diagnostics.bestFitShadowNodalForceScale
                      > 2.0
                  && comparison->diagnostics.bestFitShadowNodalForceScale
                      < 3.0
                  && comparison->diagnostics.nodalForceCosineSimilarity
                      > 1.0 - 1.0e-14
                  && comparison->diagnostics
                         .relativeBestFitNodalForceShapeResidualL2
                      < 1.0e-14
                  && comparison->shadowControlCellFingerprint
                      == endpoint->controlCells.fingerprint
                  && comparison->shadowPressureSourceFingerprint
                      == endpoint->pressureSources.fingerprint
                  && comparison->samples.size()
                      == endpoint->pressureEpoch
                             .acceptedPressureSamples.bindings.size()
                  && 2 * endpoint->pressureEpoch
                             .acceptedPressureSamples.bindings.size()
                      == endpoint->controlCells.materialWallHalfFaceCount,
              "opt-in mimetic audit publishes one complete endpoint with roundoff-equivalent graph and shadow sources after graph convergence");
        if (endpoint != nullptr && comparison != nullptr
            && ownerTransition != nullptr) {
            fsi::validateSceneFluidMimeticPressureAuditEndpointIntegrity(
                *endpoint);
            fsi::validateSceneFluidPressureShadowComparisonIntegrity(
                *comparison);
            fsi::validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
                *ownerTransition, *comparison);
            check(endpoint->usesConsecutiveWarmStart == (step != 0)
                      && endpoint->usesRegionWallPrediction == (step != 0)
                      && endpoint->pressureEpoch
                             .topologyTransitionFingerprint
                          == (step == 0
                                  ? 0
                                  : endpoint
                                        ->pressureTopologyTransitionFingerprint)
                      && endpoint->fingerprint != previousAuditFingerprint
                      && comparison->fingerprint
                          != previousComparisonFingerprint,
                  "mimetic audit bootstraps once then advances through the transported wall predictor and consecutive warm state");
            previousAuditFingerprint = endpoint->fingerprint;
            previousComparisonFingerprint = comparison->fingerprint;
        }
    }
    const auto& acceptedComparison =
        *audited.acceptedMimeticPressureComparison();
    const auto& acceptedTransition =
        *audited.acceptedMimeticPressureOwnerTransition();
    const auto expectedTransition =
        fsi::decideSceneFluidPressureOwnerTransition(acceptedComparison);
    const auto repeatedTransition =
        fsi::decideSceneFluidPressureOwnerTransition(acceptedComparison);
    fsi::validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
        acceptedTransition, acceptedComparison);
    check(acceptedTransition == expectedTransition
              && acceptedTransition == repeatedTransition
              && acceptedTransition.selectedOwner
                  == fsi::SceneFluidPressureOwner::ReferenceGraph
              && acceptedTransition.rejectionMask == 0xeb00ULL
              && acceptedTransition.rejectionCount == 6
              && fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                  acceptedTransition,
                  fsi::SceneFluidPressureOwnerTransitionRejection::
                      PressureDifferenceMismatch)
              && fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                  acceptedTransition,
                  fsi::SceneFluidPressureOwnerTransitionRejection::
                      PressureScaleMismatch)
              && fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                  acceptedTransition,
                  fsi::SceneFluidPressureOwnerTransitionRejection::
                      NodalForceScaleMismatch)
              && fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                  acceptedTransition,
                  fsi::SceneFluidPressureOwnerTransitionRejection::
                      NetForceMismatch)
              && !fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                  acceptedTransition,
                  fsi::SceneFluidPressureOwnerTransitionRejection::
                      IntegratedSourceMismatch),
          "pressure-owner transition gate retains graph loads for the live mimetic pressure and force mismatch while accepting roundoff-equivalent sources");
    auto corruptedTransition = acceptedTransition;
    corruptedTransition.selectedOwner =
        fsi::SceneFluidPressureOwner::ShadowMimetic;
    bool rejected = false;
    try {
        fsi::validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
            corruptedTransition, acceptedComparison);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "pressure-owner transition decision rejects selected-owner corruption");
    auto invalidPolicy = acceptedTransition.policy;
    invalidPolicy.maximumRelativeNetForceDelta =
        -std::numeric_limits<double>::epsilon();
    rejected = false;
    try {
        static_cast<void>(fsi::decideSceneFluidPressureOwnerTransition(
            acceptedComparison, invalidPolicy));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "pressure-owner transition rejects invalid policy before selection");
    fsi::SceneFluidMimeticPressureAuditConfiguration invalidConfiguration;
    invalidConfiguration.enabled = true;
    invalidConfiguration.ownerTransitionPolicy = invalidPolicy;
    rejected = false;
    try {
        fsi::ScenePressureCellCase invalid(invalidConfiguration);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "pressure coupling rejects an invalid owner-transition policy during construction");

    fsi::SceneFluidMimeticPressureAuditConfiguration
        permissiveConfiguration;
    permissiveConfiguration.enabled = true;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumRelativePressureDifferenceDeltaL2 = 1.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumAbsolutePressureScaleDeviation = 2.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumRelativePressureShapeResidualL2 = 1.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumAbsoluteNodalForceScaleDeviation = 2.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumRelativeNodalForceShapeResidualL2 = 1.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumRelativeNetForceDelta = 1.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumRelativeNetMomentDelta = 1.0;
    permissiveConfiguration.ownerTransitionPolicy
        .maximumAbsolutePowerDeltaWatts = 100.0;
    fsi::ScenePressureCellCase permissive(permissiveConfiguration);
    fsi::ScenePressureCellCase permissiveReference;
    const auto permissiveFrame = permissive.advance();
    const auto referenceFrame = permissiveReference.advance();
    const auto* permissiveTransition =
        permissive.acceptedMimeticPressureOwnerTransition();
    check(serialized(permissiveFrame) == serialized(referenceFrame)
              && permissiveTransition != nullptr
              && permissiveTransition->selectedOwner
                  == fsi::SceneFluidPressureOwner::ShadowMimetic
              && permissiveTransition->rejectionMask == 0,
          "a permissive zero-rejection owner decision remains diagnostic and leaves graph-owned production frames byte-identical");

    auto corrupted = *audited.acceptedMimeticPressureAudit();
    ++corrupted.pressureEpoch.diagnostics.pressureSolve
          .reducedTraceSolve.iterationCount;
    rejected = false;
    try {
        fsi::validateSceneFluidMimeticPressureAuditEndpointIntegrity(
            corrupted);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "mimetic pressure-audit endpoint fingerprints complete nested solve diagnostics");
    auto corruptedComparison =
        *audited.acceptedMimeticPressureComparison();
    corruptedComparison.samples.front().shadowMinusReferencePascals += 0.01;
    rejected = false;
    try {
        fsi::validateSceneFluidPressureShadowComparisonIntegrity(
            corruptedComparison);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "mimetic pressure comparison rejects per-sample delta corruption");
    corruptedComparison = *audited.acceptedMimeticPressureComparison();
    corruptedComparison.controlSources.front()
        .integratedSourceDeltaPascalsMeters += 0.01;
    rejected = false;
    try {
        fsi::validateSceneFluidPressureShadowComparisonIntegrity(
            corruptedComparison);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "mimetic pressure comparison rejects per-control source corruption");

    fsi::SceneFluidMimeticPressureAuditConfiguration limitedConfiguration;
    limitedConfiguration.enabled = true;
    limitedConfiguration.settings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    limitedConfiguration.settings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    limitedConfiguration.limits.maximumOwnedBytes = 0;
    fsi::ScenePressureCellCase limited(limitedConfiguration);
    const auto beforeRejection = serializedCheckpoint(limited.checkpoint());
    rejected = false;
    try {
        static_cast<void>(limited.advance());
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected
              && limited.acceptedMimeticPressureAudit() == nullptr
              && limited.acceptedMimeticPressureOwnerTransition() == nullptr
              && serializedCheckpoint(limited.checkpoint())
                  == beforeRejection,
          "mimetic pressure-audit limit failure rolls Structure and graph pressure back transactionally");

    fsi::SceneFluidMimeticPressureAuditConfiguration
        comparisonLimitedConfiguration;
    comparisonLimitedConfiguration.enabled = true;
    comparisonLimitedConfiguration.settings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    comparisonLimitedConfiguration.settings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    comparisonLimitedConfiguration.comparisonLimits.maximumOwnedBytes = 0;
    fsi::ScenePressureCellCase comparisonLimited(
        comparisonLimitedConfiguration);
    const auto beforeComparisonRejection = serializedCheckpoint(
        comparisonLimited.checkpoint());
    rejected = false;
    try {
        static_cast<void>(comparisonLimited.advance());
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected
              && comparisonLimited.acceptedMimeticPressureAudit() == nullptr
              && comparisonLimited.acceptedMimeticPressureComparison()
                  == nullptr
              && comparisonLimited
                     .acceptedMimeticPressureOwnerTransition()
                  == nullptr
              && serializedCheckpoint(comparisonLimited.checkpoint())
                  == beforeComparisonRejection,
          "mimetic pressure-comparison limit failure rolls Structure and every accepted owner back transactionally");

    auto sourceCountLimitedConfiguration = comparisonLimitedConfiguration;
    sourceCountLimitedConfiguration.comparisonLimits.maximumOwnedBytes =
        1024ULL * 1024ULL * 1024ULL;
    sourceCountLimitedConfiguration.comparisonLimits.maximumControlVolumes =
        0;
    fsi::ScenePressureCellCase sourceCountLimited(
        sourceCountLimitedConfiguration);
    const auto beforeSourceCountRejection = serializedCheckpoint(
        sourceCountLimited.checkpoint());
    rejected = false;
    try {
        static_cast<void>(sourceCountLimited.advance());
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected
              && sourceCountLimited.acceptedMimeticPressureAudit() == nullptr
              && sourceCountLimited.acceptedMimeticPressureComparison()
                  == nullptr
              && sourceCountLimited
                     .acceptedMimeticPressureOwnerTransition()
                  == nullptr
              && serializedCheckpoint(sourceCountLimited.checkpoint())
                  == beforeSourceCountRejection,
          "mimetic pressure source-comparison count failure rolls Structure and every accepted owner back transactionally");
}

void testManufacturedPressureOperatorResponses() {
    fsi::ScenePressureCellCase simulation(true);
    static_cast<void>(simulation.advance());
    const auto* endpoint = simulation.acceptedMimeticPressureAudit();
    check(endpoint != nullptr,
          "manufactured operator audit has an accepted mimetic topology");
    if (endpoint == nullptr) return;
    fsi::SceneFluidPressureOperatorResponseAuditSettings settings;
    settings.graphSolve.absoluteResidualTolerancePascalsMeters = 1.0e-10;
    settings.graphSolve.relativeResidualTolerance = 1.0e-11;
    settings.graphSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    settings.shadowSolve.absoluteResidualTolerancePascalsMeters = 1.0e-10;
    settings.shadowSolve.relativeResidualTolerance = 1.0e-11;
    settings.shadowSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    const auto acceptedSource =
        fsi::sceneFluidMimeticIntegratedCellSources(
            endpoint->pressureSources);
    const auto audit = fsi::auditSceneFluidPressureOperatorResponses(
        simulation.acceptedPressureEpoch().pressureOperator,
        endpoint->controlCells, endpoint->fullTraceSystem,
        endpoint->condensedTraceSystem, acceptedSource, settings);
    const auto repeated = fsi::auditSceneFluidPressureOperatorResponses(
        simulation.acceptedPressureEpoch().pressureOperator,
        endpoint->controlCells, endpoint->fullTraceSystem,
        endpoint->condensedTraceSystem, acceptedSource, settings);
    fsi::validateSceneFluidPressureOperatorResponseAuditIntegrity(audit);
    check(audit == repeated && audit.includesAcceptedSource
              && audit.modes.size() == 7
              && audit.responses.size()
                  == 7 * endpoint->controlCells.controlCells.size()
              && audit.modes[0].bestFitShadowPressureScale > 2.0
              && audit.modes[0].bestFitShadowPressureScale < 3.0
              && std::abs(
                     audit.modes[0].bestFitShadowPressureScale
                     - 2.53935427074188)
                  < 1.0e-9
              && std::abs(
                     audit.modes[0].relativeBestFitShapeResidualL2
                     - 0.024412292727059)
                  < 1.0e-9,
          "accepted pressure source is a special high-gain cut-cell response mode");
    constexpr double expectedGains[] = {
        0.999104345858666,
        1.00488697756112,
        1.00722994317487,
        0.99966665199924,
        1.00469290834379,
    };
    constexpr double expectedResiduals[] = {
        0.086616065926422,
        0.012315965046049,
        0.04435725726353,
        0.037803572674267,
        0.164611785658836,
    };
    for (std::size_t mode = 1; mode + 1 < audit.modes.size(); ++mode) {
        check(audit.modes[mode].finite
                  && audit.modes[mode].bestFitShadowPressureScale > 0.9
                  && audit.modes[mode].bestFitShadowPressureScale < 1.1
                  && audit.modes[mode].relativeBestFitShapeResidualL2 < 0.18
                  && std::abs(
                         audit.modes[mode].bestFitShadowPressureScale
                         - expectedGains[mode - 1])
                      < 1.0e-9
                  && std::abs(
                         audit.modes[mode].relativeBestFitShapeResidualL2
                         - expectedResiduals[mode - 1])
                      < 1.0e-9,
              "independent manufactured pressure modes remain near unit gain with bounded shape disagreement");
    }
    const auto& regionContrast = audit.modes.back();
    check(regionContrast.kind
              == fsi::SceneFluidPressureOperatorResponseModeKind::
                  RegionContrast
              && regionContrast.finite
              && regionContrast.hasTwoTerminalConductance
              && regionContrast.lowerTerminalRegionId == 1
              && regionContrast.upperTerminalRegionId == 2
              && std::abs(
                     regionContrast.bestFitShadowPressureScale
                     - 2.56237302967949)
                  < 1.0e-9
              && std::abs(
                     regionContrast.relativeBestFitShapeResidualL2
                     - 0.0198229425861661)
                  < 1.0e-9
              && std::abs(
                     regionContrast.graphTwoTerminalConductanceMeters
                     - 0.18)
                  < 1.0e-12
              && std::abs(
                     regionContrast.shadowTwoTerminalConductanceMeters
                     - 0.0700820848334838)
                  < 1.0e-12
              && std::abs(
                     regionContrast.graphToShadowTwoTerminalConductanceRatio
                     - 2.56841674199167)
                  < 1.0e-9
              && std::abs(
                     regionContrast.shadowToGraphSourceComplianceRatio
                     - regionContrast
                         .graphToShadowTwoTerminalConductanceRatio)
                  < 1.0e-12
              && regionContrast.pressureCosineSimilarity > 0.9998,
          "authored-region contrast measures the intake-only graph and shadow two-terminal conductances");
    fsi::SceneFluidMimeticRegionConductanceAuditSettings
        terminalSettings;
    terminalSettings.solve = settings.shadowSolve;
    const auto terminal =
        fsi::auditSceneFluidMimeticRegionConductance(
            endpoint->controlCells, endpoint->fullTraceSystem,
            endpoint->condensedTraceSystem, terminalSettings);
    const auto repeatedTerminal =
        fsi::auditSceneFluidMimeticRegionConductance(
            endpoint->controlCells, endpoint->fullTraceSystem,
            endpoint->condensedTraceSystem, terminalSettings);
    fsi::validateSceneFluidMimeticRegionConductanceAuditIntegrity(
        terminal);
    check(terminal == repeatedTerminal
              && terminal.componentCount == 1
              && terminal.lowerTerminalRegionId == 1
              && terminal.upperTerminalRegionId == 2
              && terminal.openings.size() == 1
              && terminal.openings.front().traceKind
                  == fsi::SceneFluidMimeticHalfFaceKind::CartesianTrace
              && std::abs(terminal.openingAreaSquareMeters - 0.18)
                  < 1.0e-14
              && terminal.openings.front()
                     .integratedTransferPascalsMeters == 1.0
              && terminal.lowerTerminalIntegratedSourcePascalsMeters
                  == 1.0
              && terminal.upperTerminalIntegratedSourcePascalsMeters
                  == -1.0
              && terminal.componentIntegratedSourcePascalsMeters == 0.0
              && std::abs(
                     terminal.sourcePressureWorkPascalsSquaredMeters
                     - 14.2689818999465)
                  < 1.0e-12
              && std::abs(
                     terminal.conductanceMeters
                     - 0.0700820848335194)
                  < 1.0e-14
              && std::abs(
                     terminal.conductanceMeters
                     - regionContrast.shadowTwoTerminalConductanceMeters)
                  < 1.0e-12,
          "graph-independent terminal audit closes against the face-aligned shadow conductance");
    auto corruptTerminal = terminal;
    corruptTerminal.responses.front().gaugeAlignedPressurePascals += 0.01;
    bool terminalRejected = false;
    try {
        fsi::validateSceneFluidMimeticRegionConductanceAuditIntegrity(
            corruptTerminal);
    } catch (const std::exception&) {
        terminalRejected = true;
    }
    check(terminalRejected,
          "terminal conductance audit rejects response corruption");
    fsi::SceneFluidMimeticRegionConductanceAuditLimits terminalLimits;
    terminalLimits.maximumControlCells =
        endpoint->controlCells.controlCells.size() - 1;
    terminalRejected = false;
    try {
        static_cast<void>(
            fsi::auditSceneFluidMimeticRegionConductance(
                endpoint->controlCells, endpoint->fullTraceSystem,
                endpoint->condensedTraceSystem, terminalSettings,
                terminalLimits));
    } catch (const std::length_error&) {
        terminalRejected = true;
    }
    check(terminalRejected,
          "terminal conductance audit enforces its control-cell limit before solving");
    terminalLimits = {};
    terminalLimits.maximumOwnedBytes = 0;
    terminalRejected = false;
    try {
        static_cast<void>(
            fsi::auditSceneFluidMimeticRegionConductance(
                endpoint->controlCells, endpoint->fullTraceSystem,
                endpoint->condensedTraceSystem, terminalSettings,
                terminalLimits));
    } catch (const std::length_error&) {
        terminalRejected = true;
    }
    check(terminalRejected,
          "terminal conductance audit enforces its aggregate byte limit before solving");
    const auto manufacturedOnly =
        fsi::auditSceneFluidPressureOperatorResponses(
            simulation.acceptedPressureEpoch().pressureOperator,
            endpoint->controlCells, endpoint->fullTraceSystem,
            endpoint->condensedTraceSystem, {}, settings);
    check(!manufacturedOnly.includesAcceptedSource
              && manufacturedOnly.modes.size() == 6
              && manufacturedOnly.modes.front().kind
                  == fsi::SceneFluidPressureOperatorResponseModeKind::
                      CoordinateX,
          "operator-response audit supports a six-mode manufactured-only oracle");
    auto corrupt = audit;
    corrupt.responses.front()
        .shadowMinusBestFitGraphPressurePascals += 0.01;
    bool rejected = false;
    try {
        fsi::validateSceneFluidPressureOperatorResponseAuditIntegrity(
            corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "operator-response audit rejects response corruption");

    fsi::SceneFluidPressureOperatorResponseAuditLimits limited;
    limited.maximumModes = 6;
    rejected = false;
    try {
        static_cast<void>(fsi::auditSceneFluidPressureOperatorResponses(
            simulation.acceptedPressureEpoch().pressureOperator,
            endpoint->controlCells, endpoint->fullTraceSystem,
            endpoint->condensedTraceSystem, acceptedSource, settings,
            limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "operator-response audit enforces its aggregate mode limit before solving");
    limited.maximumModes = 64;
    limited.maximumOwnedBytes = 0;
    rejected = false;
    try {
        static_cast<void>(fsi::auditSceneFluidPressureOperatorResponses(
            simulation.acceptedPressureEpoch().pressureOperator,
            endpoint->controlCells, endpoint->fullTraceSystem,
            endpoint->condensedTraceSystem, acceptedSource, settings,
            limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "operator-response audit enforces its aggregate byte limit before solving");
}

void testRestOperatorRefinementAudit() {
    fsi::ScenePressureCellOperatorRefinementAuditSettings settings;
    settings.response = pressureCellResponseAuditSettings();
    const std::vector<fsi::fluid::GridCellCounts> resolutions{
        {2, 2, 2}, {4, 4, 4}, {8, 8, 8},
    };
    const auto audit = fsi::auditScenePressureCellOperatorRefinement(
        resolutions, settings);
    const auto repeated = fsi::auditScenePressureCellOperatorRefinement(
        resolutions, settings);
    fsi::validateScenePressureCellOperatorRefinementAuditIntegrity(audit);
    check(audit == repeated
              && audit.samples.size() == resolutions.size()
              && audit.structureDefinitionFingerprint != 0
              && audit.ownedStorageBytes > 0,
          "rest pressure-cell refinement audit owns three validated resolutions");
    constexpr double expectedGraphConductance[] = {
        0.306328421033639,
        2.58638935233124,
        1.90639953337707,
    };
    constexpr double expectedShadowConductance[] = {
        0.0109310650855869,
        0.0768124070857098,
        0.215905072160453,
    };
    constexpr double expectedRatios[] = {
        28.0236572223457,
        33.6715050401331,
        8.82980429454805,
    };
    constexpr double expectedNormalizedGraph[] = {
        2.85744545106773,
        12.0629787869036,
        4.44574539981043,
    };
    constexpr double expectedNormalizedShadow[] = {
        0.101965472543292,
        0.358254814346008,
        0.503493084501934,
    };
    for (std::size_t index = 0; index < audit.samples.size(); ++index) {
        const auto& sample = audit.samples[index];
        check(std::abs(
                  sample.graphConductanceMeters
                  - expectedGraphConductance[index]) < 1.0e-9
                  && std::abs(
                      sample.shadowConductanceMeters
                      - expectedShadowConductance[index]) < 1.0e-9
                  && std::abs(
                      sample.graphToShadowConductanceRatio
                      - expectedRatios[index]) < 1.0e-9
                  && std::abs(
                      sample.normalizedGraphConductance
                      - expectedNormalizedGraph[index]) < 1.0e-9
                  && std::abs(
                      sample.normalizedShadowConductance
                      - expectedNormalizedShadow[index]) < 1.0e-9,
              "skew-intake refinement retains exact graph and shadow conductance spectra");
    }
    check(audit.samples[0].normalizedShadowConductance
                  < audit.samples[1].normalizedShadowConductance
              && audit.samples[1].normalizedShadowConductance
                  < audit.samples[2].normalizedShadowConductance
              && audit.samples[0].graphToShadowConductanceRatio > 20.0
              && audit.samples[1].graphToShadowConductanceRatio > 20.0
              && audit.samples[2].graphToShadowConductanceRatio > 8.0,
          "skew-intake refinement exposes a smoother shadow trend while both operators remain far apart");

    auto corrupt = audit;
    corrupt.samples.back().normalizedShadowConductance += 0.01;
    bool rejected = false;
    try {
        fsi::validateScenePressureCellOperatorRefinementAuditIntegrity(
            corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "rest pressure-cell refinement audit rejects derived-value corruption");

    const std::vector<fsi::fluid::GridCellCounts> invalid{{5, 4, 5}};
    rejected = false;
    try {
        static_cast<void>(fsi::auditScenePressureCellOperatorRefinement(
            invalid, settings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "rest pressure-cell refinement audit rejects an anisotropic grid");

    fsi::ScenePressureCellOperatorRefinementAuditLimits limited;
    limited.maximumGridCellsPerSample = 7;
    rejected = false;
    try {
        static_cast<void>(fsi::auditScenePressureCellOperatorRefinement(
            std::span(resolutions).first(1), settings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "rest pressure-cell refinement audit enforces its grid-cell limit before assembly");
}

void testRestOperatorPhaseAudit() {
    const auto& phases = pressureCellGridPhases();
    fsi::ScenePressureCellOperatorPhaseAuditSettings phaseSettings;
    phaseSettings.response = pressureCellResponseAuditSettings();
    const auto phaseAudit =
        fsi::auditScenePressureCellOperatorGridPhases(
            phases, phaseSettings);
    const auto repeated =
        fsi::auditScenePressureCellOperatorGridPhases(
            phases, phaseSettings);
    fsi::validateScenePressureCellOperatorPhaseAuditIntegrity(phaseAudit);
    check(phaseAudit == repeated
              && phaseAudit.samples.size() == phases.size()
              && phaseAudit.structureDefinitionFingerprint != 0
              && phaseAudit.ownedStorageBytes > 0,
          "fixed-grid phase audit is deterministic and self-contained");
    check(phaseAudit.statistics.acceptedSampleCount == 6
              && phaseAudit.statistics
                  .rejectedIncompleteFaceOwnershipSampleCount == 2
              && std::abs(
                  phaseAudit.statistics.minimumNormalizedGraphConductance
                  - 2.5574976589654592) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.maximumNormalizedGraphConductance
                  - 13.985353433666207) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.meanNormalizedGraphConductance
                  - 6.2524874416504694) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.graphCoefficientOfVariation
                  - 0.77175443562108292) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.minimumNormalizedShadowConductance
                  - 0.13695429534604528) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.maximumNormalizedShadowConductance
                  - 0.40878148469499564) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.meanNormalizedShadowConductance
                  - 0.2690883993497567) < 1.0e-10
              && std::abs(
                  phaseAudit.statistics.shadowCoefficientOfVariation
                  - 0.34030273798624311) < 1.0e-10,
          "fixed-grid phase audit locks topology yield and graph/shadow placement spectra");
    constexpr std::size_t rejectedIndices[] = {2, 6};
    constexpr std::size_t expectedResolvedFull[] = {188, 180};
    constexpr std::size_t expectedResolvedPartition[] = {4, 12};
    constexpr std::size_t expectedUnresolvedOpeningPatch[] = {1, 2};
    for (std::size_t index = 0; index < 2; ++index) {
        const std::size_t rejectedIndex = rejectedIndices[index];
        const auto& sample = phaseAudit.samples[rejectedIndex];
        check(sample.status
                      == fsi::ScenePressureCellOperatorPhaseSampleStatus::
                          RejectedIncompleteFaceOwnership
                  && sample.faceOwnershipRejection.has_value()
                  && !sample.acceptedAudit.has_value(),
              "phase audit preserves typed incomplete-face topology rejection");
        if (sample.faceOwnershipRejection.has_value()) {
            const auto& rejection = *sample.faceOwnershipRejection;
            check(rejection.faceCount == 192
                      && rejection.resolvedFullFaceCount
                          == expectedResolvedFull[index]
                      && rejection.resolvedPartitionFaceCount
                          == expectedResolvedPartition[index]
                      && rejection.resolvedOpeningFaceCount == 0
                      && rejection.unresolvedActiveFaceCount == 0
                      && rejection.unresolvedCappedFaceCount == 0
                      && rejection.unresolvedAmbiguousFaceCount == 0
                      && rejection.unresolvedOpeningFaceCount == 0
                      && rejection.unresolvedEmbeddedOpeningPatchCount
                          == expectedUnresolvedOpeningPatch[index],
                  "phase audit retains exact unresolved embedded-opening diagnostics");
        }
    }

    auto corrupt = phaseAudit;
    corrupt.statistics.graphCoefficientOfVariation += 0.01;
    bool rejected = false;
    try {
        fsi::validateScenePressureCellOperatorPhaseAuditIntegrity(corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "phase audit rejects aggregate-statistic corruption");

    corrupt = phaseAudit;
    corrupt.samples[rejectedIndices[0]].faceOwnershipRejection
        ->unresolvedEmbeddedOpeningPatchCount = 0;
    rejected = false;
    try {
        fsi::validateScenePressureCellOperatorPhaseAuditIntegrity(corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "phase audit rejects a falsely complete ownership rejection");

    const std::vector<fsi::fluid::Vector3> invalidPhase{{1.0, 0.0, 0.0}};
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorGridPhases(
                invalidPhase, phaseSettings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "phase audit rejects a non-canonical unit-cell phase");

    const std::vector<fsi::fluid::Vector3> duplicatePhases{
        {0.25, 0.25, 0.25}, {0.25, 0.25, 0.25},
    };
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorGridPhases(
                duplicatePhases, phaseSettings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "phase audit rejects duplicate placement samples");

    fsi::ScenePressureCellOperatorPhaseAuditLimits limited;
    limited.maximumGridCellsPerSample = 63;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorGridPhases(
                std::span(phases).first(1), phaseSettings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "phase audit enforces its grid limit before topology assembly");

    limited = {};
    limited.maximumOwnedBytes = 0;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorGridPhases(
                std::span(phases).first(1), phaseSettings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "phase audit enforces its aggregate byte limit before publication");

    limited = {};
    limited.maximumSamples = 0;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorGridPhases(
                std::span(phases).first(1), phaseSettings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "phase audit enforces its sample limit before assembly");

}

void testRestOperatorPhaseRefinementAudit() {
    const std::vector<fsi::fluid::GridCellCounts> resolutions{
        {2, 2, 2}, {4, 4, 4}, {8, 8, 8},
    };
    const auto& phases = pressureCellGridPhases();
    fsi::ScenePressureCellOperatorPhaseRefinementAuditSettings settings;
    settings.response = pressureCellResponseAuditSettings();
    const auto audit =
        fsi::auditScenePressureCellOperatorPhaseRefinement(
            resolutions, phases, settings);
    const auto repeated =
        fsi::auditScenePressureCellOperatorPhaseRefinement(
            resolutions, phases, settings);
    fsi::validateScenePressureCellOperatorPhaseRefinementAuditIntegrity(
        audit);
    check(audit == repeated
              && audit.levels.size() == resolutions.size()
              && audit.gridPhaseFractions == phases
              && audit.structureDefinitionFingerprint != 0
              && audit.ownedStorageBytes > 0,
          "phase-refinement audit is deterministic and owns every nested phase spectrum");

    constexpr std::size_t expectedAccepted[] = {4, 6, 2};
    constexpr std::size_t expectedRejected[] = {4, 2, 6};
    constexpr double expectedTopologyFraction[] = {0.5, 0.75, 0.25};
    constexpr double expectedGraphMinimum[] = {
        2.7504856987378798,
        2.5574976589654592,
        4.4457453998104262,
    };
    constexpr double expectedGraphMaximum[] = {
        5.1709987460683449,
        13.985353433666207,
        4.5916922290114295,
    };
    constexpr double expectedGraphMean[] = {
        3.9024884705941725,
        6.2524874416504694,
        4.5187188144109278,
    };
    constexpr double expectedGraphVariation[] = {
        0.28333875171817424,
        0.77175443562108292,
        0.016149138195494173,
    };
    constexpr double expectedShadowMinimum[] = {
        0.10139886521551159,
        0.13695429534604528,
        0.41157189435124092,
    };
    constexpr double expectedShadowMaximum[] = {
        0.10196547254329238,
        0.40878148469499564,
        0.50349308450193431,
    };
    constexpr double expectedShadowMean[] = {
        0.10168136399302657,
        0.2690883993497567,
        0.45753248942658764,
    };
    constexpr double expectedShadowVariation[] = {
        0.0023598291157185796,
        0.34030273798624311,
        0.10045318340769153,
    };
    for (std::size_t index = 0; index < audit.levels.size(); ++index) {
        const auto& level = audit.levels[index];
        const auto& statistics = level.phaseAudit.statistics;
        check(statistics.acceptedSampleCount == expectedAccepted[index]
                  && statistics
                      .rejectedIncompleteFaceOwnershipSampleCount
                      == expectedRejected[index]
                  && level.acceptedTopologyFraction
                      == expectedTopologyFraction[index]
                  && std::abs(
                      statistics.minimumNormalizedGraphConductance
                      - expectedGraphMinimum[index]) < 1.0e-9
                  && std::abs(
                      statistics.maximumNormalizedGraphConductance
                      - expectedGraphMaximum[index]) < 1.0e-9
                  && std::abs(
                      statistics.meanNormalizedGraphConductance
                      - expectedGraphMean[index]) < 1.0e-9
                  && std::abs(
                      statistics.graphCoefficientOfVariation
                      - expectedGraphVariation[index]) < 1.0e-9
                  && std::abs(
                      statistics.minimumNormalizedShadowConductance
                      - expectedShadowMinimum[index]) < 1.0e-9
                  && std::abs(
                      statistics.maximumNormalizedShadowConductance
                      - expectedShadowMaximum[index]) < 1.0e-9
                  && std::abs(
                      statistics.meanNormalizedShadowConductance
                      - expectedShadowMean[index]) < 1.0e-9
                  && std::abs(
                      statistics.shadowCoefficientOfVariation
                      - expectedShadowVariation[index]) < 1.0e-9,
              "phase-refinement audit locks topology yield and conductance distributions");
        for (const auto& sample : level.phaseAudit.samples) {
            if (!sample.faceOwnershipRejection.has_value()) {
                continue;
            }
            const auto& rejection = *sample.faceOwnershipRejection;
            check(rejection.resolvedFullFaceCount
                          + rejection.resolvedPartitionFaceCount
                          + rejection.resolvedOpeningFaceCount
                      == rejection.faceCount
                      && rejection.unresolvedActiveFaceCount == 0
                      && rejection.unresolvedCappedFaceCount == 0
                      && rejection.unresolvedAmbiguousFaceCount == 0
                      && rejection.unresolvedOpeningFaceCount == 0
                      && rejection.unresolvedEmbeddedOpeningPatchCount > 0,
                  "phase-refinement rejection isolates embedded intake ownership");
        }
    }
    check(audit.levels[0].phaseAudit.statistics
                  .meanNormalizedShadowConductance
                  < audit.levels[1].phaseAudit.statistics
                      .meanNormalizedShadowConductance
              && audit.levels[1].phaseAudit.statistics
                  .meanNormalizedShadowConductance
                  < audit.levels[2].phaseAudit.statistics
                      .meanNormalizedShadowConductance
              && audit.levels[0].acceptedTopologyFraction
                  > audit.levels[2].acceptedTopologyFraction,
          "phase-refinement audit exposes a monotone shadow mean but deteriorating fine-grid topology yield");

    auto corrupt = audit;
    corrupt.levels[1].acceptedTopologyFraction += 0.01;
    bool rejected = false;
    try {
        fsi::validateScenePressureCellOperatorPhaseRefinementAuditIntegrity(
            corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "phase-refinement audit rejects derived topology-yield corruption");

    const std::vector<fsi::fluid::GridCellCounts> unordered{
        {4, 4, 4}, {2, 2, 2},
    };
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorPhaseRefinement(
                unordered, phases, settings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "phase-refinement audit rejects non-increasing resolutions before assembly");

    fsi::ScenePressureCellOperatorPhaseRefinementAuditLimits limited;
    limited.maximumAggregatePhaseSamples = 23;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorPhaseRefinement(
                resolutions, phases, settings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "phase-refinement audit enforces its aggregate sample limit before assembly");

    limited = {};
    limited.maximumOwnedBytes = 0;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellOperatorPhaseRefinement(
                resolutions, phases, settings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "phase-refinement audit enforces its byte limit before assembly");
}

void testUncensoredMimeticConductancePhaseRefinementAudit() {
    const std::vector<fsi::fluid::GridCellCounts> resolutions{
        {2, 2, 2}, {4, 4, 4}, {8, 8, 8},
    };
    const auto& phases = pressureCellGridPhases();
    fsi::ScenePressureCellMimeticConductancePhaseRefinementAuditSettings
        settings;
    settings.conductance.solve =
        pressureCellResponseAuditSettings().shadowSolve;
    const auto audit =
        fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
            resolutions, phases, settings);
    const auto repeated =
        fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
            resolutions, phases, settings);
    fsi::validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        audit);
    check(audit == repeated
              && audit.levels.size() == resolutions.size()
              && audit.gridPhaseFractions == phases
              && audit.structureDefinitionFingerprint != 0
              && audit.ownedStorageBytes > 0,
          "uncensored mimetic phase-refinement audit is deterministic and self-contained");
    constexpr std::size_t expectedAccepted[]{8, 8, 8};
    constexpr std::size_t expectedRejected[]{0, 0, 0};
    constexpr double expectedMinimum[]{
        0.099347529340182233,
        0.12925278175548274,
        0.29373390290273943,
    };
    constexpr double expectedMaximum[]{
        0.10196547254313142,
        0.408781484694966,
        0.88867102845759738,
    };
    constexpr double expectedMean[]{
        0.10066047211051987,
        0.24092951956532438,
        0.52124824069717568,
    };
    constexpr double expectedVariation[]{
        0.010428734257699697,
        0.39050354377024471,
        0.38827709459290677,
    };
    constexpr double expectedNormalized[][8]{
        {
            0.10196547254313142, 0.099933389492765257,
            0.10149376949041197, 0.10186734872294091,
            0.099443375352812405, 0.099834026726442912,
            0.10139886521547191, 0.099347529340182233,
        },
        {
            0.35825481434598105, 0.408781484694966,
            0.18244616118707896, 0.24123026839889594,
            0.20221810122979983, 0.26800858093279595,
            0.12925278175548274, 0.13724396397759467,
        },
        {
            0.51940722905705417, 0.35090782886457544,
            0.88867102845759738, 0.44826142651391082,
            0.44945862708767215, 0.29373390290273943,
            0.81684064037937898, 0.40270524231447685,
        },
    };
    constexpr std::size_t expectedOpeningCounts[][8]{
        {1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 2, 2, 2, 2, 4, 4},
        {4, 6, 4, 6, 6, 8, 7, 9},
    };
    for (std::size_t levelIndex = 0;
         levelIndex < audit.levels.size(); ++levelIndex) {
        const auto& level = audit.levels[levelIndex];
        check(level.samples.size() == phases.size()
                  && level.acceptedSampleCount
                      == expectedAccepted[levelIndex]
                  && level.rejectedLocalCellLinearConsistencySampleCount
                      == expectedRejected[levelIndex]
                  && std::abs(
                      level.minimumNormalizedConductance
                      - expectedMinimum[levelIndex]) < 1.0e-12
                  && std::abs(
                      level.maximumNormalizedConductance
                      - expectedMaximum[levelIndex]) < 1.0e-12
                  && std::abs(
                      level.meanNormalizedConductance
                      - expectedMean[levelIndex]) < 1.0e-12
                  && std::abs(
                      level.normalizedConductanceCoefficientOfVariation
                      - expectedVariation[levelIndex]) < 1.0e-12,
              "uncensored mimetic audit locks each phase population and conditional spectrum");
        for (std::size_t phaseIndex = 0;
             phaseIndex < level.samples.size(); ++phaseIndex) {
            const auto& sample = level.samples[phaseIndex];
            check(sample.openingTraceCount
                          == expectedOpeningCounts[levelIndex][phaseIndex]
                      && std::abs(
                          sample.normalizedConductance
                          - expectedNormalized[levelIndex][phaseIndex])
                          < 1.0e-12,
                  "uncensored mimetic audit retains each area-weighted phase response");
            check(sample.status
                          == fsi::ScenePressureCellMimeticConductancePhaseSampleStatus::Accepted
                      && sample.conductanceAudit.has_value()
                      && !sample.localCellLinearConsistencyRejection
                              .has_value(),
                  "uncensored mimetic audit publishes every accepted terminal solve");
        }
    }
    check(audit.levels[0].acceptedSampleCount == phases.size()
              && audit.levels[1].acceptedSampleCount == phases.size()
              && audit.levels[2].acceptedSampleCount == phases.size()
              && audit.levels[0].meanNormalizedConductance
                  < audit.levels[1].meanNormalizedConductance
              && audit.levels[1].meanNormalizedConductance
                  < audit.levels[2].meanNormalizedConductance,
          "graph-independent source exposes the complete phase spectrum across refinement");

    auto translatedSettings = settings;
    translatedSettings.geometryTranslationMeters = {
        256.0, -512.0, 1024.0};
    const std::array<fsi::fluid::GridCellCounts, 1>
        translatedResolution{{{8, 8, 8}}};
    const std::array<fsi::fluid::Vector3, 1> translatedPhase{{
        phases[2],
    }};
    const auto translated =
        fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
            translatedResolution, translatedPhase, translatedSettings);
    fsi::validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        translated);
    const auto& originSample = audit.levels[2].samples[2];
    const auto& translatedSample = translated.levels[0].samples[0];
    check(translated.structureDefinitionFingerprint
                  != audit.structureDefinitionFingerprint
              && translatedSample.status
                  == fsi::ScenePressureCellMimeticConductancePhaseSampleStatus::Accepted
              && translatedSample.openingTraceCount
                  == originSample.openingTraceCount
              && translatedSample.controlVolumeCount
                  == originSample.controlVolumeCount
              && translatedSample.fullTraceCount
                  == originSample.fullTraceCount
              && translatedSample.reducedTraceCount
                  == originSample.reducedTraceCount
              && std::abs(
                  translatedSample.intakeAreaSquareMeters
                  - originSample.intakeAreaSquareMeters) < 1.0e-12
              && std::abs(
                  translatedSample.normalizedConductance
                  - originSample.normalizedConductance) < 5.0e-11,
          "graph-independent shadow response is invariant under a distant common coordinate translation");

    const std::array<fsi::fluid::GridCellCounts, 1>
        fineResolution{{{16, 16, 16}}};
    const auto fineAudit =
        fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
            fineResolution, phases, settings);
    fsi::validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        fineAudit);
    const auto& fineLevel = fineAudit.levels.front();
    constexpr std::size_t fineOpeningCounts[]{
        15, 13, 14, 11, 11, 10, 11, 9,
    };
    constexpr double fineNormalized[]{
        0.62454444526833486,
        0.9992120361521708,
        0.69710200382794307,
        0.78394637046749738,
        1.00214013973588,
        1.141234407607552,
        0.845443837466477,
        1.1269364948078211,
    };
    check(fineLevel.samples.size() == phases.size()
              && fineLevel.acceptedSampleCount == phases.size()
              && fineLevel.rejectedLocalCellLinearConsistencySampleCount == 0
              && std::abs(fineLevel.minimumNormalizedConductance
                          - 0.62454444526833486) < 1.0e-12
              && std::abs(fineLevel.maximumNormalizedConductance
                          - 1.141234407607552) < 1.0e-12
              && std::abs(fineLevel.meanNormalizedConductance
                          - 0.90256996691670954) < 1.0e-12
              && std::abs(
                  fineLevel.normalizedConductanceCoefficientOfVariation
                  - 0.20104176667262066) < 1.0e-12,
          "16-cubed mimetic audit retains the complete phase spectrum");
    for (std::size_t phaseIndex = 0;
         phaseIndex < fineLevel.samples.size(); ++phaseIndex) {
        const auto& sample = fineLevel.samples[phaseIndex];
        check(sample.status
                      == fsi::ScenePressureCellMimeticConductancePhaseSampleStatus::Accepted
                  && sample.conductanceAudit.has_value()
                  && !sample.localCellLinearConsistencyRejection.has_value()
                  && sample.openingTraceCount
                      == fineOpeningCounts[phaseIndex]
                  && std::abs(sample.normalizedConductance
                              - fineNormalized[phaseIndex]) < 1.0e-12,
              "16-cubed mimetic audit locks each uncensored phase response");
    }
    check(fineLevel.normalizedConductanceCoefficientOfVariation
                  < 0.55 * audit.levels[2]
                      .normalizedConductanceCoefficientOfVariation,
          "uncensored phase sensitivity contracts materially from 8 to 16 cubed");

    auto corrupt = audit;
    corrupt.levels[2].samples[1].normalizedConductance += 1.0e-6;
    bool rejected = false;
    try {
        fsi::validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
            corrupt);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected,
          "uncensored mimetic audit rejects accepted-sample corruption");

    fsi::ScenePressureCellMimeticConductancePhaseRefinementAuditLimits
        limited;
    limited.maximumAggregateSamples = 23;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
                resolutions, phases, settings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "uncensored mimetic audit enforces its aggregate sample limit before assembly");
    limited = {};
    limited.maximumOwnedBytes = 0;
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
                resolutions, phases, settings, limited));
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected,
          "uncensored mimetic audit enforces its aggregate byte limit before assembly");

    auto invalidSettings = settings;
    invalidSettings.geometryTranslationMeters.x =
        std::numeric_limits<double>::infinity();
    rejected = false;
    try {
        static_cast<void>(
            fsi::auditScenePressureCellMimeticConductancePhaseRefinement(
                translatedResolution, translatedPhase, invalidSettings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "uncensored mimetic audit rejects a non-finite diagnostic translation");
}

void testPersistentMimeticPressureAuditRestart() {
    fsi::ScenePressureCellCase initial(true);
    const auto initialCheckpoint = initial.checkpoint();
    const auto initialBytes = serializedCheckpoint(initialCheckpoint);
    fsi::ScenePressureCellCheckpoint decodedInitial;
    fsi::ScenePressureCellCheckpointPersistenceError error;
    check(initialCheckpoint.coupling
                  .mimeticPressureAuditSettingsFingerprint != 0
              && !initialCheckpoint.coupling
                      .mimeticPressureState.has_value()
              && fsi::deserializeScenePressureCellCheckpoint(
                  initialBytes, decodedInitial, &error),
          "initial audited checkpoint round-trips its mode without inventing pressure state");
    fsi::ScenePressureCellCase initialReplay(true);
    initialReplay.restore(decodedInitial);
    check(serialized(initialReplay.advance())
              == serialized(initial.advance()),
          "initial audited checkpoint reproduces the first shadowed frame");

    fsi::ScenePressureCellCase source(true);
    for (std::size_t step = 0; step < 8; ++step) {
        static_cast<void>(source.advance());
    }
    const auto saved = source.checkpoint();
    const auto bytes = serializedCheckpoint(saved);
    fsi::ScenePressureCellCheckpoint decoded;
    check(saved.coupling.mimeticPressureState.has_value()
              && fsi::deserializeScenePressureCellCheckpoint(
                  bytes, decoded, &error)
              && decoded.coupling.mimeticPressureState
                  == saved.coupling.mimeticPressureState
              && serializedCheckpoint(decoded) == bytes,
          "audited pressure-cell checkpoint composes canonical compact SWMP state with the graph restart");

    fsi::ScenePressureCellCase restored(true);
    restored.restore(decoded);
    check(restored.acceptedMimeticPressureAudit() == nullptr
              && restored.acceptedMimeticPressureComparison() == nullptr
              && restored.acceptedMimeticPressureOwnerTransition() == nullptr
              && serializedCheckpoint(restored.checkpoint()) == bytes,
          "audited restore retains rebuilt warm topology without fabricating transient endpoint diagnostics");
    const auto expectedFrame = source.advance();
    const auto expectedDiagnostics = source.diagnostics();
    const auto expectedEndpoint = *source.acceptedMimeticPressureAudit();
    const auto expectedComparison =
        *source.acceptedMimeticPressureComparison();
    const auto expectedOwnerTransition =
        *source.acceptedMimeticPressureOwnerTransition();
    const auto replayFrame = restored.advance();
    const auto* replayEndpoint = restored.acceptedMimeticPressureAudit();
    const auto* replayComparison =
        restored.acceptedMimeticPressureComparison();
    const auto* replayOwnerTransition =
        restored.acceptedMimeticPressureOwnerTransition();
    check(serialized(replayFrame) == serialized(expectedFrame)
              && restored.diagnostics() == expectedDiagnostics
              && replayEndpoint != nullptr
              && *replayEndpoint == expectedEndpoint
              && replayComparison != nullptr
              && *replayComparison == expectedComparison
              && replayOwnerTransition != nullptr
              && *replayOwnerTransition == expectedOwnerTransition
              && replayEndpoint->usesConsecutiveWarmStart
              && replayEndpoint->usesRegionWallPrediction,
          "restored compact SWMP state reproduces the exact next consecutive wall-predicted endpoint");

    fsi::ScenePressureCellCase production;
    const auto productionBefore = serializedCheckpoint(
        production.checkpoint());
    bool rejected = false;
    try {
        production.restore(decoded);
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected
              && serializedCheckpoint(production.checkpoint())
                  == productionBefore,
          "default pressure-cell owner rejects an audited checkpoint transactionally");
    fsi::ScenePressureCellCase auditedDestination(true);
    const auto auditedBefore = serializedCheckpoint(
        auditedDestination.checkpoint());
    rejected = false;
    try {
        auditedDestination.restore(production.checkpoint());
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected
              && serializedCheckpoint(auditedDestination.checkpoint())
                  == auditedBefore,
          "audited pressure-cell owner rejects a graph-only checkpoint transactionally");

    auto corrupt = saved;
    corrupt.coupling.mimeticPressureState->controls.front()
        .pressurePascals += 0.01;
    std::vector<std::uint8_t> unchanged{1, 2, 3};
    const auto original = unchanged;
    check(!fsi::serializeScenePressureCellCheckpoint(
              corrupt, unchanged, &error)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::
                      InvalidData
              && unchanged == original,
          "composite checkpoint rejects corrupt SWMP state without publishing bytes");

    auto limits = fsi::ScenePressureCellCheckpointPersistenceLimits{};
    limits.mimeticPressureState.maximumControlCells =
        saved.coupling.mimeticPressureState->controls.size() - 1;
    auto preserved = decoded;
    check(!fsi::deserializeScenePressureCellCheckpoint(
              bytes, preserved, &error, limits)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::
                      LimitExceeded
              && serializedCheckpoint(preserved) == bytes,
          "composite checkpoint bounds nested SWMP rows before publication");
}

} // namespace

int main() {
    try {
        testVisibleStrongPressureCellAndReplay();
        testPersistentCheckpointAndRejection();
        testOptInMimeticPressureAuditIsShadowOnly();
        testManufacturedPressureOperatorResponses();
        testRestOperatorRefinementAudit();
        testRestOperatorPhaseAudit();
        testRestOperatorPhaseRefinementAudit();
        testUncensoredMimeticConductancePhaseRefinementAudit();
        testPersistentMimeticPressureAuditRestart();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene pressure cell check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene pressure cell checks passed");
    return 0;
}
