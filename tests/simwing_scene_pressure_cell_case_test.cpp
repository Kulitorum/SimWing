#include "scene_pressure_cell_case.h"
#include "scene_pressure_cell_checkpoint_persistence.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
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
                  && audited.diagnostics().coupling
                         .mimeticPressureComparisonFingerprint
                      == comparison->fingerprint
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
        if (endpoint != nullptr && comparison != nullptr) {
            fsi::validateSceneFluidMimeticPressureAuditEndpointIntegrity(
                *endpoint);
            fsi::validateSceneFluidPressureShadowComparisonIntegrity(
                *comparison);
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
    auto corrupted = *audited.acceptedMimeticPressureAudit();
    ++corrupted.pressureEpoch.diagnostics.pressureSolve
          .reducedTraceSolve.iterationCount;
    bool rejected = false;
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
              && serializedCheckpoint(sourceCountLimited.checkpoint())
                  == beforeSourceCountRejection,
          "mimetic pressure source-comparison count failure rolls Structure and every accepted owner back transactionally");
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
              && serializedCheckpoint(restored.checkpoint()) == bytes,
          "audited restore retains rebuilt warm topology without fabricating transient endpoint diagnostics");
    const auto expectedFrame = source.advance();
    const auto expectedDiagnostics = source.diagnostics();
    const auto expectedEndpoint = *source.acceptedMimeticPressureAudit();
    const auto expectedComparison =
        *source.acceptedMimeticPressureComparison();
    const auto replayFrame = restored.advance();
    const auto* replayEndpoint = restored.acceptedMimeticPressureAudit();
    const auto* replayComparison =
        restored.acceptedMimeticPressureComparison();
    check(serialized(replayFrame) == serialized(expectedFrame)
              && restored.diagnostics() == expectedDiagnostics
              && replayEndpoint != nullptr
              && *replayEndpoint == expectedEndpoint
              && replayComparison != nullptr
              && *replayComparison == expectedComparison
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
