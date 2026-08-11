#include "frozen_scene_pressure_case.h"
#include "scene_pressure_cell_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <stdexcept>
#include <vector>

namespace {

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

const simwing::viewer::ScalarField* scalarField(
    const simwing::viewer::DiagnosticFrame& frame,
    const char* name) {
    for (const auto& field : frame.scalarFields) {
        if (field.name == name) {
            return &field;
        }
    }
    return nullptr;
}

const simwing::viewer::VectorField* vectorField(
    const simwing::viewer::DiagnosticFrame& frame,
    const char* name) {
    for (const auto& field : frame.vectorFields) {
        if (field.name == name) {
            return &field;
        }
    }
    return nullptr;
}

std::vector<std::uint8_t> serialized(
    const simwing::viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> result;
    simwing::viewer::ProtocolError error;
    if (!simwing::viewer::serializeFrame(frame, result, &error)) {
        throw std::runtime_error(error.message);
    }
    return result;
}

bool sameGeometry(const simwing::viewer::DiagnosticFrame& first,
                  const simwing::viewer::DiagnosticFrame& second) {
    if (first.vertices.size() != second.vertices.size()
        || first.triangles.size() != second.triangles.size()) {
        return false;
    }
    for (std::size_t index = 0; index < first.vertices.size(); ++index) {
        const auto& a = first.vertices[index];
        const auto& b = second.vertices[index];
        if (a.stableId != b.stableId
            || a.positionMetres.x != b.positionMetres.x
            || a.positionMetres.y != b.positionMetres.y
            || a.positionMetres.z != b.positionMetres.z) {
            return false;
        }
    }
    for (std::size_t index = 0; index < first.triangles.size(); ++index) {
        const auto& a = first.triangles[index];
        const auto& b = second.triangles[index];
        if (a.stableId != b.stableId
            || a.vertex0 != b.vertex0 || a.vertex1 != b.vertex1
            || a.vertex2 != b.vertex2
            || a.negativeRegionId != b.negativeRegionId
            || a.positiveRegionId != b.positiveRegionId) {
            return false;
        }
    }
    return true;
}

void testFrozenScenePressureCase() {
    const auto scene = simwing::fsi::makeScenePressureCellGeometry();
    simwing::fsi::FrozenScenePressureCaseSettings settings;
    settings.cellCounts = {4, 4, 4};
    settings.useExplicitDomain = true;
    settings.lowerMeters = {};
    settings.upperMeters = {4.0, 4.0, 4.0};
    check(!settings.useCorrectedTraceFlowContinuation,
          "corrected trace-flow continuation remains opt-in");
    check(!settings.useRegionalTransportFlowPrediction,
          "regional transport flow prediction remains opt-in");
    check(!settings.useMovingGeometryFsi,
          "moving geometry FSI remains opt-in");
    settings.useRegionalTransportFlowPrediction = true;
    simwing::fsi::FrozenScenePressureCase first(scene, settings);
    simwing::fsi::FrozenScenePressureCase second(scene, settings);
    const auto firstFrame = first.advance();
    const auto repeatedFrame = second.advance();
    const auto firstDiagnostics = first.diagnostics();
    const auto secondFrame = first.advance();
    const auto repeatedSecondFrame = second.advance();
    const auto* pressure = scalarField(
        firstFrame, "frozen_scene.pressure_jump");
    const auto* nodalForce = vectorField(
        firstFrame, "frozen_scene.nodal_pressure_force");
    const auto* fluidVelocity = vectorField(
        firstFrame, "frozen_scene.fluid_velocity");
    const auto* secondPressure = scalarField(
        secondFrame, "frozen_scene.pressure_jump");
    const auto& diagnostics = first.diagnostics();
    check(serialized(firstFrame) == serialized(repeatedFrame)
              && serialized(secondFrame) == serialized(repeatedSecondFrame)
              && firstFrame.step == 1
              && secondFrame.step == 2
              && sameGeometry(firstFrame, secondFrame)
              && secondPressure != nullptr
              && pressure != nullptr
              && secondPressure->values != pressure->values
              && firstDiagnostics.flowAdvanceCount == 0
              && firstDiagnostics.windRampSeconds
                  == settings.windRampSeconds
              && firstDiagnostics.windRampFraction > 0.0
              && diagnostics.windRampFraction
                  > firstDiagnostics.windRampFraction
              && diagnostics.windRampFraction <= 1.0
              && diagnostics.flowAdvanceCount == 1
              && !diagnostics.usesCorrectedTraceFlowContinuation
              && diagnostics.usesRegionalTransportFlowPrediction
              && diagnostics
                     .maximumRegionalTransportFlowDifferenceFromBulkBaselineCubicMetersPerSecond
                  > 0.0
              && diagnostics
                     .maximumCarriedTraceCorrectionCubicMetersPerSecond
                  == 0.0
              && diagnostics.maximumTraceBulkIncrementCubicMetersPerSecond
                  == 0.0,
          "frozen scene flow evolves deterministically while leaving geometry unchanged");
    check(pressure != nullptr
              && pressure->association
                  == simwing::viewer::FieldAssociation::Triangle
              && pressure->values.size() == firstFrame.triangles.size()
              && nodalForce != nullptr
              && nodalForce->association
                  == simwing::viewer::FieldAssociation::Vertex
              && nodalForce->values.size() == firstFrame.vertices.size()
              && fluidVelocity != nullptr
              && fluidVelocity->association
                  == simwing::viewer::FieldAssociation::Vertex
              && fluidVelocity->values.size() == firstFrame.vertices.size()
              && std::any_of(
                  fluidVelocity->values.begin(), fluidVelocity->values.end(),
                  [](const auto value) {
                      return value.x != 0.0 || value.y != 0.0
                          || value.z != 0.0;
                  }),
          "frozen scene frame publishes complete pressure, nodal load, and fluid velocity fields");
    check(diagnostics.finite
              && diagnostics.gridCellCounts
                  == simwing::fsi::fluid::GridCellCounts{4, 4, 4}
              && diagnostics.pressureControlCount > 0
              && diagnostics.sharedTraceCount > 0
              && diagnostics.pressureIterationCount > 0
              && diagnostics.extrapolatedZeroVolumePressureSideCount == 0
              && diagnostics.maximumPressureExtrapolationDistanceMeters
                  == 0.0
              && diagnostics.maximumAbsolutePressureDifferencePascals > 0.0
              && diagnostics.reconstructedMaterialWallPressureSideCount > 0
              && diagnostics.fallbackMaterialWallPressureSideCount == 0
              && diagnostics
                     .maximumAbsoluteMaterialWallTracePressureDifferencePascals
                  > 0.0
              && diagnostics
                     .maximumAbsoluteMaterialWallTraceDifferenceFromControlSamplePascals
                  > 0.0
              && diagnostics
                     .materialWallTracePressureTransferForceResidualNewtons
                  < 1.0e-8
              && diagnostics
                     .materialWallTracePressureTransferMomentResidualNewtonMeters
                  < 1.0e-8
              && diagnostics
                     .maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond
                  <= diagnostics.correctedContinuityToleranceCubicMetersPerSecond
              && diagnostics.maximumCollapsedMacVelocityMetersPerSecond > 0.0
              && diagnostics.maximumRegionalVelocityMetersPerSecond > 0.0
              && diagnostics
                     .maximumRegionalLinkVelocityResidualMetersPerSecond
                  >= 0.0
              && diagnostics.regionalKineticEnergyJoules > 0.0
              && diagnostics.regionalTransportSubstepCount > 0
              && diagnostics
                     .regionalTransportMaximumVelocityChangeMetersPerSecond
                  > 0.0
              && diagnostics
                     .regionalTransportMomentumResidualKilogramMetersPerSecond
                  < 1.0e-10
              && diagnostics.regionalTransportAdvectiveEnergyLossJoules
                  >= 0.0
              && diagnostics.regionalTransportViscousEnergyLossJoules
                  >= 0.0
              && diagnostics.bulkFlowSubstepCount > 0
              && diagnostics.bulkProjectionDivergenceAfterPerSecond
                  < diagnostics.bulkProjectionDivergenceBeforePerSecond
              && diagnostics.transferForceResidualNewtons < 1.0e-8
              && diagnostics.transferMomentResidualNewtonMeters < 1.0e-8,
          "frozen scene pressure solve and conservative transfer are accepted");

    settings.windRampSeconds = 0.0;
    simwing::fsi::FrozenScenePressureCase impulsive(scene, settings);
    static_cast<void>(impulsive.advance());
    check(impulsive.diagnostics().windRampFraction == 1.0,
          "zero ramp duration explicitly restores full-wind startup");

    simwing::fsi::FrozenScenePressureCaseSettings movingSettings;
    movingSettings.cellCounts = {4, 4, 4};
    movingSettings.useExplicitDomain = true;
    movingSettings.lowerMeters = {};
    movingSettings.upperMeters = {4.0, 4.0, 4.0};
    movingSettings.backgroundWindMetersPerSecond = {0.1, 0.0, 0.0};
    movingSettings.windRampSeconds = 0.0;
    movingSettings.timeStepSeconds = 0.01;
    simwing::fsi::FrozenScenePressureCase frozenReference(
        scene, movingSettings);
    const auto frozenReferenceFrame = frozenReference.advance();
    movingSettings.useMovingGeometryFsi = true;
    auto prestrainedScene = scene;
    for (auto& triangle : prestrainedScene.triangles) {
        for (auto& coordinate : triangle.materialCoordinates) {
            coordinate.x *= 0.5;
            coordinate.y *= 0.5;
        }
    }
    bool prestrainedRejected = false;
    try {
        simwing::fsi::FrozenScenePressureCase rejected(
            prestrainedScene, movingSettings);
    } catch (const std::runtime_error& exception) {
        prestrainedRejected = std::string(exception.what()).find(
            "initial Structure rest audit rejected") != std::string::npos;
    }
    check(prestrainedRejected,
          "moving scene FSI rejects a nonstationary imported rest state before CFD assembly");

    simwing::fsi::FrozenScenePressureCase movingFirst(
        scene, movingSettings);
    simwing::fsi::FrozenScenePressureCase movingSecond(
        scene, movingSettings);
    const auto movingFrame = movingFirst.advance();
    const auto repeatedMovingFrame = movingSecond.advance();
    const auto& movingDiagnostics = movingFirst.diagnostics();
    check(serialized(movingFrame) == serialized(repeatedMovingFrame)
              && movingFirst.traceHeader().solverCommit
                  == simwing::fsi::movingScenePressureSolverId
              && movingFrame.step == 1
              && movingFirst.acceptedStepCount() == 1
              && !sameGeometry(frozenReferenceFrame, movingFrame),
          "moving scene FSI deterministically publishes its first deformed Structure frame");
    check(movingDiagnostics.finite
              && movingDiagnostics.usesMovingGeometryFsi
              && movingDiagnostics.geometryAdvanceCount == 1
              && movingDiagnostics.maximumGeometryDisplacementMeters > 0.0,
          "moving scene FSI accepts one finite XPBD geometry advance");
    check(movingDiagnostics.usesConsecutivePressureWarmStart
              && movingDiagnostics.usesRegionWallPrediction,
          "moving scene FSI accepts the rebased consecutive pressure endpoint");
    check(movingDiagnostics.wallMomentumResidualKilogramMetersPerSecond
                  < 1.0e-12
              && movingDiagnostics.totalFluidTransferForceResidualNewtons
                  < 1.0e-8,
          "moving scene FSI closes wall momentum and total conservative transfer");
    check(vectorField(movingFrame, "moving_scene.nodal_wall_force")
                  != nullptr
              && vectorField(
                     movingFrame,
                     "moving_scene.nodal_total_fluid_force")
                  != nullptr,
          "moving scene FSI publishes wall and total nodal load fields");
}

} // namespace

int main() {
    try {
        testFrozenScenePressureCase();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "FAIL: %s\n", exception.what());
        ++failures;
    }
    if (failures == 0) {
        std::printf("simwing frozen scene pressure case tests passed\n");
    }
    return failures == 0 ? 0 : 1;
}
