#include "open_piston_case.h"
#include "viewer_protocol.h"

#include <cmath>
#include <cstdio>
#include <ranges>
#include <string_view>
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

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g, tolerance %.3g)\n",
                     message, actual, expected, tolerance);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "open piston: deterministic frame serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const std::string_view name) {
    const auto found = std::ranges::find(
        frame.scalarFields, name, &viewer::ScalarField::name);
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

void testVisibleOpenPistonAndDeterminism() {
    fsi::OpenPistonCase first;
    fsi::OpenPistonCase second;
    const auto firstInitialFrame = first.advance();
    const auto secondInitialFrame = second.advance();
    check(serialized(firstInitialFrame) == serialized(secondInitialFrame),
          "open piston: first accelerated frame replays bit-for-bit");

    const auto* initialPressure = scalarField(
        firstInitialFrame, "interface.pressure_traction");
    const auto* initialActuator = scalarField(
        firstInitialFrame, "actuator.step_impulse");
    const auto* initialGcl = scalarField(
        firstInitialFrame, "fluid.gcl_residual");
    check(initialPressure != nullptr
              && initialPressure->values.size() == 2
              && initialPressure->values.front() < 0.0
              && initialPressure->values[0] == initialPressure->values[1],
          "open piston: first frame exposes uniform resisting CFD pressure");
    check(initialActuator != nullptr
              && initialActuator->values.size() == 1
              && initialActuator->values.front() > 300.0,
          "open piston: explicit actuator overcomes inertia and fluid reaction");
    check(initialGcl != nullptr
              && std::abs(initialGcl->values.front()) < 1.0e-13,
          "open piston: first visible frame closes geometric conservation");
    const auto& initialBridge = first.bridgeDiagnostics();
    const auto& initialCutSurface = first.cutSurfaceDiagnostics();
    const double firstEndPlane = 3.0 + 0.05
        * first.stepSettings().timeStepSeconds;
    check(initialBridge.correspondenceMode
              == fsi::PlanarFaceCorrespondenceMode::RigidNormalTranslation
              && initialBridge.overlapPatchCount >= 6
              && initialBridge.gridPlaneCoordinateMeters == 3.0,
          "open piston: pressure uses moving face-resolved correspondence");
    checkNear(initialBridge.physicalPlaneCoordinateMeters,
              firstEndPlane, 0.0,
              "open piston: physical pressure plane follows the plate");
    check(initialBridge.maximumRigidPositionResidualMeters < 1.0e-14
              && initialBridge
                     .maximumRigidVelocityResidualMetersPerSecond < 1.0e-14
              && initialBridge.forceResidualNormNewtons < 1.0e-10
              && initialBridge.momentResidualNormNewtonMeters < 1.0e-10
              && std::abs(initialBridge.powerResidualWatts) < 1.0e-10,
          "open piston: moving correspondence ledgers meet their budgets");
    check(initialCutSurface.accepted
              && initialCutSurface.faceCount == 6
              && initialCutSurface.gridPlaneCoordinateMeters == 3.0
              && initialCutSurface.physicalPlaneCoordinateMeters
                  == firstEndPlane
              && initialCutSurface.periodicPositionResidualMeters < 1.0e-14
              && initialCutSurface.forceResidualNormNewtons < 1.0e-10
              && std::abs(initialCutSurface.powerResidualWatts) < 1.0e-10,
          "open piston: fluid accepts pressure geometry at the physical cut plane");
    const auto& initialConservation = first.conservationDiagnostics();
    check(initialConservation.accepted
              && initialConservation.version
                  == fsi::openPistonConservationVersion
              && initialConservation
                     .structureMomentumResidualNormNewtonSeconds < 1.0e-10
              && initialConservation
                     .fluidMomentumResidualNormNewtonSeconds < 1.0e-10
              && initialConservation
                     .systemMomentumResidualNormNewtonSeconds < 1.0e-10
              && std::abs(initialConservation.structureEnergyResidualJoules)
                  < 1.0e-10
              && std::abs(initialConservation.fluidEnergyResidualJoules)
                  < 1.0e-10
              && std::abs(initialConservation.systemEnergyResidualJoules)
                  < 1.0e-10,
          "open piston: first coupled momentum and energy ledgers close");
    checkNear(initialConservation.fluidMomentumChangeNewtonSeconds.x,
              1.44, 3.0e-12,
              "open piston: fluid plug gains its analytic momentum");
    checkNear(
        initialConservation.pressureImpulseToStructureNewtonSeconds.x,
        -1.44, 3.0e-12,
        "open piston: complete CFD reaction delivers opposite impulse");
    checkNear(initialConservation.actuatorImpulseNewtonSeconds.x,
              301.44, 3.0e-12,
              "open piston: actuator supplies structure plus fluid momentum");
    checkNear(initialConservation.fluidKineticEnergyChangeJoules,
              0.036, 2.0e-14,
              "open piston: fluid plug gains its analytic kinetic energy");
    checkNear(initialConservation.pressureWorkToStructureJoules,
              -0.036, 2.0e-14,
              "open piston: time-averaged CFD reaction loses fluid kinetic work");
    checkNear(initialConservation.actuatorWorkJoules,
              7.536, 2.0e-13,
              "open piston: actuator work balances both kinetic-energy gains");
    const auto* initialGridPlane = scalarField(
        firstInitialFrame, "interface.grid_plane");
    const auto* initialPhysicalPlane = scalarField(
        firstInitialFrame, "interface.physical_plane");
    const auto* initialCorrespondenceResidual = scalarField(
        firstInitialFrame, "interface.correspondence_residual");
    const auto* initialCorrespondenceVelocityResidual = scalarField(
        firstInitialFrame, "interface.correspondence_velocity_residual");
    const auto* initialCutPositionResidual = scalarField(
        firstInitialFrame, "fluid.cut_surface_periodic_residual");
    const auto* initialCutForceResidual = scalarField(
        firstInitialFrame, "fluid.cut_surface_force_residual");
    const auto* initialCutPowerResidual = scalarField(
        firstInitialFrame, "fluid.cut_surface_power_residual");
    const auto* initialSystemMomentumResidual = scalarField(
        firstInitialFrame, "conservation.system_momentum_residual");
    const auto* initialSystemEnergyResidual = scalarField(
        firstInitialFrame, "conservation.system_energy_residual");
    check(initialGridPlane != nullptr
              && initialGridPlane->values.front() == 3.0
              && initialPhysicalPlane != nullptr
              && initialPhysicalPlane->values.front() == firstEndPlane
              && initialCorrespondenceResidual != nullptr
              && initialCorrespondenceResidual->values.front() < 1.0e-14
              && initialCorrespondenceVelocityResidual != nullptr
              && initialCorrespondenceVelocityResidual->values.front()
                  < 1.0e-14
              && initialCutPositionResidual != nullptr
              && initialCutPositionResidual->values.front() < 1.0e-14
              && initialCutForceResidual != nullptr
              && initialCutForceResidual->values.front() < 1.0e-10
              && initialCutPowerResidual != nullptr
              && std::abs(initialCutPowerResidual->values.front()) < 1.0e-10
              && initialSystemMomentumResidual != nullptr
              && initialSystemMomentumResidual->values.front() < 1.0e-10
              && initialSystemEnergyResidual != nullptr
              && std::abs(initialSystemEnergyResidual->values.front())
                  < 1.0e-10,
          "open piston: trace exposes grid/physical correspondence geometry");

    viewer::DiagnosticFrame finalFrame = firstInitialFrame;
    constexpr std::uint64_t steps = 24;
    for (std::uint64_t step = 1; step < steps; ++step) {
        finalFrame = first.advance();
        const auto replay = second.advance();
        check(serialized(finalFrame) == serialized(replay),
              "open piston: continued accepted frames replay bit-for-bit");
    }

    const auto checkpoint = first.structure().checkpoint();
    check(checkpoint.acceptedStepCount == steps
              && finalFrame.step == steps
              && finalFrame.sceneChecksum == fsi::openPistonCaseChecksum
              && finalFrame.solverCommit == fsi::openPistonCaseSolverId,
          "open piston: trace provenance and accepted step remain aligned");
    check(finalFrame.triangles.size() == 2
              && finalFrame.triangles[0].negativeRegionId == 9
              && finalFrame.triangles[0].positiveRegionId == 9,
          "open piston: viewer exposes a nonseparating two-sided sheet");
    const double expectedOffset = 0.05
        * first.stepSettings().timeStepSeconds
        * static_cast<double>(steps);
    checkNear(first.surfaceOffsetMeters(), expectedOffset, 2.0e-17,
              "open piston: accepted plate displacement drives cut-cell offset");
    for (const auto& node : checkpoint.nodes) {
        checkNear(node.positionMeters.x, 3.0 + expectedOffset, 2.0e-14,
                  "open piston: prescribed actuator leaves visible translation");
        checkNear(node.velocityMetersPerSecond.x, 0.05, 2.0e-14,
                  "open piston: plate coasts at the prescribed verification speed");
    }
    checkNear(first.structure().diagnostics()
                  .linearMomentumKgMetersPerSecond.x,
              300.0, 2.0e-10,
              "open piston: structural momentum matches the driven plate");

    const auto& control = first.controlVolumeDiagnostics();
    check(control.accepted && control.finite,
          "open piston: final control-volume ledger remains accepted");
    checkNear(control.endVolumeCubicMeters,
              12.0 + 6.0 * expectedOffset, 2.0e-14,
              "open piston: visible motion accumulates analytic chamber volume");
    check(std::abs(control.surfaceGeometryResidualCubicMeters) < 2.0e-15
              && std::abs(control.continuityResidualCubicMeters) < 2.0e-15,
          "open piston: geometry sweep and opening transport remain closed");
    check(finalFrame.couplingResiduals.fluid < 2.0e-11
              && finalFrame.couplingResiduals.tractionNewtons < 1.0e-10
              && finalFrame.couplingResiduals.interfacePowerWatts < 1.0e-10,
          "open piston: visible numerical residuals meet their budgets");
}

void testAcceptedTopologyRebase() {
    fsi::OpenPistonCase first;
    fsi::OpenPistonCase second;
    viewer::DiagnosticFrame crossingFrame;
    constexpr std::uint64_t stepsToFirstRebase = 1200;
    for (std::uint64_t step = 0; step < stepsToFirstRebase; ++step) {
        crossingFrame = first.advance();
        const auto replay = second.advance();
        if (step + 2 >= stepsToFirstRebase) {
            check(serialized(crossingFrame) == serialized(replay),
                  "rebase: frames around the topology crossing replay bit-for-bit");
        }
    }

    check(first.topologyRebaseCount() == 1
              && second.topologyRebaseCount() == 1
              && first.movingPlaneCoordinate() == 7
              && first.surfaceOffsetMeters() == 0.0,
          "rebase: the accepted worker advances exactly one MAC plane");
    const auto& control = first.controlVolumeDiagnostics();
    const auto& rebase = first.lastRebaseDiagnostics();
    check(control.accepted
              && control.movingPlaneCoordinate == 6
              && control.endCutCellVolumeFraction == 1.0
              && rebase.accepted
              && rebase.previousMovingPlaneCoordinate == 6
              && rebase.rebasedMovingPlaneCoordinate == 7,
          "rebase: terminal and candidate epochs remain explicit");
    checkNear(control.endVolumeCubicMeters, 15.0, 2.0e-13,
              "rebase: the terminal old epoch reaches five full layers");
    checkNear(rebase.rebasedReferenceVolumeCubicMeters, 15.0, 0.0,
              "rebase: the new epoch starts at the same chamber volume");
    checkNear(rebase.volumeContinuityResidualCubicMeters, 0.0, 0.0,
              "rebase: topology transition preserves volume exactly");
    check(first.lastRebaseVelocityResidualMetersPerSecond() < 2.0e-12,
          "rebase: constraint remap preserves material velocity within budget");
    for (const auto& node : first.structure().checkpoint().nodes) {
        checkNear(node.positionMeters.x, 3.5, 2.0e-13,
                  "rebase: structural plate reaches the crossed grid face");
    }

    const auto* crossingCount = scalarField(
        crossingFrame, "fluid.topology_rebase_count");
    const auto* crossingVolume = scalarField(
        crossingFrame, "fluid.chamber_volume");
    const auto* crossingFraction = scalarField(
        crossingFrame, "fluid.cut_cell_fraction");
    const auto* crossingVolumeResidual = scalarField(
        crossingFrame, "fluid.rebase_volume_residual");
    const auto* crossingVelocityResidual = scalarField(
        crossingFrame, "fluid.rebase_velocity_residual");
    check(crossingCount != nullptr
              && crossingCount->values.front() == 1.0
              && crossingVolume != nullptr
              && crossingVolume->values.front() == 15.0
              && crossingFraction != nullptr
              && crossingFraction->values.front() == 1.0,
          "rebase: crossing frame publishes its epoch and terminal geometry");
    check(crossingVolumeResidual != nullptr
              && crossingVolumeResidual->values.front() == 0.0
              && crossingVelocityResidual != nullptr
              && crossingVelocityResidual->values.front() < 2.0e-12,
          "rebase: crossing frame exposes volume and velocity residuals");

    const auto continuedFrame = first.advance();
    const auto replayContinuedFrame = second.advance();
    check(serialized(continuedFrame) == serialized(replayContinuedFrame),
          "rebase: the first partial-cell step in the new epoch replays");
    const double nextOffset = 0.05
        * first.stepSettings().timeStepSeconds;
    check(first.topologyRebaseCount() == 1
              && first.movingPlaneCoordinate() == 7,
          "rebase: ordinary continuation does not create another epoch");
    checkNear(first.surfaceOffsetMeters(), nextOffset, 2.0e-17,
              "rebase: partial-cell offset restarts from zero");
    checkNear(first.controlVolumeDiagnostics().startVolumeCubicMeters,
              15.0, 0.0,
              "rebase: continued geometry starts at the terminal old volume");
    checkNear(first.controlVolumeDiagnostics().endVolumeCubicMeters,
              15.0 + 6.0 * nextOffset, 2.0e-13,
              "rebase: chamber growth continues in the new epoch");
    const auto& continuedBridge = first.bridgeDiagnostics();
    check(continuedBridge.correspondenceMode
              == fsi::PlanarFaceCorrespondenceMode::RigidNormalTranslation
              && continuedBridge.gridPlaneCoordinateMeters == 3.5
              && continuedBridge.overlapPatchCount >= 6,
          "rebase: face-resolved material patches survive the new grid epoch");
    checkNear(continuedBridge.physicalPlaneCoordinateMeters,
              3.5 + nextOffset, 2.0e-13,
              "rebase: physical correspondence continues beyond the grid plane");
    check(continuedBridge.maximumRigidPositionResidualMeters < 1.0e-13
              && continuedBridge
                     .maximumRigidVelocityResidualMetersPerSecond < 1.0e-13
              && continuedBridge.forceResidualNormNewtons < 1.0e-10
              && continuedBridge.momentResidualNormNewtonMeters < 1.0e-10
              && std::abs(continuedBridge.powerResidualWatts) < 1.0e-10,
          "rebase: continued face-resolved ledgers remain closed");
    check(first.cutSurfaceDiagnostics().accepted
              && first.cutSurfaceDiagnostics().gridPlaneCoordinateMeters
                  == 3.5
              && first.cutSurfaceDiagnostics()
                     .periodicPositionResidualMeters < 1.0e-13,
          "rebase: fluid-side cut pressure follows the new physical epoch");
    check(first.conservationDiagnostics().accepted
              && first.conservationDiagnostics()
                     .systemMomentumResidualNormNewtonSeconds < 1.0e-8
              && std::abs(first.conservationDiagnostics()
                              .systemEnergyResidualJoules) < 2.0e-9,
          "rebase: continued coupled conservation remains accepted");
    checkNear(first.cutSurfaceDiagnostics().physicalPlaneCoordinateMeters,
              3.5 + nextOffset, 2.0e-13,
              "rebase: fluid-side cut pressure uses the physical plate position");

    viewer::DiagnosticFrame wrappedCrossingFrame = continuedFrame;
    while (first.structure().checkpoint().acceptedStepCount < 2400) {
        wrappedCrossingFrame = first.advance();
        const auto replay = second.advance();
        if (first.structure().checkpoint().acceptedStepCount >= 2399) {
            check(serialized(wrappedCrossingFrame) == serialized(replay),
                  "rebase: periodic crossing frames replay bit-for-bit");
        }
    }
    check(first.topologyRebaseCount() == 2
              && first.movingPlaneCoordinate() == 0
              && first.surfaceOffsetMeters() == 0.0,
          "rebase: second crossing wraps to periodic grid plane zero");
    check(first.lastRebaseDiagnostics().accepted
              && first.lastRebaseDiagnostics()
                     .previousMovingPlaneCoordinate == 7
              && first.lastRebaseDiagnostics()
                     .rebasedMovingPlaneCoordinate == 0,
          "rebase: periodic old and new topology epochs remain explicit");
    checkNear(first.controlVolumeDiagnostics().endVolumeCubicMeters,
              18.0, 3.0e-13,
              "rebase: periodic crossing preserves the analytic chamber volume");
    checkNear(first.bridgeDiagnostics().gridPlaneCoordinateMeters,
              3.5, 0.0,
              "rebase: terminal transfer remains on the old Eulerian plane");
    checkNear(first.bridgeDiagnostics().physicalPlaneCoordinateMeters,
              4.0, 3.0e-13,
              "rebase: terminal physical plane remains unwrapped");
    check(first.cutSurfaceDiagnostics().accepted
              && first.cutSurfaceDiagnostics().gridPlaneCoordinateMeters
                  == 3.5
              && first.cutSurfaceDiagnostics()
                     .periodicPositionResidualMeters < 3.0e-13,
          "rebase: terminal fluid reaction is placed on the unwrapped cut plane");
    checkNear(first.cutSurfaceDiagnostics().physicalPlaneCoordinateMeters,
              4.0, 3.0e-13,
              "rebase: terminal cut pressure position remains unwrapped");

    const auto wrappedContinued = first.advance();
    const auto wrappedReplay = second.advance();
    check(serialized(wrappedContinued) == serialized(wrappedReplay),
          "rebase: first step after periodic wrap replays bit-for-bit");
    checkNear(first.bridgeDiagnostics().gridPlaneCoordinateMeters,
              0.0, 0.0,
              "rebase: continued transfer uses wrapped Eulerian plane zero");
    checkNear(first.bridgeDiagnostics().physicalPlaneCoordinateMeters,
              4.0 + nextOffset, 3.0e-13,
              "rebase: continued transfer retains unwrapped physical position");
    check(first.bridgeDiagnostics().maximumRigidPositionResidualMeters
              < 2.0e-13
              && first.bridgeDiagnostics()
                     .maximumRigidVelocityResidualMetersPerSecond < 1.0e-13,
          "rebase: wrapped moving correspondence remains kinematically bound");
    check(first.cutSurfaceDiagnostics().accepted
              && first.cutSurfaceDiagnostics().gridPlaneCoordinateMeters
                  == 0.0
              && first.cutSurfaceDiagnostics()
                     .periodicPositionResidualMeters < 3.0e-13,
          "rebase: cut pressure accepts the periodic physical image");
    check(first.conservationDiagnostics().accepted
              && first.conservationDiagnostics()
                     .systemMomentumResidualNormNewtonSeconds < 1.0e-8
              && std::abs(first.conservationDiagnostics()
                              .systemEnergyResidualJoules) < 2.0e-9,
          "rebase: wrapped coupled conservation remains accepted");
    checkNear(first.cutSurfaceDiagnostics().physicalPlaneCoordinateMeters,
              4.0 + nextOffset, 3.0e-13,
              "rebase: wrapped cut pressure retains the unwrapped position");
}

} // namespace

int main() {
    testVisibleOpenPistonAndDeterminism();
    testAcceptedTopologyRebase();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing open piston case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing open piston case checks passed");
    return 0;
}
