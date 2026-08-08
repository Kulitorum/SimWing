#include "fluid/moving_control_volume.h"
#include "fluid/planar_cut_surface.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::GridFacePressureJump;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::MovingPorousProjectionSettings;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PlanarControlVolumeStep;
using simwing::fsi::fluid::PlanarMovingControlVolume;
using simwing::fsi::fluid::PorousGridFaceCrossing;
using simwing::fsi::fluid::SharpPressureJumpField;
using simwing::fsi::fluid::evaluatePlanarCutSurfacePressure;
using simwing::fsi::fluid::projectVelocityWithMovingAndPorousInterfaces;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;
using simwing::fsi::fluid::rebasePlanarMovingControlVolume;

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

double maximumVelocityDifference(
    const MacVelocityField& first,
    const MacVelocityField& second) {
    double result = 0.0;
    const auto accumulate = [&](const auto firstValues,
                                const auto secondValues) {
        for (std::size_t index = 0; index < firstValues.size(); ++index) {
            result = std::max(
                result, std::abs(firstValues[index] - secondValues[index]));
        }
    };
    accumulate(first.xFaces(), second.xFaces());
    accumulate(first.yFaces(), second.yFaces());
    accumulate(first.zFaces(), second.zFaces());
    return result;
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

PeriodicCartesianGrid pistonGrid() {
    return PeriodicCartesianGrid({8, 2, 3}, {}, {4.0, 2.0, 3.0});
}

std::vector<GridFaceMovingInterface> completePlaneFaces(
    const PeriodicCartesianGrid& grid,
    const GridFaceAxis axis,
    const std::size_t planeCoordinate,
    const std::uint64_t surfaceStableId,
    const std::uint64_t regionStableId,
    const double speedMetersPerSecond) {
    std::vector<GridFaceMovingInterface> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t coordinate = axis == GridFaceAxis::X
                    ? i : (axis == GridFaceAxis::Y ? j : k);
                if (coordinate == planeCoordinate) {
                    result.push_back({
                        surfaceStableId,
                        regionStableId,
                        regionStableId,
                        axis,
                        i,
                        j,
                        k,
                        speedMetersPerSecond,
                    });
                }
            }
        }
    }
    return result;
}

std::vector<GridFaceMovingInterface> pistonFaces(
    const PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond = 0.25,
    const std::size_t planeCoordinate = 6,
    const std::uint64_t regionStableId = 9) {
    return completePlaneFaces(
        grid, GridFaceAxis::X, planeCoordinate, 300,
        regionStableId, speedMetersPerSecond);
}

MovingInterfaceProjectionSettings projectionSettings() {
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.4;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return settings;
}

std::vector<PorousGridFaceCrossing> porousOpeningFaces(
    const PeriodicCartesianGrid& grid) {
    std::vector<PorousGridFaceCrossing> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                600, 10, 11, GridFaceAxis::X, 2, j, k,
                0.35, 0.0, {10.0, 0.0}});
        }
    }
    return result;
}

SharpPressureJumpField porousOpeningBalance(
    const PeriodicCartesianGrid& grid) {
    std::vector<GridFacePressureJump> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                601, 11, 10, GridFaceAxis::X, 4, j, k,
                2.5, 0.65});
        }
    }
    return SharpPressureJumpField(grid, std::move(result));
}

MovingPorousProjectionSettings porousControlVolumeSettings() {
    MovingPorousProjectionSettings result;
    result.movingProjection = projectionSettings();
    result.iteration.absoluteNormalVelocityToleranceMetersPerSecond =
        1.0e-12;
    result.iteration.relativeNormalVelocityTolerance = 1.0e-12;
    result.iteration.absolutePressureJumpTolerancePascals = 1.0e-11;
    result.iteration.relativePressureJumpTolerance = 1.0e-12;
    result.iteration.relaxation = 0.5;
    result.iteration.maximumNonlinearIterations = 100;
    return result;
}

void testAcceleratedOpenPistonAndGcl() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(
        grid, pistonFaces(grid));
    const PlanarMovingControlVolume controlVolume(
        grid, interfaces, 300, 2);
    MacVelocityField firstVelocity(grid);
    MacVelocityField secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure = firstPressure;
    const auto settings = projectionSettings();
    const auto firstFluid = projectVelocityWithMovingInterfaces(
        grid, firstVelocity, firstPressure, interfaces, settings);
    const auto secondFluid = projectVelocityWithMovingInterfaces(
        grid, secondVelocity, secondPressure, interfaces, settings);

    check(firstFluid.projection.converged
              && firstFluid.projection.iterationCount > 0,
          "open piston: accelerating the connected fluid requires pressure projection");
    check(firstVelocity == secondVelocity
              && firstPressure == secondPressure
              && firstFluid == secondFluid,
          "open piston: pressure, velocity, and diagnostics replay bit-for-bit");
    check(firstFluid.fluidRegionCount == 1
              && firstFluid.regions.size() == 1
              && firstFluid.regions.front().stableId == 9,
          "open piston: both membrane sides remain one connected fluid region");
    checkNear(
        firstFluid.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond,
        0.0, 2.0e-15,
        "open piston: the connected periodic region remains globally compatible");
    check(firstFluid.maximumNormalVelocityErrorMetersPerSecond == 0.0,
          "open piston: moving-surface velocity remains bit-exact");
    check(firstFluid.projection.divergenceL2AfterPerSecond < 2.0e-11,
          "open piston: projected opening flow is discretely divergence-free");
    for (const double velocity : firstVelocity.xFaces()) {
        checkNear(velocity, 0.25, 2.0e-12,
                  "open piston: projection routes plug flow around the opening");
    }
    check(firstFluid.surfaces.size() == 1
              && firstFluid.surfaces.front().stableId == 300
              && firstFluid.surfaces.front().pressureForceNewtons.x < 0.0
              && firstFluid.surfaces.front().pressurePowerWatts < 0.0,
          "open piston: accelerating fluid exerts a resisting pressure load");
    checkNear(firstFluid.surfaces.front()
                  .constraintReactionImpulseNewtonSeconds.x,
              -7.2, 3.0e-12,
              "open piston: complete constraint reaction balances fluid momentum");
    checkNear(firstFluid.surfaces.front()
                  .constraintReactionWorkJoules,
              -1.8, 1.0e-12,
              "open piston: end-velocity reaction power retains its rectangular work");
    check(firstFluid.surfaces.front().directConstraintForceNewtons.x < 0.0
              && firstFluid.surfaces.front()
                     .constraintReactionForceNewtons.x
                  < firstFluid.surfaces.front().pressureForceNewtons.x,
          "open piston: direct MAC enforcement completes adjacent pressure traction");

    const PlanarControlVolumeStep step{0.0, 0.1, 0.4};
    const auto first = controlVolume.evaluate(
        grid, firstVelocity, firstFluid, step);
    const auto second = controlVolume.evaluate(
        grid, secondVelocity, secondFluid, step);
    check(first == second,
          "open piston: control-volume ledgers replay bit-for-bit");
    check(first.version
              == simwing::fsi::fluid::planarMovingControlVolumeVersion
              && first.movingSurfaceStableId == 300
              && first.fluidRegionStableId == 9
              && first.movingSurfaceFaceCount == 6
              && first.openingFaceCount == 6,
          "open piston: stable surface, region, and complete planes are explicit");
    checkNear(controlVolume.crossSectionAreaSquareMeters(), 6.0, 0.0,
              "open piston: cross-section area is exact");
    checkNear(controlVolume.referenceVolumeCubicMeters(), 12.0, 0.0,
              "open piston: four full chamber layers have analytic volume");
    checkNear(first.startVolumeCubicMeters, 12.0, 0.0,
              "open piston: initial chamber volume matches reference geometry");
    checkNear(first.endVolumeCubicMeters, 12.6, 2.0e-15,
              "open piston: moving partial cell increases chamber volume");
    checkNear(first.endCutCellVolumeCubicMeters, 0.6, 2.0e-16,
              "open piston: partial-cell volume is area times offset");
    checkNear(first.endCutCellVolumeFraction, 0.2, 0.0,
              "open piston: partial-cell fraction uses normal grid spacing");
    checkNear(first.geometryVolumeChangeCubicMeters, 0.6, 2.0e-15,
              "open piston: geometry ledger reports analytic swept volume");
    checkNear(first.surfaceSweptVolumeCubicMeters, 0.6, 2.0e-15,
              "open piston: moving-face sweep independently recovers volume");
    checkNear(first.openingTransportVolumeCubicMeters, 0.6, 3.0e-12,
              "open piston: projected opening flux fills the swept chamber");
    check(std::abs(first.surfaceGeometryResidualCubicMeters) < 3.0e-15
              && std::abs(first.continuityResidualCubicMeters) < 3.0e-12
              && first.maximumSurfaceVelocityErrorMetersPerSecond < 2.0e-12
              && first.finite && first.accepted,
          "open piston: surface geometry and opening-flux GCL ledgers close");
}

void testPorousOpeningAcrossControlAndCutSurface() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(
        grid, pistonFaces(grid));
    const PlanarMovingControlVolume controlVolume(
        grid, interfaces, 300, 2);
    const auto porous = porousOpeningFaces(grid);
    const auto balance = porousOpeningBalance(grid);
    const auto settings = porousControlVolumeSettings();
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto fluid = projectVelocityWithMovingAndPorousInterfaces(
        grid, velocity, pressure, interfaces,
        porous, balance, settings);
    const PlanarControlVolumeStep step{0.0, 0.1, 0.4};
    const auto control = controlVolume.evaluate(
        grid, velocity, fluid, step);
    const auto cut = evaluatePlanarCutSurfacePressure(
        grid, controlVolume, fluid, 0.1, 3.1);

    double porousTransportCubicMeters = 0.0;
    for (const auto& sample : fluid.porous.samples) {
        porousTransportCubicMeters +=
            sample.volumeFlowRateCubicMetersPerSecond
            * step.durationSeconds;
    }
    check(fluid.accepted && control.accepted && cut.accepted,
          "porous control volume: projection, GCL, and cut reaction all accept");
    checkNear(porousTransportCubicMeters,
              control.geometryVolumeChangeCubicMeters, 8.0e-12,
              "porous control volume: resolved porous opening flux fills the swept chamber");
    checkNear(control.openingTransportVolumeCubicMeters,
              porousTransportCubicMeters, 8.0e-12,
              "porous control volume: GCL opening transport is the porous tile sum");
    checkNear(fluid.porous.samples.front()
                  .pressureJump.pressureJumpPascals,
              -2.5, 3.0e-11,
              "porous control volume: opening slip closes the analytic Darcy loss");
    checkNear(fluid.porous.totalPorousDissipationJoules,
              1.5, 3.0e-10,
              "porous control volume: opening dissipation integrates over the epoch");
    check(cut.pressureForceNewtons
              == fluid.movingInterface.surfaces.front()
                     .constraintReactionForceNewtons
              && cut.pressurePowerWatts
                  == fluid.movingInterface.surfaces.front()
                         .constraintReactionPowerWatts,
          "porous control volume: physical cut geometry retains the final moving reaction");

    auto rejected = fluid;
    rejected.accepted = false;
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, rejected, step)); },
        "porous control volume: an unaccepted outer iteration cannot expose its inner GCL");
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, rejected, 0.1, 3.1)); },
        "porous control volume: an unaccepted outer iteration cannot expose a cut reaction");

    auto inconsistent = fluid;
    ++inconsistent.porous.projection.iterationCount;
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, inconsistent, step)); },
        "porous control volume: divergent nested projection diagnostics are rejected");
}

void testAllAxisControlVolumes() {
    const PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    for (std::size_t axisIndex = 0; axisIndex < 3; ++axisIndex) {
        const auto axis = static_cast<GridFaceAxis>(axisIndex);
        const FaceAlignedMovingInterface interfaces(
            grid, completePlaneFaces(
                grid, axis, 3, 400 + axisIndex,
                20 + axisIndex, 0.125));
        const PlanarMovingControlVolume controlVolume(
            grid, interfaces, 400 + axisIndex, 1);
        MacVelocityField velocity(grid);
        CellScalarField pressure(grid);
        auto settings = projectionSettings();
        settings.projection.timeStepSeconds = 0.1;
        const auto fluid = projectVelocityWithMovingInterfaces(
            grid, velocity, pressure, interfaces, settings);
        const auto control = controlVolume.evaluate(
            grid, velocity, fluid, {0.0, 0.0125, 0.1});
        check(fluid.projection.converged && control.accepted,
              "axes: every normal component accepts the open piston projection");
        check(control.axis == axis
                  && control.movingSurfaceFaceCount == 16
                  && control.openingFaceCount == 16,
              "axes: complete moving and opening planes retain their orientation");
        checkNear(control.crossSectionAreaSquareMeters, 16.0, 0.0,
                  "axes: cross-section area is orientation independent");
        checkNear(control.referenceVolumeCubicMeters, 32.0, 0.0,
                  "axes: two-layer reference volume is orientation independent");
        checkNear(control.geometryVolumeChangeCubicMeters, 0.2, 4.0e-15,
                  "axes: partial-cell geometry has analytic volume change");
        checkNear(control.openingTransportVolumeCubicMeters, 0.2, 2.0e-12,
                  "axes: projected opening transport matches every orientation");

        const auto terminal = controlVolume.evaluate(
            grid, velocity, fluid, {0.9, 1.0, 0.8, true});
        const FaceAlignedMovingInterface rebasedInterfaces(
            grid, completePlaneFaces(
                grid, axis, 0, 400 + axisIndex,
                20 + axisIndex, 0.125));
        const auto rebased = rebasePlanarMovingControlVolume(
            grid, controlVolume, rebasedInterfaces, terminal);
        check(terminal.accepted && rebased.diagnostics.accepted
                  && rebased.controlVolume.movingPlaneCoordinate() == 0,
              "axes: terminal steps rebase across the periodic boundary");
        checkNear(rebased.diagnostics.previousTerminalVolumeCubicMeters,
                  48.0, 0.0,
                  "axes: terminal volume includes the completed partial cell");
        checkNear(rebased.controlVolume.referenceVolumeCubicMeters(),
                  48.0, 0.0,
                  "axes: wrapped reference volume preserves continuity");
        checkNear(rebased.diagnostics.volumeContinuityResidualCubicMeters,
                  0.0, 0.0,
                  "axes: periodic topology rebase has exact volume continuity");
    }
}

void testExplicitTopologyRebase() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(
        grid, pistonFaces(grid));
    const PlanarMovingControlVolume controlVolume(
        grid, interfaces, 300, 2);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto fluid = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());
    const auto terminal = controlVolume.evaluate(
        grid, velocity, fluid, {0.4, 0.5, 0.4, true});
    check(terminal.accepted
              && terminal.endCutCellVolumeFraction == 1.0
              && terminal.movingPlaneCoordinate == 6
              && terminal.openingPlaneCoordinate == 2,
          "rebase: an explicitly terminal old epoch reaches one full cell");

    const FaceAlignedMovingInterface rebasedInterfaces(
        grid, pistonFaces(grid, 0.25, 7));
    const auto first = rebasePlanarMovingControlVolume(
        grid, controlVolume, rebasedInterfaces, terminal);
    const auto second = rebasePlanarMovingControlVolume(
        grid, controlVolume, rebasedInterfaces, terminal);
    check(first.diagnostics == second.diagnostics
              && first.diagnostics.accepted
              && first.diagnostics.previousMovingPlaneCoordinate == 6
              && first.diagnostics.rebasedMovingPlaneCoordinate == 7,
          "rebase: candidate topology and ledger replay deterministically");
    checkNear(first.diagnostics.completedCellOffsetMeters, 0.5, 0.0,
              "rebase: completed offset equals one normal cell spacing");
    checkNear(first.diagnostics.previousTerminalVolumeCubicMeters, 15.0, 0.0,
              "rebase: old terminal chamber has five full layers");
    checkNear(first.diagnostics.rebasedReferenceVolumeCubicMeters, 15.0, 0.0,
              "rebase: next topology converts the partial cell to a full layer");
    checkNear(first.diagnostics.volumeContinuityResidualCubicMeters, 0.0, 0.0,
              "rebase: old and new topology volumes meet exactly");

    auto rebasedVelocity = velocity;
    auto rebasedPressure = pressure;
    const auto rebasedFluid = projectVelocityWithMovingInterfaces(
        grid, rebasedVelocity, rebasedPressure, rebasedInterfaces,
        projectionSettings());
    const auto continued = first.controlVolume.evaluate(
        grid, rebasedVelocity, rebasedFluid, {0.0, 0.1, 0.4});
    check(maximumVelocityDifference(rebasedVelocity, velocity) < 2.0e-12
              && continued.accepted,
          "rebase: compatible plug flow remaps within projection tolerance");
    checkNear(continued.startVolumeCubicMeters, 15.0, 0.0,
              "rebase: the next epoch starts at the old terminal volume");
    checkNear(continued.endVolumeCubicMeters, 15.6, 4.0e-15,
              "rebase: partial-cell growth continues after the crossing");

    auto corruptedTerminal = terminal;
    corruptedTerminal.endVolumeCubicMeters += 1.0;
    const auto corrupted = rebasePlanarMovingControlVolume(
        grid, controlVolume, rebasedInterfaces, corruptedTerminal);
    check(!corrupted.diagnostics.accepted
              && std::abs(corrupted.diagnostics
                              .volumeContinuityResidualCubicMeters) > 0.9,
          "rebase: a broken terminal volume ledger is not accepted");
    auto corruptedArea = terminal;
    corruptedArea.crossSectionAreaSquareMeters += 1.0;
    const auto areaMismatch = rebasePlanarMovingControlVolume(
        grid, controlVolume, rebasedInterfaces, corruptedArea);
    check(!areaMismatch.diagnostics.accepted,
          "rebase: a broken terminal area ledger is not accepted");

    const auto nonterminal = controlVolume.evaluate(
        grid, velocity, fluid, {0.0, 0.1, 0.4});
    const auto early = rebasePlanarMovingControlVolume(
        grid, controlVolume, rebasedInterfaces, nonterminal);
    check(!early.diagnostics.accepted,
          "rebase: an incomplete partial cell cannot advance topology");

    const FaceAlignedMovingInterface skippedInterfaces(
        grid, pistonFaces(grid, 0.25, 0));
    expectRejected(
        [&] { static_cast<void>(rebasePlanarMovingControlVolume(
            grid, controlVolume, skippedInterfaces, terminal)); },
        "rebase: candidate topology cannot skip a MAC plane");
    const FaceAlignedMovingInterface changedRegionInterfaces(
        grid, pistonFaces(grid, 0.25, 7, 10));
    expectRejected(
        [&] { static_cast<void>(rebasePlanarMovingControlVolume(
            grid, controlVolume, changedRegionInterfaces, terminal)); },
        "rebase: candidate topology cannot change stable fluid identity");

    const FaceAlignedMovingInterface finalEpochInterfaces(
        grid, pistonFaces(grid, 0.25, 1));
    const PlanarMovingControlVolume finalEpoch(
        grid, finalEpochInterfaces, 300, 2);
    MacVelocityField finalVelocity(grid);
    CellScalarField finalPressure(grid);
    const auto finalFluid = projectVelocityWithMovingInterfaces(
        grid, finalVelocity, finalPressure, finalEpochInterfaces,
        projectionSettings());
    const auto fullDomainTerminal = finalEpoch.evaluate(
        grid, finalVelocity, finalFluid, {0.4, 0.5, 0.4, true});
    const FaceAlignedMovingInterface collidingInterfaces(
        grid, pistonFaces(grid, 0.25, 2));
    expectRejected(
        [&] { static_cast<void>(rebasePlanarMovingControlVolume(
            grid, finalEpoch, collidingInterfaces,
            fullDomainTerminal)); },
        "rebase: moving surface cannot consume its resolved opening");
}

void testValidationAndFailedLedgers() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(
        grid, pistonFaces(grid));
    const PlanarMovingControlVolume controlVolume(
        grid, interfaces, 300, 2);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto fluid = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());

    expectRejected(
        [&] { static_cast<void>(PlanarMovingControlVolume(
            grid, interfaces, 300, 6)); },
        "validation: opening cannot alias the moving surface plane");
    expectRejected(
        [&] { static_cast<void>(PlanarMovingControlVolume(
            grid, interfaces, 999, 2)); },
        "validation: absent moving surface stable ID is rejected");

    const PlanarControlVolumeStep mismatchedGeometry{0.0, 0.05, 0.4};
    const auto mismatched = controlVolume.evaluate(
        grid, velocity, fluid, mismatchedGeometry);
    check(!mismatched.accepted
              && mismatched.maximumSurfaceVelocityErrorMetersPerSecond > 0.0
              && std::abs(mismatched.surfaceGeometryResidualCubicMeters) > 0.0,
          "validation: geometry that disagrees with surface speed fails acceptance");

    auto missingOpeningFlux = velocity;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            missingOpeningFlux.xFaces()[grid.cellIndex(2, j, k)] = 0.0;
        }
    }
    const auto missingFlux = controlVolume.evaluate(
        grid, missingOpeningFlux, fluid, {0.0, 0.1, 0.4});
    check(!missingFlux.accepted
              && std::abs(missingFlux.continuityResidualCubicMeters) > 0.5,
          "validation: absent resolved opening transport fails GCL acceptance");

    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, fluid, {0.0, 0.5, 0.4})); },
        "validation: motion reaching the next grid face requires topology rebase");
    auto failedFluid = fluid;
    failedFluid.projection.converged = false;
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, failedFluid, {0.0, 0.1, 0.4})); },
        "validation: unaccepted fluid projection is rejected");

    auto inconsistentPower = fluid;
    inconsistentPower.surfaces.front()
        .constraintReactionPowerWatts += 1.0;
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, inconsistentPower, {0.0, 0.1, 0.4})); },
        "validation: face and surface reaction-power ledgers must agree");

    auto invalidSettings = simwing::fsi::fluid::PlanarControlVolumeSettings{};
    invalidSettings.relativeVolumeTolerance =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, fluid, {0.0, 0.1, 0.4},
            invalidSettings)); },
        "validation: non-finite conservation tolerances are rejected");
}

} // namespace

int main() {
    testAcceleratedOpenPistonAndGcl();
    testPorousOpeningAcrossControlAndCutSurface();
    testAllAxisControlVolumes();
    testExplicitTopologyRebase();
    testValidationAndFailedLedgers();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing moving control-volume check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing moving control-volume checks passed");
    return 0;
}
