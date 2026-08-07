#include "fluid/moving_control_volume.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PlanarControlVolumeStep;
using simwing::fsi::fluid::PlanarMovingControlVolume;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;

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

std::vector<GridFaceMovingInterface> pistonFaces(
    const PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond = 0.25) {
    std::vector<GridFaceMovingInterface> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                300, 9, 9, GridFaceAxis::X, 6, j, k,
                speedMetersPerSecond,
            });
        }
    }
    return result;
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

void testAllAxisControlVolumes() {
    const PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    const auto counts = grid.cellCounts();
    for (std::size_t axisIndex = 0; axisIndex < 3; ++axisIndex) {
        const auto axis = static_cast<GridFaceAxis>(axisIndex);
        std::vector<GridFaceMovingInterface> faces;
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    const std::size_t coordinate = axis == GridFaceAxis::X
                        ? i : (axis == GridFaceAxis::Y ? j : k);
                    if (coordinate == 3) {
                        faces.push_back({
                            400 + axisIndex,
                            20 + axisIndex,
                            20 + axisIndex,
                            axis,
                            i,
                            j,
                            k,
                            0.125,
                        });
                    }
                }
            }
        }
        const FaceAlignedMovingInterface interfaces(
            grid, std::move(faces));
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
    }
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
    inconsistentPower.surfaces.front().pressurePowerWatts += 1.0;
    expectRejected(
        [&] { static_cast<void>(controlVolume.evaluate(
            grid, velocity, inconsistentPower, {0.0, 0.1, 0.4})); },
        "validation: face and surface pressure-power ledgers must agree");

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
    testAllAxisControlVolumes();
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
