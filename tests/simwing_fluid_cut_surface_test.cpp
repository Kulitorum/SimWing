#include "fluid/planar_cut_surface.h"

#include <algorithm>
#include <array>
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
using simwing::fsi::fluid::PlanarCutSurfacePressureSettings;
using simwing::fsi::fluid::PlanarMovingControlVolume;
using simwing::fsi::fluid::Vector3;
using simwing::fsi::fluid::evaluatePlanarCutSurfacePressure;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;
using simwing::fsi::fluid::rebasePlanarMovingControlVolume;
using simwing::fsi::fluid::resamplePlanarCutSurfaceReaction;

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

PeriodicCartesianGrid makeGrid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

std::vector<GridFaceMovingInterface> completePlane(
    const PeriodicCartesianGrid& grid,
    const GridFaceAxis axis,
    const std::size_t planeCoordinate,
    const std::uint64_t surfaceStableId,
    const std::uint64_t regionStableId,
    const double speedMetersPerSecond = 0.125) {
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

MovingInterfaceProjectionSettings projectionSettings() {
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.1;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return settings;
}

double normalCoordinate(const Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return value.x;
    case GridFaceAxis::Y:
        return value.y;
    case GridFaceAxis::Z:
        return value.z;
    }
    throw std::invalid_argument("unknown test axis");
}

void testAllAxesAndPeriodicImages() {
    const auto grid = makeGrid();
    for (std::size_t axisIndex = 0; axisIndex < 3; ++axisIndex) {
        const auto axis = static_cast<GridFaceAxis>(axisIndex);
        const std::uint64_t surfaceStableId = 400 + axisIndex;
        const std::uint64_t regionStableId = 20 + axisIndex;
        const FaceAlignedMovingInterface interfaces(
            grid, completePlane(
                grid, axis, 3, surfaceStableId, regionStableId));
        const PlanarMovingControlVolume controlVolume(
            grid, interfaces, surfaceStableId, 1);
        MacVelocityField velocity(grid);
        CellScalarField pressure(grid);
        const auto fluid = projectVelocityWithMovingInterfaces(
            grid, velocity, pressure, interfaces, projectionSettings());
        PlanarCutSurfacePressureSettings settings;
        settings.momentReferenceMeters = {0.25, -0.5, 0.75};
        const auto first = evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 0.25, 3.25, settings);
        const auto replay = evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 0.25, 3.25, settings);
        check(first == replay,
              "cut surface: every axis replays bit-for-bit");
        check(first.version
                  == simwing::fsi::fluid::planarCutSurfacePressureVersion
                  && first.sourceInterfaceVersion
                      == simwing::fsi::fluid::faceAlignedMovingInterfaceVersion
                  && first.surfaceStableId == surfaceStableId
                  && first.fluidRegionStableId == regionStableId
                  && first.axis == axis
                  && first.movingPlaneCoordinate == 3
                  && first.faceCount == 16
                  && first.faces.size() == 16,
              "cut surface: stable topology and source version are explicit");
        checkNear(first.gridPlaneCoordinateMeters, 3.0, 0.0,
                  "cut surface: Eulerian plane is retained");
        checkNear(first.physicalPlaneCoordinateMeters, 3.25, 0.0,
                  "cut surface: physical plane is explicit");
        checkNear(first.periodicPositionResidualMeters, 0.0, 0.0,
                  "cut surface: physical offset matches the current epoch");
        checkNear(first.areaSquareMeters, 16.0, 0.0,
                  "cut surface: complete plane has analytic area");
        check(first.forceResidualNormNewtons < 1.0e-12
                  && std::abs(first.powerResidualWatts) < 1.0e-12
                  && first.finite && first.accepted,
              "cut surface: source and physical reaction ledgers close");
        checkNear(
            normalCoordinate(first.pressureForceNewtons, axis)
                * projectionSettings().projection.timeStepSeconds,
            -9.6, 5.0e-12,
            "cut surface: complete reaction balances plug momentum on every axis");
        const auto startSample = resamplePlanarCutSurfaceReaction(
            grid, controlVolume, first, 0.0, 3.0, 0.0, settings);
        check(startSample.accepted && startSample.kinematicsResampled
                  && startSample.pressureForceNewtons
                      == first.pressureForceNewtons
                  && startSample.normalVelocityMetersPerSecond == 0.0
                  && startSample.pressurePowerWatts == 0.0
                  && startSample.reactionSourcePhysicalPlaneCoordinateMeters
                      == first.physicalPlaneCoordinateMeters
                  && startSample.reactionSourceNormalVelocityMetersPerSecond
                      == first.normalVelocityMetersPerSecond,
              "cut surface: average reaction resamples start kinematics without changing force");
        checkNear(
            fluid.projection.kineticEnergyAfterJoules
                + 0.5 * (startSample.pressurePowerWatts
                         + first.pressurePowerWatts)
                    * projectionSettings().projection.timeStepSeconds,
            0.0, 2.0e-14,
            "cut surface: average reaction work balances plug-flow kinetic energy");
        for (const auto& face : first.faces) {
            checkNear(normalCoordinate(
                          face.gridLowerCornerMeters, axis),
                      3.0, 0.0,
                      "cut surface: grid tiles remain on the MAC plane");
            checkNear(normalCoordinate(
                          face.physicalLowerCornerMeters, axis),
                      3.25, 0.0,
                      "cut surface: pressure geometry lies on the cut plane");
        }

        const auto terminalControl = controlVolume.evaluate(
            grid, velocity, fluid, {0.0, 1.0, 8.0, true});
        const auto terminal = evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 1.0, 4.0, settings);
        const FaceAlignedMovingInterface wrappedInterfaces(
            grid, completePlane(
                grid, axis, 0, surfaceStableId, regionStableId));
        const auto rebase = rebasePlanarMovingControlVolume(
            grid, controlVolume, wrappedInterfaces, terminalControl);
        auto wrappedVelocity = velocity;
        auto wrappedPressure = pressure;
        const auto wrappedFluid = projectVelocityWithMovingInterfaces(
            grid, wrappedVelocity, wrappedPressure,
            wrappedInterfaces, projectionSettings());
        const auto wrapped = evaluatePlanarCutSurfacePressure(
            grid, rebase.controlVolume, wrappedFluid,
            0.0, 4.0, settings);
        check(terminalControl.accepted
                  && rebase.diagnostics.accepted
                  && terminal.accepted && wrapped.accepted,
              "cut surface: terminal and rebased pressure geometry are accepted");
        checkNear(terminal.gridPlaneCoordinateMeters, 3.0, 0.0,
                  "cut surface: terminal pressure retains the old grid plane");
        checkNear(wrapped.gridPlaneCoordinateMeters, 0.0, 0.0,
                  "cut surface: periodic rebase uses the wrapped grid plane");
        checkNear(wrapped.physicalPlaneCoordinateMeters, 4.0, 0.0,
                  "cut surface: unwrapped physical position survives rebase");
        checkNear(wrapped.periodicPositionResidualMeters, 0.0, 0.0,
                  "cut surface: periodic image is geometrically congruent");
    }
}

void testStrictValidation() {
    const auto grid = makeGrid();
    const FaceAlignedMovingInterface interfaces(
        grid, completePlane(grid, GridFaceAxis::X, 3, 500, 30));
    const PlanarMovingControlVolume controlVolume(
        grid, interfaces, 500, 1);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto fluid = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, projectionSettings());

    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, -0.01, 2.99)); },
        "cut surface validation: negative partial-cell offsets are rejected");
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 1.01, 4.01)); },
        "cut surface validation: offsets beyond the partial cell are rejected");
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 0.25,
            std::numeric_limits<double>::quiet_NaN())); },
        "cut surface validation: the physical plane must be finite");

    const auto mismatchedPosition = evaluatePlanarCutSurfacePressure(
        grid, controlVolume, fluid, 0.25, 3.3);
    check(mismatchedPosition.finite && !mismatchedPosition.accepted
              && mismatchedPosition.periodicPositionResidualMeters > 0.04,
          "cut surface validation: a noncongruent physical plane is not accepted");

    auto failedFluid = fluid;
    failedFluid.projection.converged = false;
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, failedFluid, 0.25, 3.25)); },
        "cut surface validation: failed fluid projection is rejected");

    auto duplicateTile = fluid;
    duplicateTile.faces[1] = duplicateTile.faces[0];
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, duplicateTile, 0.25, 3.25)); },
        "cut surface validation: duplicate pressure tiles are rejected");

    auto corruptedFaceForce = fluid;
    corruptedFaceForce.faces.front()
        .constraintReactionForceNewtons.x += 1.0;
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, corruptedFaceForce, 0.25, 3.25)); },
        "cut surface validation: local traction and force must agree");

    auto corruptedFacePower = fluid;
    corruptedFacePower.faces.front()
        .constraintReactionPowerWatts += 1.0;
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, corruptedFacePower, 0.25, 3.25)); },
        "cut surface validation: local force and power must agree");

    auto nonrigidVelocity = fluid;
    auto& changedFace = nonrigidVelocity.faces.front();
    const double previousPower = changedFace.constraintReactionPowerWatts;
    changedFace.normalVelocityMetersPerSecond += 0.01;
    changedFace.constraintReactionPowerWatts =
        changedFace.constraintReactionForceNewtons.x
        * changedFace.normalVelocityMetersPerSecond;
    nonrigidVelocity.surfaces.front().constraintReactionPowerWatts +=
        changedFace.constraintReactionPowerWatts - previousPower;
    const auto nonrigid = evaluatePlanarCutSurfacePressure(
        grid, controlVolume, nonrigidVelocity, 0.25, 3.25);
    check(nonrigid.finite && !nonrigid.accepted
              && nonrigid.maximumNormalVelocitySpreadMetersPerSecond > 0.009,
          "cut surface validation: nonrigid face velocity is not accepted");

    auto corruptedAggregate = fluid;
    corruptedAggregate.surfaces.front()
        .constraintReactionForceNewtons.x += 1.0;
    const auto corrupted = evaluatePlanarCutSurfacePressure(
        grid, controlVolume, corruptedAggregate, 0.25, 3.25);
    check(corrupted.finite && !corrupted.accepted
              && corrupted.forceResidualNormNewtons > 0.9,
          "cut surface validation: corrupted force aggregates fail acceptance");

    auto unacceptedReaction = evaluatePlanarCutSurfacePressure(
        grid, controlVolume, fluid, 0.25, 3.25);
    unacceptedReaction.accepted = false;
    expectRejected(
        [&] { static_cast<void>(resamplePlanarCutSurfaceReaction(
            grid, controlVolume, unacceptedReaction,
            0.0, 3.0, 0.0)); },
        "cut surface validation: unaccepted reactions cannot be resampled");

    PlanarCutSurfacePressureSettings invalidSettings;
    invalidSettings.relativePowerTolerance = -1.0;
    expectRejected(
        [&] { static_cast<void>(evaluatePlanarCutSurfacePressure(
            grid, controlVolume, fluid, 0.25, 3.25,
            invalidSettings)); },
        "cut surface validation: invalid numerical settings are rejected");
}

} // namespace

int main() {
    testAllAxesAndPeriodicImages();
    testStrictValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d planar cut-surface check(s) failed\n",
                     failures);
        return 1;
    }
    return 0;
}
