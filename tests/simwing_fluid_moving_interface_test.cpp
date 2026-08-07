#include "fluid/moving_interface.h"

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
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::Vector3;
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

void checkVectorNear(const Vector3& actual,
                     const Vector3& expected,
                     const double tolerance,
                     const char* message) {
    if (!std::isfinite(actual.x) || !std::isfinite(actual.y)
        || !std::isfinite(actual.z)
        || std::abs(actual.x - expected.x) > tolerance
        || std::abs(actual.y - expected.y) > tolerance
        || std::abs(actual.z - expected.z) > tolerance) {
        std::fprintf(
            stderr,
            "FAIL: %s (actual [%.17g %.17g %.17g], expected [%.17g %.17g %.17g])\n",
            message, actual.x, actual.y, actual.z,
            expected.x, expected.y, expected.z);
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

std::vector<GridFaceMovingInterface> slabFaces(
    const PeriodicCartesianGrid& grid,
    const double leftVelocityMetersPerSecond = 0.25,
    const double rightVelocityMetersPerSecond = 0.25) {
    std::vector<GridFaceMovingInterface> result;
    const auto counts = grid.cellCounts();
    result.reserve(2 * counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                100, 1, 2, GridFaceAxis::X, 2, j, k,
                leftVelocityMetersPerSecond});
            result.push_back({
                200, 2, 1, GridFaceAxis::X, 6, j, k,
                rightVelocityMetersPerSecond});
        }
    }
    return result;
}

MovingInterfaceProjectionSettings pistonSettings() {
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.4;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond = 1.0e-12;
    return settings;
}

void fillSlabPressure(const PeriodicCartesianGrid& grid,
                      CellScalarField& pressure,
                      const double insidePascals = 110.0) {
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                pressure.values()[grid.cellIndex(i, j, k)] =
                    i >= 2 && i < 6 ? insidePascals : 0.0;
            }
        }
    }
}

void fillDeterministicVelocity(MacVelocityField& velocity) {
    for (std::size_t index = 0; index < velocity.xFaces().size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        velocity.xFaces()[index] = 0.2 * std::sin(0.31 * sample);
        velocity.yFaces()[index] = 0.15 * std::cos(0.23 * sample);
        velocity.zFaces()[index] = 0.1 * std::sin(0.17 * sample + 0.2);
    }
}

void testTopologyCanonicalizationAndValidation() {
    const auto grid = pistonGrid();
    auto authored = slabFaces(grid);
    auto reversed = authored;
    std::reverse(reversed.begin(), reversed.end());
    const FaceAlignedMovingInterface first(grid, authored);
    const FaceAlignedMovingInterface second(grid, reversed);
    check(first == second,
          "topology: authored face order canonicalizes bit-for-bit");
    check(first.version()
              == simwing::fsi::fluid::faceAlignedMovingInterfaceVersion
              && first.faceCount() == 12 && first.regionCount() == 2,
          "topology: version, complete piston faces, and two regions are explicit");
    check(first.regionStableIds()[0] == 1
              && first.regionStableIds()[1] == 2,
          "topology: stable fluid regions are canonicalized by ID");
    check(first.cellRegionStableIds()[grid.cellIndex(0, 0, 0)] == 1
              && first.cellRegionStableIds()[grid.cellIndex(3, 0, 0)] == 2,
          "topology: disconnected cells retain their stable region identity");

    auto invalid = authored;
    invalid[0].surfaceStableId = 0;
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: zero surface stable ID is rejected");
    invalid = authored;
    invalid[0].plusRegionStableId = invalid[0].minusRegionStableId;
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: one equal-side face cannot contradict separating topology");
    invalid = authored;
    invalid[0].normalVelocityMetersPerSecond =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: non-finite interface velocity is rejected");
    invalid = authored;
    invalid[0].i = grid.cellCounts().x;
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: out-of-range face coordinate is rejected");
    invalid = authored;
    invalid[0].axis = static_cast<GridFaceAxis>(99);
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: unknown face axis is rejected");
    invalid = authored;
    invalid.push_back(invalid.front());
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: multiple interfaces on one grid face are rejected");
    invalid = authored;
    invalid.back().surfaceStableId = invalid.front().surfaceStableId;
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, invalid)); },
        "validation: one surface ID cannot alias different oriented regions");

    std::vector<GridFaceMovingInterface> onePlane;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            onePlane.push_back(
                {100, 1, 2, GridFaceAxis::X, 2, j, k, 0.0});
        }
    }
    expectRejected(
        [&] { static_cast<void>(FaceAlignedMovingInterface(grid, onePlane)); },
        "validation: a periodic plane that does not partition regions is rejected");

    for (auto& face : onePlane) {
        face.minusRegionStableId = 9;
        face.plusRegionStableId = 9;
    }
    const FaceAlignedMovingInterface nonseparating(
        grid, std::move(onePlane));
    check(nonseparating.regionCount() == 1
              && nonseparating.regionStableIds().front() == 9
              && std::ranges::all_of(
                  nonseparating.cellRegionStableIds(),
                  [](const std::uint64_t stableId) {
                      return stableId == 9;
                  }),
          "topology: equal side IDs retain one region around a resolved opening");
}

void testAnalyticTranslatingSlabPressureWork() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(grid, slabFaces(grid));
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 0.25);
    CellScalarField pressure(grid);
    fillSlabPressure(grid, pressure);
    const auto originalVelocity = velocity;
    const auto originalPressure = pressure;
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, pistonSettings());

    check(diagnostics.projection.converged
              && diagnostics.projection.iterationCount == 0,
          "piston: already compatible translation takes the exact zero-correction path");
    check(velocity == originalVelocity && pressure == originalPressure,
          "piston: compatible velocity and two-region pressure remain bit-identical");
    check(diagnostics.interfaceFaceCount == 12
              && diagnostics.fluidRegionCount == 2
              && diagnostics.faces.size() == 12
              && diagnostics.surfaces.size() == 2,
          "piston: every face, fluid region, and membrane surface is diagnosed");
    checkNear(diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond,
              0.0, 0.0,
              "piston: equal wall velocities preserve each sealed region volume");
    checkNear(diagnostics.maximumNormalVelocityErrorMetersPerSecond,
              0.0, 0.0,
              "piston: no-penetration velocity is exact on every constrained face");
    check(diagnostics.regions[0].stableId == 1
              && diagnostics.regions[1].stableId == 2,
          "piston: region diagnostics retain canonical stable IDs");
    checkNear(diagnostics.regions[0].pressureMeanAfterPascals,
              0.0, 0.0,
              "piston: outside pressure gauge is retained");
    checkNear(diagnostics.regions[1].pressureMeanAfterPascals,
              110.0, 0.0,
              "piston: inside pressure gauge is retained");

    const auto& left = diagnostics.surfaces[0];
    const auto& right = diagnostics.surfaces[1];
    const auto& firstRightFace = diagnostics.faces[1];
    check(firstRightFace.surfaceStableId == 200
              && firstRightFace.minusRegionStableId == 2
              && firstRightFace.plusRegionStableId == 1
              && firstRightFace.axis == GridFaceAxis::X
              && firstRightFace.i == 6
              && firstRightFace.j == 0
              && firstRightFace.k == 0,
          "piston: face diagnostics preserve canonical grid and region identity");
    checkVectorNear(firstRightFace.lowerCornerMeters, {3.0, 0.0, 0.0}, 0.0,
                    "piston: face tile lower corner is exact grid geometry");
    checkVectorNear(firstRightFace.upperCornerMeters, {3.0, 1.0, 1.0}, 0.0,
                    "piston: face tile upper corner is exact grid geometry");
    checkNear(firstRightFace.areaSquareMeters, 1.0, 0.0,
              "piston: one face tile carries its exact area");
    checkVectorNear(firstRightFace.pressureTractionPascals,
                    {110.0, 0.0, 0.0}, 0.0,
                    "piston: one face tile carries its pressure traction");
    checkVectorNear(firstRightFace.pressureForceNewtons,
                    {110.0, 0.0, 0.0}, 0.0,
                    "piston: one face tile carries its integrated force");
    checkNear(firstRightFace.pressurePowerWatts, 27.5, 0.0,
              "piston: one face tile carries its pressure power");
    check(left.stableId == 100 && right.stableId == 200,
          "piston: surface diagnostics are canonicalized by stable ID");
    checkNear(left.areaSquareMeters, 6.0, 0.0,
              "piston: left grid-face area integrates analytically");
    checkNear(right.areaSquareMeters, 6.0, 0.0,
              "piston: right grid-face area integrates analytically");
    checkVectorNear(left.pressureForceNewtons, {-660.0, 0.0, 0.0}, 0.0,
                    "piston: left pressure force is analytic and outward");
    checkVectorNear(right.pressureForceNewtons, {660.0, 0.0, 0.0}, 0.0,
                    "piston: right pressure force is analytic and outward");
    checkNear(left.maximumPressureTractionDeviationPascals, 0.0, 0.0,
              "piston: left pressure traction is exactly uniform");
    checkNear(right.maximumPressureTractionDeviationPascals, 0.0, 0.0,
              "piston: right pressure traction is exactly uniform");
    checkVectorNear(left.pressureImpulseNewtonSeconds,
                    {-264.0, 0.0, 0.0}, 2.0e-13,
                    "piston: left pressure impulse matches the temporal canonical");
    checkVectorNear(right.pressureImpulseNewtonSeconds,
                    {264.0, 0.0, 0.0}, 2.0e-13,
                    "piston: right pressure impulse matches the Structure exchange");
    checkNear(left.pressureWorkJoules, -66.0, 5.0e-14,
              "piston: left wall pressure work is negative during translation");
    checkNear(right.pressureWorkJoules, 66.0, 5.0e-14,
              "piston: right wall pressure work equals pressure times swept volume");
    checkVectorNear(diagnostics.totalPressureForceNewtons, {}, 0.0,
                    "piston: closed slab has zero net pressure force");
    checkVectorNear(diagnostics.totalPressureImpulseNewtonSeconds, {}, 0.0,
                    "piston: closed slab has zero net pressure impulse");
    checkNear(diagnostics.totalPressureWorkJoules, 0.0, 0.0,
              "piston: rigid translation cancels the two wall work ledgers");
    check(diagnostics.finite,
          "piston: all moving-interface ledgers remain finite");
}

void testAllAxisConstraintStencils() {
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
                    if (coordinate == 1) {
                        faces.push_back({
                            100 + axisIndex, 1, 2, axis, i, j, k, 0.125});
                    } else if (coordinate == 3) {
                        faces.push_back({
                            200 + axisIndex, 2, 1, axis, i, j, k, 0.125});
                    }
                }
            }
        }
        const FaceAlignedMovingInterface interfaces(grid, std::move(faces));
        MacVelocityField velocity(grid);
        switch (axis) {
        case GridFaceAxis::X:
            std::ranges::fill(velocity.xFaces(), 0.125);
            break;
        case GridFaceAxis::Y:
            std::ranges::fill(velocity.yFaces(), 0.125);
            break;
        case GridFaceAxis::Z:
            std::ranges::fill(velocity.zFaces(), 0.125);
            break;
        }
        CellScalarField pressure(grid);
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    const std::size_t coordinate = axis == GridFaceAxis::X
                        ? i : (axis == GridFaceAxis::Y ? j : k);
                    pressure.values()[grid.cellIndex(i, j, k)] =
                        coordinate >= 1 && coordinate < 3 ? 7.0 : 0.0;
                }
            }
        }
        const auto originalVelocity = velocity;
        const auto originalPressure = pressure;
        auto settings = pistonSettings();
        settings.projection.timeStepSeconds = 0.1;
        const auto diagnostics = projectVelocityWithMovingInterfaces(
            grid, velocity, pressure, interfaces, settings);
        check(diagnostics.projection.converged
                  && diagnostics.projection.iterationCount == 0
                  && diagnostics.maximumNormalVelocityErrorMetersPerSecond
                      == 0.0,
              "axes: every MAC component preserves its exact interface velocity");
        check(velocity == originalVelocity && pressure == originalPressure,
              "axes: compatible pressure and translation need no correction");
        Vector3 expectedRightForce;
        if (axis == GridFaceAxis::X) {
            expectedRightForce.x = 112.0;
        } else if (axis == GridFaceAxis::Y) {
            expectedRightForce.y = 112.0;
        } else {
            expectedRightForce.z = 112.0;
        }
        checkVectorNear(
            diagnostics.surfaces[1].pressureForceNewtons,
            expectedRightForce, 0.0,
            "axes: oriented pressure force uses the selected face component");
    }
}

void testDisturbedProjectionAndDeterminism() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(grid, slabFaces(grid));
    MacVelocityField firstVelocity(grid);
    fillDeterministicVelocity(firstVelocity);
    MacVelocityField secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    fillSlabPressure(grid, firstPressure);
    CellScalarField secondPressure = firstPressure;
    const auto settings = pistonSettings();
    const auto firstDiagnostics = projectVelocityWithMovingInterfaces(
        grid, firstVelocity, firstPressure, interfaces, settings);
    const auto secondDiagnostics = projectVelocityWithMovingInterfaces(
        grid, secondVelocity, secondPressure, interfaces, settings);

    check(firstDiagnostics.projection.converged
              && firstDiagnostics.projection.iterationCount > 0,
          "projection: disturbed regional flow converges through disconnected CG");
    check(firstVelocity == secondVelocity
              && firstPressure == secondPressure
              && firstDiagnostics == secondDiagnostics,
          "projection: identical moving-interface solves replay bit-for-bit");
    check(firstDiagnostics.projection.divergenceL2AfterPerSecond < 2.0e-11,
          "projection: regional divergence is reduced to the solver budget");
    check(firstDiagnostics.maximumNormalVelocityErrorMetersPerSecond == 0.0,
          "projection: constrained normal velocities remain bit-exact");
    check(firstDiagnostics.surfaces[0]
                  .maximumPressureTractionDeviationPascals > 0.0
              || firstDiagnostics.surfaces[1]
                  .maximumPressureTractionDeviationPascals > 0.0,
          "projection: disturbed pressure reports its nonuniform surface traction");
    checkNear(
        firstDiagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond,
        0.0, 2.0e-15,
        "projection: each disturbed region remains volume compatible");
    checkNear(firstDiagnostics.regions[0].pressureMeanAfterPascals,
              firstDiagnostics.regions[0].pressureMeanBeforePascals, 2.0e-14,
              "projection: outside pressure gauge survives the correction");
    checkNear(firstDiagnostics.regions[1].pressureMeanAfterPascals,
              firstDiagnostics.regions[1].pressureMeanBeforePascals, 2.0e-14,
              "projection: inside pressure gauge survives the correction");
    check(firstDiagnostics.finite,
          "projection: disturbed moving-interface diagnostics are finite");
}

void testIncompatibilityValidationAndRollback() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface incompatible(
        grid, slabFaces(grid, 0.25, 0.0));
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid, 3.0);
    const auto originalVelocity = velocity;
    const auto originalPressure = pressure;
    const auto diagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, incompatible, pistonSettings());
    check(!diagnostics.projection.converged,
          "compatibility: one-sided fixed-grid volume change is rejected");
    checkNear(diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond,
              1.5, 2.0e-16,
              "compatibility: rejected piston reports its exact missing volume rate");
    check(velocity == originalVelocity && pressure == originalPressure,
          "compatibility: rejected volume change commits neither field");

    const FaceAlignedMovingInterface compatible(grid, slabFaces(grid));
    fillDeterministicVelocity(velocity);
    pressure = CellScalarField(grid, 2.0);
    const auto beforeTruncatedVelocity = velocity;
    const auto beforeTruncatedPressure = pressure;
    auto truncated = pistonSettings();
    truncated.projection.absoluteResidualTolerance = 1.0e-30;
    truncated.projection.relativeResidualTolerance = 0.0;
    truncated.projection.maximumIterations = 1;
    const auto truncatedDiagnostics = projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, compatible, truncated);
    check(!truncatedDiagnostics.projection.converged,
          "rollback: intentionally truncated regional CG does not converge");
    check(velocity == beforeTruncatedVelocity
              && pressure == beforeTruncatedPressure,
          "rollback: failed regional CG commits neither pressure nor velocity");

    auto invalidSettings = pistonSettings();
    invalidSettings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(projectVelocityWithMovingInterfaces(
            grid, velocity, pressure, compatible, invalidSettings)); },
        "validation: non-finite compatibility tolerance is rejected");

    const PeriodicCartesianGrid otherGrid(
        {8, 2, 3}, {}, {8.0, 2.0, 3.0});
    const FaceAlignedMovingInterface foreign(
        otherGrid, slabFaces(otherGrid));
    const auto beforeForeignVelocity = velocity;
    expectRejected(
        [&] { static_cast<void>(projectVelocityWithMovingInterfaces(
            grid, velocity, pressure, foreign, pistonSettings())); },
        "validation: a moving interface from another physical grid is rejected");
    check(velocity == beforeForeignVelocity,
          "validation: foreign topology is rejected before field mutation");
}

} // namespace

int main() {
    testTopologyCanonicalizationAndValidation();
    testAnalyticTranslatingSlabPressureWork();
    testAllAxisConstraintStencils();
    testDisturbedProjectionAndDeterminism();
    testIncompatibilityValidationAndRollback();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing moving-interface check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing moving-interface checks passed");
    return 0;
}
