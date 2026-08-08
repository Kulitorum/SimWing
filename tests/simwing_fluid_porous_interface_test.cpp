#include "fluid/porous_interface.h"
#include "fluid/porous_flow.h"
#include "fluid/projection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::fsi::fluid;

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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void testConstitutiveLawAndInverse() {
    const DarcyForchheimerResistance combined{100.0, 25.0};
    checkNear(porousPressureJumpPascals(combined, 2.0),
              -300.0, 0.0,
              "porous law opposes positive normal flow");
    checkNear(porousPressureJumpPascals(combined, -2.0),
              300.0, 0.0,
              "porous law opposes negative normal flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  combined, -300.0),
              2.0, 5.0e-16,
              "porous inverse recovers positive combined-law flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  combined, 300.0),
              -2.0, 5.0e-16,
              "porous inverse recovers negative combined-law flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {50.0, 0.0}, -100.0),
              2.0, 0.0,
              "porous inverse supports a pure Darcy fit");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {0.0, 25.0}, -100.0),
              2.0, 0.0,
              "porous inverse supports a pure quadratic fit");
    check(porousPressureJumpPascals(combined, 0.0) == 0.0
              && porousRelativeNormalVelocityMetersPerSecond(
                     combined, 0.0) == 0.0,
          "porous law preserves the exact zero state");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {1.0e308, 1.0e308}, -1.0e308),
              0.5 * (std::sqrt(5.0) - 1.0), 2.0e-16,
              "porous inverse remains stable at large matched scales");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {0.0, 1.0e-308}, -1.0e308),
              1.0e308, 2.0e292,
              "porous quadratic inverse avoids an intermediate overflow");

    for (const double velocity : {-20.0, -1.0, -1.0e-6,
                                  1.0e-6, 1.0, 20.0}) {
        const double jump = porousPressureJumpPascals(combined, velocity);
        checkNear(porousRelativeNormalVelocityMetersPerSecond(
                      combined, jump),
                  velocity, 3.0e-15 * std::max(1.0, std::abs(velocity)),
                  "porous combined law round-trips across flow scales");
        check(-jump * velocity >= 0.0,
              "porous pressure work is dissipative");
    }
}

std::vector<PorousGridFaceCrossing> porousPlane(
    const PeriodicCartesianGrid& grid) {
    std::vector<PorousGridFaceCrossing> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                10, 1, 2, GridFaceAxis::X,
                1, j, k, 0.35, 0.5, {100.0, 20.0}});
        }
    }
    return result;
}

void testFluxDrivenProjectionCanonical() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 2.0);
    auto definitions = porousPlane(grid);
    auto reversedDefinitions = definitions;
    std::reverse(reversedDefinitions.begin(), reversedDefinitions.end());
    const PorousPressureJumpField porous(
        grid, velocity, definitions);
    const PorousPressureJumpField reversed(
        grid, velocity, reversedDefinitions);

    check(porous.pressureJumps() == reversed.pressureJumps()
              && std::ranges::equal(porous.samples(), reversed.samples()),
          "porous field canonicalizes reversed authored crossings");
    check(porous.samples().size() == 6
              && porous.pressureJumps().faceCount() == 6,
          "porous field retains every resolved plane tile");
    for (const auto& sample : porous.samples()) {
        checkNear(sample.fluidNormalVelocityMetersPerSecond,
                  2.0, 0.0,
                  "porous sample reads the owning MAC normal velocity");
        checkNear(sample.relativeNormalVelocityMetersPerSecond,
                  1.5, 0.0,
                  "porous sample subtracts sheet normal velocity");
        checkNear(sample.pressureJump.pressureJumpPascals,
                  -195.0, 0.0,
                  "porous sample applies the calibrated signed jump");
        checkNear(sample.volumeFlowRateCubicMetersPerSecond,
                  1.5, 0.0,
                  "porous sample retains resolved tile volume flux");
        checkNear(sample.dissipationWatts,
                  292.5, 0.0,
                  "porous sample closes positive pressure dissipation");
    }
    checkNear(porous.totalDissipationWatts(),
              1755.0, 0.0,
              "porous field aggregates tile dissipation");

    std::vector<GridFacePressureJump> drivenJumps(
        porous.pressureJumps().faces().begin(),
        porous.pressureJumps().faces().end());
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            drivenJumps.push_back({
                20, 2, 1, GridFaceAxis::X,
                3, j, k, 195.0, 0.5});
        }
    }
    const SharpPressureJumpField completeJumps(
        grid, std::move(drivenJumps));
    CellScalarField pressure(grid);
    ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.absoluteResidualTolerance = 1.0e-11;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    const auto diagnostics = projectVelocityWithPressureJumps(
        grid, velocity, pressure, completeJumps, settings);
    check(diagnostics.converged
              && diagnostics.pressureJumpFaceCount == 12
              && diagnostics.divergenceL2AfterPerSecond < 1.0e-13,
          "porous plane plus pressure source projects a compatible flow");
    check(std::ranges::all_of(
              velocity.xFaces(),
              [](const double sample) { return sample == 2.0; })
              && std::ranges::all_of(
                  velocity.yFaces(),
                  [](const double sample) { return sample == 0.0; })
              && std::ranges::all_of(
                  velocity.zFaces(),
                  [](const double sample) { return sample == 0.0; }),
          "flux-driven porous canonical creates no spurious velocity");
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i == 0 || i == 3
                    ? 97.5 : -97.5;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(pressure.values()[grid.cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 2.0e-12,
          "flux-driven porous canonical preserves its analytic pressure loss");
}

PorousProjectionSettings coupledSettings() {
    PorousProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.0;
    settings.projection.timeStepSeconds = 0.1;
    settings.projection.absoluteResidualTolerance = 1.0e-12;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    settings.absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-12;
    settings.relativeNormalVelocityTolerance = 1.0e-12;
    settings.absolutePressureJumpTolerancePascals = 1.0e-11;
    settings.relativePressureJumpTolerance = 1.0e-12;
    settings.relaxation = 0.5;
    settings.maximumNonlinearIterations = 100;
    return settings;
}

std::vector<PorousGridFaceCrossing> coupledPorousPlane(
    const PeriodicCartesianGrid& grid,
    const double surfaceVelocityMetersPerSecond = 0.0) {
    std::vector<PorousGridFaceCrossing> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                30, 1, 2, GridFaceAxis::X, 1, j, k, 0.35,
                surfaceVelocityMetersPerSecond, {10.0, 0.0}});
        }
    }
    return result;
}

SharpPressureJumpField coupledDrivingPlane(
    const PeriodicCartesianGrid& grid,
    const double pressureRisePascals = 20.0) {
    std::vector<GridFacePressureJump> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                40, 2, 1, GridFaceAxis::X, 3, j, k,
                pressureRisePascals, 0.5});
        }
    }
    return SharpPressureJumpField(grid, std::move(faces));
}

void testImplicitGridPressureFluxCoupling() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    const auto porous = coupledPorousPlane(grid);
    const auto driving = coupledDrivingPlane(grid);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithPorousInterfaces(
        grid, velocity, pressure, porous, driving, coupledSettings());
    MacVelocityField replayVelocity(grid);
    CellScalarField replayPressure(grid);
    const auto replay = projectVelocityWithPorousInterfaces(
        grid, replayVelocity, replayPressure, porous, driving,
        coupledSettings());

    constexpr double expectedVelocity = 0.4;
    double maximumVelocityError = 0.0;
    for (std::size_t face = 0; face < grid.cellCount(); ++face) {
        maximumVelocityError = std::max({
            maximumVelocityError,
            std::abs(velocity.xFaces()[face] - expectedVelocity),
            std::abs(velocity.yFaces()[face]),
            std::abs(velocity.zFaces()[face]),
        });
    }
    check(diagnostics.accepted && diagnostics.finite
              && diagnostics.nonlinearIterationCount > 1,
          "coupled porous: endpoint pressure/flux iteration converges transactionally");
    check(velocity == replayVelocity
              && pressure == replayPressure
              && diagnostics == replay,
          "coupled porous: identical nonlinear solves replay bit-for-bit");
    check(diagnostics.porousCrossingCount == 6
              && diagnostics.projection.pressureJumpFaceCount == 12
              && diagnostics.samples.size() == 6,
          "coupled porous: porous and prescribed crossings remain separately counted");
    check(maximumVelocityError < 2.0e-12,
          "coupled porous: uniform Darcy flow matches the analytic implicit endpoint");
    check(diagnostics.projection.divergenceMaximumAfterPerSecond < 1.0e-13,
          "coupled porous: accepted grid flow remains divergence-free");
    check(diagnostics.finalMaximumNormalVelocityResidualMetersPerSecond
              < 2.0e-12
              && diagnostics.finalMaximumPressureJumpResidualPascals
                  < 2.0e-11,
          "coupled porous: independent velocity and constitutive jump residuals close");
    for (const auto& sample : diagnostics.samples) {
        checkNear(sample.fluidNormalVelocityMetersPerSecond,
                  expectedVelocity, 2.0e-12,
                  "coupled porous: endpoint sample owns the committed MAC velocity");
        checkNear(sample.pressureJump.pressureJumpPascals,
                  -4.0, 2.0e-11,
                  "coupled porous: endpoint sample owns the analytic Darcy loss");
    }
    checkNear(diagnostics.totalDissipationWatts,
              9.6, 1.0e-10,
              "coupled porous: endpoint tile dissipation integrates analytically");
    const double fluidMomentum = 1.0 * 4.0 * 6.0 * expectedVelocity;
    const double pressureImpulse = (20.0 - 4.0) * 6.0 * 0.1;
    checkNear(fluidMomentum, pressureImpulse, 2.0e-15,
              "coupled porous: pressure impulse matches the uniform fluid momentum");

    MacVelocityField movingSheetVelocity(grid);
    CellScalarField movingSheetPressure(grid);
    const auto movingSheet = projectVelocityWithPorousInterfaces(
        grid, movingSheetVelocity, movingSheetPressure,
        coupledPorousPlane(grid, 0.1), driving, coupledSettings());
    check(movingSheet.accepted,
          "coupled porous: prescribed sheet-normal motion remains solvable");
    checkNear(movingSheet.samples.front()
                  .fluidNormalVelocityMetersPerSecond,
              0.42, 2.0e-12,
              "coupled porous: moving-sheet endpoint matches the analytic pressure balance");
    checkNear(movingSheet.samples.front()
                  .relativeNormalVelocityMetersPerSecond,
              0.32, 2.0e-12,
              "coupled porous: constitutive flow remains relative to sheet motion");

    MacVelocityField reverseVelocity(grid);
    CellScalarField reversePressure(grid);
    const auto reverse = projectVelocityWithPorousInterfaces(
        grid, reverseVelocity, reversePressure, porous,
        coupledDrivingPlane(grid, -20.0), coupledSettings());
    check(reverse.accepted,
          "coupled porous: reversed pressure drive converges");
    checkNear(reverse.samples.front().fluidNormalVelocityMetersPerSecond,
              -expectedVelocity, 2.0e-12,
              "coupled porous: endpoint coupling is orientation symmetric");

    auto nonlinearPorous = porous;
    for (auto& crossing : nonlinearPorous) {
        crossing.resistance = {10.0, 5.0};
    }
    MacVelocityField nonlinearVelocity(grid);
    CellScalarField nonlinearPressure(grid);
    const auto nonlinear = projectVelocityWithPorousInterfaces(
        grid, nonlinearVelocity, nonlinearPressure, nonlinearPorous,
        driving, coupledSettings());
    const double expectedNonlinearVelocity =
        (-1.25 + std::sqrt(1.8125)) / 0.25;
    check(nonlinear.accepted,
          "coupled porous: the nonlinear Darcy-Forchheimer iteration converges");
    checkNear(nonlinear.samples.front()
                  .fluidNormalVelocityMetersPerSecond,
              expectedNonlinearVelocity, 2.0e-12,
              "coupled porous: nonlinear endpoint flow matches the analytic root");

    const PeriodicCartesianGrid heterogeneousGrid(
        {4, 2, 2}, {}, {4.0, 2.0, 2.0});
    const std::vector<PorousGridFaceCrossing> heterogeneousPorous = {
        {50, 1, 2, GridFaceAxis::X, 1, 0, 0, 0.35,
         0.0, {5.0, 0.0}},
        {50, 1, 2, GridFaceAxis::X, 1, 1, 0, 0.35,
         0.0, {20.0, 0.0}},
        {50, 1, 2, GridFaceAxis::X, 1, 0, 1, 0.35,
         0.0, {5.0, 0.0}},
        {50, 1, 2, GridFaceAxis::X, 1, 1, 1, 0.35,
         0.0, {20.0, 0.0}},
    };
    MacVelocityField heterogeneousVelocity(heterogeneousGrid);
    CellScalarField heterogeneousPressure(heterogeneousGrid);
    const auto heterogeneous = projectVelocityWithPorousInterfaces(
        heterogeneousGrid, heterogeneousVelocity, heterogeneousPressure,
        heterogeneousPorous, coupledDrivingPlane(heterogeneousGrid),
        coupledSettings());
    check(heterogeneous.accepted
              && heterogeneous.samples.size() == 4,
          "coupled porous: heterogeneous resolved tiles converge in one grid solve");
    if (heterogeneous.samples.size() == 4) {
        check(heterogeneous.samples[0].fluidNormalVelocityMetersPerSecond
                  > heterogeneous.samples[1]
                        .fluidNormalVelocityMetersPerSecond
                  && heterogeneous.samples[1]
                         .fluidNormalVelocityMetersPerSecond > 0.0
                  && heterogeneous.samples[2]
                         .fluidNormalVelocityMetersPerSecond
                      > heterogeneous.samples[3]
                            .fluidNormalVelocityMetersPerSecond,
              "coupled porous: the lower-resistance tile carries the larger forward flux");
    }
    for (const auto& sample : heterogeneous.samples) {
        checkNear(
            sample.pressureJump.pressureJumpPascals,
            porousPressureJumpPascals(
                heterogeneousPorous[sample.pressureJump.j].resistance,
                sample.relativeNormalVelocityMetersPerSecond),
            0.0,
            "coupled porous: every heterogeneous endpoint sample obeys its own law");
    }
    check(heterogeneous.projection.divergenceMaximumAfterPerSecond
              < 2.0e-12,
          "coupled porous: heterogeneous accepted flow remains divergence-free");
}

void testCoupledPorousDelegationValidationAndRollback() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    MacVelocityField directVelocity(grid);
    for (std::size_t face = 0; face < grid.cellCount(); ++face) {
        directVelocity.xFaces()[face] =
            0.1 * std::sin(static_cast<double>(face + 1));
        directVelocity.yFaces()[face] =
            0.1 * std::cos(static_cast<double>(face + 1));
    }
    MacVelocityField delegatedVelocity = directVelocity;
    CellScalarField directPressure(grid, 3.0);
    CellScalarField delegatedPressure = directPressure;
    const auto settings = coupledSettings();
    const auto direct = projectVelocity(
        grid, directVelocity, directPressure, settings.projection);
    const auto delegated = projectVelocityWithPorousInterfaces(
        grid, delegatedVelocity, delegatedPressure, {}, settings);
    check(directVelocity == delegatedVelocity
              && directPressure == delegatedPressure
              && direct == delegated.projection,
          "coupled porous: empty topology delegates to the exact base projection");

    MacVelocityField failedVelocity(grid);
    CellScalarField failedPressure(grid, 7.0);
    const auto originalVelocity = failedVelocity;
    const auto originalPressure = failedPressure;
    auto truncated = settings;
    truncated.maximumNonlinearIterations = 1;
    truncated.absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-30;
    truncated.relativeNormalVelocityTolerance = 0.0;
    truncated.absolutePressureJumpTolerancePascals = 1.0e-30;
    truncated.relativePressureJumpTolerance = 0.0;
    const auto failed = projectVelocityWithPorousInterfaces(
        grid, failedVelocity, failedPressure, coupledPorousPlane(grid),
        coupledDrivingPlane(grid), truncated);
    check(!failed.accepted && failed.nonlinearIterationCount == 1,
          "coupled porous: an exhausted nonlinear iteration is rejected");
    check(failedVelocity == originalVelocity
              && failedPressure == originalPressure,
          "coupled porous: nonlinear failure commits neither field");

    auto invalid = settings;
    invalid.relaxation = 0.0;
    expectInvalid(
        [&] {
            static_cast<void>(projectVelocityWithPorousInterfaces(
                grid, failedVelocity, failedPressure,
                coupledPorousPlane(grid), invalid));
        },
        "coupled porous: invalid relaxation is rejected");
    check(failedVelocity == originalVelocity
              && failedPressure == originalPressure,
          "coupled porous: invalid settings are rejected before mutation");

    const PeriodicCartesianGrid foreignGrid(
        {5, 3, 2}, {}, {5.0, 3.0, 2.0});
    const SharpPressureJumpField foreign(foreignGrid);
    expectInvalid(
        [&] {
            static_cast<void>(projectVelocityWithPorousInterfaces(
                grid, failedVelocity, failedPressure,
                coupledPorousPlane(grid), foreign, settings));
        },
        "coupled porous: a prescribed jump field from another grid is rejected");
}

void testAxisAreaAndValidation() {
    const PeriodicCartesianGrid grid(
        {2, 2, 2}, {}, {4.0, 6.0, 8.0});
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 3.0);
    std::ranges::fill(velocity.yFaces(), 4.0);
    std::ranges::fill(velocity.zFaces(), 5.0);
    const DarcyForchheimerResistance linear{10.0, 0.0};
    const PorousPressureJumpField axes(grid, velocity, {
        {10, 1, 2, GridFaceAxis::X, 0, 0, 0, 0.2, 0.0, linear},
        {20, 3, 4, GridFaceAxis::Y, 1, 0, 0, 0.4, 0.0, linear},
        {30, 5, 6, GridFaceAxis::Z, 1, 1, 0, 0.6, 0.0, linear},
    });
    check(axes.samples().size() == 3
              && axes.samples()[0].faceAreaSquareMeters == 12.0
              && axes.samples()[1].faceAreaSquareMeters == 8.0
              && axes.samples()[2].faceAreaSquareMeters == 6.0
              && axes.samples()[0].volumeFlowRateCubicMetersPerSecond == 36.0
              && axes.samples()[1].volumeFlowRateCubicMetersPerSecond == 32.0
              && axes.samples()[2].volumeFlowRateCubicMetersPerSecond == 30.0,
          "porous sampling uses the correct X/Y/Z MAC tile areas");

    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals({}, 1.0)); },
        "porous law rejects zero resistance");
    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals(
            {-1.0, 1.0}, 1.0)); },
        "porous law rejects negative resistance");
    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals(
            linear, std::numeric_limits<double>::infinity())); },
        "porous law rejects non-finite velocity");
    expectInvalid(
        [&] { static_cast<void>(
            porousRelativeNormalVelocityMetersPerSecond(
                linear, std::numeric_limits<double>::quiet_NaN())); },
        "porous inverse rejects non-finite pressure");

    const PeriodicCartesianGrid otherGrid(
        {3, 2, 2}, {}, {6.0, 6.0, 8.0});
    const MacVelocityField foreignVelocity(otherGrid);
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, foreignVelocity, {})); },
        "porous field rejects foreign MAC topology");
    auto invalid = PorousGridFaceCrossing{
        10, 1, 2, GridFaceAxis::X, 2, 0, 0,
        0.5, 0.0, linear};
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field rejects an out-of-range face");
    invalid.i = 0;
    invalid.surfaceNormalVelocityMetersPerSecond =
        std::numeric_limits<double>::infinity();
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field rejects non-finite sheet velocity");
    invalid.surfaceNormalVelocityMetersPerSecond = 0.0;
    invalid.surfaceStableId = 0;
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field delegates stable topology validation");
}

PorousPlugFlowSettings plugSettings() {
    PorousPlugFlowSettings settings;
    settings.resistance = {100.0, 25.0};
    settings.densityKgPerCubicMeter = 1.2;
    settings.flowLengthMeters = 4.0;
    settings.crossSectionAreaSquareMeters = 3.0;
    settings.drivingPressureRisePascals = 250.0;
    settings.timeStepSeconds = 1.0 / 60.0;
    return settings;
}

void testPressureDrivenPlugFlow() {
    const auto settings = plugSettings();
    double velocity = 0.0;
    const auto first = advancePorousPlugFlow(velocity, settings);
    check(first.accepted
              && velocity == first.velocityAfterMetersPerSecond
              && first.midpointVelocityMetersPerSecond
                  == 0.5 * first.velocityAfterMetersPerSecond,
          "porous plug flow commits its implicit midpoint candidate");
    checkNear(first.momentumResidualNewtonSeconds,
              0.0, 4.0e-15,
              "porous plug flow closes pressure impulse");
    checkNear(first.energyResidualJoules,
              0.0, 4.0e-15,
              "porous plug flow closes pressure work and dissipation");
    check(first.porousDissipationJoules > 0.0
              && first.drivingPressureWorkJoules
                  > first.porousDissipationJoules,
          "accelerating porous plug retains positive dissipation");

    for (std::size_t step = 1; step < 600; ++step) {
        static_cast<void>(advancePorousPlugFlow(velocity, settings));
    }
    const double steady = 0.5 * (std::sqrt(56.0) - 4.0);
    checkNear(velocity, steady, 2.0e-15,
              "porous plug flow converges to the analytic pressure-driven speed");
    const auto steadyStep = advancePorousPlugFlow(velocity, settings);
    checkNear(steadyStep.endpointPressureDropPascals,
              settings.drivingPressureRisePascals, 3.0e-13,
              "porous steady pressure loss balances the driving rise");
    checkNear(steadyStep.drivingPressureWorkJoules,
              steadyStep.porousDissipationJoules, 4.0e-14,
              "porous steady pressure work becomes material dissipation");

    auto reverseSettings = settings;
    reverseSettings.drivingPressureRisePascals *= -1.0;
    double reverseVelocity = 0.0;
    for (std::size_t step = 0; step < 600; ++step) {
        static_cast<void>(advancePorousPlugFlow(
            reverseVelocity, reverseSettings));
    }
    checkNear(reverseVelocity, -steady, 2.0e-15,
              "porous plug flow is orientation symmetric");

    auto decaySettings = settings;
    decaySettings.drivingPressureRisePascals = 0.0;
    double decayingVelocity = 2.0;
    const auto decay = advancePorousPlugFlow(
        decayingVelocity, decaySettings);
    check(decayingVelocity > 0.0 && decayingVelocity < 2.0
              && decay.drivingPressureWorkJoules == 0.0
              && decay.porousDissipationJoules > 0.0,
          "unforced porous plug flow decelerates dissipatively");
    checkNear(decay.kineticEnergyAfterJoules
                  - decay.kineticEnergyBeforeJoules,
              -decay.porousDissipationJoules, 4.0e-14,
              "unforced porous loss becomes the exact kinetic-energy change");

    auto invalidSettings = settings;
    invalidSettings.timeStepSeconds = 0.0;
    const double preserved = velocity;
    expectInvalid(
        [&] { static_cast<void>(advancePorousPlugFlow(
            velocity, invalidSettings)); },
        "porous plug flow rejects a nonpositive time step");
    check(velocity == preserved,
          "rejected porous plug flow leaves caller state unchanged");
}

} // namespace

int main() {
    testConstitutiveLawAndInverse();
    testFluxDrivenProjectionCanonical();
    testImplicitGridPressureFluxCoupling();
    testCoupledPorousDelegationValidationAndRollback();
    testAxisAreaAndValidation();
    testPressureDrivenPlugFlow();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing porous-interface check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing porous-interface checks passed");
    return 0;
}
