#include "fluid/advection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::UniformMacAdvectionSettings;
using simwing::fsi::fluid::advectVelocityByUniformFlow;
using simwing::fsi::fluid::computeDivergence;
using simwing::fsi::fluid::maximumAbsoluteValue;

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

double componentSum(const std::span<const double> values) {
    double result = 0.0;
    for (const double value : values) {
        result += value;
    }
    return result;
}

MacVelocityField deterministicVelocity(
    const PeriodicCartesianGrid& grid) {
    MacVelocityField result(grid);
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        const double sample = static_cast<double>(index + 1);
        result.xFaces()[index] = 0.5 + std::sin(0.37 * sample);
        result.yFaces()[index] = -0.3 + std::cos(0.29 * sample + 0.2);
        result.zFaces()[index] = 0.1 + std::sin(0.19 * sample - 0.4);
    }
    return result;
}

UniformMacAdvectionSettings settings() {
    UniformMacAdvectionSettings result;
    result.densityKgPerCubicMeter = 1.2;
    result.transportVelocityMetersPerSecond = {0.2, -0.1, 0.05};
    result.timeStepSeconds = 0.1;
    result.absoluteMomentumToleranceNewtonSeconds = 3.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 3.0e-12;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

void testZeroTransportAndUniformModes() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto arbitrary = deterministicVelocity(grid);
    const auto original = arbitrary;
    auto zeroSettings = settings();
    zeroSettings.transportVelocityMetersPerSecond = {};
    const auto zero = advectVelocityByUniformFlow(
        grid, arbitrary, zeroSettings);
    check(zero.accepted && zero.stable && zero.bounded && zero.finite
              && zero.totalAbsoluteCourantNumber == 0.0
              && zero.maximumVelocityChangeMetersPerSecond == 0.0
              && arbitrary == original,
          "zero transport: exact no-op commits bit-identically");

    MacVelocityField uniform(grid);
    std::ranges::fill(uniform.xFaces(), 1.25);
    std::ranges::fill(uniform.yFaces(), -0.5);
    std::ranges::fill(uniform.zFaces(), 0.125);
    const auto uniformBefore = uniform;
    const auto transported = advectVelocityByUniformFlow(
        grid, uniform, settings());
    check(transported.accepted && transported.bounded
              && uniform == uniformBefore
              && transported.maximumVelocityChangeMetersPerSecond == 0.0
              && transported.numericalKineticEnergyLossJoules == 0.0,
          "uniform mode: conservative transport preserves constant velocity exactly");
}

void testExactOneCellPeriodicShift() {
    const PeriodicCartesianGrid grid(
        {8, 4, 3}, {}, {4.0, 2.0, 1.5});
    auto velocity = deterministicVelocity(grid);
    const auto before = velocity;
    auto shiftSettings = settings();
    shiftSettings.timeStepSeconds = 0.25;
    shiftSettings.transportVelocityMetersPerSecond = {2.0, 0.0, 0.0};
    const auto diagnostics = advectVelocityByUniformFlow(
        grid, velocity, shiftSettings);
    const auto counts = grid.cellCounts();
    double maximumError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t previousI = i == 0 ? counts.x - 1 : i - 1;
                const std::size_t destination = grid.cellIndex(i, j, k);
                const std::size_t source = grid.cellIndex(previousI, j, k);
                maximumError = std::max({
                    maximumError,
                    std::abs(velocity.xFaces()[destination]
                             - before.xFaces()[source]),
                    std::abs(velocity.yFaces()[destination]
                             - before.yFaces()[source]),
                    std::abs(velocity.zFaces()[destination]
                             - before.zFaces()[source]),
                });
            }
        }
    }
    check(diagnostics.accepted && diagnostics.stable
              && diagnostics.bounded
              && diagnostics.totalAbsoluteCourantNumber == 1.0,
          "exact shift: the sharp convex CFL boundary is accepted");
    check(maximumError == 0.0
              && diagnostics.maximumBoundViolationMetersPerSecond == 0.0
              && std::abs(diagnostics.numericalKineticEnergyLossJoules)
                  < 2.0e-12,
          "exact shift: CFL-one translation wraps every MAC component exactly");
}

void testDeterminismConservationAndBounds() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto first = deterministicVelocity(grid);
    auto second = first;
    const double xSum = componentSum(first.xFaces());
    const double ySum = componentSum(first.yFaces());
    const double zSum = componentSum(first.zFaces());
    const auto firstDiagnostics = advectVelocityByUniformFlow(
        grid, first, settings());
    const auto secondDiagnostics = advectVelocityByUniformFlow(
        grid, second, settings());
    check(first == second && firstDiagnostics == secondDiagnostics,
          "determinism: identical transport inputs replay bit-for-bit");
    checkNear(componentSum(first.xFaces()), xSum, 1.0e-13,
              "conservation: periodic transport preserves x momentum sum");
    checkNear(componentSum(first.yFaces()), ySum, 1.0e-13,
              "conservation: periodic transport preserves y momentum sum");
    checkNear(componentSum(first.zFaces()), zSum, 1.0e-13,
              "conservation: periodic transport preserves z momentum sum");
    check(firstDiagnostics.accepted && firstDiagnostics.bounded
              && firstDiagnostics.maximumBoundViolationMetersPerSecond == 0.0
              && firstDiagnostics.momentumResidualNormNewtonSeconds
                  < 3.0e-12
              && firstDiagnostics.kineticEnergyAfterJoules
                  <= firstDiagnostics.kineticEnergyBeforeJoules
              && firstDiagnostics.numericalKineticEnergyLossJoules > 0.0,
          "conservation: mixed-sign transport closes bounds, momentum, and energy");
}

void testSolenoidalTransport() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {24, 24, 4}, {}, {twoPi, twoPi, twoPi});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                velocity.xFaces()[index] =
                    std::sin(xFace.x) * std::cos(xFace.y);
                velocity.yFaces()[index] =
                    -std::cos(yFace.x) * std::sin(yFace.y);
            }
        }
    }
    CellScalarField before(grid);
    computeDivergence(grid, velocity, before);
    auto transportSettings = settings();
    transportSettings.transportVelocityMetersPerSecond = {0.3, -0.2, 0.1};
    transportSettings.timeStepSeconds = 0.1;
    const auto diagnostics = advectVelocityByUniformFlow(
        grid, velocity, transportSettings);
    CellScalarField after(grid);
    computeDivergence(grid, velocity, after);
    check(diagnostics.accepted,
          "solenoidal mode: transport step is accepted");
    checkNear(maximumAbsoluteValue(before), 0.0, 2.0e-14,
              "solenoidal mode: initial MAC field is divergence-free");
    checkNear(maximumAbsoluteValue(after), 0.0, 3.0e-14,
              "solenoidal mode: uniform transport commutes with divergence");
}

double fullPeriodError(const std::size_t resolution) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, 2, 2}, {}, {twoPi, 1.0, 1.0});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                velocity.xFaces()[grid.cellIndex(i, j, k)] = std::sin(
                    grid.xFaceCenterMeters(i, j, k).x);
            }
        }
    }
    const auto expected = velocity;
    auto transportSettings = settings();
    transportSettings.transportVelocityMetersPerSecond = {1.0, 0.0, 0.0};
    transportSettings.timeStepSeconds =
        0.5 * grid.cellSpacingMeters().x;
    const std::size_t steps = 2 * resolution;
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advectVelocityByUniformFlow(
            grid, velocity, transportSettings);
        check(diagnostics.accepted,
              "convergence: every full-period transport step is accepted");
    }
    double squaredError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        const double error = velocity.xFaces()[index]
            - expected.xFaces()[index];
        squaredError += error * error;
    }
    return std::sqrt(
        squaredError / static_cast<double>(grid.cellCount()));
}

void testFirstOrderConvergence() {
    const double coarseError = fullPeriodError(16);
    const double mediumError = fullPeriodError(32);
    const double fineError = fullPeriodError(64);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 1.7 && coarseRatio < 2.2,
          "convergence: first donor-cell refinement approaches first order");
    check(fineRatio > 1.8 && fineRatio < 2.1,
          "convergence: second donor-cell refinement approaches first order");
}

void testTransactionalRejection() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    const auto original = deterministicVelocity(grid);
    auto unstable = original;
    auto unstableSettings = settings();
    unstableSettings.transportVelocityMetersPerSecond = {10.0, 10.0, 10.0};
    unstableSettings.timeStepSeconds = 1.0;
    const auto rejected = advectVelocityByUniformFlow(
        grid, unstable, unstableSettings);
    check(!rejected.stable && !rejected.accepted && unstable == original,
          "stability: excessive total CFL is rejected transactionally");

    auto invalid = original;
    auto invalidSettings = settings();
    invalidSettings.maximumTotalCourantNumber = 1.01;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByUniformFlow(
            grid, invalid, invalidSettings)); },
        "validation: an unsafe configured CFL ceiling is rejected");
    check(invalid == original,
          "validation: invalid settings leave velocity unchanged");

    auto nonfinite = original;
    nonfinite.zFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(advectVelocityByUniformFlow(
            grid, nonfinite, settings())); },
        "validation: non-finite velocity is rejected");
    check(std::isnan(nonfinite.zFaces().front()),
          "validation: non-finite rejection does not rewrite its input");

    auto nonfiniteTransport = settings();
    nonfiniteTransport.transportVelocityMetersPerSecond.x =
        std::numeric_limits<double>::infinity();
    auto finiteInput = original;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByUniformFlow(
            grid, finiteInput, nonfiniteTransport)); },
        "validation: non-finite transport velocity is rejected");
    check(finiteInput == original,
          "validation: invalid transport is transactional");

    const PeriodicCartesianGrid foreignGrid(
        {8, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto mismatched = original;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByUniformFlow(
            foreignGrid, mismatched, settings())); },
        "validation: a field from another grid is rejected");
    check(mismatched == original,
          "validation: grid mismatch is transactional");
}

} // namespace

int main() {
    testZeroTransportAndUniformModes();
    testExactOneCellPeriodicShift();
    testDeterminismConservationAndBounds();
    testSolenoidalTransport();
    testFirstOrderConvergence();
    testTransactionalRejection();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid advection check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid advection checks passed");
    return 0;
}
