#include "fluid/diffusion.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::GridCellCounts;
using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PeriodicMacDiffusionSettings;
using simwing::fsi::fluid::diffuseVelocityExplicit;
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
        result.xFaces()[index] = 0.4 + std::sin(0.31 * sample);
        result.yFaces()[index] = -0.2 + std::cos(0.23 * sample + 0.1);
        result.zFaces()[index] = 0.1 + std::sin(0.17 * sample - 0.2);
    }
    return result;
}

PeriodicMacDiffusionSettings settings() {
    PeriodicMacDiffusionSettings result;
    result.densityKgPerCubicMeter = 1.2;
    result.kinematicViscositySquareMetersPerSecond = 0.02;
    result.timeStepSeconds = 0.01;
    result.absoluteMomentumToleranceNewtonSeconds = 2.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 2.0e-12;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

void testZeroViscosityAndUniformModes() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto arbitrary = deterministicVelocity(grid);
    const auto original = arbitrary;
    auto zeroSettings = settings();
    zeroSettings.kinematicViscositySquareMetersPerSecond = 0.0;
    const auto zero = diffuseVelocityExplicit(
        grid, arbitrary, zeroSettings);
    check(zero.accepted && zero.stable && zero.finite
              && zero.totalDiffusionNumber == 0.0
              && zero.maximumVelocityChangeMetersPerSecond == 0.0
              && arbitrary == original,
          "zero viscosity: exact no-op commits bit-identically");

    MacVelocityField uniform(grid);
    std::ranges::fill(uniform.xFaces(), 1.25);
    std::ranges::fill(uniform.yFaces(), -0.5);
    std::ranges::fill(uniform.zFaces(), 0.125);
    const auto uniformBefore = uniform;
    const auto diffused = diffuseVelocityExplicit(
        grid, uniform, settings());
    check(diffused.accepted && uniform == uniformBefore
              && diffused.maximumVelocityChangeMetersPerSecond == 0.0
              && diffused.dissipatedKineticEnergyJoules == 0.0,
          "uniform mode: periodic viscosity preserves constant velocity exactly");
}

void testAnalyticDiscreteFourierDecay() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {16, 12, 8}, {}, {twoPi, twoPi, twoPi});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                const auto zFace = grid.zFaceCenterMeters(i, j, k);
                velocity.xFaces()[index] =
                    0.4 + std::sin(xFace.y);
                velocity.yFaces()[index] =
                    -0.2 + 0.5 * std::cos(2.0 * yFace.x);
                velocity.zFaces()[index] =
                    0.1 + 0.25 * std::sin(zFace.x)
                        * std::cos(zFace.y);
            }
        }
    }
    const auto before = velocity;
    CellScalarField divergenceBefore(grid);
    computeDivergence(grid, velocity, divergenceBefore);
    const auto diffusionSettings = settings();
    const auto diagnostics = diffuseVelocityExplicit(
        grid, velocity, diffusionSettings);
    const auto spacing = grid.cellSpacingMeters();
    const double lambdaX1 = 4.0
        * std::pow(std::sin(std::numbers::pi
                           / static_cast<double>(counts.x)), 2)
        / (spacing.x * spacing.x);
    const double lambdaX2 = 4.0
        * std::pow(std::sin(2.0 * std::numbers::pi
                           / static_cast<double>(counts.x)), 2)
        / (spacing.x * spacing.x);
    const double lambdaY1 = 4.0
        * std::pow(std::sin(std::numbers::pi
                           / static_cast<double>(counts.y)), 2)
        / (spacing.y * spacing.y);
    const double viscosityTime = diffusionSettings
        .kinematicViscositySquareMetersPerSecond
        * diffusionSettings.timeStepSeconds;
    const double xFactor = 1.0 - viscosityTime * lambdaY1;
    const double yFactor = 1.0 - viscosityTime * lambdaX2;
    const double zFactor = 1.0
        - viscosityTime * (lambdaX1 + lambdaY1);
    double maximumError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        maximumError = std::max({
            maximumError,
            std::abs(velocity.xFaces()[index]
                     - (0.4 + xFactor
                         * (before.xFaces()[index] - 0.4))),
            std::abs(velocity.yFaces()[index]
                     - (-0.2 + yFactor
                         * (before.yFaces()[index] + 0.2))),
            std::abs(velocity.zFaces()[index]
                     - (0.1 + zFactor
                         * (before.zFaces()[index] - 0.1))),
        });
    }
    check(diagnostics.accepted && diagnostics.stable
              && diagnostics.totalDiffusionNumber < 0.01,
          "Fourier decay: stable explicit step is accepted");
    check(maximumError < 8.0e-16,
          "Fourier decay: every staggered mode follows its discrete eigenvalue");
    check(diagnostics.momentumResidualNormNewtonSeconds < 2.0e-12
              && diagnostics.kineticEnergyAfterJoules
                  < diagnostics.kineticEnergyBeforeJoules
              && diagnostics.dissipatedKineticEnergyJoules > 0.0,
          "Fourier decay: momentum closes and kinetic energy dissipates");
    CellScalarField divergenceAfter(grid);
    computeDivergence(grid, velocity, divergenceAfter);
    check(maximumAbsoluteValue(divergenceBefore) < 2.0e-14
              && maximumAbsoluteValue(divergenceAfter) < 2.0e-14,
          "Fourier decay: componentwise viscosity preserves the solenoidal mode");
}

double discreteEigenvalueError(const std::size_t resolution) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, 2, 2}, {}, {twoPi, 1.0, 1.0});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                velocity.xFaces()[index] = std::sin(
                    grid.xFaceCenterMeters(i, j, k).x);
            }
        }
    }
    const double before = velocity.xFaces()[grid.cellIndex(1, 0, 0)];
    auto diffusionSettings = settings();
    diffusionSettings.kinematicViscositySquareMetersPerSecond = 0.1;
    diffusionSettings.timeStepSeconds = 1.0e-4;
    const auto diagnostics = diffuseVelocityExplicit(
        grid, velocity, diffusionSettings);
    check(diagnostics.accepted,
          "convergence: manufactured Fourier step is accepted");
    const double after = velocity.xFaces()[grid.cellIndex(1, 0, 0)];
    const double measuredEigenvalue = (1.0 - after / before)
        / (diffusionSettings.kinematicViscositySquareMetersPerSecond
           * diffusionSettings.timeStepSeconds);
    return std::abs(measuredEigenvalue - 1.0);
}

void testSecondOrderSpatialConvergence() {
    const double coarseError = discreteEigenvalueError(16);
    const double mediumError = discreteEigenvalueError(32);
    const double fineError = discreteEigenvalueError(64);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 3.95 && coarseRatio < 4.05,
          "convergence: first viscous eigenvalue refinement is second order");
    check(fineRatio > 3.98 && fineRatio < 4.02,
          "convergence: second viscous eigenvalue refinement is second order");
}

void testConservationDeterminismAndRejection() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto first = deterministicVelocity(grid);
    auto second = first;
    const auto original = first;
    const double xSum = componentSum(first.xFaces());
    const double ySum = componentSum(first.yFaces());
    const double zSum = componentSum(first.zFaces());
    const auto firstDiagnostics = diffuseVelocityExplicit(
        grid, first, settings());
    const auto secondDiagnostics = diffuseVelocityExplicit(
        grid, second, settings());
    check(first == second && firstDiagnostics == secondDiagnostics,
          "determinism: identical viscous inputs replay bit-for-bit");
    checkNear(componentSum(first.xFaces()), xSum, 1.0e-13,
              "conservation: periodic viscosity preserves x momentum sum");
    checkNear(componentSum(first.yFaces()), ySum, 1.0e-13,
              "conservation: periodic viscosity preserves y momentum sum");
    checkNear(componentSum(first.zFaces()), zSum, 1.0e-13,
              "conservation: periodic viscosity preserves z momentum sum");
    check(firstDiagnostics.accepted
              && firstDiagnostics.momentumResidualNormNewtonSeconds
                  < 2.0e-12
              && firstDiagnostics.dissipatedKineticEnergyJoules > 0.0,
          "conservation: accepted viscous step closes physical ledgers");

    auto unstable = original;
    auto unstableSettings = settings();
    unstableSettings.kinematicViscositySquareMetersPerSecond = 10.0;
    unstableSettings.timeStepSeconds = 1.0;
    const auto rejected = diffuseVelocityExplicit(
        grid, unstable, unstableSettings);
    check(!rejected.stable && !rejected.accepted && unstable == original,
          "stability: excessive explicit diffusion is rejected transactionally");

    auto invalid = original;
    auto invalidSettings = settings();
    invalidSettings.maximumDiffusionNumber = 0.51;
    expectRejected(
        [&] { static_cast<void>(diffuseVelocityExplicit(
            grid, invalid, invalidSettings)); },
        "validation: an unsafe configured stability ceiling is rejected");
    check(invalid == original,
          "validation: invalid settings leave velocity unchanged");

    auto nonfinite = original;
    nonfinite.yFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(diffuseVelocityExplicit(
            grid, nonfinite, settings())); },
        "validation: non-finite velocity is rejected");
    check(std::isnan(nonfinite.yFaces().front()),
          "validation: non-finite rejection does not rewrite its input");

    const PeriodicCartesianGrid foreignGrid(
        {8, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto mismatched = original;
    expectRejected(
        [&] { static_cast<void>(diffuseVelocityExplicit(
            foreignGrid, mismatched, settings())); },
        "validation: a field from another grid is rejected");
    check(mismatched == original,
          "validation: grid mismatch is transactional");
}

void testSharpStabilityBoundary() {
    const PeriodicCartesianGrid grid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double value = (i + j + k) % 2 == 0 ? 1.0 : -1.0;
                const std::size_t index = grid.cellIndex(i, j, k);
                velocity.xFaces()[index] = value;
                velocity.yFaces()[index] = value;
                velocity.zFaces()[index] = value;
            }
        }
    }
    const auto before = velocity;
    auto boundarySettings = settings();
    boundarySettings.kinematicViscositySquareMetersPerSecond = 1.0 / 6.0;
    boundarySettings.timeStepSeconds = 1.0;
    const auto diagnostics = diffuseVelocityExplicit(
        grid, velocity, boundarySettings);
    double maximumSignFlipError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        maximumSignFlipError = std::max({
            maximumSignFlipError,
            std::abs(velocity.xFaces()[index] + before.xFaces()[index]),
            std::abs(velocity.yFaces()[index] + before.yFaces()[index]),
            std::abs(velocity.zFaces()[index] + before.zFaces()[index]),
        });
    }
    check(diagnostics.accepted && diagnostics.stable
              && diagnostics.totalDiffusionNumber == 0.5,
          "stability boundary: the non-amplifying CFL limit is accepted");
    check(maximumSignFlipError < 3.0e-16
              && std::abs(diagnostics.dissipatedKineticEnergyJoules)
                  < 2.0e-12,
          "stability boundary: Nyquist mode flips sign without energy growth");
}

} // namespace

int main() {
    testZeroViscosityAndUniformModes();
    testAnalyticDiscreteFourierDecay();
    testSecondOrderSpatialConvergence();
    testConservationDeterminismAndRejection();
    testSharpStabilityBoundary();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid diffusion check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid diffusion checks passed");
    return 0;
}
