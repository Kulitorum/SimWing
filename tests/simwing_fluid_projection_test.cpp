#include "fluid/grid.h"
#include "fluid/projection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::GridCellCounts;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::ProjectionSettings;
using simwing::fsi::fluid::Vector3;
using simwing::fsi::fluid::applyNegativeLaplacian;
using simwing::fsi::fluid::computeDivergence;
using simwing::fsi::fluid::computePressureGradient;
using simwing::fsi::fluid::mean;
using simwing::fsi::fluid::projectVelocity;

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

PeriodicCartesianGrid unitGrid(const GridCellCounts counts) {
    return PeriodicCartesianGrid(counts, {}, {2.0, 3.0, 4.0});
}

double componentSum(const std::span<const double> values) {
    double result = 0.0;
    for (const double value : values) {
        result += value;
    }
    return result;
}

void testGridContractAndValidation() {
    const PeriodicCartesianGrid grid({4, 3, 2}, {-1.0, 2.0, 4.0}, {1.0, 5.0, 8.0});
    check(grid.cellCount() == 24, "grid: cell count is the product of dimensions");
    checkNear(grid.cellSpacingMeters().x, 0.5, 0.0,
              "grid: x spacing follows physical bounds");
    checkNear(grid.cellSpacingMeters().y, 1.0, 0.0,
              "grid: y spacing follows physical bounds");
    checkNear(grid.cellSpacingMeters().z, 2.0, 0.0,
              "grid: z spacing follows physical bounds");
    check(grid.cellIndex(3, 2, 1) == 23,
          "grid: indexing is deterministic with x as the fast axis");
    checkNear(grid.xFaceCenterMeters(0, 0, 0).x, -1.0, 0.0,
              "grid: x velocity is stored on the lower x face");
    checkNear(grid.yFaceCenterMeters(0, 0, 0).x, -0.75, 0.0,
              "grid: y velocity is centred in x");
    checkNear(grid.cellCenterMeters(0, 0, 0).z, 5.0, 0.0,
              "grid: pressure is cell centred");

    bool rejected = false;
    try {
        static_cast<void>(PeriodicCartesianGrid({1, 2, 2}, {}, {1.0, 1.0, 1.0}));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, "grid: periodic axes require at least two cells");

    rejected = false;
    try {
        static_cast<void>(PeriodicCartesianGrid(
            {2, 2, 2}, {}, {1.0, 0.0, 1.0}));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, "grid: degenerate physical bounds are rejected");
}

void testDiscreteOperatorIdentity() {
    const auto grid = unitGrid({7, 5, 4});
    CellScalarField scalar(grid);
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                scalar.values()[index] = std::sin(0.7 * static_cast<double>(i + 1))
                    + 0.3 * std::cos(0.4 * static_cast<double>(j + 2))
                    - 0.2 * std::sin(0.9 * static_cast<double>(k + 3));
                velocity.xFaces()[index] =
                    std::cos(0.31 * static_cast<double>(index + 1));
                velocity.yFaces()[index] =
                    std::sin(0.23 * static_cast<double>(index + 2));
                velocity.zFaces()[index] =
                    std::cos(0.17 * static_cast<double>(index + 5));
            }
        }
    }

    CellScalarField divergence(grid);
    MacVelocityField gradient(grid);
    computeDivergence(grid, velocity, divergence);
    computePressureGradient(grid, scalar, gradient);
    double scalarDivergenceIntegral = 0.0;
    double gradientVelocityIntegral = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        scalarDivergenceIntegral +=
            scalar.values()[index] * divergence.values()[index];
        gradientVelocityIntegral +=
            gradient.xFaces()[index] * velocity.xFaces()[index]
            + gradient.yFaces()[index] * velocity.yFaces()[index]
            + gradient.zFaces()[index] * velocity.zFaces()[index];
    }
    const double volume = grid.cellVolumeCubicMeters();
    checkNear(volume * scalarDivergenceIntegral,
              -volume * gradientVelocityIntegral,
              2.0e-14,
              "operators: gradient is the negative adjoint of divergence");

    CellScalarField negativeLaplacian(grid);
    CellScalarField divergenceOfGradient(grid);
    applyNegativeLaplacian(grid, scalar, negativeLaplacian);
    computeDivergence(grid, gradient, divergenceOfGradient);
    double maximumStencilDifference = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        maximumStencilDifference = std::max(
            maximumStencilDifference,
            std::abs(negativeLaplacian.values()[index]
                     + divergenceOfGradient.values()[index]));
    }
    check(maximumStencilDifference < 2.0e-14,
          "operators: negative Laplacian exactly composes divergence and gradient");
}

ProjectionSettings strictSettings() {
    ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.025;
    settings.absoluteResidualTolerance = 1.0e-11;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    return settings;
}

void testDiscretelyManufacturedPressureProjection() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 10, 8}, {}, {twoPi, twoPi, twoPi});
    CellScalarField expectedPressure(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const auto center = grid.cellCenterMeters(i, j, k);
                expectedPressure.values()[index] =
                    2.0 * std::sin(center.x)
                    + 0.5 * std::cos(2.0 * center.y)
                    - 0.25 * std::sin(3.0 * center.z);
            }
        }
    }
    MacVelocityField pressureGradient(grid);
    computePressureGradient(grid, expectedPressure, pressureGradient);
    const auto settings = strictSettings();
    MacVelocityField predictedVelocity(grid);
    const double scale = settings.timeStepSeconds
        / settings.densityKgPerCubicMeter;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        predictedVelocity.xFaces()[index] =
            1.25 + scale * pressureGradient.xFaces()[index];
        predictedVelocity.yFaces()[index] =
            -0.5 + scale * pressureGradient.yFaces()[index];
        predictedVelocity.zFaces()[index] =
            0.75 + scale * pressureGradient.zFaces()[index];
    }
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocity(
        grid, predictedVelocity, pressure, settings);

    check(diagnostics.converged,
          "manufactured projection: conjugate gradient converges");
    check(diagnostics.iterationCount > 0,
          "manufactured projection: a nonzero pressure correction is solved");
    check(diagnostics.divergenceL2AfterPerSecond < 2.0e-12,
          "manufactured projection: corrected divergence reaches roundoff");
    check(diagnostics.divergenceL2AfterPerSecond
              < 1.0e-10 * diagnostics.divergenceL2BeforePerSecond,
          "manufactured projection: divergence falls by ten orders of magnitude");
    checkNear(mean(pressure), 0.0, 2.0e-16,
              "manufactured projection: periodic pressure uses a zero-mean gauge");
    double maximumPressureError = 0.0;
    double maximumVelocityError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        maximumPressureError = std::max(
            maximumPressureError,
            std::abs(pressure.values()[index]
                     - expectedPressure.values()[index]));
        maximumVelocityError = std::max({
            maximumVelocityError,
            std::abs(predictedVelocity.xFaces()[index] - 1.25),
            std::abs(predictedVelocity.yFaces()[index] + 0.5),
            std::abs(predictedVelocity.zFaces()[index] - 0.75),
        });
    }
    check(maximumPressureError < 2.0e-12,
          "manufactured projection: discrete pressure is recovered up to its gauge");
    check(maximumVelocityError < 2.0e-14,
          "manufactured projection: only the solenoidal velocity remains");
}

void testTaylorGreenIsPreserved() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {32, 32, 2}, {}, {twoPi, twoPi, twoPi});
    MacVelocityField velocity(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                velocity.xFaces()[index] =
                    std::sin(xFace.x) * std::cos(xFace.y);
                velocity.yFaces()[index] =
                    -std::cos(yFace.x) * std::sin(yFace.y);
            }
        }
    }
    const auto original = velocity;
    CellScalarField pressure(grid);
    auto settings = strictSettings();
    settings.absoluteResidualTolerance = 1.0e-10;
    const auto diagnostics = projectVelocity(grid, velocity, pressure, settings);
    check(diagnostics.converged && diagnostics.iterationCount == 0,
          "Taylor-Green: the discrete solenoidal field needs no correction");
    check(velocity == original,
          "Taylor-Green: projection preserves a solenoidal field bit-for-bit");
    check(diagnostics.divergenceMaximumAfterPerSecond < 2.0e-14,
          "Taylor-Green: staggered analytic velocity has roundoff divergence");
    checkNear(diagnostics.kineticEnergyAfterJoules,
              diagnostics.kineticEnergyBeforeJoules,
              0.0,
              "Taylor-Green: a zero correction preserves kinetic energy exactly");
}

double continuousManufacturedPressureError(const std::size_t resolution) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, 2, 2}, {}, {twoPi, 1.0, 1.0});
    auto settings = strictSettings();
    settings.absoluteResidualTolerance = 1.0e-12;
    MacVelocityField predictedVelocity(grid);
    const double scale = settings.timeStepSeconds
        / settings.densityKgPerCubicMeter;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const auto face = grid.xFaceCenterMeters(i, j, k);
                predictedVelocity.xFaces()[index] = scale * std::cos(face.x);
            }
        }
    }
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocity(
        grid, predictedVelocity, pressure, settings);
    check(diagnostics.converged,
          "continuous manufactured solution: pressure solve converges");
    double sumSquaredError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const double expected = std::sin(
                    grid.cellCenterMeters(i, j, k).x);
                const double error = pressure.values()[index] - expected;
                sumSquaredError += error * error;
            }
        }
    }
    return std::sqrt(sumSquaredError / static_cast<double>(grid.cellCount()));
}

void testManufacturedSecondOrderConvergence() {
    const double coarseError = continuousManufacturedPressureError(12);
    const double mediumError = continuousManufacturedPressureError(24);
    const double fineError = continuousManufacturedPressureError(48);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 3.9 && coarseRatio < 4.2,
          "convergence: first pressure refinement is second order");
    check(fineRatio > 3.9 && fineRatio < 4.2,
          "convergence: second pressure refinement is second order");
}

MacVelocityField deterministicPredictedVelocity(
    const PeriodicCartesianGrid& grid) {
    MacVelocityField velocity(grid);
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        const double sample = static_cast<double>(index + 1);
        velocity.xFaces()[index] = 0.3 + std::sin(0.37 * sample);
        velocity.yFaces()[index] = -0.2 + std::cos(0.23 * sample);
        velocity.zFaces()[index] = 0.1 + std::sin(0.19 * sample + 0.4);
    }
    return velocity;
}

void testProjectionConservationDeterminismAndRollback() {
    const auto grid = unitGrid({7, 6, 5});
    auto firstVelocity = deterministicPredictedVelocity(grid);
    auto secondVelocity = firstVelocity;
    const double xSum = componentSum(firstVelocity.xFaces());
    const double ySum = componentSum(firstVelocity.yFaces());
    const double zSum = componentSum(firstVelocity.zFaces());
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    const auto settings = strictSettings();
    const auto firstDiagnostics = projectVelocity(
        grid, firstVelocity, firstPressure, settings);
    const auto secondDiagnostics = projectVelocity(
        grid, secondVelocity, secondPressure, settings);
    check(firstDiagnostics.converged,
          "projection: deterministic mixed-mode solve converges");
    check(firstVelocity == secondVelocity
              && firstPressure == secondPressure
              && firstDiagnostics == secondDiagnostics,
          "projection: identical inputs replay bit-for-bit");
    checkNear(componentSum(firstVelocity.xFaces()), xSum, 5.0e-14,
              "projection: periodic correction preserves x momentum");
    checkNear(componentSum(firstVelocity.yFaces()), ySum, 5.0e-14,
              "projection: periodic correction preserves y momentum");
    checkNear(componentSum(firstVelocity.zFaces()), zSum, 5.0e-14,
              "projection: periodic correction preserves z momentum");
    check(firstDiagnostics.kineticEnergyAfterJoules
              <= firstDiagnostics.kineticEnergyBeforeJoules + 1.0e-13,
          "projection: the orthogonal correction cannot add kinetic energy");
    check(firstDiagnostics.divergenceL2AfterPerSecond < 2.0e-11,
          "projection: mixed-mode divergence satisfies the declared solve tolerance");

    auto failedVelocity = deterministicPredictedVelocity(grid);
    CellScalarField failedPressure(grid, 0.125);
    const auto originalVelocity = failedVelocity;
    const auto originalPressure = failedPressure;
    auto failingSettings = settings;
    failingSettings.absoluteResidualTolerance = 1.0e-30;
    failingSettings.relativeResidualTolerance = 0.0;
    failingSettings.maximumIterations = 1;
    const auto failedDiagnostics = projectVelocity(
        grid, failedVelocity, failedPressure, failingSettings);
    check(!failedDiagnostics.converged,
          "rollback: an intentionally truncated solve reports non-convergence");
    check(failedVelocity == originalVelocity && failedPressure == originalPressure,
          "rollback: a failed projection commits neither pressure nor velocity");
}

void testInvalidInputsAreRejectedTransactionally() {
    const auto grid = unitGrid({4, 4, 4});
    auto velocity = deterministicPredictedVelocity(grid);
    CellScalarField pressure(grid);
    const auto originalVelocity = velocity;
    auto invalidSettings = strictSettings();
    invalidSettings.timeStepSeconds = 0.0;
    bool rejected = false;
    try {
        static_cast<void>(projectVelocity(
            grid, velocity, pressure, invalidSettings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && velocity == originalVelocity,
          "validation: invalid settings are rejected before mutation");

    velocity.xFaces()[0] = std::numeric_limits<double>::quiet_NaN();
    rejected = false;
    try {
        static_cast<void>(projectVelocity(
            grid, velocity, pressure, strictSettings()));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected && std::isnan(velocity.xFaces()[0]),
          "validation: non-finite input is rejected without rewriting it");
}

} // namespace

int main() {
    testGridContractAndValidation();
    testDiscreteOperatorIdentity();
    testDiscretelyManufacturedPressureProjection();
    testTaylorGreenIsPreserved();
    testManufacturedSecondOrderConvergence();
    testProjectionConservationDeterminismAndRollback();
    testInvalidInputsAreRejectedTransactionally();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing fluid projection check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid projection checks passed");
    return 0;
}
