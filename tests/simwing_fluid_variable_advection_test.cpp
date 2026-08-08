#include "fluid/advection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::UniformMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacReconstruction;
using simwing::fsi::fluid::advectVelocityByMacFlow;
using simwing::fsi::fluid::advectVelocityByMacFlowSspRk2;
using simwing::fsi::fluid::advectVelocityByUniformFlow;

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
        result.xFaces()[index] = 0.4 + std::sin(0.37 * sample);
        result.yFaces()[index] = -0.2 + std::cos(0.23 * sample + 0.2);
        result.zFaces()[index] = 0.1 + std::sin(0.19 * sample - 0.3);
    }
    return result;
}

MacVelocityField uniformAdvector(
    const PeriodicCartesianGrid& grid,
    const double x,
    const double y,
    const double z) {
    MacVelocityField result(grid);
    std::ranges::fill(result.xFaces(), x);
    std::ranges::fill(result.yFaces(), y);
    std::ranges::fill(result.zFaces(), z);
    return result;
}

MacVelocityField shearAdvector(const PeriodicCartesianGrid& grid) {
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const double y = grid.xFaceCenterMeters(i, j, k).y;
                result.xFaces()[index] = 0.75 + 0.25 * std::sin(y);
            }
        }
    }
    return result;
}

MacVelocityField taylorGreen(const PeriodicCartesianGrid& grid) {
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                result.xFaces()[index] =
                    std::sin(xFace.x) * std::cos(xFace.y);
                result.yFaces()[index] =
                    -std::cos(yFace.x) * std::sin(yFace.y);
            }
        }
    }
    return result;
}

VariableMacAdvectionSettings settings() {
    VariableMacAdvectionSettings result;
    result.densityKgPerCubicMeter = 1.2;
    result.timeStepSeconds = 0.05;
    result.absoluteDivergenceTolerancePerSecond = 3.0e-12;
    result.relativeDivergenceTolerance = 1.0e-12;
    result.absoluteMomentumToleranceNewtonSeconds = 4.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 4.0e-12;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

void testExactUniformOracleSubset() {
    const PeriodicCartesianGrid grid(
        {7, 6, 5}, {}, {2.0, 3.0, 4.0});
    auto expected = deterministicVelocity(grid);
    auto actual = expected;
    const auto advector = uniformAdvector(grid, 0.2, -0.1, 0.05);
    auto variableSettings = settings();
    variableSettings.timeStepSeconds = 0.1;
    UniformMacAdvectionSettings uniformSettings;
    uniformSettings.densityKgPerCubicMeter =
        variableSettings.densityKgPerCubicMeter;
    uniformSettings.transportVelocityMetersPerSecond = {0.2, -0.1, 0.05};
    uniformSettings.timeStepSeconds = variableSettings.timeStepSeconds;
    uniformSettings.maximumTotalCourantNumber =
        variableSettings.maximumLocalOutgoingCourantNumber;
    uniformSettings.absoluteMomentumToleranceNewtonSeconds =
        variableSettings.absoluteMomentumToleranceNewtonSeconds;
    uniformSettings.relativeMomentumTolerance =
        variableSettings.relativeMomentumTolerance;
    uniformSettings.absoluteEnergyToleranceJoules =
        variableSettings.absoluteEnergyToleranceJoules;
    uniformSettings.relativeEnergyTolerance =
        variableSettings.relativeEnergyTolerance;
    uniformSettings.absoluteBoundToleranceMetersPerSecond =
        variableSettings.absoluteBoundToleranceMetersPerSecond;
    uniformSettings.relativeBoundTolerance =
        variableSettings.relativeBoundTolerance;
    const auto uniform = advectVelocityByUniformFlow(
        grid, expected, uniformSettings);
    const auto variable = advectVelocityByMacFlow(
        grid, actual, advector, variableSettings);
    check(uniform.accepted && variable.accepted
              && variable.uniformAdvector
              && variable.divergenceCompatible
              && actual == expected,
          "uniform subset: variable-flow transport delegates bit-exactly to its oracle");
    check(variable.maximumLocalOutgoingCourantNumber
              == uniform.totalAbsoluteCourantNumber
              && variable.momentumAfterNewtonSeconds
                  == uniform.momentumAfterNewtonSeconds
              && variable.kineticEnergyAfterJoules
                  == uniform.kineticEnergyAfterJoules
              && variable.maximumVelocityChangeMetersPerSecond
                  == uniform.maximumVelocityChangeMetersPerSecond,
          "uniform subset: delegated diagnostics retain the oracle ledger exactly");
}

void testVariableFlowConservationBoundsAndDeterminism() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {18, 15, 3}, {}, {twoPi, twoPi, 1.0});
    const auto advector = shearAdvector(grid);
    auto first = deterministicVelocity(grid);
    auto second = first;
    const double xSum = componentSum(first.xFaces());
    const double ySum = componentSum(first.yFaces());
    const double zSum = componentSum(first.zFaces());
    auto flowSettings = settings();
    flowSettings.timeStepSeconds = 0.08;
    const auto firstDiagnostics = advectVelocityByMacFlow(
        grid, first, advector, flowSettings);
    const auto secondDiagnostics = advectVelocityByMacFlow(
        grid, second, advector, flowSettings);
    check(firstDiagnostics.accepted
              && firstDiagnostics.divergenceCompatible
              && firstDiagnostics.stable
              && firstDiagnostics.bounded
              && !firstDiagnostics.uniformAdvector,
          "variable flow: a divergence-free shear transport is accepted");
    check(first == second && firstDiagnostics == secondDiagnostics,
          "variable flow: identical staggered flux updates replay bit-for-bit");
    checkNear(componentSum(first.xFaces()), xSum, 2.0e-13,
              "variable flow: shared periodic fluxes preserve x momentum sum");
    checkNear(componentSum(first.yFaces()), ySum, 2.0e-13,
              "variable flow: shared periodic fluxes preserve y momentum sum");
    checkNear(componentSum(first.zFaces()), zSum, 2.0e-13,
              "variable flow: shared periodic fluxes preserve z momentum sum");
    check(firstDiagnostics.maximumAdvectingDivergencePerSecond == 0.0
              && firstDiagnostics.maximumControlVolumeDivergencePerSecond
                  == 0.0
              && firstDiagnostics.maximumBoundViolationMetersPerSecond
                  == 0.0
              && firstDiagnostics.momentumResidualNormNewtonSeconds
                  < 4.0e-12
              && firstDiagnostics.kineticEnergyAfterJoules
                  < firstDiagnostics.kineticEnergyBeforeJoules,
          "variable flow: divergence, bounds, momentum, and energy ledgers close");
}

void testNonlinearSelfAdvectionAlias() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {24, 24, 3}, {}, {twoPi, twoPi, 1.0});
    auto first = taylorGreen(grid);
    auto second = first;
    const auto before = first;
    auto flowSettings = settings();
    flowSettings.timeStepSeconds = 0.02;
    const auto firstDiagnostics = advectVelocityByMacFlow(
        grid, first, first, flowSettings);
    const auto secondDiagnostics = advectVelocityByMacFlow(
        grid, second, second, flowSettings);
    check(firstDiagnostics.accepted && first != before,
          "self advection: the transported field may safely alias its advector");
    check(first == second && firstDiagnostics == secondDiagnostics,
          "self advection: aliased nonlinear updates replay bit-for-bit");
    check(firstDiagnostics.maximumAdvectingDivergencePerSecond < 3.0e-14
              && firstDiagnostics.maximumControlVolumeDivergencePerSecond
                  < 3.0e-14
              && firstDiagnostics.momentumResidualNormNewtonSeconds
                  < 4.0e-12
              && firstDiagnostics.numericalKineticEnergyLossJoules > 0.0,
          "self advection: Taylor-Green closes divergence, momentum, and energy ledgers");
}

double shearTransportError(const std::size_t resolution) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, resolution, 2}, {}, {twoPi, twoPi, 1.0});
    const auto advector = shearAdvector(grid);
    MacVelocityField transported(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto center = grid.xFaceCenterMeters(i, j, k);
                transported.xFaces()[grid.cellIndex(i, j, k)] =
                    std::sin(center.x);
            }
        }
    }
    constexpr double finalTime = 1.0;
    const double requestedTimeStep =
        0.4 * grid.cellSpacingMeters().x;
    const std::size_t steps = static_cast<std::size_t>(
        std::ceil(finalTime / requestedTimeStep));
    auto flowSettings = settings();
    flowSettings.timeStepSeconds =
        finalTime / static_cast<double>(steps);
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advectVelocityByMacFlow(
            grid, transported, advector, flowSettings);
        check(diagnostics.accepted,
              "convergence: every variable-shear transport step is accepted");
    }
    double squaredError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto center = grid.xFaceCenterMeters(i, j, k);
                const double shear = 0.75 + 0.25 * std::sin(center.y);
                const double expected = std::sin(
                    center.x - shear * finalTime);
                const double error = transported.xFaces()[
                    grid.cellIndex(i, j, k)] - expected;
                squaredError += error * error;
            }
        }
    }
    return std::sqrt(
        squaredError / static_cast<double>(grid.cellCount()));
}

void testFirstOrderVariableFlowConvergence() {
    const double coarseError = shearTransportError(16);
    const double mediumError = shearTransportError(32);
    const double fineError = shearTransportError(64);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 1.7 && coarseRatio < 2.2,
          "convergence: variable-shear donor transport approaches first order");
    check(fineRatio > 1.8 && fineRatio < 2.1,
          "convergence: refined variable-shear transport remains first order");
}

void testMusclSspRk2ExactCompositionAndBounds() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {18, 15, 3}, {}, {twoPi, twoPi, 1.0});
    const auto advector = shearAdvector(grid);
    const auto original = deterministicVelocity(grid);
    auto reconstructedSettings = settings();
    reconstructedSettings.reconstruction =
        VariableMacReconstruction::MonotonizedCentral;
    reconstructedSettings.timeStepSeconds = 0.02;
    auto eulerSettings = reconstructedSettings;
    eulerSettings.enforceEulerEnergyNonIncrease = false;

    auto twiceAdvanced = original;
    const auto expectedFirst = advectVelocityByMacFlow(
        grid, twiceAdvanced, advector, eulerSettings);
    const auto expectedSecond = advectVelocityByMacFlow(
        grid, twiceAdvanced, advector, eulerSettings);
    auto expected = original;
    const auto average = [](const std::span<const double> before,
                            const std::span<const double> twice,
                            const std::span<double> destination) {
        for (std::size_t index = 0; index < before.size(); ++index) {
            destination[index] = before[index]
                + 0.5 * (twice[index] - before[index]);
        }
    };
    average(original.xFaces(), twiceAdvanced.xFaces(), expected.xFaces());
    average(original.yFaces(), twiceAdvanced.yFaces(), expected.yFaces());
    average(original.zFaces(), twiceAdvanced.zFaces(), expected.zFaces());

    auto first = original;
    auto second = original;
    const auto firstDiagnostics = advectVelocityByMacFlowSspRk2(
        grid, first, advector, reconstructedSettings);
    const auto secondDiagnostics = advectVelocityByMacFlowSspRk2(
        grid, second, advector, reconstructedSettings);
    check(expectedFirst.accepted && expectedSecond.accepted
              && firstDiagnostics.accepted
              && firstDiagnostics.reconstruction
                  == VariableMacReconstruction::MonotonizedCentral
              && firstDiagnostics.firstEulerStage == expectedFirst
              && firstDiagnostics.secondEulerStage == expectedSecond
              && first == expected,
          "MUSCL SSPRK2: result and diagnostics equal two Euler stages exactly");
    check(first == second && firstDiagnostics == secondDiagnostics,
          "MUSCL SSPRK2: reconstructed transport replays bit-for-bit");
    check(firstDiagnostics.bounded
              && firstDiagnostics.maximumBoundViolationMetersPerSecond
                  == 0.0
              && firstDiagnostics.momentumResidualNormNewtonSeconds
                  < 4.0e-12
              && firstDiagnostics.kineticEnergyAfterJoules
                  <= firstDiagnostics.kineticEnergyBeforeJoules,
          "MUSCL SSPRK2: bounds, momentum, and energy ledgers close");

    const PeriodicCartesianGrid pulseGrid(
        {32, 2, 2}, {}, {1.0, 1.0, 1.0});
    auto pulse = MacVelocityField(pulseGrid);
    for (std::size_t index = 0; index < pulseGrid.cellCount(); ++index) {
        const std::size_t i = index % pulseGrid.cellCounts().x;
        pulse.xFaces()[index] = i >= 8 && i < 20 ? 1.0 : 0.0;
    }
    const auto pulseAdvector = uniformAdvector(pulseGrid, 1.0, 0.0, 0.0);
    auto pulseSettings = reconstructedSettings;
    pulseSettings.timeStepSeconds = 0.4
        * pulseGrid.cellSpacingMeters().x;
    const auto pulseDiagnostics = advectVelocityByMacFlowSspRk2(
        pulseGrid, pulse, pulseAdvector, pulseSettings);
    check(pulseDiagnostics.accepted && pulseDiagnostics.bounded
              && pulseDiagnostics.componentMinimumAfterMetersPerSecond.x
                  >= 0.0
              && pulseDiagnostics.componentMaximumAfterMetersPerSecond.x
                  <= 1.0,
          "MUSCL SSPRK2: a discontinuous pulse stays inside its exact old-time bounds");
}

double musclFullPeriodError(const std::size_t resolution) {
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
    const auto advector = uniformAdvector(grid, 1.0, 0.0, 0.0);
    const double requestedTimeStep =
        0.35 * grid.cellSpacingMeters().x;
    const std::size_t steps = static_cast<std::size_t>(
        std::ceil(twoPi / requestedTimeStep));
    auto transportSettings = settings();
    transportSettings.reconstruction =
        VariableMacReconstruction::MonotonizedCentral;
    transportSettings.timeStepSeconds =
        twoPi / static_cast<double>(steps);
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advectVelocityByMacFlowSspRk2(
            grid, velocity, advector, transportSettings);
        check(diagnostics.accepted,
              "MUSCL convergence: every full-period SSPRK2 step is accepted");
        if (!diagnostics.accepted) {
            std::fprintf(
                stderr,
                "MUSCL rejection: first=%d stable=%d bounded=%d bound=%.17g energyLoss=%.17g cfl=%.17g\n",
                diagnostics.firstEulerStage.accepted ? 1 : 0,
                diagnostics.firstEulerStage.stable ? 1 : 0,
                diagnostics.firstEulerStage.bounded ? 1 : 0,
                diagnostics.firstEulerStage.maximumBoundViolationMetersPerSecond,
                diagnostics.firstEulerStage.numericalKineticEnergyLossJoules,
                diagnostics.firstEulerStage.maximumLocalOutgoingCourantNumber);
            break;
        }
    }
    double absoluteError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        const double error = velocity.xFaces()[index]
            - expected.xFaces()[index];
        absoluteError += std::abs(error);
    }
    return absoluteError / static_cast<double>(grid.cellCount());
}

void testMusclObservedSecondOrderSpatialConvergence() {
    const double coarseError = musclFullPeriodError(32);
    const double mediumError = musclFullPeriodError(64);
    const double fineError = musclFullPeriodError(128);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    if (!(coarseRatio > 3.1 && coarseRatio < 4.8)
        || !(fineRatio > 3.1 && fineRatio < 4.8)) {
        std::fprintf(stderr,
                     "MUSCL ratios: %.17g %.17g errors %.17g %.17g %.17g\n",
                     coarseRatio, fineRatio,
                     coarseError, mediumError, fineError);
    }
    check(coarseRatio > 3.1 && coarseRatio < 4.8,
          "MUSCL convergence: first smooth-wave L1 refinement approaches second order");
    check(fineRatio > 3.1 && fineRatio < 4.8,
          "MUSCL convergence: refined smooth-wave L1 transport remains near second order");
}

void testTransactionalRejection() {
    const PeriodicCartesianGrid grid(
        {8, 7, 6}, {}, {2.0, 3.0, 4.0});
    const auto original = deterministicVelocity(grid);
    auto divergentAdvector = uniformAdvector(grid, 0.0, 0.0, 0.0);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                divergentAdvector.xFaces()[grid.cellIndex(i, j, k)] =
                    std::sin(2.0 * std::numbers::pi
                             * static_cast<double>(i)
                             / static_cast<double>(counts.x));
            }
        }
    }
    auto divergent = original;
    const auto rejectedDivergence = advectVelocityByMacFlow(
        grid, divergent, divergentAdvector, settings());
    check(!rejectedDivergence.divergenceCompatible
              && !rejectedDivergence.accepted
              && divergent == original,
          "rollback: a divergent advector is rejected without mutation");

    auto unstable = original;
    const auto shear = shearAdvector(grid);
    auto unstableSettings = settings();
    unstableSettings.timeStepSeconds = 10.0;
    const auto rejectedCfl = advectVelocityByMacFlow(
        grid, unstable, shear, unstableSettings);
    check(!rejectedCfl.stable && !rejectedCfl.accepted
              && unstable == original,
          "rollback: excessive local outgoing CFL is rejected without mutation");

    auto invalid = original;
    auto invalidSettings = settings();
    invalidSettings.maximumLocalOutgoingCourantNumber = 1.01;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByMacFlow(
            grid, invalid, shear, invalidSettings)); },
        "validation: an unsafe variable-flow CFL ceiling is rejected");
    check(invalid == original,
          "validation: invalid settings are transactional");

    invalid = original;
    invalidSettings = settings();
    invalidSettings.reconstruction =
        static_cast<VariableMacReconstruction>(255);
    expectRejected(
        [&] { static_cast<void>(advectVelocityByMacFlow(
            grid, invalid, shear, invalidSettings)); },
        "validation: an unknown reconstruction is rejected");
    check(invalid == original,
          "validation: rejected reconstruction is transactional");

    auto nonfiniteAdvector = shear;
    nonfiniteAdvector.yFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    auto finiteInput = original;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByMacFlow(
            grid, finiteInput, nonfiniteAdvector, settings())); },
        "validation: a non-finite advector is rejected");
    check(finiteInput == original,
          "validation: non-finite advector rejection is transactional");

    const PeriodicCartesianGrid foreignGrid(
        {9, 7, 6}, {}, {2.0, 3.0, 4.0});
    auto mismatched = original;
    expectRejected(
        [&] { static_cast<void>(advectVelocityByMacFlow(
            foreignGrid, mismatched, shear, settings())); },
        "validation: mismatched variable-flow fields are rejected");
    check(mismatched == original,
          "validation: grid mismatch leaves the transported field unchanged");
}

} // namespace

int main() {
    testExactUniformOracleSubset();
    testVariableFlowConservationBoundsAndDeterminism();
    testNonlinearSelfAdvectionAlias();
    testFirstOrderVariableFlowConvergence();
    testMusclSspRk2ExactCompositionAndBounds();
    testMusclObservedSecondOrderSpatialConvergence();
    testTransactionalRejection();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing variable advection check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing variable advection checks passed");
    return 0;
}
