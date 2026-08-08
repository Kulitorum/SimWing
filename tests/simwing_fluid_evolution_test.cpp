#include "fluid/evolution.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PeriodicFlowAdvectionMode;
using simwing::fsi::fluid::PeriodicFlowDiffusionMode;
using simwing::fsi::fluid::PeriodicFlowFailureStage;
using simwing::fsi::fluid::PeriodicFlowSettings;
using simwing::fsi::fluid::PeriodicMacDiffusionSettings;
using simwing::fsi::fluid::ProjectionSettings;
using simwing::fsi::fluid::UniformMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacAdvectionSettings;
using simwing::fsi::fluid::advancePeriodicFlow;
using simwing::fsi::fluid::advectVelocityByMacFlow;
using simwing::fsi::fluid::advectVelocityByUniformFlow;
using simwing::fsi::fluid::diffuseVelocityExplicit;
using simwing::fsi::fluid::diffuseVelocitySspRk2;
using simwing::fsi::fluid::projectVelocity;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
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

PeriodicFlowSettings settings() {
    PeriodicFlowSettings result;
    result.densityKgPerCubicMeter = 1.2;
    result.transportVelocityMetersPerSecond = {0.3, -0.2, 0.0};
    result.kinematicViscositySquareMetersPerSecond = 0.02;
    result.timeStepSeconds = 0.05;
    result.projectionAbsoluteResidualTolerance = 1.0e-11;
    result.projectionRelativeResidualTolerance = 1.0e-13;
    result.projectionMaximumIterations = 1000;
    result.absoluteMomentumToleranceNewtonSeconds = 3.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 3.0e-12;
    result.relativeEnergyTolerance = 1.0e-12;
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

void testExactStageCompositionAndDeterminism() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {24, 24, 4}, {}, {twoPi, twoPi, twoPi});
    const auto flowSettings = settings();
    auto expectedVelocity = taylorGreen(grid);
    CellScalarField expectedPressure(grid);

    UniformMacAdvectionSettings advectionSettings;
    advectionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    advectionSettings.transportVelocityMetersPerSecond =
        flowSettings.transportVelocityMetersPerSecond;
    advectionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    advectionSettings.absoluteMomentumToleranceNewtonSeconds =
        flowSettings.absoluteMomentumToleranceNewtonSeconds;
    advectionSettings.relativeMomentumTolerance =
        flowSettings.relativeMomentumTolerance;
    advectionSettings.absoluteEnergyToleranceJoules =
        flowSettings.absoluteEnergyToleranceJoules;
    advectionSettings.relativeEnergyTolerance =
        flowSettings.relativeEnergyTolerance;
    const auto expectedAdvection = advectVelocityByUniformFlow(
        grid, expectedVelocity, advectionSettings);

    PeriodicMacDiffusionSettings diffusionSettings;
    diffusionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    diffusionSettings.kinematicViscositySquareMetersPerSecond =
        flowSettings.kinematicViscositySquareMetersPerSecond;
    diffusionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    diffusionSettings.absoluteMomentumToleranceNewtonSeconds =
        flowSettings.absoluteMomentumToleranceNewtonSeconds;
    diffusionSettings.relativeMomentumTolerance =
        flowSettings.relativeMomentumTolerance;
    diffusionSettings.absoluteEnergyToleranceJoules =
        flowSettings.absoluteEnergyToleranceJoules;
    diffusionSettings.relativeEnergyTolerance =
        flowSettings.relativeEnergyTolerance;
    const auto expectedDiffusion = diffuseVelocityExplicit(
        grid, expectedVelocity, diffusionSettings);

    ProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    projectionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    projectionSettings.absoluteResidualTolerance =
        flowSettings.projectionAbsoluteResidualTolerance;
    projectionSettings.relativeResidualTolerance =
        flowSettings.projectionRelativeResidualTolerance;
    projectionSettings.maximumIterations =
        flowSettings.projectionMaximumIterations;
    const auto expectedProjection = projectVelocity(
        grid, expectedVelocity, expectedPressure, projectionSettings);

    auto firstVelocity = taylorGreen(grid);
    auto secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    const auto first = advancePeriodicFlow(
        grid, firstVelocity, firstPressure, flowSettings);
    const auto second = advancePeriodicFlow(
        grid, secondVelocity, secondPressure, flowSettings);
    check(first.accepted
              && first.failureStage == PeriodicFlowFailureStage::None
              && firstVelocity == expectedVelocity
              && firstPressure == expectedPressure,
          "composition: accepted macro-step equals its three verified stages exactly");
    check(first.uniformAdvection == expectedAdvection
              && first.explicitDiffusion == expectedDiffusion
              && first.projection == expectedProjection,
          "composition: stage diagnostics retain their exact standalone contracts");
    check(first == second
              && firstVelocity == secondVelocity
              && firstPressure == secondPressure,
          "composition: identical macro-steps replay bit-for-bit");
    check(first.projection.iterationCount == 0
              && first.finalDivergenceL2PerSecond < 3.0e-14,
          "composition: transported and diffused solenoidal mode needs no correction");
    check(first.momentumResidualNormNewtonSeconds < 3.0e-12
              && first.kineticEnergyAfterJoules
                  < first.kineticEnergyBeforeJoules
              && first.totalEnergyLossJoules > 0.0,
          "composition: final momentum and kinetic-energy ledgers close");
    check(std::abs(
              first.totalEnergyLossJoules
              - first.advectionNumericalEnergyLossJoules
              - first.viscousEnergyLossJoules
              - first.projectionEnergyLossJoules) < 5.0e-12,
          "composition: stage energy losses sum to the macro-step loss");
}

void testExactNoOp() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {16, 16, 2}, {}, {twoPi, twoPi, twoPi});
    auto velocity = taylorGreen(grid);
    CellScalarField pressure(grid);
    const auto originalVelocity = velocity;
    const auto originalPressure = pressure;
    auto noOpSettings = settings();
    noOpSettings.transportVelocityMetersPerSecond = {};
    noOpSettings.kinematicViscositySquareMetersPerSecond = 0.0;
    const auto diagnostics = advancePeriodicFlow(
        grid, velocity, pressure, noOpSettings);
    check(diagnostics.accepted
              && diagnostics.projection.iterationCount == 0
              && velocity == originalVelocity
              && pressure == originalPressure
              && diagnostics.totalEnergyLossJoules == 0.0,
          "no-op: solenoidal zero-transport inviscid step is bit-identical");
}

void testExactNonlinearStageComposition() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {24, 24, 3}, {}, {twoPi, twoPi, 1.0});
    auto flowSettings = settings();
    flowSettings.advectionMode =
        PeriodicFlowAdvectionMode::SelfAdvectingMac;
    flowSettings.diffusionMode = PeriodicFlowDiffusionMode::SspRk2;
    flowSettings.timeStepSeconds = 0.02;
    flowSettings.advectionAbsoluteDivergenceTolerancePerSecond = 3.0e-12;

    auto expectedVelocity = taylorGreen(grid);
    CellScalarField expectedPressure(grid);
    VariableMacAdvectionSettings advectionSettings;
    advectionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    advectionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    advectionSettings.maximumLocalOutgoingCourantNumber =
        flowSettings.maximumTotalCourantNumber;
    advectionSettings.absoluteDivergenceTolerancePerSecond =
        flowSettings.advectionAbsoluteDivergenceTolerancePerSecond;
    advectionSettings.relativeDivergenceTolerance =
        flowSettings.advectionRelativeDivergenceTolerance;
    advectionSettings.absoluteMomentumToleranceNewtonSeconds =
        flowSettings.absoluteMomentumToleranceNewtonSeconds;
    advectionSettings.relativeMomentumTolerance =
        flowSettings.relativeMomentumTolerance;
    advectionSettings.absoluteEnergyToleranceJoules =
        flowSettings.absoluteEnergyToleranceJoules;
    advectionSettings.relativeEnergyTolerance =
        flowSettings.relativeEnergyTolerance;
    const auto expectedAdvection = advectVelocityByMacFlow(
        grid, expectedVelocity, expectedVelocity, advectionSettings);

    PeriodicMacDiffusionSettings diffusionSettings;
    diffusionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    diffusionSettings.kinematicViscositySquareMetersPerSecond =
        flowSettings.kinematicViscositySquareMetersPerSecond;
    diffusionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    diffusionSettings.maximumDiffusionNumber =
        flowSettings.maximumDiffusionNumber;
    diffusionSettings.absoluteMomentumToleranceNewtonSeconds =
        flowSettings.absoluteMomentumToleranceNewtonSeconds;
    diffusionSettings.relativeMomentumTolerance =
        flowSettings.relativeMomentumTolerance;
    diffusionSettings.absoluteEnergyToleranceJoules =
        flowSettings.absoluteEnergyToleranceJoules;
    diffusionSettings.relativeEnergyTolerance =
        flowSettings.relativeEnergyTolerance;
    const auto expectedDiffusion = diffuseVelocitySspRk2(
        grid, expectedVelocity, diffusionSettings);

    ProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    projectionSettings.timeStepSeconds = flowSettings.timeStepSeconds;
    projectionSettings.absoluteResidualTolerance =
        flowSettings.projectionAbsoluteResidualTolerance;
    projectionSettings.relativeResidualTolerance =
        flowSettings.projectionRelativeResidualTolerance;
    projectionSettings.maximumIterations =
        flowSettings.projectionMaximumIterations;
    const auto expectedProjection = projectVelocity(
        grid, expectedVelocity, expectedPressure, projectionSettings);

    auto actualVelocity = taylorGreen(grid);
    CellScalarField actualPressure(grid);
    const auto actual = advancePeriodicFlow(
        grid, actualVelocity, actualPressure, flowSettings);
    check(expectedAdvection.accepted && expectedDiffusion.accepted
              && expectedProjection.converged && actual.accepted
              && actual.advectionMode
                  == PeriodicFlowAdvectionMode::SelfAdvectingMac
              && actualVelocity == expectedVelocity
              && actualPressure == expectedPressure,
          "nonlinear composition: macro-step equals its standalone stages exactly");
    check(actual.variableAdvection == expectedAdvection
              && actual.sspRk2Diffusion == expectedDiffusion
              && actual.projection == expectedProjection,
          "nonlinear composition: variable-advection diagnostics remain exact");
    check(actual.variableAdvection.numericalKineticEnergyLossJoules > 0.0
              && actual.projection.iterationCount > 0
              && actual.finalDivergenceL2PerSecond < 1.0e-10,
          "nonlinear composition: self-advection is dissipative and projection restores continuity");
}

void testDeterministicNonlinearSequence() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {20, 20, 3}, {}, {twoPi, twoPi, 1.0});
    auto firstVelocity = taylorGreen(grid);
    auto secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    auto flowSettings = settings();
    flowSettings.advectionMode =
        PeriodicFlowAdvectionMode::SelfAdvectingMac;
    flowSettings.diffusionMode = PeriodicFlowDiffusionMode::SspRk2;
    flowSettings.timeStepSeconds = 0.01;
    flowSettings.kinematicViscositySquareMetersPerSecond = 0.01;
    for (std::size_t step = 0; step < 12; ++step) {
        const auto first = advancePeriodicFlow(
            grid, firstVelocity, firstPressure, flowSettings);
        const auto second = advancePeriodicFlow(
            grid, secondVelocity, secondPressure, flowSettings);
        check(first.accepted && second.accepted,
              "nonlinear sequence: every projected self-advection step is accepted");
        check(first == second
                  && firstVelocity == secondVelocity
                  && firstPressure == secondPressure,
              "nonlinear sequence: repeated macro-steps replay bit-for-bit");
        check(first.variableAdvection.divergenceCompatible
                  && first.finalDivergenceL2PerSecond < 1.0e-10,
              "nonlinear sequence: each committed field remains eligible for its next step");
    }
}

void testStageFailureRollback() {
    const PeriodicCartesianGrid grid(
        {8, 7, 6}, {}, {2.0, 3.0, 4.0});
    const auto originalVelocity = deterministicVelocity(grid);
    const CellScalarField originalPressure(grid, 0.125);

    auto advectionVelocity = originalVelocity;
    auto advectionPressure = originalPressure;
    auto advectionFailure = settings();
    advectionFailure.transportVelocityMetersPerSecond = {20.0, 20.0, 20.0};
    advectionFailure.timeStepSeconds = 1.0;
    const auto rejectedAdvection = advancePeriodicFlow(
        grid, advectionVelocity, advectionPressure, advectionFailure);
    check(!rejectedAdvection.accepted
              && rejectedAdvection.failureStage
                  == PeriodicFlowFailureStage::Advection
              && advectionVelocity == originalVelocity
              && advectionPressure == originalPressure,
          "rollback: unstable advection commits neither velocity nor pressure");

    auto diffusionVelocity = originalVelocity;
    auto diffusionPressure = originalPressure;
    auto diffusionFailure = settings();
    diffusionFailure.transportVelocityMetersPerSecond = {};
    diffusionFailure.kinematicViscositySquareMetersPerSecond = 20.0;
    diffusionFailure.timeStepSeconds = 1.0;
    const auto rejectedDiffusion = advancePeriodicFlow(
        grid, diffusionVelocity, diffusionPressure, diffusionFailure);
    check(!rejectedDiffusion.accepted
              && rejectedDiffusion.failureStage
                  == PeriodicFlowFailureStage::Diffusion
              && diffusionVelocity == originalVelocity
              && diffusionPressure == originalPressure,
          "rollback: unstable diffusion discards an accepted advection candidate");

    auto projectionVelocity = originalVelocity;
    auto projectionPressure = originalPressure;
    auto projectionFailure = settings();
    projectionFailure.transportVelocityMetersPerSecond = {};
    projectionFailure.kinematicViscositySquareMetersPerSecond = 0.0;
    projectionFailure.projectionAbsoluteResidualTolerance = 1.0e-30;
    projectionFailure.projectionRelativeResidualTolerance = 0.0;
    projectionFailure.projectionMaximumIterations = 1;
    const auto rejectedProjection = advancePeriodicFlow(
        grid, projectionVelocity, projectionPressure, projectionFailure);
    check(!rejectedProjection.accepted
              && rejectedProjection.failureStage
                  == PeriodicFlowFailureStage::Projection
              && projectionVelocity == originalVelocity
              && projectionPressure == originalPressure,
          "rollback: failed projection discards all earlier stage candidates");
}

void testInvalidInputsAreTransactional() {
    const PeriodicCartesianGrid grid(
        {8, 7, 6}, {}, {2.0, 3.0, 4.0});
    const auto originalVelocity = deterministicVelocity(grid);
    const CellScalarField originalPressure(grid);
    auto velocity = originalVelocity;
    auto pressure = originalPressure;
    auto invalid = settings();
    invalid.maximumDiffusionNumber = 0.51;
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlow(
            grid, velocity, pressure, invalid)); },
        "validation: unsafe composed settings are rejected");
    check(velocity == originalVelocity && pressure == originalPressure,
          "validation: rejected settings mutate neither field");

    velocity = originalVelocity;
    pressure = originalPressure;
    invalid = settings();
    invalid.advectionMode = static_cast<PeriodicFlowAdvectionMode>(255);
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlow(
            grid, velocity, pressure, invalid)); },
        "validation: an unknown composed advection mode is rejected");
    check(velocity == originalVelocity && pressure == originalPressure,
          "validation: rejected advection mode mutates neither field");

    velocity = originalVelocity;
    pressure = originalPressure;
    invalid = settings();
    invalid.diffusionMode = static_cast<PeriodicFlowDiffusionMode>(255);
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlow(
            grid, velocity, pressure, invalid)); },
        "validation: an unknown composed diffusion mode is rejected");
    check(velocity == originalVelocity && pressure == originalPressure,
          "validation: rejected diffusion mode mutates neither field");

    velocity = originalVelocity;
    pressure = originalPressure;
    velocity.xFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlow(
            grid, velocity, pressure, settings())); },
        "validation: non-finite composed input is rejected");
    check(std::isnan(velocity.xFaces().front())
              && pressure == originalPressure,
          "validation: non-finite rejection does not rewrite either field");

    const PeriodicCartesianGrid foreignGrid(
        {9, 7, 6}, {}, {2.0, 3.0, 4.0});
    velocity = originalVelocity;
    pressure = originalPressure;
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlow(
            foreignGrid, velocity, pressure, settings())); },
        "validation: mismatched composed fields are rejected");
    check(velocity == originalVelocity && pressure == originalPressure,
          "validation: grid mismatch is transactional");
}

} // namespace

int main() {
    testExactStageCompositionAndDeterminism();
    testExactNoOp();
    testExactNonlinearStageComposition();
    testDeterministicNonlinearSequence();
    testStageFailureRollback();
    testInvalidInputsAreTransactional();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid evolution check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid evolution checks passed");
    return 0;
}
