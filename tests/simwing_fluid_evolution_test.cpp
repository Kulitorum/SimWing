#include "fluid/evolution.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFacePressureJump;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PeriodicFlowAdvectionMode;
using simwing::fsi::fluid::PeriodicFlowDiffusionMode;
using simwing::fsi::fluid::PeriodicFlowFailureStage;
using simwing::fsi::fluid::PeriodicFlowSettings;
using simwing::fsi::fluid::PeriodicFlowStrangFailureStage;
using simwing::fsi::fluid::PeriodicFlowStrangSspRk2Settings;
using simwing::fsi::fluid::PeriodicFlowStrangSubcyclingFailureStage;
using simwing::fsi::fluid::PeriodicFlowStrangSubcyclingSettings;
using simwing::fsi::fluid::PeriodicMacDiffusionSettings;
using simwing::fsi::fluid::ProjectedMacAdvectionSspRk2Settings;
using simwing::fsi::fluid::ProjectionSettings;
using simwing::fsi::fluid::PorousConstitutiveEvaluation;
using simwing::fsi::fluid::PorousGridFaceCrossing;
using simwing::fsi::fluid::PorousIterationSettings;
using simwing::fsi::fluid::PorousProjectionSettings;
using simwing::fsi::fluid::SharpPressureJumpField;
using simwing::fsi::fluid::UniformMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacReconstruction;
using simwing::fsi::fluid::advancePeriodicFlow;
using simwing::fsi::fluid::advancePeriodicFlowStrangSspRk2;
using simwing::fsi::fluid::advancePeriodicFlowStrangSspRk2Subcycled;
using simwing::fsi::fluid::advectVelocityProjectedSspRk2;
using simwing::fsi::fluid::advectVelocityByMacFlow;
using simwing::fsi::fluid::advectVelocityByUniformFlow;
using simwing::fsi::fluid::diffuseVelocityExplicit;
using simwing::fsi::fluid::diffuseVelocitySspRk2;
using simwing::fsi::fluid::projectVelocity;
using simwing::fsi::fluid::projectVelocityWithPorousInterfaces;

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

PeriodicFlowStrangSspRk2Settings strangSettings() {
    PeriodicFlowStrangSspRk2Settings result;
    result.densityKgPerCubicMeter = 1.2;
    result.kinematicViscositySquareMetersPerSecond = 0.02;
    result.timeStepSeconds = 0.01;
    result.advectionReconstruction =
        VariableMacReconstruction::MonotonizedCentral;
    result.advectionAbsoluteDivergenceTolerancePerSecond = 2.0e-10;
    result.advectionRelativeDivergenceTolerance = 1.0e-12;
    result.projectionAbsoluteResidualTolerance = 1.0e-12;
    result.projectionRelativeResidualTolerance = 1.0e-12;
    result.projectionMaximumIterations = 1000;
    result.absoluteMomentumToleranceNewtonSeconds = 4.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 4.0e-12;
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

double streamFunction(const double x, const double y) {
    return std::sin(x) * std::sin(y)
        + 0.18 * std::sin(2.0 * x + y)
        - 0.11 * std::cos(x - 2.0 * y);
}

MacVelocityField vorticalVelocity(
    const PeriodicCartesianGrid& grid) {
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const double x = lower.x
                    + static_cast<double>(i) * spacing.x;
                const double nextX = lower.x
                    + static_cast<double>(nextI) * spacing.x;
                const double y = lower.y
                    + static_cast<double>(j) * spacing.y;
                const double nextY = lower.y
                    + static_cast<double>(nextJ) * spacing.y;
                const std::size_t index = grid.cellIndex(i, j, k);
                result.xFaces()[index] =
                    (streamFunction(x, nextY) - streamFunction(x, y))
                    / spacing.y;
                result.yFaces()[index] =
                    -(streamFunction(nextX, y) - streamFunction(x, y))
                    / spacing.x;
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

PeriodicMacDiffusionSettings halfDiffusionSettings(
    const PeriodicFlowStrangSspRk2Settings& source) {
    PeriodicMacDiffusionSettings result;
    result.densityKgPerCubicMeter = source.densityKgPerCubicMeter;
    result.kinematicViscositySquareMetersPerSecond =
        source.kinematicViscositySquareMetersPerSecond;
    result.timeStepSeconds = 0.5 * source.timeStepSeconds;
    result.maximumDiffusionNumber = source.maximumDiffusionNumber;
    result.absoluteMomentumToleranceNewtonSeconds =
        source.absoluteMomentumToleranceNewtonSeconds;
    result.relativeMomentumTolerance = source.relativeMomentumTolerance;
    result.absoluteEnergyToleranceJoules =
        source.absoluteEnergyToleranceJoules;
    result.relativeEnergyTolerance = source.relativeEnergyTolerance;
    return result;
}

ProjectedMacAdvectionSspRk2Settings projectedSettings(
    const PeriodicFlowStrangSspRk2Settings& source) {
    ProjectedMacAdvectionSspRk2Settings result;
    result.densityKgPerCubicMeter = source.densityKgPerCubicMeter;
    result.timeStepSeconds = source.timeStepSeconds;
    result.reconstruction = source.advectionReconstruction;
    result.maximumLocalOutgoingCourantNumber =
        source.maximumLocalOutgoingCourantNumber;
    result.absoluteDivergenceTolerancePerSecond =
        source.advectionAbsoluteDivergenceTolerancePerSecond;
    result.relativeDivergenceTolerance =
        source.advectionRelativeDivergenceTolerance;
    result.projectionAbsoluteResidualTolerance =
        source.projectionAbsoluteResidualTolerance;
    result.projectionRelativeResidualTolerance =
        source.projectionRelativeResidualTolerance;
    result.projectionMaximumIterations = source.projectionMaximumIterations;
    result.absoluteMomentumToleranceNewtonSeconds =
        source.absoluteMomentumToleranceNewtonSeconds;
    result.relativeMomentumTolerance = source.relativeMomentumTolerance;
    result.absoluteEnergyToleranceJoules =
        source.absoluteEnergyToleranceJoules;
    result.relativeEnergyTolerance = source.relativeEnergyTolerance;
    return result;
}

void testStrangExactCompositionAndDeterminism() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {14, 12, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid);
    const auto integrationSettings = strangSettings();
    const auto diffusionSettings = halfDiffusionSettings(
        integrationSettings);
    const auto advectionSettings = projectedSettings(
        integrationSettings);

    auto expectedVelocity = originalVelocity;
    auto expectedPressure = originalPressure;
    const auto expectedFirstDiffusion = diffuseVelocitySspRk2(
        grid, expectedVelocity, diffusionSettings);
    const auto expectedAdvection = advectVelocityProjectedSspRk2(
        grid, expectedVelocity, expectedPressure, advectionSettings);
    const auto expectedSecondDiffusion = diffuseVelocitySspRk2(
        grid, expectedVelocity, diffusionSettings);

    auto firstVelocity = originalVelocity;
    auto secondVelocity = originalVelocity;
    auto firstPressure = originalPressure;
    auto secondPressure = originalPressure;
    const auto first = advancePeriodicFlowStrangSspRk2(
        grid, firstVelocity, firstPressure, integrationSettings);
    const auto second = advancePeriodicFlowStrangSspRk2(
        grid, secondVelocity, secondPressure, integrationSettings);
    check(expectedFirstDiffusion.accepted
              && expectedAdvection.accepted
              && expectedSecondDiffusion.accepted
              && first.accepted
              && first.failureStage
                  == PeriodicFlowStrangFailureStage::None
              && firstVelocity == expectedVelocity
              && firstPressure == expectedPressure,
          "Strang composition: result equals the three standalone sub-integrators exactly");
    check(first.firstHalfDiffusion == expectedFirstDiffusion
              && first.projectedAdvection == expectedAdvection
              && first.secondHalfDiffusion == expectedSecondDiffusion,
          "Strang composition: every sub-integrator diagnostic remains exact");
    check(first == second && firstVelocity == secondVelocity
              && firstPressure == secondPressure,
          "Strang composition: identical split steps replay bit-for-bit");
    check(first.momentumResidualNormNewtonSeconds < 4.0e-12
              && first.kineticEnergyAfterJoules
                  < first.kineticEnergyBeforeJoules
              && first.finalDivergenceL2PerSecond < 2.0e-10,
          "Strang composition: aggregate momentum, energy, and continuity close");
    check(std::abs(
              first.totalEnergyLossJoules
              - first.firstHalfViscousEnergyLossJoules
              - first.transportProjectionEnergyLossJoules
              - first.secondHalfViscousEnergyLossJoules) < 8.0e-12,
          "Strang composition: sub-integrator energy losses sum to the full step");

    const PeriodicCartesianGrid uniformGrid(
        {8, 7, 3}, {}, {2.0, 3.0, 4.0});
    MacVelocityField uniform(uniformGrid);
    std::ranges::fill(uniform.xFaces(), 0.5);
    std::ranges::fill(uniform.yFaces(), -0.25);
    std::ranges::fill(uniform.zFaces(), 0.125);
    const auto uniformBefore = uniform;
    CellScalarField uniformPressure(uniformGrid);
    const auto noOp = advancePeriodicFlowStrangSspRk2(
        uniformGrid, uniform, uniformPressure, integrationSettings);
    check(noOp.accepted && uniform == uniformBefore
              && noOp.totalEnergyLossJoules == 0.0
              && noOp.maximumVelocityChangeMetersPerSecond == 0.0,
          "Strang composition: uniform flow is a bit-exact no-op");
}

void testStrangSubcyclingExactComposition() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 10, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid);
    PeriodicFlowStrangSubcyclingSettings subcyclingSettings;
    subcyclingSettings.flow = strangSettings();
    subcyclingSettings.flow.kinematicViscositySquareMetersPerSecond = 1.2;
    subcyclingSettings.flow.timeStepSeconds = 0.1;
    subcyclingSettings.maximumSubsteps = 8;

    auto expectedVelocity = originalVelocity;
    auto expectedPressure = originalPressure;
    auto manualSettings = subcyclingSettings.flow;
    manualSettings.timeStepSeconds = 0.05;
    const auto expectedFirst = advancePeriodicFlowStrangSspRk2(
        grid, expectedVelocity, expectedPressure, manualSettings);
    const auto expectedSecond = advancePeriodicFlowStrangSspRk2(
        grid, expectedVelocity, expectedPressure, manualSettings);

    auto firstVelocity = originalVelocity;
    auto secondVelocity = originalVelocity;
    auto firstPressure = originalPressure;
    auto secondPressure = originalPressure;
    const auto first = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, firstVelocity, firstPressure, subcyclingSettings);
    const auto second = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, secondVelocity, secondPressure, subcyclingSettings);
    check(expectedFirst.accepted && expectedSecond.accepted
              && first.accepted
              && first.failureStage
                  == PeriodicFlowStrangSubcyclingFailureStage::None
              && first.plannedSubstepCount == 2
              && first.completedSubstepCount == 2
              && first.stabilityRetryCount == 0
              && first.substepSeconds == 0.05
              && first.substeps.size() == 2
              && first.substeps[0] == expectedFirst
              && first.substeps[1] == expectedSecond
              && firstVelocity == expectedVelocity
              && firstPressure == expectedPressure,
          "Strang subcycling: viscosity sizing equals two manual split steps exactly");
    check(first == second && firstVelocity == secondVelocity
              && firstPressure == secondPressure,
          "Strang subcycling: identical outer intervals replay bit-for-bit");
    check(first.maximumObservedDiffusionNumber <= 0.5
              && first.momentumResidualNormNewtonSeconds < 4.0e-12
              && first.kineticEnergyAfterJoules
                  <= first.kineticEnergyBeforeJoules
              && first.finalDivergenceL2PerSecond < 2.0e-10,
          "Strang subcycling: aggregate stability and conservation ledgers close");
}

void testStrangSubcyclingStabilityRetryAndRollback() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid, 0.125);
    PeriodicFlowStrangSubcyclingSettings retrySettings;
    retrySettings.flow = strangSettings();
    retrySettings.flow.kinematicViscositySquareMetersPerSecond = 0.0;
    retrySettings.flow.timeStepSeconds = 1.0;
    retrySettings.maximumSubsteps = 16;

    auto velocity = originalVelocity;
    auto pressure = originalPressure;
    const auto retried = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, velocity, pressure, retrySettings);
    auto replayVelocity = originalVelocity;
    auto replayPressure = originalPressure;
    const auto replayed = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, replayVelocity, replayPressure, retrySettings);
    if (!(retried.accepted && retried.stabilityRetryCount > 0
          && retried.plannedSubstepCount > 1)) {
        std::fprintf(
            stderr,
            "subcycling retry: accepted=%d stage=%u retries=%zu "
            "planned=%zu completed=%zu maxCfl=%.17g failed=%u "
            "nested=%u finite=%d currentCfl=%.17g stable=%d "
            "bounded=%d div=%d momentum=%.17g energy=%d\n",
            retried.accepted ? 1 : 0,
            static_cast<unsigned int>(retried.failureStage),
            retried.stabilityRetryCount,
            retried.plannedSubstepCount,
            retried.completedSubstepCount,
            retried.maximumObservedOutgoingCourantNumber,
            static_cast<unsigned int>(retried.failedSubstep.failureStage),
            static_cast<unsigned int>(
                retried.failedSubstep.projectedAdvection.failureStage),
            retried.finite ? 1 : 0,
            retried.failedSubstep.projectedAdvection.firstAdvection
                .maximumLocalOutgoingCourantNumber,
            retried.failedSubstep.projectedAdvection.firstAdvection.stable
                ? 1 : 0,
            retried.failedSubstep.projectedAdvection.firstAdvection.bounded
                ? 1 : 0,
            retried.failedSubstep.projectedAdvection.firstAdvection
                .divergenceCompatible ? 1 : 0,
            retried.failedSubstep.projectedAdvection.firstAdvection
                .momentumResidualNormNewtonSeconds,
            retried.failedSubstep.projectedAdvection.firstAdvection
                .energyNonIncreasing ? 1 : 0);
    }
    check(retried.accepted && retried.stabilityRetryCount > 0
              && retried.plannedSubstepCount > 1
              && retried.completedSubstepCount
                  == retried.plannedSubstepCount
              && retried.maximumObservedOutgoingCourantNumber > 1.0
              && retried.failedSubstep.failureStage
                  == PeriodicFlowStrangFailureStage::ProjectedAdvection,
          "Strang subcycling: an unstable CFL attempt restarts at a safe subdivision");
    check(retried == replayed && velocity == replayVelocity
              && pressure == replayPressure,
          "Strang subcycling: the stability retry replays bit-for-bit");
    auto expectedVelocity = originalVelocity;
    auto expectedPressure = originalPressure;
    auto manualSettings = retrySettings.flow;
    bool manualAccepted = retried.accepted
        && retried.plannedSubstepCount != 0
        && retried.substeps.size() == retried.plannedSubstepCount;
    if (manualAccepted) {
        manualSettings.timeStepSeconds = retrySettings.flow.timeStepSeconds
            / static_cast<double>(retried.plannedSubstepCount);
        for (std::size_t index = 0;
             index < retried.plannedSubstepCount; ++index) {
            const auto manual = advancePeriodicFlowStrangSspRk2(
                grid, expectedVelocity, expectedPressure, manualSettings);
            manualAccepted = manualAccepted && manual.accepted
                && manual == retried.substeps[index];
        }
    }
    check(manualAccepted && velocity == expectedVelocity
              && pressure == expectedPressure,
          "Strang subcycling: retried result equals its final equal-step schedule exactly");

    auto limitedVelocity = originalVelocity;
    auto limitedPressure = originalPressure;
    auto limitedSettings = retrySettings;
    limitedSettings.maximumSubsteps = 1;
    const auto limited = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, limitedVelocity, limitedPressure, limitedSettings);
    check(!limited.accepted
              && limited.failureStage
                  == PeriodicFlowStrangSubcyclingFailureStage::SubstepLimit
              && limited.plannedSubstepCount > 1
              && limitedVelocity == originalVelocity
              && limitedPressure == originalPressure,
          "Strang subcycling: a bounded retry limit rolls back the whole interval");

    auto failedVelocity = originalVelocity;
    auto failedPressure = originalPressure;
    auto failedSettings = retrySettings;
    failedSettings.flow.timeStepSeconds = 0.01;
    failedSettings.flow.projectionAbsoluteResidualTolerance = 1.0e-30;
    failedSettings.flow.projectionRelativeResidualTolerance = 0.0;
    failedSettings.flow.projectionMaximumIterations = 1;
    const auto failed = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, failedVelocity, failedPressure, failedSettings);
    check(!failed.accepted && failed.stabilityRetryCount == 0
              && failed.failureStage
                  == PeriodicFlowStrangSubcyclingFailureStage::Substep
              && failed.failedSubstep.failureStage
                  == PeriodicFlowStrangFailureStage::ProjectedAdvection
              && failedVelocity == originalVelocity
              && failedPressure == originalPressure,
          "Strang subcycling: projection failure is fatal and transactional");
}

struct StrangState {
    MacVelocityField velocity;
    CellScalarField pressure;

    explicit StrangState(const PeriodicCartesianGrid& grid)
        : velocity(vorticalVelocity(grid)), pressure(grid) {}
};

StrangState integrateStrang(const PeriodicCartesianGrid& grid,
                            const std::size_t steps,
                            const double finalTime) {
    StrangState result(grid);
    auto integrationSettings = strangSettings();
    integrationSettings.kinematicViscositySquareMetersPerSecond = 0.03;
    integrationSettings.timeStepSeconds =
        finalTime / static_cast<double>(steps);
    integrationSettings.projectionAbsoluteResidualTolerance = 1.0e-13;
    integrationSettings.projectionRelativeResidualTolerance = 1.0e-13;
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advancePeriodicFlowStrangSspRk2(
            grid, result.velocity, result.pressure, integrationSettings);
        check(diagnostics.accepted,
              "Strang refinement: every second-order split step is accepted");
    }
    return result;
}

double velocityError(const MacVelocityField& actual,
                     const MacVelocityField& expected) {
    double squaredError = 0.0;
    const auto accumulate = [&squaredError](
                                const std::span<const double> first,
                                const std::span<const double> second) {
        for (std::size_t index = 0; index < first.size(); ++index) {
            const double error = first[index] - second[index];
            squaredError += error * error;
        }
    };
    accumulate(actual.xFaces(), expected.xFaces());
    accumulate(actual.yFaces(), expected.yFaces());
    accumulate(actual.zFaces(), expected.zFaces());
    return std::sqrt(
        squaredError
        / static_cast<double>(3 * actual.xFaces().size()));
}

MacVelocityField viscousTranslatingTaylorGreen(
    const PeriodicCartesianGrid& grid,
    const double timeSeconds,
    const double kinematicViscositySquareMetersPerSecond) {
    constexpr double backgroundX = 0.35;
    constexpr double backgroundY = -0.2;
    const double amplitude = std::exp(
        -2.0 * kinematicViscositySquareMetersPerSecond * timeSeconds);
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                const double xFacePhaseX =
                    xFace.x - backgroundX * timeSeconds;
                const double xFacePhaseY =
                    xFace.y - backgroundY * timeSeconds;
                const double yFacePhaseX =
                    yFace.x - backgroundX * timeSeconds;
                const double yFacePhaseY =
                    yFace.y - backgroundY * timeSeconds;
                result.xFaces()[index] = backgroundX
                    + amplitude
                        * std::sin(xFacePhaseX) * std::cos(xFacePhaseY);
                result.yFaces()[index] = backgroundY
                    - amplitude
                        * std::cos(yFacePhaseX) * std::sin(yFacePhaseY);
            }
        }
    }
    return result;
}

double viscousTranslatingTaylorGreenError(const std::size_t resolution) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, resolution, 2}, {}, {twoPi, twoPi, 1.0});
    constexpr double viscosity = 0.02;
    constexpr double finalTime = 0.08;
    auto velocity = viscousTranslatingTaylorGreen(grid, 0.0, viscosity);
    CellScalarField pressure(grid);
    PeriodicFlowStrangSubcyclingSettings integrationSettings;
    integrationSettings.flow = strangSettings();
    integrationSettings.flow.densityKgPerCubicMeter = 1.225;
    integrationSettings.flow.kinematicViscositySquareMetersPerSecond =
        viscosity;
    integrationSettings.flow.advectionReconstruction =
        VariableMacReconstruction::MonotonizedCentral;
    integrationSettings.flow.advectionAbsoluteDivergenceTolerancePerSecond =
        2.0e-10;
    integrationSettings.flow.projectionAbsoluteResidualTolerance = 1.0e-11;
    integrationSettings.flow.projectionRelativeResidualTolerance = 1.0e-12;
    integrationSettings.flow.projectionMaximumIterations = 2000;
    integrationSettings.flow.absoluteMomentumToleranceNewtonSeconds = 1.0e-10;
    integrationSettings.flow.absoluteEnergyToleranceJoules = 1.0e-10;
    integrationSettings.maximumSubsteps = 64;
    const double spacing = grid.cellSpacingMeters().x;
    const double requestedTimeStep = 0.12 * spacing * spacing;
    const std::size_t steps = static_cast<std::size_t>(
        std::ceil(finalTime / requestedTimeStep));
    integrationSettings.flow.timeStepSeconds =
        finalTime / static_cast<double>(steps);
    bool accepted = true;
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advancePeriodicFlowStrangSspRk2Subcycled(
            grid, velocity, pressure, integrationSettings);
        accepted = accepted && diagnostics.accepted
            && diagnostics.completedSubstepCount
                == diagnostics.plannedSubstepCount
            && diagnostics.finalDivergenceL2PerSecond < 2.0e-10;
    }
    check(accepted,
          "analytic viscous Taylor-Green: every outer interval is accepted");

    const auto expected = viscousTranslatingTaylorGreen(
        grid, finalTime, viscosity);
    double absoluteError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        absoluteError += std::abs(
            velocity.xFaces()[index] - expected.xFaces()[index]);
        absoluteError += std::abs(
            velocity.yFaces()[index] - expected.yFaces()[index]);
    }
    return absoluteError
        / static_cast<double>(2 * grid.cellCount());
}

void testSubcycledStrangAnalyticViscousTaylorGreenRefinement() {
    const double coarseError = viscousTranslatingTaylorGreenError(12);
    const double mediumError = viscousTranslatingTaylorGreenError(24);
    const double fineError = viscousTranslatingTaylorGreenError(48);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    if (!(coarseRatio > 3.0 && coarseRatio < 5.0
          && fineRatio > 3.0 && fineRatio < 5.0
          && fineError < mediumError && mediumError < coarseError)) {
        std::fprintf(
            stderr,
            "analytic viscous Taylor-Green refinement: errors %.17g %.17g "
            "%.17g, ratios %.17g %.17g\n",
            coarseError, mediumError, fineError, coarseRatio, fineRatio);
    }
    check(coarseRatio > 3.0 && coarseRatio < 5.0,
          "analytic viscous Taylor-Green: first full-flow ratio is "
          "second order");
    check(fineRatio > 3.0 && fineRatio < 5.0,
          "analytic viscous Taylor-Green: refined ratio remains second order");
    check(fineError < mediumError && mediumError < coarseError,
          "analytic viscous Taylor-Green: L1 error decreases monotonically");
}

void testStrangObservedSecondOrderTemporalRefinement() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    constexpr double finalTime = 0.08;
    const auto reference = integrateStrang(grid, 256, finalTime);
    const auto coarse = integrateStrang(grid, 8, finalTime);
    const auto medium = integrateStrang(grid, 16, finalTime);
    const auto fine = integrateStrang(grid, 32, finalTime);
    const double coarseError = velocityError(
        coarse.velocity, reference.velocity);
    const double mediumError = velocityError(
        medium.velocity, reference.velocity);
    const double fineError = velocityError(
        fine.velocity, reference.velocity);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 3.5 && coarseRatio < 4.6,
          "Strang refinement: first full-flow ratio is second order");
    check(fineRatio > 3.3 && fineRatio < 4.7,
          "Strang refinement: refined full-flow ratio remains second order");
}

void testStrangRollback() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid, 0.125);

    auto diffusionVelocity = originalVelocity;
    auto diffusionPressure = originalPressure;
    auto diffusionFailure = strangSettings();
    diffusionFailure.kinematicViscositySquareMetersPerSecond = 100.0;
    diffusionFailure.timeStepSeconds = 1.0;
    const auto rejectedDiffusion = advancePeriodicFlowStrangSspRk2(
        grid, diffusionVelocity, diffusionPressure, diffusionFailure);
    check(!rejectedDiffusion.accepted
              && rejectedDiffusion.failureStage
                  == PeriodicFlowStrangFailureStage::FirstHalfDiffusion
              && diffusionVelocity == originalVelocity
              && diffusionPressure == originalPressure,
          "Strang rollback: unstable first half diffusion commits neither field");

    auto advectionVelocity = originalVelocity;
    auto advectionPressure = originalPressure;
    auto advectionFailure = strangSettings();
    advectionFailure.kinematicViscositySquareMetersPerSecond = 0.0;
    advectionFailure.timeStepSeconds = 10.0;
    const auto rejectedAdvection = advancePeriodicFlowStrangSspRk2(
        grid, advectionVelocity, advectionPressure, advectionFailure);
    check(!rejectedAdvection.accepted
              && rejectedAdvection.failureStage
                  == PeriodicFlowStrangFailureStage::ProjectedAdvection
              && advectionVelocity == originalVelocity
              && advectionPressure == originalPressure,
          "Strang rollback: failed projected transport discards the first half step");

    auto projectionVelocity = originalVelocity;
    auto projectionPressure = originalPressure;
    auto projectionFailure = strangSettings();
    projectionFailure.kinematicViscositySquareMetersPerSecond = 0.0;
    projectionFailure.projectionAbsoluteResidualTolerance = 1.0e-30;
    projectionFailure.projectionRelativeResidualTolerance = 0.0;
    projectionFailure.projectionMaximumIterations = 1;
    const auto rejectedProjection = advancePeriodicFlowStrangSspRk2(
        grid, projectionVelocity, projectionPressure, projectionFailure);
    check(!rejectedProjection.accepted
              && rejectedProjection.failureStage
                  == PeriodicFlowStrangFailureStage::ProjectedAdvection
              && projectionVelocity == originalVelocity
              && projectionPressure == originalPressure,
          "Strang rollback: failed internal projection commits neither field");

    auto invalidVelocity = originalVelocity;
    auto invalidPressure = originalPressure;
    auto invalidSettings = strangSettings();
    invalidSettings.maximumDiffusionNumber = 0.51;
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlowStrangSspRk2(
            grid, invalidVelocity, invalidPressure, invalidSettings)); },
        "Strang validation: unsafe settings are rejected");
    check(invalidVelocity == originalVelocity
              && invalidPressure == originalPressure,
          "Strang validation: invalid settings mutate neither field");

    invalidVelocity = originalVelocity;
    invalidPressure = originalPressure;
    invalidSettings = strangSettings();
    invalidSettings.advectionReconstruction =
        static_cast<VariableMacReconstruction>(255);
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlowStrangSspRk2(
            grid, invalidVelocity, invalidPressure, invalidSettings)); },
        "Strang validation: unknown reconstruction is rejected");
    check(invalidVelocity == originalVelocity
              && invalidPressure == originalPressure,
          "Strang validation: rejected reconstruction mutates neither field");

    PeriodicFlowStrangSubcyclingSettings invalidSubcycling;
    invalidSubcycling.flow = strangSettings();
    invalidSubcycling.maximumSubsteps = 0;
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlowStrangSspRk2Subcycled(
            grid, invalidVelocity, invalidPressure,
            invalidSubcycling)); },
        "Strang subcycling validation: a zero retry bound is rejected");
    check(invalidVelocity == originalVelocity
              && invalidPressure == originalPressure,
          "Strang subcycling validation: invalid settings mutate neither field");
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

SharpPressureJumpField staticSlabJumps(
    const PeriodicCartesianGrid& grid) {
    std::vector<GridFacePressureJump> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                10, 1, 2, GridFaceAxis::X,
                2, j, k, 100.0, 0.5});
            faces.push_back({
                20, 2, 1, GridFaceAxis::X,
                6, j, k, -100.0, 0.5});
        }
    }
    return SharpPressureJumpField(grid, std::move(faces));
}

void checkStaticSlabState(const PeriodicCartesianGrid& grid,
                          const MacVelocityField& velocity,
                          const CellScalarField& pressure,
                          const char* velocityMessage,
                          const char* pressureMessage) {
    double maximumVelocityError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        maximumVelocityError = std::max({
            maximumVelocityError,
            std::abs(velocity.xFaces()[index] - 0.4),
            std::abs(velocity.yFaces()[index]),
            std::abs(velocity.zFaces()[index]),
        });
    }
    check(maximumVelocityError < 3.0e-13, velocityMessage);
    const auto counts = grid.cellCounts();
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i >= 2 && i < 6
                    ? 50.0 : -50.0;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(pressure.values()[grid.cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 3.0e-12, pressureMessage);
}

void testStaticJumpsAcrossCompleteFlowPaths() {
    const PeriodicCartesianGrid grid(
        {8, 3, 2}, {}, {4.0, 1.5, 1.0});
    const auto jumps = staticSlabJumps(grid);
    MacVelocityField initialVelocity(grid);
    std::ranges::fill(initialVelocity.xFaces(), 0.4);

    auto firstOrderVelocity = initialVelocity;
    CellScalarField firstOrderPressure(grid);
    const auto firstOrder = advancePeriodicFlow(
        grid, firstOrderVelocity, firstOrderPressure,
        jumps, settings());
    check(firstOrder.accepted
              && firstOrder.projection.pressureJumpFaceCount == 12,
          "sharp full flow: first-order projection retains every crossing");
    checkStaticSlabState(
        grid, firstOrderVelocity, firstOrderPressure,
        "sharp full flow: first-order step creates no spurious velocity",
        "sharp full flow: first-order step retains analytic pressure");

    auto strangVelocity = initialVelocity;
    CellScalarField strangPressure(grid);
    const auto strang = advancePeriodicFlowStrangSspRk2(
        grid, strangVelocity, strangPressure,
        jumps, strangSettings());
    check(strang.accepted
              && strang.projectedAdvection.firstProjection
                     .pressureJumpFaceCount == 12
              && strang.projectedAdvection.secondProjection
                     .pressureJumpFaceCount == 12,
          "sharp full flow: Strang transport retains crossings at both stages");
    checkStaticSlabState(
        grid, strangVelocity, strangPressure,
        "sharp full flow: Strang step creates no spurious velocity",
        "sharp full flow: Strang step retains analytic pressure");

    PeriodicFlowStrangSubcyclingSettings subcyclingSettings;
    subcyclingSettings.flow = strangSettings();
    subcyclingSettings.flow.timeStepSeconds = 0.02;
    subcyclingSettings.maximumSubsteps = 64;
    auto subcycledVelocity = initialVelocity;
    CellScalarField subcycledPressure(grid);
    const auto subcycled = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, subcycledVelocity, subcycledPressure,
        jumps, subcyclingSettings);
    check(subcycled.accepted
              && subcycled.completedSubstepCount
                  == subcycled.plannedSubstepCount
              && std::ranges::all_of(
                  subcycled.substeps,
                  [](const auto& substep) {
                      return substep.projectedAdvection.firstProjection
                                     .pressureJumpFaceCount == 12
                          && substep.projectedAdvection.secondProjection
                                     .pressureJumpFaceCount == 12;
                  }),
          "sharp full flow: every private substep retains the immutable topology");
    checkStaticSlabState(
        grid, subcycledVelocity, subcycledPressure,
        "sharp full flow: subcycling creates no spurious velocity",
        "sharp full flow: subcycling retains analytic pressure");

    const SharpPressureJumpField empty(grid);
    auto legacyVelocity = initialVelocity;
    auto emptyVelocity = initialVelocity;
    CellScalarField legacyPressure(grid);
    CellScalarField emptyPressure(grid);
    const auto legacy = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, legacyVelocity, legacyPressure, subcyclingSettings);
    const auto withEmpty = advancePeriodicFlowStrangSspRk2Subcycled(
        grid, emptyVelocity, emptyPressure, empty, subcyclingSettings);
    check(legacy == withEmpty
              && legacyVelocity == emptyVelocity
              && legacyPressure == emptyPressure,
          "sharp full flow: empty subcycled overload is bit-exact to legacy path");

    const PeriodicCartesianGrid foreignGrid(
        {7, 3, 2}, {}, {3.5, 1.5, 1.0});
    const SharpPressureJumpField foreign(foreignGrid);
    auto rejectedVelocity = initialVelocity;
    CellScalarField rejectedPressure(grid);
    expectRejected(
        [&] { static_cast<void>(advancePeriodicFlowStrangSspRk2Subcycled(
            grid, rejectedVelocity, rejectedPressure,
            foreign, subcyclingSettings)); },
        "sharp full flow: foreign static topology is rejected");
    check(rejectedVelocity == initialVelocity
              && rejectedPressure == CellScalarField(grid),
          "sharp full flow: rejected topology mutates neither field");
}

std::vector<PorousGridFaceCrossing> composedPorousPlane(
    const PeriodicCartesianGrid& grid) {
    std::vector<PorousGridFaceCrossing> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                100, 1, 2, GridFaceAxis::X, 1, j, k,
                0.35, 0.0, {10.0, 0.0}});
        }
    }
    return result;
}

SharpPressureJumpField composedDrivingPlane(
    const PeriodicCartesianGrid& grid) {
    std::vector<GridFacePressureJump> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                200, 2, 1, GridFaceAxis::X, 3, j, k,
                20.0, 0.5});
        }
    }
    return SharpPressureJumpField(grid, std::move(result));
}

PorousIterationSettings composedPorousIteration() {
    PorousIterationSettings result;
    result.constitutiveEvaluation =
        PorousConstitutiveEvaluation::Midpoint;
    result.absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-13;
    result.relativeNormalVelocityTolerance = 1.0e-13;
    result.absolutePressureJumpTolerancePascals = 1.0e-12;
    result.relativePressureJumpTolerance = 1.0e-13;
    result.relaxation = 0.5;
    result.maximumNonlinearIterations = 200;
    return result;
}

PeriodicFlowSettings composedPorousFlowSettings() {
    PeriodicFlowSettings result;
    result.densityKgPerCubicMeter = 1.0;
    result.transportVelocityMetersPerSecond = {};
    result.kinematicViscositySquareMetersPerSecond = 0.0;
    result.timeStepSeconds = 0.1;
    result.projectionAbsoluteResidualTolerance = 1.0e-12;
    result.projectionRelativeResidualTolerance = 1.0e-13;
    result.projectionMaximumIterations = 1000;
    result.absoluteMomentumToleranceNewtonSeconds = 2.0e-10;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 2.0e-10;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

void testPorousProjectionAcrossCompleteFlow() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    const auto porous = composedPorousPlane(grid);
    const auto driving = composedDrivingPlane(grid);
    const auto iteration = composedPorousIteration();
    const auto flowSettings = composedPorousFlowSettings();
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = advancePeriodicFlow(
        grid, velocity, pressure, porous, driving,
        iteration, flowSettings);
    MacVelocityField replayVelocity(grid);
    CellScalarField replayPressure(grid);
    const auto replay = advancePeriodicFlow(
        grid, replayVelocity, replayPressure, porous, driving,
        iteration, flowSettings);

    constexpr double expectedVelocity = 4.0 / 9.0;
    double maximumVelocityError = 0.0;
    for (std::size_t face = 0; face < grid.cellCount(); ++face) {
        maximumVelocityError = std::max({
            maximumVelocityError,
            std::abs(velocity.xFaces()[face] - expectedVelocity),
            std::abs(velocity.yFaces()[face]),
            std::abs(velocity.zFaces()[face]),
        });
    }
    check(diagnostics.version
              == simwing::fsi::fluid::periodicFlowStepVersion
              && diagnostics.accepted
              && diagnostics.failureStage
                  == PeriodicFlowFailureStage::None
              && diagnostics.porousProjection.accepted,
          "porous full flow: all stages and the interface-aware ledger accept");
    check(velocity == replayVelocity
              && pressure == replayPressure
              && diagnostics == replay,
          "porous full flow: identical composed steps replay bit-for-bit");
    check(maximumVelocityError < 2.0e-13,
          "porous full flow: midpoint pressure stage reaches the analytic endpoint");
    check(diagnostics.projection.pressureJumpFaceCount == 12
              && diagnostics.porousProjection.samples.size() == 6,
          "porous full flow: composed diagnostics retain all fixed and porous crossings");
    check(diagnostics.momentumResidualNormNewtonSeconds < 2.0e-11,
          "porous full flow: jump impulse closes aggregate momentum");
    check(std::abs(diagnostics.totalEnergyLossJoules) < 2.0e-11,
          "porous full flow: midpoint jump work closes aggregate kinetic energy");
    check(std::abs(
              diagnostics.pressureJumpImpulseOnFluidNewtonSeconds.x
              - diagnostics.momentumAfterNewtonSeconds.x) < 2.0e-11,
          "porous full flow: diagnosed jump impulse produces bulk momentum");
    check(std::abs(
              diagnostics.pressureJumpWorkToFluidJoules
              - diagnostics.kineticEnergyAfterJoules) < 2.0e-11,
          "porous full flow: diagnosed jump work produces bulk kinetic energy");

    MacVelocityField directVelocity(grid);
    CellScalarField directPressure(grid);
    PorousProjectionSettings directSettings;
    directSettings.iteration = iteration;
    directSettings.projection.densityKgPerCubicMeter =
        flowSettings.densityKgPerCubicMeter;
    directSettings.projection.timeStepSeconds =
        flowSettings.timeStepSeconds;
    directSettings.projection.absoluteResidualTolerance =
        flowSettings.projectionAbsoluteResidualTolerance;
    directSettings.projection.relativeResidualTolerance =
        flowSettings.projectionRelativeResidualTolerance;
    directSettings.projection.maximumIterations =
        flowSettings.projectionMaximumIterations;
    const auto direct = projectVelocityWithPorousInterfaces(
        grid, directVelocity, directPressure, porous, driving,
        directSettings);
    check(velocity == directVelocity
              && pressure == directPressure
              && diagnostics.porousProjection == direct,
          "porous full flow: zero transport and viscosity compose the exact standalone pressure step");

    MacVelocityField decayingVelocity(grid);
    std::ranges::fill(decayingVelocity.xFaces(), 1.0);
    CellScalarField decayingPressure(grid);
    const auto decay = advancePeriodicFlow(
        grid, decayingVelocity, decayingPressure,
        porous, iteration, flowSettings);
    check(decay.accepted,
          "porous full flow: an unforced dissipative sheet passes the interface-aware ledger");
    check(std::abs(decayingVelocity.xFaces().front() - 7.0 / 9.0)
              < 2.0e-13,
          "porous full flow: unforced midpoint drag reaches its analytic endpoint");
    check(std::abs(decay.pressureJumpWorkToFluidJoules
                   + decay.porousProjection
                         .totalPorousDissipationJoules) < 2.0e-11,
          "porous full flow: stationary porous pressure work equals negative dissipation");
    check(std::abs(decay.totalEnergyLossJoules) < 2.0e-11,
          "porous full flow: unforced midpoint decay closes its energy ledger");

    MacVelocityField legacyVelocity(grid);
    MacVelocityField emptyVelocity(grid);
    CellScalarField legacyPressure(grid);
    CellScalarField emptyPressure(grid);
    const std::vector<PorousGridFaceCrossing> empty;
    const auto legacy = advancePeriodicFlow(
        grid, legacyVelocity, legacyPressure, flowSettings);
    const auto withEmpty = advancePeriodicFlow(
        grid, emptyVelocity, emptyPressure, empty,
        iteration, flowSettings);
    check(legacy == withEmpty
              && legacyVelocity == emptyVelocity
              && legacyPressure == emptyPressure,
          "porous full flow: empty topology delegates bit-exactly to the original path");

    MacVelocityField failedVelocity(grid);
    CellScalarField failedPressure(grid, 3.0);
    const auto originalFailedVelocity = failedVelocity;
    const auto originalFailedPressure = failedPressure;
    auto truncated = iteration;
    truncated.maximumNonlinearIterations = 1;
    truncated.absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-30;
    truncated.relativeNormalVelocityTolerance = 0.0;
    truncated.absolutePressureJumpTolerancePascals = 1.0e-30;
    truncated.relativePressureJumpTolerance = 0.0;
    const auto failed = advancePeriodicFlow(
        grid, failedVelocity, failedPressure, porous, driving,
        truncated, flowSettings);
    check(!failed.accepted
              && failed.failureStage
                  == PeriodicFlowFailureStage::Projection,
          "porous full flow: exhausted nonlinear projection is reported at its stage");
    check(failedVelocity == originalFailedVelocity
              && failedPressure == originalFailedPressure,
          "porous full flow: nonlinear failure rolls back the complete composed step");
}

} // namespace

int main() {
    testExactStageCompositionAndDeterminism();
    testExactNoOp();
    testExactNonlinearStageComposition();
    testDeterministicNonlinearSequence();
    testStrangExactCompositionAndDeterminism();
    testStrangSubcyclingExactComposition();
    testStrangSubcyclingStabilityRetryAndRollback();
    testSubcycledStrangAnalyticViscousTaylorGreenRefinement();
    testStrangObservedSecondOrderTemporalRefinement();
    testStrangRollback();
    testStageFailureRollback();
    testInvalidInputsAreTransactional();
    testStaticJumpsAcrossCompleteFlowPaths();
    testPorousProjectionAcrossCompleteFlow();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid evolution check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid evolution checks passed");
    return 0;
}
