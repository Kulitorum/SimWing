#include "fluid/projection.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

double dot(const CellScalarField& first, const CellScalarField& second) {
    const auto firstValues = first.values();
    const auto secondValues = second.values();
    double result = 0.0;
    for (std::size_t index = 0; index < firstValues.size(); ++index) {
        result += firstValues[index] * secondValues[index];
    }
    return result;
}

void subtractMean(CellScalarField& field) {
    const double fieldMean = mean(field);
    for (double& value : field.values()) {
        value -= fieldMean;
    }
}

void validateSettings(const ProjectionSettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument("projection density must be finite and positive");
    }
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "projection time step must be finite and positive");
    }
    if (!std::isfinite(settings.absoluteResidualTolerance)
        || settings.absoluteResidualTolerance < 0.0
        || !std::isfinite(settings.relativeResidualTolerance)
        || settings.relativeResidualTolerance < 0.0
        || (settings.absoluteResidualTolerance == 0.0
            && settings.relativeResidualTolerance == 0.0)) {
        throw std::invalid_argument(
            "projection residual tolerances must be finite, non-negative, and not both zero");
    }
}

} // namespace

ProjectionDiagnostics projectVelocity(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const ProjectionSettings& settings) {
    validateSettings(settings);
    if (!predictedVelocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)) {
        throw std::invalid_argument("projection fields do not match their grid");
    }
    if (!isFinite(predictedVelocityMetersPerSecond)
        || !isFinite(pressurePascals)) {
        throw std::invalid_argument("projection fields must contain finite values");
    }

    ProjectionDiagnostics diagnostics;
    CellScalarField divergence(grid);
    computeDivergence(grid, predictedVelocityMetersPerSecond, divergence);
    diagnostics.compatibilityDivergencePerSecond = mean(divergence);
    diagnostics.divergenceL2BeforePerSecond = l2Norm(divergence);
    diagnostics.divergenceMaximumBeforePerSecond =
        maximumAbsoluteValue(divergence);
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, predictedVelocityMetersPerSecond,
        settings.densityKgPerCubicMeter);

    CellScalarField rightHandSide(grid);
    const double rightHandSideScale =
        -settings.densityKgPerCubicMeter / settings.timeStepSeconds;
    const double compatibilityRightHandSide = rightHandSideScale
        * diagnostics.compatibilityDivergencePerSecond;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        rightHandSide.values()[index] =
            rightHandSideScale * divergence.values()[index]
            - compatibilityRightHandSide;
    }

    CellScalarField candidatePressure = pressurePascals;
    subtractMean(candidatePressure);
    CellScalarField operatorPressure(grid);
    applyNegativeLaplacian(grid, candidatePressure, operatorPressure);
    CellScalarField residual(grid);
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        residual.values()[index] = rightHandSide.values()[index]
            - operatorPressure.values()[index];
    }
    subtractMean(residual);
    diagnostics.initialResidualPascalsPerSquareMeter = l2Norm(residual);
    diagnostics.finalResidualPascalsPerSquareMeter =
        diagnostics.initialResidualPascalsPerSquareMeter;
    const double convergenceThreshold = std::max(
        settings.absoluteResidualTolerance,
        settings.relativeResidualTolerance
            * diagnostics.initialResidualPascalsPerSquareMeter);

    CellScalarField direction = residual;
    CellScalarField operatorDirection(grid);
    double residualSquared = dot(residual, residual);
    diagnostics.converged =
        diagnostics.finalResidualPascalsPerSquareMeter <= convergenceThreshold;

    while (!diagnostics.converged
           && diagnostics.iterationCount < settings.maximumIterations) {
        applyNegativeLaplacian(grid, direction, operatorDirection);
        const double denominator = dot(direction, operatorDirection);
        if (!std::isfinite(denominator) || !(denominator > 0.0)) {
            break;
        }
        const double alpha = residualSquared / denominator;
        for (std::size_t index = 0; index < grid.cellCount(); ++index) {
            candidatePressure.values()[index] +=
                alpha * direction.values()[index];
            residual.values()[index] -=
                alpha * operatorDirection.values()[index];
        }
        ++diagnostics.iterationCount;
        subtractMean(candidatePressure);
        subtractMean(residual);
        const double nextResidualSquared = dot(residual, residual);
        diagnostics.finalResidualPascalsPerSquareMeter = std::sqrt(
            nextResidualSquared / static_cast<double>(grid.cellCount()));
        diagnostics.converged =
            diagnostics.finalResidualPascalsPerSquareMeter <= convergenceThreshold;
        if (diagnostics.converged) {
            residualSquared = nextResidualSquared;
            break;
        }
        if (!std::isfinite(nextResidualSquared) || !(residualSquared > 0.0)) {
            break;
        }
        const double beta = nextResidualSquared / residualSquared;
        for (std::size_t index = 0; index < grid.cellCount(); ++index) {
            direction.values()[index] = residual.values()[index]
                + beta * direction.values()[index];
        }
        subtractMean(direction);
        residualSquared = nextResidualSquared;
    }

    // Recursive CG residuals can drift from b-Ap on difficult systems. The
    // transactional acceptance decision is based on the explicitly
    // recomputed residual, never only on the recurrence's estimate.
    if (diagnostics.converged) {
        applyNegativeLaplacian(grid, candidatePressure, operatorPressure);
        for (std::size_t index = 0; index < grid.cellCount(); ++index) {
            residual.values()[index] = rightHandSide.values()[index]
                - operatorPressure.values()[index];
        }
        subtractMean(residual);
        diagnostics.finalResidualPascalsPerSquareMeter = l2Norm(residual);
        diagnostics.converged = std::isfinite(
            diagnostics.finalResidualPascalsPerSquareMeter)
            && diagnostics.finalResidualPascalsPerSquareMeter
                <= convergenceThreshold;
    }

    if (!diagnostics.converged) {
        diagnostics.divergenceL2AfterPerSecond =
            diagnostics.divergenceL2BeforePerSecond;
        diagnostics.divergenceMaximumAfterPerSecond =
            diagnostics.divergenceMaximumBeforePerSecond;
        diagnostics.kineticEnergyAfterJoules =
            diagnostics.kineticEnergyBeforeJoules;
        diagnostics.pressureMeanPascals = mean(pressurePascals);
        return diagnostics;
    }

    MacVelocityField pressureGradient(grid);
    computePressureGradient(grid, candidatePressure, pressureGradient);
    MacVelocityField candidateVelocity = predictedVelocityMetersPerSecond;
    const double correctionScale =
        settings.timeStepSeconds / settings.densityKgPerCubicMeter;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        candidateVelocity.xFaces()[index] -=
            correctionScale * pressureGradient.xFaces()[index];
        candidateVelocity.yFaces()[index] -=
            correctionScale * pressureGradient.yFaces()[index];
        candidateVelocity.zFaces()[index] -=
            correctionScale * pressureGradient.zFaces()[index];
    }
    CellScalarField correctedDivergence(grid);
    computeDivergence(grid, candidateVelocity, correctedDivergence);
    diagnostics.divergenceL2AfterPerSecond = l2Norm(correctedDivergence);
    diagnostics.divergenceMaximumAfterPerSecond =
        maximumAbsoluteValue(correctedDivergence);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.pressureMeanPascals = mean(candidatePressure);

    predictedVelocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
