#include "fluid/evolution.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

bool finite(const Vector3& value) noexcept {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 subtract(const Vector3& first, const Vector3& second) noexcept {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

double length(const Vector3& value) noexcept {
    return std::hypot(value.x, value.y, value.z);
}

double combinedTolerance(const double absoluteTolerance,
                         const double relativeTolerance,
                         const double firstScale,
                         const double secondScale) noexcept {
    return absoluteTolerance
        + relativeTolerance * std::max(firstScale, secondScale);
}

void validateSettings(const PeriodicFlowSettings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.transportVelocityMetersPerSecond.x,
        settings.transportVelocityMetersPerSecond.y,
        settings.transportVelocityMetersPerSecond.z,
        settings.kinematicViscositySquareMetersPerSecond,
        settings.timeStepSeconds,
        settings.maximumTotalCourantNumber,
        settings.advectionAbsoluteDivergenceTolerancePerSecond,
        settings.advectionRelativeDivergenceTolerance,
        settings.maximumDiffusionNumber,
        settings.projectionAbsoluteResidualTolerance,
        settings.projectionRelativeResidualTolerance,
        settings.absoluteMomentumToleranceNewtonSeconds,
        settings.relativeMomentumTolerance,
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
    };
    if (!std::ranges::all_of(finiteValues, [](const double value) {
            return std::isfinite(value);
        })
        || (settings.advectionMode
                != PeriodicFlowAdvectionMode::PrescribedUniform
            && settings.advectionMode
                != PeriodicFlowAdvectionMode::SelfAdvectingMac)
        || (settings.diffusionMode
                != PeriodicFlowDiffusionMode::ForwardEuler
            && settings.diffusionMode
                != PeriodicFlowDiffusionMode::SspRk2)
        || settings.densityKgPerCubicMeter <= 0.0
        || settings.kinematicViscositySquareMetersPerSecond < 0.0
        || settings.timeStepSeconds <= 0.0
        || settings.maximumTotalCourantNumber <= 0.0
        || settings.maximumTotalCourantNumber > 1.0
        || settings.advectionAbsoluteDivergenceTolerancePerSecond < 0.0
        || settings.advectionRelativeDivergenceTolerance < 0.0
        || settings.maximumDiffusionNumber <= 0.0
        || settings.maximumDiffusionNumber > 0.5
        || settings.projectionAbsoluteResidualTolerance < 0.0
        || settings.projectionRelativeResidualTolerance < 0.0
        || settings.projectionMaximumIterations == 0
        || settings.absoluteMomentumToleranceNewtonSeconds < 0.0
        || settings.relativeMomentumTolerance < 0.0
        || settings.absoluteEnergyToleranceJoules < 0.0
        || settings.relativeEnergyTolerance < 0.0) {
        throw std::invalid_argument(
            "periodic flow settings are invalid");
    }
}

Vector3 momentumNewtonSeconds(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocity,
    const double densityKgPerCubicMeter) {
    Vector3 result;
    for (const double value : velocity.xFaces()) {
        result.x += value;
    }
    for (const double value : velocity.yFaces()) {
        result.y += value;
    }
    for (const double value : velocity.zFaces()) {
        result.z += value;
    }
    const double sampleMass = densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters();
    result.x *= sampleMass;
    result.y *= sampleMass;
    result.z *= sampleMass;
    return result;
}

} // namespace

PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)) {
        throw std::invalid_argument(
            "periodic flow fields do not match their grid");
    }
    if (!isFinite(velocityMetersPerSecond) || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "periodic flow fields must be finite");
    }

    PeriodicFlowDiagnostics diagnostics;
    diagnostics.advectionMode = settings.advectionMode;
    diagnostics.diffusionMode = settings.diffusionMode;
    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    CellScalarField divergence(grid);
    computeDivergence(grid, velocityMetersPerSecond, divergence);
    diagnostics.initialDivergenceL2PerSecond = l2Norm(divergence);

    auto candidateVelocity = velocityMetersPerSecond;
    auto candidatePressure = pressurePascals;
    switch (settings.advectionMode) {
    case PeriodicFlowAdvectionMode::PrescribedUniform: {
        UniformMacAdvectionSettings advectionSettings;
        advectionSettings.densityKgPerCubicMeter =
            settings.densityKgPerCubicMeter;
        advectionSettings.transportVelocityMetersPerSecond =
            settings.transportVelocityMetersPerSecond;
        advectionSettings.timeStepSeconds = settings.timeStepSeconds;
        advectionSettings.maximumTotalCourantNumber =
            settings.maximumTotalCourantNumber;
        advectionSettings.absoluteMomentumToleranceNewtonSeconds =
            settings.absoluteMomentumToleranceNewtonSeconds;
        advectionSettings.relativeMomentumTolerance =
            settings.relativeMomentumTolerance;
        advectionSettings.absoluteEnergyToleranceJoules =
            settings.absoluteEnergyToleranceJoules;
        advectionSettings.relativeEnergyTolerance =
            settings.relativeEnergyTolerance;
        diagnostics.uniformAdvection = advectVelocityByUniformFlow(
            grid, candidateVelocity, advectionSettings);
        if (!diagnostics.uniformAdvection.accepted) {
            diagnostics.failureStage =
                PeriodicFlowFailureStage::Advection;
            diagnostics.finite = diagnostics.uniformAdvection.finite;
            return diagnostics;
        }
        diagnostics.advectionNumericalEnergyLossJoules =
            diagnostics.uniformAdvection.numericalKineticEnergyLossJoules;
        break;
    }
    case PeriodicFlowAdvectionMode::SelfAdvectingMac: {
        VariableMacAdvectionSettings advectionSettings;
        advectionSettings.densityKgPerCubicMeter =
            settings.densityKgPerCubicMeter;
        advectionSettings.timeStepSeconds = settings.timeStepSeconds;
        advectionSettings.maximumLocalOutgoingCourantNumber =
            settings.maximumTotalCourantNumber;
        advectionSettings.absoluteDivergenceTolerancePerSecond =
            settings.advectionAbsoluteDivergenceTolerancePerSecond;
        advectionSettings.relativeDivergenceTolerance =
            settings.advectionRelativeDivergenceTolerance;
        advectionSettings.absoluteMomentumToleranceNewtonSeconds =
            settings.absoluteMomentumToleranceNewtonSeconds;
        advectionSettings.relativeMomentumTolerance =
            settings.relativeMomentumTolerance;
        advectionSettings.absoluteEnergyToleranceJoules =
            settings.absoluteEnergyToleranceJoules;
        advectionSettings.relativeEnergyTolerance =
            settings.relativeEnergyTolerance;
        diagnostics.variableAdvection = advectVelocityByMacFlow(
            grid, candidateVelocity, candidateVelocity,
            advectionSettings);
        if (!diagnostics.variableAdvection.accepted) {
            diagnostics.failureStage =
                PeriodicFlowFailureStage::Advection;
            diagnostics.finite = diagnostics.variableAdvection.finite;
            return diagnostics;
        }
        diagnostics.advectionNumericalEnergyLossJoules =
            diagnostics.variableAdvection.numericalKineticEnergyLossJoules;
        break;
    }
    default:
        throw std::invalid_argument(
            "periodic flow advection mode is invalid");
    }

    PeriodicMacDiffusionSettings diffusionSettings;
    diffusionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    diffusionSettings.kinematicViscositySquareMetersPerSecond =
        settings.kinematicViscositySquareMetersPerSecond;
    diffusionSettings.timeStepSeconds = settings.timeStepSeconds;
    diffusionSettings.maximumDiffusionNumber =
        settings.maximumDiffusionNumber;
    diffusionSettings.absoluteMomentumToleranceNewtonSeconds =
        settings.absoluteMomentumToleranceNewtonSeconds;
    diffusionSettings.relativeMomentumTolerance =
        settings.relativeMomentumTolerance;
    diffusionSettings.absoluteEnergyToleranceJoules =
        settings.absoluteEnergyToleranceJoules;
    diffusionSettings.relativeEnergyTolerance =
        settings.relativeEnergyTolerance;
    bool diffusionAccepted = false;
    bool diffusionFinite = true;
    switch (settings.diffusionMode) {
    case PeriodicFlowDiffusionMode::ForwardEuler:
        diagnostics.explicitDiffusion = diffuseVelocityExplicit(
            grid, candidateVelocity, diffusionSettings);
        diffusionAccepted = diagnostics.explicitDiffusion.accepted;
        diffusionFinite = diagnostics.explicitDiffusion.finite;
        diagnostics.viscousEnergyLossJoules =
            diagnostics.explicitDiffusion.dissipatedKineticEnergyJoules;
        break;
    case PeriodicFlowDiffusionMode::SspRk2:
        diagnostics.sspRk2Diffusion = diffuseVelocitySspRk2(
            grid, candidateVelocity, diffusionSettings);
        diffusionAccepted = diagnostics.sspRk2Diffusion.accepted;
        diffusionFinite = diagnostics.sspRk2Diffusion.finite;
        diagnostics.viscousEnergyLossJoules =
            diagnostics.sspRk2Diffusion.dissipatedKineticEnergyJoules;
        break;
    default:
        throw std::invalid_argument(
            "periodic flow diffusion mode is invalid");
    }
    if (!diffusionAccepted) {
        diagnostics.failureStage =
            PeriodicFlowFailureStage::Diffusion;
        const bool advectionFinite = settings.advectionMode
                == PeriodicFlowAdvectionMode::PrescribedUniform
            ? diagnostics.uniformAdvection.finite
            : diagnostics.variableAdvection.finite;
        diagnostics.finite = advectionFinite
            && diffusionFinite;
        return diagnostics;
    }

    ProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    projectionSettings.timeStepSeconds = settings.timeStepSeconds;
    projectionSettings.absoluteResidualTolerance =
        settings.projectionAbsoluteResidualTolerance;
    projectionSettings.relativeResidualTolerance =
        settings.projectionRelativeResidualTolerance;
    projectionSettings.maximumIterations =
        settings.projectionMaximumIterations;
    diagnostics.projection = projectVelocity(
        grid, candidateVelocity, candidatePressure, projectionSettings);
    if (!diagnostics.projection.converged) {
        diagnostics.failureStage =
            PeriodicFlowFailureStage::Projection;
        const bool advectionFinite = settings.advectionMode
                == PeriodicFlowAdvectionMode::PrescribedUniform
            ? diagnostics.uniformAdvection.finite
            : diagnostics.variableAdvection.finite;
        diagnostics.finite = advectionFinite
            && diffusionFinite;
        return diagnostics;
    }

    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        diagnostics.momentumAfterNewtonSeconds,
        diagnostics.momentumBeforeNewtonSeconds);
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.projectionEnergyLossJoules =
        diagnostics.projection.kineticEnergyBeforeJoules
        - diagnostics.projection.kineticEnergyAfterJoules;
    diagnostics.totalEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.finalDivergenceL2PerSecond =
        diagnostics.projection.divergenceL2AfterPerSecond;
    const bool advectionFinite = settings.advectionMode
            == PeriodicFlowAdvectionMode::PrescribedUniform
        ? diagnostics.uniformAdvection.finite
        : diagnostics.variableAdvection.finite;
    diagnostics.finite = advectionFinite
        && diffusionFinite
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(
            diagnostics.advectionNumericalEnergyLossJoules)
        && std::isfinite(diagnostics.viscousEnergyLossJoules)
        && std::isfinite(diagnostics.projectionEnergyLossJoules)
        && std::isfinite(diagnostics.totalEnergyLossJoules)
        && std::isfinite(diagnostics.initialDivergenceL2PerSecond)
        && std::isfinite(diagnostics.finalDivergenceL2PerSecond)
        && isFinite(candidateVelocity) && isFinite(candidatePressure);
    const double momentumTolerance = combinedTolerance(
        settings.absoluteMomentumToleranceNewtonSeconds,
        settings.relativeMomentumTolerance,
        length(diagnostics.momentumBeforeNewtonSeconds),
        length(diagnostics.momentumAfterNewtonSeconds));
    const double energyTolerance = combinedTolerance(
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
        std::abs(diagnostics.kineticEnergyBeforeJoules),
        std::abs(diagnostics.kineticEnergyAfterJoules));
    diagnostics.accepted = diagnostics.finite
        && diagnostics.momentumResidualNormNewtonSeconds
            <= momentumTolerance
        && diagnostics.kineticEnergyAfterJoules
            <= diagnostics.kineticEnergyBeforeJoules + energyTolerance;
    if (!diagnostics.accepted) {
        diagnostics.failureStage =
            PeriodicFlowFailureStage::Conservation;
        return diagnostics;
    }

    velocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
