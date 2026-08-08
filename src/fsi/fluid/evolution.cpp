#include "fluid/evolution.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
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

Vector3 add(const Vector3& first, const Vector3& second) noexcept {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
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

void validateSettings(
    const PeriodicFlowStrangSspRk2Settings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.kinematicViscositySquareMetersPerSecond,
        settings.timeStepSeconds,
        settings.maximumLocalOutgoingCourantNumber,
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
        || (settings.advectionReconstruction
                != VariableMacReconstruction::DonorCell
            && settings.advectionReconstruction
                != VariableMacReconstruction::MonotonizedCentral)
        || settings.densityKgPerCubicMeter <= 0.0
        || settings.kinematicViscositySquareMetersPerSecond < 0.0
        || settings.timeStepSeconds <= 0.0
        || settings.maximumLocalOutgoingCourantNumber <= 0.0
        || settings.maximumLocalOutgoingCourantNumber > 1.0
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
            "periodic Strang-SSPRK2 flow settings are invalid");
    }
}

void validateSettings(
    const PeriodicFlowStrangSubcyclingSettings& settings) {
    validateSettings(settings.flow);
    if (settings.maximumSubsteps == 0
        || settings.maximumSubsteps
            > periodicFlowStrangMaximumSubsteps) {
        throw std::invalid_argument(
            "periodic Strang subcycling settings are invalid");
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

double maximumDifference(const std::span<const double> first,
                         const std::span<const double> second) noexcept {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result = std::max(
            result, std::abs(first[index] - second[index]));
    }
    return result;
}

double maximumOutgoingCourantNumber(
    const PeriodicFlowStrangSspRk2Diagnostics& diagnostics) noexcept {
    return std::max(
        diagnostics.projectedAdvection
            .firstAdvection.maximumLocalOutgoingCourantNumber,
        diagnostics.projectedAdvection
            .secondAdvection.maximumLocalOutgoingCourantNumber);
}

double maximumDiffusionNumber(
    const PeriodicFlowStrangSspRk2Diagnostics& diagnostics) noexcept {
    return std::max({
        diagnostics.firstHalfDiffusion
            .firstEulerStage.totalDiffusionNumber,
        diagnostics.firstHalfDiffusion
            .secondEulerStage.totalDiffusionNumber,
        diagnostics.secondHalfDiffusion
            .firstEulerStage.totalDiffusionNumber,
        diagnostics.secondHalfDiffusion
            .secondEulerStage.totalDiffusionNumber,
    });
}

double rejectedStabilityRatio(
    const PeriodicFlowStrangSspRk2Diagnostics& diagnostics) noexcept {
    if (diagnostics.failureStage
        == PeriodicFlowStrangFailureStage::FirstHalfDiffusion) {
        const auto& diffusion = diagnostics.firstHalfDiffusion;
        if (!diffusion.firstEulerStage.stable
            || (!diffusion.accepted
                && !diffusion.secondEulerStage.stable
                && diffusion.secondEulerStage.timeStepSeconds != 0.0)) {
            return maximumDiffusionNumber(diagnostics)
                / diffusion.firstEulerStage
                    .maximumAcceptedDiffusionNumber;
        }
    }
    if (diagnostics.failureStage
        == PeriodicFlowStrangFailureStage::SecondHalfDiffusion) {
        const auto& diffusion = diagnostics.secondHalfDiffusion;
        if (!diffusion.firstEulerStage.stable
            || (!diffusion.accepted
                && !diffusion.secondEulerStage.stable
                && diffusion.secondEulerStage.timeStepSeconds != 0.0)) {
            return maximumDiffusionNumber(diagnostics)
                / diffusion.firstEulerStage
                    .maximumAcceptedDiffusionNumber;
        }
    }
    if (diagnostics.failureStage
        == PeriodicFlowStrangFailureStage::ProjectedAdvection) {
        const auto& transport = diagnostics.projectedAdvection;
        if (transport.failureStage
                == ProjectedMacAdvectionFailureStage::FirstAdvection
            && !transport.firstAdvection.stable) {
            return transport.firstAdvection
                .maximumLocalOutgoingCourantNumber
                / transport.firstAdvection
                    .maximumAcceptedLocalOutgoingCourantNumber;
        }
        if (transport.failureStage
                == ProjectedMacAdvectionFailureStage::SecondAdvection
            && !transport.secondAdvection.stable) {
            return transport.secondAdvection
                .maximumLocalOutgoingCourantNumber
                / transport.secondAdvection
                    .maximumAcceptedLocalOutgoingCourantNumber;
        }
        const VariableMacAdvectionDiagnostics* rejectedAdvection = nullptr;
        if (transport.failureStage
            == ProjectedMacAdvectionFailureStage::FirstAdvection) {
            rejectedAdvection = &transport.firstAdvection;
        } else if (transport.failureStage
                   == ProjectedMacAdvectionFailureStage::SecondAdvection) {
            rejectedAdvection = &transport.secondAdvection;
        }
        if (rejectedAdvection != nullptr
            && rejectedAdvection->finite
            && rejectedAdvection->stable
            && rejectedAdvection->divergenceCompatible
            && !rejectedAdvection->bounded) {
            // Limited reconstruction has a stricter multidimensional
            // maximum-principle CFL than the donor outgoing-flux ceiling.
            // Halving the step is the bounded deterministic fallback.
            return 2.0;
        }
    }
    return 0.0;
}

std::size_t boundedSubstepCount(
    const double requestedCount,
    const std::size_t maximumSubsteps) noexcept {
    if (!std::isfinite(requestedCount)
        || requestedCount
            > static_cast<double>(maximumSubsteps)) {
        return maximumSubsteps + 1;
    }
    return std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(requestedCount)));
}

struct PorousPeriodicProjectionInput {
    const std::vector<PorousGridFaceCrossing>* crossings = nullptr;
    const SharpPressureJumpField* prescribedPressureJumps = nullptr;
    const PorousIterationSettings* iteration = nullptr;
};

PeriodicFlowDiagnostics advancePeriodicFlowImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PorousPeriodicProjectionInput* porousProjection,
    const PeriodicFlowSettings& settings);

PeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2Impl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PeriodicFlowStrangSspRk2Settings& settings);

PorousPeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithPorousInterfacesImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& firstHalfPorousCrossings,
    const std::vector<PorousGridFaceCrossing>& secondHalfPorousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings);

PeriodicFlowStrangSubcyclingDiagnostics
advancePeriodicFlowStrangSspRk2SubcycledImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PeriodicFlowStrangSubcyclingSettings& settings);

} // namespace

PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowSettings& settings) {
    return advancePeriodicFlowImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        nullptr, nullptr, settings);
}

PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    const PeriodicFlowSettings& settings) {
    return advancePeriodicFlowImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        &pressureJumps, nullptr, settings);
}

PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowSettings& settings) {
    if (porousCrossings.empty()) {
        return advancePeriodicFlow(
            grid, velocityMetersPerSecond, pressurePascals, settings);
    }
    const PorousPeriodicProjectionInput input{
        &porousCrossings, nullptr, &porousIteration};
    return advancePeriodicFlowImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        nullptr, &input, settings);
}

PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowSettings& settings) {
    if (porousCrossings.empty() && prescribedPressureJumps.empty()) {
        if (!prescribedPressureJumps.matches(grid)) {
            throw std::invalid_argument(
                "periodic porous-flow jumps do not match their grid");
        }
        return advancePeriodicFlow(
            grid, velocityMetersPerSecond, pressurePascals, settings);
    }
    const PorousPeriodicProjectionInput input{
        &porousCrossings, &prescribedPressureJumps, &porousIteration};
    return advancePeriodicFlowImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        nullptr, &input, settings);
}

namespace {

PeriodicFlowDiagnostics advancePeriodicFlowImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PorousPeriodicProjectionInput* porousProjection,
    const PeriodicFlowSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || (pressureJumps != nullptr
            && !pressureJumps->matches(grid))
        || (porousProjection != nullptr
            && porousProjection->prescribedPressureJumps != nullptr
            && !porousProjection->prescribedPressureJumps->matches(grid))) {
        throw std::invalid_argument(
            "periodic flow fields or pressure jumps do not match their grid");
    }
    if (pressureJumps != nullptr && porousProjection != nullptr) {
        throw std::invalid_argument(
            "periodic flow cannot select two projection owners");
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
    bool projectionAccepted = false;
    if (porousProjection != nullptr) {
        PorousProjectionSettings coupledSettings;
        coupledSettings.projection = projectionSettings;
        coupledSettings.iteration = *porousProjection->iteration;
        if (porousProjection->prescribedPressureJumps == nullptr) {
            diagnostics.porousProjection =
                projectVelocityWithPorousInterfaces(
                    grid, candidateVelocity, candidatePressure,
                    *porousProjection->crossings, coupledSettings);
        } else {
            diagnostics.porousProjection =
                projectVelocityWithPorousInterfaces(
                    grid, candidateVelocity, candidatePressure,
                    *porousProjection->crossings,
                    *porousProjection->prescribedPressureJumps,
                    coupledSettings);
        }
        diagnostics.projection = diagnostics.porousProjection.projection;
        diagnostics.pressureJumpImpulseOnFluidNewtonSeconds =
            diagnostics.porousProjection
                .totalPressureJumpImpulseOnFluidNewtonSeconds;
        diagnostics.pressureJumpWorkToFluidJoules =
            diagnostics.porousProjection
                .totalPressureJumpWorkToFluidJoules;
        projectionAccepted = diagnostics.porousProjection.accepted;
    } else {
        diagnostics.projection = pressureJumps == nullptr
            ? projectVelocity(
                grid, candidateVelocity, candidatePressure,
                projectionSettings)
            : projectVelocityWithPressureJumps(
                grid, candidateVelocity, candidatePressure,
                *pressureJumps, projectionSettings);
        projectionAccepted = diagnostics.projection.converged;
    }
    if (!projectionAccepted) {
        diagnostics.failureStage =
            PeriodicFlowFailureStage::Projection;
        const bool advectionFinite = settings.advectionMode
                == PeriodicFlowAdvectionMode::PrescribedUniform
            ? diagnostics.uniformAdvection.finite
            : diagnostics.variableAdvection.finite;
        diagnostics.finite = advectionFinite
            && diffusionFinite
            && (porousProjection == nullptr
                || diagnostics.porousProjection.finite);
        return diagnostics;
    }

    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        subtract(
            diagnostics.momentumAfterNewtonSeconds,
            diagnostics.momentumBeforeNewtonSeconds),
        diagnostics.pressureJumpImpulseOnFluidNewtonSeconds);
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.projectionEnergyLossJoules =
        diagnostics.projection.kineticEnergyBeforeJoules
        + diagnostics.pressureJumpWorkToFluidJoules
        - diagnostics.projection.kineticEnergyAfterJoules;
    diagnostics.totalEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        + diagnostics.pressureJumpWorkToFluidJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.finalDivergenceL2PerSecond =
        diagnostics.projection.divergenceL2AfterPerSecond;
    const bool advectionFinite = settings.advectionMode
            == PeriodicFlowAdvectionMode::PrescribedUniform
        ? diagnostics.uniformAdvection.finite
        : diagnostics.variableAdvection.finite;
    diagnostics.finite = advectionFinite
        && diffusionFinite
        && (porousProjection == nullptr
            || diagnostics.porousProjection.finite)
        && finite(diagnostics.pressureJumpImpulseOnFluidNewtonSeconds)
        && std::isfinite(diagnostics.pressureJumpWorkToFluidJoules)
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
        length(add(
            diagnostics.momentumBeforeNewtonSeconds,
            diagnostics.pressureJumpImpulseOnFluidNewtonSeconds)),
        length(diagnostics.momentumAfterNewtonSeconds));
    const double energyTolerance = combinedTolerance(
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
        std::abs(diagnostics.kineticEnergyBeforeJoules
                 + diagnostics.pressureJumpWorkToFluidJoules),
        std::abs(diagnostics.kineticEnergyAfterJoules));
    diagnostics.accepted = diagnostics.finite
        && diagnostics.momentumResidualNormNewtonSeconds
            <= momentumTolerance
        && diagnostics.kineticEnergyAfterJoules
            <= diagnostics.kineticEnergyBeforeJoules
                + diagnostics.pressureJumpWorkToFluidJoules
                + energyTolerance;
    if (!diagnostics.accepted) {
        diagnostics.failureStage =
            PeriodicFlowFailureStage::Conservation;
        return diagnostics;
    }

    velocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace

PeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return advancePeriodicFlowStrangSspRk2Impl(
        grid, velocityMetersPerSecond, pressurePascals,
        nullptr, settings);
}

PeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return advancePeriodicFlowStrangSspRk2Impl(
        grid, velocityMetersPerSecond, pressurePascals,
        &pressureJumps, settings);
}

namespace {

PeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2Impl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || (pressureJumps != nullptr
            && !pressureJumps->matches(grid))) {
        throw std::invalid_argument(
            "periodic Strang-SSPRK2 fields or pressure jumps do not match their grid");
    }
    if (!isFinite(velocityMetersPerSecond) || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "periodic Strang-SSPRK2 fields must be finite");
    }

    PeriodicFlowStrangSspRk2Diagnostics diagnostics;
    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    CellScalarField initialDivergence(grid);
    computeDivergence(
        grid, velocityMetersPerSecond, initialDivergence);
    diagnostics.initialDivergenceL2PerSecond = l2Norm(initialDivergence);

    PeriodicMacDiffusionSettings halfDiffusionSettings;
    halfDiffusionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    halfDiffusionSettings.kinematicViscositySquareMetersPerSecond =
        settings.kinematicViscositySquareMetersPerSecond;
    halfDiffusionSettings.timeStepSeconds =
        0.5 * settings.timeStepSeconds;
    halfDiffusionSettings.maximumDiffusionNumber =
        settings.maximumDiffusionNumber;
    halfDiffusionSettings.absoluteMomentumToleranceNewtonSeconds =
        settings.absoluteMomentumToleranceNewtonSeconds;
    halfDiffusionSettings.relativeMomentumTolerance =
        settings.relativeMomentumTolerance;
    halfDiffusionSettings.absoluteEnergyToleranceJoules =
        settings.absoluteEnergyToleranceJoules;
    halfDiffusionSettings.relativeEnergyTolerance =
        settings.relativeEnergyTolerance;

    ProjectedMacAdvectionSspRk2Settings advectionSettings;
    advectionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    advectionSettings.timeStepSeconds = settings.timeStepSeconds;
    advectionSettings.reconstruction = settings.advectionReconstruction;
    advectionSettings.maximumLocalOutgoingCourantNumber =
        settings.maximumLocalOutgoingCourantNumber;
    advectionSettings.absoluteDivergenceTolerancePerSecond =
        settings.advectionAbsoluteDivergenceTolerancePerSecond;
    advectionSettings.relativeDivergenceTolerance =
        settings.advectionRelativeDivergenceTolerance;
    advectionSettings.projectionAbsoluteResidualTolerance =
        settings.projectionAbsoluteResidualTolerance;
    advectionSettings.projectionRelativeResidualTolerance =
        settings.projectionRelativeResidualTolerance;
    advectionSettings.projectionMaximumIterations =
        settings.projectionMaximumIterations;
    advectionSettings.absoluteMomentumToleranceNewtonSeconds =
        settings.absoluteMomentumToleranceNewtonSeconds;
    advectionSettings.relativeMomentumTolerance =
        settings.relativeMomentumTolerance;
    advectionSettings.absoluteEnergyToleranceJoules =
        settings.absoluteEnergyToleranceJoules;
    advectionSettings.relativeEnergyTolerance =
        settings.relativeEnergyTolerance;

    auto candidateVelocity = velocityMetersPerSecond;
    auto candidatePressure = pressurePascals;
    diagnostics.firstHalfDiffusion = diffuseVelocitySspRk2(
        grid, candidateVelocity, halfDiffusionSettings);
    if (!diagnostics.firstHalfDiffusion.accepted) {
        diagnostics.failureStage =
            PeriodicFlowStrangFailureStage::FirstHalfDiffusion;
        diagnostics.finite = diagnostics.firstHalfDiffusion.finite;
        return diagnostics;
    }
    diagnostics.firstHalfViscousEnergyLossJoules =
        diagnostics.firstHalfDiffusion.dissipatedKineticEnergyJoules;

    diagnostics.projectedAdvection = pressureJumps == nullptr
        ? advectVelocityProjectedSspRk2(
            grid, candidateVelocity, candidatePressure,
            advectionSettings)
        : advectVelocityProjectedSspRk2(
            grid, candidateVelocity, candidatePressure,
            *pressureJumps, advectionSettings);
    if (!diagnostics.projectedAdvection.accepted) {
        diagnostics.failureStage =
            PeriodicFlowStrangFailureStage::ProjectedAdvection;
        diagnostics.finite = diagnostics.firstHalfDiffusion.finite
            && diagnostics.projectedAdvection.finite;
        return diagnostics;
    }
    diagnostics.transportProjectionEnergyLossJoules =
        diagnostics.projectedAdvection.totalKineticEnergyLossJoules;

    diagnostics.secondHalfDiffusion = diffuseVelocitySspRk2(
        grid, candidateVelocity, halfDiffusionSettings);
    if (!diagnostics.secondHalfDiffusion.accepted) {
        diagnostics.failureStage =
            PeriodicFlowStrangFailureStage::SecondHalfDiffusion;
        diagnostics.finite = diagnostics.firstHalfDiffusion.finite
            && diagnostics.projectedAdvection.finite
            && diagnostics.secondHalfDiffusion.finite;
        return diagnostics;
    }
    diagnostics.secondHalfViscousEnergyLossJoules =
        diagnostics.secondHalfDiffusion.dissipatedKineticEnergyJoules;

    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        diagnostics.momentumAfterNewtonSeconds,
        diagnostics.momentumBeforeNewtonSeconds);
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.totalEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.maximumVelocityChangeMetersPerSecond = std::max({
        maximumDifference(
            velocityMetersPerSecond.xFaces(), candidateVelocity.xFaces()),
        maximumDifference(
            velocityMetersPerSecond.yFaces(), candidateVelocity.yFaces()),
        maximumDifference(
            velocityMetersPerSecond.zFaces(), candidateVelocity.zFaces()),
    });
    CellScalarField finalDivergence(grid);
    computeDivergence(grid, candidateVelocity, finalDivergence);
    diagnostics.finalDivergenceL2PerSecond = l2Norm(finalDivergence);
    diagnostics.finite = diagnostics.firstHalfDiffusion.finite
        && diagnostics.projectedAdvection.finite
        && diagnostics.secondHalfDiffusion.finite
        && isFinite(candidateVelocity) && isFinite(candidatePressure)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(
            diagnostics.firstHalfViscousEnergyLossJoules)
        && std::isfinite(
            diagnostics.transportProjectionEnergyLossJoules)
        && std::isfinite(
            diagnostics.secondHalfViscousEnergyLossJoules)
        && std::isfinite(diagnostics.totalEnergyLossJoules)
        && std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        && std::isfinite(diagnostics.initialDivergenceL2PerSecond)
        && std::isfinite(diagnostics.finalDivergenceL2PerSecond);
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
            PeriodicFlowStrangFailureStage::Conservation;
        return diagnostics;
    }

    velocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

PorousPeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithPorousInterfacesImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& firstHalfPorousCrossings,
    const std::vector<PorousGridFaceCrossing>& secondHalfPorousCrossings,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || (prescribedPressureJumps != nullptr
            && !prescribedPressureJumps->matches(grid))) {
        throw std::invalid_argument(
            "periodic porous Strang fields or pressure jumps do not match their grid");
    }
    if (!isFinite(velocityMetersPerSecond) || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "periodic porous Strang fields must be finite");
    }

    PorousPeriodicFlowStrangSspRk2Diagnostics diagnostics;
    if (firstHalfPorousCrossings.empty()
        && secondHalfPorousCrossings.empty()
        && (prescribedPressureJumps == nullptr
            || prescribedPressureJumps->empty())) {
        diagnostics.bulkFlow = advancePeriodicFlowStrangSspRk2(
            grid, velocityMetersPerSecond, pressurePascals, settings);
        diagnostics.momentumBeforeNewtonSeconds =
            diagnostics.bulkFlow.momentumBeforeNewtonSeconds;
        diagnostics.momentumAfterNewtonSeconds =
            diagnostics.bulkFlow.momentumAfterNewtonSeconds;
        diagnostics.momentumResidualNewtonSeconds =
            diagnostics.bulkFlow.momentumResidualNewtonSeconds;
        diagnostics.momentumResidualNormNewtonSeconds =
            diagnostics.bulkFlow.momentumResidualNormNewtonSeconds;
        diagnostics.kineticEnergyBeforeJoules =
            diagnostics.bulkFlow.kineticEnergyBeforeJoules;
        diagnostics.kineticEnergyAfterJoules =
            diagnostics.bulkFlow.kineticEnergyAfterJoules;
        diagnostics.totalEnergyLossJoules =
            diagnostics.bulkFlow.totalEnergyLossJoules;
        diagnostics.maximumVelocityChangeMetersPerSecond =
            diagnostics.bulkFlow.maximumVelocityChangeMetersPerSecond;
        diagnostics.initialDivergenceL2PerSecond =
            diagnostics.bulkFlow.initialDivergenceL2PerSecond;
        diagnostics.finalDivergenceL2PerSecond =
            diagnostics.bulkFlow.finalDivergenceL2PerSecond;
        diagnostics.finite = diagnostics.bulkFlow.finite;
        diagnostics.accepted = diagnostics.bulkFlow.accepted;
        if (!diagnostics.accepted) {
            diagnostics.failureStage =
                PorousPeriodicFlowStrangFailureStage::BulkFlow;
        }
        return diagnostics;
    }
    if (porousIteration.constitutiveEvaluation
        != PorousConstitutiveEvaluation::Midpoint) {
        throw std::invalid_argument(
            "periodic porous Strang coupling requires midpoint constitutive evaluation");
    }

    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    CellScalarField initialDivergence(grid);
    computeDivergence(
        grid, velocityMetersPerSecond, initialDivergence);
    diagnostics.initialDivergenceL2PerSecond = l2Norm(initialDivergence);

    PorousProjectionSettings halfPorousSettings;
    halfPorousSettings.iteration = porousIteration;
    halfPorousSettings.projection.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    halfPorousSettings.projection.timeStepSeconds =
        0.5 * settings.timeStepSeconds;
    halfPorousSettings.projection.absoluteResidualTolerance =
        settings.projectionAbsoluteResidualTolerance;
    halfPorousSettings.projection.relativeResidualTolerance =
        settings.projectionRelativeResidualTolerance;
    halfPorousSettings.projection.maximumIterations =
        settings.projectionMaximumIterations;

    auto candidateVelocity = velocityMetersPerSecond;
    auto candidatePressure = pressurePascals;
    diagnostics.firstHalfPorous = prescribedPressureJumps == nullptr
        ? projectVelocityWithPorousInterfaces(
            grid, candidateVelocity, candidatePressure,
            firstHalfPorousCrossings, halfPorousSettings)
        : projectVelocityWithPorousInterfaces(
            grid, candidateVelocity, candidatePressure,
            firstHalfPorousCrossings, *prescribedPressureJumps,
            halfPorousSettings);
    if (!diagnostics.firstHalfPorous.accepted) {
        diagnostics.failureStage =
            PorousPeriodicFlowStrangFailureStage::FirstHalfPorous;
        diagnostics.finite = diagnostics.firstHalfPorous.finite;
        return diagnostics;
    }

    diagnostics.bulkFlow = advancePeriodicFlowStrangSspRk2(
        grid, candidateVelocity, candidatePressure, settings);
    if (!diagnostics.bulkFlow.accepted) {
        diagnostics.failureStage =
            PorousPeriodicFlowStrangFailureStage::BulkFlow;
        diagnostics.finite = diagnostics.firstHalfPorous.finite
            && diagnostics.bulkFlow.finite;
        return diagnostics;
    }

    diagnostics.secondHalfPorous = prescribedPressureJumps == nullptr
        ? projectVelocityWithPorousInterfaces(
            grid, candidateVelocity, candidatePressure,
            secondHalfPorousCrossings, halfPorousSettings)
        : projectVelocityWithPorousInterfaces(
            grid, candidateVelocity, candidatePressure,
            secondHalfPorousCrossings, *prescribedPressureJumps,
            halfPorousSettings);
    if (!diagnostics.secondHalfPorous.accepted) {
        diagnostics.failureStage =
            PorousPeriodicFlowStrangFailureStage::SecondHalfPorous;
        diagnostics.finite = diagnostics.firstHalfPorous.finite
            && diagnostics.bulkFlow.finite
            && diagnostics.secondHalfPorous.finite;
        return diagnostics;
    }

    diagnostics.pressureJumpImpulseOnFluidNewtonSeconds = add(
        diagnostics.firstHalfPorous
            .totalPressureJumpImpulseOnFluidNewtonSeconds,
        diagnostics.secondHalfPorous
            .totalPressureJumpImpulseOnFluidNewtonSeconds);
    diagnostics.pressureJumpWorkToFluidJoules =
        diagnostics.firstHalfPorous.totalPressureJumpWorkToFluidJoules
        + diagnostics.secondHalfPorous.totalPressureJumpWorkToFluidJoules;
    diagnostics.porousDissipationJoules =
        diagnostics.firstHalfPorous.totalPorousDissipationJoules
        + diagnostics.secondHalfPorous.totalPorousDissipationJoules;
    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        diagnostics.momentumAfterNewtonSeconds,
        add(diagnostics.momentumBeforeNewtonSeconds,
            diagnostics.pressureJumpImpulseOnFluidNewtonSeconds));
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity, settings.densityKgPerCubicMeter);
    diagnostics.totalEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        + diagnostics.pressureJumpWorkToFluidJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.maximumVelocityChangeMetersPerSecond = std::max({
        maximumDifference(
            velocityMetersPerSecond.xFaces(), candidateVelocity.xFaces()),
        maximumDifference(
            velocityMetersPerSecond.yFaces(), candidateVelocity.yFaces()),
        maximumDifference(
            velocityMetersPerSecond.zFaces(), candidateVelocity.zFaces()),
    });
    CellScalarField finalDivergence(grid);
    computeDivergence(grid, candidateVelocity, finalDivergence);
    diagnostics.finalDivergenceL2PerSecond = l2Norm(finalDivergence);

    diagnostics.finite = diagnostics.firstHalfPorous.finite
        && diagnostics.bulkFlow.finite
        && diagnostics.secondHalfPorous.finite
        && isFinite(candidateVelocity) && isFinite(candidatePressure)
        && finite(diagnostics.pressureJumpImpulseOnFluidNewtonSeconds)
        && std::isfinite(diagnostics.pressureJumpWorkToFluidJoules)
        && std::isfinite(diagnostics.porousDissipationJoules)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.totalEnergyLossJoules)
        && std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        && std::isfinite(diagnostics.initialDivergenceL2PerSecond)
        && std::isfinite(diagnostics.finalDivergenceL2PerSecond);
    const double momentumTolerance = combinedTolerance(
        settings.absoluteMomentumToleranceNewtonSeconds,
        settings.relativeMomentumTolerance,
        length(add(
            diagnostics.momentumBeforeNewtonSeconds,
            diagnostics.pressureJumpImpulseOnFluidNewtonSeconds)),
        length(diagnostics.momentumAfterNewtonSeconds));
    const double energyTolerance = combinedTolerance(
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
        std::abs(diagnostics.kineticEnergyBeforeJoules
                 + diagnostics.pressureJumpWorkToFluidJoules),
        std::abs(diagnostics.kineticEnergyAfterJoules));
    diagnostics.accepted = diagnostics.finite
        && diagnostics.momentumResidualNormNewtonSeconds
            <= momentumTolerance
        && diagnostics.kineticEnergyAfterJoules
            <= diagnostics.kineticEnergyBeforeJoules
                + diagnostics.pressureJumpWorkToFluidJoules
                + energyTolerance;
    if (!diagnostics.accepted) {
        diagnostics.failureStage =
            PorousPeriodicFlowStrangFailureStage::Conservation;
        return diagnostics;
    }

    velocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

MovingPlanarPorousFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheetImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const MovingPlanarPorousSheetStrangStages& sheetStages,
    const SharpPressureJumpField* prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    const auto& first = sheetStages.firstHalf;
    const auto& second = sheetStages.secondHalf;
    if (first.surfaceStableId != second.surfaceStableId
        || first.minusRegionStableId != second.minusRegionStableId
        || first.plusRegionStableId != second.plusRegionStableId
        || first.resistance != second.resistance) {
        throw std::invalid_argument(
            "moving planar porous Strang stages change immutable sheet identity or material");
    }
    std::vector<PorousGridFaceCrossing> firstCrossings;
    std::vector<PorousGridFaceCrossing> secondCrossings;
    try {
        firstCrossings = makePlanarPorousSheetCrossings(
            grid, first);
        secondCrossings = makePlanarPorousSheetCrossings(
            grid, second);
    } catch (const std::runtime_error&) {
        throw std::invalid_argument(
            "moving planar porous Strang stage placement is outside its topology");
    }
    MovingPorousTopologySelection transition;
    try {
        transition = selectMovingPorousTopology(
            grid, first.topology,
            second.physicalPlaneCoordinateMeters);
    } catch (const std::runtime_error&) {
        throw std::invalid_argument(
            "moving planar porous Strang stages skip a topology segment");
    }
    if (transition.topology != second.topology) {
        throw std::invalid_argument(
            "moving planar porous Strang stages have discontinuous topology");
    }

    MovingPlanarPorousFlowStrangSspRk2Diagnostics diagnostics;
    diagnostics.firstHalfSheet = first;
    diagnostics.secondHalfSheet = second;
    diagnostics.flow =
        advancePeriodicFlowStrangSspRk2WithPorousInterfacesImpl(
            grid, velocityMetersPerSecond, pressurePascals,
            firstCrossings, secondCrossings,
            prescribedPressureJumps, porousIteration, settings);
    diagnostics.finite = diagnostics.flow.finite;
    diagnostics.accepted = diagnostics.flow.accepted;
    return diagnostics;
}

} // namespace

PorousPeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return advancePeriodicFlowStrangSspRk2WithPorousInterfacesImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        porousCrossings, porousCrossings,
        nullptr, porousIteration, settings);
}

PorousPeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return advancePeriodicFlowStrangSspRk2WithPorousInterfacesImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        porousCrossings, porousCrossings,
        &prescribedPressureJumps,
        porousIteration, settings);
}

MovingPlanarPorousFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheet(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const MovingPlanarPorousSheetStrangStages& sheetStages,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return
        advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheetImpl(
            grid, velocityMetersPerSecond, pressurePascals,
            sheetStages, nullptr, porousIteration, settings);
}

MovingPlanarPorousFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheet(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const MovingPlanarPorousSheetStrangStages& sheetStages,
    const SharpPressureJumpField& prescribedPressureJumps,
    const PorousIterationSettings& porousIteration,
    const PeriodicFlowStrangSspRk2Settings& settings) {
    return
        advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheetImpl(
            grid, velocityMetersPerSecond, pressurePascals,
            sheetStages, &prescribedPressureJumps,
            porousIteration, settings);
}

PeriodicFlowStrangSubcyclingDiagnostics
advancePeriodicFlowStrangSspRk2Subcycled(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowStrangSubcyclingSettings& settings) {
    return advancePeriodicFlowStrangSspRk2SubcycledImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        nullptr, settings);
}

PeriodicFlowStrangSubcyclingDiagnostics
advancePeriodicFlowStrangSspRk2Subcycled(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    const PeriodicFlowStrangSubcyclingSettings& settings) {
    return advancePeriodicFlowStrangSspRk2SubcycledImpl(
        grid, velocityMetersPerSecond, pressurePascals,
        &pressureJumps, settings);
}

namespace {

PeriodicFlowStrangSubcyclingDiagnostics
advancePeriodicFlowStrangSspRk2SubcycledImpl(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField* pressureJumps,
    const PeriodicFlowStrangSubcyclingSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || (pressureJumps != nullptr
            && !pressureJumps->matches(grid))) {
        throw std::invalid_argument(
            "periodic Strang subcycling fields or pressure jumps do not match their grid");
    }
    if (!isFinite(velocityMetersPerSecond) || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "periodic Strang subcycling fields must be finite");
    }

    PeriodicFlowStrangSubcyclingDiagnostics diagnostics;
    diagnostics.requestedIntervalSeconds =
        settings.flow.timeStepSeconds;
    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond,
        settings.flow.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond,
        settings.flow.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    CellScalarField initialDivergence(grid);
    computeDivergence(
        grid, velocityMetersPerSecond, initialDivergence);
    diagnostics.initialDivergenceL2PerSecond = l2Norm(initialDivergence);
    diagnostics.finalDivergenceL2PerSecond =
        diagnostics.initialDivergenceL2PerSecond;

    const Vector3 spacing = grid.cellSpacingMeters();
    const double inverseSpacingSquaredSum =
        1.0 / (spacing.x * spacing.x)
        + 1.0 / (spacing.y * spacing.y)
        + 1.0 / (spacing.z * spacing.z);
    const double oneStepHalfDiffusionNumber =
        0.5 * settings.flow.kinematicViscositySquareMetersPerSecond
        * settings.flow.timeStepSeconds
        * inverseSpacingSquaredSum;
    std::size_t substepCount = boundedSubstepCount(
        oneStepHalfDiffusionNumber
            / settings.flow.maximumDiffusionNumber,
        settings.maximumSubsteps);

    while (true) {
        diagnostics.plannedSubstepCount = substepCount;
        diagnostics.substepSeconds = settings.flow.timeStepSeconds
            / static_cast<double>(substepCount);
        if (substepCount > settings.maximumSubsteps) {
            diagnostics.failureStage =
                PeriodicFlowStrangSubcyclingFailureStage::SubstepLimit;
            return diagnostics;
        }

        auto candidateVelocity = velocityMetersPerSecond;
        auto candidatePressure = pressurePascals;
        diagnostics.substeps.clear();
        diagnostics.substeps.reserve(substepCount);
        diagnostics.completedSubstepCount = 0;
        bool restart = false;
        for (std::size_t substepIndex = 0;
             substepIndex < substepCount; ++substepIndex) {
            auto substepSettings = settings.flow;
            substepSettings.timeStepSeconds = diagnostics.substepSeconds;
            const auto substep = pressureJumps == nullptr
                ? advancePeriodicFlowStrangSspRk2(
                    grid, candidateVelocity, candidatePressure,
                    substepSettings)
                : advancePeriodicFlowStrangSspRk2(
                    grid, candidateVelocity, candidatePressure,
                    *pressureJumps, substepSettings);
            diagnostics.maximumObservedOutgoingCourantNumber = std::max(
                diagnostics.maximumObservedOutgoingCourantNumber,
                maximumOutgoingCourantNumber(substep));
            diagnostics.maximumObservedDiffusionNumber = std::max(
                diagnostics.maximumObservedDiffusionNumber,
                maximumDiffusionNumber(substep));
            if (!substep.accepted) {
                diagnostics.failedSubstepIndex = substepIndex;
                diagnostics.failedSubstep = substep;
                const double stabilityRatio =
                    rejectedStabilityRatio(substep);
                if (stabilityRatio > 1.0) {
                    std::size_t requiredCount = boundedSubstepCount(
                        static_cast<double>(substepCount)
                            * stabilityRatio,
                        settings.maximumSubsteps);
                    if (requiredCount <= substepCount) {
                        requiredCount = substepCount + 1;
                    }
                    if (requiredCount > settings.maximumSubsteps) {
                        diagnostics.plannedSubstepCount = requiredCount;
                        diagnostics.substepSeconds =
                            settings.flow.timeStepSeconds
                            / static_cast<double>(requiredCount);
                        diagnostics.failureStage =
                            PeriodicFlowStrangSubcyclingFailureStage::
                                SubstepLimit;
                        diagnostics.finite = substep.finite;
                        return diagnostics;
                    }
                    ++diagnostics.stabilityRetryCount;
                    substepCount = requiredCount;
                    restart = true;
                    break;
                }
                diagnostics.failureStage =
                    PeriodicFlowStrangSubcyclingFailureStage::Substep;
                diagnostics.finite = substep.finite;
                return diagnostics;
            }
            diagnostics.substeps.push_back(substep);
            diagnostics.completedSubstepCount =
                diagnostics.substeps.size();
        }
        if (restart) {
            continue;
        }

        diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
            grid, candidateVelocity,
            settings.flow.densityKgPerCubicMeter);
        diagnostics.momentumResidualNewtonSeconds = subtract(
            diagnostics.momentumAfterNewtonSeconds,
            diagnostics.momentumBeforeNewtonSeconds);
        diagnostics.momentumResidualNormNewtonSeconds = length(
            diagnostics.momentumResidualNewtonSeconds);
        diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
            grid, candidateVelocity,
            settings.flow.densityKgPerCubicMeter);
        diagnostics.totalEnergyLossJoules =
            diagnostics.kineticEnergyBeforeJoules
            - diagnostics.kineticEnergyAfterJoules;
        diagnostics.maximumVelocityChangeMetersPerSecond = std::max({
            maximumDifference(
                velocityMetersPerSecond.xFaces(),
                candidateVelocity.xFaces()),
            maximumDifference(
                velocityMetersPerSecond.yFaces(),
                candidateVelocity.yFaces()),
            maximumDifference(
                velocityMetersPerSecond.zFaces(),
                candidateVelocity.zFaces()),
        });
        CellScalarField finalDivergence(grid);
        computeDivergence(grid, candidateVelocity, finalDivergence);
        diagnostics.finalDivergenceL2PerSecond = l2Norm(finalDivergence);
        diagnostics.finite = isFinite(candidateVelocity)
            && isFinite(candidatePressure)
            && finite(diagnostics.momentumBeforeNewtonSeconds)
            && finite(diagnostics.momentumAfterNewtonSeconds)
            && finite(diagnostics.momentumResidualNewtonSeconds)
            && std::isfinite(
                diagnostics.momentumResidualNormNewtonSeconds)
            && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
            && std::isfinite(diagnostics.kineticEnergyAfterJoules)
            && std::isfinite(diagnostics.totalEnergyLossJoules)
            && std::isfinite(
                diagnostics.maximumVelocityChangeMetersPerSecond)
            && std::isfinite(diagnostics.initialDivergenceL2PerSecond)
            && std::isfinite(diagnostics.finalDivergenceL2PerSecond)
            && std::isfinite(
                diagnostics.maximumObservedOutgoingCourantNumber)
            && std::isfinite(
                diagnostics.maximumObservedDiffusionNumber);
        const double momentumTolerance = combinedTolerance(
            settings.flow.absoluteMomentumToleranceNewtonSeconds,
            settings.flow.relativeMomentumTolerance,
            length(diagnostics.momentumBeforeNewtonSeconds),
            length(diagnostics.momentumAfterNewtonSeconds));
        const double energyTolerance = combinedTolerance(
            settings.flow.absoluteEnergyToleranceJoules,
            settings.flow.relativeEnergyTolerance,
            std::abs(diagnostics.kineticEnergyBeforeJoules),
            std::abs(diagnostics.kineticEnergyAfterJoules));
        diagnostics.accepted = diagnostics.finite
            && diagnostics.momentumResidualNormNewtonSeconds
                <= momentumTolerance
            && diagnostics.kineticEnergyAfterJoules
                <= diagnostics.kineticEnergyBeforeJoules
                    + energyTolerance;
        if (!diagnostics.accepted) {
            diagnostics.failureStage =
                PeriodicFlowStrangSubcyclingFailureStage::Conservation;
            return diagnostics;
        }

        velocityMetersPerSecond = std::move(candidateVelocity);
        pressurePascals = std::move(candidatePressure);
        return diagnostics;
    }
}

} // namespace

} // namespace simwing::fsi::fluid
