#include "fluid/projected_advection.h"

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

void validateSettings(
    const ProjectedMacAdvectionSspRk2Settings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.timeStepSeconds,
        settings.maximumLocalOutgoingCourantNumber,
        settings.absoluteDivergenceTolerancePerSecond,
        settings.relativeDivergenceTolerance,
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
        || settings.densityKgPerCubicMeter <= 0.0
        || settings.timeStepSeconds <= 0.0
        || settings.maximumLocalOutgoingCourantNumber <= 0.0
        || settings.maximumLocalOutgoingCourantNumber > 1.0
        || settings.absoluteDivergenceTolerancePerSecond < 0.0
        || settings.relativeDivergenceTolerance < 0.0
        || settings.projectionAbsoluteResidualTolerance < 0.0
        || settings.projectionRelativeResidualTolerance < 0.0
        || settings.projectionMaximumIterations == 0
        || settings.absoluteMomentumToleranceNewtonSeconds < 0.0
        || settings.relativeMomentumTolerance < 0.0
        || settings.absoluteEnergyToleranceJoules < 0.0
        || settings.relativeEnergyTolerance < 0.0) {
        throw std::invalid_argument(
            "projected MAC SSPRK2 advection settings are invalid");
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

void averageSspRk2Component(
    const std::span<const double> original,
    const std::span<const double> twiceAdvanced,
    const std::span<double> destination,
    double& maximumChange) {
    for (std::size_t index = 0; index < original.size(); ++index) {
        destination[index] = original[index]
            + 0.5 * (twiceAdvanced[index] - original[index]);
        maximumChange = std::max(
            maximumChange,
            std::abs(destination[index] - original[index]));
    }
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

VariableMacAdvectionSettings advectionSettings(
    const ProjectedMacAdvectionSspRk2Settings& settings) {
    VariableMacAdvectionSettings result;
    result.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    result.timeStepSeconds = settings.timeStepSeconds;
    result.maximumLocalOutgoingCourantNumber =
        settings.maximumLocalOutgoingCourantNumber;
    result.absoluteDivergenceTolerancePerSecond =
        settings.absoluteDivergenceTolerancePerSecond;
    result.relativeDivergenceTolerance =
        settings.relativeDivergenceTolerance;
    result.absoluteMomentumToleranceNewtonSeconds =
        settings.absoluteMomentumToleranceNewtonSeconds;
    result.relativeMomentumTolerance = settings.relativeMomentumTolerance;
    result.absoluteEnergyToleranceJoules =
        settings.absoluteEnergyToleranceJoules;
    result.relativeEnergyTolerance = settings.relativeEnergyTolerance;
    return result;
}

ProjectionSettings projectionSettings(
    const ProjectedMacAdvectionSspRk2Settings& settings) {
    ProjectionSettings result;
    result.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    result.timeStepSeconds = settings.timeStepSeconds;
    result.absoluteResidualTolerance =
        settings.projectionAbsoluteResidualTolerance;
    result.relativeResidualTolerance =
        settings.projectionRelativeResidualTolerance;
    result.maximumIterations = settings.projectionMaximumIterations;
    return result;
}

} // namespace

ProjectedMacAdvectionSspRk2Diagnostics
advectVelocityProjectedSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const ProjectedMacAdvectionSspRk2Settings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)) {
        throw std::invalid_argument(
            "projected MAC SSPRK2 fields do not match their grid");
    }
    if (!isFinite(velocityMetersPerSecond) || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "projected MAC SSPRK2 fields must be finite");
    }

    ProjectedMacAdvectionSspRk2Diagnostics diagnostics;
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

    const auto donorSettings = advectionSettings(settings);
    const auto pressureSettings = projectionSettings(settings);
    auto firstVelocity = velocityMetersPerSecond;
    auto firstPressure = pressurePascals;
    diagnostics.firstAdvection = advectVelocityByMacFlow(
        grid, firstVelocity, firstVelocity, donorSettings);
    if (!diagnostics.firstAdvection.accepted) {
        diagnostics.failureStage =
            ProjectedMacAdvectionFailureStage::FirstAdvection;
        diagnostics.finite = diagnostics.firstAdvection.finite;
        return diagnostics;
    }
    diagnostics.firstProjection = projectVelocity(
        grid, firstVelocity, firstPressure, pressureSettings);
    if (!diagnostics.firstProjection.converged) {
        diagnostics.failureStage =
            ProjectedMacAdvectionFailureStage::FirstProjection;
        diagnostics.finite = diagnostics.firstAdvection.finite;
        return diagnostics;
    }

    auto twiceAdvancedVelocity = firstVelocity;
    auto finalPressure = firstPressure;
    diagnostics.secondAdvection = advectVelocityByMacFlow(
        grid, twiceAdvancedVelocity, twiceAdvancedVelocity,
        donorSettings);
    if (!diagnostics.secondAdvection.accepted) {
        diagnostics.failureStage =
            ProjectedMacAdvectionFailureStage::SecondAdvection;
        diagnostics.finite = diagnostics.firstAdvection.finite
            && diagnostics.secondAdvection.finite;
        return diagnostics;
    }

    MacVelocityField candidate = velocityMetersPerSecond;
    averageSspRk2Component(
        velocityMetersPerSecond.xFaces(),
        twiceAdvancedVelocity.xFaces(), candidate.xFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    averageSspRk2Component(
        velocityMetersPerSecond.yFaces(),
        twiceAdvancedVelocity.yFaces(), candidate.yFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    averageSspRk2Component(
        velocityMetersPerSecond.zFaces(),
        twiceAdvancedVelocity.zFaces(), candidate.zFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    diagnostics.secondProjection = projectVelocity(
        grid, candidate, finalPressure, pressureSettings);
    if (!diagnostics.secondProjection.converged) {
        diagnostics.failureStage =
            ProjectedMacAdvectionFailureStage::SecondProjection;
        diagnostics.finite = diagnostics.firstAdvection.finite
            && diagnostics.secondAdvection.finite;
        return diagnostics;
    }
    diagnostics.maximumVelocityChangeMetersPerSecond = std::max({
        maximumDifference(
            velocityMetersPerSecond.xFaces(), candidate.xFaces()),
        maximumDifference(
            velocityMetersPerSecond.yFaces(), candidate.yFaces()),
        maximumDifference(
            velocityMetersPerSecond.zFaces(), candidate.zFaces()),
    });

    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidate, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        diagnostics.momentumAfterNewtonSeconds,
        diagnostics.momentumBeforeNewtonSeconds);
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidate, settings.densityKgPerCubicMeter);
    diagnostics.totalKineticEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.finalDivergenceL2PerSecond =
        diagnostics.secondProjection.divergenceL2AfterPerSecond;
    diagnostics.finite = diagnostics.firstAdvection.finite
        && diagnostics.secondAdvection.finite
        && isFinite(candidate) && isFinite(finalPressure)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.totalKineticEnergyLossJoules)
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
            ProjectedMacAdvectionFailureStage::Conservation;
        return diagnostics;
    }

    velocityMetersPerSecond = std::move(candidate);
    pressurePascals = std::move(finalPressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
