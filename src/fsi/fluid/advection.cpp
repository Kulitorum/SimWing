#include "fluid/advection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

std::size_t uncheckedIndex(const GridCellCounts counts,
                           const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) noexcept {
    return (k * counts.y + j) * counts.x + i;
}

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

void validateSettings(const UniformMacAdvectionSettings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.transportVelocityMetersPerSecond.x,
        settings.transportVelocityMetersPerSecond.y,
        settings.transportVelocityMetersPerSecond.z,
        settings.timeStepSeconds,
        settings.maximumTotalCourantNumber,
        settings.absoluteMomentumToleranceNewtonSeconds,
        settings.relativeMomentumTolerance,
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
        settings.absoluteBoundToleranceMetersPerSecond,
        settings.relativeBoundTolerance,
    };
    if (!std::ranges::all_of(finiteValues, [](const double value) {
            return std::isfinite(value);
        })
        || settings.densityKgPerCubicMeter <= 0.0
        || settings.timeStepSeconds <= 0.0
        || settings.maximumTotalCourantNumber <= 0.0
        || settings.maximumTotalCourantNumber > 1.0
        || settings.absoluteMomentumToleranceNewtonSeconds < 0.0
        || settings.relativeMomentumTolerance < 0.0
        || settings.absoluteEnergyToleranceJoules < 0.0
        || settings.relativeEnergyTolerance < 0.0
        || settings.absoluteBoundToleranceMetersPerSecond < 0.0
        || settings.relativeBoundTolerance < 0.0) {
        throw std::invalid_argument(
            "uniform MAC advection settings are invalid");
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

std::pair<double, double> extrema(const std::span<const double> values) {
    const auto [minimum, maximum] = std::ranges::minmax_element(values);
    return {*minimum, *maximum};
}

void advectComponent(
    const GridCellCounts counts,
    const Vector3& courantMagnitudes,
    const Vector3& transportVelocity,
    const std::span<const double> source,
    const std::span<double> destination,
    double& maximumChange) {
    const double centerWeight = 1.0
        - courantMagnitudes.x
        - courantMagnitudes.y
        - courantMagnitudes.z;
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t upwindK = transportVelocity.z >= 0.0
            ? (k == 0 ? counts.z - 1 : k - 1)
            : (k + 1 == counts.z ? 0 : k + 1);
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t upwindJ = transportVelocity.y >= 0.0
                ? (j == 0 ? counts.y - 1 : j - 1)
                : (j + 1 == counts.y ? 0 : j + 1);
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t upwindI = transportVelocity.x >= 0.0
                    ? (i == 0 ? counts.x - 1 : i - 1)
                    : (i + 1 == counts.x ? 0 : i + 1);
                const std::size_t center = uncheckedIndex(counts, i, j, k);
                destination[center] =
                    centerWeight * source[center]
                    + courantMagnitudes.x
                        * source[uncheckedIndex(counts, upwindI, j, k)]
                    + courantMagnitudes.y
                        * source[uncheckedIndex(counts, i, upwindJ, k)]
                    + courantMagnitudes.z
                        * source[uncheckedIndex(counts, i, j, upwindK)];
                maximumChange = std::max(
                    maximumChange,
                    std::abs(destination[center] - source[center]));
            }
        }
    }
}

} // namespace

UniformMacAdvectionDiagnostics advectVelocityByUniformFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const UniformMacAdvectionSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)) {
        throw std::invalid_argument(
            "uniform MAC advection velocity does not match its grid");
    }
    if (!isFinite(velocityMetersPerSecond)) {
        throw std::invalid_argument(
            "uniform MAC advection velocity must be finite");
    }

    UniformMacAdvectionDiagnostics diagnostics;
    diagnostics.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    diagnostics.transportVelocityMetersPerSecond =
        settings.transportVelocityMetersPerSecond;
    diagnostics.timeStepSeconds = settings.timeStepSeconds;
    diagnostics.maximumAcceptedTotalCourantNumber =
        settings.maximumTotalCourantNumber;
    const Vector3 spacing = grid.cellSpacingMeters();
    diagnostics.directionalCourantNumbers = {
        settings.transportVelocityMetersPerSecond.x
            * settings.timeStepSeconds / spacing.x,
        settings.transportVelocityMetersPerSecond.y
            * settings.timeStepSeconds / spacing.y,
        settings.transportVelocityMetersPerSecond.z
            * settings.timeStepSeconds / spacing.z,
    };
    const Vector3 courantMagnitudes{
        std::abs(diagnostics.directionalCourantNumbers.x),
        std::abs(diagnostics.directionalCourantNumbers.y),
        std::abs(diagnostics.directionalCourantNumbers.z),
    };
    diagnostics.totalAbsoluteCourantNumber =
        courantMagnitudes.x + courantMagnitudes.y + courantMagnitudes.z;
    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    const auto [minimumX, maximumX] = extrema(
        velocityMetersPerSecond.xFaces());
    const auto [minimumY, maximumY] = extrema(
        velocityMetersPerSecond.yFaces());
    const auto [minimumZ, maximumZ] = extrema(
        velocityMetersPerSecond.zFaces());
    diagnostics.componentMinimumBeforeMetersPerSecond = {
        minimumX, minimumY, minimumZ};
    diagnostics.componentMaximumBeforeMetersPerSecond = {
        maximumX, maximumY, maximumZ};
    diagnostics.componentMinimumAfterMetersPerSecond =
        diagnostics.componentMinimumBeforeMetersPerSecond;
    diagnostics.componentMaximumAfterMetersPerSecond =
        diagnostics.componentMaximumBeforeMetersPerSecond;
    diagnostics.stable = diagnostics.totalAbsoluteCourantNumber
        <= settings.maximumTotalCourantNumber;
    diagnostics.finite = finite(diagnostics.directionalCourantNumbers)
        && std::isfinite(diagnostics.totalAbsoluteCourantNumber)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && finite(diagnostics.componentMinimumBeforeMetersPerSecond)
        && finite(diagnostics.componentMaximumBeforeMetersPerSecond);
    if (!diagnostics.stable || !diagnostics.finite) {
        return diagnostics;
    }
    if (diagnostics.totalAbsoluteCourantNumber == 0.0) {
        diagnostics.bounded = true;
        diagnostics.accepted = true;
        return diagnostics;
    }
    if (minimumX == maximumX
        && minimumY == maximumY
        && minimumZ == maximumZ) {
        diagnostics.bounded = true;
        diagnostics.accepted = true;
        return diagnostics;
    }

    MacVelocityField candidate = velocityMetersPerSecond;
    advectComponent(
        grid.cellCounts(), courantMagnitudes,
        settings.transportVelocityMetersPerSecond,
        velocityMetersPerSecond.xFaces(), candidate.xFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    advectComponent(
        grid.cellCounts(), courantMagnitudes,
        settings.transportVelocityMetersPerSecond,
        velocityMetersPerSecond.yFaces(), candidate.yFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    advectComponent(
        grid.cellCounts(), courantMagnitudes,
        settings.transportVelocityMetersPerSecond,
        velocityMetersPerSecond.zFaces(), candidate.zFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    const auto [candidateMinimumX, candidateMaximumX] = extrema(
        candidate.xFaces());
    const auto [candidateMinimumY, candidateMaximumY] = extrema(
        candidate.yFaces());
    const auto [candidateMinimumZ, candidateMaximumZ] = extrema(
        candidate.zFaces());
    diagnostics.componentMinimumAfterMetersPerSecond = {
        candidateMinimumX, candidateMinimumY, candidateMinimumZ};
    diagnostics.componentMaximumAfterMetersPerSecond = {
        candidateMaximumX, candidateMaximumY, candidateMaximumZ};
    diagnostics.maximumBoundViolationMetersPerSecond = std::max({
        0.0,
        minimumX - candidateMinimumX,
        minimumY - candidateMinimumY,
        minimumZ - candidateMinimumZ,
        candidateMaximumX - maximumX,
        candidateMaximumY - maximumY,
        candidateMaximumZ - maximumZ,
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
    diagnostics.numericalKineticEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.finite = diagnostics.finite
        && isFinite(candidate)
        && finite(diagnostics.componentMinimumAfterMetersPerSecond)
        && finite(diagnostics.componentMaximumAfterMetersPerSecond)
        && std::isfinite(diagnostics.maximumBoundViolationMetersPerSecond)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.numericalKineticEnergyLossJoules)
        && std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond);
    const double boundScale = std::max({
        std::abs(minimumX), std::abs(maximumX),
        std::abs(minimumY), std::abs(maximumY),
        std::abs(minimumZ), std::abs(maximumZ),
    });
    const double boundTolerance = combinedTolerance(
        settings.absoluteBoundToleranceMetersPerSecond,
        settings.relativeBoundTolerance, boundScale, boundScale);
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
    diagnostics.bounded = diagnostics.finite
        && diagnostics.maximumBoundViolationMetersPerSecond
            <= boundTolerance;
    diagnostics.accepted = diagnostics.bounded
        && diagnostics.momentumResidualNormNewtonSeconds
            <= momentumTolerance
        && diagnostics.kineticEnergyAfterJoules
            <= diagnostics.kineticEnergyBeforeJoules + energyTolerance;
    if (diagnostics.accepted) {
        velocityMetersPerSecond = std::move(candidate);
    }
    return diagnostics;
}

} // namespace simwing::fsi::fluid
