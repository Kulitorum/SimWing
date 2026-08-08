#include "fluid/diffusion.h"

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

void validateSettings(const PeriodicMacDiffusionSettings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.kinematicViscositySquareMetersPerSecond,
        settings.timeStepSeconds,
        settings.maximumDiffusionNumber,
        settings.absoluteMomentumToleranceNewtonSeconds,
        settings.relativeMomentumTolerance,
        settings.absoluteEnergyToleranceJoules,
        settings.relativeEnergyTolerance,
    };
    if (!std::ranges::all_of(finiteValues, [](const double value) {
            return std::isfinite(value);
        })
        || settings.densityKgPerCubicMeter <= 0.0
        || settings.kinematicViscositySquareMetersPerSecond < 0.0
        || settings.timeStepSeconds <= 0.0
        || settings.maximumDiffusionNumber <= 0.0
        || settings.maximumDiffusionNumber > 0.5
        || settings.absoluteMomentumToleranceNewtonSeconds < 0.0
        || settings.relativeMomentumTolerance < 0.0
        || settings.absoluteEnergyToleranceJoules < 0.0
        || settings.relativeEnergyTolerance < 0.0) {
        throw std::invalid_argument(
            "periodic MAC diffusion settings are invalid");
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

void diffuseComponent(const GridCellCounts counts,
                      const Vector3 inverseSpacingSquared,
                      const double viscosityTime,
                      const std::span<const double> source,
                      const std::span<double> destination,
                      double& maximumChange) {
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t previousK = k == 0 ? counts.z - 1 : k - 1;
        const std::size_t nextK = k + 1 == counts.z ? 0 : k + 1;
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t previousJ = j == 0 ? counts.y - 1 : j - 1;
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t previousI = i == 0 ? counts.x - 1 : i - 1;
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const std::size_t center = uncheckedIndex(counts, i, j, k);
                const double laplacian =
                    (source[uncheckedIndex(counts, previousI, j, k)]
                     - 2.0 * source[center]
                     + source[uncheckedIndex(counts, nextI, j, k)])
                        * inverseSpacingSquared.x
                    + (source[uncheckedIndex(counts, i, previousJ, k)]
                       - 2.0 * source[center]
                       + source[uncheckedIndex(counts, i, nextJ, k)])
                        * inverseSpacingSquared.y
                    + (source[uncheckedIndex(counts, i, j, previousK)]
                       - 2.0 * source[center]
                       + source[uncheckedIndex(counts, i, j, nextK)])
                        * inverseSpacingSquared.z;
                destination[center] = source[center]
                    + viscosityTime * laplacian;
                maximumChange = std::max(
                    maximumChange,
                    std::abs(destination[center] - source[center]));
            }
        }
    }
}

} // namespace

PeriodicMacDiffusionDiagnostics diffuseVelocityExplicit(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const PeriodicMacDiffusionSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)) {
        throw std::invalid_argument(
            "periodic MAC diffusion velocity does not match its grid");
    }
    if (!isFinite(velocityMetersPerSecond)) {
        throw std::invalid_argument(
            "periodic MAC diffusion velocity must be finite");
    }

    PeriodicMacDiffusionDiagnostics diagnostics;
    diagnostics.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    diagnostics.kinematicViscositySquareMetersPerSecond =
        settings.kinematicViscositySquareMetersPerSecond;
    diagnostics.timeStepSeconds = settings.timeStepSeconds;
    diagnostics.maximumAcceptedDiffusionNumber =
        settings.maximumDiffusionNumber;
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vector3 inverseSpacingSquared{
        1.0 / (spacing.x * spacing.x),
        1.0 / (spacing.y * spacing.y),
        1.0 / (spacing.z * spacing.z),
    };
    const double viscosityTime =
        settings.kinematicViscositySquareMetersPerSecond
        * settings.timeStepSeconds;
    diagnostics.directionalDiffusionNumbers = {
        viscosityTime * inverseSpacingSquared.x,
        viscosityTime * inverseSpacingSquared.y,
        viscosityTime * inverseSpacingSquared.z,
    };
    diagnostics.totalDiffusionNumber =
        diagnostics.directionalDiffusionNumbers.x
        + diagnostics.directionalDiffusionNumbers.y
        + diagnostics.directionalDiffusionNumbers.z;
    diagnostics.momentumBeforeNewtonSeconds = momentumNewtonSeconds(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.momentumAfterNewtonSeconds =
        diagnostics.momentumBeforeNewtonSeconds;
    diagnostics.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, velocityMetersPerSecond, settings.densityKgPerCubicMeter);
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyBeforeJoules;
    diagnostics.stable = diagnostics.totalDiffusionNumber
        <= settings.maximumDiffusionNumber;
    diagnostics.finite = finite(diagnostics.directionalDiffusionNumbers)
        && std::isfinite(diagnostics.totalDiffusionNumber)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules);
    if (!diagnostics.stable || !diagnostics.finite) {
        return diagnostics;
    }
    if (settings.kinematicViscositySquareMetersPerSecond == 0.0) {
        diagnostics.accepted = true;
        return diagnostics;
    }

    MacVelocityField candidate = velocityMetersPerSecond;
    diffuseComponent(
        grid.cellCounts(), inverseSpacingSquared, viscosityTime,
        velocityMetersPerSecond.xFaces(), candidate.xFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    diffuseComponent(
        grid.cellCounts(), inverseSpacingSquared, viscosityTime,
        velocityMetersPerSecond.yFaces(), candidate.yFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    diffuseComponent(
        grid.cellCounts(), inverseSpacingSquared, viscosityTime,
        velocityMetersPerSecond.zFaces(), candidate.zFaces(),
        diagnostics.maximumVelocityChangeMetersPerSecond);
    diagnostics.momentumAfterNewtonSeconds = momentumNewtonSeconds(
        grid, candidate, settings.densityKgPerCubicMeter);
    diagnostics.momentumResidualNewtonSeconds = subtract(
        diagnostics.momentumAfterNewtonSeconds,
        diagnostics.momentumBeforeNewtonSeconds);
    diagnostics.momentumResidualNormNewtonSeconds = length(
        diagnostics.momentumResidualNewtonSeconds);
    diagnostics.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidate, settings.densityKgPerCubicMeter);
    diagnostics.dissipatedKineticEnergyJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    diagnostics.finite = diagnostics.finite
        && isFinite(candidate)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.dissipatedKineticEnergyJoules)
        && std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond);
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
    if (diagnostics.accepted) {
        velocityMetersPerSecond = std::move(candidate);
    }
    return diagnostics;
}

} // namespace simwing::fsi::fluid
