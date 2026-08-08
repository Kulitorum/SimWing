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

void validateSettings(const VariableMacAdvectionSettings& settings) {
    const std::array finiteValues{
        settings.densityKgPerCubicMeter,
        settings.timeStepSeconds,
        settings.maximumLocalOutgoingCourantNumber,
        settings.absoluteDivergenceTolerancePerSecond,
        settings.relativeDivergenceTolerance,
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
        || settings.maximumLocalOutgoingCourantNumber <= 0.0
        || settings.maximumLocalOutgoingCourantNumber > 1.0
        || settings.absoluteDivergenceTolerancePerSecond < 0.0
        || settings.relativeDivergenceTolerance < 0.0
        || settings.absoluteMomentumToleranceNewtonSeconds < 0.0
        || settings.relativeMomentumTolerance < 0.0
        || settings.absoluteEnergyToleranceJoules < 0.0
        || settings.relativeEnergyTolerance < 0.0
        || settings.absoluteBoundToleranceMetersPerSecond < 0.0
        || settings.relativeBoundTolerance < 0.0) {
        throw std::invalid_argument(
            "variable MAC advection settings are invalid");
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

enum class Axis : std::uint8_t {
    X = 0,
    Y = 1,
    Z = 2,
};

std::span<const double> component(const MacVelocityField& field,
                                  const Axis axis) noexcept {
    switch (axis) {
    case Axis::X:
        return field.xFaces();
    case Axis::Y:
        return field.yFaces();
    case Axis::Z:
        return field.zFaces();
    }
    return {};
}

std::span<double> component(MacVelocityField& field,
                            const Axis axis) noexcept {
    switch (axis) {
    case Axis::X:
        return field.xFaces();
    case Axis::Y:
        return field.yFaces();
    case Axis::Z:
        return field.zFaces();
    }
    return {};
}

std::size_t previous(const std::size_t value,
                     const std::size_t count) noexcept {
    return value == 0 ? count - 1 : value - 1;
}

std::size_t next(const std::size_t value,
                 const std::size_t count) noexcept {
    return value + 1 == count ? 0 : value + 1;
}

void shiftPositive(std::array<std::size_t, 3>& coordinates,
                   const GridCellCounts counts,
                   const Axis axis) noexcept {
    const std::array sizes{counts.x, counts.y, counts.z};
    const auto index = static_cast<std::size_t>(axis);
    coordinates[index] = next(coordinates[index], sizes[index]);
}

void shiftNegative(std::array<std::size_t, 3>& coordinates,
                   const GridCellCounts counts,
                   const Axis axis) noexcept {
    const std::array sizes{counts.x, counts.y, counts.z};
    const auto index = static_cast<std::size_t>(axis);
    coordinates[index] = previous(coordinates[index], sizes[index]);
}

std::size_t uncheckedIndex(
    const GridCellCounts counts,
    const std::array<std::size_t, 3>& coordinates) noexcept {
    return uncheckedIndex(
        counts, coordinates[0], coordinates[1], coordinates[2]);
}

double positiveFluxVelocity(
    const GridCellCounts counts,
    const MacVelocityField& advectingVelocity,
    const Axis transportedComponent,
    const Axis direction,
    const std::array<std::size_t, 3>& coordinates) noexcept {
    const auto advectingComponent = component(
        advectingVelocity, direction);
    if (transportedComponent == direction) {
        auto positive = coordinates;
        shiftPositive(positive, counts, direction);
        return 0.5 * (
            advectingComponent[uncheckedIndex(counts, coordinates)]
            + advectingComponent[uncheckedIndex(counts, positive)]);
    }

    // The transverse advecting component is native half a cell away from
    // this component control-volume face. Average the two native MAC samples
    // that straddle the transported component's coordinate.
    auto first = coordinates;
    shiftPositive(first, counts, direction);
    auto second = first;
    shiftNegative(second, counts, transportedComponent);
    return 0.5 * (
        advectingComponent[uncheckedIndex(counts, first)]
        + advectingComponent[uncheckedIndex(counts, second)]);
}

struct VariableComponentMetrics {
    double maximumOutgoingCourantNumber = 0.0;
    double maximumAbsoluteDivergencePerSecond = 0.0;
    double maximumChangeMetersPerSecond = 0.0;
};

void advectComponentByMacFlow(
    const GridCellCounts counts,
    const Vector3& spacing,
    const double timeStepSeconds,
    const MacVelocityField& advectingVelocity,
    const Axis transportedComponent,
    const std::span<const double> source,
    const std::span<double> destination,
    VariableComponentMetrics& metrics) {
    const std::array axes{Axis::X, Axis::Y, Axis::Z};
    const std::array inverseSpacing{
        1.0 / spacing.x, 1.0 / spacing.y, 1.0 / spacing.z};
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::array coordinates{i, j, k};
                const std::size_t center = uncheckedIndex(
                    counts, coordinates);
                double fluxDivergence = 0.0;
                double velocityDivergence = 0.0;
                double outgoingCourantNumber = 0.0;
                for (std::size_t directionIndex = 0;
                     directionIndex < axes.size(); ++directionIndex) {
                    const Axis direction = axes[directionIndex];
                    auto positiveCoordinates = coordinates;
                    shiftPositive(
                        positiveCoordinates, counts, direction);
                    auto negativeCoordinates = coordinates;
                    shiftNegative(
                        negativeCoordinates, counts, direction);
                    const double positiveVelocity = positiveFluxVelocity(
                        counts, advectingVelocity, transportedComponent,
                        direction, coordinates);
                    const double negativeVelocity = positiveFluxVelocity(
                        counts, advectingVelocity, transportedComponent,
                        direction, negativeCoordinates);
                    const double positiveValue = positiveVelocity >= 0.0
                        ? source[center]
                        : source[uncheckedIndex(
                            counts, positiveCoordinates)];
                    const double negativeValue = negativeVelocity >= 0.0
                        ? source[uncheckedIndex(
                            counts, negativeCoordinates)]
                        : source[center];
                    fluxDivergence += inverseSpacing[directionIndex]
                        * (positiveVelocity * positiveValue
                           - negativeVelocity * negativeValue);
                    velocityDivergence += inverseSpacing[directionIndex]
                        * (positiveVelocity - negativeVelocity);
                    outgoingCourantNumber += timeStepSeconds
                        * inverseSpacing[directionIndex]
                        * (std::max(positiveVelocity, 0.0)
                           + std::max(-negativeVelocity, 0.0));
                }
                destination[center] = source[center]
                    - timeStepSeconds * fluxDivergence;
                metrics.maximumOutgoingCourantNumber = std::max(
                    metrics.maximumOutgoingCourantNumber,
                    outgoingCourantNumber);
                metrics.maximumAbsoluteDivergencePerSecond = std::max(
                    metrics.maximumAbsoluteDivergencePerSecond,
                    std::abs(velocityDivergence));
                metrics.maximumChangeMetersPerSecond = std::max(
                    metrics.maximumChangeMetersPerSecond,
                    std::abs(destination[center] - source[center]));
            }
        }
    }
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

VariableMacAdvectionDiagnostics advectVelocityByMacFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const MacVelocityField& advectingVelocityMetersPerSecond,
    const VariableMacAdvectionSettings& settings) {
    validateSettings(settings);
    if (!velocityMetersPerSecond.matches(grid)
        || !advectingVelocityMetersPerSecond.matches(grid)) {
        throw std::invalid_argument(
            "variable MAC advection velocity does not match its grid");
    }
    if (!isFinite(velocityMetersPerSecond)
        || !isFinite(advectingVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "variable MAC advection velocity must be finite");
    }

    VariableMacAdvectionDiagnostics diagnostics;
    diagnostics.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    diagnostics.timeStepSeconds = settings.timeStepSeconds;
    diagnostics.maximumAcceptedLocalOutgoingCourantNumber =
        settings.maximumLocalOutgoingCourantNumber;
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

    const auto [advectingMinimumX, advectingMaximumX] = extrema(
        advectingVelocityMetersPerSecond.xFaces());
    const auto [advectingMinimumY, advectingMaximumY] = extrema(
        advectingVelocityMetersPerSecond.yFaces());
    const auto [advectingMinimumZ, advectingMaximumZ] = extrema(
        advectingVelocityMetersPerSecond.zFaces());
    const Vector3 spacing = grid.cellSpacingMeters();
    const double divergenceScale =
        std::max(std::abs(advectingMinimumX),
                 std::abs(advectingMaximumX)) / spacing.x
        + std::max(std::abs(advectingMinimumY),
                   std::abs(advectingMaximumY)) / spacing.y
        + std::max(std::abs(advectingMinimumZ),
                   std::abs(advectingMaximumZ)) / spacing.z;
    diagnostics.acceptedDivergenceTolerancePerSecond = combinedTolerance(
        settings.absoluteDivergenceTolerancePerSecond,
        settings.relativeDivergenceTolerance,
        divergenceScale, divergenceScale);
    CellScalarField advectingDivergence(grid);
    computeDivergence(
        grid, advectingVelocityMetersPerSecond, advectingDivergence);
    diagnostics.maximumAdvectingDivergencePerSecond =
        maximumAbsoluteValue(advectingDivergence);

    diagnostics.uniformAdvector =
        advectingMinimumX == advectingMaximumX
        && advectingMinimumY == advectingMaximumY
        && advectingMinimumZ == advectingMaximumZ;
    if (diagnostics.uniformAdvector) {
        UniformMacAdvectionSettings uniformSettings;
        uniformSettings.densityKgPerCubicMeter =
            settings.densityKgPerCubicMeter;
        uniformSettings.transportVelocityMetersPerSecond = {
            advectingMinimumX, advectingMinimumY, advectingMinimumZ};
        uniformSettings.timeStepSeconds = settings.timeStepSeconds;
        uniformSettings.maximumTotalCourantNumber =
            settings.maximumLocalOutgoingCourantNumber;
        uniformSettings.absoluteMomentumToleranceNewtonSeconds =
            settings.absoluteMomentumToleranceNewtonSeconds;
        uniformSettings.relativeMomentumTolerance =
            settings.relativeMomentumTolerance;
        uniformSettings.absoluteEnergyToleranceJoules =
            settings.absoluteEnergyToleranceJoules;
        uniformSettings.relativeEnergyTolerance =
            settings.relativeEnergyTolerance;
        uniformSettings.absoluteBoundToleranceMetersPerSecond =
            settings.absoluteBoundToleranceMetersPerSecond;
        uniformSettings.relativeBoundTolerance =
            settings.relativeBoundTolerance;
        const auto uniform = advectVelocityByUniformFlow(
            grid, velocityMetersPerSecond, uniformSettings);
        diagnostics.maximumLocalOutgoingCourantNumber =
            uniform.totalAbsoluteCourantNumber;
        diagnostics.maximumControlVolumeDivergencePerSecond = 0.0;
        diagnostics.momentumBeforeNewtonSeconds =
            uniform.momentumBeforeNewtonSeconds;
        diagnostics.momentumAfterNewtonSeconds =
            uniform.momentumAfterNewtonSeconds;
        diagnostics.momentumResidualNewtonSeconds =
            uniform.momentumResidualNewtonSeconds;
        diagnostics.momentumResidualNormNewtonSeconds =
            uniform.momentumResidualNormNewtonSeconds;
        diagnostics.kineticEnergyBeforeJoules =
            uniform.kineticEnergyBeforeJoules;
        diagnostics.kineticEnergyAfterJoules =
            uniform.kineticEnergyAfterJoules;
        diagnostics.numericalKineticEnergyLossJoules =
            uniform.numericalKineticEnergyLossJoules;
        diagnostics.componentMinimumBeforeMetersPerSecond =
            uniform.componentMinimumBeforeMetersPerSecond;
        diagnostics.componentMaximumBeforeMetersPerSecond =
            uniform.componentMaximumBeforeMetersPerSecond;
        diagnostics.componentMinimumAfterMetersPerSecond =
            uniform.componentMinimumAfterMetersPerSecond;
        diagnostics.componentMaximumAfterMetersPerSecond =
            uniform.componentMaximumAfterMetersPerSecond;
        diagnostics.maximumBoundViolationMetersPerSecond =
            uniform.maximumBoundViolationMetersPerSecond;
        diagnostics.maximumVelocityChangeMetersPerSecond =
            uniform.maximumVelocityChangeMetersPerSecond;
        diagnostics.divergenceCompatible = true;
        diagnostics.stable = uniform.stable;
        diagnostics.bounded = uniform.bounded;
        diagnostics.finite = uniform.finite;
        diagnostics.accepted = uniform.accepted;
        return diagnostics;
    }

    MacVelocityField candidate = velocityMetersPerSecond;
    VariableComponentMetrics componentMetrics;
    advectComponentByMacFlow(
        grid.cellCounts(), spacing, settings.timeStepSeconds,
        advectingVelocityMetersPerSecond, Axis::X,
        velocityMetersPerSecond.xFaces(), candidate.xFaces(),
        componentMetrics);
    advectComponentByMacFlow(
        grid.cellCounts(), spacing, settings.timeStepSeconds,
        advectingVelocityMetersPerSecond, Axis::Y,
        velocityMetersPerSecond.yFaces(), candidate.yFaces(),
        componentMetrics);
    advectComponentByMacFlow(
        grid.cellCounts(), spacing, settings.timeStepSeconds,
        advectingVelocityMetersPerSecond, Axis::Z,
        velocityMetersPerSecond.zFaces(), candidate.zFaces(),
        componentMetrics);
    diagnostics.maximumLocalOutgoingCourantNumber =
        componentMetrics.maximumOutgoingCourantNumber;
    diagnostics.maximumControlVolumeDivergencePerSecond =
        componentMetrics.maximumAbsoluteDivergencePerSecond;
    diagnostics.maximumVelocityChangeMetersPerSecond =
        componentMetrics.maximumChangeMetersPerSecond;
    diagnostics.divergenceCompatible =
        diagnostics.maximumAdvectingDivergencePerSecond
            <= diagnostics.acceptedDivergenceTolerancePerSecond
        && diagnostics.maximumControlVolumeDivergencePerSecond
            <= diagnostics.acceptedDivergenceTolerancePerSecond;
    diagnostics.stable =
        diagnostics.maximumLocalOutgoingCourantNumber
        <= settings.maximumLocalOutgoingCourantNumber;

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
    diagnostics.finite =
        std::isfinite(diagnostics.maximumLocalOutgoingCourantNumber)
        && std::isfinite(
            diagnostics.maximumAdvectingDivergencePerSecond)
        && std::isfinite(
            diagnostics.maximumControlVolumeDivergencePerSecond)
        && std::isfinite(
            diagnostics.acceptedDivergenceTolerancePerSecond)
        && isFinite(candidate)
        && finite(diagnostics.componentMinimumAfterMetersPerSecond)
        && finite(diagnostics.componentMaximumAfterMetersPerSecond)
        && std::isfinite(diagnostics.maximumBoundViolationMetersPerSecond)
        && finite(diagnostics.momentumBeforeNewtonSeconds)
        && finite(diagnostics.momentumAfterNewtonSeconds)
        && finite(diagnostics.momentumResidualNewtonSeconds)
        && std::isfinite(
            diagnostics.momentumResidualNormNewtonSeconds)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
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
    diagnostics.accepted = diagnostics.finite
        && diagnostics.divergenceCompatible
        && diagnostics.stable
        && diagnostics.bounded
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
