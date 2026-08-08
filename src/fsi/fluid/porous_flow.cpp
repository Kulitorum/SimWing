#include "fluid/porous_flow.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

double pressureDropPascals(
    const DarcyForchheimerResistance& resistance,
    const double velocityMetersPerSecond) {
    return -porousPressureJumpPascals(
        resistance, velocityMetersPerSecond);
}

double ledgerTolerance(const std::initializer_list<double> values) {
    double scale = 1.0;
    for (const double value : values) {
        scale = std::max(scale, std::abs(value));
    }
    return 512.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

PorousPlugFlowDiagnostics advancePorousPlugFlow(
    double& velocityMetersPerSecond,
    const PorousPlugFlowSettings& settings) {
    if (!std::isfinite(velocityMetersPerSecond)
        || !std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.flowLengthMeters)
        || !(settings.flowLengthMeters > 0.0)
        || !std::isfinite(settings.crossSectionAreaSquareMeters)
        || !(settings.crossSectionAreaSquareMeters > 0.0)
        || !std::isfinite(settings.drivingPressureRisePascals)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "porous plug flow requires finite positive physical settings");
    }
    // Validates the resistance without changing the exact zero state.
    static_cast<void>(porousPressureJumpPascals(
        settings.resistance, 0.0));

    const double inverseArealInertia = settings.timeStepSeconds
        / (settings.densityKgPerCubicMeter * settings.flowLengthMeters);
    const double midpointRightHandSide = velocityMetersPerSecond
        + 0.5 * inverseArealInertia
            * settings.drivingPressureRisePascals;
    const DarcyForchheimerResistance midpointResistance{
        1.0 + 0.5 * inverseArealInertia
            * settings.resistance.linearPascalSecondsPerMeter,
        0.5 * inverseArealInertia
            * settings.resistance
                .quadraticPascalSecondsSquaredPerSquareMeter,
    };
    if (!std::isfinite(inverseArealInertia)
        || !std::isfinite(midpointRightHandSide)
        || !std::isfinite(
            midpointResistance.linearPascalSecondsPerMeter)
        || !std::isfinite(
            midpointResistance
                .quadraticPascalSecondsSquaredPerSquareMeter)) {
        throw std::overflow_error(
            "porous plug-flow midpoint equation is not finite");
    }
    const double midpointVelocity =
        porousRelativeNormalVelocityMetersPerSecond(
            midpointResistance, -midpointRightHandSide);
    const double candidateVelocity =
        2.0 * midpointVelocity - velocityMetersPerSecond;
    if (!std::isfinite(candidateVelocity)) {
        throw std::overflow_error(
            "porous plug-flow candidate velocity is not finite");
    }

    PorousPlugFlowDiagnostics result;
    result.velocityBeforeMetersPerSecond = velocityMetersPerSecond;
    result.midpointVelocityMetersPerSecond = midpointVelocity;
    result.velocityAfterMetersPerSecond = candidateVelocity;
    result.midpointPressureDropPascals = pressureDropPascals(
        settings.resistance, midpointVelocity);
    result.endpointPressureDropPascals = pressureDropPascals(
        settings.resistance, candidateVelocity);
    result.fluidMassKilograms = settings.densityKgPerCubicMeter
        * settings.flowLengthMeters
        * settings.crossSectionAreaSquareMeters;
    result.momentumBeforeNewtonSeconds = result.fluidMassKilograms
        * velocityMetersPerSecond;
    result.momentumAfterNewtonSeconds = result.fluidMassKilograms
        * candidateVelocity;
    result.netPressureImpulseNewtonSeconds =
        settings.crossSectionAreaSquareMeters
        * (settings.drivingPressureRisePascals
           - result.midpointPressureDropPascals)
        * settings.timeStepSeconds;
    result.momentumResidualNewtonSeconds =
        result.momentumAfterNewtonSeconds
        - result.momentumBeforeNewtonSeconds
        - result.netPressureImpulseNewtonSeconds;
    result.kineticEnergyBeforeJoules = 0.5 * result.fluidMassKilograms
        * velocityMetersPerSecond * velocityMetersPerSecond;
    result.kineticEnergyAfterJoules = 0.5 * result.fluidMassKilograms
        * candidateVelocity * candidateVelocity;
    result.drivingPressureWorkJoules =
        settings.crossSectionAreaSquareMeters
        * settings.drivingPressureRisePascals
        * midpointVelocity * settings.timeStepSeconds;
    result.porousDissipationJoules =
        settings.crossSectionAreaSquareMeters
        * result.midpointPressureDropPascals
        * midpointVelocity * settings.timeStepSeconds;
    result.energyResidualJoules =
        result.kineticEnergyAfterJoules
        - result.kineticEnergyBeforeJoules
        - result.drivingPressureWorkJoules
        + result.porousDissipationJoules;

    if (!std::ranges::all_of(
            std::initializer_list<double>{
                result.midpointPressureDropPascals,
                result.endpointPressureDropPascals,
                result.fluidMassKilograms,
                result.momentumBeforeNewtonSeconds,
                result.momentumAfterNewtonSeconds,
                result.netPressureImpulseNewtonSeconds,
                result.momentumResidualNewtonSeconds,
                result.kineticEnergyBeforeJoules,
                result.kineticEnergyAfterJoules,
                result.drivingPressureWorkJoules,
                result.porousDissipationJoules,
                result.energyResidualJoules,
            },
            [](const double value) { return std::isfinite(value); })
        || result.porousDissipationJoules < 0.0) {
        throw std::overflow_error(
            "porous plug-flow diagnostics are not finite and dissipative");
    }
    const double momentumTolerance = ledgerTolerance({
        result.momentumBeforeNewtonSeconds,
        result.momentumAfterNewtonSeconds,
        result.netPressureImpulseNewtonSeconds,
    });
    const double energyTolerance = ledgerTolerance({
        result.kineticEnergyBeforeJoules,
        result.kineticEnergyAfterJoules,
        result.drivingPressureWorkJoules,
        result.porousDissipationJoules,
    });
    if (std::abs(result.momentumResidualNewtonSeconds)
            > momentumTolerance
        || std::abs(result.energyResidualJoules) > energyTolerance) {
        throw std::runtime_error(
            "porous plug-flow conservation ledger did not close");
    }

    result.accepted = true;
    velocityMetersPerSecond = candidateVelocity;
    return result;
}

} // namespace simwing::fsi::fluid
