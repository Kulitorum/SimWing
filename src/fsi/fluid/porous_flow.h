#pragma once

#include "fluid/porous_interface.h"

namespace simwing::fsi::fluid {

struct PorousPlugFlowSettings {
    DarcyForchheimerResistance resistance;
    double densityKgPerCubicMeter = 1.225;
    double flowLengthMeters = 1.0;
    double crossSectionAreaSquareMeters = 1.0;
    double drivingPressureRisePascals = 0.0;
    double timeStepSeconds = 1.0 / 60.0;

    bool operator==(const PorousPlugFlowSettings&) const = default;
};

struct PorousPlugFlowDiagnostics {
    bool accepted = false;
    double velocityBeforeMetersPerSecond = 0.0;
    double midpointVelocityMetersPerSecond = 0.0;
    double velocityAfterMetersPerSecond = 0.0;
    double midpointPressureDropPascals = 0.0;
    double endpointPressureDropPascals = 0.0;
    double fluidMassKilograms = 0.0;
    double momentumBeforeNewtonSeconds = 0.0;
    double momentumAfterNewtonSeconds = 0.0;
    double netPressureImpulseNewtonSeconds = 0.0;
    double momentumResidualNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double drivingPressureWorkJoules = 0.0;
    double porousDissipationJoules = 0.0;
    double energyResidualJoules = 0.0;

    bool operator==(const PorousPlugFlowDiagnostics&) const = default;
};

// Advances the uniform normal velocity of a fluid plug through one calibrated
// porous plane under a prescribed pressure rise. The Darcy-Forchheimer loss is
// evaluated at the implicit midpoint, giving exact discrete pressure-impulse
// and kinetic-energy identities up to independently checked roundoff. Failure
// leaves velocityMetersPerSecond unchanged. This is a one-degree-of-freedom
// verification oracle, not a replacement for a general grid coupling solve.
[[nodiscard]] PorousPlugFlowDiagnostics advancePorousPlugFlow(
    double& velocityMetersPerSecond,
    const PorousPlugFlowSettings& settings);

} // namespace simwing::fsi::fluid
