#pragma once

#include "fluid/grid.h"

#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t uniformMacAdvectionVersion = 1;

struct UniformMacAdvectionSettings {
    double densityKgPerCubicMeter = 1.225;
    Vector3 transportVelocityMetersPerSecond;
    double timeStepSeconds = 1.0 / 60.0;
    // The unsplit donor-cell update is a convex combination only while the
    // sum of absolute directional Courant numbers does not exceed one.
    double maximumTotalCourantNumber = 1.0;
    double absoluteMomentumToleranceNewtonSeconds = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
    double absoluteBoundToleranceMetersPerSecond = 1.0e-14;
    double relativeBoundTolerance = 1.0e-12;
};

struct UniformMacAdvectionDiagnostics {
    std::uint32_t version = uniformMacAdvectionVersion;
    double densityKgPerCubicMeter = 0.0;
    Vector3 transportVelocityMetersPerSecond;
    double timeStepSeconds = 0.0;
    Vector3 directionalCourantNumbers;
    double totalAbsoluteCourantNumber = 0.0;
    double maximumAcceptedTotalCourantNumber = 0.0;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double numericalKineticEnergyLossJoules = 0.0;
    Vector3 componentMinimumBeforeMetersPerSecond;
    Vector3 componentMaximumBeforeMetersPerSecond;
    Vector3 componentMinimumAfterMetersPerSecond;
    Vector3 componentMaximumAfterMetersPerSecond;
    double maximumBoundViolationMetersPerSecond = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    bool stable = false;
    bool bounded = false;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const UniformMacAdvectionDiagnostics&) const = default;
};

// Conservatively transports every periodic MAC component by one prescribed
// uniform world-space velocity using an unsplit donor-cell update. This is a
// bounded first-order verification oracle, not the later higher-order
// nonlinear self-advection operator. Unstable requests return accepted=false;
// invalid or non-finite input throws. No rejected request mutates velocity.
[[nodiscard]] UniformMacAdvectionDiagnostics advectVelocityByUniformFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const UniformMacAdvectionSettings& settings = {});

} // namespace simwing::fsi::fluid
