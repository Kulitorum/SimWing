#pragma once

#include "fluid/grid.h"

#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t uniformMacAdvectionVersion = 1;
inline constexpr std::uint32_t variableMacAdvectionVersion = 2;
inline constexpr std::uint32_t variableMacAdvectionSspRk2Version = 1;

enum class VariableMacReconstruction : std::uint8_t {
    DonorCell = 0,
    MonotonizedCentral = 1,
};

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

struct VariableMacAdvectionSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    VariableMacReconstruction reconstruction =
        VariableMacReconstruction::DonorCell;
    // Each staggered component control volume must export no more than its
    // complete old-time value during one donor-cell step.
    double maximumLocalOutgoingCourantNumber = 1.0;
    double absoluteDivergenceTolerancePerSecond = 1.0e-11;
    double relativeDivergenceTolerance = 1.0e-12;
    double absoluteMomentumToleranceNewtonSeconds = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
    // Donor-cell Euler stages are dissipative. A limited MUSCL Euler stage may
    // add O(dt^2) energy even though its enclosing SSPRK2 update does not; only
    // that enclosing integrator disables this intermediate-stage criterion.
    bool enforceEulerEnergyNonIncrease = true;
    double absoluteBoundToleranceMetersPerSecond = 1.0e-14;
    double relativeBoundTolerance = 1.0e-12;
};

struct VariableMacAdvectionDiagnostics {
    std::uint32_t version = variableMacAdvectionVersion;
    VariableMacReconstruction reconstruction =
        VariableMacReconstruction::DonorCell;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    double maximumLocalOutgoingCourantNumber = 0.0;
    double maximumAcceptedLocalOutgoingCourantNumber = 0.0;
    double maximumAdvectingDivergencePerSecond = 0.0;
    double maximumControlVolumeDivergencePerSecond = 0.0;
    double acceptedDivergenceTolerancePerSecond = 0.0;
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
    bool energyCriterionEnabled = true;
    bool energyNonIncreasing = false;
    bool uniformAdvector = false;
    bool divergenceCompatible = false;
    bool stable = false;
    bool bounded = false;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const VariableMacAdvectionDiagnostics&) const = default;
};

// Conservatively transports every periodic MAC component through the
// staggered control-volume fluxes induced by another divergence-free MAC
// field. Passing the transported field as advectingVelocity gives the first
// nonlinear self-advection path. Donor-cell reconstruction is first order; the
// monotonized-central option is a conservative limited second-order spatial
// stage intended for an SSPRK2 enclosure. A uniform donor-cell advector
// delegates to the exact uniform-flow oracle. No rejected request mutates
// velocity.
[[nodiscard]] VariableMacAdvectionDiagnostics advectVelocityByMacFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const MacVelocityField& advectingVelocityMetersPerSecond,
    const VariableMacAdvectionSettings& settings = {});

struct VariableMacAdvectionSspRk2Diagnostics {
    std::uint32_t version = variableMacAdvectionSspRk2Version;
    VariableMacReconstruction reconstruction =
        VariableMacReconstruction::DonorCell;
    VariableMacAdvectionDiagnostics firstEulerStage;
    VariableMacAdvectionDiagnostics secondEulerStage;
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
    bool bounded = false;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const VariableMacAdvectionSspRk2Diagnostics&) const = default;
};

// Applies two conservative variable-flow Euler stages with one fixed,
// divergence-free MAC advector and convexly averages the twice-advanced field
// with the old field. With monotonized-central reconstruction this supplies
// the bounded second-order temporal companion needed to measure smooth spatial
// refinement. Nonlinear self-advection instead requires pressure projection
// between stages and belongs to projected_advection.*. Failure is
// transactional.
[[nodiscard]] VariableMacAdvectionSspRk2Diagnostics
advectVelocityByMacFlowSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const MacVelocityField& advectingVelocityMetersPerSecond,
    const VariableMacAdvectionSettings& settings = {});

} // namespace simwing::fsi::fluid
