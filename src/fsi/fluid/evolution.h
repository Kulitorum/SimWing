#pragma once

#include "fluid/advection.h"
#include "fluid/diffusion.h"
#include "fluid/projection.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t periodicLinearFlowStepVersion = 1;

enum class PeriodicLinearFlowFailureStage : std::uint8_t {
    None = 0,
    Advection = 1,
    Diffusion = 2,
    Projection = 3,
    Conservation = 4,
};

struct PeriodicLinearFlowSettings {
    double densityKgPerCubicMeter = 1.225;
    Vector3 transportVelocityMetersPerSecond;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    double timeStepSeconds = 1.0 / 60.0;
    double maximumTotalCourantNumber = 1.0;
    double maximumDiffusionNumber = 0.5;
    double projectionAbsoluteResidualTolerance = 1.0e-12;
    double projectionRelativeResidualTolerance = 1.0e-10;
    std::size_t projectionMaximumIterations = 1000;
    double absoluteMomentumToleranceNewtonSeconds = 2.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 2.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
};

struct PeriodicLinearFlowDiagnostics {
    std::uint32_t version = periodicLinearFlowStepVersion;
    UniformMacAdvectionDiagnostics advection;
    PeriodicMacDiffusionDiagnostics diffusion;
    ProjectionDiagnostics projection;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double advectionNumericalEnergyLossJoules = 0.0;
    double viscousEnergyLossJoules = 0.0;
    double projectionEnergyLossJoules = 0.0;
    double totalEnergyLossJoules = 0.0;
    double initialDivergenceL2PerSecond = 0.0;
    double finalDivergenceL2PerSecond = 0.0;
    PeriodicLinearFlowFailureStage failureStage =
        PeriodicLinearFlowFailureStage::None;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicLinearFlowDiagnostics&) const = default;
};

// First complete periodic verification step: bounded uniform-flow transport,
// explicit laminar viscosity, then the zero-mean pressure projection. It is a
// linear transport canonical, not nonlinear Navier-Stokes convection. All
// stages run on candidates; velocity and pressure commit together only after
// every stage and the final momentum/energy ledger are accepted.
[[nodiscard]] PeriodicLinearFlowDiagnostics advancePeriodicLinearFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicLinearFlowSettings& settings = {});

} // namespace simwing::fsi::fluid
