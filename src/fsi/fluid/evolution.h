#pragma once

#include "fluid/advection.h"
#include "fluid/diffusion.h"
#include "fluid/projection.h"
#include "fluid/projected_advection.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t periodicFlowStepVersion = 3;
inline constexpr std::uint32_t periodicFlowStrangSspRk2Version = 1;
inline constexpr std::uint32_t periodicFlowStrangSubcyclingVersion = 1;
inline constexpr std::size_t periodicFlowStrangMaximumSubsteps = 4096;

enum class PeriodicFlowAdvectionMode : std::uint8_t {
    PrescribedUniform = 0,
    SelfAdvectingMac = 1,
};

enum class PeriodicFlowDiffusionMode : std::uint8_t {
    ForwardEuler = 0,
    SspRk2 = 1,
};

enum class PeriodicFlowFailureStage : std::uint8_t {
    None = 0,
    Advection = 1,
    Diffusion = 2,
    Projection = 3,
    Conservation = 4,
};

struct PeriodicFlowSettings {
    double densityKgPerCubicMeter = 1.225;
    PeriodicFlowAdvectionMode advectionMode =
        PeriodicFlowAdvectionMode::PrescribedUniform;
    Vector3 transportVelocityMetersPerSecond;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    PeriodicFlowDiffusionMode diffusionMode =
        PeriodicFlowDiffusionMode::ForwardEuler;
    double timeStepSeconds = 1.0 / 60.0;
    double maximumTotalCourantNumber = 1.0;
    double advectionAbsoluteDivergenceTolerancePerSecond = 1.0e-11;
    double advectionRelativeDivergenceTolerance = 1.0e-12;
    double maximumDiffusionNumber = 0.5;
    double projectionAbsoluteResidualTolerance = 1.0e-12;
    double projectionRelativeResidualTolerance = 1.0e-10;
    std::size_t projectionMaximumIterations = 1000;
    double absoluteMomentumToleranceNewtonSeconds = 2.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 2.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
};

struct PeriodicFlowDiagnostics {
    std::uint32_t version = periodicFlowStepVersion;
    PeriodicFlowAdvectionMode advectionMode =
        PeriodicFlowAdvectionMode::PrescribedUniform;
    UniformMacAdvectionDiagnostics uniformAdvection;
    VariableMacAdvectionDiagnostics variableAdvection;
    PeriodicFlowDiffusionMode diffusionMode =
        PeriodicFlowDiffusionMode::ForwardEuler;
    PeriodicMacDiffusionDiagnostics explicitDiffusion;
    PeriodicMacDiffusionSspRk2Diagnostics sspRk2Diffusion;
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
    PeriodicFlowFailureStage failureStage =
        PeriodicFlowFailureStage::None;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicFlowDiagnostics&) const = default;
};

// Complete periodic verification step: bounded donor-cell transport, laminar
// viscosity, then the zero-mean pressure projection. Transport may be
// the exact prescribed-uniform oracle or first-order nonlinear MAC
// self-advection; viscosity may use forward Euler or second-order SSPRK2. All
// stages run on candidates; velocity and pressure commit together only after
// every stage and the final momentum/energy ledger pass.
[[nodiscard]] PeriodicFlowDiagnostics advancePeriodicFlow(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowSettings& settings = {});

enum class PeriodicFlowStrangFailureStage : std::uint8_t {
    None = 0,
    FirstHalfDiffusion = 1,
    ProjectedAdvection = 2,
    SecondHalfDiffusion = 3,
    Conservation = 4,
};

struct PeriodicFlowStrangSspRk2Settings {
    double densityKgPerCubicMeter = 1.225;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    double timeStepSeconds = 1.0 / 60.0;
    VariableMacReconstruction advectionReconstruction =
        VariableMacReconstruction::DonorCell;
    double maximumLocalOutgoingCourantNumber = 1.0;
    double advectionAbsoluteDivergenceTolerancePerSecond = 1.0e-11;
    double advectionRelativeDivergenceTolerance = 1.0e-12;
    double maximumDiffusionNumber = 0.5;
    double projectionAbsoluteResidualTolerance = 1.0e-12;
    double projectionRelativeResidualTolerance = 1.0e-10;
    std::size_t projectionMaximumIterations = 1000;
    double absoluteMomentumToleranceNewtonSeconds = 2.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 2.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
};

struct PeriodicFlowStrangSspRk2Diagnostics {
    std::uint32_t version = periodicFlowStrangSspRk2Version;
    PeriodicMacDiffusionSspRk2Diagnostics firstHalfDiffusion;
    ProjectedMacAdvectionSspRk2Diagnostics projectedAdvection;
    PeriodicMacDiffusionSspRk2Diagnostics secondHalfDiffusion;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double firstHalfViscousEnergyLossJoules = 0.0;
    double transportProjectionEnergyLossJoules = 0.0;
    double secondHalfViscousEnergyLossJoules = 0.0;
    double totalEnergyLossJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    double initialDivergenceL2PerSecond = 0.0;
    double finalDivergenceL2PerSecond = 0.0;
    PeriodicFlowStrangFailureStage failureStage =
        PeriodicFlowStrangFailureStage::None;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicFlowStrangSspRk2Diagnostics&) const = default;
};

// Second-order temporal verification path for periodic incompressible flow:
// SSPRK2 viscosity over half a step, pressure-projected nonlinear SSPRK2
// transport over the full step, and the symmetric viscous half step. All
// candidates remain private and pressure/velocity commit together only after
// every sub-integrator and the independent aggregate ledger pass. Donor-cell
// convection remains the default; limited monotonized-central reconstruction
// is selectable for the higher-order smooth-flow path.
[[nodiscard]] PeriodicFlowStrangSspRk2Diagnostics
advancePeriodicFlowStrangSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowStrangSspRk2Settings& settings = {});

enum class PeriodicFlowStrangSubcyclingFailureStage : std::uint8_t {
    None = 0,
    SubstepLimit = 1,
    Substep = 2,
    Conservation = 3,
};

struct PeriodicFlowStrangSubcyclingSettings {
    // flow.timeStepSeconds is the requested outer interval. Each private
    // Strang step receives that interval divided by the current subdivision.
    PeriodicFlowStrangSspRk2Settings flow;
    std::size_t maximumSubsteps = 1024;
};

struct PeriodicFlowStrangSubcyclingDiagnostics {
    std::uint32_t version = periodicFlowStrangSubcyclingVersion;
    double requestedIntervalSeconds = 0.0;
    double substepSeconds = 0.0;
    std::size_t plannedSubstepCount = 0;
    std::size_t completedSubstepCount = 0;
    std::size_t stabilityRetryCount = 0;
    // substeps contains only the final equal-step attempt. The last rejected
    // retry trigger remains available even after a later attempt succeeds.
    std::size_t failedSubstepIndex = 0;
    std::vector<PeriodicFlowStrangSspRk2Diagnostics> substeps;
    PeriodicFlowStrangSspRk2Diagnostics failedSubstep;
    double maximumObservedOutgoingCourantNumber = 0.0;
    double maximumObservedDiffusionNumber = 0.0;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double totalEnergyLossJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    double initialDivergenceL2PerSecond = 0.0;
    double finalDivergenceL2PerSecond = 0.0;
    PeriodicFlowStrangSubcyclingFailureStage failureStage =
        PeriodicFlowStrangSubcyclingFailureStage::None;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicFlowStrangSubcyclingDiagnostics&) const = default;
};

// Advances one requested outer interval through equal Strang/SSPRK2
// substeps. The initial subdivision satisfies the explicit viscous limit. A
// rejected advective CFL, limited maximum-principle, or viscous stability stage
// increases the subdivision and restarts the complete interval from private
// original candidates; projection, ledger, and other numerical failures are
// not hidden by retry. Velocity and pressure commit together only after every
// substep and an independent outer momentum/energy ledger pass.
[[nodiscard]] PeriodicFlowStrangSubcyclingDiagnostics
advancePeriodicFlowStrangSspRk2Subcycled(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const PeriodicFlowStrangSubcyclingSettings& settings = {});

} // namespace simwing::fsi::fluid
