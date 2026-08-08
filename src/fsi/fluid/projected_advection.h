#pragma once

#include "fluid/advection.h"
#include "fluid/projection.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t projectedMacAdvectionSspRk2Version = 2;

enum class ProjectedMacAdvectionFailureStage : std::uint8_t {
    None = 0,
    FirstAdvection = 1,
    FirstProjection = 2,
    SecondAdvection = 3,
    SecondProjection = 4,
    Conservation = 5,
};

struct ProjectedMacAdvectionSspRk2Settings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    VariableMacReconstruction reconstruction =
        VariableMacReconstruction::DonorCell;
    double maximumLocalOutgoingCourantNumber = 1.0;
    double absoluteDivergenceTolerancePerSecond = 1.0e-11;
    double relativeDivergenceTolerance = 1.0e-12;
    double projectionAbsoluteResidualTolerance = 1.0e-12;
    double projectionRelativeResidualTolerance = 1.0e-10;
    std::size_t projectionMaximumIterations = 1000;
    double absoluteMomentumToleranceNewtonSeconds = 2.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 2.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
};

struct ProjectedMacAdvectionSspRk2Diagnostics {
    std::uint32_t version = projectedMacAdvectionSspRk2Version;
    VariableMacReconstruction reconstruction =
        VariableMacReconstruction::DonorCell;
    VariableMacAdvectionDiagnostics firstAdvection;
    ProjectionDiagnostics firstProjection;
    VariableMacAdvectionDiagnostics secondAdvection;
    ProjectionDiagnostics secondProjection;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double totalKineticEnergyLossJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    double initialDivergenceL2PerSecond = 0.0;
    double finalDivergenceL2PerSecond = 0.0;
    ProjectedMacAdvectionFailureStage failureStage =
        ProjectedMacAdvectionFailureStage::None;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const ProjectedMacAdvectionSspRk2Diagnostics&) const = default;
};

// Advances the periodic incompressible nonlinear transport equation with the
// two-stage SSPRK2 method. Each donor-cell Euler stage self-advects from a
// divergence-free MAC state; pressure projection makes stage one eligible to
// act as stage two's advector and projects the final convex combination. Both
// fields commit together only after all four stages and the aggregate ledger
// pass. Donor-cell reconstruction remains the default; the selectable
// monotonized-central reconstruction is limited and conservative, with its
// intermediate Euler energy criterion deferred to this aggregate SSPRK2
// acceptance.
[[nodiscard]] ProjectedMacAdvectionSspRk2Diagnostics
advectVelocityProjectedSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const ProjectedMacAdvectionSspRk2Settings& settings = {});

// Static sharp crossings are applied at both SSPRK2 projection stages. The
// authored field is immutable over this interval; moving topology and
// pressure-dependent porous updates belong outside this operator. An empty
// field delegates bit-exactly to the original no-jump path.
[[nodiscard]] ProjectedMacAdvectionSspRk2Diagnostics
advectVelocityProjectedSspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    const ProjectedMacAdvectionSspRk2Settings& settings = {});

} // namespace simwing::fsi::fluid
