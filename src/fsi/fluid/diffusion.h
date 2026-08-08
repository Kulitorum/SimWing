#pragma once

#include "fluid/grid.h"

#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t periodicMacDiffusionVersion = 1;
inline constexpr std::uint32_t periodicMacDiffusionSspRk2Version = 1;

struct PeriodicMacDiffusionSettings {
    double densityKgPerCubicMeter = 1.225;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    double timeStepSeconds = 1.0 / 60.0;
    // Forward Euler with the centred seven-point stencil is non-amplifying
    // only while nu*dt*(1/dx^2+1/dy^2+1/dz^2) <= 1/2.
    double maximumDiffusionNumber = 0.5;
    double absoluteMomentumToleranceNewtonSeconds = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-12;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-12;
};

struct PeriodicMacDiffusionDiagnostics {
    std::uint32_t version = periodicMacDiffusionVersion;
    double densityKgPerCubicMeter = 0.0;
    double kinematicViscositySquareMetersPerSecond = 0.0;
    double timeStepSeconds = 0.0;
    Vector3 directionalDiffusionNumbers;
    double totalDiffusionNumber = 0.0;
    double maximumAcceptedDiffusionNumber = 0.0;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double dissipatedKineticEnergyJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    bool stable = false;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicMacDiffusionDiagnostics&) const = default;
};

// Applies one explicit centred viscous step to every unique periodic MAC-face
// component. Each staggered component uses the same translated cell lattice,
// so the scalar periodic seven-point Laplacian applies without interpolation.
// Unstable requests return accepted=false; invalid or non-finite input throws.
// In every failure mode the caller's velocity remains bit-for-bit unchanged.
[[nodiscard]] PeriodicMacDiffusionDiagnostics diffuseVelocityExplicit(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const PeriodicMacDiffusionSettings& settings = {});

struct PeriodicMacDiffusionSspRk2Diagnostics {
    std::uint32_t version = periodicMacDiffusionSspRk2Version;
    PeriodicMacDiffusionDiagnostics firstEulerStage;
    PeriodicMacDiffusionDiagnostics secondEulerStage;
    Vector3 momentumBeforeNewtonSeconds;
    Vector3 momentumAfterNewtonSeconds;
    Vector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double dissipatedKineticEnergyJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PeriodicMacDiffusionSspRk2Diagnostics&) const = default;
};

// Applies the two-stage second-order strong-stability-preserving Runge-Kutta
// method to the same periodic centred viscous operator. Both full-step Euler
// candidates must satisfy the existing sharp diffusion-number contract; the
// final state is the exact convex average of the old field and stage two.
// Stage or aggregate failure leaves the caller's velocity unchanged.
[[nodiscard]] PeriodicMacDiffusionSspRk2Diagnostics
diffuseVelocitySspRk2(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& velocityMetersPerSecond,
    const PeriodicMacDiffusionSettings& settings = {});

} // namespace simwing::fsi::fluid
