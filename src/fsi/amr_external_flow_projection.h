#pragma once

#include "amr_external_flow.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::amr {

inline constexpr std::uint32_t amrWindTunnelProjectionVersion = 1;

struct WindTunnelProjectionSettings {
    WindTunnelGridSettings grid;
    double timeStepSeconds = 0.01;
    double airDensityKilogramsPerCubicMeter = 1.225;
    double relativeTolerance = 1.0e-10;
    std::size_t maximumIterations = 200;

    bool operator==(const WindTunnelProjectionSettings&) const = default;
};

struct WindTunnelProjectionDiagnostics {
    std::uint32_t version = amrWindTunnelProjectionVersion;
    WindTunnelBoundaryDiagnostics hierarchy;
    std::size_t activeCompositeCellCount = 0;
    std::size_t solverIterations = 0;
    double solverFinalResidual = 0.0;
    double initialMaximumDivergencePerSecond = 0.0;
    double projectedMaximumDivergencePerSecond = 0.0;
    double maximumDivergenceReductionRatio = 0.0;
    double maximumPressureCorrectionPascals = 0.0;
    double lowerYInflowNormalVelocityErrorMetersPerSecond = 0.0;
    double upperYOutletNormalVelocityChangeMetersPerSecond = 0.0;
    bool pressureOutletReferenceOwned = false;
    bool finite = false;
    bool accepted = false;

    bool operator==(const WindTunnelProjectionDiagnostics&) const = default;
};

// Projects a deliberately divergent, nonuniform positive-Y wake on the
// two-level non-periodic hierarchy. Homogeneous pressure-correction Neumann
// conditions preserve prescribed normal velocity at -Y and the X/Z far
// field; homogeneous Dirichlet at +Y owns the pressure reference and allows
// the outlet flux to adjust. This is a projection canonical, not yet an
// advective Navier-Stokes step or a wing-aerodynamics result.
[[nodiscard]] WindTunnelProjectionDiagnostics
evaluateWindTunnelPressureProjection(
    const WindTunnelProjectionSettings& settings = {});

} // namespace simwing::fsi::amr
