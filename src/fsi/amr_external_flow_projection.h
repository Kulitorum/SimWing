#pragma once

#include "amr_static_wing_interface.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace simwing::fsi::amr {

inline constexpr std::uint32_t amrWindTunnelProjectionVersion = 1;
inline constexpr std::uint32_t amrWindTunnelMomentumStepVersion = 1;
inline constexpr std::uint32_t amrStaticWingDirectForcingVersion = 1;

struct WindTunnelProjectionSettings {
    WindTunnelGridSettings grid;
    double timeStepSeconds = 0.005;
    double airDensityKilogramsPerCubicMeter = 1.225;
    double relativeTolerance = 1.0e-10;
    std::size_t maximumIterations = 200;
    std::size_t staticWingForcingProjectionIterations = 12;
    double staticWingDirectForcingRelaxation = 1.8;

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

// Owning X-fast diagnostic copy of the projected coarse level. The complete
// two-level solve remains represented by diagnostics; covered coarse pressure
// cells and faces are averaged down from the refined patch before publication.
// Face storage is intentionally not exposed through this visualization DTO.
struct WindTunnelProjectedCoarseGrid {
    WindTunnelProjectionDiagnostics diagnostics;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::vector<fluid::Vector3> velocityMetersPerSecond;
    std::vector<double> pressurePascals;
    std::vector<double> divergencePerSecond;
};

struct WindTunnelMomentumStepDiagnostics {
    std::uint32_t version = amrWindTunnelMomentumStepVersion;
    WindTunnelProjectionDiagnostics initialProjection;
    WindTunnelProjectionDiagnostics correctedProjection;
    double maximumOutgoingCourantNumber = 0.0;
    double maximumCellVelocityChangeMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    struct StaticWingDirectForcingDiagnostics {
        std::uint32_t version = amrStaticWingDirectForcingVersion;
        StaticWingInterfaceDiagnostics binding;
        double maximumSurfaceNormalSpeedBeforeMetersPerSecond = 0.0;
        double maximumSurfaceNormalSpeedAfterForcingMetersPerSecond = 0.0;
        double maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond = 0.0;
        std::size_t forcingProjectionIterations = 0;
        bool active = false;
        bool finite = true;
        bool accepted = true;

        bool operator==(const StaticWingDirectForcingDiagnostics&) const =
            default;
    } staticWing;
    bool finite = false;
    bool accepted = false;

    bool operator==(const WindTunnelMomentumStepDiagnostics&) const = default;
};

struct WindTunnelMomentumStepResult {
    WindTunnelMomentumStepDiagnostics diagnostics;
    WindTunnelProjectedCoarseGrid projectedCoarseGrid;
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

[[nodiscard]] WindTunnelProjectedCoarseGrid
evaluateWindTunnelProjectedCoarseGrid(
    const WindTunnelProjectionSettings& settings = {});

// Advances all face-centred velocity components through one first-order
// donor-cell convective predictor on both AMR levels, reapplies the physical
// velocity contract, synchronizes fine/coarse faces, and performs a fresh
// composite pressure projection. This remains an empty wind-tunnel canonical:
// no wing/interface force or turbulence model is present.
[[nodiscard]] WindTunnelMomentumStepResult
evaluateWindTunnelMomentumAdvance(
    const WindTunnelProjectionSettings& settings = {});

class WindTunnelMomentumState final {
public:
    explicit WindTunnelMomentumState(
        WindTunnelProjectionSettings settings = {});
    WindTunnelMomentumState(
        WindTunnelProjectionSettings settings,
        const Scene& staticWingScene);
    ~WindTunnelMomentumState();

    WindTunnelMomentumState(const WindTunnelMomentumState&) = delete;
    WindTunnelMomentumState& operator=(const WindTunnelMomentumState&) =
        delete;
    WindTunnelMomentumState(WindTunnelMomentumState&&) = delete;
    WindTunnelMomentumState& operator=(WindTunnelMomentumState&&) = delete;

    [[nodiscard]] WindTunnelMomentumStepResult advance();
    [[nodiscard]] const WindTunnelProjectedCoarseGrid&
    projectedCoarseGrid() const noexcept;
    [[nodiscard]] const StaticWingInterface*
    staticWingInterface() const noexcept;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace simwing::fsi::amr
