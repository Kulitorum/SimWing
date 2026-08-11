#pragma once

#include "fluid/grid.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi::amr {

inline constexpr std::uint32_t amrWindTunnelBoundaryVersion = 1;

enum class PhysicalBoundaryKind : std::uint8_t {
    FarField = 0,
    PrescribedInflow = 1,
    PressureOutlet = 2,
};

// Physical-face order is -X, +X, -Y, +Y, -Z, +Z. The imported wing uses
// +Y chordwise flow, so the first wind-tunnel contract owns -Y as inflow and
// +Y as the pressure outlet. No physical face is periodic.
inline constexpr std::array<PhysicalBoundaryKind, 6>
    positiveYWindTunnelBoundaries{
        PhysicalBoundaryKind::FarField,
        PhysicalBoundaryKind::FarField,
        PhysicalBoundaryKind::PrescribedInflow,
        PhysicalBoundaryKind::PressureOutlet,
        PhysicalBoundaryKind::FarField,
        PhysicalBoundaryKind::FarField,
    };

struct WindTunnelGridSettings {
    fluid::GridCellCounts coarseCellCounts{32, 48, 16};
    fluid::Vector3 lowerMeters{-6.0, -4.0, -4.0};
    fluid::Vector3 upperMeters{6.0, 8.0, 2.0};
    fluid::Vector3 freestreamMetersPerSecond{0.0, 10.0, 0.0};
    std::size_t maximumGridSize = 16;
    std::size_t refinementRatio = 2;
    // Half-open coarse-cell box refined by refinementRatio. The default
    // follows the complete authoritative gnuC2 canopy plus a local wake band,
    // rather than refining an arbitrary cube around the domain centre.
    std::array<std::size_t, 3> refinedCoarseCellLower{1, 12, 1};
    std::array<std::size_t, 3> refinedCoarseCellUpperExclusive{31, 32, 12};

    bool operator==(const WindTunnelGridSettings&) const = default;
};

struct WindTunnelBoundaryDiagnostics {
    std::uint32_t version = amrWindTunnelBoundaryVersion;
    fluid::GridCellCounts coarseCellCounts;
    fluid::GridCellCounts refinedDomainCellCounts;
    fluid::Vector3 coarseCellSpacingMeters;
    std::array<PhysicalBoundaryKind, 6> physicalBoundaries{};
    std::size_t coarseBlockCount = 0;
    std::size_t refinedBlockCount = 0;
    std::size_t coarseValidCellCount = 0;
    std::size_t refinedValidCellCount = 0;
    std::size_t allocatedVelocityBytes = 0;
    double maximumInteriorWakePerturbationMetersPerSecond = 0.0;
    double lowerYInflowMaximumErrorMetersPerSecond = 0.0;
    double upperYOutflowGradientMaximumErrorMetersPerSecond = 0.0;
    double farFieldMaximumErrorMetersPerSecond = 0.0;
    bool allPhysicalDirectionsNonPeriodic = false;
    bool finite = false;
    bool accepted = false;

    bool operator==(const WindTunnelBoundaryDiagnostics&) const = default;
};

// Process-scoped AMReX runtime boundary. The wrapper keeps AMReX types out of
// public SimWing headers and finalizes only the instance it initialized.
class Runtime final {
public:
    Runtime(int& argc, char**& argv);
    ~Runtime();

    Runtime(const Runtime&) = delete;
    Runtime& operator=(const Runtime&) = delete;
    Runtime(Runtime&&) = delete;
    Runtime& operator=(Runtime&&) = delete;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

// Builds a two-level block-structured Cartesian wind-tunnel hierarchy,
// allocates face-centred velocity, initializes a nonuniform interior wake, and
// fills the first explicit physical ghost-cell contract. This is the AMR and
// boundary-ownership gate before advection/projection is admitted; it does not
// claim a Navier-Stokes step or wing aerodynamics.
[[nodiscard]] WindTunnelBoundaryDiagnostics
evaluateWindTunnelBoundaryInitialization(
    const WindTunnelGridSettings& settings = {});

} // namespace simwing::fsi::amr
