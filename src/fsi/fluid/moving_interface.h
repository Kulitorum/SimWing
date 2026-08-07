#pragma once

#include "fluid/interface_jump.h"
#include "fluid/projection.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t faceAlignedMovingInterfaceVersion = 1;

// A zero-thickness interface exactly coincident with one periodic MAC face.
// normalVelocityMetersPerSecond is signed in the positive coordinate-axis
// direction. The face separates the previous axis cell from (i,j,k), using the
// same orientation convention as GridFacePressureJump.
struct GridFaceMovingInterface {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    double normalVelocityMetersPerSecond = 0.0;

    bool operator==(const GridFaceMovingInterface&) const = default;
};

// Immutable fixed-topology verification boundary. Removing its constrained
// faces from the periodic cell graph must produce exactly one connected
// component for every stable fluid region. Multiple crossings on one grid face
// and disconnected aliases of one region ID are rejected.
class FaceAlignedMovingInterface final {
public:
    FaceAlignedMovingInterface(
        const PeriodicCartesianGrid& grid,
        std::vector<GridFaceMovingInterface> faces);

    [[nodiscard]] std::uint32_t version() const noexcept;
    [[nodiscard]] GridCellCounts cellCounts() const noexcept;
    [[nodiscard]] bool matches(const PeriodicCartesianGrid& grid) const noexcept;
    [[nodiscard]] std::size_t faceCount() const noexcept;
    [[nodiscard]] std::size_t regionCount() const noexcept;
    [[nodiscard]] std::span<const GridFaceMovingInterface> faces() const noexcept;
    [[nodiscard]] std::span<const std::uint64_t>
    regionStableIds() const noexcept;
    [[nodiscard]] std::span<const std::uint64_t>
    cellRegionStableIds() const noexcept;
    [[nodiscard]] std::span<const std::uint8_t>
    xFaceConstraints() const noexcept;
    [[nodiscard]] std::span<const std::uint8_t>
    yFaceConstraints() const noexcept;
    [[nodiscard]] std::span<const std::uint8_t>
    zFaceConstraints() const noexcept;
    [[nodiscard]] std::span<const double>
    xFaceNormalVelocitiesMetersPerSecond() const noexcept;
    [[nodiscard]] std::span<const double>
    yFaceNormalVelocitiesMetersPerSecond() const noexcept;
    [[nodiscard]] std::span<const double>
    zFaceNormalVelocitiesMetersPerSecond() const noexcept;

    bool operator==(const FaceAlignedMovingInterface&) const = default;

private:
    std::uint32_t version_ = faceAlignedMovingInterfaceVersion;
    GridCellCounts cellCounts_;
    Vector3 lowerMeters_;
    Vector3 upperMeters_;
    std::vector<GridFaceMovingInterface> faces_;
    std::vector<std::uint64_t> regionStableIds_;
    std::vector<std::uint64_t> cellRegionStableIds_;
    std::vector<std::uint8_t> xFaceConstraints_;
    std::vector<std::uint8_t> yFaceConstraints_;
    std::vector<std::uint8_t> zFaceConstraints_;
    std::vector<double> xFaceNormalVelocitiesMetersPerSecond_;
    std::vector<double> yFaceNormalVelocitiesMetersPerSecond_;
    std::vector<double> zFaceNormalVelocitiesMetersPerSecond_;
};

struct MovingInterfaceProjectionSettings {
    ProjectionSettings projection;
    // Each fixed-topology region must have zero net prescribed volume rate.
    // Nonzero volume change requires a future cut-cell geometric-conservation
    // term and is rejected rather than silently subtracting compatibility.
    double absoluteRegionVolumeRateToleranceCubicMetersPerSecond = 1.0e-12;
};

struct MovingFluidRegionDiagnostics {
    std::uint64_t stableId = 0;
    std::size_t cellCount = 0;
    double compatibilityVolumeRateCubicMetersPerSecond = 0.0;
    double pressureMeanBeforePascals = 0.0;
    double pressureMeanAfterPascals = 0.0;

    bool operator==(const MovingFluidRegionDiagnostics&) const = default;
};

struct MovingInterfaceSurfaceDiagnostics {
    std::uint64_t stableId = 0;
    std::size_t faceCount = 0;
    double areaSquareMeters = 0.0;
    // Adjacent cell-centre pressures are treated as the two one-sided values.
    // This is exact for the piecewise-constant slab canonical, but is not a
    // replacement for future higher-order interface reconstruction. Force is
    // fluid-on-interface: (p_minus - p_plus) times the positive-axis area
    // vector. Positive power/work is therefore energy delivered to the moving
    // interface.
    Vector3 pressureForceNewtons;
    Vector3 pressureImpulseNewtonSeconds;
    double pressurePowerWatts = 0.0;
    double pressureWorkJoules = 0.0;

    bool operator==(const MovingInterfaceSurfaceDiagnostics&) const = default;
};

struct MovingInterfaceProjectionDiagnostics {
    ProjectionDiagnostics projection;
    std::uint32_t interfaceVersion = faceAlignedMovingInterfaceVersion;
    std::size_t interfaceFaceCount = 0;
    std::size_t fluidRegionCount = 0;
    double maximumAbsoluteRegionVolumeRateCubicMetersPerSecond = 0.0;
    double maximumNormalVelocityErrorMetersPerSecond = 0.0;
    Vector3 totalPressureForceNewtons;
    Vector3 totalPressureImpulseNewtonSeconds;
    double totalPressurePowerWatts = 0.0;
    double totalPressureWorkJoules = 0.0;
    std::vector<MovingFluidRegionDiagnostics> regions;
    std::vector<MovingInterfaceSurfaceDiagnostics> surfaces;
    bool finite = true;

    bool operator==(const MovingInterfaceProjectionDiagnostics&) const = default;
};

// Transactionally projects a predicted velocity while holding every interface
// normal velocity exactly fixed. Pressure corrections use a disconnected
// finite-volume Laplacian and a zero-mean gauge in each stable fluid region;
// the caller's prior mean pressure in each region is retained. Those separate
// means are null modes of this incompressible solve, so their difference is an
// explicit scenario/coupling input rather than a pressure level discovered by
// this kernel. A failed or region-incompatible solve commits neither velocity
// nor pressure.
[[nodiscard]] MovingInterfaceProjectionDiagnostics
projectVelocityWithMovingInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionSettings& settings = {});

} // namespace simwing::fsi::fluid
