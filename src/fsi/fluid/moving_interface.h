#pragma once

#include "fluid/interface_jump.h"
#include "fluid/projection.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t faceAlignedMovingInterfaceVersion = 2;

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
// component for every stable fluid region. Distinct side IDs describe a
// separating closed sheet; equal side IDs describe a nonseparating sheet whose
// two sides remain connected around a resolved opening. Multiple crossings on
// one grid face and disconnected aliases of one region ID are rejected.
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

// Canonical face-resolved pressure sample. lower/upperCornerMeters describe
// the complete axis-aligned MAC-face tile; their coordinate on axis is equal.
// Adjacent cell-centre pressure is constant on this first-order tile, so its
// force acts exactly at the rectangle centroid. The tuple (axis,i,j,k) is its
// fixed-grid identity within the bound interface topology.
struct MovingInterfaceFaceDiagnostics {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    Vector3 lowerCornerMeters;
    Vector3 upperCornerMeters;
    double areaSquareMeters = 0.0;
    double normalVelocityMetersPerSecond = 0.0;
    Vector3 pressureTractionPascals;
    Vector3 pressureForceNewtons;
    double pressurePowerWatts = 0.0;
    // Full fluid-on-interface constraint reaction. The adjacent-cell pressure
    // contribution above is insufficient when imposing the face velocity also
    // changes the predicted MAC degree of freedom directly. The direct term is
    // zero for an already compatible prescribed velocity.
    Vector3 directConstraintForceNewtons;
    Vector3 constraintReactionTractionPascals;
    Vector3 constraintReactionForceNewtons;
    double constraintReactionPowerWatts = 0.0;

    bool operator==(const MovingInterfaceFaceDiagnostics&) const = default;
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
    // Maximum facewise pressure-traction departure from the surface mean.
    // A downstream uniform-traction bridge may accept the surface only when
    // this value meets its explicit reconstruction tolerance.
    double maximumPressureTractionDeviationPascals = 0.0;
    double pressurePowerWatts = 0.0;
    double pressureWorkJoules = 0.0;
    // The complete constraint reaction adds the impulse required to replace
    // the predicted normal MAC velocity with the prescribed value. This is the
    // load that must cross a coupled fluid/structure boundary. For an already
    // compatible predicted face velocity it reduces exactly to pressureForce.
    Vector3 directConstraintForceNewtons;
    Vector3 constraintReactionForceNewtons;
    Vector3 constraintReactionImpulseNewtonSeconds;
    double maximumConstraintReactionTractionDeviationPascals = 0.0;
    double constraintReactionPowerWatts = 0.0;
    double constraintReactionWorkJoules = 0.0;

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
    Vector3 totalConstraintReactionForceNewtons;
    Vector3 totalConstraintReactionImpulseNewtonSeconds;
    double totalConstraintReactionPowerWatts = 0.0;
    double totalConstraintReactionWorkJoules = 0.0;
    std::vector<MovingFluidRegionDiagnostics> regions;
    std::vector<MovingInterfaceFaceDiagnostics> faces;
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
// this kernel. Face diagnostics separate adjacent pressure traction from the
// complete constraint reaction, which also accounts for direct replacement of
// a predicted constrained-face velocity. A failed or region-incompatible solve
// commits neither velocity nor pressure.
[[nodiscard]] MovingInterfaceProjectionDiagnostics
projectVelocityWithMovingInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionSettings& settings = {});

// Applies the same disconnected moving-interface projection while retaining
// prescribed sharp pressure jumps on unconstrained faces. A grid face cannot
// be owned by both boundaries: an impermeable moving constraint and a
// pressure-jump crossing on the same MAC degree of freedom are rejected.
// Empty jump fields take the exact moving-interface-only path.
[[nodiscard]] MovingInterfaceProjectionDiagnostics
projectVelocityWithMovingInterfacesAndPressureJumps(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const SharpPressureJumpField& pressureJumps,
    const MovingInterfaceProjectionSettings& settings = {});

} // namespace simwing::fsi::fluid
