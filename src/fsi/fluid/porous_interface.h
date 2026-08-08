#pragma once

#include "fluid/moving_interface.h"
#include "fluid/projection.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

// Calibrated normal-flow resistance per unit fabric area:
//
//   pressureDrop = linearResistance * q + quadraticResistance * |q| q
//
// where q is fluid velocity relative to the oriented sheet normal. The
// coefficients may come from a material coupon fit; this boundary deliberately
// does not invent permeability, thickness, or Forchheimer constants that are
// absent from the authoritative material data.
struct DarcyForchheimerResistance {
    double linearPascalSecondsPerMeter = 0.0;
    double quadraticPascalSecondsSquaredPerSquareMeter = 0.0;

    bool operator==(const DarcyForchheimerResistance&) const = default;
};

// Returns p_plus - p_minus. Positive relative normal flow travels from the
// minus region to the plus region and therefore creates a negative jump.
[[nodiscard]] double porousPressureJumpPascals(
    const DarcyForchheimerResistance& resistance,
    double relativeNormalVelocityMetersPerSecond);

// Exact monotone inverse of porousPressureJumpPascals().
[[nodiscard]] double porousRelativeNormalVelocityMetersPerSecond(
    const DarcyForchheimerResistance& resistance,
    double pressureJumpPascals);

struct PorousGridFaceCrossing {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    double crossingFraction = 0.5;
    double surfaceNormalVelocityMetersPerSecond = 0.0;
    DarcyForchheimerResistance resistance;

    bool operator==(const PorousGridFaceCrossing&) const = default;
};

struct PorousGridFaceSample {
    GridFacePressureJump pressureJump;
    double fluidNormalVelocityMetersPerSecond = 0.0;
    double surfaceNormalVelocityMetersPerSecond = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;
    double faceAreaSquareMeters = 0.0;
    double volumeFlowRateCubicMetersPerSecond = 0.0;
    double dissipationWatts = 0.0;

    bool operator==(const PorousGridFaceSample&) const = default;
};

// Samples a finite MAC field at authored porous crossings, applies the
// calibrated law, and exposes the result through the same canonical sharp-jump
// field used by projection. Several crossings may share one face when their
// fractions and region chain are valid. This is the explicit constitutive
// evaluation used by the coupled iteration below.
class PorousPressureJumpField final {
public:
    PorousPressureJumpField(
        const PeriodicCartesianGrid& grid,
        const MacVelocityField& fluidVelocityMetersPerSecond,
        std::vector<PorousGridFaceCrossing> crossings = {});

    [[nodiscard]] const SharpPressureJumpField&
    pressureJumps() const noexcept;
    [[nodiscard]] std::span<const PorousGridFaceSample>
    samples() const noexcept;
    [[nodiscard]] double totalDissipationWatts() const noexcept;

private:
    SharpPressureJumpField pressureJumps_;
    std::vector<PorousGridFaceSample> samples_;
    double totalDissipationWatts_ = 0.0;
};

enum class PorousConstitutiveEvaluation : std::uint8_t {
    Endpoint = 0,
    Midpoint = 1,
};

struct PorousIterationSettings {
    PorousConstitutiveEvaluation constitutiveEvaluation =
        PorousConstitutiveEvaluation::Endpoint;
    double absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-10;
    double relativeNormalVelocityTolerance = 1.0e-8;
    double absolutePressureJumpTolerancePascals = 1.0e-8;
    double relativePressureJumpTolerance = 1.0e-8;
    double relaxation = 0.5;
    std::size_t maximumNonlinearIterations = 100;

    bool operator==(const PorousIterationSettings&) const = default;
};

struct PorousProjectionSettings {
    ProjectionSettings projection;
    PorousIterationSettings iteration;

    bool operator==(const PorousProjectionSettings&) const = default;
};

struct PorousProjectionDiagnostics {
    bool accepted = false;
    bool finite = true;
    PorousConstitutiveEvaluation constitutiveEvaluation =
        PorousConstitutiveEvaluation::Endpoint;
    std::size_t nonlinearIterationCount = 0;
    std::size_t porousCrossingCount = 0;
    double initialMaximumNormalVelocityResidualMetersPerSecond = 0.0;
    double finalMaximumNormalVelocityResidualMetersPerSecond = 0.0;
    double finalMaximumPressureJumpResidualPascals = 0.0;
    double totalDissipationWatts = 0.0;
    double totalPorousDissipationJoules = 0.0;
    // Sum of every porous and separately prescribed oriented jump acting on
    // the fluid at the selected constitutive time. For one crossing,
    // F_fluid = (p_plus - p_minus) * area * positiveAxis.
    Vector3 totalPressureJumpForceOnFluidNewtons;
    Vector3 totalPressureJumpImpulseOnFluidNewtonSeconds;
    double totalPressureJumpPowerToFluidWatts = 0.0;
    double totalPressureJumpWorkToFluidJoules = 0.0;
    ProjectionDiagnostics projection;
    // Samples at constitutiveEvaluation; midpoint mode therefore deliberately
    // reports midpoint fluid/slip velocity rather than the committed endpoint.
    std::vector<PorousGridFaceSample> samples;

    bool operator==(const PorousProjectionDiagnostics&) const = default;
};

// Transactionally closes the Darcy-Forchheimer law with the periodic sharp
// projection. Each nonlinear iterate resamples the porous jump from either its
// relaxed endpoint field or the midpoint between the original prediction and
// that endpoint, then reprojects the original predicted velocity. The authored
// surfaceNormalVelocity is understood at the selected constitutive time. An
// accepted result must satisfy independent normal-velocity and jump residual
// tolerances. Optional prescribed jumps represent separately owned pressure
// sources or other static interfaces and participate in every trial. Failure
// commits neither velocity nor pressure. This fixed-grid Picard boundary
// supports nonuniform porous tiles and prescribed sheet-normal velocity, but it
// is not moving cut-cell topology. Complete fixed-grid macro-step composition
// is owned by fluid/evolution.h.
[[nodiscard]] PorousProjectionDiagnostics projectVelocityWithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const PorousProjectionSettings& settings = {});

[[nodiscard]] PorousProjectionDiagnostics projectVelocityWithPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const PorousProjectionSettings& settings = {});

struct MovingPorousProjectionSettings {
    MovingInterfaceProjectionSettings movingProjection;
    PorousIterationSettings iteration;

    bool operator==(const MovingPorousProjectionSettings&) const = default;
};

struct MovingPorousProjectionDiagnostics {
    PorousProjectionDiagnostics porous;
    MovingInterfaceProjectionDiagnostics movingInterface;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const MovingPorousProjectionDiagnostics&) const = default;
};

inline constexpr std::uint32_t porousSurfaceTractionVersion = 1;

struct PorousFaceTractionDiagnostics {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    double crossingFraction = 0.5;
    Vector3 lowerCornerMeters;
    Vector3 upperCornerMeters;
    double areaSquareMeters = 0.0;
    double pressureJumpPascals = 0.0;
    double fluidNormalVelocityMetersPerSecond = 0.0;
    double surfaceNormalVelocityMetersPerSecond = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;
    Vector3 pressureForceOnFluidNewtons;
    Vector3 pressureImpulseOnFluidNewtonSeconds;
    double pressurePowerToFluidWatts = 0.0;
    double pressureWorkToFluidJoules = 0.0;
    Vector3 pressureForceOnSurfaceNewtons;
    Vector3 pressureImpulseOnSurfaceNewtonSeconds;
    double pressurePowerToSurfaceWatts = 0.0;
    double pressureWorkToSurfaceJoules = 0.0;
    double dissipationWatts = 0.0;
    double dissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;

    bool operator==(const PorousFaceTractionDiagnostics&) const = default;
};

struct PorousSurfaceTractionAggregate {
    std::uint64_t stableId = 0;
    std::size_t faceCount = 0;
    double areaSquareMeters = 0.0;
    Vector3 pressureForceOnFluidNewtons;
    Vector3 pressureImpulseOnFluidNewtonSeconds;
    double pressurePowerToFluidWatts = 0.0;
    double pressureWorkToFluidJoules = 0.0;
    Vector3 pressureForceOnSurfaceNewtons;
    Vector3 pressureImpulseOnSurfaceNewtonSeconds;
    double pressurePowerToSurfaceWatts = 0.0;
    double pressureWorkToSurfaceJoules = 0.0;
    double dissipationWatts = 0.0;
    double dissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;

    bool operator==(const PorousSurfaceTractionAggregate&) const = default;
};

struct PorousSurfaceTractionDiagnostics {
    std::uint32_t version = porousSurfaceTractionVersion;
    std::vector<PorousFaceTractionDiagnostics> faces;
    std::vector<PorousSurfaceTractionAggregate> surfaces;
    Vector3 totalPressureForceOnFluidNewtons;
    Vector3 totalPressureImpulseOnFluidNewtonSeconds;
    double totalPressurePowerToFluidWatts = 0.0;
    double totalPressureWorkToFluidJoules = 0.0;
    Vector3 totalPressureForceOnSurfaceNewtons;
    Vector3 totalPressureImpulseOnSurfaceNewtonSeconds;
    double totalPressurePowerToSurfaceWatts = 0.0;
    double totalPressureWorkToSurfaceJoules = 0.0;
    double totalDissipationWatts = 0.0;
    double totalDissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;
    double maximumAbsoluteFaceEnergyResidualJoules = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PorousSurfaceTractionDiagnostics&) const = default;
};

// Closes the same endpoint/midpoint porous law through the disconnected
// moving-interface projector. Impermeable moving faces and porous/prescribed
// jump faces remain separately owned; any overlap is rejected by the combined
// projector before caller fields can change. The retained moving diagnostic is
// the last inner projection, so its complete face/surface reaction ledger
// corresponds to the accepted endpoint or the reported failure. Empty porous
// and prescribed topology delegates to the exact moving-only projection.
[[nodiscard]] MovingPorousProjectionDiagnostics
projectVelocityWithMovingAndPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& movingInterfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const MovingPorousProjectionSettings& settings = {});

[[nodiscard]] MovingPorousProjectionDiagnostics
projectVelocityWithMovingAndPorousInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& movingInterfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const MovingPorousProjectionSettings& settings = {});

// Reconstructs the equal-and-opposite pressure load carried by only the
// calibrated porous sheets. Separately prescribed jumps remain owned by their
// own source and are intentionally excluded. For every face and aggregate,
// fluid pressure work + sheet pressure work + porous dissipation closes to
// zero. The moving overload additionally requires complete outer and nested
// acceptance before exposing traction to a downstream structure adapter.
[[nodiscard]] PorousSurfaceTractionDiagnostics
evaluatePorousSurfaceTraction(
    const PeriodicCartesianGrid& grid,
    const PorousProjectionDiagnostics& diagnostics,
    double timeStepSeconds);

[[nodiscard]] PorousSurfaceTractionDiagnostics
evaluatePorousSurfaceTraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousProjectionDiagnostics& diagnostics,
    double timeStepSeconds);

} // namespace simwing::fsi::fluid
