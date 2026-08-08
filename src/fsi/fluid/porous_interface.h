#pragma once

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

struct PorousProjectionSettings {
    ProjectionSettings projection;
    PorousConstitutiveEvaluation constitutiveEvaluation =
        PorousConstitutiveEvaluation::Endpoint;
    double absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-10;
    double relativeNormalVelocityTolerance = 1.0e-8;
    double absolutePressureJumpTolerancePascals = 1.0e-8;
    double relativePressureJumpTolerance = 1.0e-8;
    double relaxation = 0.5;
    std::size_t maximumNonlinearIterations = 100;

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
// is not yet moving cut-cell topology or a complete coupled flow integrator.
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

} // namespace simwing::fsi::fluid
