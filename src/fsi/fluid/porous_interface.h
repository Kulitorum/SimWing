#pragma once

#include "fluid/interface_jump.h"

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
// fractions and region chain are valid. This is an explicit flux-driven
// constitutive evaluation, not yet an implicit pressure/flux coupling solve.
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

} // namespace simwing::fsi::fluid
