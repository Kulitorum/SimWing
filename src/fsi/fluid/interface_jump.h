#pragma once

#include "fluid/grid.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

enum class GridFaceAxis : std::uint8_t {
    X = 0,
    Y = 1,
    Z = 2,
};

// One sharp crossing along the cell-centre segment normal to a periodic MAC
// face. The indexed face separates the previous cell in axis direction from
// the cell at (i,j,k). crossingFraction is the open-interval position on that
// -axis to +axis segment. Region IDs and the signed pressure jump follow the
// same traversal:
//
//     pressureJumpPascals = p(plusRegion) - p(minusRegion)
//
// surfaceStableId identifies the oriented discrete-surface entity responsible
// for the crossing. Multiple crossings on one grid face remain individually
// visible, are ordered by crossingFraction, and must form a continuous region
// chain. The dense stencil stores their deterministic signed sum. The fraction
// establishes topology order only; this first-order operator does not claim
// subcell pressure interpolation or cut-cell geometry.
struct GridFacePressureJump {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    double pressureJumpPascals = 0.0;
    double crossingFraction = 0.5;

    bool operator==(const GridFacePressureJump&) const = default;
};

// Validated, immutable sharp-interface source data. Dense component arrays are
// retained for the projection stencil while canonical crossing metadata
// preserves stable surface, two-sided region identity, and normal ordering for
// diagnostics and later surface reconstruction.
class SharpPressureJumpField final {
public:
    SharpPressureJumpField(
        const PeriodicCartesianGrid& grid,
        std::vector<GridFacePressureJump> faces = {});

    [[nodiscard]] GridCellCounts cellCounts() const noexcept;
    [[nodiscard]] bool matches(const PeriodicCartesianGrid& grid) const noexcept;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t faceCount() const noexcept;
    [[nodiscard]] std::span<const GridFacePressureJump> faces() const noexcept;
    [[nodiscard]] std::span<const double> xFaceJumpsPascals() const noexcept;
    [[nodiscard]] std::span<const double> yFaceJumpsPascals() const noexcept;
    [[nodiscard]] std::span<const double> zFaceJumpsPascals() const noexcept;

    bool operator==(const SharpPressureJumpField&) const = default;

private:
    GridCellCounts cellCounts_;
    std::vector<GridFacePressureJump> faces_;
    std::vector<double> xFaceJumpsPascals_;
    std::vector<double> yFaceJumpsPascals_;
    std::vector<double> zFaceJumpsPascals_;
};

// Computes grad(p) with the prescribed discontinuity removed on every crossed
// face. A piecewise-constant pressure field satisfying every jump therefore
// has exactly zero sharp gradient away from floating-point roundoff.
void computePressureGradientWithJumps(
    const PeriodicCartesianGrid& grid,
    const CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    MacVelocityField& gradientPascalsPerMeter);

// Computes div(jump/distance), the finite-volume source paired with the sharp
// gradient. Projection subtracts this source from its Poisson right-hand side.
void computePressureJumpSource(
    const PeriodicCartesianGrid& grid,
    const SharpPressureJumpField& pressureJumps,
    CellScalarField& sourcePascalsPerSquareMeter);

} // namespace simwing::fsi::fluid
