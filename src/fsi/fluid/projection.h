#pragma once

#include "fluid/grid.h"
#include "fluid/interface_jump.h"

#include <cstddef>

namespace simwing::fsi::fluid {

struct ProjectionSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    double absoluteResidualTolerance = 1.0e-12;
    double relativeResidualTolerance = 1.0e-10;
    std::size_t maximumIterations = 1000;
};

struct ProjectionDiagnostics {
    bool converged = false;
    std::size_t iterationCount = 0;
    double compatibilityDivergencePerSecond = 0.0;
    double initialResidualPascalsPerSquareMeter = 0.0;
    double finalResidualPascalsPerSquareMeter = 0.0;
    double divergenceL2BeforePerSecond = 0.0;
    double divergenceL2AfterPerSecond = 0.0;
    double divergenceMaximumBeforePerSecond = 0.0;
    double divergenceMaximumAfterPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double pressureMeanPascals = 0.0;
    // Number of authored sharp crossings; several may share one grid face.
    std::size_t pressureJumpFaceCount = 0;
    double pressureJumpSourceCompatibilityPascalsPerSquareMeter = 0.0;

    bool operator==(const ProjectionDiagnostics&) const = default;
};

// Projects a periodic predicted velocity onto the discretely divergence-free
// subspace. Pressure is a zero-mean warm start and is overwritten only when
// the conjugate-gradient solve converges. A failed solve leaves both fields
// bit-for-bit unchanged so a caller can reduce its step and retry.
[[nodiscard]] ProjectionDiagnostics projectVelocity(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const ProjectionSettings& settings = {});

// Applies the same transactional projection with a prescribed sharp,
// two-sided pressure discontinuity. Empty jump fields take the exact no-jump
// path. Ordered multiple crossings on one face contribute their signed sum.
[[nodiscard]] ProjectionDiagnostics projectVelocityWithPressureJumps(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    const ProjectionSettings& settings = {});

} // namespace simwing::fsi::fluid
