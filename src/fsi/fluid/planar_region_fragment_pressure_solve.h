#pragma once

#include "fluid/planar_region_fragment_opening_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_operator.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentPressureSolveSettings {
    double absoluteResidualTolerancePascalsMeters = 1.0e-12;
    double relativeResidualTolerance = 1.0e-10;
    double absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-12;
    std::size_t maximumIterations = 4000;

    bool operator==(
        const PlanarPressureRegionFragmentPressureSolveSettings&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureSolveComponentDiagnostics {
    std::size_t componentIndex = 0;
    std::size_t fragmentCount = 0;
    double totalVolumeCubicMeters = 0.0;
    double rightHandSideSumPascalsMeters = 0.0;
    double compatibilityCorrectionPascalsMeters = 0.0;
    double correctionVolumeMeanBeforePascals = 0.0;
    double correctionVolumeMeanAfterPascals = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureSolveComponentDiagnostics&)
        const = default;
};

struct PlanarPressureRegionFragmentPressureSolveDiagnostics {
    bool compatible = false;
    bool converged = false;
    bool finite = false;
    bool usesOpeningPressureOperator = false;
    std::uint64_t pressureOperatorFingerprint = 0;
    std::uint64_t basePressureOperatorFingerprint = 0;
    std::uint64_t openingFingerprint = 0;
    std::uint64_t fragmentFingerprint = 0;
    std::size_t rowCount = 0;
    std::size_t componentCount = 0;
    std::size_t iterationCount = 0;
    double maximumAbsoluteComponentCompatibilityPascalsMeters = 0.0;
    double initialResidualL2PascalsMeters = 0.0;
    double finalResidualL2PascalsMeters = 0.0;
    double finalResidualMaximumPascalsMeters = 0.0;
    double maximumAbsoluteCorrectionVolumeMeanPascals = 0.0;
    std::vector<
        PlanarPressureRegionFragmentPressureSolveComponentDiagnostics>
        components;

    bool operator==(
        const PlanarPressureRegionFragmentPressureSolveDiagnostics&) const =
        default;
};

// Solves A*deltaP=b for a correction only. The static regional pressure
// potential retained by the fragment set is not part of this vector and is
// never shifted or overwritten. Integrated RHS compatibility is required per
// disconnected component; only admitted roundoff is removed. A converged
// correction is committed with roundoff-zero volume-weighted mean in each
// component. Any incompatible, non-finite, or non-converged attempt leaves the
// caller's warm start bit-for-bit unchanged.
[[nodiscard]] PlanarPressureRegionFragmentPressureSolveDiagnostics
solvePlanarPressureRegionFragmentPressureCorrection(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& correctionPascals,
    const PlanarPressureRegionFragmentPressureSolveSettings& settings = {});

// Opt-in counterpart over the exact opening-augmented graph. The same
// compatibility, constant-nullspace, volume-gauge, residual, and rollback
// contract applies, but compatibility is measured over opening-connected
// components. This solves only the supplied integrated RHS; it does not
// invent an opening source, aperture velocity, inertia, or resistance law.
[[nodiscard]] PlanarPressureRegionFragmentPressureSolveDiagnostics
solvePlanarPressureRegionFragmentOpeningPressureCorrection(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& correctionPascals,
    const PlanarPressureRegionFragmentPressureSolveSettings& settings = {});

} // namespace simwing::fsi::fluid
