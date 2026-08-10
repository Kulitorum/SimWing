#pragma once

#include "fluid/planar_region_fragment_opening_flux.h"
#include "fluid/planar_region_fragment_pressure_solve.h"
#include "fluid/planar_region_fragment_volume_rate.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentOpeningPressureProjectionSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    double absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeContinuityTolerance = 1.0e-10;
    double absoluteMomentumResidualToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeMomentumResidualTolerance = 1.0e-10;
    double absoluteEnergyResidualToleranceJoules = 1.0e-12;
    double relativeEnergyResidualTolerance = 1.0e-10;
    PlanarPressureRegionFragmentPressureSolveSettings pressureSolve;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureProjectionSettings&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningPressureProjectionLimits {
    PlanarPressureRegionFragmentOpeningPressureOperatorLimits
        pressureOperatorLimits;
    PlanarPressureRegionFragmentVolumeRateLimits volumeRateLimits;
    PlanarPressureRegionFragmentOpeningFluxLimits openingFluxLimits;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics {
    bool accepted = false;
    bool finite = false;
    bool energyAccepted = false;
    std::uint64_t pressureOperatorFingerprint = 0;
    std::uint64_t basePressureOperatorFingerprint = 0;
    std::uint64_t topologyFingerprint = 0;
    std::uint64_t fragmentFingerprint = 0;
    std::uint64_t volumeRateFingerprint = 0;
    std::uint64_t openingFingerprint = 0;
    std::uint64_t predictedOpeningFluxFingerprint = 0;
    std::uint64_t correctedOpeningFluxFingerprint = 0;
    std::size_t fragmentCount = 0;
    std::size_t topologyLinkCount = 0;
    std::size_t openingPatchCount = 0;
    std::size_t correctedSameRegionGridLinkCount = 0;
    std::size_t retainedPressureLayerWallLinkCount = 0;
    std::size_t workingStorageBytes = 0;
    double predictedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double predictedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double maximumAbsolutePredictedComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double correctedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double correctedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double continuityToleranceCubicMetersPerSecond = 0.0;
    double maximumAbsoluteGridVelocityCorrectionMetersPerSecond = 0.0;
    double maximumAbsoluteOpeningVelocityCorrectionMetersPerSecond = 0.0;
    double maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        0.0;
    double momentumResidualToleranceKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double midpointPressureWorkJoules = 0.0;
    double workEnergyResidualJoules = 0.0;
    double maximumAbsoluteWorkEnergyResidualJoules = 0.0;
    double correctionKineticEnergyJoules = 0.0;
    double geometryPressureWorkJoules = 0.0;
    double affineEnergyResidualJoules = 0.0;
    double energyResidualToleranceJoules = 0.0;
    PlanarPressureRegionFragmentPressureSolveDiagnostics pressureSolve;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics&)
        const = default;
};

// Opt-in moving projection over the exact aperture-augmented graph. Predicted
// same-region Cartesian velocities and material-relative aperture velocities
// jointly enter dV/dt + div(Q). A compatible pressure correction updates both
// kinds of degree with dt/rho * (p_minus-p_plus)/centerDistance, while the
// remaining pressure-layer wall velocity stays exactly zero.
//
// The diagonal mass rho*area*centerDistance supplies a per-degree pressure-
// impulse/work identity. Acceptance also requires the global moving affine
// identity deltaK = geometryPressureWork - correctionKineticEnergy. Link
// velocities, opening samples, immutable opening-flux state, and the pressure
// warm start commit together only after continuity and energy close. This is
// an inviscid projection; it owns no aperture resistance/loss law, authored
// static-pressure relaxation, traction subtraction, or production state.
[[nodiscard]]
PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics
projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::vector<double>& orientedTopologyLinkVelocityMetersPerSecond,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureProjectionSettings&
        settings = {},
    const PlanarPressureRegionFragmentOpeningPressureProjectionLimits& limits =
        {});

} // namespace simwing::fsi::fluid
