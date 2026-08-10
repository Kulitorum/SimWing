#pragma once

#include "fluid/planar_region_fragment_opening_pressure_step.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningAcceptedStateVersion = 1;

struct PlanarPressureRegionFragmentOpeningAcceptedStateLimits {
    PlanarPressureRegionFragmentOpeningPressureOperatorLimits
        pressureOperatorLimits;
    PlanarPressureRegionFragmentVolumeRateLimits volumeRateLimits;
    PlanarPressureRegionFragmentOpeningFluxLimits openingFluxLimits;
    std::size_t maximumTopologyLinkVelocities = 140'000'000;
    std::size_t maximumOpeningSamples = 20'000'000;
    std::size_t maximumPressureCorrections = 20'000'000;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Immutable continuation endpoint for one accepted active-aperture pressure
// step. It owns the corrected same-region Cartesian link velocities, canonical
// material-relative opening samples, rebuilt opening-flux state, and connected-
// component-gauged correction pressure as one fingerprinted artifact. Source
// fingerprints bind the exact pressure graph, topology, volume-rate epoch,
// opening overlay, resistance coefficients, and pre-step flux.
//
// Capture independently reconstructs final continuity, pressure gauges, and
// diagonal kinetic energy before accepting the nested step evidence. This is
// aperture-flow continuation only: it does not compose total regional pressure,
// subtract open-area sheet loads, apply Structure traction, rebase topology, or
// select a production worker path.
struct PlanarPressureRegionFragmentOpeningAcceptedState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningAcceptedStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourcePressureOperatorFingerprint = 0;
    std::uint64_t sourceBasePressureOperatorFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    std::uint64_t sourceOpeningFluxFingerprint = 0;
    std::uint64_t resultOpeningFluxFingerprint = 0;
    std::uint64_t resistanceDefinitionFingerprint = 0;
    PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
    std::vector<double>
        orientedTopologyLinkVelocityMetersPerSecond;
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        openingVelocitySamples;
    PlanarPressureRegionFragmentOpeningFluxState openingFlux;
    std::vector<double> pressureCorrectionPascals;
    double correctedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double correctedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double continuityToleranceCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectionVolumeMeanPascals = 0.0;
    double maximumAbsolutePressureCorrectionPascals = 0.0;
    double maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        0.0;
    double momentumResidualToleranceKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double authoredPressureWorkJoules = 0.0;
    double geometryPressureWorkJoules = 0.0;
    double correctionKineticEnergyJoules = 0.0;
    double dissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;
    double energyToleranceJoules = 0.0;
    std::size_t pressureSolveIterationCount = 0;
    double pressureSolveFinalResidualL2PascalsMeters = 0.0;
    double pressureSolveFinalResidualMaximumPascalsMeters = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningAcceptedState&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningAcceptedState
capturePlanarPressureRegionFragmentOpeningAcceptedState(
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
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepDiagnostics&
        diagnostics,
    std::span<const double>
        orientedTopologyLinkVelocityMetersPerSecond,
    std::span<const PlanarPressureRegionFragmentOpeningVelocitySample>
        openingVelocitySamples,
    const PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    std::span<const double> pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state);

void validatePlanarPressureRegionFragmentOpeningAcceptedState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
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
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits = {});

} // namespace simwing::fsi::fluid
