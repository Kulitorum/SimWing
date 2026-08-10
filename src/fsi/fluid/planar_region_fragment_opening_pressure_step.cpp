#include "fluid/planar_region_fragment_opening_pressure_step.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error("opening pressure-step storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error("opening pressure-step storage overflows");
    }
    return first + second;
}

double roundoffTolerance(const std::initializer_list<double> values) {
    double scale = 1.0;
    for (const double value : values)
        scale = std::max(scale, std::abs(value));
    return 2048.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureStepDiagnostics
advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    std::vector<double>& orientedTopologyLinkVelocityMetersPerSecond,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureStepLimits& limits) {
    if (limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening pressure-step limits are invalid");
    }
    PlanarPressureRegionFragmentOpeningPressureStepDiagnostics diagnostics;
    diagnostics.sourceOpeningFluxFingerprint = openingFluxState.fingerprint;
    diagnostics.workingStorageBytes = checkedAdd(
        checkedAdd(
            checkedMultiply(
                orientedTopologyLinkVelocityMetersPerSecond.size(),
                sizeof(double)),
            checkedMultiply(
                pressureCorrectionPascals.size(), sizeof(double))),
        checkedAdd(
            checkedMultiply(
                openingVelocitySamples.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocitySample)),
            openingFluxState.ownedStorageBytes));
    if (diagnostics.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening pressure-step working storage exceeds its limit");
    }

    auto candidateTopologyVelocity =
        orientedTopologyLinkVelocityMetersPerSecond;
    auto candidateOpeningSamples = openingVelocitySamples;
    auto candidateOpeningFlux = openingFluxState;
    auto candidatePressure = pressureCorrectionPascals;
    PlanarPressureRegionFragmentOpeningResistanceSettings resistanceSettings;
    resistanceSettings.densityKgPerCubicMeter =
        settings.projection.densityKgPerCubicMeter;
    resistanceSettings.timeStepSeconds =
        settings.projection.timeStepSeconds;
    resistanceSettings.useAuthoredPressureDrive =
        settings.useAuthoredPressureDrive;
    diagnostics.resistance =
        advancePlanarPressureRegionFragmentOpeningResistance(
            grid, sweep, fragments, topology, openingDefinitions, openings,
            resistanceDefinitions, candidateOpeningSamples,
            candidateOpeningFlux, resistanceSettings,
            limits.resistanceLimits);
    diagnostics.finite = diagnostics.resistance.finite;
    if (!diagnostics.resistance.accepted
        || !diagnostics.resistance.finite) {
        return diagnostics;
    }

    diagnostics.projection =
        projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, volumeRates, openingDefinitions, openings,
            candidateTopologyVelocity, candidateOpeningSamples,
            candidateOpeningFlux, candidatePressure, settings.projection,
            limits.projectionLimits);
    diagnostics.finite = diagnostics.projection.finite;
    if (!diagnostics.projection.accepted
        || !diagnostics.projection.finite) {
        return diagnostics;
    }

    diagnostics.resultOpeningFluxFingerprint =
        candidateOpeningFlux.fingerprint;
    diagnostics.kineticEnergyBeforeJoules =
        diagnostics.projection.kineticEnergyBeforeJoules
        + diagnostics.resistance.kineticEnergyBeforeJoules
        - diagnostics.resistance.kineticEnergyAfterJoules;
    diagnostics.kineticEnergyAfterJoules =
        diagnostics.projection.kineticEnergyAfterJoules;
    diagnostics.kineticEnergyChangeJoules =
        diagnostics.kineticEnergyAfterJoules
        - diagnostics.kineticEnergyBeforeJoules;
    diagnostics.geometryPressureWorkJoules =
        diagnostics.projection.geometryPressureWorkJoules;
    diagnostics.authoredPressureWorkJoules =
        diagnostics.resistance.authoredPressureWorkJoules;
    diagnostics.correctionKineticEnergyJoules =
        diagnostics.projection.correctionKineticEnergyJoules;
    diagnostics.dissipatedEnergyJoules =
        diagnostics.resistance.dissipatedEnergyJoules;
    diagnostics.energyResidualJoules =
        diagnostics.kineticEnergyChangeJoules
        - diagnostics.authoredPressureWorkJoules
        - diagnostics.geometryPressureWorkJoules
        + diagnostics.correctionKineticEnergyJoules
        + diagnostics.dissipatedEnergyJoules;
    diagnostics.energyToleranceJoules = std::max({
        diagnostics.resistance.energyToleranceJoules,
        diagnostics.projection.energyResidualToleranceJoules,
        roundoffTolerance({
            diagnostics.kineticEnergyBeforeJoules,
            diagnostics.kineticEnergyAfterJoules,
            diagnostics.kineticEnergyChangeJoules,
            diagnostics.authoredPressureWorkJoules,
            diagnostics.geometryPressureWorkJoules,
            diagnostics.correctionKineticEnergyJoules,
            diagnostics.dissipatedEnergyJoules,
        }),
    });
    diagnostics.finite = std::ranges::all_of(
        std::initializer_list<double>{
            diagnostics.kineticEnergyBeforeJoules,
            diagnostics.kineticEnergyAfterJoules,
            diagnostics.kineticEnergyChangeJoules,
            diagnostics.authoredPressureWorkJoules,
            diagnostics.geometryPressureWorkJoules,
            diagnostics.correctionKineticEnergyJoules,
            diagnostics.dissipatedEnergyJoules,
            diagnostics.energyResidualJoules,
            diagnostics.energyToleranceJoules,
        }, [](const double value) { return std::isfinite(value); });
    diagnostics.energyAccepted = diagnostics.finite
        && diagnostics.dissipatedEnergyJoules
            >= -diagnostics.energyToleranceJoules
        && std::abs(diagnostics.energyResidualJoules)
            <= diagnostics.energyToleranceJoules;
    diagnostics.accepted = diagnostics.energyAccepted;
    if (!diagnostics.accepted) return diagnostics;

    orientedTopologyLinkVelocityMetersPerSecond =
        std::move(candidateTopologyVelocity);
    openingVelocitySamples = std::move(candidateOpeningSamples);
    openingFluxState = std::move(candidateOpeningFlux);
    pressureCorrectionPascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
