#include "fluid/planar_region_fragment_opening_pressure_projection.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

void validateTolerancePair(const double absolute,
                           const double relative,
                           const char* message) {
    if (!std::isfinite(absolute) || absolute < 0.0
        || !std::isfinite(relative) || relative < 0.0
        || (absolute == 0.0 && relative == 0.0)) {
        throw std::invalid_argument(message);
    }
}

void validateSettings(
    const PlanarPressureRegionFragmentOpeningPressureProjectionSettings&
        settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "opening pressure-projection physical settings are invalid");
    }
    validateTolerancePair(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance,
        "opening pressure-projection continuity tolerances are invalid");
    validateTolerancePair(
        settings.absoluteMomentumResidualToleranceKilogramMetersPerSecond,
        settings.relativeMomentumResidualTolerance,
        "opening pressure-projection momentum tolerances are invalid");
    validateTolerancePair(
        settings.absoluteEnergyResidualToleranceJoules,
        settings.relativeEnergyResidualTolerance,
        "opening pressure-projection energy tolerances are invalid");
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening pressure-projection storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening pressure-projection storage overflows");
    }
    return first + second;
}

std::size_t workingStorageBytes(const std::size_t fragmentCount,
                                const std::size_t topologyLinkCount,
                                const std::size_t openingPatchCount) {
    return checkedAdd(
        checkedAdd(
            checkedMultiply(10 * sizeof(double), fragmentCount),
            checkedMultiply(sizeof(double), topologyLinkCount)),
        checkedMultiply(
            2 * sizeof(
                PlanarPressureRegionFragmentOpeningVelocitySample)
                + sizeof(std::size_t),
            openingPatchCount));
}

bool allFinite(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

double rootMeanSquare(const std::span<const double> values) {
    double squared = 0.0;
    for (const double value : values) squared += value * value;
    return std::sqrt(squared / static_cast<double>(values.size()));
}

double maximumAbsolute(const std::span<const double> values) {
    double maximum = 0.0;
    for (const double value : values)
        maximum = std::max(maximum, std::abs(value));
    return maximum;
}

double scaledTolerance(const double absolute,
                       const double relative,
                       const std::initializer_list<double> scales) {
    double scale = 0.0;
    for (const double value : scales)
        scale = std::max(scale, std::abs(value));
    return std::max(absolute, relative * scale);
}

void accumulateGridFlow(
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<const double> velocities,
    std::vector<double>& flow) {
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        const double volumeFlow =
            link.areaSquareMeters * velocities[link.linkIndex];
        flow[link.minusFragmentIndex] += volumeFlow;
        flow[link.plusFragmentIndex] -= volumeFlow;
    }
}

void accumulateOpeningFlow(
    const PlanarPressureRegionFragmentOpeningFluxState& openingFlux,
    std::vector<double>& flow) {
    for (const auto& fragment : openingFlux.fragments) {
        flow[fragment.fragmentIndex] +=
            fragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond;
    }
}

void accumulateGeometryRate(
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::vector<double>& flow) {
    for (const auto& fragment : volumeRates.fragments) {
        flow[fragment.fragmentIndex] +=
            fragment.geometryVolumeChangeRateCubicMetersPerSecond;
    }
}

double maximumComponentResidual(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const std::span<const double> residual) {
    double maximum = 0.0;
    for (const auto& component : pressureOperator.components) {
        double sum = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            sum += residual[pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset]];
        }
        maximum = std::max(maximum, std::abs(sum));
    }
    return maximum;
}

void accumulateDegreeEnergy(
    const double density,
    const double timeStep,
    const double area,
    const double distance,
    const double pressureDifference,
    const double before,
    const double after,
    PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics&
        diagnostics,
    double& maximumMomentumScale) {
    const double mass = density * area * distance;
    const double velocityChange = after - before;
    const double momentumChange = mass * velocityChange;
    const double pressureImpulse =
        timeStep * area * pressureDifference;
    const double momentumResidual = momentumChange - pressureImpulse;
    const double energyBefore = 0.5 * mass * before * before;
    const double energyAfter = 0.5 * mass * after * after;
    const double energyChange = energyAfter - energyBefore;
    const double pressureWork =
        pressureImpulse * 0.5 * (before + after);
    const double energyResidual = energyChange - pressureWork;
    const double correctionEnergy =
        0.5 * mass * velocityChange * velocityChange;
    diagnostics.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        std::max(
            diagnostics
                .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
            std::abs(momentumResidual));
    maximumMomentumScale = std::max(
        maximumMomentumScale,
        std::max(std::abs(momentumChange), std::abs(pressureImpulse)));
    diagnostics.kineticEnergyBeforeJoules += energyBefore;
    diagnostics.kineticEnergyAfterJoules += energyAfter;
    diagnostics.kineticEnergyChangeJoules += energyChange;
    diagnostics.midpointPressureWorkJoules += pressureWork;
    diagnostics.correctionKineticEnergyJoules += correctionEnergy;
    diagnostics.workEnergyResidualJoules += energyResidual;
    diagnostics.maximumAbsoluteWorkEnergyResidualJoules = std::max(
        diagnostics.maximumAbsoluteWorkEnergyResidualJoules,
        std::abs(energyResidual));
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics
projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
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
    std::vector<double>& orientedTopologyLinkVelocityMetersPerSecond,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureProjectionSettings&
        settings,
    const PlanarPressureRegionFragmentOpeningPressureProjectionLimits& limits) {
    validateSettings(settings);
    if (limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening pressure-projection limits are invalid");
    }
    validatePlanarPressureRegionFragmentOpeningPressureOperator(
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, openingDefinitions, openings,
        limits.pressureOperatorLimits);
    validatePlanarPressureRegionFragmentVolumeRates(
        volumeRates, grid, sweep, fragments, topology,
        limits.volumeRateLimits);
    if (volumeRates.durationSeconds != settings.timeStepSeconds) {
        throw std::invalid_argument(
            "opening pressure-projection duration does not match volume rates");
    }
    validatePlanarPressureRegionFragmentOpeningFluxState(
        openingFluxState, grid, sweep, fragments, topology,
        openingDefinitions, openings, openingVelocitySamples,
        limits.openingFluxLimits);
    if (orientedTopologyLinkVelocityMetersPerSecond.size()
            != topology.links.size()
        || pressureCorrectionPascals.size() != pressureOperator.rows.size()
        || !allFinite(orientedTopologyLinkVelocityMetersPerSecond)
        || !allFinite(pressureCorrectionPascals)) {
        throw std::invalid_argument(
            "opening pressure-projection fields are invalid");
    }
    for (const auto& link : topology.links) {
        if (link.kind
                == PlanarPressureRegionFragmentFaceKind::PressureLayerWall
            && orientedTopologyLinkVelocityMetersPerSecond[link.linkIndex]
                != 0.0) {
            throw std::invalid_argument(
                "opening pressure-projection solid-wall velocity is nonzero");
        }
    }

    PlanarPressureRegionFragmentOpeningPressureProjectionDiagnostics
        diagnostics;
    diagnostics.finite = true;
    diagnostics.pressureOperatorFingerprint = pressureOperator.fingerprint;
    diagnostics.basePressureOperatorFingerprint =
        basePressureOperator.fingerprint;
    diagnostics.topologyFingerprint = topology.fingerprint;
    diagnostics.fragmentFingerprint = fragments.fingerprint;
    diagnostics.volumeRateFingerprint = volumeRates.fingerprint;
    diagnostics.openingFingerprint = openings.fingerprint;
    diagnostics.predictedOpeningFluxFingerprint = openingFluxState.fingerprint;
    diagnostics.fragmentCount = fragments.fragments.size();
    diagnostics.topologyLinkCount = topology.links.size();
    diagnostics.openingPatchCount = openings.patches.size();
    diagnostics.workingStorageBytes = workingStorageBytes(
        fragments.fragments.size(), topology.links.size(),
        openings.patches.size());
    if (diagnostics.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening pressure-projection working storage exceeds its limit");
    }
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            ++diagnostics.correctedSameRegionGridLinkCount;
        } else {
            ++diagnostics.retainedPressureLayerWallLinkCount;
        }
    }

    std::vector<double> predictedResidual(
        fragments.fragments.size(), 0.0);
    accumulateGridFlow(
        topology, orientedTopologyLinkVelocityMetersPerSecond,
        predictedResidual);
    accumulateOpeningFlow(openingFluxState, predictedResidual);
    accumulateGeometryRate(volumeRates, predictedResidual);
    diagnostics.predictedContinuityResidualL2CubicMetersPerSecond =
        rootMeanSquare(predictedResidual);
    diagnostics.predictedContinuityResidualMaximumCubicMetersPerSecond =
        maximumAbsolute(predictedResidual);
    diagnostics
        .maximumAbsolutePredictedComponentContinuityResidualCubicMetersPerSecond =
        maximumComponentResidual(pressureOperator, predictedResidual);

    std::vector<double> integratedRightHandSide = predictedResidual;
    const double rightHandSideScale =
        -settings.densityKgPerCubicMeter / settings.timeStepSeconds;
    for (double& value : integratedRightHandSide)
        value *= rightHandSideScale;
    if (!allFinite(integratedRightHandSide)) {
        diagnostics.finite = false;
        return diagnostics;
    }

    std::vector<double> candidatePressure = pressureCorrectionPascals;
    diagnostics.pressureSolve =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, openingDefinitions, openings, integratedRightHandSide,
            candidatePressure, settings.pressureSolve);
    diagnostics.finite = diagnostics.pressureSolve.finite;
    if (!diagnostics.pressureSolve.compatible
        || !diagnostics.pressureSolve.converged
        || !diagnostics.pressureSolve.finite) {
        return diagnostics;
    }

    std::vector<double> candidateTopologyVelocity =
        orientedTopologyLinkVelocityMetersPerSecond;
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        candidateOpeningSamples = openingVelocitySamples;
    std::vector<std::size_t> sampleOrder(candidateOpeningSamples.size());
    std::iota(sampleOrder.begin(), sampleOrder.end(), 0);
    std::ranges::sort(sampleOrder, [&](const std::size_t first,
                                      const std::size_t second) {
        return candidateOpeningSamples[first].patchStableId
            < candidateOpeningSamples[second].patchStableId;
    });
    const double velocityCorrectionScale =
        settings.timeStepSeconds / settings.densityKgPerCubicMeter;
    double maximumMomentumScale = 0.0;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        const double pressureDifference =
            candidatePressure[link.minusFragmentIndex]
            - candidatePressure[link.plusFragmentIndex];
        const double correction = velocityCorrectionScale
            * pressureDifference / link.centerDistanceMeters;
        const double before = candidateTopologyVelocity[link.linkIndex];
        const double after = before + correction;
        if (!std::isfinite(after)) {
            diagnostics.finite = false;
            return diagnostics;
        }
        candidateTopologyVelocity[link.linkIndex] = after;
        diagnostics.maximumAbsoluteGridVelocityCorrectionMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteGridVelocityCorrectionMetersPerSecond,
                std::abs(correction));
        accumulateDegreeEnergy(
            settings.densityKgPerCubicMeter, settings.timeStepSeconds,
            link.areaSquareMeters, link.centerDistanceMeters,
            pressureDifference, before, after, diagnostics,
            maximumMomentumScale);
    }
    for (const auto& patch : openings.patches) {
        const auto ordered = std::ranges::lower_bound(
            sampleOrder, patch.patchStableId, {},
            [&](const std::size_t index) {
                return candidateOpeningSamples[index].patchStableId;
            });
        if (ordered == sampleOrder.end()
            || candidateOpeningSamples[*ordered].patchStableId
                != patch.patchStableId) {
            throw std::logic_error(
                "opening pressure-projection sample mapping is incomplete");
        }
        auto& sample = candidateOpeningSamples[*ordered];
        const double pressureDifference =
            candidatePressure[patch.minusFragmentIndex]
            - candidatePressure[patch.plusFragmentIndex];
        const double correction = velocityCorrectionScale
            * pressureDifference / patch.centerDistanceMeters;
        const double before = sample.relativeNormalVelocityMetersPerSecond;
        const double after = before + correction;
        if (!std::isfinite(after)) {
            diagnostics.finite = false;
            return diagnostics;
        }
        sample.relativeNormalVelocityMetersPerSecond = after;
        diagnostics.maximumAbsoluteOpeningVelocityCorrectionMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteOpeningVelocityCorrectionMetersPerSecond,
                std::abs(correction));
        accumulateDegreeEnergy(
            settings.densityKgPerCubicMeter, settings.timeStepSeconds,
            patch.areaSquareMeters, patch.centerDistanceMeters,
            pressureDifference, before, after, diagnostics,
            maximumMomentumScale);
    }

    auto candidateOpeningFlux =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            grid, sweep, fragments, topology, openingDefinitions, openings,
            candidateOpeningSamples, limits.openingFluxLimits);
    diagnostics.correctedOpeningFluxFingerprint =
        candidateOpeningFlux.fingerprint;
    std::vector<double> correctedResidual(
        fragments.fragments.size(), 0.0);
    accumulateGridFlow(topology, candidateTopologyVelocity, correctedResidual);
    accumulateOpeningFlow(candidateOpeningFlux, correctedResidual);
    accumulateGeometryRate(volumeRates, correctedResidual);
    diagnostics.correctedContinuityResidualL2CubicMetersPerSecond =
        rootMeanSquare(correctedResidual);
    diagnostics.correctedContinuityResidualMaximumCubicMetersPerSecond =
        maximumAbsolute(correctedResidual);
    diagnostics
        .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond =
        maximumComponentResidual(pressureOperator, correctedResidual);
    diagnostics.continuityToleranceCubicMetersPerSecond = scaledTolerance(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance,
        {diagnostics.predictedContinuityResidualMaximumCubicMetersPerSecond});

    diagnostics.kineticEnergyChangeJoules =
        diagnostics.kineticEnergyAfterJoules
        - diagnostics.kineticEnergyBeforeJoules;
    diagnostics.workEnergyResidualJoules =
        diagnostics.kineticEnergyChangeJoules
        - diagnostics.midpointPressureWorkJoules;
    for (std::size_t index = 0; index < candidatePressure.size(); ++index) {
        diagnostics.geometryPressureWorkJoules -=
            settings.timeStepSeconds * candidatePressure[index]
            * volumeRates.fragments[index]
                  .geometryVolumeChangeRateCubicMetersPerSecond;
    }
    diagnostics.affineEnergyResidualJoules =
        diagnostics.kineticEnergyChangeJoules
        - diagnostics.geometryPressureWorkJoules
        + diagnostics.correctionKineticEnergyJoules;
    diagnostics.momentumResidualToleranceKilogramMetersPerSecond =
        scaledTolerance(
            settings.absoluteMomentumResidualToleranceKilogramMetersPerSecond,
            settings.relativeMomentumResidualTolerance,
            {maximumMomentumScale});
    diagnostics.energyResidualToleranceJoules = scaledTolerance(
        settings.absoluteEnergyResidualToleranceJoules,
        settings.relativeEnergyResidualTolerance,
        {diagnostics.kineticEnergyBeforeJoules,
         diagnostics.kineticEnergyAfterJoules,
         diagnostics.kineticEnergyChangeJoules,
         diagnostics.midpointPressureWorkJoules,
         diagnostics.correctionKineticEnergyJoules,
         diagnostics.geometryPressureWorkJoules});
    diagnostics.finite = allFinite(candidateTopologyVelocity)
        && allFinite(candidatePressure)
        && std::isfinite(
            diagnostics.correctedContinuityResidualL2CubicMetersPerSecond)
        && std::isfinite(
            diagnostics.correctedContinuityResidualMaximumCubicMetersPerSecond)
        && std::isfinite(
            diagnostics
                .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond)
        && std::isfinite(
            diagnostics
                .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond)
        && std::isfinite(diagnostics.workEnergyResidualJoules)
        && std::isfinite(
            diagnostics.maximumAbsoluteWorkEnergyResidualJoules)
        && std::isfinite(diagnostics.affineEnergyResidualJoules);
    diagnostics.energyAccepted = diagnostics.finite
        && diagnostics
                .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond
            <= diagnostics
                .momentumResidualToleranceKilogramMetersPerSecond
        && std::abs(diagnostics.workEnergyResidualJoules)
            <= diagnostics.energyResidualToleranceJoules
        && diagnostics.maximumAbsoluteWorkEnergyResidualJoules
            <= diagnostics.energyResidualToleranceJoules
        && std::abs(diagnostics.affineEnergyResidualJoules)
            <= diagnostics.energyResidualToleranceJoules;
    diagnostics.accepted = diagnostics.energyAccepted
        && diagnostics.correctedContinuityResidualMaximumCubicMetersPerSecond
            <= diagnostics.continuityToleranceCubicMetersPerSecond
        && diagnostics
               .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
            <= diagnostics.continuityToleranceCubicMetersPerSecond;
    if (!diagnostics.accepted) return diagnostics;

    orientedTopologyLinkVelocityMetersPerSecond =
        std::move(candidateTopologyVelocity);
    openingVelocitySamples = std::move(candidateOpeningSamples);
    openingFluxState = std::move(candidateOpeningFlux);
    pressureCorrectionPascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
