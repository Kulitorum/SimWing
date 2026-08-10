#include "fluid/planar_region_fragment_pressure_projection.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

void validateSettings(
    const PlanarPressureRegionFragmentPressureProjectionSettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(
            settings.absoluteContinuityToleranceCubicMetersPerSecond)
        || settings.absoluteContinuityToleranceCubicMetersPerSecond < 0.0
        || !std::isfinite(settings.relativeContinuityTolerance)
        || settings.relativeContinuityTolerance < 0.0
        || (settings.absoluteContinuityToleranceCubicMetersPerSecond == 0.0
            && settings.relativeContinuityTolerance == 0.0)) {
        throw std::invalid_argument(
            "planar regional pressure-projection settings are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional pressure-projection storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional pressure-projection storage overflows");
    }
    return first + second;
}

std::size_t workingStorageBytes(const std::size_t fragmentCount,
                                const std::size_t linkCount,
                                const std::size_t componentCount) {
    return checkedSum(
        checkedSum(
            checkedProduct(9 * sizeof(double), fragmentCount),
            checkedProduct(sizeof(double), linkCount)),
        checkedProduct(
            sizeof(
                PlanarPressureRegionFragmentPressureSolveComponentDiagnostics),
            componentCount));
}

bool allFinite(const std::vector<double>& values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

bool isStaticGeometry(const PlanarPressureRegionSweepLedger& sweep) {
    if (sweep.previousProfile != sweep.currentProfile
        || sweep.globalGeometryVolumeChangeCubicMeters != 0.0
        || sweep.globalBoundarySweptVolumeCubicMeters != 0.0
        || sweep.globalSurfaceGeometryResidualCubicMeters != 0.0) {
        return false;
    }
    return std::ranges::all_of(
        sweep.intervals,
        [](const PlanarPressureRegionIntervalSweep& interval) {
            return interval.previousVolumeCubicMeters
                    == interval.currentVolumeCubicMeters
                && interval.geometryVolumeChangeCubicMeters == 0.0
                && interval.lowerSurfaceDisplacementMeters == 0.0
                && interval.upperSurfaceDisplacementMeters == 0.0
                && interval.lowerSurfaceVelocityMetersPerSecond == 0.0
                && interval.upperSurfaceVelocityMetersPerSecond == 0.0
                && interval.boundarySweptVolumeCubicMeters == 0.0
                && interval.surfaceGeometryResidualCubicMeters == 0.0;
        });
}

void accumulateNetOutwardFlow(
    const PlanarPressureRegionFragmentTopology& topology,
    const std::vector<double>& orientedNormalVelocityMetersPerSecond,
    std::vector<double>& netOutwardFlowCubicMetersPerSecond) {
    std::ranges::fill(netOutwardFlowCubicMetersPerSecond, 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        const double flow = link.areaSquareMeters
            * orientedNormalVelocityMetersPerSecond[link.linkIndex];
        netOutwardFlowCubicMetersPerSecond[link.minusFragmentIndex] += flow;
        netOutwardFlowCubicMetersPerSecond[link.plusFragmentIndex] -= flow;
    }
}

double rootMeanSquare(const std::vector<double>& values) {
    double squaredSum = 0.0;
    for (const double value : values) squaredSum += value * value;
    return std::sqrt(squaredSum / static_cast<double>(values.size()));
}

double maximumAbsolute(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double maximumComponentBalance(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const std::vector<double>& netOutwardFlowCubicMetersPerSecond) {
    double maximum = 0.0;
    for (const auto& component : pressureOperator.components) {
        double balance = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            balance += netOutwardFlowCubicMetersPerSecond[
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset]];
        }
        maximum = std::max(maximum, std::abs(balance));
    }
    return maximum;
}

} // namespace

PlanarPressureRegionFragmentPressureProjectionDiagnostics
projectStaticPlanarPressureRegionFragmentFaceVelocities(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::vector<double>& orientedNormalVelocityMetersPerSecond,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentPressureProjectionSettings& settings,
    const PlanarPressureRegionFragmentPressureProjectionLimits& limits) {
    validateSettings(settings);
    if (limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional pressure-projection limits are invalid");
    }
    validatePlanarPressureRegionFragmentPressureOperator(
        pressureOperator, grid, sweep, fragments, topology,
        limits.pressureOperatorLimits);
    if (orientedNormalVelocityMetersPerSecond.size()
            != topology.links.size()
        || pressureCorrectionPascals.size()
            != pressureOperator.rows.size()
        || !allFinite(orientedNormalVelocityMetersPerSecond)
        || !allFinite(pressureCorrectionPascals)) {
        throw std::invalid_argument(
            "planar regional pressure-projection fields are invalid");
    }

    PlanarPressureRegionFragmentPressureProjectionDiagnostics diagnostics;
    diagnostics.finite = true;
    diagnostics.pressureOperatorFingerprint = pressureOperator.fingerprint;
    diagnostics.topologyFingerprint = topology.fingerprint;
    diagnostics.fragmentFingerprint = fragments.fingerprint;
    diagnostics.fragmentCount = fragments.fragments.size();
    diagnostics.linkCount = topology.links.size();
    diagnostics.workingStorageBytes = workingStorageBytes(
        fragments.fragments.size(), topology.links.size(),
        pressureOperator.components.size());
    if (diagnostics.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional pressure-projection working storage exceeds its limit");
    }
    diagnostics.staticGeometry = isStaticGeometry(sweep);
    if (!diagnostics.staticGeometry) {
        throw std::invalid_argument(
            "planar regional pressure projection requires static geometry");
    }
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            ++diagnostics.sealedPressureLayerWallLinkCount;
            if (orientedNormalVelocityMetersPerSecond[link.linkIndex]
                != 0.0) {
                throw std::invalid_argument(
                    "planar regional pressure-layer wall velocity is nonzero");
            }
        } else {
            ++diagnostics.projectedSameRegionGridLinkCount;
        }
    }

    std::vector<double> predictedNetOutwardFlow(
        fragments.fragments.size(), 0.0);
    accumulateNetOutwardFlow(
        topology, orientedNormalVelocityMetersPerSecond,
        predictedNetOutwardFlow);
    diagnostics.predictedNetOutwardFlowL2CubicMetersPerSecond =
        rootMeanSquare(predictedNetOutwardFlow);
    diagnostics.predictedNetOutwardFlowMaximumCubicMetersPerSecond =
        maximumAbsolute(predictedNetOutwardFlow);
    diagnostics.maximumAbsolutePredictedComponentBalanceCubicMetersPerSecond =
        maximumComponentBalance(pressureOperator, predictedNetOutwardFlow);

    std::vector<double> integratedRightHandSide(
        pressureOperator.rows.size(), 0.0);
    const double rightHandSideScale =
        -settings.densityKgPerCubicMeter / settings.timeStepSeconds;
    for (std::size_t row = 0; row < integratedRightHandSide.size(); ++row) {
        integratedRightHandSide[row] =
            rightHandSideScale * predictedNetOutwardFlow[row];
    }
    if (!allFinite(integratedRightHandSide)) {
        diagnostics.finite = false;
        return diagnostics;
    }

    std::vector<double> candidatePressure = pressureCorrectionPascals;
    diagnostics.pressureSolve =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, grid, sweep, fragments, topology,
            integratedRightHandSide, candidatePressure,
            settings.pressureSolve);
    diagnostics.finite = diagnostics.pressureSolve.finite;
    if (!diagnostics.pressureSolve.compatible
        || !diagnostics.pressureSolve.converged
        || !diagnostics.pressureSolve.finite) {
        return diagnostics;
    }

    std::vector<double> candidateVelocity =
        orientedNormalVelocityMetersPerSecond;
    const double velocityCorrectionScale =
        settings.timeStepSeconds / settings.densityKgPerCubicMeter;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        const double correction = velocityCorrectionScale
            * (candidatePressure[link.minusFragmentIndex]
               - candidatePressure[link.plusFragmentIndex])
            / link.centerDistanceMeters;
        candidateVelocity[link.linkIndex] += correction;
        diagnostics.maximumAbsoluteVelocityCorrectionMetersPerSecond =
            std::max(
                diagnostics.maximumAbsoluteVelocityCorrectionMetersPerSecond,
                std::abs(correction));
    }
    if (!allFinite(candidateVelocity)) {
        diagnostics.finite = false;
        return diagnostics;
    }

    std::vector<double> correctedNetOutwardFlow(
        fragments.fragments.size(), 0.0);
    accumulateNetOutwardFlow(
        topology, candidateVelocity, correctedNetOutwardFlow);
    diagnostics.correctedNetOutwardFlowL2CubicMetersPerSecond =
        rootMeanSquare(correctedNetOutwardFlow);
    diagnostics.correctedNetOutwardFlowMaximumCubicMetersPerSecond =
        maximumAbsolute(correctedNetOutwardFlow);
    diagnostics.maximumAbsoluteCorrectedComponentBalanceCubicMetersPerSecond =
        maximumComponentBalance(pressureOperator, correctedNetOutwardFlow);
    diagnostics.continuityToleranceCubicMetersPerSecond = std::max(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance
            * diagnostics
                .predictedNetOutwardFlowMaximumCubicMetersPerSecond);
    diagnostics.finite = std::isfinite(
            diagnostics.correctedNetOutwardFlowL2CubicMetersPerSecond)
        && std::isfinite(
            diagnostics.correctedNetOutwardFlowMaximumCubicMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsoluteVelocityCorrectionMetersPerSecond)
        && std::isfinite(
            diagnostics.continuityToleranceCubicMetersPerSecond);
    diagnostics.accepted = diagnostics.finite
        && diagnostics.correctedNetOutwardFlowMaximumCubicMetersPerSecond
            <= diagnostics.continuityToleranceCubicMetersPerSecond;
    if (!diagnostics.accepted) return diagnostics;

    orientedNormalVelocityMetersPerSecond = std::move(candidateVelocity);
    pressureCorrectionPascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
