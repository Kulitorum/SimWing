#include "fluid/planar_region_fragment_opening_momentum_prediction.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening momentum-prediction storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening momentum-prediction storage overflows");
    }
    return first * second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits&
        limits) {
    if (limits.maximumFragments == 0 || limits.maximumDofs == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening momentum-prediction limits are invalid");
    }
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 add(const Vector3& first, const Vector3& second) {
    return {first.x + second.x, first.y + second.y,
            first.z + second.z};
}

Vector3 subtract(const Vector3& first, const Vector3& second) {
    return {first.x - second.x, first.y - second.y,
            first.z - second.z};
}

double norm(const Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

double component(const Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "opening momentum-prediction axis is invalid");
}

bool equalWithinRoundoff(const double first, const double second) {
    const double scaleValue = std::max(
        {1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second)
        <= 128.0 * std::numeric_limits<double>::epsilon() * scaleValue;
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t predictionFingerprint(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction) {
    Fingerprint fingerprint;
    fingerprint.integer(prediction.version);
    fingerprint.integer(prediction.sourceTransportFingerprint);
    fingerprint.integer(prediction.sourceTransportMetricFingerprint);
    fingerprint.integer(prediction.currentMetricFingerprint);
    fingerprint.integer(prediction.currentVolumeRateFingerprint);
    fingerprint.real(prediction.densityKgPerCubicMeter);
    fingerprint.real(prediction.timeStepSeconds);
    const auto& diagnostics = prediction.diagnostics;
    for (const std::size_t value : {
             diagnostics.fragmentCount,
             diagnostics.dofCount,
             diagnostics.sharedRegionGridDofCount,
             diagnostics.solidWallTraceDofCount,
             diagnostics.openingDofCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprintVector(
        fingerprint, diagnostics.sourceMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, diagnostics.remappedMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint,
        diagnostics.geometricMomentumChangeKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint,
        diagnostics.predictedStateMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint,
        diagnostics.reconstructionMomentumChangeKilogramMetersPerSecond);
    for (const double value : {
             diagnostics.sourceKineticEnergyJoules,
             diagnostics.remappedKineticEnergyJoules,
             diagnostics.predictedStateKineticEnergyJoules,
             diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
             diagnostics.maximumSourceSpeedMetersPerSecond,
             diagnostics.maximumEndpointNormalVelocityJumpMetersPerSecond,
             diagnostics.maximumAbsolutePredictedRelativeVelocityMetersPerSecond,
             diagnostics.maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.finite ? 1U : 0U));
    fingerprint.integer(prediction.predictedVelocityState.fingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.workingStorageBytes));
    return fingerprint.value();
}

std::size_t expectedStateOwnedBytes(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric) {
    return checkedAdd(
        checkedMultiply(
            metric.dofs.size(),
            sizeof(
                PlanarPressureRegionFragmentOpeningVelocityStateSample)),
        checkedAdd(
            checkedMultiply(
                metric.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateFragment)),
            checkedMultiply(
                metric.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateComponent))));
}

void validateDegreeIdentity(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric) {
    if (sourceMetric.dofs.size() != currentMetric.dofs.size()
        || sourceMetric.components.size() != currentMetric.components.size()
        || sourceMetric.profileAxis != currentMetric.profileAxis) {
        throw std::invalid_argument(
            "opening momentum-prediction degree topology changed");
    }
    for (std::size_t index = 0;
         index < sourceMetric.components.size(); ++index) {
        const auto& source = sourceMetric.components[index];
        const auto& current = currentMetric.components[index];
        if (source.componentIndex != index
            || current.componentIndex != index
            || source.stableId == 0
            || source.stableId != current.stableId
            || source.baseComponentCount != current.baseComponentCount) {
            throw std::invalid_argument(
                "opening momentum-prediction component topology changed");
        }
    }
    for (std::size_t index = 0; index < sourceMetric.dofs.size(); ++index) {
        const auto& source = sourceMetric.dofs[index];
        const auto& current = currentMetric.dofs[index];
        if (source.dofIndex != index || current.dofIndex != index
            || source.stableId == 0
            || source.stableId != current.stableId
            || source.kind != current.kind
            || source.sourceFaceLinkStableId
                != current.sourceFaceLinkStableId
            || source.sourceOpeningPatchStableId
                != current.sourceOpeningPatchStableId
            || source.axis != current.axis
            || source.surfaceStableId != current.surfaceStableId
            || source.ownerFragmentStableId
                != current.ownerFragmentStableId
            || source.oppositeFragmentStableId
                != current.oppositeFragmentStableId
            || source.connectedComponentIndex
                != current.connectedComponentIndex) {
            throw std::invalid_argument(
                "opening momentum-prediction degree identity changed");
        }
    }
}

PlanarPressureRegionFragmentOpeningMomentumPrediction buildPrediction(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits&
        limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        transport);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        transportTargetMetric);
    validatePlanarPressureRegionFragmentVolumeRates(
        currentVolumeRates, grid, currentSweep, currentFragments,
        currentTopology, limits.volumeRateLimits);
    validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        currentMetric, grid, currentSweep, currentFragments,
        currentTopology, currentBaseMetric, currentOpeningDefinitions,
        currentOpenings, limits.stateLimits.metricLimits);
    if (!transport.diagnostics.accepted
        || transport.targetMetricFingerprint
            != transportTargetMetric.fingerprint) {
        throw std::invalid_argument(
            "opening momentum-prediction source transport is invalid");
    }
    validateDegreeIdentity(transportTargetMetric, currentMetric);

    const std::size_t fragmentCount = currentMetric.fragments.size();
    const std::size_t dofCount = currentMetric.dofs.size();
    if (fragmentCount == 0
        || transport.controls.size() != fragmentCount
        || transportTargetMetric.fragments.size() != fragmentCount
        || currentVolumeRates.fragments.size() != fragmentCount
        || currentMetric.sourceFragmentFingerprint
            != currentFragments.fingerprint
        || currentMetric.sourceTopologyFingerprint
            != currentTopology.fingerprint
        || transport.densityKgPerCubicMeter <= 0.0) {
        throw std::invalid_argument(
            "opening momentum-prediction endpoint identity is invalid");
    }
    if (fragmentCount > limits.maximumFragments
        || dofCount > limits.maximumDofs) {
        throw std::length_error(
            "opening momentum-prediction entity limit exceeded");
    }
    const std::size_t expectedOwnedBytes = expectedStateOwnedBytes(
        currentMetric);
    const std::size_t expectedStateWorkingBytes = checkedMultiply(
        checkedMultiply(dofCount, 3), sizeof(double));
    std::size_t expectedWorkingBytes = checkedMultiply(
        checkedMultiply(fragmentCount, 3), sizeof(std::size_t));
    expectedWorkingBytes = checkedAdd(
        expectedWorkingBytes, expectedStateWorkingBytes);
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || expectedWorkingBytes > limits.maximumWorkingBytes
        || dofCount > limits.stateLimits.maximumSamples
        || fragmentCount > limits.stateLimits.maximumFragments
        || currentMetric.components.size()
            > limits.stateLimits.maximumComponents
        || expectedOwnedBytes > limits.stateLimits.maximumOwnedBytes
        || expectedStateWorkingBytes
            > limits.stateLimits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum-prediction storage limit exceeded");
    }

    std::vector<std::size_t> sourceOrder(fragmentCount);
    std::vector<std::size_t> currentOrder(fragmentCount);
    std::iota(sourceOrder.begin(), sourceOrder.end(), 0);
    std::iota(currentOrder.begin(), currentOrder.end(), 0);
    std::ranges::sort(sourceOrder, {}, [&](const std::size_t index) {
        return transport.controls[index].stableId;
    });
    std::ranges::sort(currentOrder, {}, [&](const std::size_t index) {
        return currentMetric.fragments[index].stableId;
    });
    std::vector<std::size_t> sourceByCurrent(fragmentCount);

    PlanarPressureRegionFragmentOpeningMomentumPrediction result;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.sourceTransportMetricFingerprint =
        transportTargetMetric.fingerprint;
    result.currentMetricFingerprint = currentMetric.fingerprint;
    result.currentVolumeRateFingerprint = currentVolumeRates.fingerprint;
    result.densityKgPerCubicMeter = transport.densityKgPerCubicMeter;
    result.timeStepSeconds = currentVolumeRates.durationSeconds;
    result.workingStorageBytes = expectedWorkingBytes;
    auto& diagnostics = result.diagnostics;
    diagnostics.fragmentCount = fragmentCount;
    diagnostics.dofCount = dofCount;
    diagnostics.sourceMomentumKilogramMetersPerSecond =
        transport.diagnostics.momentumAfterKilogramMetersPerSecond;
    diagnostics.sourceKineticEnergyJoules =
        transport.diagnostics.kineticEnergyAfterJoules;

    for (std::size_t offset = 0; offset < fragmentCount; ++offset) {
        const std::size_t sourceIndex = sourceOrder[offset];
        const std::size_t currentIndex = currentOrder[offset];
        const auto& source = transport.controls[sourceIndex];
        const auto& current = currentMetric.fragments[currentIndex];
        const auto& rate = currentVolumeRates.fragments[currentIndex];
        if (source.stableId == 0 || source.stableId != current.stableId
            || source.regionStableId != current.regionStableId
            || source.connectedComponentIndex
                != current.connectedComponentIndex
            || rate.fragmentIndex != currentIndex
            || rate.stableId != current.stableId
            || rate.regionStableId != current.regionStableId
            || !equalWithinRoundoff(
                source.volumeCubicMeters,
                rate.previousVolumeCubicMeters)
            || !equalWithinRoundoff(
                current.sourceVolumeCubicMeters,
                rate.currentVolumeCubicMeters)) {
            throw std::invalid_argument(
                "opening momentum-prediction fragment mapping is invalid");
        }
        sourceByCurrent[currentIndex] = sourceIndex;
        diagnostics.maximumAbsoluteVolumeChangeCubicMeters = std::max(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
            std::abs(
                current.sourceVolumeCubicMeters
                - source.volumeCubicMeters));
        diagnostics.maximumSourceSpeedMetersPerSecond = std::max(
            diagnostics.maximumSourceSpeedMetersPerSecond,
            norm(source.velocityMetersPerSecond));
        const double currentMass = result.densityKgPerCubicMeter
            * current.sourceVolumeCubicMeters;
        diagnostics.remappedMomentumKilogramMetersPerSecond = add(
            diagnostics.remappedMomentumKilogramMetersPerSecond,
            {
                currentMass * source.velocityMetersPerSecond.x,
                currentMass * source.velocityMetersPerSecond.y,
                currentMass * source.velocityMetersPerSecond.z,
            });
        const double speed = norm(source.velocityMetersPerSecond);
        diagnostics.remappedKineticEnergyJoules +=
            0.5 * currentMass * speed * speed;
    }
    diagnostics.geometricMomentumChangeKilogramMetersPerSecond = subtract(
        diagnostics.remappedMomentumKilogramMetersPerSecond,
        diagnostics.sourceMomentumKilogramMetersPerSecond);

    std::vector<double> normal(dofCount, 0.0);
    std::vector<double> material(dofCount, 0.0);
    std::vector<double> relative(dofCount, 0.0);
    for (const auto& dof : currentMetric.dofs) {
        const auto& link = currentTopology.links.at(
            dof.sourceFaceLinkIndex);
        if (link.stableId != dof.sourceFaceLinkStableId) {
            throw std::logic_error(
                "opening momentum-prediction link binding is invalid");
        }
        const double ownerNormal = component(
            transport.controls[sourceByCurrent[dof.ownerFragmentIndex]]
                .velocityMetersPerSecond,
            dof.axis);
        const double oppositeNormal = component(
            transport.controls[sourceByCurrent[dof.oppositeFragmentIndex]]
                .velocityMetersPerSecond,
            dof.axis);
        diagnostics.maximumEndpointNormalVelocityJumpMetersPerSecond =
            std::max(
                diagnostics
                    .maximumEndpointNormalVelocityJumpMetersPerSecond,
                std::abs(ownerNormal - oppositeNormal));
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            relative[dof.dofIndex] =
                0.5 * (ownerNormal + oppositeNormal);
            normal[dof.dofIndex] = relative[dof.dofIndex];
            ++diagnostics.sharedRegionGridDofCount;
        } else {
            if (link.kind
                != PlanarPressureRegionFragmentFaceKind::
                    PressureLayerWall) {
                throw std::logic_error(
                    "opening momentum-prediction wall degree is invalid");
            }
            const double minusMaterial = currentVolumeRates
                .fragments.at(link.minusFragmentIndex)
                .upperBoundaryVelocityMetersPerSecond;
            const double plusMaterial = currentVolumeRates
                .fragments.at(link.plusFragmentIndex)
                .lowerBoundaryVelocityMetersPerSecond;
            if (minusMaterial != plusMaterial) {
                throw std::invalid_argument(
                    "opening momentum-prediction wall motion disagrees");
            }
            material[dof.dofIndex] = minusMaterial;
            if (dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch) {
                relative[dof.dofIndex] =
                    0.5 * (ownerNormal + oppositeNormal)
                    - material[dof.dofIndex];
                ++diagnostics.openingDofCount;
                diagnostics
                    .maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond =
                    std::max(
                        diagnostics
                            .maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond,
                        std::abs(relative[dof.dofIndex]));
            } else {
                ++diagnostics.solidWallTraceDofCount;
            }
            normal[dof.dofIndex] = material[dof.dofIndex]
                + relative[dof.dofIndex];
        }
        diagnostics
            .maximumAbsolutePredictedRelativeVelocityMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsolutePredictedRelativeVelocityMetersPerSecond,
                std::abs(relative[dof.dofIndex]));
    }
    result.predictedVelocityState =
        buildPlanarPressureRegionFragmentOpeningVelocityState(
            currentMetric, normal, material, relative,
            result.densityKgPerCubicMeter, limits.stateLimits);
    diagnostics.predictedStateMomentumKilogramMetersPerSecond =
        result.predictedVelocityState.momentumKilogramMetersPerSecond;
    diagnostics.reconstructionMomentumChangeKilogramMetersPerSecond =
        subtract(
            diagnostics.predictedStateMomentumKilogramMetersPerSecond,
            diagnostics.remappedMomentumKilogramMetersPerSecond);
    diagnostics.predictedStateKineticEnergyJoules =
        result.predictedVelocityState.collocatedKineticEnergyJoules;
    diagnostics.finite = finiteVector(
            diagnostics.sourceMomentumKilogramMetersPerSecond)
        && finiteVector(
            diagnostics.remappedMomentumKilogramMetersPerSecond)
        && finiteVector(
            diagnostics.geometricMomentumChangeKilogramMetersPerSecond)
        && finiteVector(
            diagnostics.predictedStateMomentumKilogramMetersPerSecond)
        && finiteVector(
            diagnostics.reconstructionMomentumChangeKilogramMetersPerSecond)
        && std::isfinite(diagnostics.sourceKineticEnergyJoules)
        && std::isfinite(diagnostics.remappedKineticEnergyJoules)
        && std::isfinite(diagnostics.predictedStateKineticEnergyJoules)
        && std::isfinite(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters)
        && std::isfinite(diagnostics.maximumSourceSpeedMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumEndpointNormalVelocityJumpMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsolutePredictedRelativeVelocityMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "opening momentum-prediction diagnostics are non-finite");
    }
    result.ownedStorageBytes = result.predictedVelocityState.ownedStorageBytes;
    if (result.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "opening momentum-prediction owned storage disagrees");
    }
    result.fingerprint = predictionFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningMomentumPrediction
predictPlanarPressureRegionFragmentOpeningMomentum(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits& limits) {
    return buildPrediction(
        transport, transportTargetMetric, grid, currentSweep,
        currentFragments, currentTopology, currentVolumeRates,
        currentOpeningDefinitions, currentOpenings, currentBaseMetric,
        currentMetric, limits);
}

void validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction) {
    validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(
        prediction.predictedVelocityState);
    const auto& diagnostics = prediction.diagnostics;
    std::size_t sharedRegionGridDofCount = 0;
    std::size_t solidWallTraceDofCount = 0;
    std::size_t openingDofCount = 0;
    for (const auto& sample : prediction.predictedVelocityState.samples) {
        if (sample.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            ++sharedRegionGridDofCount;
        } else if (sample.kind
                   == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                       OpeningPatch) {
            ++openingDofCount;
        } else {
            ++solidWallTraceDofCount;
        }
    }
    if (prediction.version
            != planarPressureRegionFragmentOpeningMomentumPredictionVersion
        || prediction.fingerprint == 0
        || prediction.sourceTransportFingerprint == 0
        || prediction.sourceTransportMetricFingerprint == 0
        || prediction.currentMetricFingerprint == 0
        || prediction.currentVolumeRateFingerprint == 0
        || !std::isfinite(prediction.densityKgPerCubicMeter)
        || !(prediction.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(prediction.timeStepSeconds)
        || !(prediction.timeStepSeconds > 0.0)
        || prediction.predictedVelocityState.sourceMetricFingerprint
            != prediction.currentMetricFingerprint
        || prediction.predictedVelocityState.mappedFromAcceptedEndpoint
        || prediction.predictedVelocityState.sourceAcceptedStateFingerprint
            != 0
        || prediction.predictedVelocityState.sourceVolumeRateFingerprint
            != 0
        || prediction.predictedVelocityState.densityKgPerCubicMeter
            != prediction.densityKgPerCubicMeter
        || diagnostics.fragmentCount == 0
        || diagnostics.dofCount == 0
        || diagnostics.fragmentCount
            != prediction.predictedVelocityState.fragments.size()
        || diagnostics.dofCount
            != prediction.predictedVelocityState.samples.size()
        || diagnostics.sharedRegionGridDofCount
            != sharedRegionGridDofCount
        || diagnostics.solidWallTraceDofCount
            != solidWallTraceDofCount
        || diagnostics.openingDofCount != openingDofCount
        || diagnostics.sharedRegionGridDofCount
                + diagnostics.solidWallTraceDofCount
                + diagnostics.openingDofCount
            != diagnostics.dofCount
        || !finiteVector(
            diagnostics.sourceMomentumKilogramMetersPerSecond)
        || !finiteVector(
            diagnostics.remappedMomentumKilogramMetersPerSecond)
        || !finiteVector(
            diagnostics.geometricMomentumChangeKilogramMetersPerSecond)
        || diagnostics.geometricMomentumChangeKilogramMetersPerSecond
            != subtract(
                diagnostics.remappedMomentumKilogramMetersPerSecond,
                diagnostics.sourceMomentumKilogramMetersPerSecond)
        || !finiteVector(
            diagnostics.predictedStateMomentumKilogramMetersPerSecond)
        || diagnostics.predictedStateMomentumKilogramMetersPerSecond
            != prediction.predictedVelocityState
                   .momentumKilogramMetersPerSecond
        || !finiteVector(
            diagnostics.reconstructionMomentumChangeKilogramMetersPerSecond)
        || diagnostics.reconstructionMomentumChangeKilogramMetersPerSecond
            != subtract(
                diagnostics.predictedStateMomentumKilogramMetersPerSecond,
                diagnostics.remappedMomentumKilogramMetersPerSecond)
        || !std::isfinite(diagnostics.sourceKineticEnergyJoules)
        || diagnostics.sourceKineticEnergyJoules < 0.0
        || !std::isfinite(diagnostics.remappedKineticEnergyJoules)
        || diagnostics.remappedKineticEnergyJoules < 0.0
        || !std::isfinite(diagnostics.predictedStateKineticEnergyJoules)
        || diagnostics.predictedStateKineticEnergyJoules < 0.0
        || diagnostics.predictedStateKineticEnergyJoules
            != prediction.predictedVelocityState
                   .collocatedKineticEnergyJoules
        || !std::isfinite(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters)
        || diagnostics.maximumAbsoluteVolumeChangeCubicMeters < 0.0
        || !std::isfinite(diagnostics.maximumSourceSpeedMetersPerSecond)
        || diagnostics.maximumSourceSpeedMetersPerSecond < 0.0
        || !std::isfinite(
            diagnostics.maximumEndpointNormalVelocityJumpMetersPerSecond)
        || diagnostics.maximumEndpointNormalVelocityJumpMetersPerSecond < 0.0
        || !std::isfinite(
            diagnostics.maximumAbsolutePredictedRelativeVelocityMetersPerSecond)
        || diagnostics.maximumAbsolutePredictedRelativeVelocityMetersPerSecond
            < 0.0
        || !std::isfinite(
            diagnostics.maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond)
        || diagnostics
               .maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond
            < 0.0
        || !diagnostics.finite
        || prediction.ownedStorageBytes
            != prediction.predictedVelocityState.ownedStorageBytes
        || prediction.workingStorageBytes == 0
        || prediction.fingerprint != predictionFingerprint(prediction)) {
        throw std::invalid_argument(
            "opening momentum-prediction integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
        prediction);
    if (prediction.diagnostics.fragmentCount > limits.maximumFragments
        || prediction.diagnostics.dofCount > limits.maximumDofs
        || prediction.ownedStorageBytes > limits.maximumOwnedBytes
        || prediction.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum-prediction validation limit exceeded");
    }
    if (prediction != buildPrediction(
            transport, transportTargetMetric, grid, currentSweep,
            currentFragments, currentTopology, currentVolumeRates,
            currentOpeningDefinitions, currentOpenings, currentBaseMetric,
            currentMetric, limits)) {
        throw std::invalid_argument(
            "opening momentum prediction is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
