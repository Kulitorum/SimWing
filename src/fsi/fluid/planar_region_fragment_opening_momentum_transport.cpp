#include "fluid/planar_region_fragment_opening_momentum_transport.h"

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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
            "opening momentum-transport storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening momentum-transport storage overflows");
    }
    return first * second;
}

bool validSettings(
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings) {
    return std::isfinite(settings.maximumOutgoingCourantNumber)
        && settings.maximumOutgoingCourantNumber > 0.0
        && settings.maximumOutgoingCourantNumber <= 1.0
        && settings.maximumSubsteps > 0
        && std::isfinite(
            settings.absoluteContinuityToleranceCubicMetersPerSecond)
        && settings.absoluteContinuityToleranceCubicMetersPerSecond >= 0.0
        && std::isfinite(settings.relativeContinuityTolerance)
        && settings.relativeContinuityTolerance >= 0.0
        && std::isfinite(
            settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && settings.absoluteMomentumToleranceKilogramMetersPerSecond >= 0.0
        && std::isfinite(settings.relativeMomentumTolerance)
        && settings.relativeMomentumTolerance >= 0.0
        && std::isfinite(settings.absoluteEnergyToleranceJoules)
        && settings.absoluteEnergyToleranceJoules >= 0.0
        && std::isfinite(settings.relativeEnergyTolerance)
        && settings.relativeEnergyTolerance >= 0.0;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits) {
    if (limits.maximumFragments == 0 || limits.maximumDofs == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening momentum-transport limits are invalid");
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

Vector3 scale(const Vector3& value, const double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

double norm(const Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

double tolerance(const double absolute, const double relative,
                 const double scaleValue) {
    return std::max(absolute, relative * std::abs(scaleValue));
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

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport) {
    return checkedMultiply(
        transport.controls.size(),
        sizeof(
            PlanarPressureRegionFragmentOpeningMomentumTransportControl));
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings) {
    fingerprint.real(settings.maximumOutgoingCourantNumber);
    fingerprint.integer(static_cast<std::uint64_t>(settings.maximumSubsteps));
    fingerprint.real(
        settings.absoluteContinuityToleranceCubicMetersPerSecond);
    fingerprint.real(settings.relativeContinuityTolerance);
    fingerprint.real(
        settings.absoluteMomentumToleranceKilogramMetersPerSecond);
    fingerprint.real(settings.relativeMomentumTolerance);
    fingerprint.real(settings.absoluteEnergyToleranceJoules);
    fingerprint.real(settings.relativeEnergyTolerance);
}

std::uint64_t transportFingerprint(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport) {
    Fingerprint fingerprint;
    fingerprint.integer(transport.version);
    fingerprint.integer(transport.sourceStateFingerprint);
    fingerprint.integer(transport.sourceMetricFingerprint);
    fingerprint.integer(transport.targetFlowStateFingerprint);
    fingerprint.integer(transport.targetMetricFingerprint);
    fingerprint.integer(transport.targetVolumeRateFingerprint);
    fingerprint.real(transport.densityKgPerCubicMeter);
    fingerprint.real(transport.timeStepSeconds);
    fingerprintSettings(fingerprint, transport.settings);
    const auto& diagnostics = transport.diagnostics;
    for (const std::size_t value : {
             diagnostics.fragmentCount,
             diagnostics.transportDofCount,
             diagnostics.openingDofCount,
             diagnostics.substepCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    for (const double value : {
             diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
             diagnostics.maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond,
             diagnostics.maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond,
             diagnostics.maximumContinuityResidualCubicMetersPerSecond,
             diagnostics.continuityToleranceCubicMetersPerSecond,
             diagnostics.maximumFullStepOutgoingCourantNumber,
             diagnostics.maximumAcceptedSubstepOutgoingCourantNumber}) {
        fingerprint.real(value);
    }
    fingerprintVector(
        fingerprint, diagnostics.momentumBeforeKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, diagnostics.momentumAfterKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, diagnostics.momentumResidualKilogramMetersPerSecond);
    for (const double value : {
             diagnostics.momentumResidualNormKilogramMetersPerSecond,
             diagnostics.kineticEnergyBeforeJoules,
             diagnostics.kineticEnergyAfterJoules,
             diagnostics.advectiveKineticEnergyLossJoules,
             diagnostics.maximumVelocityChangeMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.enumeration(diagnostics.failureStage);
    fingerprint.boolean(diagnostics.finite);
    fingerprint.boolean(diagnostics.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(transport.controls.size()));
    for (const auto& control : transport.controls) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.fragmentIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(control.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.connectedComponentIndex));
        fingerprint.real(control.volumeCubicMeters);
        fingerprintVector(fingerprint, control.velocityMetersPerSecond);
        fingerprintVector(fingerprint, control.momentumKilogramMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        transport.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        transport.workingStorageBytes));
    return fingerprint.value();
}

Vector3 totalMomentum(
    const std::vector<
        PlanarPressureRegionFragmentOpeningMomentumTransportControl>&
        controls) {
    Vector3 result;
    for (const auto& control : controls)
        result = add(result, control.momentumKilogramMetersPerSecond);
    return result;
}

double kineticEnergy(
    const std::vector<
        PlanarPressureRegionFragmentOpeningMomentumTransportControl>& controls,
    const double density) {
    double result = 0.0;
    for (const auto& control : controls) {
        const double speed = norm(control.velocityMetersPerSecond);
        result += 0.5 * density * control.volumeCubicMeters
            * speed * speed;
    }
    return result;
}

std::size_t requiredSubsteps(const double fullStepCourant,
                             const double maximumCourant) {
    if (!(fullStepCourant > maximumCourant)) return 1;
    const double required = std::ceil(fullStepCourant / maximumCourant);
    if (!std::isfinite(required)
        || required > static_cast<double>(
            std::numeric_limits<std::size_t>::max())) {
        return std::numeric_limits<std::size_t>::max();
    }
    return static_cast<std::size_t>(required);
}

PlanarPressureRegionFragmentOpeningMomentumTransport buildTransport(
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityState& targetFlowState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& targetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& targetSweep,
    const PlanarPressureRegionFragmentSet& targetFragments,
    const PlanarPressureRegionFragmentTopology& targetTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& targetVolumeRates,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings,
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits) {
    validateLimits(limits);
    if (!validSettings(settings)) {
        throw std::invalid_argument(
            "opening momentum-transport settings are invalid");
    }
    validatePlanarPressureRegionFragmentOpeningVelocityState(
        sourceState, sourceMetric, limits.stateLimits);
    validatePlanarPressureRegionFragmentOpeningVelocityState(
        targetFlowState, targetMetric, limits.stateLimits);
    validatePlanarPressureRegionFragmentVolumeRates(
        targetVolumeRates, grid, targetSweep, targetFragments,
        targetTopology, limits.volumeRateLimits);
    const std::size_t fragmentCount = targetMetric.fragments.size();
    const std::size_t dofCount = targetMetric.dofs.size();
    if (fragmentCount == 0
        || sourceMetric.fragments.size() != fragmentCount
        || sourceMetric.dofs.size() != dofCount
        || sourceState.fragments.size() != fragmentCount
        || targetFlowState.fragments.size() != fragmentCount
        || targetVolumeRates.fragments.size() != fragmentCount
        || sourceMetric.components.size() != targetMetric.components.size()
        || sourceState.densityKgPerCubicMeter
            != targetFlowState.densityKgPerCubicMeter
        || sourceMetric.profileAxis != targetMetric.profileAxis
        || targetMetric.profileAxis != targetVolumeRates.axis
        || targetMetric.sourceFragmentFingerprint
            != targetFragments.fingerprint
        || targetMetric.sourceTopologyFingerprint
            != targetTopology.fingerprint) {
        throw std::invalid_argument(
            "opening momentum-transport endpoint identity is invalid");
    }
    if (fragmentCount > limits.maximumFragments
        || dofCount > limits.maximumDofs) {
        throw std::length_error(
            "opening momentum-transport entity limit exceeded");
    }
    const std::size_t expectedOwnedBytes = checkedMultiply(
        fragmentCount,
        sizeof(
            PlanarPressureRegionFragmentOpeningMomentumTransportControl));
    std::size_t workingBytes = checkedMultiply(
        checkedMultiply(fragmentCount, 2), sizeof(std::size_t));
    workingBytes = checkedAdd(workingBytes, checkedMultiply(
        fragmentCount, sizeof(std::size_t)));
    workingBytes = checkedAdd(workingBytes, checkedMultiply(
        dofCount, sizeof(double)));
    workingBytes = checkedAdd(workingBytes, checkedMultiply(
        fragmentCount, 3 * sizeof(double)));
    workingBytes = checkedAdd(workingBytes, checkedMultiply(
        fragmentCount,
        sizeof(
            PlanarPressureRegionFragmentOpeningMomentumTransportControl)));
    workingBytes = checkedAdd(workingBytes, checkedMultiply(
        fragmentCount, sizeof(Vector3)));
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum-transport storage limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningMomentumTransport result;
    result.sourceStateFingerprint = sourceState.fingerprint;
    result.sourceMetricFingerprint = sourceMetric.fingerprint;
    result.targetFlowStateFingerprint = targetFlowState.fingerprint;
    result.targetMetricFingerprint = targetMetric.fingerprint;
    result.targetVolumeRateFingerprint = targetVolumeRates.fingerprint;
    result.densityKgPerCubicMeter = sourceState.densityKgPerCubicMeter;
    result.timeStepSeconds = targetVolumeRates.durationSeconds;
    result.settings = settings;
    result.workingStorageBytes = workingBytes;
    auto& diagnostics = result.diagnostics;
    diagnostics.fragmentCount = fragmentCount;

    std::vector<std::size_t> sourceOrder(fragmentCount);
    std::vector<std::size_t> targetOrder(fragmentCount);
    std::iota(sourceOrder.begin(), sourceOrder.end(), 0);
    std::iota(targetOrder.begin(), targetOrder.end(), 0);
    std::ranges::sort(sourceOrder, {}, [&](const std::size_t index) {
        return sourceState.fragments[index].stableId;
    });
    std::ranges::sort(targetOrder, {}, [&](const std::size_t index) {
        return targetFlowState.fragments[index].stableId;
    });
    std::vector<std::size_t> sourceByTarget(fragmentCount);
    for (std::size_t offset = 0; offset < fragmentCount; ++offset) {
        const std::size_t sourceIndex = sourceOrder[offset];
        const std::size_t targetIndex = targetOrder[offset];
        const auto& source = sourceState.fragments[sourceIndex];
        const auto& target = targetFlowState.fragments[targetIndex];
        const auto& rate = targetVolumeRates.fragments[targetIndex];
        if (source.stableId == 0 || source.stableId != target.stableId
            || source.regionStableId != target.regionStableId
            || source.connectedComponentIndex
                != target.connectedComponentIndex) {
            throw std::invalid_argument(
                "opening momentum-transport stable fragment mapping is invalid");
        }
        if (rate.fragmentIndex != targetIndex
            || rate.stableId != target.stableId
            || rate.regionStableId != target.regionStableId) {
            throw std::invalid_argument(
                "opening momentum-transport volume-rate mapping is invalid");
        }
        if (!equalWithinRoundoff(
                rate.previousVolumeCubicMeters,
                source.volumeCubicMeters)
            || !equalWithinRoundoff(
                rate.currentVolumeCubicMeters,
                target.volumeCubicMeters)) {
            throw std::invalid_argument(
                "opening momentum-transport endpoint volume mapping is invalid");
        }
        sourceByTarget[targetIndex] = sourceIndex;
    }
    for (std::size_t index = 0;
         index < sourceMetric.components.size(); ++index) {
        if (sourceMetric.components[index].stableId
                != targetMetric.components[index].stableId
            || sourceMetric.components[index].baseComponentCount
                != targetMetric.components[index].baseComponentCount) {
            throw std::invalid_argument(
                "opening momentum-transport component mapping is invalid");
        }
    }
    for (std::size_t index = 0; index < dofCount; ++index) {
        const auto& source = sourceMetric.dofs[index];
        const auto& target = targetMetric.dofs[index];
        if (source.dofIndex != index || target.dofIndex != index
            || source.stableId == 0
            || source.stableId != target.stableId
            || source.kind != target.kind
            || source.sourceFaceLinkStableId
                != target.sourceFaceLinkStableId
            || source.sourceOpeningPatchStableId
                != target.sourceOpeningPatchStableId
            || source.axis != target.axis
            || source.surfaceStableId != target.surfaceStableId
            || source.ownerFragmentStableId
                != target.ownerFragmentStableId
            || source.oppositeFragmentStableId
                != target.oppositeFragmentStableId
            || source.connectedComponentIndex
                != target.connectedComponentIndex) {
            throw std::invalid_argument(
                "opening momentum-transport degree mapping is invalid");
        }
    }

    std::vector<double> relativeFlows(dofCount, 0.0);
    std::vector<double> netOutward(fragmentCount, 0.0);
    std::vector<double> outward(fragmentCount, 0.0);
    std::vector<double> geometryRates(fragmentCount, 0.0);
    double continuityScale = 0.0;
    for (std::size_t index = 0; index < fragmentCount; ++index) {
        const double rate = targetVolumeRates.fragments[index]
            .geometryVolumeChangeRateCubicMetersPerSecond;
        if (!std::isfinite(rate)) {
            throw std::invalid_argument(
                "opening momentum-transport geometry rate is non-finite");
        }
        geometryRates[index] = rate;
        diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
                std::abs(rate));
        continuityScale = std::max(continuityScale, std::abs(rate));
    }
    for (const auto& dof : targetMetric.dofs) {
        const bool transports = dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid
            || dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch;
        if (!transports) continue;
        const double flow = dof.areaSquareMeters
            * targetFlowState.samples[dof.dofIndex]
                .relativeNormalVelocityMetersPerSecond;
        if (!std::isfinite(flow)) {
            throw std::invalid_argument(
                "opening momentum-transport relative flow is non-finite");
        }
        relativeFlows[dof.dofIndex] = flow;
        netOutward[dof.ownerFragmentIndex] += flow;
        netOutward[dof.oppositeFragmentIndex] -= flow;
        outward[dof.ownerFragmentIndex] += std::max(0.0, flow);
        outward[dof.oppositeFragmentIndex] += std::max(0.0, -flow);
        ++diagnostics.transportDofCount;
        diagnostics.openingDofCount += dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                OpeningPatch ? 1 : 0;
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                OpeningPatch) {
            diagnostics
                .maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond =
                std::max(
                    diagnostics
                        .maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond,
                    std::abs(flow));
        }
        diagnostics.maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(flow));
        continuityScale = std::max(continuityScale, std::abs(flow));
    }
    diagnostics.continuityToleranceCubicMetersPerSecond = tolerance(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance, continuityScale);
    for (std::size_t index = 0; index < fragmentCount; ++index) {
        const double residual = geometryRates[index] + netOutward[index];
        diagnostics.maximumContinuityResidualCubicMetersPerSecond = std::max(
            diagnostics.maximumContinuityResidualCubicMetersPerSecond,
            std::abs(residual));
        const double previousVolume = sourceState.fragments[
            sourceByTarget[index]].volumeCubicMeters;
        const double currentVolume =
            targetFlowState.fragments[index].volumeCubicMeters;
        if (!(previousVolume > 0.0) || !(currentVolume > 0.0)) {
            diagnostics.failureStage =
                PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                    GeometryVolume;
        }
        diagnostics.maximumFullStepOutgoingCourantNumber = std::max(
            diagnostics.maximumFullStepOutgoingCourantNumber,
            result.timeStepSeconds * outward[index]
                / std::min(previousVolume, currentVolume));
    }
    if (diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && diagnostics.maximumContinuityResidualCubicMetersPerSecond
            > diagnostics.continuityToleranceCubicMetersPerSecond) {
        diagnostics.failureStage =
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                FlowContinuity;
    }
    diagnostics.substepCount = requiredSubsteps(
        diagnostics.maximumFullStepOutgoingCourantNumber,
        settings.maximumOutgoingCourantNumber);
    if (diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && diagnostics.substepCount > settings.maximumSubsteps) {
        diagnostics.failureStage =
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                SubstepLimit;
    }
    if (diagnostics.substepCount == 0) diagnostics.substepCount = 1;
    diagnostics.maximumAcceptedSubstepOutgoingCourantNumber =
        diagnostics.maximumFullStepOutgoingCourantNumber
        / static_cast<double>(diagnostics.substepCount);

    std::vector<
        PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        candidate;
    candidate.reserve(fragmentCount);
    for (std::size_t targetIndex = 0;
         targetIndex < fragmentCount; ++targetIndex) {
        const auto& source =
            sourceState.fragments[sourceByTarget[targetIndex]];
        candidate.push_back({
            targetIndex,
            targetFlowState.fragments[targetIndex].stableId,
            targetFlowState.fragments[targetIndex].regionStableId,
            targetFlowState.fragments[targetIndex].connectedComponentIndex,
            source.volumeCubicMeters,
            source.collocatedVelocityMetersPerSecond,
            source.momentumKilogramMetersPerSecond,
        });
    }
    diagnostics.momentumBeforeKilogramMetersPerSecond =
        totalMomentum(candidate);
    diagnostics.kineticEnergyBeforeJoules = kineticEnergy(
        candidate, result.densityKgPerCubicMeter);

    if (diagnostics.failureStage
        == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
            None) {
        std::vector<Vector3> impulse(fragmentCount);
        const double substepSeconds = result.timeStepSeconds
            / static_cast<double>(diagnostics.substepCount);
        double previousEnergy = diagnostics.kineticEnergyBeforeJoules;
        for (std::size_t substep = 0;
             substep < diagnostics.substepCount; ++substep) {
            std::ranges::fill(impulse, Vector3{});
            for (const auto& dof : targetMetric.dofs) {
                const double flow = relativeFlows[dof.dofIndex];
                if (flow == 0.0) continue;
                const auto& donor = flow >= 0.0
                    ? candidate[dof.ownerFragmentIndex]
                    : candidate[dof.oppositeFragmentIndex];
                const Vector3 transported = scale(
                    donor.velocityMetersPerSecond,
                    result.densityKgPerCubicMeter * flow * substepSeconds);
                impulse[dof.ownerFragmentIndex] = subtract(
                    impulse[dof.ownerFragmentIndex], transported);
                impulse[dof.oppositeFragmentIndex] = add(
                    impulse[dof.oppositeFragmentIndex], transported);
            }
            const double elapsed = substep + 1 == diagnostics.substepCount
                ? result.timeStepSeconds
                : substepSeconds * static_cast<double>(substep + 1);
            for (std::size_t index = 0; index < fragmentCount; ++index) {
                candidate[index].momentumKilogramMetersPerSecond = add(
                    candidate[index].momentumKilogramMetersPerSecond,
                    impulse[index]);
                const double sourceVolume = sourceState.fragments[
                    sourceByTarget[index]].volumeCubicMeters;
                const double targetVolume =
                    targetFlowState.fragments[index].volumeCubicMeters;
                candidate[index].volumeCubicMeters =
                    substep + 1 == diagnostics.substepCount
                    ? targetVolume
                    : sourceVolume
                        + (targetVolume - sourceVolume)
                            * (elapsed / result.timeStepSeconds);
                candidate[index].velocityMetersPerSecond = scale(
                    candidate[index].momentumKilogramMetersPerSecond,
                    1.0 / (result.densityKgPerCubicMeter
                           * candidate[index].volumeCubicMeters));
            }
            const double energy = kineticEnergy(
                candidate, result.densityKgPerCubicMeter);
            if (!std::isfinite(energy)
                || energy > previousEnergy
                    + tolerance(
                        settings.absoluteEnergyToleranceJoules,
                        settings.relativeEnergyTolerance,
                        previousEnergy)) {
                diagnostics.failureStage =
                    PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                        AdvectionEnergy;
                break;
            }
            previousEnergy = energy;
        }
    }
    diagnostics.momentumAfterKilogramMetersPerSecond =
        totalMomentum(candidate);
    diagnostics.momentumResidualKilogramMetersPerSecond = subtract(
        diagnostics.momentumAfterKilogramMetersPerSecond,
        diagnostics.momentumBeforeKilogramMetersPerSecond);
    diagnostics.momentumResidualNormKilogramMetersPerSecond = norm(
        diagnostics.momentumResidualKilogramMetersPerSecond);
    diagnostics.kineticEnergyAfterJoules = kineticEnergy(
        candidate, result.densityKgPerCubicMeter);
    diagnostics.advectiveKineticEnergyLossJoules =
        diagnostics.kineticEnergyBeforeJoules
        - diagnostics.kineticEnergyAfterJoules;
    const double momentumScale = std::max(
        norm(diagnostics.momentumBeforeKilogramMetersPerSecond),
        norm(diagnostics.momentumAfterKilogramMetersPerSecond));
    if (diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && diagnostics.momentumResidualNormKilogramMetersPerSecond
            > tolerance(
                settings.absoluteMomentumToleranceKilogramMetersPerSecond,
                settings.relativeMomentumTolerance, momentumScale)) {
        diagnostics.failureStage =
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                Conservation;
    }
    for (std::size_t index = 0; index < fragmentCount; ++index) {
        diagnostics.maximumVelocityChangeMetersPerSecond = std::max(
            diagnostics.maximumVelocityChangeMetersPerSecond,
            norm(subtract(
                candidate[index].velocityMetersPerSecond,
                sourceState.fragments[sourceByTarget[index]]
                    .collocatedVelocityMetersPerSecond)));
    }
    diagnostics.finite = finiteVector(
            diagnostics.momentumBeforeKilogramMetersPerSecond)
        && finiteVector(diagnostics.momentumAfterKilogramMetersPerSecond)
        && finiteVector(diagnostics.momentumResidualKilogramMetersPerSecond)
        && std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.advectiveKineticEnergyLossJoules)
        && std::isfinite(diagnostics.maximumVelocityChangeMetersPerSecond)
        && std::ranges::all_of(candidate, [](const auto& control) {
            return std::isfinite(control.volumeCubicMeters)
                && control.volumeCubicMeters > 0.0
                && finiteVector(control.velocityMetersPerSecond)
                && finiteVector(control.momentumKilogramMetersPerSecond);
        });
    if (!diagnostics.finite
        && diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None) {
        diagnostics.failureStage =
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                NonFinite;
    }
    if (diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && diagnostics.finite) {
        diagnostics.accepted = true;
        result.controls = std::move(candidate);
        result.ownedStorageBytes = expectedOwnedBytes;
    }
    result.fingerprint = transportFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningMomentumTransport
advancePlanarPressureRegionFragmentOpeningMomentum(
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityState& targetFlowState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& targetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& targetSweep,
    const PlanarPressureRegionFragmentSet& targetFragments,
    const PlanarPressureRegionFragmentTopology& targetTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& targetVolumeRates,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings,
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits) {
    return buildTransport(
        sourceState, sourceMetric, targetFlowState, targetMetric, grid,
        targetSweep, targetFragments, targetTopology, targetVolumeRates,
        settings, limits);
}

void validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport) {
    const auto& diagnostics = transport.diagnostics;
    bool controlsValid = true;
    for (std::size_t index = 0;
         index < transport.controls.size(); ++index) {
        const auto& control = transport.controls[index];
        const Vector3 reconstructedMomentum = scale(
            control.velocityMetersPerSecond,
            transport.densityKgPerCubicMeter
                * control.volumeCubicMeters);
        controlsValid = controlsValid
            && control.fragmentIndex == index
            && control.stableId != 0
            && control.regionStableId != 0
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && finiteVector(control.velocityMetersPerSecond)
            && finiteVector(control.momentumKilogramMetersPerSecond)
            && norm(subtract(
                   control.momentumKilogramMetersPerSecond,
                   reconstructedMomentum))
                <= tolerance(
                    transport.settings
                        .absoluteMomentumToleranceKilogramMetersPerSecond,
                    transport.settings.relativeMomentumTolerance,
                    norm(control.momentumKilogramMetersPerSecond));
    }
    const bool acceptedShape = diagnostics.accepted
        && diagnostics.failureStage
            == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && transport.controls.size() == diagnostics.fragmentCount
        && transport.ownedStorageBytes == ownedStorageBytes(transport);
    const bool rejectedShape = !diagnostics.accepted
        && diagnostics.failureStage
            != PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None
        && transport.controls.empty() && transport.ownedStorageBytes == 0;
    if (transport.version
            != planarPressureRegionFragmentOpeningMomentumTransportVersion
        || transport.fingerprint == 0
        || transport.sourceStateFingerprint == 0
        || transport.sourceMetricFingerprint == 0
        || transport.targetFlowStateFingerprint == 0
        || transport.targetMetricFingerprint == 0
        || transport.targetVolumeRateFingerprint == 0
        || !std::isfinite(transport.densityKgPerCubicMeter)
        || !(transport.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(transport.timeStepSeconds)
        || !(transport.timeStepSeconds > 0.0)
        || !validSettings(transport.settings)
        || diagnostics.fragmentCount == 0
        || diagnostics.transportDofCount == 0
        || diagnostics.openingDofCount > diagnostics.transportDofCount
        || diagnostics.substepCount == 0
        || !std::isfinite(
            diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumContinuityResidualCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.continuityToleranceCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumFullStepOutgoingCourantNumber)
        || !std::isfinite(
            diagnostics.maximumAcceptedSubstepOutgoingCourantNumber)
        || !finiteVector(diagnostics.momentumBeforeKilogramMetersPerSecond)
        || !finiteVector(diagnostics.momentumAfterKilogramMetersPerSecond)
        || !finiteVector(diagnostics.momentumResidualKilogramMetersPerSecond)
        || !std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        || !std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        || !std::isfinite(diagnostics.kineticEnergyAfterJoules)
        || !std::isfinite(diagnostics.advectiveKineticEnergyLossJoules)
        || !std::isfinite(diagnostics.maximumVelocityChangeMetersPerSecond)
        || !diagnostics.finite
        || !controlsValid
        || (diagnostics.accepted
            && diagnostics.substepCount > transport.settings.maximumSubsteps)
        || (diagnostics.accepted
            && diagnostics.maximumAcceptedSubstepOutgoingCourantNumber
                > transport.settings.maximumOutgoingCourantNumber
                    + 32.0 * std::numeric_limits<double>::epsilon()
                        * std::max(
                            1.0,
                            diagnostics
                                .maximumAcceptedSubstepOutgoingCourantNumber))
        || (diagnostics.accepted
            && totalMomentum(transport.controls)
                != diagnostics.momentumAfterKilogramMetersPerSecond)
        || (diagnostics.accepted
            && kineticEnergy(
                   transport.controls,
                   transport.densityKgPerCubicMeter)
                != diagnostics.kineticEnergyAfterJoules)
        || (!acceptedShape && !rejectedShape)
        || transport.workingStorageBytes == 0
        || transport.fingerprint != transportFingerprint(transport)) {
        throw std::invalid_argument(
            "opening momentum-transport integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumTransport(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityState& targetFlowState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& targetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& targetSweep,
    const PlanarPressureRegionFragmentSet& targetFragments,
    const PlanarPressureRegionFragmentTopology& targetTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& targetVolumeRates,
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        transport);
    if (transport.controls.size() > limits.maximumFragments
        || transport.ownedStorageBytes > limits.maximumOwnedBytes
        || transport.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening momentum-transport validation limit exceeded");
    }
    if (transport != buildTransport(
            sourceState, sourceMetric, targetFlowState, targetMetric, grid,
            targetSweep, targetFragments, targetTopology, targetVolumeRates,
            transport.settings, limits)) {
        throw std::invalid_argument(
            "opening momentum transport is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
