#include "fluid/planar_region_fragment_opening_momentum_adjustment_state.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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
        for (std::size_t index = 0; index < sizeof(value); ++index) {
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

Vector3 add(const Vector3& first, const Vector3& second) {
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

Vector3 subtract(const Vector3& first, const Vector3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

Vector3 scale(const Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double norm(const Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

bool finite(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double tolerance(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings,
    const double reference) {
    return settings.absoluteMomentumToleranceKilogramMetersPerSecond
        + settings.relativeMomentumTolerance
            * std::max(1.0, std::abs(reference));
}

bool validSettings(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings) {
    return std::isfinite(
               settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && settings.absoluteMomentumToleranceKilogramMetersPerSecond >= 0.0
        && std::isfinite(settings.relativeMomentumTolerance)
        && settings.relativeMomentumTolerance >= 0.0;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    if (limits.maximumFragments == 0 || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "opening momentum-adjustment state limits are invalid");
    }
}

std::size_t ownedStorageBytes(const std::size_t controls) {
    if (controls > std::numeric_limits<std::size_t>::max()
                       / sizeof(
                           PlanarPressureRegionFragmentOpeningMomentumTransportControl)) {
        throw std::length_error(
            "opening momentum-adjustment state storage overflows");
    }
    return controls
        * sizeof(
            PlanarPressureRegionFragmentOpeningMomentumTransportControl);
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceTransportFingerprint);
    fingerprint.integer(state.sourceAdjustmentFingerprint);
    fingerprint.integer(state.sourceMetricFingerprint);
    fingerprint.real(state.densityKgPerCubicMeter);
    fingerprint.real(state.timeStepSeconds);
    fingerprint.real(
        state.settings.absoluteMomentumToleranceKilogramMetersPerSecond);
    fingerprint.real(state.settings.relativeMomentumTolerance);
    fingerprintVector(
        fingerprint, state.sourceMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, state.adjustedMomentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, state.adjustmentImpulseKilogramMetersPerSecond);
    fingerprint.real(state.sourceKineticEnergyJoules);
    fingerprint.real(state.adjustedKineticEnergyJoules);
    fingerprint.real(state.kineticEnergyChangeJoules);
    fingerprint.real(
        state.maximumMomentumVelocityResidualKilogramMetersPerSecond);
    fingerprint.integer(static_cast<std::uint64_t>(state.controls.size()));
    for (const auto& control : state.controls) {
        fingerprint.integer(
            static_cast<std::uint64_t>(control.fragmentIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(control.regionStableId);
        fingerprint.integer(
            static_cast<std::uint64_t>(control.connectedComponentIndex));
        fingerprint.real(control.volumeCubicMeters);
        fingerprintVector(fingerprint, control.velocityMetersPerSecond);
        fingerprintVector(fingerprint, control.momentumKilogramMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentOpeningMomentumAdjustmentState buildState(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const std::uint64_t sourceAdjustmentFingerprint,
    const std::span<
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    validateLimits(limits);
    if (!validSettings(settings) || sourceAdjustmentFingerprint == 0) {
        throw std::invalid_argument(
            "opening momentum-adjustment state settings or source are invalid");
    }
    validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        transport);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(metric);
    if (!transport.diagnostics.accepted
        || transport.targetMetricFingerprint != metric.fingerprint
        || transport.controls.size() != metric.fragments.size()
        || adjustedControls.size() != transport.controls.size()
        || adjustedControls.empty()) {
        throw std::invalid_argument(
            "opening momentum-adjustment state endpoint is invalid");
    }
    const std::size_t storage = ownedStorageBytes(adjustedControls.size());
    if (adjustedControls.size() > limits.maximumFragments
        || storage > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening momentum-adjustment state limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningMomentumAdjustmentState result;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.sourceAdjustmentFingerprint = sourceAdjustmentFingerprint;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.densityKgPerCubicMeter = transport.densityKgPerCubicMeter;
    result.timeStepSeconds = transport.timeStepSeconds;
    result.settings = settings;
    result.sourceMomentumKilogramMetersPerSecond =
        transport.diagnostics.momentumAfterKilogramMetersPerSecond;
    result.sourceKineticEnergyJoules =
        transport.diagnostics.kineticEnergyAfterJoules;
    result.controls.reserve(adjustedControls.size());
    for (std::size_t index = 0; index < adjustedControls.size(); ++index) {
        const auto& source = transport.controls[index];
        const auto& adjusted = adjustedControls[index];
        const auto& metricFragment = metric.fragments[index];
        if (source.fragmentIndex != index
            || adjusted.fragmentIndex != index
            || metricFragment.fragmentIndex != index
            || source.stableId == 0
            || adjusted.stableId != source.stableId
            || metricFragment.stableId != source.stableId
            || adjusted.regionStableId != source.regionStableId
            || metricFragment.regionStableId != source.regionStableId
            || adjusted.connectedComponentIndex
                != source.connectedComponentIndex
            || metricFragment.connectedComponentIndex
                != source.connectedComponentIndex
            || adjusted.volumeCubicMeters != source.volumeCubicMeters
            || metricFragment.sourceVolumeCubicMeters
                != source.volumeCubicMeters
            || !finite(adjusted.velocityMetersPerSecond)
            || !finite(adjusted.momentumKilogramMetersPerSecond)) {
            throw std::invalid_argument(
                "opening momentum-adjustment control identity is invalid");
        }
        const Vector3 reconstructed = scale(
            adjusted.velocityMetersPerSecond,
            result.densityKgPerCubicMeter * adjusted.volumeCubicMeters);
        const double residual = norm(subtract(
            adjusted.momentumKilogramMetersPerSecond, reconstructed));
        result.maximumMomentumVelocityResidualKilogramMetersPerSecond =
            std::max(
                result.maximumMomentumVelocityResidualKilogramMetersPerSecond,
                residual);
        if (residual > tolerance(
                settings, norm(adjusted.momentumKilogramMetersPerSecond))) {
            throw std::invalid_argument(
                "opening momentum-adjustment control momentum is inconsistent");
        }
        result.adjustedMomentumKilogramMetersPerSecond = add(
            result.adjustedMomentumKilogramMetersPerSecond,
            adjusted.momentumKilogramMetersPerSecond);
        const double speed = norm(adjusted.velocityMetersPerSecond);
        result.adjustedKineticEnergyJoules +=
            0.5 * result.densityKgPerCubicMeter
            * adjusted.volumeCubicMeters * speed * speed;
        result.controls.push_back(adjusted);
    }
    result.adjustmentImpulseKilogramMetersPerSecond = subtract(
        result.adjustedMomentumKilogramMetersPerSecond,
        result.sourceMomentumKilogramMetersPerSecond);
    result.kineticEnergyChangeJoules =
        result.adjustedKineticEnergyJoules
        - result.sourceKineticEnergyJoules;
    result.ownedStorageBytes = storage;
    result.fingerprint = stateFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
        result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningMomentumAdjustmentState
capturePlanarPressureRegionFragmentOpeningMomentumAdjustmentState(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const std::uint64_t sourceAdjustmentFingerprint,
    const std::span<
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    return buildState(
        transport, metric, sourceAdjustmentFingerprint, adjustedControls,
        settings, limits);
}

void validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state) {
    Vector3 adjustedMomentum;
    double adjustedEnergy = 0.0;
    double maximumResidual = 0.0;
    bool controlsValid = true;
    for (std::size_t index = 0; index < state.controls.size(); ++index) {
        const auto& control = state.controls[index];
        const Vector3 reconstructed = scale(
            control.velocityMetersPerSecond,
            state.densityKgPerCubicMeter * control.volumeCubicMeters);
        const double residual = norm(subtract(
            control.momentumKilogramMetersPerSecond, reconstructed));
        controlsValid = controlsValid
            && control.fragmentIndex == index && control.stableId != 0
            && control.regionStableId != 0
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && finite(control.velocityMetersPerSecond)
            && finite(control.momentumKilogramMetersPerSecond)
            && residual <= tolerance(
                state.settings,
                norm(control.momentumKilogramMetersPerSecond));
        maximumResidual = std::max(maximumResidual, residual);
        adjustedMomentum = add(
            adjustedMomentum, control.momentumKilogramMetersPerSecond);
        const double speed = norm(control.velocityMetersPerSecond);
        adjustedEnergy += 0.5 * state.densityKgPerCubicMeter
            * control.volumeCubicMeters * speed * speed;
    }
    if (state.version
            != planarPressureRegionFragmentOpeningMomentumAdjustmentStateVersion
        || state.fingerprint == 0
        || state.sourceTransportFingerprint == 0
        || state.sourceAdjustmentFingerprint == 0
        || state.sourceMetricFingerprint == 0
        || !(state.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(state.densityKgPerCubicMeter)
        || !(state.timeStepSeconds > 0.0)
        || !std::isfinite(state.timeStepSeconds)
        || !validSettings(state.settings) || state.controls.empty()
        || !finite(state.sourceMomentumKilogramMetersPerSecond)
        || !finite(state.adjustedMomentumKilogramMetersPerSecond)
        || !finite(state.adjustmentImpulseKilogramMetersPerSecond)
        || state.adjustedMomentumKilogramMetersPerSecond != adjustedMomentum
        || state.adjustmentImpulseKilogramMetersPerSecond
            != subtract(
                state.adjustedMomentumKilogramMetersPerSecond,
                state.sourceMomentumKilogramMetersPerSecond)
        || !std::isfinite(state.sourceKineticEnergyJoules)
        || state.sourceKineticEnergyJoules < 0.0
        || !std::isfinite(state.adjustedKineticEnergyJoules)
        || state.adjustedKineticEnergyJoules < 0.0
        || state.adjustedKineticEnergyJoules != adjustedEnergy
        || state.kineticEnergyChangeJoules
            != state.adjustedKineticEnergyJoules
                - state.sourceKineticEnergyJoules
        || state.maximumMomentumVelocityResidualKilogramMetersPerSecond
            != maximumResidual
        || !controlsValid
        || state.ownedStorageBytes != ownedStorageBytes(state.controls.size())
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening momentum-adjustment state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentState(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const std::uint64_t sourceAdjustmentFingerprint,
    const std::span<
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
        state);
    if (state.controls.size() > limits.maximumFragments
        || state.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening momentum-adjustment state validation limit exceeded");
    }
    if (state != buildState(
            transport, metric, sourceAdjustmentFingerprint,
            adjustedControls, state.settings, limits)) {
        throw std::invalid_argument(
            "opening momentum-adjustment state is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
