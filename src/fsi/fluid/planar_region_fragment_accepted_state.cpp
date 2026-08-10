#include "fluid/planar_region_fragment_accepted_state.h"

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
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
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

void validateLimits(
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar regional accepted-state limits are invalid");
    }
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional accepted-state storage overflows");
    }
    return first + second;
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentAcceptedState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceVelocityStateFingerprint);
    fingerprint.integer(state.sourcePressureStateFingerprint);
    fingerprint.integer(state.sourceSurfaceLoadFingerprint);
    fingerprint.integer(state.sourceFragmentFingerprint);
    fingerprint.integer(state.sourceTopologyFingerprint);
    fingerprint.integer(state.sourceMetricFingerprint);
    fingerprint.boolean(state.staticGeometry);
    fingerprint.boolean(state.usesMovingVolumeRates);
    fingerprint.real(state.timeStepSeconds);
    fingerprint.real(state.densityKgPerCubicMeter);
    fingerprintVector(
        fingerprint, state.fluidMomentumKilogramMetersPerSecond);
    fingerprint.real(state.fluidKineticEnergyJoules);
    fingerprintVector(
        fingerprint, state.pressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, state.pressureImpulseOnSheetNewtonSeconds);
    fingerprint.real(state.pressureWorkToSheetJoules);
    fingerprint.integer(state.velocity.fingerprint);
    fingerprint.integer(state.pressure.fingerprint);
    fingerprint.integer(state.surfaceLoads.fingerprint);
    fingerprint.boolean(state.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        state.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentAcceptedState buildState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocity,
    const PlanarPressureRegionFragmentPressureState& pressure,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentVelocityState(
        velocity, grid, sweep, fragments, topology, metric,
        limits.velocityStateLimits);
    validatePlanarPressureRegionFragmentPressureStateIntegrity(pressure);
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressure, limits.surfaceLoadLimits);
    if (pressure.sourceAfterVelocityStateFingerprint != velocity.fingerprint
        || pressure.sourceFragmentFingerprint != fragments.fingerprint
        || pressure.sourceTopologyFingerprint != topology.fingerprint
        || pressure.sourceMetricFingerprint != metric.fingerprint
        || surfaceLoads.sourcePressureStateFingerprint
            != pressure.fingerprint
        || surfaceLoads.sourceTopologyFingerprint != topology.fingerprint
        || pressure.staticGeometry != surfaceLoads.staticGeometry
        || pressure.usesMovingVolumeRates
            != surfaceLoads.usesMovingVolumeRates
        || pressure.timeStepSeconds != surfaceLoads.timeStepSeconds) {
        throw std::invalid_argument(
            "planar regional accepted-state sources are incompatible");
    }

    PlanarPressureRegionFragmentAcceptedState result;
    result.sourceVelocityStateFingerprint = velocity.fingerprint;
    result.sourcePressureStateFingerprint = pressure.fingerprint;
    result.sourceSurfaceLoadFingerprint = surfaceLoads.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.staticGeometry = pressure.staticGeometry;
    result.usesMovingVolumeRates = pressure.usesMovingVolumeRates;
    result.timeStepSeconds = pressure.timeStepSeconds;
    result.densityKgPerCubicMeter = velocity.densityKgPerCubicMeter;
    result.fluidMomentumKilogramMetersPerSecond =
        velocity.momentumKilogramMetersPerSecond;
    result.fluidKineticEnergyJoules = velocity.kineticEnergyJoules;
    result.pressureForceOnSheetNewtons =
        surfaceLoads.totalPressureForceOnSheetNewtons;
    result.pressureImpulseOnSheetNewtonSeconds =
        surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds;
    result.pressureWorkToSheetJoules =
        surfaceLoads.totalPressureWorkToSheetJoules;
    result.ownedStorageBytes = checkedSum(
        checkedSum(
            velocity.ownedStorageBytes, pressure.ownedStorageBytes),
        surfaceLoads.ownedStorageBytes);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional accepted-state storage limit exceeded");
    }
    result.velocity = velocity;
    result.pressure = pressure;
    result.surfaceLoads = surfaceLoads;
    result.accepted = true;
    result.fingerprint = stateFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentAcceptedState
capturePlanarPressureRegionFragmentAcceptedState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocity,
    const PlanarPressureRegionFragmentPressureState& pressure,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits) {
    return buildState(
        grid, sweep, fragments, topology, metric, velocity, pressure,
        surfaceLoads, limits);
}

void validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
    const PlanarPressureRegionFragmentAcceptedState& state) {
    validatePlanarPressureRegionFragmentVelocityStateIntegrity(state.velocity);
    validatePlanarPressureRegionFragmentPressureStateIntegrity(state.pressure);
    validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
        state.surfaceLoads);
    if (state.version != planarPressureRegionFragmentAcceptedStateVersion
        || state.fingerprint == 0 || !state.accepted
        || state.sourceVelocityStateFingerprint != state.velocity.fingerprint
        || state.sourcePressureStateFingerprint != state.pressure.fingerprint
        || state.sourceSurfaceLoadFingerprint
            != state.surfaceLoads.fingerprint
        || state.sourceFragmentFingerprint
            != state.velocity.sourceFragmentFingerprint
        || state.sourceFragmentFingerprint
            != state.pressure.sourceFragmentFingerprint
        || state.sourceTopologyFingerprint
            != state.velocity.sourceTopologyFingerprint
        || state.sourceTopologyFingerprint
            != state.pressure.sourceTopologyFingerprint
        || state.sourceMetricFingerprint
            != state.velocity.sourceMetricFingerprint
        || state.sourceMetricFingerprint
            != state.pressure.sourceMetricFingerprint
        || state.pressure.sourceAfterVelocityStateFingerprint
            != state.velocity.fingerprint
        || state.surfaceLoads.sourcePressureStateFingerprint
            != state.pressure.fingerprint
        || state.staticGeometry != state.pressure.staticGeometry
        || state.staticGeometry != state.surfaceLoads.staticGeometry
        || state.usesMovingVolumeRates
            != state.pressure.usesMovingVolumeRates
        || state.usesMovingVolumeRates
            != state.surfaceLoads.usesMovingVolumeRates
        || state.timeStepSeconds != state.pressure.timeStepSeconds
        || state.timeStepSeconds != state.surfaceLoads.timeStepSeconds
        || state.densityKgPerCubicMeter
            != state.velocity.densityKgPerCubicMeter
        || state.fluidMomentumKilogramMetersPerSecond
            != state.velocity.momentumKilogramMetersPerSecond
        || state.fluidKineticEnergyJoules
            != state.velocity.kineticEnergyJoules
        || state.pressureForceOnSheetNewtons
            != state.surfaceLoads.totalPressureForceOnSheetNewtons
        || state.pressureImpulseOnSheetNewtonSeconds
            != state.surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds
        || state.pressureWorkToSheetJoules
            != state.surfaceLoads.totalPressureWorkToSheetJoules
        || !std::isfinite(state.timeStepSeconds)
        || !(state.timeStepSeconds > 0.0)
        || !std::isfinite(state.densityKgPerCubicMeter)
        || !(state.densityKgPerCubicMeter > 0.0)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "planar regional accepted-state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentAcceptedState(
    const PlanarPressureRegionFragmentAcceptedState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits) {
    validateLimits(limits);
    if (state != buildState(
                     grid, sweep, fragments, topology, metric,
                     state.velocity, state.pressure, state.surfaceLoads,
                     limits)) {
        throw std::invalid_argument(
            "planar regional accepted state is corrupted");
    }
}

} // namespace simwing::fsi::fluid
