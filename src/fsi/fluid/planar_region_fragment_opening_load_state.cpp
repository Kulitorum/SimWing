#include "fluid/planar_region_fragment_opening_load_state.h"

#include <bit>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned> void integer(Unsigned value) {
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
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "opening load-state limits are invalid");
    }
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error("opening load-state storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger&
        openingSurfaceLoads) {
    return checkedAdd(
        checkedAdd(
            acceptedState.ownedStorageBytes, pressureState.ownedStorageBytes),
        checkedAdd(
            surfaceLoads.ownedStorageBytes,
            openingSurfaceLoads.ownedStorageBytes));
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningLoadState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    for (const std::uint64_t value : {
             state.sourceAcceptedStateFingerprint,
             state.sourcePressureStateFingerprint,
             state.sourceSurfaceLoadFingerprint,
             state.sourceOpeningSurfaceLoadFingerprint,
             state.sourceOpeningFingerprint,
             state.sourceFragmentFingerprint,
             state.sourceTopologyFingerprint,
             state.sourceVolumeRateFingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.boolean(state.staticGeometry);
    fingerprint.boolean(state.usesMovingVolumeRates);
    for (const double value : {
             state.timeStepSeconds,
             state.fullWallAreaSquareMeters,
             state.openingAreaSquareMeters,
             state.solidAreaSquareMeters}) {
        fingerprint.real(value);
    }
    for (const Vector3& value : {
             state.fullWallPressureForceOnSheetNewtons,
             state.openingRemovedPressureForceOnSheetNewtons,
             state.solidPressureForceOnSheetNewtons,
             state.fullWallPressureImpulseOnSheetNewtonSeconds,
             state.openingRemovedPressureImpulseOnSheetNewtonSeconds,
             state.solidPressureImpulseOnSheetNewtonSeconds,
             state.fullWallPressureMomentOnSheetNewtonMeters,
             state.openingRemovedPressureMomentOnSheetNewtonMeters,
             state.solidPressureMomentOnSheetNewtonMeters}) {
        fingerprintVector(fingerprint, value);
    }
    fingerprint.real(state.fullWallPressureWorkToSheetJoules);
    fingerprint.real(state.openingRemovedPressureWorkToSheetJoules);
    fingerprint.real(state.solidPressureWorkToSheetJoules);
    fingerprint.integer(state.acceptedFlow.fingerprint);
    fingerprint.integer(state.pressure.fingerprint);
    fingerprint.integer(state.surfaceLoads.fingerprint);
    fingerprint.integer(state.openingSurfaceLoads.fingerprint);
    fingerprint.boolean(state.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentOpeningLoadState buildState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger&
        openingSurfaceLoads,
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
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions,
        limits.pressureStateLimits.acceptedStateLimits);
    validatePlanarPressureRegionFragmentOpeningPressureState(
        pressureState, acceptedState, pressureOperator, basePressureOperator,
        grid, sweep, fragments, topology, volumeRates, openingDefinitions,
        openings, resistanceDefinitions, limits.pressureStateLimits);
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState,
        limits.surfaceLoadLimits.surfaceLoadLimits);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        openingSurfaceLoads, surfaceLoads, pressureState, grid, sweep,
        fragments, topology, openingDefinitions, openings,
        limits.surfaceLoadLimits);

    if (pressureState.sourceAcceptedStateFingerprint
            != acceptedState.fingerprint
        || pressureState.sourceOpeningFingerprint != openings.fingerprint
        || pressureState.sourceFragmentFingerprint != fragments.fingerprint
        || pressureState.sourceTopologyFingerprint != topology.fingerprint
        || pressureState.sourceVolumeRateFingerprint != volumeRates.fingerprint
        || surfaceLoads.sourcePressureStateFingerprint
            != pressureState.fingerprint
        || surfaceLoads.sourceTopologyFingerprint != topology.fingerprint
        || openingSurfaceLoads.sourceSurfaceLoadFingerprint
            != surfaceLoads.fingerprint
        || openingSurfaceLoads.sourcePressureStateFingerprint
            != pressureState.fingerprint
        || openingSurfaceLoads.sourceOpeningFingerprint
            != openings.fingerprint
        || openingSurfaceLoads.sourceTopologyFingerprint
            != topology.fingerprint
        || pressureState.staticGeometry != surfaceLoads.staticGeometry
        || pressureState.staticGeometry
            != openingSurfaceLoads.staticGeometry
        || pressureState.usesMovingVolumeRates
            != surfaceLoads.usesMovingVolumeRates
        || pressureState.usesMovingVolumeRates
            != openingSurfaceLoads.usesMovingVolumeRates
        || pressureState.timeStepSeconds != surfaceLoads.timeStepSeconds
        || pressureState.timeStepSeconds
            != openingSurfaceLoads.timeStepSeconds) {
        throw std::invalid_argument(
            "opening load-state sources are incompatible");
    }

    PlanarPressureRegionFragmentOpeningLoadState result;
    result.sourceAcceptedStateFingerprint = acceptedState.fingerprint;
    result.sourcePressureStateFingerprint = pressureState.fingerprint;
    result.sourceSurfaceLoadFingerprint = surfaceLoads.fingerprint;
    result.sourceOpeningSurfaceLoadFingerprint =
        openingSurfaceLoads.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceVolumeRateFingerprint = volumeRates.fingerprint;
    result.staticGeometry = pressureState.staticGeometry;
    result.usesMovingVolumeRates = pressureState.usesMovingVolumeRates;
    result.timeStepSeconds = pressureState.timeStepSeconds;
    result.fullWallAreaSquareMeters =
        openingSurfaceLoads.totalWallAreaSquareMeters;
    result.openingAreaSquareMeters =
        openingSurfaceLoads.totalOpeningAreaSquareMeters;
    result.solidAreaSquareMeters =
        openingSurfaceLoads.totalSolidAreaSquareMeters;
    result.fullWallPressureForceOnSheetNewtons =
        openingSurfaceLoads.sourceTotalPressureForceOnSheetNewtons;
    result.openingRemovedPressureForceOnSheetNewtons =
        openingSurfaceLoads.openingRemovedTotalPressureForceOnSheetNewtons;
    result.solidPressureForceOnSheetNewtons =
        openingSurfaceLoads.solidTotalPressureForceOnSheetNewtons;
    result.fullWallPressureImpulseOnSheetNewtonSeconds =
        openingSurfaceLoads.sourceTotalPressureImpulseOnSheetNewtonSeconds;
    result.openingRemovedPressureImpulseOnSheetNewtonSeconds =
        openingSurfaceLoads
            .openingRemovedTotalPressureImpulseOnSheetNewtonSeconds;
    result.solidPressureImpulseOnSheetNewtonSeconds =
        openingSurfaceLoads.solidTotalPressureImpulseOnSheetNewtonSeconds;
    result.fullWallPressureMomentOnSheetNewtonMeters =
        openingSurfaceLoads.sourceTotalPressureMomentOnSheetNewtonMeters;
    result.openingRemovedPressureMomentOnSheetNewtonMeters =
        openingSurfaceLoads
            .openingRemovedTotalPressureMomentOnSheetNewtonMeters;
    result.solidPressureMomentOnSheetNewtonMeters =
        openingSurfaceLoads.solidTotalPressureMomentOnSheetNewtonMeters;
    result.fullWallPressureWorkToSheetJoules =
        openingSurfaceLoads.sourceTotalPressureWorkToSheetJoules;
    result.openingRemovedPressureWorkToSheetJoules =
        openingSurfaceLoads.openingRemovedTotalPressureWorkToSheetJoules;
    result.solidPressureWorkToSheetJoules =
        openingSurfaceLoads.solidTotalPressureWorkToSheetJoules;
    result.ownedStorageBytes = ownedStorageBytes(
        acceptedState, pressureState, surfaceLoads, openingSurfaceLoads);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening load-state storage limit exceeded");
    }
    result.acceptedFlow = acceptedState;
    result.pressure = pressureState;
    result.surfaceLoads = surfaceLoads;
    result.openingSurfaceLoads = openingSurfaceLoads;
    result.accepted = true;
    result.fingerprint = stateFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningLoadState
capturePlanarPressureRegionFragmentOpeningLoadState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger&
        openingSurfaceLoads,
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
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits) {
    return buildState(
        acceptedState, pressureState, surfaceLoads, openingSurfaceLoads,
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits);
}

void validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
    const PlanarPressureRegionFragmentOpeningLoadState& state) {
    validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
        state.acceptedFlow);
    validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(
        state.pressure);
    validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
        state.surfaceLoads);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
        state.openingSurfaceLoads);
    if (state.version
            != planarPressureRegionFragmentOpeningLoadStateVersion
        || state.fingerprint == 0 || !state.accepted
        || state.sourceAcceptedStateFingerprint
            != state.acceptedFlow.fingerprint
        || state.sourcePressureStateFingerprint != state.pressure.fingerprint
        || state.sourceSurfaceLoadFingerprint
            != state.surfaceLoads.fingerprint
        || state.sourceOpeningSurfaceLoadFingerprint
            != state.openingSurfaceLoads.fingerprint
        || state.sourceOpeningFingerprint
            != state.acceptedFlow.sourceOpeningFingerprint
        || state.sourceOpeningFingerprint
            != state.pressure.sourceOpeningFingerprint
        || state.sourceOpeningFingerprint
            != state.openingSurfaceLoads.sourceOpeningFingerprint
        || state.sourceFragmentFingerprint
            != state.acceptedFlow.sourceFragmentFingerprint
        || state.sourceFragmentFingerprint
            != state.pressure.sourceFragmentFingerprint
        || state.sourceTopologyFingerprint
            != state.acceptedFlow.sourceTopologyFingerprint
        || state.sourceTopologyFingerprint
            != state.pressure.sourceTopologyFingerprint
        || state.sourceTopologyFingerprint
            != state.surfaceLoads.sourceTopologyFingerprint
        || state.sourceTopologyFingerprint
            != state.openingSurfaceLoads.sourceTopologyFingerprint
        || state.sourceVolumeRateFingerprint
            != state.acceptedFlow.sourceVolumeRateFingerprint
        || state.sourceVolumeRateFingerprint
            != state.pressure.sourceVolumeRateFingerprint
        || state.pressure.sourceAcceptedStateFingerprint
            != state.acceptedFlow.fingerprint
        || state.surfaceLoads.sourcePressureStateFingerprint
            != state.pressure.fingerprint
        || state.openingSurfaceLoads.sourceSurfaceLoadFingerprint
            != state.surfaceLoads.fingerprint
        || state.openingSurfaceLoads.sourcePressureStateFingerprint
            != state.pressure.fingerprint
        || state.staticGeometry != state.pressure.staticGeometry
        || state.staticGeometry != state.surfaceLoads.staticGeometry
        || state.staticGeometry != state.openingSurfaceLoads.staticGeometry
        || state.usesMovingVolumeRates
            != state.pressure.usesMovingVolumeRates
        || state.usesMovingVolumeRates
            != state.surfaceLoads.usesMovingVolumeRates
        || state.usesMovingVolumeRates
            != state.openingSurfaceLoads.usesMovingVolumeRates
        || state.timeStepSeconds != state.pressure.timeStepSeconds
        || state.timeStepSeconds != state.surfaceLoads.timeStepSeconds
        || state.timeStepSeconds
            != state.openingSurfaceLoads.timeStepSeconds
        || state.fullWallAreaSquareMeters
            != state.openingSurfaceLoads.totalWallAreaSquareMeters
        || state.openingAreaSquareMeters
            != state.openingSurfaceLoads.totalOpeningAreaSquareMeters
        || state.solidAreaSquareMeters
            != state.openingSurfaceLoads.totalSolidAreaSquareMeters
        || state.fullWallPressureForceOnSheetNewtons
            != state.openingSurfaceLoads
                   .sourceTotalPressureForceOnSheetNewtons
        || state.openingRemovedPressureForceOnSheetNewtons
            != state.openingSurfaceLoads
                   .openingRemovedTotalPressureForceOnSheetNewtons
        || state.solidPressureForceOnSheetNewtons
            != state.openingSurfaceLoads.solidTotalPressureForceOnSheetNewtons
        || state.fullWallPressureImpulseOnSheetNewtonSeconds
            != state.openingSurfaceLoads
                   .sourceTotalPressureImpulseOnSheetNewtonSeconds
        || state.openingRemovedPressureImpulseOnSheetNewtonSeconds
            != state.openingSurfaceLoads
                   .openingRemovedTotalPressureImpulseOnSheetNewtonSeconds
        || state.solidPressureImpulseOnSheetNewtonSeconds
            != state.openingSurfaceLoads
                   .solidTotalPressureImpulseOnSheetNewtonSeconds
        || state.fullWallPressureMomentOnSheetNewtonMeters
            != state.openingSurfaceLoads
                   .sourceTotalPressureMomentOnSheetNewtonMeters
        || state.openingRemovedPressureMomentOnSheetNewtonMeters
            != state.openingSurfaceLoads
                   .openingRemovedTotalPressureMomentOnSheetNewtonMeters
        || state.solidPressureMomentOnSheetNewtonMeters
            != state.openingSurfaceLoads
                   .solidTotalPressureMomentOnSheetNewtonMeters
        || state.fullWallPressureWorkToSheetJoules
            != state.openingSurfaceLoads.sourceTotalPressureWorkToSheetJoules
        || state.openingRemovedPressureWorkToSheetJoules
            != state.openingSurfaceLoads
                   .openingRemovedTotalPressureWorkToSheetJoules
        || state.solidPressureWorkToSheetJoules
            != state.openingSurfaceLoads.solidTotalPressureWorkToSheetJoules
        || !std::isfinite(state.timeStepSeconds)
        || !(state.timeStepSeconds > 0.0)
        || !std::isfinite(state.fullWallAreaSquareMeters)
        || !(state.fullWallAreaSquareMeters > 0.0)
        || !std::isfinite(state.openingAreaSquareMeters)
        || !(state.openingAreaSquareMeters > 0.0)
        || !std::isfinite(state.solidAreaSquareMeters)
        || state.solidAreaSquareMeters < 0.0
        || !finiteVector(state.fullWallPressureForceOnSheetNewtons)
        || !finiteVector(state.openingRemovedPressureForceOnSheetNewtons)
        || !finiteVector(state.solidPressureForceOnSheetNewtons)
        || !finiteVector(state.fullWallPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(
            state.openingRemovedPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(state.solidPressureImpulseOnSheetNewtonSeconds)
        || !finiteVector(state.fullWallPressureMomentOnSheetNewtonMeters)
        || !finiteVector(
            state.openingRemovedPressureMomentOnSheetNewtonMeters)
        || !finiteVector(state.solidPressureMomentOnSheetNewtonMeters)
        || !std::isfinite(state.fullWallPressureWorkToSheetJoules)
        || !std::isfinite(state.openingRemovedPressureWorkToSheetJoules)
        || !std::isfinite(state.solidPressureWorkToSheetJoules)
        || state.ownedStorageBytes
            != ownedStorageBytes(
                state.acceptedFlow, state.pressure, state.surfaceLoads,
                state.openingSurfaceLoads)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening load-state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningLoadState(
    const PlanarPressureRegionFragmentOpeningLoadState& state,
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
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits) {
    validateLimits(limits);
    if (state != buildState(
                     state.acceptedFlow, state.pressure, state.surfaceLoads,
                     state.openingSurfaceLoads, pressureOperator,
                     basePressureOperator, grid, sweep, fragments, topology,
                     volumeRates, openingDefinitions, openings,
                     resistanceDefinitions, limits)) {
        throw std::invalid_argument("opening load state is corrupted");
    }
}

} // namespace simwing::fsi::fluid
