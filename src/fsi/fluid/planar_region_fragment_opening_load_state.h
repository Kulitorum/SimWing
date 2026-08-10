#pragma once

#include "fluid/planar_region_fragment_opening_pressure_state.h"
#include "fluid/planar_region_fragment_opening_surface_load.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningLoadStateVersion = 1;

struct PlanarPressureRegionFragmentOpeningLoadStateLimits {
    PlanarPressureRegionFragmentOpeningPressureStateLimits pressureStateLimits;
    PlanarPressureRegionFragmentOpeningSurfaceLoadLimits surfaceLoadLimits;
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
};

// Atomic immutable owner for one accepted active-aperture pressure/load
// endpoint. It recursively retains the aperture-flow continuation, composed
// opening-connected pressure, full-wall load source, and exact opening/solid
// area partition. Aggregate force, impulse, origin moment, and material work
// make the retained-solid handoff visible without weakening nested provenance.
//
// Capture applies no Structure load, advances no fluid state, and owns no
// topology rebase, scene mapping, pressure relaxation, or worker selection.
struct PlanarPressureRegionFragmentOpeningLoadState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningLoadStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourcePressureStateFingerprint = 0;
    std::uint64_t sourceSurfaceLoadFingerprint = 0;
    std::uint64_t sourceOpeningSurfaceLoadFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    double fullWallAreaSquareMeters = 0.0;
    double openingAreaSquareMeters = 0.0;
    double solidAreaSquareMeters = 0.0;
    Vector3 fullWallPressureForceOnSheetNewtons;
    Vector3 openingRemovedPressureForceOnSheetNewtons;
    Vector3 solidPressureForceOnSheetNewtons;
    Vector3 fullWallPressureImpulseOnSheetNewtonSeconds;
    Vector3 openingRemovedPressureImpulseOnSheetNewtonSeconds;
    Vector3 solidPressureImpulseOnSheetNewtonSeconds;
    Vector3 fullWallPressureMomentOnSheetNewtonMeters;
    Vector3 openingRemovedPressureMomentOnSheetNewtonMeters;
    Vector3 solidPressureMomentOnSheetNewtonMeters;
    double fullWallPressureWorkToSheetJoules = 0.0;
    double openingRemovedPressureWorkToSheetJoules = 0.0;
    double solidPressureWorkToSheetJoules = 0.0;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedFlow;
    PlanarPressureRegionFragmentOpeningPressureState pressure;
    PlanarPressureRegionFragmentSurfaceLoadLedger surfaceLoads;
    PlanarPressureRegionFragmentOpeningSurfaceLoadLedger openingSurfaceLoads;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningLoadState&) const = default;
};

// Deterministically composes accepted opening flow through connected-gauge
// pressure, full-wall loads, and retained-solid aperture partition before
// capturing the same immutable load state exposed by the lower-level overload.
[[nodiscard]] PlanarPressureRegionFragmentOpeningLoadState
composePlanarPressureRegionFragmentOpeningLoadState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
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
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings =
        {},
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits = {});

[[nodiscard]] PlanarPressureRegionFragmentOpeningLoadState
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
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
    const PlanarPressureRegionFragmentOpeningLoadState& state);

void validatePlanarPressureRegionFragmentOpeningLoadState(
    const PlanarPressureRegionFragmentOpeningLoadState& state,
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
    const PlanarPressureRegionFragmentOpeningLoadStateLimits& limits = {});

} // namespace simwing::fsi::fluid
