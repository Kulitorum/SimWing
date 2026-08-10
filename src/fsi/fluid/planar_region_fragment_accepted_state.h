#pragma once

#include "fluid/planar_region_fragment_surface_load.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentAcceptedStateVersion = 1;

struct PlanarPressureRegionFragmentAcceptedStateLimits {
    PlanarPressureRegionFragmentVelocityStateLimits velocityStateLimits;
    PlanarPressureRegionFragmentSurfaceLoadLimits surfaceLoadLimits;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
};

// Rollback-safe accepted endpoint for one regional pressure-projection epoch.
// It atomically retains the projected diagonal velocity, composed total
// pressure, and provenance-rich sheet-load ledger. Every nested product is
// integrity checked and must share exact fragment/topology/metric/after-state
// identity before capture.
//
// The endpoint owns continuation data only. It performs no momentum transport,
// topology rebase, structural load application, or production worker commit.
struct PlanarPressureRegionFragmentAcceptedState {
    std::uint32_t version =
        planarPressureRegionFragmentAcceptedStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceVelocityStateFingerprint = 0;
    std::uint64_t sourcePressureStateFingerprint = 0;
    std::uint64_t sourceSurfaceLoadFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    Vector3 fluidMomentumKilogramMetersPerSecond;
    double fluidKineticEnergyJoules = 0.0;
    Vector3 pressureForceOnSheetNewtons;
    Vector3 pressureImpulseOnSheetNewtonSeconds;
    double pressureWorkToSheetJoules = 0.0;
    PlanarPressureRegionFragmentVelocityState velocity;
    PlanarPressureRegionFragmentPressureState pressure;
    PlanarPressureRegionFragmentSurfaceLoadLedger surfaceLoads;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentAcceptedState&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentAcceptedState
capturePlanarPressureRegionFragmentAcceptedState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocity,
    const PlanarPressureRegionFragmentPressureState& pressure,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits = {});

void validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
    const PlanarPressureRegionFragmentAcceptedState& state);

void validatePlanarPressureRegionFragmentAcceptedState(
    const PlanarPressureRegionFragmentAcceptedState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentAcceptedStateLimits& limits = {});

} // namespace simwing::fsi::fluid
