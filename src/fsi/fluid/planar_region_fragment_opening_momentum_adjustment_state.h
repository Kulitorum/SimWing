#pragma once

#include "fluid/planar_region_fragment_opening_momentum_transport.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumAdjustmentStateVersion = 1;

struct PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings {
    double absoluteMomentumToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-11;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits {
    PlanarPressureRegionFragmentOpeningVelocityMetricLimits metricLimits;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Fluid-owned identity-preserving collocated momentum endpoint after one
// externally proven same-epoch adjustment. Stable fragment, region, component,
// volume, density, time-step, accepted transport, and transport-target metric
// ownership remain exact; only velocity and momentum may differ. The opaque
// nonzero adjustment fingerprint lets a higher-level adapter bind its own
// physical conservation and energy proof without introducing scene types into
// the fluid library.
//
// This state performs no pressure prediction, wall law, Structure mutation,
// topology rebase, persistence, or production-worker selection.
struct PlanarPressureRegionFragmentOpeningMomentumAdjustmentState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumAdjustmentStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceAdjustmentFingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings
        settings;
    Vector3 sourceMomentumKilogramMetersPerSecond;
    Vector3 adjustedMomentumKilogramMetersPerSecond;
    Vector3 adjustmentImpulseKilogramMetersPerSecond;
    double sourceKineticEnergyJoules = 0.0;
    double adjustedKineticEnergyJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double maximumMomentumVelocityResidualKilogramMetersPerSecond = 0.0;
    std::vector<PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        controls;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState&)
        const = default;
};

[[nodiscard]]
PlanarPressureRegionFragmentOpeningMomentumAdjustmentState
capturePlanarPressureRegionFragmentOpeningMomentumAdjustmentState(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    std::uint64_t sourceAdjustmentFingerprint,
    std::span<
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings = {},
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits = {});

void validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state);

void validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentState(
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    std::uint64_t sourceAdjustmentFingerprint,
    std::span<
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls,
    const PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits = {});

} // namespace simwing::fsi::fluid
