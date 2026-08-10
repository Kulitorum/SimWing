#pragma once

#include "fluid/planar_region_fragment_opening_momentum_cycle.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumCycleStateVersion = 1;

struct PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits {
    PlanarPressureRegionFragmentOpeningMomentumTransportLimits transport;
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits acceptedState;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// Immutable restart endpoint for the accepted staggered opening-momentum
// cycle. Transported collocated momentum belongs to transportMetric, while the
// accepted pressure endpoint belongs to acceptedMetric. The intervening
// prediction, pressure-only warm start, and predicted aperture flux bind the
// two endpoints to one successfully committed cycle.
//
// Capture accepts only a complete cycle result. Full validation additionally
// binds the transport controls to their live current volume-rate epoch and the
// pressure endpoint to trusted next-geometry sources. No rejected attempt,
// temporary predictor, worker state, or topology-rebase state is retained.
struct PlanarPressureRegionFragmentOpeningMomentumCycleState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumCycleStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t transportMetricFingerprint = 0;
    std::uint64_t acceptedMetricFingerprint = 0;
    std::uint64_t currentAcceptedFlowStateFingerprint = 0;
    std::uint64_t predictionFingerprint = 0;
    std::uint64_t pressureWarmStartFingerprint = 0;
    std::uint64_t predictedOpeningFluxFingerprint = 0;
    PlanarPressureRegionFragmentOpeningMomentumTransport transport;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedState;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumCycleState&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningMomentumCycleState
capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state);

void validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    const PlanarPressureRegionFragmentVolumeRateSet& transportVolumeRates,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& acceptedSweep,
    const PlanarPressureRegionFragmentSet& acceptedFragments,
    const PlanarPressureRegionFragmentTopology& acceptedTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits =
        {});

} // namespace simwing::fsi::fluid
