#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state_persistence.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint16_t
    planarPressureRegionFragmentOpeningMomentumCycleStateProtocolVersion = 1;

struct PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits {
    PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits state;
    std::size_t maximumEncodedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumTransportControls = 20'000'000;
};

enum class
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
    SourceMismatch,
};

struct PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError {
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode
        code =
            PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                None;
    std::string message;
};

// Deterministic bounded little-endian SWRM persistence for one accepted
// staggered opening-momentum endpoint pair. The transport controls are stored
// directly. The accepted pressure endpoint remains an embedded SWRO artifact,
// so decode rebuilds its derived opening-flux rows through the existing trusted
// source validator. Publication is transactional.
[[nodiscard]] bool
serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
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
    std::vector<std::uint8_t>& bytes,
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError*
        error = nullptr,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits&
        limits = {});

[[nodiscard]] bool
deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
    std::span<const std::uint8_t> bytes,
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
    PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError*
        error = nullptr,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits&
        limits = {});

} // namespace simwing::fsi::fluid
