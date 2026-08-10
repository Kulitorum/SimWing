#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint16_t
    planarPressureRegionFragmentOpeningAcceptedStateProtocolVersion = 1;

struct PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits {
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits stateLimits;
    std::size_t maximumEncodedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumTopologyLinkVelocities = 140'000'000;
    std::size_t maximumOpeningSamples = 20'000'000;
    std::size_t maximumPressureCorrections = 20'000'000;
};

enum class
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode {
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

struct PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError {
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode code =
        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
            None;
    std::string message;
};

// Deterministic bounded little-endian persistence for one accepted regional
// active-aperture endpoint. The wire owns the primary link/aperture/pressure
// fields and acceptance certificates, but not redundant opening-flux rows.
// Decode rebuilds that ledger from trusted current sources, then publishes
// only after the original state fingerprint and full source validation pass.
[[nodiscard]] bool
serializePlanarPressureRegionFragmentOpeningAcceptedState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
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
    std::vector<std::uint8_t>& bytes,
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError* error =
        nullptr,
    const PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits&
        limits = {});

[[nodiscard]] bool
deserializePlanarPressureRegionFragmentOpeningAcceptedState(
    std::span<const std::uint8_t> bytes,
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
    PlanarPressureRegionFragmentOpeningAcceptedState& state,
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError* error =
        nullptr,
    const PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits&
        limits = {});

} // namespace simwing::fsi::fluid
