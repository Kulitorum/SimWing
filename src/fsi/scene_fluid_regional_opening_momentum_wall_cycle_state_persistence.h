#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state_persistence.h"
#include "scene_fluid_regional_opening_momentum_wall_cycle_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    sceneFluidRegionalOpeningMomentumWallCycleStateProtocolVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits {
    SceneFluidRegionalOpeningMomentumWallCycleStateLimits state;
    std::size_t maximumEncodedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumAdjustmentControls = 20'000'000;
    std::size_t maximumWallTractions = 100'000'000;
};

enum class
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode {
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

struct SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError {
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode code =
        SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode::
            None;
    std::string message;
};

// Deterministic bounded little-endian SWRW persistence for the opt-in
// adjusted-momentum/accepted-pressure/wall-traction endpoint. Adjusted
// collocated controls and conservative quadrature tractions are primary wire
// data. The pressure endpoint remains an embedded SWRO artifact whose derived
// opening-flux rows are rebuilt from trusted accepted-geometry sources.
// Encode and decode publish transactionally and do not modify SWRM.
[[nodiscard]] bool
serializeSceneFluidRegionalOpeningMomentumWallCycleState(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    std::vector<std::uint8_t>& bytes,
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError* error =
        nullptr,
    const SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits&
        limits = {});

[[nodiscard]] bool
deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
    std::span<const std::uint8_t> bytes,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    SceneFluidRegionalOpeningMomentumWallCycleState& state,
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError* error =
        nullptr,
    const SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits&
        limits = {});

} // namespace simwing::fsi
