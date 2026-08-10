#pragma once

#include "scene_fluid_regional_opening_momentum_wall_coupled_state.h"
#include "scene_fluid_regional_opening_momentum_wall_cycle_state_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    sceneFluidRegionalOpeningMomentumWallCoupledStateProtocolVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits {
    SceneFluidRegionalOpeningMomentumWallCoupledStateLimits state;
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits
        cycleState;
    std::size_t maximumEncodedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

enum class
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
    SourceMismatch,
    TopologyMismatch,
    ReplayMismatch,
};

struct SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError {
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode code =
        SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::
            None;
    std::string message;
};

// Deterministic bounded little-endian SWRC persistence for one accepted fluid/
// Structure pair. The compact fluid endpoint remains a nested SWRW artifact.
// Only the existing complete pre/post Structure checkpoint encodings are
// retained from the structural receipt. Decode restores the pre-step state and
// deterministically replays the one load/XPBD epoch to rebuild the full
// in-memory receipt before publishing the output.
[[nodiscard]] bool
serializeSceneFluidRegionalOpeningMomentumWallCoupledState(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& structureDefinitionOwner,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    std::vector<std::uint8_t>& bytes,
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError* error =
        nullptr,
    const SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits&
        limits = {});

[[nodiscard]] bool
deserializeSceneFluidRegionalOpeningMomentumWallCoupledState(
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& structureDefinitionOwner,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError* error =
        nullptr,
    const SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits&
        limits = {});

} // namespace simwing::fsi
