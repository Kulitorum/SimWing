#pragma once

#include "scene_fluid_regional_opening_momentum_wall_cycle_state.h"

#include <optional>

namespace simwing::fsi {

// Transactional in-memory owner for the opt-in adjusted-momentum/accepted-
// pressure/wall-traction endpoint. Rejected pressure receipts retain the
// prior state. Accepted commit and restore validate a complete candidate
// before one no-throw publication swap.
//
// The owner does not run wall exchange, transport, pressure, Structure, byte
// persistence, retries, or production selection.
class SceneFluidRegionalOpeningMomentumWallCycleOwner final {
public:
    SceneFluidRegionalOpeningMomentumWallCycleOwner() = default;

    [[nodiscard]] bool hasState() const noexcept;
    [[nodiscard]] const SceneFluidRegionalOpeningMomentumWallCycleState&
    state() const;
    [[nodiscard]] SceneFluidRegionalOpeningMomentumWallCycleState
    checkpoint() const;

    [[nodiscard]] bool tryCommit(
        const SceneFluidRegionalOpeningMomentumWallPressureEpoch& pressureEpoch,
        const SceneFluidRegionalOpeningMomentumWallExchange& wallExchange,
        const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
            transportMetric,
        const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
            acceptedMetric,
        const SceneFluidQuadratureDefinition& quadrature,
        const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits =
            {});

    void restore(
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
        const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits =
            {});

private:
    std::optional<SceneFluidRegionalOpeningMomentumWallCycleState> state_;
};

} // namespace simwing::fsi
