#include "scene_fluid_regional_opening_momentum_wall_cycle_owner.h"

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {

bool SceneFluidRegionalOpeningMomentumWallCycleOwner::hasState()
    const noexcept {
    return state_.has_value();
}

const SceneFluidRegionalOpeningMomentumWallCycleState&
SceneFluidRegionalOpeningMomentumWallCycleOwner::state() const {
    if (!state_) {
        throw std::logic_error(
            "regional opening wall-cycle owner has no accepted state");
    }
    return *state_;
}

SceneFluidRegionalOpeningMomentumWallCycleState
SceneFluidRegionalOpeningMomentumWallCycleOwner::checkpoint() const {
    const auto& current = state();
    validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(current);
    return current;
}

bool SceneFluidRegionalOpeningMomentumWallCycleOwner::tryCommit(
    const SceneFluidRegionalOpeningMomentumWallPressureEpoch& pressureEpoch,
    const SceneFluidRegionalOpeningMomentumWallExchange& wallExchange,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
        pressureEpoch);
    if (!pressureEpoch.pressureEpoch.diagnostics.accepted) return false;

    std::optional<SceneFluidRegionalOpeningMomentumWallCycleState> candidate{
        captureSceneFluidRegionalOpeningMomentumWallCycleState(
            pressureEpoch, wallExchange, transportMetric, acceptedMetric,
            quadrature, limits)};
    static_assert(std::is_nothrow_move_constructible_v<
        SceneFluidRegionalOpeningMomentumWallCycleState>);
    static_assert(std::is_nothrow_swappable_v<
        SceneFluidRegionalOpeningMomentumWallCycleState>);
    state_.swap(candidate);
    return true;
}

void SceneFluidRegionalOpeningMomentumWallCycleOwner::restore(
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
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    std::optional<SceneFluidRegionalOpeningMomentumWallCycleState> candidate{
        state};
    validateSceneFluidRegionalOpeningMomentumWallCycleState(
        *candidate, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric, acceptedMetric,
        quadrature, limits);
    state_.swap(candidate);
}

} // namespace simwing::fsi
