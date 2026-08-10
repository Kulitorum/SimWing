#include "fluid/planar_region_fragment_opening_momentum_cycle_owner.h"

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi::fluid {

bool PlanarPressureRegionFragmentOpeningMomentumCycleOwner::hasState()
    const noexcept {
    return state_.has_value();
}

const PlanarPressureRegionFragmentOpeningMomentumCycleState&
PlanarPressureRegionFragmentOpeningMomentumCycleOwner::state() const {
    if (!state_) {
        throw std::logic_error(
            "opening momentum-cycle owner has no accepted state");
    }
    return *state_;
}

PlanarPressureRegionFragmentOpeningMomentumCycleState
PlanarPressureRegionFragmentOpeningMomentumCycleOwner::checkpoint() const {
    const auto& current = state();
    validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
        current);
    return current;
}

bool PlanarPressureRegionFragmentOpeningMomentumCycleOwner::tryCommit(
    const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
        result);
    if (!result.diagnostics.accepted) return false;

    std::optional<PlanarPressureRegionFragmentOpeningMomentumCycleState>
        candidate{capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
            result, transportMetric, acceptedMetric, limits)};
    static_assert(std::is_nothrow_move_constructible_v<
        PlanarPressureRegionFragmentOpeningMomentumCycleState>);
    static_assert(std::is_nothrow_swappable_v<
        PlanarPressureRegionFragmentOpeningMomentumCycleState>);
    state_.swap(candidate);
    return true;
}

void PlanarPressureRegionFragmentOpeningMomentumCycleOwner::restore(
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
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits& limits) {
    std::optional<PlanarPressureRegionFragmentOpeningMomentumCycleState>
        candidate{state};
    validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
        *candidate, transportVolumeRates, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedMetric, limits);
    state_.swap(candidate);
}

} // namespace simwing::fsi::fluid
