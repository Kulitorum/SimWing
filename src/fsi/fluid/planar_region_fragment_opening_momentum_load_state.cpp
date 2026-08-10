#include "fluid/planar_region_fragment_opening_momentum_load_state.h"

namespace simwing::fsi::fluid {

PlanarPressureRegionFragmentOpeningLoadState
composePlanarPressureRegionFragmentOpeningMomentumLoadState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& cycleState,
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
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings,
    const PlanarPressureRegionFragmentOpeningMomentumLoadStateLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
        cycleState, transportVolumeRates, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedMetric, limits.cycleState);
    return composePlanarPressureRegionFragmentOpeningLoadState(
        cycleState.acceptedState, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, settings, limits.loadState);
}

} // namespace simwing::fsi::fluid
