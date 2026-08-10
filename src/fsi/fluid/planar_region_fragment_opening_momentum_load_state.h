#pragma once

#include "fluid/planar_region_fragment_opening_load_state.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentOpeningMomentumLoadStateLimits {
    PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits cycleState;
    PlanarPressureRegionFragmentOpeningLoadStateLimits loadState;
};

// Structure-facing pressure-load composition for one accepted transported
// opening cycle. Full validation first binds collocated transport to its exact
// current metric/volume epoch and accepted pressure to the consecutive next
// epoch. The established connected-pressure/full-wall/retained-solid pipeline
// then produces the immutable load state.
//
// The source cycle state remains immutable. This adapter applies no Structure
// load, adds no viscosity or wall shear, performs no topology rebase, and does
// not select a production worker path.
[[nodiscard]] PlanarPressureRegionFragmentOpeningLoadState
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
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings =
        {},
    const PlanarPressureRegionFragmentOpeningMomentumLoadStateLimits& limits =
        {});

} // namespace simwing::fsi::fluid
