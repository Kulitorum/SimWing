#pragma once

#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"

#include <optional>

namespace simwing::fsi::fluid {

// Opt-in transactional owner for one accepted staggered opening-momentum
// endpoint pair. A rejected numerical cycle returns false and retains the
// prior state. Invalid input throws before mutation. Accepted commit and
// restore construct and validate a complete candidate before one no-throw
// publication swap.
//
// This owner deliberately does not build geometry, choose bootstrap/re-entrant
// source policy, retry with modified settings, persist bytes, rebase topology,
// or select a production worker path.
class PlanarPressureRegionFragmentOpeningMomentumCycleOwner final {
public:
    PlanarPressureRegionFragmentOpeningMomentumCycleOwner() = default;

    [[nodiscard]] bool hasState() const noexcept;
    [[nodiscard]] const PlanarPressureRegionFragmentOpeningMomentumCycleState&
    state() const;
    [[nodiscard]] PlanarPressureRegionFragmentOpeningMomentumCycleState
    checkpoint() const;

    [[nodiscard]] bool tryCommit(
        const PlanarPressureRegionFragmentOpeningMomentumCycleResult& result,
        const PlanarPressureRegionFragmentOpeningVelocityMetric&
            transportMetric,
        const PlanarPressureRegionFragmentOpeningVelocityMetric&
            acceptedMetric,
        const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits&
            limits = {});

    void restore(
        const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
        const PlanarPressureRegionFragmentVolumeRateSet& transportVolumeRates,
        const PlanarPressureRegionFragmentOpeningVelocityMetric&
            transportMetric,
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
        std::span<
            const PlanarPressureRegionFragmentOpeningResistanceDefinition>
            acceptedResistanceDefinitions,
        const PlanarPressureRegionFragmentOpeningVelocityMetric&
            acceptedMetric,
        const PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits&
            limits = {});

private:
    std::optional<
        PlanarPressureRegionFragmentOpeningMomentumCycleState> state_;
};

} // namespace simwing::fsi::fluid
