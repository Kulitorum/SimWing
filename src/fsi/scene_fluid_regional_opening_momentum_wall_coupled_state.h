#pragma once

#include "scene_fluid_regional_opening_momentum_wall_structure_step_epoch.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallCoupledStateVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallCoupledStateLimits {
    fluid::PlanarPressureRegionFragmentOpeningVelocityStateLimits
        consecutiveFlow;
    fluid::PlanarPressureRegionFragmentOpeningMomentumTransportLimits
        consecutiveTransport;
    SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits
        structureStep;
    std::size_t maximumOwnedBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One accepted opt-in coupled endpoint. The fluid side owns the adjusted
// collocated momentum, consecutive pressure, and conservative wall traction;
// the structural side owns the exact one-step receipt and its complete
// post-step Structure checkpoint encoding.
//
// This is an in-memory ownership boundary, not a new persistence protocol or
// a production-worker selection.
struct SceneFluidRegionalOpeningMomentumWallCoupledState {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallCoupledStateVersion;
    std::uint64_t fingerprint = 0;
    SceneFluidRegionalOpeningMomentumWallCycleState cycleState;
    SceneFluidRegionalOpeningMomentumWallStructureStepEpoch structureStep;
    std::size_t ownedStorageBytes = 0;
};

void validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state);

void validateSceneFluidRegionalOpeningMomentumWallCoupledState(
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
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits =
        {});

// Owns Structure and its matching accepted regional wall-cycle endpoint as one
// transaction. Bootstrap advance requires an empty owner; consecutive
// replacement requires explicit fixed-metric fluid lineage. Both leave the
// prior owner state untouched until the structural step has fully
// replay-validated. Restore validates and allocates a complete replacement
// Structure and state before publishing either.
class SceneFluidRegionalOpeningMomentumWallCoupledStateOwner final {
public:
    explicit SceneFluidRegionalOpeningMomentumWallCoupledStateOwner(
        Structure structure);

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner(
        const SceneFluidRegionalOpeningMomentumWallCoupledStateOwner&) =
        delete;
    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& operator=(
        const SceneFluidRegionalOpeningMomentumWallCoupledStateOwner&) =
        delete;
    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner(
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner&&) noexcept =
        default;
    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& operator=(
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner&&) noexcept =
        default;

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] bool hasState() const noexcept;
    [[nodiscard]] const SceneFluidRegionalOpeningMomentumWallCoupledState&
    state() const;
    [[nodiscard]] SceneFluidRegionalOpeningMomentumWallCoupledState
    checkpoint() const;

    void advance(
        const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
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
        const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
            settings,
        const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits =
            {});

    // Repeats the opt-in coupled transaction while the regional CFD metric
    // remains fixed. The supplied transport must replay exactly from this
    // owner's accepted adjusted momentum into its accepted pressure flow;
    // the new cycle must in turn descend from that transport. This is the
    // first explicit multi-step ownership boundary, not a moving-topology
    // rebase or a production-worker selection.
    void advanceFixedMetricConsecutive(
        const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
        const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
            consecutiveTransport,
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
        const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
            settings,
        const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits =
            {});

    void restore(
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
        const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
            settings,
        const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits =
            {});

private:
    Structure structure_;
    std::optional<SceneFluidRegionalOpeningMomentumWallCoupledState> state_;
};

} // namespace simwing::fsi
