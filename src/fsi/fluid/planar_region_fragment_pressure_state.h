#pragma once

#include "fluid/planar_region_fragment_pressure_jump_energy.h"
#include "fluid/planar_region_fragment_projection_energy.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentPressureStateVersion = 1;

struct PlanarPressureRegionFragmentPressureStateControl {
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::size_t componentIndex = 0;
    std::uint64_t regionStableId = 0;
    double volumeCubicMeters = 0.0;
    double authoredPressurePascals = 0.0;
    double correctionPressurePascals = 0.0;
    double totalPressurePascals = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureStateControl&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureStateWall {
    std::size_t wallIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t minusFragmentIndex = 0;
    std::size_t plusFragmentIndex = 0;
    std::size_t minusComponentIndex = 0;
    std::size_t plusComponentIndex = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    double areaSquareMeters = 0.0;
    Vector3 unitNormalMinusToPlus;
    double minusTotalPressurePascals = 0.0;
    double plusTotalPressurePascals = 0.0;
    double authoredPressureJumpPascals = 0.0;
    double correctionPressureJumpPascals = 0.0;
    double totalPressureJumpPascals = 0.0;
    double pressureSplitResidualPascals = 0.0;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 pressureForceSplitResidualNewtons;
    double materialWallVelocityMetersPerSecond = 0.0;
    double authoredPressureWorkToFluidJoules = 0.0;
    double correctionPressureWorkToFluidJoules = 0.0;
    double totalPressureWorkToFluidJoules = 0.0;
    double pressureWorkSplitResidualJoules = 0.0;
    double totalPressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureStateWall&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureStateComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    double volumeCubicMeters = 0.0;
    double authoredVolumeMeanPressurePascals = 0.0;
    double correctionVolumeMeanPressurePascals = 0.0;
    double totalVolumeMeanPressurePascals = 0.0;
    double pressureMeanSplitResidualPascals = 0.0;
    double authoredGeometryPressureWorkJoules = 0.0;
    double correctionGeometryPressureWorkJoules = 0.0;
    double totalGeometryPressureWorkJoules = 0.0;
    double totalWallPressureWorkToFluidJoules = 0.0;
    double wallGeometryWorkResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureStateComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureStateLimits {
    PlanarPressureRegionFragmentProjectionEnergyLimits projectionEnergyLimits;
    PlanarPressureRegionFragmentPressureJumpEnergyLimits pressureJumpLimits;
    std::size_t maximumControls = 20'000'000;
    std::size_t maximumWalls = 60'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Immutable composition of the authored regional pressure and one accepted
// correction-pressure audit. Every wall publishes the corresponding authored,
// correction and total sheet force without applying it. Moving composition
// also proves that the total material-wall work equals the sum of the two
// independently accepted geometry-work ledgers. The source jump audit must be
// bound to the projection's after-state and both audits must describe the same
// epoch and time step.
//
// This state is an opt-in diagnostic endpoint. It owns no momentum update,
// structural load transfer, transport, topology rebase, or production state.
struct PlanarPressureRegionFragmentPressureState {
    std::uint32_t version =
        planarPressureRegionFragmentPressureStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t sourceProjectionEnergyFingerprint = 0;
    std::uint64_t sourcePressureJumpEnergyFingerprint = 0;
    std::uint64_t volumeRateFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    std::vector<PlanarPressureRegionFragmentPressureStateControl> controls;
    std::vector<PlanarPressureRegionFragmentPressureStateWall> walls;
    std::vector<PlanarPressureRegionFragmentPressureStateComponent>
        components;
    double maximumAbsoluteCorrectionGaugePascals = 0.0;
    double maximumAbsolutePressureSplitResidualPascals = 0.0;
    double maximumAbsoluteForceSplitResidualNewtons = 0.0;
    double maximumAbsoluteWorkResidualJoules = 0.0;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 pressureForceSplitResidualNewtons;
    double authoredPressureWorkToFluidJoules = 0.0;
    double correctionPressureWorkToFluidJoules = 0.0;
    double totalPressureWorkToFluidJoules = 0.0;
    double pressureWorkSplitResidualJoules = 0.0;
    double totalPressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;
    double totalGeometryPressureWorkJoules = 0.0;
    double wallGeometryWorkResidualJoules = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureState&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentPressureState
composeStaticPlanarPressureRegionFragmentPressureState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits = {});

void validateStaticPlanarPressureRegionFragmentPressureState(
    const PlanarPressureRegionFragmentPressureState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits = {});

[[nodiscard]] PlanarPressureRegionFragmentPressureState
composeMovingPlanarPressureRegionFragmentPressureState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits = {});

void validateMovingPlanarPressureRegionFragmentPressureState(
    const PlanarPressureRegionFragmentPressureState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits = {});

} // namespace simwing::fsi::fluid
