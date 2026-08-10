#pragma once

#include "fluid/planar_region_fragment_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentPressureJumpEnergyVersion = 1;

struct PlanarPressureRegionFragmentPressureJumpEnergySettings {
    double timeStepSeconds = 1.0 / 60.0;
    double absolutePressureResidualTolerancePascals = 1.0e-12;
    double relativePressureResidualTolerance = 1.0e-10;
    double absoluteVelocityResidualToleranceMetersPerSecond = 1.0e-12;
    double relativeVelocityResidualTolerance = 1.0e-10;
    double absoluteForceResidualToleranceNewtons = 1.0e-12;
    double relativeForceResidualTolerance = 1.0e-10;
    double absoluteWorkResidualToleranceJoules = 1.0e-12;
    double relativeWorkResidualTolerance = 1.0e-10;

    bool operator==(
        const PlanarPressureRegionFragmentPressureJumpEnergySettings&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureJumpEnergyLayer {
    std::size_t layerIndex = 0;
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
    std::size_t minusTraceDofIndex = 0;
    std::size_t plusTraceDofIndex = 0;
    double areaSquareMeters = 0.0;
    Vector3 unitNormalMinusToPlus;
    double minusPressurePascals = 0.0;
    double plusPressurePascals = 0.0;
    double authoredPressureJumpPascals = 0.0;
    double reconstructedPressureJumpPascals = 0.0;
    double pressureJumpResidualPascals = 0.0;
    double materialWallVelocityMetersPerSecond = 0.0;
    double minusTraceVelocityMetersPerSecond = 0.0;
    double plusTraceVelocityMetersPerSecond = 0.0;
    double maximumAbsoluteWallVelocityResidualMetersPerSecond = 0.0;
    Vector3 resolvedPressureForceOnMinusFluidNewtons;
    Vector3 resolvedPressureForceOnPlusFluidNewtons;
    Vector3 resolvedPressureForceOnFluidNewtons;
    Vector3 authoredPressureJumpForceOnFluidNewtons;
    Vector3 pressureForceClosureResidualNewtons;
    Vector3 pressureForceOnSheetNewtons;
    Vector3 actionReactionForceResidualNewtons;
    Vector3 pressureJumpImpulseOnFluidNewtonSeconds;
    Vector3 pressureImpulseOnSheetNewtonSeconds;
    Vector3 actionReactionImpulseResidualNewtonSeconds;
    double resolvedPressurePowerToFluidWatts = 0.0;
    double authoredPressureJumpPowerToFluidWatts = 0.0;
    double pressurePowerClosureResidualWatts = 0.0;
    double pressureJumpWorkToFluidJoules = 0.0;
    double pressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureJumpEnergyLayer&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureJumpEnergyComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t pressureLayerSideCount = 0;
    Vector3 resolvedPressureForceOnFluidNewtons;
    Vector3 closedBoundaryForceResidualNewtons;
    Vector3 pressureImpulseOnFluidNewtonSeconds;
    double resolvedPressurePowerToFluidWatts = 0.0;
    double pressureWorkToFluidJoules = 0.0;
    double geometryPressureWorkToFluidJoules = 0.0;
    double workGeometryResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureJumpEnergyComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureJumpEnergyLimits {
    PlanarPressureRegionFragmentVelocityStateLimits velocityStateLimits;
    PlanarPressureRegionFragmentVolumeRateLimits volumeRateLimits;
    std::size_t maximumLayers = 60'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Read-only certificate for the authored pressure discontinuity across every
// regional pressure-layer wall. Absolute side pressures independently produce
// the two fluid forces; their sum must equal area times the authored signed
// jump. The published sheet force is the opposite reaction. Static wall traces
// remain zero. The moving overload binds each trace to material velocity and
// closes fluid pressure work to -dt*sum(p*dV/dt) per component and globally.
//
// This ledger applies no force or impulse, changes no velocity, and does not
// compose projection, momentum transport, topology rebase, or production
// worker ownership.
struct PlanarPressureRegionFragmentPressureJumpEnergyAudit {
    std::uint32_t version =
        planarPressureRegionFragmentPressureJumpEnergyVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t sourceVelocityStateFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t volumeRateFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    PlanarPressureRegionFragmentPressureJumpEnergySettings settings;
    std::vector<PlanarPressureRegionFragmentPressureJumpEnergyLayer> layers;
    std::vector<PlanarPressureRegionFragmentPressureJumpEnergyComponent>
        components;
    std::size_t pressureLayerTraceCount = 0;
    double maximumAbsolutePressureJumpResidualPascals = 0.0;
    double maximumAbsoluteWallVelocityResidualMetersPerSecond = 0.0;
    double maximumAbsoluteForceClosureResidualNewtons = 0.0;
    double maximumAbsoluteComponentClosedBoundaryForceResidualNewtons = 0.0;
    double maximumAbsoluteWorkClosureResidualJoules = 0.0;
    Vector3 resolvedPressureForceOnFluidNewtons;
    Vector3 authoredPressureJumpForceOnFluidNewtons;
    Vector3 pressureForceClosureResidualNewtons;
    Vector3 pressureForceOnSheetNewtons;
    Vector3 actionReactionForceResidualNewtons;
    Vector3 pressureJumpImpulseOnFluidNewtonSeconds;
    Vector3 pressureImpulseOnSheetNewtonSeconds;
    Vector3 actionReactionImpulseResidualNewtonSeconds;
    double resolvedPressurePowerToFluidWatts = 0.0;
    double authoredPressureJumpPowerToFluidWatts = 0.0;
    double pressurePowerClosureResidualWatts = 0.0;
    double pressureJumpWorkToFluidJoules = 0.0;
    double pressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;
    double geometryPressureWorkToFluidJoules = 0.0;
    double workGeometryResidualJoules = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureJumpEnergyAudit&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentPressureJumpEnergyAudit
auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings =
        {},
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits = {});

void validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits = {});

[[nodiscard]] PlanarPressureRegionFragmentPressureJumpEnergyAudit
auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings =
        {},
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits = {});

void validateMovingPlanarPressureRegionFragmentPressureJumpEnergyAudit(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits = {});

} // namespace simwing::fsi::fluid
