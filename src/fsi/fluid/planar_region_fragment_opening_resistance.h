#pragma once

#include "fluid/planar_region_fragment_opening_flux.h"
#include "fluid/porous_flow.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentOpeningResistanceDefinition {
    std::uint64_t patchStableId = 0;
    DarcyForchheimerResistance resistance;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningResistanceDefinition&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningResistanceSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningResistanceSettings&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningResistancePatchDiagnostics {
    std::size_t patchIndex = 0;
    std::uint64_t patchStableId = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    double areaSquareMeters = 0.0;
    double centerDistanceMeters = 0.0;
    DarcyForchheimerResistance resistance;
    bool zeroResistanceIdentity = false;
    PorousPlugFlowDiagnostics plugFlow;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningResistancePatchDiagnostics&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningResistanceLimits {
    PlanarPressureRegionFragmentOpeningFluxLimits openingFluxLimits;
    std::size_t maximumPatches = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningResistanceDiagnostics {
    bool accepted = false;
    bool finite = false;
    bool nonIncreasingKineticEnergy = false;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceOpeningFluxFingerprint = 0;
    std::uint64_t resultOpeningFluxFingerprint = 0;
    std::uint64_t resistanceDefinitionFingerprint = 0;
    PlanarPressureRegionFragmentOpeningResistanceSettings settings;
    std::vector<
        PlanarPressureRegionFragmentOpeningResistancePatchDiagnostics>
        patches;
    std::size_t zeroResistancePatchCount = 0;
    double maximumAbsoluteMidpointPressureDropPascals = 0.0;
    double maximumAbsoluteMomentumResidualKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double dissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;
    double energyToleranceJoules = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningResistanceDiagnostics&) const =
        default;
};

// Applies one implicit-midpoint Darcy-Forchheimer decay to every exact
// aperture patch using fluid slug length equal to that patch's fragment-center
// distance. Coefficients are calibrated pressure-drop-per-relative-velocity
// values keyed by stable patch identity. A zero/zero pair is an explicit
// bit-exact inviscid identity; active pairs reuse the canonical porous plug
// oracle and close pressure impulse, momentum, midpoint work, and kinetic-
// energy dissipation independently per patch and in aggregate.
//
// Opening samples and their immutable flux state commit together only after
// the aggregate passive-energy check succeeds. This operator supplies no
// driving pressure, geometry projection, authored static-pressure evolution,
// traction, or production ownership.
[[nodiscard]] PlanarPressureRegionFragmentOpeningResistanceDiagnostics
advancePlanarPressureRegionFragmentOpeningResistance(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    const PlanarPressureRegionFragmentOpeningResistanceSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningResistanceLimits& limits = {});

} // namespace simwing::fsi::fluid
