#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state.h"
#include "fluid/planar_region_fragment_opening_velocity_metric.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningVelocityStateVersion = 1;

struct PlanarPressureRegionFragmentOpeningVelocityStateSample {
    std::size_t dofIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentOpeningVelocityDofKind kind =
        PlanarPressureRegionFragmentOpeningVelocityDofKind::
            SharedRegionGrid;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::size_t sourceOpeningPatchIndex = 0;
    std::uint64_t sourceOpeningPatchStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t connectedComponentIndex = 0;
    double dualVolumeCubicMeters = 0.0;
    double materialNormalVelocityMetersPerSecond = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;
    double normalVelocityMetersPerSecond = 0.0;
    double diagonalMassKilograms = 0.0;
    double normalMomentumKilogramMetersPerSecond = 0.0;
    double diagonalKineticEnergyJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityStateSample&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningVelocityStateFragment {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t baseComponentIndex = 0;
    std::size_t connectedComponentIndex = 0;
    double volumeCubicMeters = 0.0;
    double massKilograms = 0.0;
    Vector3 collocatedVelocityMetersPerSecond;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 diagonalKineticEnergyByAxisJoules;
    double diagonalKineticEnergyJoules = 0.0;
    double collocatedKineticEnergyJoules = 0.0;
    double staggeringKineticEnergyJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityStateFragment&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningVelocityStateComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t fragmentCount = 0;
    double volumeCubicMeters = 0.0;
    double massKilograms = 0.0;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 diagonalKineticEnergyByAxisJoules;
    double diagonalKineticEnergyJoules = 0.0;
    double collocatedKineticEnergyJoules = 0.0;
    double staggeringKineticEnergyJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityStateComponent&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningVelocityStateLimits {
    PlanarPressureRegionFragmentOpeningVelocityMetricLimits metricLimits;
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits acceptedStateLimits;
    std::size_t maximumSamples = 140'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Immutable opening-aware regional velocity and momentum ownership. Every
// metric degree stores an absolute normal velocity as material plus relative
// velocity: fixed-grid links use zero material velocity, retained solid traces
// use zero relative velocity, and apertures retain both terms. Diagonal mass,
// normal momentum, and kinetic energy follow rho*dualVolume exactly.
//
// Each half-volume then contributes to its owning fragment. The resulting
// face-volume-weighted vector is one collocated momentum/velocity per physical
// fragment, ready for a later conservative ALE transport. Diagonal energy and
// the nonnegative staggering excess over collocated energy remain explicit.
// The accepted-endpoint capture overload derives all three velocity fields
// from validated pressure flow and material motion. Neither overload advances,
// projects, rebases, or selects production worker arithmetic.
struct PlanarPressureRegionFragmentOpeningVelocityState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningVelocityStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    bool mappedFromAcceptedEndpoint = false;
    double densityKgPerCubicMeter = 0.0;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityStateSample>
        samples;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityStateFragment>
        fragments;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityStateComponent>
        components;
    std::size_t sharedRegionGridSampleCount = 0;
    std::size_t solidWallTraceSampleCount = 0;
    std::size_t openingPatchSampleCount = 0;
    Vector3 diagonalMassByAxisKilograms;
    double physicalMassKilograms = 0.0;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 diagonalKineticEnergyByAxisJoules;
    double diagonalKineticEnergyJoules = 0.0;
    double collocatedKineticEnergyJoules = 0.0;
    double staggeringKineticEnergyJoules = 0.0;
    double maximumAbsoluteNormalVelocityMetersPerSecond = 0.0;
    double maximumCollocatedSpeedMetersPerSecond = 0.0;
    double maximumAbsoluteVelocityCompositionResidualMetersPerSecond = 0.0;
    double minimumFragmentStaggeringKineticEnergyJoules = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityState&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningVelocityState
buildPlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    std::span<const double> normalVelocityMetersPerSecond,
    std::span<const double> materialNormalVelocityMetersPerSecond,
    std::span<const double> relativeNormalVelocityMetersPerSecond,
    double densityKgPerCubicMeter,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits =
        {});

[[nodiscard]] PlanarPressureRegionFragmentOpeningVelocityState
capturePlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(
    const PlanarPressureRegionFragmentOpeningVelocityState& state);

void validatePlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityState& state,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityState& state,
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits =
        {});

} // namespace simwing::fsi::fluid
