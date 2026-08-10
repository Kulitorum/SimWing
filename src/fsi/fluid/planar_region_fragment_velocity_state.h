#pragma once

#include "fluid/planar_region_fragment_velocity_metric.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentVelocityStateVersion = 1;

struct PlanarPressureRegionFragmentVelocitySample {
    std::size_t dofIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentVelocityDofKind kind =
        PlanarPressureRegionFragmentVelocityDofKind::SharedRegionGrid;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t surfaceStableId = 0;
    std::size_t componentIndex = 0;
    std::uint64_t regionStableId = 0;
    double dualVolumeCubicMeters = 0.0;
    double normalVelocityMetersPerSecond = 0.0;
    double normalMomentumKilogramMetersPerSecond = 0.0;
    double kineticEnergyJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentVelocitySample&) const = default;
};

struct PlanarPressureRegionFragmentVelocityStateFragment {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t componentIndex = 0;
    Vector3 massByAxisKilograms;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 kineticEnergyByAxisJoules;
    double kineticEnergyJoules = 0.0;
    Vector3 massClosureResidualByAxisKilograms;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityStateFragment&) const =
        default;
};

struct PlanarPressureRegionFragmentVelocityStateComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    Vector3 massByAxisKilograms;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 kineticEnergyByAxisJoules;
    double kineticEnergyJoules = 0.0;
    Vector3 massClosureResidualByAxisKilograms;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityStateComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentVelocityStateLimits {
    PlanarPressureRegionFragmentVelocityMetricLimits metricLimits;
    std::size_t maximumSamples = 120'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

// Immutable diagonal-inertia state over one regional velocity metric. Input
// samples are ordered by metric DOF and point in the topology link's
// minus-to-plus Cartesian direction. Shared grid samples contribute their two
// half masses to both adjacent fragments. Pressure-layer traces remain two
// independent one-sided samples and contribute only to their owning fluid
// fragment, so no wall velocity or momentum is averaged across fabric.
//
// This state owns density, scalar normal momentum, and kinetic-energy ledgers.
// It does not prescribe wall velocity, prove that a field was pressure
// projected, transport momentum, apply pressure work, or enter production.
struct PlanarPressureRegionFragmentVelocityState {
    std::uint32_t version =
        planarPressureRegionFragmentVelocityStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    double densityKgPerCubicMeter = 0.0;
    std::vector<PlanarPressureRegionFragmentVelocitySample> samples;
    std::vector<PlanarPressureRegionFragmentVelocityStateFragment> fragments;
    std::vector<PlanarPressureRegionFragmentVelocityStateComponent>
        components;
    std::size_t sharedRegionGridSampleCount = 0;
    std::size_t pressureLayerTraceSampleCount = 0;
    Vector3 massByAxisKilograms;
    Vector3 momentumKilogramMetersPerSecond;
    Vector3 kineticEnergyByAxisJoules;
    double kineticEnergyJoules = 0.0;
    double maximumAbsoluteVelocityMetersPerSecond = 0.0;
    double maximumAbsoluteFragmentMassClosureResidualKilograms = 0.0;
    double maximumAbsoluteComponentMassClosureResidualKilograms = 0.0;
    Vector3 domainMassClosureResidualByAxisKilograms;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityState&) const = default;
};

void validatePlanarPressureRegionFragmentVelocityStateIntegrity(
    const PlanarPressureRegionFragmentVelocityState& state);

[[nodiscard]] PlanarPressureRegionFragmentVelocityState
buildPlanarPressureRegionFragmentVelocityState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::vector<double>& normalVelocityMetersPerSecond,
    double densityKgPerCubicMeter,
    const PlanarPressureRegionFragmentVelocityStateLimits& limits = {});

void validatePlanarPressureRegionFragmentVelocityState(
    const PlanarPressureRegionFragmentVelocityState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityStateLimits& limits = {});

} // namespace simwing::fsi::fluid
