#pragma once

#include "fluid/planar_region_fragment_topology.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentVolumeRateVersion = 1;

struct PlanarPressureRegionFragmentVolumeRate {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t cellIndex = 0;
    double lowerBoundaryDisplacementMeters = 0.0;
    double upperBoundaryDisplacementMeters = 0.0;
    double lowerBoundaryVelocityMetersPerSecond = 0.0;
    double upperBoundaryVelocityMetersPerSecond = 0.0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentVolumeRate&) const = default;
};

struct PlanarPressureRegionFragmentCellVolumeRate {
    std::size_t cellIndex = 0;
    std::size_t fragmentCount = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;
    double fixedCellVolumeClosureResidualCubicMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentCellVolumeRate&) const = default;
};

struct PlanarPressureRegionFragmentRegionVolumeRate {
    std::uint64_t regionStableId = 0;
    std::size_t fragmentCount = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double sourceSweepGeometryVolumeChangeCubicMeters = 0.0;
    double geometryVolumeChangeClosureResidualCubicMeters = 0.0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentRegionVolumeRate&) const = default;
};

struct PlanarPressureRegionFragmentComponentVolumeRate {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t fragmentCount = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentComponentVolumeRate&) const =
        default;
};

struct PlanarPressureRegionFragmentVolumeRateLimits {
    PlanarPressureRegionFragmentTopologyLimits topologyLimits;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumCells = 10'000'000;
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Reconstructs a previous volume and constant geometry-volume rate for every
// current regional fragment from its exact pressure-layer boundary
// displacements. Grid boundaries remain fixed. This first local moving-volume
// map requires every authored layer to stay inside the same unwrapped
// Cartesian segment, so current fragment identity is also valid at the
// previous endpoint. Cell, region, component, and global ledgers close the
// independently reconstructed rates back to the source sweep.
//
// The immutable result supplies geometry dV/dt only. It owns no face flow,
// pressure RHS, topology-transition retirement, or production state.
struct PlanarPressureRegionFragmentVolumeRateSet {
    std::uint32_t version =
        planarPressureRegionFragmentVolumeRateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint32_t sourceSweepVersion = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    double durationSeconds = 0.0;
    bool topologyStable = false;
    std::vector<PlanarPressureRegionFragmentVolumeRate> fragments;
    std::vector<PlanarPressureRegionFragmentCellVolumeRate> cells;
    std::vector<PlanarPressureRegionFragmentRegionVolumeRate> regions;
    std::vector<PlanarPressureRegionFragmentComponentVolumeRate> components;
    double maximumAbsoluteFragmentVolumeChangeCubicMeters = 0.0;
    double maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCellClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteRegionClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteComponentVolumeRateCubicMetersPerSecond = 0.0;
    double globalGeometryVolumeChangeCubicMeters = 0.0;
    double globalGeometryVolumeChangeRateCubicMetersPerSecond = 0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentVolumeRateSet&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentVolumeRateSet
buildPlanarPressureRegionFragmentVolumeRates(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateLimits& limits = {});

void validatePlanarPressureRegionFragmentVolumeRates(
    const PlanarPressureRegionFragmentVolumeRateSet& rates,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateLimits& limits = {});

} // namespace simwing::fsi::fluid
