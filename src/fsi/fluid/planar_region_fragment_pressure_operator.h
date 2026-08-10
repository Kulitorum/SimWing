#pragma once

#include "fluid/planar_region_fragment_topology.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentPressureOperatorVersion = 1;

struct PlanarPressureRegionFragmentPressureOperatorRow {
    std::size_t rowIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::size_t componentIndex = 0;
    bool isGauge = false;
    std::size_t firstEntry = 0;
    std::size_t entryCount = 0;
    double diagonalGeometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureOperatorRow&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureOperatorEntry {
    std::size_t entryIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::size_t columnFragmentIndex = 0;
    std::uint64_t columnFragmentStableId = 0;
    double geometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureOperatorEntry&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureOperatorComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t gaugeFragmentIndex = 0;
    std::size_t firstFragmentMember = 0;
    std::size_t fragmentCount = 0;
    double totalVolumeCubicMeters = 0.0;
    double totalGeometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureOperatorComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureOperatorLimits {
    PlanarPressureRegionFragmentTopologyLimits topologyLimits;
    std::size_t maximumRows = 20'000'000;
    std::size_t maximumEntries = 120'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Ungauged integrated finite-volume graph Laplacian over the immutable
// regional fragments. Every same-region grid link contributes +w to both row
// diagonals and one directed -w entry to each endpoint, with w=area/distance.
// Pressure-layer walls contribute no entries and remain disconnected pressure
// boundaries. The product owns neither a right-hand side, velocity field,
// linear solve, nor production pressure state.
struct PlanarPressureRegionFragmentPressureOperator {
    std::uint32_t version =
        planarPressureRegionFragmentPressureOperatorVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::vector<PlanarPressureRegionFragmentPressureOperatorRow> rows;
    std::vector<PlanarPressureRegionFragmentPressureOperatorEntry> entries;
    std::vector<PlanarPressureRegionFragmentPressureOperatorComponent>
        components;
    std::vector<std::size_t> componentFragmentIndices;
    std::size_t includedSameRegionGridLinkCount = 0;
    std::size_t excludedPressureLayerWallLinkCount = 0;
    double totalGeometryWeightMeters = 0.0;
    double totalDiagonalGeometryWeightMeters = 0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentPressureOperator&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentPressureOperator
buildPlanarPressureRegionFragmentPressureOperator(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits = {});

void validatePlanarPressureRegionFragmentPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentPressureOperatorLimits& limits = {});

[[nodiscard]] std::vector<double>
applyPlanarPressureRegionFragmentPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    std::span<const double> pressurePascals);

} // namespace simwing::fsi::fluid
