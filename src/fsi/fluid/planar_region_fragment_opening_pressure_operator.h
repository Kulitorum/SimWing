#pragma once

#include "fluid/planar_region_fragment_opening.h"
#include "fluid/planar_region_fragment_pressure_operator.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningPressureOperatorVersion = 1;

enum class PlanarPressureRegionFragmentOpeningPressureEntryKind
    : std::uint8_t {
    SameRegionGrid = 1,
    OpeningPatch = 2,
};

struct PlanarPressureRegionFragmentOpeningPressureOperatorRow {
    std::size_t rowIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::size_t connectedComponentIndex = 0;
    bool isGauge = false;
    std::size_t firstEntry = 0;
    std::size_t entryCount = 0;
    double diagonalGeometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureOperatorRow&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureOperatorEntry {
    std::size_t entryIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentOpeningPressureEntryKind kind =
        PlanarPressureRegionFragmentOpeningPressureEntryKind::SameRegionGrid;
    std::size_t sourceIndex = 0;
    std::uint64_t sourceStableId = 0;
    std::size_t columnFragmentIndex = 0;
    std::uint64_t columnFragmentStableId = 0;
    double geometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureOperatorEntry&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureOperatorComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t gaugeFragmentIndex = 0;
    std::size_t firstFragmentMember = 0;
    std::size_t fragmentCount = 0;
    double totalVolumeCubicMeters = 0.0;
    double totalGeometryWeightMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureOperatorComponent&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningPressureOperatorLimits {
    PlanarPressureRegionFragmentPressureOperatorLimits baseOperatorLimits;
    PlanarPressureRegionFragmentOpeningLimits openingLimits;
    std::size_t maximumRows = 20'000'000;
    std::size_t maximumEntries = 140'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Ungauged integrated graph Laplacian over the exact opening overlay. The
// sealed base operator contributes every same-region Cartesian edge. Every
// aperture patch additionally contributes w=area/center-distance between its
// two source fragments; the remaining solid pressure-wall area contributes
// nothing. Opening-connected components own the gauges.
//
// This is projection geometry only. It owns no pressure RHS or solve, aperture
// velocity/momentum, constitutive resistance, authored-jump evolution,
// traction subtraction, or worker state.
struct PlanarPressureRegionFragmentOpeningPressureOperator {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningPressureOperatorVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceBaseOperatorFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::vector<PlanarPressureRegionFragmentOpeningPressureOperatorRow> rows;
    std::vector<PlanarPressureRegionFragmentOpeningPressureOperatorEntry>
        entries;
    std::vector<PlanarPressureRegionFragmentOpeningPressureOperatorComponent>
        components;
    std::vector<std::size_t> componentFragmentIndices;
    std::size_t includedSameRegionGridLinkCount = 0;
    std::size_t includedOpeningPatchCount = 0;
    double totalPressureLayerWallAreaSquareMeters = 0.0;
    double totalOpeningAreaSquareMeters = 0.0;
    double totalSolidPressureLayerWallAreaSquareMeters = 0.0;
    double wallAreaPartitionResidualSquareMeters = 0.0;
    double sameRegionGeometryWeightMeters = 0.0;
    double openingGeometryWeightMeters = 0.0;
    double totalGeometryWeightMeters = 0.0;
    double totalDiagonalGeometryWeightMeters = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureOperator&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningPressureOperator
buildPlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentPressureOperator& baseOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& baseOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningPressureOperatorLimits& limits =
        {});

[[nodiscard]] std::vector<double>
applyPlanarPressureRegionFragmentOpeningPressureOperator(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    std::span<const double> pressurePascals);

} // namespace simwing::fsi::fluid
