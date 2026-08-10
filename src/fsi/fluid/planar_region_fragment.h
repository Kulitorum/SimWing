#pragma once

#include "fluid/planar_region_sweep.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarPressureRegionFragmentVersion = 1;

enum class PlanarPressureRegionFragmentBoundaryKind : std::uint8_t {
    GridFace = 0,
    PressureLayer = 1,
};

struct PlanarPressureRegionFragmentBoundary {
    PlanarPressureRegionFragmentBoundaryKind kind =
        PlanarPressureRegionFragmentBoundaryKind::GridFace;
    std::uint64_t surfaceStableId = 0;
    std::size_t faceCoordinate = 0;
    std::int64_t periodicImage = 0;
    double unwrappedCoordinateMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentBoundary&) const = default;
};

struct PlanarPressureRegionFragment {
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::int64_t axisCellPeriodicImage = 0;
    PlanarPressureRegionFragmentBoundary lowerBoundary;
    PlanarPressureRegionFragmentBoundary upperBoundary;
    double unwrappedLowerCoordinateMeters = 0.0;
    double unwrappedUpperCoordinateMeters = 0.0;
    double unwrappedAxisCentroidMeters = 0.0;
    double transverseAreaSquareMeters = 0.0;
    double volumeCubicMeters = 0.0;
    Vector3 wrappedCentroidMeters;
    double pressurePascals = 0.0;

    bool operator==(const PlanarPressureRegionFragment&) const = default;
};

struct PlanarPressureRegionFragmentRegionSummary {
    std::uint64_t regionStableId = 0;
    std::size_t fragmentCount = 0;
    double volumeCubicMeters = 0.0;
    double sourceProfileVolumeCubicMeters = 0.0;
    double volumeClosureResidualCubicMeters = 0.0;
    double pressurePascals = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentRegionSummary&) const = default;
};

struct PlanarPressureRegionFragmentCellSummary {
    std::size_t cellIndex = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t fragmentCount = 0;
    double fragmentVolumeCubicMeters = 0.0;
    double volumeClosureResidualCubicMeters = 0.0;
    Vector3 fragmentFirstMomentCubicMetersSquared;
    Vector3 firstMomentClosureResidualCubicMetersSquared;

    bool operator==(
        const PlanarPressureRegionFragmentCellSummary&) const = default;
};

struct PlanarPressureRegionFragmentLimits {
    std::size_t maximumIntervals = 1'000'000;
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumCells = 10'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
};

// Splits the current static regional profile at every Cartesian cell face and
// publishes one full-dimensional control fragment per transverse cell tile.
// Multiple pressure-layer intervals may therefore coexist in one ordinary
// Cartesian cell with distinct region, pressure, volume, centroid, and
// boundary identity. Grid-face boundaries retain periodic image provenance;
// layer boundaries retain authored surface identity. This is immutable
// one-epoch geometry only, not a velocity basis or pressure solve.
struct PlanarPressureRegionFragmentSet {
    std::uint32_t version = planarPressureRegionFragmentVersion;
    std::uint64_t fingerprint = 0;
    std::uint32_t sourceSweepVersion = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    Vector3 spacingMeters;
    double profileWindowLowerCoordinateMeters = 0.0;
    double profileWindowUpperCoordinateMeters = 0.0;
    std::vector<PlanarPressureRegionFragment> fragments;
    std::vector<PlanarPressureRegionFragmentRegionSummary> regions;
    std::vector<PlanarPressureRegionFragmentCellSummary> cells;
    std::size_t maximumFragmentsPerCell = 0;
    double fragmentVolumeCubicMeters = 0.0;
    double domainVolumeClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteRegionVolumeClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteCellVolumeClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared =
        0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentSet&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentSet
buildPlanarPressureRegionFragments(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentLimits& limits = {});

void validatePlanarPressureRegionFragments(
    const PlanarPressureRegionFragmentSet& fragments,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentLimits& limits = {});

} // namespace simwing::fsi::fluid
