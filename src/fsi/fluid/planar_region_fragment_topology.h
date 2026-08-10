#pragma once

#include "fluid/planar_region_fragment.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentTopologyVersion = 1;

enum class PlanarPressureRegionFragmentFaceKind : std::uint8_t {
    SameRegionGrid = 1,
    PressureLayerWall = 2,
};

struct PlanarPressureRegionFragmentFaceLink {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentFaceKind kind =
        PlanarPressureRegionFragmentFaceKind::SameRegionGrid;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::int64_t facePeriodicImage = 0;
    std::uint64_t surfaceStableId = 0;
    std::size_t minusFragmentIndex = 0;
    std::size_t plusFragmentIndex = 0;
    std::uint64_t minusFragmentStableId = 0;
    std::uint64_t plusFragmentStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    std::size_t minusComponentIndex = 0;
    std::size_t plusComponentIndex = 0;
    double areaSquareMeters = 0.0;
    Vector3 wrappedCentroidMeters;
    double centerDistanceMeters = 0.0;
    // Nonzero only for a same-region grid link. A pressure-layer wall has no
    // implied conductance even though both geometric sides are retained.
    double sameRegionGeometryWeightMeters = 0.0;
    double pressureJumpPascals = 0.0;
    Vector3 unitNormalMinusToPlus;
    bool crossesPeriodicBoundary = false;

    bool operator==(
        const PlanarPressureRegionFragmentFaceLink&) const = default;
};

struct PlanarPressureRegionFragmentTopologySummary {
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t incidentFaceCount = 0;
    std::size_t sameRegionGridFaceCount = 0;
    std::size_t pressureLayerWallFaceCount = 0;
    double incidentFaceAreaSquareMeters = 0.0;
    double expectedBoundaryAreaSquareMeters = 0.0;
    double boundaryAreaClosureResidualSquareMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentTopologySummary&) const = default;
};

struct PlanarPressureRegionFragmentComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t fragmentCount = 0;
    std::size_t sameRegionGridLinkCount = 0;
    double volumeCubicMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentComponent&) const = default;
};

struct PlanarPressureRegionFragmentTopologyLimits {
    PlanarPressureRegionFragmentLimits fragmentLimits;
    std::size_t maximumLinks = 60'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Pairs all six rectangular half-faces of every one-epoch regional fragment.
// Cartesian faces connect only matching regional controls; authored pressure
// layers instead publish a two-sided wall with its exact signed static pressure
// jump and zero conductance. The resulting same-region graph owns deterministic
// component identity but no velocity, momentum, pressure solve, or opening
// connection.
struct PlanarPressureRegionFragmentTopology {
    std::uint32_t version =
        planarPressureRegionFragmentTopologyVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    Vector3 spacingMeters;
    std::vector<PlanarPressureRegionFragmentFaceLink> links;
    std::vector<PlanarPressureRegionFragmentTopologySummary> fragments;
    std::vector<PlanarPressureRegionFragmentComponent> components;
    std::size_t sameRegionGridLinkCount = 0;
    std::size_t pressureLayerWallLinkCount = 0;
    std::size_t periodicGridLinkCount = 0;
    double totalUniqueFaceAreaSquareMeters = 0.0;
    double totalIncidentFaceAreaSquareMeters = 0.0;
    double totalExpectedFragmentBoundaryAreaSquareMeters = 0.0;
    double maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters =
        0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentTopology&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentTopology
buildPlanarPressureRegionFragmentTopology(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopologyLimits& limits = {});

void validatePlanarPressureRegionFragmentTopology(
    const PlanarPressureRegionFragmentTopology& topology,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopologyLimits& limits = {});

} // namespace simwing::fsi::fluid
