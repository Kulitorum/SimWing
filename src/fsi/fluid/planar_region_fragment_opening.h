#pragma once

#include "fluid/planar_region_fragment_topology.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningVersion = 2;

struct PlanarPressureRegionFragmentOpeningPatchDefinition {
    std::uint64_t patchStableId = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    double areaSquareMeters = 0.0;
    // Exact sub-tile centroid from the authoritative opening partition when
    // available. An omitted centroid preserves the original area-only model.
    std::optional<Vector3> authoredWrappedCentroidMeters;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPatchDefinition&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPatch {
    std::size_t patchIndex = 0;
    std::uint64_t patchStableId = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::size_t minusFragmentIndex = 0;
    std::size_t plusFragmentIndex = 0;
    std::uint64_t minusFragmentStableId = 0;
    std::uint64_t plusFragmentStableId = 0;
    std::size_t minusBaseComponentIndex = 0;
    std::size_t plusBaseComponentIndex = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    double areaSquareMeters = 0.0;
    double sourceWallAreaSquareMeters = 0.0;
    double sourceWallAreaFraction = 0.0;
    double centerDistanceMeters = 0.0;
    bool usesAuthoredCentroid = false;
    Vector3 wrappedCentroidMeters;
    Vector3 unitNormalNegativeToPositive;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPatch&) const = default;
};

struct PlanarPressureRegionFragmentOpeningWallPartition {
    std::size_t partitionIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    std::size_t openingPatchCount = 0;
    double wallAreaSquareMeters = 0.0;
    double openingAreaSquareMeters = 0.0;
    double solidAreaSquareMeters = 0.0;
    double openingAreaFraction = 0.0;
    bool hasExactSubtileCentroids = false;
    Vector3 openingAreaWeightedCentroidMeters;
    Vector3 solidAreaWeightedCentroidMeters;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningWallPartition&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningSummary {
    std::size_t openingIndex = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    std::size_t minusBaseComponentIndex = 0;
    std::size_t plusBaseComponentIndex = 0;
    std::size_t patchCount = 0;
    double areaSquareMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningSummary&) const = default;
};

struct PlanarPressureRegionFragmentOpeningBaseComponentMapping {
    std::size_t baseComponentIndex = 0;
    std::uint64_t baseComponentStableId = 0;
    std::size_t connectedComponentIndex = 0;
    std::uint64_t connectedComponentStableId = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningBaseComponentMapping&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningConnectedComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t fragmentCount = 0;
    std::size_t openingPatchCount = 0;
    double volumeCubicMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningConnectedComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningLimits {
    PlanarPressureRegionFragmentTopologyLimits topologyLimits;
    std::size_t maximumPatches = 20'000'000;
    std::size_t maximumPartitions = 20'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumConnectedComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Immutable aperture overlay on one accepted regional fragment topology.
// Every authored patch names one exact pressure-layer wall tile and partitions
// part of its area from solid fabric into an open negative-to-positive fluid
// connection. Multiple patches may share a tile or one opening, but one
// opening must keep one surface, orientation, side-region pair, and pair of
// base pressure components. Open area may never exceed the source wall area.
// Optional authoritative wrapped patch centroids must lie on and inside their
// wall tile. A fully centroid-authored partition derives the retained-solid
// centroid by exact first-moment subtraction; area-only input retains the wall
// centroid as an explicit compatibility model.
//
// The overlay also publishes the pressure-component union induced by those
// apertures. It changes no source topology and owns no conductance, flux,
// velocity degree, pressure solve, fabric-load subtraction, or worker state.
struct PlanarPressureRegionFragmentOpeningSet {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    std::vector<PlanarPressureRegionFragmentOpeningPatch> patches;
    std::vector<PlanarPressureRegionFragmentOpeningWallPartition> partitions;
    std::vector<PlanarPressureRegionFragmentOpeningSummary> openings;
    std::vector<PlanarPressureRegionFragmentOpeningBaseComponentMapping>
        baseComponents;
    std::vector<PlanarPressureRegionFragmentOpeningConnectedComponent>
        connectedComponents;
    double totalOpeningAreaSquareMeters = 0.0;
    double totalTouchedWallAreaSquareMeters = 0.0;
    double totalSolidAreaOnTouchedWallsSquareMeters = 0.0;
    double wallAreaPartitionResidualSquareMeters = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningSet&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningSet
buildPlanarPressureRegionFragmentOpenings(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningLimits& limits = {});

void validatePlanarPressureRegionFragmentOpenings(
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningLimits& limits = {});

} // namespace simwing::fsi::fluid
