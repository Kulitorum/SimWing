#pragma once

#include "fluid/planar_region_fragment_opening.h"
#include "fluid/planar_region_fragment_surface_load.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningSurfaceLoadVersion = 1;

struct PlanarPressureRegionFragmentOpeningSurfaceLoadTile {
    std::size_t tileIndex = 0;
    std::size_t sourceSurfaceLoadTileIndex = 0;
    std::size_t sourcePressureWallIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    bool touchedByOpening = false;
    std::size_t sourceOpeningPartitionIndex = 0;
    std::size_t openingPatchCount = 0;
    double wallAreaSquareMeters = 0.0;
    double openingAreaSquareMeters = 0.0;
    double solidAreaSquareMeters = 0.0;
    double openingAreaFraction = 0.0;
    Vector3 wrappedCentroidMeters;
    Vector3 unitNormalMinusToPlus;
    Vector3 authoredPressureTractionOnSheetPascals;
    Vector3 correctionPressureTractionOnSheetPascals;
    Vector3 totalPressureTractionOnSheetPascals;
    Vector3 openingRemovedAuthoredPressureForceOnSheetNewtons;
    Vector3 openingRemovedCorrectionPressureForceOnSheetNewtons;
    Vector3 openingRemovedTotalPressureForceOnSheetNewtons;
    Vector3 solidAuthoredPressureForceOnSheetNewtons;
    Vector3 solidCorrectionPressureForceOnSheetNewtons;
    Vector3 solidTotalPressureForceOnSheetNewtons;
    Vector3 openingRemovedTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 solidTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 openingRemovedTotalPressureMomentOnSheetNewtonMeters;
    Vector3 solidTotalPressureMomentOnSheetNewtonMeters;
    double openingRemovedTotalPressureWorkToSheetJoules = 0.0;
    double solidTotalPressureWorkToSheetJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningSurfaceLoadTile&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningSurfaceLoadSummary {
    std::size_t surfaceIndex = 0;
    std::size_t sourceSurfaceLoadSurfaceIndex = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    std::size_t tileCount = 0;
    std::size_t openingTouchedTileCount = 0;
    std::size_t fullyOpenTileCount = 0;
    double wallAreaSquareMeters = 0.0;
    double openingAreaSquareMeters = 0.0;
    double solidAreaSquareMeters = 0.0;
    Vector3 openingAreaWeightedCentroidMeters;
    Vector3 solidAreaWeightedCentroidMeters;
    Vector3 openingRemovedAuthoredPressureForceOnSheetNewtons;
    Vector3 openingRemovedCorrectionPressureForceOnSheetNewtons;
    Vector3 openingRemovedTotalPressureForceOnSheetNewtons;
    Vector3 solidAuthoredPressureForceOnSheetNewtons;
    Vector3 solidCorrectionPressureForceOnSheetNewtons;
    Vector3 solidTotalPressureForceOnSheetNewtons;
    Vector3 openingRemovedTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 solidTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 openingRemovedTotalPressureMomentOnSheetNewtonMeters;
    Vector3 solidTotalPressureMomentOnSheetNewtonMeters;
    double openingRemovedTotalPressureWorkToSheetJoules = 0.0;
    double solidTotalPressureWorkToSheetJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningSurfaceLoadSummary&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningSurfaceLoadLimits {
    PlanarPressureRegionFragmentSurfaceLoadLimits surfaceLoadLimits;
    PlanarPressureRegionFragmentOpeningLimits openingLimits;
    std::size_t maximumTiles = 60'000'000;
    std::size_t maximumSurfaces = 1'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Opt-in conservative pressure-load area partition for an exact aperture
// overlay. The source pressure traction is uniform on each regional wall tile.
// Every opening patch currently inherits that tile's wrapped centroid, so the
// opening-removed and retained-solid force, moment, impulse, and material work
// are the exact area fractions of that finite-volume tile. A future sub-tile
// opening centroid must be authored before this contract may claim a different
// moment arm.
//
// This ledger mutates neither fluid nor Structure state and is not selected by
// a production worker. It preserves the complete source tile set, including a
// fully open tile with zero retained solid load, so area and load ownership
// remain explicit and auditable.
struct PlanarPressureRegionFragmentOpeningSurfaceLoadLedger {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningSurfaceLoadVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceSurfaceLoadFingerprint = 0;
    std::uint64_t sourcePressureStateFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    std::vector<PlanarPressureRegionFragmentOpeningSurfaceLoadTile> tiles;
    std::vector<PlanarPressureRegionFragmentOpeningSurfaceLoadSummary>
        surfaces;
    double totalWallAreaSquareMeters = 0.0;
    double totalOpeningAreaSquareMeters = 0.0;
    double totalSolidAreaSquareMeters = 0.0;
    double wallAreaPartitionResidualSquareMeters = 0.0;
    Vector3 sourceAuthoredPressureForceOnSheetNewtons;
    Vector3 sourceCorrectionPressureForceOnSheetNewtons;
    Vector3 sourceTotalPressureForceOnSheetNewtons;
    Vector3 openingRemovedAuthoredPressureForceOnSheetNewtons;
    Vector3 openingRemovedCorrectionPressureForceOnSheetNewtons;
    Vector3 openingRemovedTotalPressureForceOnSheetNewtons;
    Vector3 solidAuthoredPressureForceOnSheetNewtons;
    Vector3 solidCorrectionPressureForceOnSheetNewtons;
    Vector3 solidTotalPressureForceOnSheetNewtons;
    Vector3 sourceTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 openingRemovedTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 solidTotalPressureImpulseOnSheetNewtonSeconds;
    Vector3 sourceTotalPressureMomentOnSheetNewtonMeters;
    Vector3 openingRemovedTotalPressureMomentOnSheetNewtonMeters;
    Vector3 solidTotalPressureMomentOnSheetNewtonMeters;
    double sourceTotalPressureWorkToSheetJoules = 0.0;
    double openingRemovedTotalPressureWorkToSheetJoules = 0.0;
    double solidTotalPressureWorkToSheetJoules = 0.0;
    double maximumAbsoluteAreaPartitionResidualSquareMeters = 0.0;
    double maximumAbsoluteForcePartitionResidualNewtons = 0.0;
    double maximumAbsoluteForceSplitResidualNewtons = 0.0;
    double maximumAbsoluteImpulsePartitionResidualNewtonSeconds = 0.0;
    double maximumAbsoluteMomentPartitionResidualNewtonMeters = 0.0;
    double maximumAbsoluteSurfaceAggregationResidualNewtons = 0.0;
    double workPartitionResidualJoules = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningSurfaceLoadLedger
capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger);

void validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLedger& ledger,
    const PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningSurfaceLoadLimits& limits = {});

} // namespace simwing::fsi::fluid
