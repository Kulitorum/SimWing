#pragma once

#include "fluid/planar_region_fragment_pressure_state.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentSurfaceLoadVersion = 1;

struct PlanarPressureRegionFragmentSurfaceLoadTile {
    std::size_t tileIndex = 0;
    std::size_t sourcePressureWallIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    double areaSquareMeters = 0.0;
    Vector3 wrappedCentroidMeters;
    Vector3 unitNormalMinusToPlus;
    Vector3 authoredPressureTractionOnSheetPascals;
    Vector3 correctionPressureTractionOnSheetPascals;
    Vector3 totalPressureTractionOnSheetPascals;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 totalPressureImpulseOnSheetNewtonSeconds;
    double totalPressureWorkToSheetJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentSurfaceLoadTile&) const = default;
};

struct PlanarPressureRegionFragmentSurfaceLoadSummary {
    std::size_t surfaceIndex = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    std::size_t tileCount = 0;
    double areaSquareMeters = 0.0;
    Vector3 areaWeightedCentroidMeters;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 totalPressureImpulseOnSheetNewtonSeconds;
    double totalPressureWorkToSheetJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentSurfaceLoadSummary&) const =
        default;
};

struct PlanarPressureRegionFragmentSurfaceLoadLimits {
    std::size_t maximumTiles = 60'000'000;
    std::size_t maximumSurfaces = 1'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

// Minimal immutable pressure-load handoff grouped by authored surface. Each
// tile retains the exact regional pressure-wall identity and geometry needed
// by a later conservative structural adapter. Force, time-integrated impulse,
// and material work remain split into authored/correction/total ownership and
// close back to the accepted composed pressure state.
//
// Capturing this ledger applies no Structure load, changes no fluid momentum,
// and owns no transport, topology rebase, or production continuation.
struct PlanarPressureRegionFragmentSurfaceLoadLedger {
    std::uint32_t version =
        planarPressureRegionFragmentSurfaceLoadVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourcePressureStateFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    std::vector<PlanarPressureRegionFragmentSurfaceLoadTile> tiles;
    std::vector<PlanarPressureRegionFragmentSurfaceLoadSummary> surfaces;
    double totalAreaSquareMeters = 0.0;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 pressureForceSplitResidualNewtons;
    Vector3 totalPressureImpulseOnSheetNewtonSeconds;
    Vector3 impulseTimeIntegrationResidualNewtonSeconds;
    double totalPressureWorkToSheetJoules = 0.0;
    double sourceWorkResidualJoules = 0.0;
    double maximumAbsoluteForceReconstructionResidualNewtons = 0.0;
    double maximumAbsoluteSurfaceAggregationResidualNewtons = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentSurfaceLoadLedger&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentSurfaceLoadLedger
capturePlanarPressureRegionFragmentSurfaceLoads(
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits = {});

void validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& ledger);

void validatePlanarPressureRegionFragmentSurfaceLoads(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& ledger,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits = {});

} // namespace simwing::fsi::fluid
