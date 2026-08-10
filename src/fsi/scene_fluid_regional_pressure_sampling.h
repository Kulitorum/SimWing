#pragma once

#include "fluid/planar_region_fragment_accepted_state.h"
#include "fluid/planar_region_fragment_opening_load_state.h"
#include "scene_fluid_pressure_traction.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalPressureSamplingVersion = 2;
inline constexpr std::uint32_t
    sceneFluidRegionalPressureLoadApplicationVersion = 2;

struct SceneFluidRegionalPressureSamplingLimits {
    std::size_t maximumSamples = 20'000'000;
    std::size_t maximumTiles = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionalPressureSampleBinding {
    std::size_t sampleIndex = 0;
    std::uint64_t stableId = 0;
    StableId triangleId = invalidStableId;
    std::size_t surfaceLoadTileIndex = 0;
    std::size_t sourcePressureWallIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    StableId surfaceStableId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    double areaSquareMeters = 0.0;
    double negativeSidePressurePascals = 0.0;
    double positiveSidePressurePascals = 0.0;
    double pressureDifferencePascals = 0.0;

    bool operator==(
        const SceneFluidRegionalPressureSampleBinding&) const = default;
};

struct SceneFluidRegionalPressureTileCoverage {
    std::size_t tileIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    StableId surfaceStableId = invalidStableId;
    std::size_t sampleCount = 0;
    double sourceAreaSquareMeters = 0.0;
    double sampledAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;

    bool operator==(
        const SceneFluidRegionalPressureTileCoverage&) const = default;
};

// Immutable bridge from one accepted sealed or opening-aware planar regional
// pressure endpoint to the authoritative scene quadrature. A regional surface
// ID must equal the scene sheet ID. Every clipped quadrature patch must resolve
// to exactly one pressure-wall tile with matching side regions, current plane,
// normal and material velocity. Sealed tiles must receive their complete wall
// area; opening-aware tiles receive only their retained-solid area, and a fully
// open tile remains explicit with zero samples and zero load.
//
// The sampled one-sided total pressures deliberately retain the regional
// correction gauge chosen by the accepted pressure state. Capture and
// evaluation are read-only. The separate explicit transactional application
// below accepts either sealed or opening-aware samples only after this complete
// conservative closure. None of these operations advances either solver or
// selects production pressure ownership.
struct SceneFluidRegionalPressureSampleSet {
    std::uint32_t version =
        sceneFluidRegionalPressureSamplingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t regionalAcceptedStateFingerprint = 0;
    std::uint64_t regionalOpeningLoadStateFingerprint = 0;
    std::uint64_t regionalPressureStateFingerprint = 0;
    std::uint64_t regionalSurfaceLoadFingerprint = 0;
    std::uint64_t regionalTopologyFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    bool openingAware = false;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    std::vector<SceneFluidRegionalPressureSampleBinding> bindings;
    std::vector<SceneFluidQuadraturePressure> pressures;
    std::vector<SceneFluidRegionalPressureTileCoverage> tiles;
    double sampledAreaSquareMeters = 0.0;
    double maximumAbsoluteTileAreaResidualSquareMeters = 0.0;
    StructureVector3 sampledPressureForceOnSheetNewtons;
    StructureVector3 sampledPressureMomentOnSheetNewtonMeters;
    double sampledPressurePowerToSheetWatts = 0.0;
    StructureVector3 sourceForceResidualNewtons;
    StructureVector3 sourceMomentResidualNewtonMeters;
    double sourcePowerResidualWatts = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalPressureSampleSet&) const = default;
};

struct SceneFluidRegionalPressureLoadApplicationLimits {
    std::size_t maximumNodeLoads = 20'000'000;
    std::size_t maximumStructureNodes = 20'000'000;
    std::size_t maximumOwnedBytes = 2048ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionalPressureAppliedNodeLoad {
    std::size_t loadIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;
    StructureVector3 priorPendingForceNewtons;
    StructureVector3 appliedPressureForceNewtons;
    StructureVector3 resultingPendingForceNewtons;
    StructureVector3 applicationResidualNewtons;

    bool operator==(
        const SceneFluidRegionalPressureAppliedNodeLoad&) const = default;
};

// Immutable receipt for one explicit pending-load mutation. The target must
// still be at the exact Structure epoch and nodal kinematics represented by
// the sampled scene state. All node bindings, arithmetic and storage limits
// are validated before the first load changes; every postcondition is checked
// against a saved Structure checkpoint and any failure restores that checkpoint.
// Existing non-pressure pending loads are retained exactly.
//
// This does not step Structure, consume the pending loads, persist regional
// state, or enable a production worker path.
struct SceneFluidRegionalPressureLoadApplication {
    std::uint32_t version =
        sceneFluidRegionalPressureLoadApplicationVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceSamplingFingerprint = 0;
    std::uint64_t sourceSurfaceStateFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t targetDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t structureNodeCount = 0;
    std::vector<SceneFluidRegionalPressureAppliedNodeLoad> nodeLoads;
    StructureVector3 priorPendingForceNewtons;
    StructureVector3 appliedPressureForceNewtons;
    StructureVector3 resultingPendingForceNewtons;
    StructureVector3 applicationResidualNewtons;
    bool applied = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalPressureLoadApplication&) const = default;
};

[[nodiscard]] SceneFluidRegionalPressureSampleSet
sampleSceneFluidRegionalAcceptedPressure(
    const fluid::PlanarPressureRegionFragmentAcceptedState& acceptedState,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSamplingLimits& limits = {});

[[nodiscard]] SceneFluidRegionalPressureSampleSet
sampleSceneFluidRegionalOpeningPressure(
    const fluid::PlanarPressureRegionFragmentOpeningLoadState& loadState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits&
        loadStateLimits = {},
    const SceneFluidRegionalPressureSamplingLimits& limits = {});

void validateSceneFluidRegionalPressureSampleIntegrity(
    const SceneFluidRegionalPressureSampleSet& samples);

void validateSceneFluidRegionalAcceptedPressureSamples(
    const SceneFluidRegionalPressureSampleSet& samples,
    const fluid::PlanarPressureRegionFragmentAcceptedState& acceptedState,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSamplingLimits& limits = {});

void validateSceneFluidRegionalOpeningPressureSamples(
    const SceneFluidRegionalPressureSampleSet& samples,
    const fluid::PlanarPressureRegionFragmentOpeningLoadState& loadState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits&
        loadStateLimits = {},
    const SceneFluidRegionalPressureSamplingLimits& limits = {});

[[nodiscard]] ConservativeTransferResult
evaluateSceneFluidRegionalAcceptedPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    const ConservativeTransferSettings& settings = {});

[[nodiscard]] SceneFluidRegionalPressureLoadApplication
applySceneFluidRegionalAcceptedPressureLoads(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    Structure& target,
    const ConservativeTransferSettings& transferSettings = {},
    const SceneFluidRegionalPressureLoadApplicationLimits& limits = {});

void validateSceneFluidRegionalPressureLoadApplicationIntegrity(
    const SceneFluidRegionalPressureLoadApplication& application);

} // namespace simwing::fsi
