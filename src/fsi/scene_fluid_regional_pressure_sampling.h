#pragma once

#include "fluid/planar_region_fragment_accepted_state.h"
#include "scene_fluid_pressure_traction.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalPressureSamplingVersion = 1;

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

// Immutable bridge from one accepted planar regional pressure endpoint to the
// authoritative scene quadrature. A regional surface ID must equal the scene
// sheet ID. Every clipped quadrature patch must resolve to exactly one
// pressure-wall tile with matching side regions, current plane, normal and
// material velocity; every tile must receive its complete area exactly once.
//
// The sampled one-sided total pressures deliberately retain the regional
// correction gauge chosen by the accepted pressure state. This adapter only
// evaluates the existing conservative Structure transfer; it does not mutate
// Structure, advance either solver, or select production pressure ownership.
struct SceneFluidRegionalPressureSampleSet {
    std::uint32_t version =
        sceneFluidRegionalPressureSamplingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t regionalAcceptedStateFingerprint = 0;
    std::uint64_t regionalPressureStateFingerprint = 0;
    std::uint64_t regionalSurfaceLoadFingerprint = 0;
    std::uint64_t regionalTopologyFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
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

[[nodiscard]] ConservativeTransferResult
evaluateSceneFluidRegionalAcceptedPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    const ConservativeTransferSettings& settings = {});

} // namespace simwing::fsi
