#pragma once

#include "scene_fluid_mimetic_pressure_epoch.h"
#include "scene_fluid_mimetic_trace_flow.h"
#include "scene_fluid_pressure_epoch.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticPressureAuditVersion = 1;

struct SceneFluidMimeticPressureAuditSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    SceneFluidMimeticControlCellSettings controlCells;
    SceneFluidMimeticTraceSystemSettings traceSystem;
    SceneFluidMimeticTraceSolveSettings pressureSolve;

    bool operator==(
        const SceneFluidMimeticPressureAuditSettings&) const = default;
};

struct SceneFluidMimeticPressureAuditLimits {
    SceneFluidMimeticControlCellLimits controlCells;
    SceneFluidMimeticTraceSystemLimits traceSystem;
    SceneFluidMimeticCondensedTraceSystemLimits condensedSystem;
    SceneFluidMimeticTraceFlowLimits traceFlow;
    SceneFluidMimeticPressureSourceLimits pressureSource;
    SceneFluidMimeticPressureEpochLimits pressureEpoch;
    std::size_t maximumOwnedBytes =
        16ULL * 1024ULL * 1024ULL * 1024ULL;
};

// Complete immutable endpoint for one opt-in shadow pressure epoch. It keeps
// the rebuilt mixed-hybrid topology, the exact physical predictor/source
// chain, and the atomically accepted pressure state plus material samples
// together. Nothing in this product applies loads or changes the production
// graph-pressure path.
struct SceneFluidMimeticPressureAuditEndpoint {
    std::uint32_t version = sceneFluidMimeticPressureAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t scenePressureEpochFingerprint = 0;
    std::uint64_t pressureTopologyTransitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    bool usesRegionWallPrediction = false;
    bool usesConsecutiveWarmStart = false;
    std::size_t ownedStorageBytes = 0;
    SceneFluidMimeticControlCellSet controlCells;
    SceneFluidMimeticTraceSystem fullTraceSystem;
    SceneFluidMimeticCondensedTraceSystem condensedTraceSystem;
    SceneFluidMimeticTraceFlowPrediction predictedTraceFlows;
    SceneFluidMimeticPressureSourceSet pressureSources;
    SceneFluidMimeticPressureEpochResult pressureEpoch;

    bool operator==(
        const SceneFluidMimeticPressureAuditEndpoint&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& gridEpoch,
    const SceneFluidOpeningCapSet& openingCaps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& pressureFaceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditEndpoint& previousEndpoint,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
buildSceneFluidMimeticPressureAuditEndpoint(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditEndpoint& previousEndpoint,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

void validateSceneFluidMimeticPressureAuditEndpointIntegrity(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint);

} // namespace simwing::fsi
