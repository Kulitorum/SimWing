#pragma once

#include "scene_fluid_mimetic_pressure_epoch.h"
#include "scene_fluid_mimetic_trace_flow.h"
#include "scene_fluid_pressure_epoch.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticPressureAuditVersion = 1;
inline constexpr std::uint32_t
    sceneFluidMimeticPressureAuditTopologyVersion = 1;
inline constexpr std::uint32_t
    sceneFluidMimeticPressureAuditWarmStateVersion = 1;

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

// Rebuilt, non-persistent topology used to validate and resume one compact
// accepted SWMP state. Persistence stores only the pressure state; every
// control shell and trace operator is reconstructed from the trusted scene
// and restored Structure checkpoint before that state is published.
struct SceneFluidMimeticPressureAuditTopology {
    std::uint32_t version =
        sceneFluidMimeticPressureAuditTopologyVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t scenePressureEpochFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    SceneFluidMimeticControlCellSet controlCells;
    SceneFluidMimeticTraceSystem fullTraceSystem;
    SceneFluidMimeticCondensedTraceSystem condensedTraceSystem;

    bool operator==(
        const SceneFluidMimeticPressureAuditTopology&) const = default;
};

// In-memory restart owner. The topology is deliberately rebuilt rather than
// decoded, while the accepted pressure rows retain their original source and
// epoch provenance for the next consecutive warm remap.
struct SceneFluidMimeticPressureAuditWarmState {
    std::uint32_t version =
        sceneFluidMimeticPressureAuditWarmStateVersion;
    std::uint64_t fingerprint = 0;
    std::size_t ownedStorageBytes = 0;
    SceneFluidMimeticPressureAuditTopology topology;
    SceneFluidMimeticPressureState acceptedPressureState;

    bool operator==(
        const SceneFluidMimeticPressureAuditWarmState&) const = default;
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

[[nodiscard]] std::uint64_t
sceneFluidMimeticPressureAuditSettingsFingerprint(
    const SceneFluidMimeticPressureAuditSettings& settings);

[[nodiscard]] SceneFluidMimeticPressureAuditTopology
buildSceneFluidMimeticPressureAuditTopology(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureEpoch& pressureEpoch,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureAuditWarmState
bindSceneFluidMimeticPressureAuditWarmState(
    SceneFluidMimeticPressureAuditTopology topology,
    SceneFluidMimeticPressureState acceptedPressureState,
    const SceneFluidMimeticPressureAuditLimits& limits = {});

void validateSceneFluidMimeticPressureAuditTopologyIntegrity(
    const SceneFluidMimeticPressureAuditTopology& topology);

void validateSceneFluidMimeticPressureAuditWarmStateIntegrity(
    const SceneFluidMimeticPressureAuditWarmState& warmState);

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

// Reuses one already accepted fixed-geometry mixed-hybrid topology with a new
// complete MAC predictor/opening-flux ledger. Only predicted trace flow,
// physical sources, and the atomic pressure epoch are rebuilt. This is a
// zero-warm-start fixed-topology continuation; moving geometry and topology
// transitions must use the consecutive-epoch overloads below.
[[nodiscard]] SceneFluidMimeticPressureAuditEndpoint
advanceSceneFluidMimeticPressureAuditFixedTopology(
    const SceneFluidMimeticPressureAuditEndpoint& acceptedTopology,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& pressureFaceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
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
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidPressureVolumeRateSet& geometryVolumeRates,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureAuditWarmState& previousWarmState,
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
    const SceneFluidMimeticPressureAuditWarmState& previousWarmState,
    const SceneFluidMimeticPressureAuditSettings& settings = {},
    const SceneFluidMimeticPressureAuditLimits& limits = {});

void validateSceneFluidMimeticPressureAuditEndpointIntegrity(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint);

} // namespace simwing::fsi
