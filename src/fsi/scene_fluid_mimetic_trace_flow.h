#pragma once

#include "scene_fluid_mimetic_trace_system.h"
#include "scene_fluid_opening_flux.h"
#include "scene_fluid_region_wall.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticTraceFlowVersion = 2;
inline constexpr std::uint32_t
    sceneFluidMimeticTraceFlowContinuationVersion = 1;
inline constexpr std::uint32_t
    sceneFluidMimeticRegionTransportFlowPredictionVersion = 1;

struct SceneFluidMimeticCorrectedTraceFlow;
struct SceneFluidRegionTransport;

struct SceneFluidMimeticTraceFlowLimits {
    std::size_t maximumSharedTraces = 200'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticPredictedTraceFlow {
    std::size_t sharedTraceOrdinal = 0;
    std::size_t traceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::uint64_t sourceStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t minusControlCellIndex = 0;
    std::size_t plusControlCellIndex = 0;
    double predictedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidMimeticPredictedTraceFlow&) const = default;
};

// Immutable predictor over every shared mimetic trace. A fixed-epoch bootstrap
// samples exact partitioned MAC faces and accepted relative opening flux. The
// transported overload instead projects accepted material-wall-adjusted region
// velocities and subtracts the same accepted cap sweep. Positive flow is
// oriented from the trace's MinusOrNegative control to its PlusOrPositive
// control. Material-wall traces remain impermeable and are intentionally
// absent. This product samples predictor flow only; it does not assemble a
// pressure source or apply a correction.
struct SceneFluidMimeticTraceFlowPrediction {
    std::uint32_t version = sceneFluidMimeticTraceFlowVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t mimeticTraceSystemFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t openingFluxFingerprint = 0;
    std::uint64_t velocityFingerprint = 0;
    std::uint64_t regionWallExchangeFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double sourceDensityKgPerCubicMeter = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t componentCount = 0;
    std::size_t cartesianTraceCount = 0;
    std::size_t authoredOpeningTraceCount = 0;
    double maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond = 0.0;
    std::vector<SceneFluidMimeticPredictedTraceFlow> traces;
    std::vector<double> componentBalanceResidualsCubicMetersPerSecond;

    bool operator==(
        const SceneFluidMimeticTraceFlowPrediction&) const = default;
};

// Fixed-topology continuation which preserves the last accepted corrected
// shared-trace flow and applies only the change between two complete bulk-MAC
// predictor samples. This keeps Cartesian subface differences and embedded
// opening flow explicit instead of reconstructing the next predictor from a
// lossy one-velocity-per-face collapse. It is a diagnostic continuation
// model: transport still owns the bulk increment and this product neither
// advects subface residuals nor supplies a pressure warm start.
struct SceneFluidMimeticTraceFlowContinuation {
    std::uint32_t version =
        sceneFluidMimeticTraceFlowContinuationVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t previousCorrectedTraceFlowFingerprint = 0;
    std::uint64_t previousBaselinePredictionFingerprint = 0;
    std::uint64_t currentBaselinePredictionFingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t mimeticTraceSystemFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    double maximumAbsoluteCarriedCorrectionCubicMetersPerSecond = 0.0;
    double maximumAbsoluteBulkIncrementCubicMetersPerSecond = 0.0;
    bool finite = false;
    SceneFluidMimeticTraceFlowPrediction prediction;

    bool operator==(
        const SceneFluidMimeticTraceFlowContinuation&) const = default;
};

// Fixed-topology predictor which projects one accepted regional momentum
// transport state onto the existing shared traces. The nested prediction keeps
// the current bulk-MAC sample's exact opening/cap-sweep provenance, while this
// wrapper records that transported per-region velocities replaced its bulk
// link velocities. It is deliberately separate from the consecutive-epoch
// region-link predictor: no geometry epoch advances and pressure still starts
// from zero.
struct SceneFluidMimeticRegionTransportFlowPrediction {
    std::uint32_t version =
        sceneFluidMimeticRegionTransportFlowPredictionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t regionTransportFingerprint = 0;
    std::uint64_t sourceCorrectedTraceFlowFingerprint = 0;
    std::uint64_t currentBulkBaselinePredictionFingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t mimeticTraceSystemFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double transportTargetSimulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    double maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond =
        0.0;
    bool finite = false;
    SceneFluidMimeticTraceFlowPrediction prediction;

    bool operator==(
        const SceneFluidMimeticRegionTransportFlowPrediction&) const =
        default;
};

[[nodiscard]] SceneFluidMimeticTraceFlowPrediction
sampleSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticTraceFlowLimits& limits = {});

[[nodiscard]] SceneFluidMimeticTraceFlowContinuation
continueSceneFluidMimeticTraceFlowsFixedTopology(
    const SceneFluidMimeticCorrectedTraceFlow& previousCorrectedFlow,
    const SceneFluidMimeticTraceFlowPrediction& previousBulkBaseline,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline,
    const SceneFluidMimeticTraceFlowLimits& limits = {});

[[nodiscard]] SceneFluidMimeticRegionTransportFlowPrediction
predictSceneFluidMimeticTraceFlowsFromRegionTransportFixedTopology(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidMimeticCorrectedTraceFlow& sourceCorrectedFlow,
    const SceneFluidRegionTransport& regionTransport,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline,
    const SceneFluidMimeticTraceFlowLimits& limits = {});

[[nodiscard]] SceneFluidMimeticTraceFlowPrediction
sampleSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidMimeticTraceFlowLimits& limits = {});

void validateSceneFluidMimeticTraceFlowPredictionIntegrity(
    const SceneFluidMimeticTraceFlowPrediction& prediction);

void validateSceneFluidMimeticTraceFlowContinuationIntegrity(
    const SceneFluidMimeticTraceFlowContinuation& continuation);

void validateSceneFluidMimeticRegionTransportFlowPredictionIntegrity(
    const SceneFluidMimeticRegionTransportFlowPrediction& prediction);

void validateSceneFluidMimeticTraceFlowContinuation(
    const SceneFluidMimeticTraceFlowContinuation& continuation,
    const SceneFluidMimeticCorrectedTraceFlow& previousCorrectedFlow,
    const SceneFluidMimeticTraceFlowPrediction& previousBulkBaseline,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline);

void validateSceneFluidMimeticRegionTransportFlowPrediction(
    const SceneFluidMimeticRegionTransportFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidMimeticCorrectedTraceFlow& sourceCorrectedFlow,
    const SceneFluidRegionTransport& regionTransport,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline);

void validateSceneFluidMimeticTraceFlowPrediction(
    const SceneFluidMimeticTraceFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond);

void validateSceneFluidMimeticTraceFlowPrediction(
    const SceneFluidMimeticTraceFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange);

} // namespace simwing::fsi
