#pragma once

#include "fluid/mimetic_local_cell.h"
#include "scene_fluid_mimetic_control_cell.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticTraceSystemVersion = 1;
inline constexpr std::size_t invalidSceneFluidMimeticTraceIndex =
    std::numeric_limits<std::size_t>::max();

struct SceneFluidMimeticTraceSystemSettings {
    fluid::MimeticLocalCellSettings localCell;

    bool operator==(
        const SceneFluidMimeticTraceSystemSettings&) const = default;
};

struct SceneFluidMimeticTraceSystemLimits {
    std::size_t maximumTraces = 200'000'000;
    std::size_t maximumIncidences = 400'000'000;
    std::size_t maximumLocalOperators = 50'000'000;
    std::size_t maximumLocalOperatorBytes =
        4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticTrace {
    std::size_t traceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::uint64_t sourceStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t firstIncidence = 0;
    std::size_t incidenceCount = 0;
    double operatorDiagonal = 0.0;
    bool isGauge = false;

    bool operator==(const SceneFluidMimeticTrace&) const = default;
};

struct SceneFluidMimeticTraceIncidence {
    std::size_t incidenceIndex = 0;
    std::size_t traceIndex = 0;
    std::size_t halfFaceIndex = 0;
    std::size_t controlCellIndex = 0;
    std::size_t localHalfFaceIndex = 0;

    bool operator==(
        const SceneFluidMimeticTraceIncidence&) const = default;
};

// Matrix-free mixed-hybrid trace topology. Shared Cartesian and opening
// half-faces identify one two-sided trace unknown; every impermeable material
// half-face owns one local wall-trace unknown whose equation is zero normal
// flux. Each control cell retains the compact local SPD factorization. Cell
// scalars are eliminated by exact integrated conservation during application,
// so the global positive-semidefinite operator is minus the sum of incident
// outward fluxes. One deterministic trace gauge is retained per authored
// pressure component. This audit product does not replace the production graph
// pressure operator.
struct SceneFluidMimeticTraceSystem {
    std::uint32_t version = sceneFluidMimeticTraceSystemVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidMimeticTraceSystemSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t localOperatorStorageBytes = 0;
    std::size_t componentCount = 0;
    std::size_t sharedTraceCount = 0;
    std::size_t materialWallTraceCount = 0;
    double minimumPositiveOperatorDiagonal = 0.0;
    double maximumOperatorDiagonal = 0.0;
    std::vector<SceneFluidMimeticTrace> traces;
    std::vector<SceneFluidMimeticTraceIncidence> incidences;
    std::vector<std::size_t> halfFaceTraceIndices;
    std::vector<std::size_t> componentGaugeTraceIndices;
    std::vector<fluid::MimeticLocalCellOperator> localOperators;

    bool operator==(const SceneFluidMimeticTraceSystem&) const = default;
};

struct SceneFluidMimeticTraceEvaluation {
    std::vector<double> cellScalars;
    std::vector<double> halfFaceIntegratedOutwardFluxes;
    std::vector<double> traceIntegratedOutwardFluxSums;
    double maximumCellConservationResidual = 0.0;
    double maximumTraceFluxImbalance = 0.0;
};

[[nodiscard]] SceneFluidMimeticTraceSystem
buildSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystemSettings& settings = {},
    const SceneFluidMimeticTraceSystemLimits& limits = {});

void validateSceneFluidMimeticTraceSystemIntegrity(
    const SceneFluidMimeticTraceSystem& system);

void validateSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    const SceneFluidMimeticControlCellSet& controlCells);

[[nodiscard]] SceneFluidMimeticTraceEvaluation
evaluateSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    std::span<const double> traceScalars,
    std::span<const double> integratedCellSources);

// Applies the positive-semidefinite trace operator with zero cell source.
[[nodiscard]] std::vector<double> applySceneFluidMimeticTraceOperator(
    const SceneFluidMimeticTraceSystem& system,
    std::span<const double> traceScalars);

// Builds the trace right-hand side induced by integrated cell sources.
[[nodiscard]] std::vector<double> buildSceneFluidMimeticTraceRightHandSide(
    const SceneFluidMimeticTraceSystem& system,
    std::span<const double> integratedCellSources);

} // namespace simwing::fsi
