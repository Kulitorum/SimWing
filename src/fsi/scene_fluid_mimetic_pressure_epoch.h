#pragma once

#include "scene_fluid_mimetic_pressure_sampling.h"
#include "scene_fluid_mimetic_pressure_warm_start.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticPressureEpochVersion = 1;

enum class SceneFluidMimeticPressureEpochFailureStage : std::uint8_t {
    None = 0,
    PressureSolve = 1,
};

struct SceneFluidMimeticPressureEpochLimits {
    std::size_t maximumBootstrapWarmStartBytes =
        2048ULL * 1024ULL * 1024ULL;
    SceneFluidMimeticPressureStateLimits state;
    SceneFluidPressureSamplingLimits sampling;
    SceneFluidMimeticPressureWarmStartLimits warmStart;
};

struct SceneFluidMimeticPressureEpochDiagnostics {
    bool accepted = false;
    bool usedConsecutiveWarmStart = false;
    SceneFluidMimeticPressureEpochFailureStage failureStage =
        SceneFluidMimeticPressureEpochFailureStage::None;
    SceneFluidMimeticPressureSolveDiagnostics pressureSolve;

    bool operator==(
        const SceneFluidMimeticPressureEpochDiagnostics&) const = default;
};

// Atomic source-to-surface-pressure transaction. A rejected pressure solve
// publishes diagnostics only. An accepted solve is independently captured as
// immutable pressure state and sampled onto the exact material quadrature
// before either payload is returned. The consecutive overload derives its
// read-only reduced warm field from the prior accepted state and topology
// transition; bootstrap uses an explicit bounded zero field.
struct SceneFluidMimeticPressureEpochResult {
    std::uint32_t version = sceneFluidMimeticPressureEpochVersion;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t pressureSourceFingerprint = 0;
    std::uint64_t warmStartFingerprint = 0;
    std::uint64_t topologyTransitionFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidMimeticPressureEpochDiagnostics diagnostics;
    SceneFluidMimeticPressureState acceptedPressureState;
    SceneFluidMimeticPressureSampleSet acceptedPressureSamples;

    bool operator==(
        const SceneFluidMimeticPressureEpochResult&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureEpochResult
acceptSceneFluidMimeticPressureEpoch(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticTraceSolveSettings& solveSettings = {},
    const SceneFluidMimeticPressureEpochLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureEpochResult
acceptSceneFluidMimeticPressureEpoch(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticTraceSolveSettings& solveSettings = {},
    const SceneFluidMimeticPressureEpochLimits& limits = {});

void validateSceneFluidMimeticPressureEpochResultIntegrity(
    const SceneFluidMimeticPressureEpochResult& result);

void validateSceneFluidMimeticPressureEpochResult(
    const SceneFluidMimeticPressureEpochResult& result,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources);

} // namespace simwing::fsi
