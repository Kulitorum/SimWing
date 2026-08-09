#pragma once

#include "scene_fluid_mimetic_pressure_state.h"
#include "scene_fluid_pressure_topology_transition.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticPressureWarmStartVersion = 1;

struct SceneFluidMimeticPressureWarmStartLimits {
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumPreviousReducedTraces = 200'000'000;
    std::size_t maximumCurrentReducedTraces = 200'000'000;
    std::size_t maximumComponents = 50'000'000;
    std::size_t maximumWorkingBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Immutable solver initialization for one consecutive accepted grid epoch.
// Control pressures follow the shared topology-transition policy. Reduced
// traces with retained stable IDs keep their previous value before the whole
// current component is shifted to its deterministic zero gauge. A new trace
// starts at the arithmetic mean of its rebased endpoint control pressures.
// Retired traces disappear. This product cannot alter either accepted state.
struct SceneFluidMimeticPressureWarmStart {
    std::uint32_t version =
        sceneFluidMimeticPressureWarmStartVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourcePressureStateFingerprint = 0;
    std::uint64_t sourceTopologyTransitionFingerprint = 0;
    std::uint64_t previousMimeticControlCellFingerprint = 0;
    std::uint64_t previousFullTraceSystemFingerprint = 0;
    std::uint64_t previousCondensedTraceSystemFingerprint = 0;
    std::uint64_t currentMimeticControlCellFingerprint = 0;
    std::uint64_t currentFullTraceSystemFingerprint = 0;
    std::uint64_t currentCondensedTraceSystemFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;
    std::size_t componentCount = 0;
    std::size_t previousReducedTraceCount = 0;
    std::size_t currentReducedTraceCount = 0;
    std::size_t retainedTraceCount = 0;
    std::size_t appearedTraceCount = 0;
    std::size_t disappearedTraceCount = 0;
    double maximumAbsoluteGaugeShiftPascals = 0.0;
    double maximumAbsolutePressurePascals = 0.0;
    std::vector<double> componentGaugeShiftsPascals;
    std::vector<double> reducedTracePascals;

    bool operator==(
        const SceneFluidMimeticPressureWarmStart&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureWarmStart
buildSceneFluidMimeticPressureWarmStart(
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidMimeticPressureWarmStartLimits& limits = {});

void validateSceneFluidMimeticPressureWarmStartIntegrity(
    const SceneFluidMimeticPressureWarmStart& warmStart);

void validateSceneFluidMimeticPressureWarmStart(
    const SceneFluidMimeticPressureWarmStart& warmStart,
    const SceneFluidMimeticPressureState& previousState,
    const SceneFluidMimeticControlCellSet& previousControlCells,
    const SceneFluidMimeticTraceSystem& previousFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& previousCondensedSystem,
    const SceneFluidMimeticControlCellSet& currentControlCells,
    const SceneFluidMimeticTraceSystem& currentFullSystem,
    const SceneFluidMimeticCondensedTraceSystem& currentCondensedSystem,
    const SceneFluidPressureTopologyTransition& topologyTransition);

} // namespace simwing::fsi
