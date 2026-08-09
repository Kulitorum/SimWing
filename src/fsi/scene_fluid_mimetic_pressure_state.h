#pragma once

#include "scene_fluid_mimetic_pressure_solve.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticPressureStateVersion = 1;

struct SceneFluidMimeticPressureStateLimits {
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumReducedTraces = 200'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticAcceptedControlPressure {
    std::size_t controlCellIndex = 0;
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    double pressurePascals = 0.0;

    bool operator==(
        const SceneFluidMimeticAcceptedControlPressure&) const = default;
};

struct SceneFluidMimeticAcceptedTracePressure {
    std::size_t reducedTraceIndex = 0;
    std::size_t fullTraceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::size_t componentIndex = 0;
    bool isGauge = false;
    double pressurePascals = 0.0;

    bool operator==(
        const SceneFluidMimeticAcceptedTracePressure&) const = default;
};

// Immutable accepted endpoint for mimetic pressure continuation. Capture is
// allowed only from a source-bound atomic solve whose reduced field exactly
// reconstructs the published full traces and cell scalars. Material-wall
// traces are derivable and are not persisted; shared gauges must already be
// exactly zero. This product owns no solver mutation and is safe to retain
// across a failed subsequent attempt.
struct SceneFluidMimeticPressureState {
    std::uint32_t version = sceneFluidMimeticPressureStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t pressureSourceFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t componentCount = 0;
    double maximumAbsoluteControlPressurePascals = 0.0;
    double maximumAbsoluteTracePressurePascals = 0.0;
    std::vector<SceneFluidMimeticAcceptedControlPressure> controls;
    std::vector<SceneFluidMimeticAcceptedTracePressure> traces;

    bool operator==(const SceneFluidMimeticPressureState&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureState
captureSceneFluidMimeticPressureState(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticPressureSolveResult& acceptedPressure,
    const SceneFluidMimeticPressureStateLimits& limits = {});

void validateSceneFluidMimeticPressureStateIntegrity(
    const SceneFluidMimeticPressureState& state);

void validateSceneFluidMimeticPressureState(
    const SceneFluidMimeticPressureState& state,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem);

} // namespace simwing::fsi
