#pragma once

#include "scene_fluid_mimetic_pressure_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    sceneFluidMimeticPressureStateProtocolVersion = 1;

struct SceneFluidMimeticPressureStatePersistenceLimits {
    std::size_t maximumEncodedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumReducedTraces = 200'000'000;
};

enum class SceneFluidMimeticPressureStatePersistenceErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
    TopologyMismatch,
};

struct SceneFluidMimeticPressureStatePersistenceError {
    SceneFluidMimeticPressureStatePersistenceErrorCode code =
        SceneFluidMimeticPressureStatePersistenceErrorCode::None;
    std::string message;
};

// Deterministic bounded little-endian persistence for one accepted pressure
// endpoint. Decode receives the trusted rebuilt control/full/condensed
// topology and publishes only after the state's fingerprint, stable rows,
// gauges, epoch, and complete topology provenance validate against it.
[[nodiscard]] bool serializeSceneFluidMimeticPressureState(
    const SceneFluidMimeticPressureState& state,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    std::vector<std::uint8_t>& bytes,
    SceneFluidMimeticPressureStatePersistenceError* error = nullptr,
    const SceneFluidMimeticPressureStatePersistenceLimits& limits = {});

[[nodiscard]] bool deserializeSceneFluidMimeticPressureState(
    std::span<const std::uint8_t> bytes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    SceneFluidMimeticPressureState& state,
    SceneFluidMimeticPressureStatePersistenceError* error = nullptr,
    const SceneFluidMimeticPressureStatePersistenceLimits& limits = {});

} // namespace simwing::fsi
