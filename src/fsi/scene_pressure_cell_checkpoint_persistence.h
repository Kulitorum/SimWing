#pragma once

#include "scene_pressure_cell_case.h"
#include "scene_fluid_mimetic_pressure_state_persistence.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    scenePressureCellCheckpointProtocolVersion = 10;

struct ScenePressureCellCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 256u * 1024u * 1024u;
    std::size_t maximumControlVolumes = 5'000'000;
    std::size_t maximumLinks = 10'000'000;
    std::size_t maximumSolveComponents = 1'000'000;
    std::size_t maximumProjectionStorageBytes =
        192u * 1024u * 1024u;
    std::size_t maximumMomentumStorageBytes =
        192u * 1024u * 1024u;
    std::size_t maximumWallTractions = 10'000'000;
    std::size_t maximumWallTractionStorageBytes =
        192u * 1024u * 1024u;
    SceneFluidMimeticPressureStatePersistenceLimits mimeticPressureState;
    StructureCheckpointPersistenceLimits structure;
};

enum class ScenePressureCellCheckpointPersistenceErrorCode {
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

struct ScenePressureCellCheckpointPersistenceError {
    ScenePressureCellCheckpointPersistenceErrorCode code =
        ScenePressureCellCheckpointPersistenceErrorCode::None;
    std::string message;
};

// Canonical pressure-cell persistence. Identity and solver settings are taken
// from a freshly rebuilt case rather than trusted from the wire. The nested
// Structure codec plus pressure-projection, material-wall traction, and
// region-momentum and compact SWMP integrity validators all run before
// serialization or
// decoded state is published. The accepted bulk MAC
// predictor is deterministically derived from that projection after restore;
// the private per-step bulk pressure is transient and neither field is
// duplicated in the envelope. An initial checkpoint reconstructs the
// canonical prescribed wind. Accepted checkpoints also persist the immutable
// region-momentum state used by the next conservative transport step. Audit
// checkpoints persist only accepted pressure rows; trusted control/trace
// topology is rebuilt from the Structure payload before SWMP decode.
[[nodiscard]] bool serializeScenePressureCellCheckpoint(
    const ScenePressureCellCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    ScenePressureCellCheckpointPersistenceError* error = nullptr,
    const ScenePressureCellCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeScenePressureCellCheckpoint(
    std::span<const std::uint8_t> bytes,
    ScenePressureCellCheckpoint& checkpoint,
    ScenePressureCellCheckpointPersistenceError* error = nullptr,
    const ScenePressureCellCheckpointPersistenceLimits& limits = {});

} // namespace simwing::fsi
