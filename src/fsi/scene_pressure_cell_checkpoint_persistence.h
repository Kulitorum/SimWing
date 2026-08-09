#pragma once

#include "scene_pressure_cell_case.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    scenePressureCellCheckpointProtocolVersion = 2;

struct ScenePressureCellCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 256u * 1024u * 1024u;
    std::size_t maximumControlVolumes = 5'000'000;
    std::size_t maximumLinks = 10'000'000;
    std::size_t maximumSolveComponents = 1'000'000;
    std::size_t maximumProjectionStorageBytes =
        192u * 1024u * 1024u;
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
// Structure codec and the pressure-projection integrity validator both run
// before serialization or decoded state is published. The bulk MAC predictor
// is deterministically derived from that accepted projection after restore;
// it is intentionally not duplicated in the envelope.
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
