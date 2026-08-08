#pragma once

#include "structure.h"

#include <softwing/checkpoint_persistence.h>
#include <softwing/suspension_checkpoint_persistence.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t structureCheckpointProtocolVersion = 1;

struct StructureCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 192u * 1024u * 1024u;
    std::size_t maximumNodes = 5'000'000;
    softwing::SoftBodyCheckpointPersistenceLimits body;
    softwing::SuspensionCheckpointPersistenceLimits suspension;
};

enum class StructureCheckpointPersistenceErrorCode {
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

struct StructureCheckpointPersistenceError {
    StructureCheckpointPersistenceErrorCode code =
        StructureCheckpointPersistenceErrorCode::None;
    std::string message;
};

// The owner supplies the trusted definition. Serialization and decoding both
// restore through an equivalent rebuilt Structure before publishing bytes or
// output, preserving the adapter's existing transactional semantic validator.
[[nodiscard]] bool serializeStructureCheckpoint(
    const Structure& owner,
    const StructureCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    StructureCheckpointPersistenceError* error = nullptr,
    const StructureCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeStructureCheckpoint(
    std::span<const std::uint8_t> bytes,
    const Structure& owner,
    StructureCheckpoint& checkpoint,
    StructureCheckpointPersistenceError* error = nullptr,
    const StructureCheckpointPersistenceLimits& limits = {});

} // namespace simwing::fsi
