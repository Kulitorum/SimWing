#pragma once

#include "softwing/suspension.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace softwing {

inline constexpr std::uint16_t suspensionCheckpointProtocolVersion = 1;

struct SuspensionCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 64u * 1024u * 1024u;
    std::size_t maximumRecords = 1'000'000;
    std::size_t maximumStringBytes = 1u * 1024u * 1024u;
    std::size_t maximumTotalTextBytes = 8u * 1024u * 1024u;
};

enum class SuspensionCheckpointPersistenceErrorCode {
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

struct SuspensionCheckpointPersistenceError {
    SuspensionCheckpointPersistenceErrorCode code =
        SuspensionCheckpointPersistenceErrorCode::None;
    std::string message;
};

[[nodiscard]] bool serializeSuspensionCheckpoint(
    const SuspensionCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    SuspensionCheckpointPersistenceError* error = nullptr,
    const SuspensionCheckpointPersistenceLimits& limits = {});

// Counts, stable identities, and the topology fingerprint must match the
// supplied checkpoint from an equivalent rebuilt SuspensionSystem. The live
// owner's restore() remains the final semantic validator before state commit.
[[nodiscard]] bool deserializeSuspensionCheckpoint(
    std::span<const std::uint8_t> bytes,
    const SuspensionCheckpoint& topologyTemplate,
    SuspensionCheckpoint& checkpoint,
    SuspensionCheckpointPersistenceError* error = nullptr,
    const SuspensionCheckpointPersistenceLimits& limits = {});

} // namespace softwing
