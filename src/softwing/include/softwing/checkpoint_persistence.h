#pragma once

#include "softwing/soft_body.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace softwing {

inline constexpr std::uint16_t softBodyCheckpointProtocolVersion = 1;

struct SoftBodyCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 64u * 1024u * 1024u;
    std::size_t maximumNodes = 5'000'000;
    std::size_t maximumTriangles = 5'000'000;
    std::size_t maximumConstraints = 10'000'000;
    std::size_t maximumMembranes = 5'000'000;
    std::size_t maximumDihedrals = 10'000'000;
    std::size_t maximumContactPairs = 1'000'000;
    std::size_t maximumContactMultipliers = 5'000'000;
    std::size_t maximumContactRecords = 1'000'000;
    std::size_t maximumAuditKeys = 5'000'000;
};

enum class SoftBodyCheckpointPersistenceErrorCode {
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

struct SoftBodyCheckpointPersistenceError {
    SoftBodyCheckpointPersistenceErrorCode code =
        SoftBodyCheckpointPersistenceErrorCode::None;
    std::string message;
};

// The wire payload contains only committed mutable state. Deserialization
// overlays it onto an immutable topology template captured from an equivalent
// rebuilt body, so persisted bytes cannot redefine masses, connectivity,
// materials, contact registration, or any other solver topology.
[[nodiscard]] bool serializeSoftBodyCheckpoint(
    const SoftBodyCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    SoftBodyCheckpointPersistenceError* error = nullptr,
    const SoftBodyCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeSoftBodyCheckpoint(
    std::span<const std::uint8_t> bytes,
    const SoftBodyCheckpoint& topologyTemplate,
    SoftBodyCheckpoint& checkpoint,
    SoftBodyCheckpointPersistenceError* error = nullptr,
    const SoftBodyCheckpointPersistenceLimits& limits = {});

} // namespace softwing
