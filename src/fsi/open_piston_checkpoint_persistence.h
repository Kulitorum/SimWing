#pragma once

#include "open_piston_case.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t openPistonCheckpointProtocolVersion = 1;

struct OpenPistonCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 256u * 1024u * 1024u;
    std::size_t maximumCutSurfaceFaces = 1'000'000;
    StructureCheckpointPersistenceLimits structure;
    fluid::MovingInterfaceFluidCheckpointLimits fluid;
};

enum class OpenPistonCheckpointPersistenceErrorCode {
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

struct OpenPistonCheckpointPersistenceError {
    OpenPistonCheckpointPersistenceErrorCode code =
        OpenPistonCheckpointPersistenceErrorCode::None;
    std::string message;
};

[[nodiscard]] bool serializeOpenPistonCheckpoint(
    const OpenPistonCaseCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    OpenPistonCheckpointPersistenceError* error = nullptr,
    const OpenPistonCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeOpenPistonCheckpoint(
    std::span<const std::uint8_t> bytes,
    OpenPistonCaseCheckpoint& checkpoint,
    OpenPistonCheckpointPersistenceError* error = nullptr,
    const OpenPistonCheckpointPersistenceLimits& limits = {});

} // namespace simwing::fsi
