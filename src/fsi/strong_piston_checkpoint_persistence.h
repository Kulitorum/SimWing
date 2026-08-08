#pragma once

#include "piston_case.h"
#include "structure_checkpoint_persistence.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    strongPistonCheckpointProtocolVersion = 1;

struct StrongPistonCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 256u * 1024u * 1024u;
    StructureCheckpointPersistenceLimits structure;
    fluid::MovingInterfaceFluidCheckpointLimits fluid;
};

enum class StrongPistonCheckpointPersistenceErrorCode {
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

struct StrongPistonCheckpointPersistenceError {
    StrongPistonCheckpointPersistenceErrorCode code =
        StrongPistonCheckpointPersistenceErrorCode::None;
    std::string message;
};

// Deterministic bounded envelope around the existing Structure and accepted
// moving-interface fluid codecs. A fresh canonical owner validates their
// coupled identity and velocity closure before bytes or decoded state commit.
[[nodiscard]] bool serializeStrongPistonCheckpoint(
    const StrongCoupledPistonCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    StrongPistonCheckpointPersistenceError* error = nullptr,
    const StrongPistonCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeStrongPistonCheckpoint(
    std::span<const std::uint8_t> bytes,
    StrongCoupledPistonCheckpoint& checkpoint,
    StrongPistonCheckpointPersistenceError* error = nullptr,
    const StrongPistonCheckpointPersistenceLimits& limits = {});

} // namespace simwing::fsi
