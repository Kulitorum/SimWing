#pragma once

#include "moving_porous_flow_case.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t
    movingPorousFlowCaseCheckpointProtocolVersion = 1;

struct MovingPorousFlowCaseCheckpointPersistenceLimits {
    std::size_t maximumEncodedBytes = 64u * 1024u * 1024u;
    std::size_t maximumScalarSamples = 5'000'000;
    std::size_t maximumPressureJumpFaces = 1'000'000;
    std::uint64_t maximumReplaySteps = 10'000;
};

enum class MovingPorousFlowCaseCheckpointPersistenceErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
};

struct MovingPorousFlowCaseCheckpointPersistenceError {
    MovingPorousFlowCaseCheckpointPersistenceErrorCode code =
        MovingPorousFlowCaseCheckpointPersistenceErrorCode::None;
    std::string message;

    [[nodiscard]] explicit operator bool() const noexcept {
        return code
            != MovingPorousFlowCaseCheckpointPersistenceErrorCode::None;
    }
};

// Bounded deterministic little-endian persistence for the canonical moving
// porous worker. Fields and sharp crossings are stored explicitly. Decode
// replays the bounded canonical history and publishes only when every stored
// field and public topology/kinematic epoch matches bit-for-bit.
[[nodiscard]] bool serializeMovingPorousFlowCaseCheckpoint(
    const MovingPorousFlowCaseCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    MovingPorousFlowCaseCheckpointPersistenceError* error = nullptr,
    const MovingPorousFlowCaseCheckpointPersistenceLimits& limits = {});

[[nodiscard]] bool deserializeMovingPorousFlowCaseCheckpoint(
    std::span<const std::uint8_t> bytes,
    MovingPorousFlowCaseCheckpoint& checkpoint,
    MovingPorousFlowCaseCheckpointPersistenceError* error = nullptr,
    const MovingPorousFlowCaseCheckpointPersistenceLimits& limits = {});

} // namespace simwing::fsi
