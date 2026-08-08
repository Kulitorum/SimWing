#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint16_t workerControlProtocolVersion = 1;

struct WorkerControlProtocolLimits {
    std::uint64_t maximumMessageBytes = 64ULL * 1024ULL;
    std::uint64_t maximumAdvanceSteps = 10'000'000;
    std::uint32_t maximumErrorMessageBytes = 4096;
};

enum class WorkerControlProtocolErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
};

struct WorkerControlProtocolError {
    WorkerControlProtocolErrorCode code =
        WorkerControlProtocolErrorCode::None;
    std::string message;

    [[nodiscard]] explicit operator bool() const noexcept {
        return code != WorkerControlProtocolErrorCode::None;
    }
};

enum class WorkerControlCommandKind : std::uint8_t {
    Advance = 1,
    Checkpoint = 2,
    Stop = 3,
};

struct WorkerControlCommand {
    WorkerControlCommandKind kind = WorkerControlCommandKind::Advance;
    std::uint64_t requestId = 0;
    std::uint64_t advanceStepCount = 0;

    bool operator==(const WorkerControlCommand&) const = default;
};

enum class WorkerControlResponseKind : std::uint8_t {
    Ready = 1,
    Advanced = 2,
    Checkpointed = 3,
    Stopped = 4,
    Error = 5,
};

enum class WorkerControlFailureCode : std::uint16_t {
    None = 0,
    InvalidCommand = 1,
    NumericalFailure = 2,
    CheckpointFailure = 3,
    InternalFailure = 4,
};

struct WorkerControlResponse {
    WorkerControlResponseKind kind = WorkerControlResponseKind::Ready;
    std::uint64_t requestId = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::uint64_t producedFrameCount = 0;
    WorkerControlFailureCode failureCode = WorkerControlFailureCode::None;
    std::string errorMessage;

    bool operator==(const WorkerControlResponse&) const = default;
};

[[nodiscard]] bool validateWorkerControlCommand(
    const WorkerControlCommand& command,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool validateWorkerControlResponse(
    const WorkerControlResponse& response,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool serializeWorkerControlCommand(
    const WorkerControlCommand& command,
    std::vector<std::uint8_t>& bytes,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool deserializeWorkerControlCommand(
    std::span<const std::uint8_t> bytes,
    WorkerControlCommand& command,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool serializeWorkerControlResponse(
    const WorkerControlResponse& response,
    std::vector<std::uint8_t>& bytes,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool deserializeWorkerControlResponse(
    std::span<const std::uint8_t> bytes,
    WorkerControlResponse& response,
    WorkerControlProtocolError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

} // namespace simwing::fsi
