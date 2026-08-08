#pragma once

#include "worker_control_protocol.h"

#include <iosfwd>
#include <string>

namespace simwing::fsi {

enum class WorkerControlStreamResult {
    Message,
    EndOfStream,
    Error,
};

enum class WorkerControlStreamErrorCode {
    None,
    IoFailure,
    Truncated,
    LimitExceeded,
    ProtocolFailure,
};

struct WorkerControlStreamError {
    WorkerControlStreamErrorCode code = WorkerControlStreamErrorCode::None;
    WorkerControlProtocolErrorCode protocolCode =
        WorkerControlProtocolErrorCode::None;
    std::string message;

    [[nodiscard]] explicit operator bool() const noexcept {
        return code != WorkerControlStreamErrorCode::None;
    }
};

// Messages are self-framed by the protocol envelope's bounded payload length;
// the stream adds no host-endian or platform-sized prefix. A physical EOF is
// clean only before the first byte of a new envelope. All decoded outputs are
// transactional, and every successful write is flushed for pipe operation.
[[nodiscard]] WorkerControlStreamResult readWorkerControlCommand(
    std::istream& input,
    WorkerControlCommand& command,
    WorkerControlStreamError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] WorkerControlStreamResult readWorkerControlResponse(
    std::istream& input,
    WorkerControlResponse& response,
    WorkerControlStreamError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool writeWorkerControlCommand(
    std::ostream& output,
    const WorkerControlCommand& command,
    WorkerControlStreamError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

[[nodiscard]] bool writeWorkerControlResponse(
    std::ostream& output,
    const WorkerControlResponse& response,
    WorkerControlStreamError* error = nullptr,
    const WorkerControlProtocolLimits& limits = {});

} // namespace simwing::fsi
