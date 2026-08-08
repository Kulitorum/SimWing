#pragma once

#include "viewer_protocol.h"
#include "worker_control_protocol.h"

#include <cstdint>
#include <functional>
#include <string>

namespace simwing::fsi {

using WorkerControlAdvance =
    std::function<viewer::DiagnosticFrame()>;
using WorkerControlFrameSink = std::function<bool(
    const viewer::DiagnosticFrame& frame,
    std::string& error)>;
using WorkerControlCheckpointAction =
    std::function<bool(std::string& error)>;
using WorkerControlStepQuery = std::function<std::uint64_t()>;
using WorkerControlTimeQuery = std::function<double()>;

struct WorkerControlSessionHooks {
    WorkerControlAdvance advance;
    WorkerControlFrameSink publishFrame;
    WorkerControlCheckpointAction writeCheckpoint;
    WorkerControlStepQuery acceptedStepCount;
    WorkerControlTimeQuery simulationTimeSeconds;
};

// Case-neutral, synchronous execution of already-decoded control commands.
// The owning case adapter supplies numerical advance/checkpoint callbacks and
// absolute accepted state queries. Output failures never roll back an accepted
// numerical step; malformed commands fail before any callback is invoked.
class WorkerControlSession {
public:
    explicit WorkerControlSession(WorkerControlSessionHooks hooks);

    WorkerControlSession(const WorkerControlSession&) = delete;
    WorkerControlSession& operator=(const WorkerControlSession&) = delete;
    WorkerControlSession(WorkerControlSession&&) = delete;
    WorkerControlSession& operator=(WorkerControlSession&&) = delete;
    virtual ~WorkerControlSession() = default;

    [[nodiscard]] WorkerControlResponse readyResponse() const;

    [[nodiscard]] bool execute(
        const WorkerControlCommand& command,
        WorkerControlResponse& response,
        WorkerControlProtocolError* protocolError = nullptr);

    [[nodiscard]] bool stopped() const noexcept;

private:
    [[nodiscard]] WorkerControlResponse response(
        WorkerControlResponseKind kind,
        std::uint64_t requestId) const;
    [[nodiscard]] WorkerControlResponse failureResponse(
        std::uint64_t requestId,
        WorkerControlFailureCode failureCode,
        std::string message) const;

    WorkerControlSessionHooks hooks_;
    bool stopped_ = false;
};

} // namespace simwing::fsi
