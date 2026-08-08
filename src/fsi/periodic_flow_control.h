#pragma once

#include "periodic_flow_case.h"
#include "worker_control_protocol.h"

#include <functional>
#include <string>

namespace simwing::fsi {

// Called only after PeriodicFlowCase has committed one accepted step. Returning
// false reports an output failure at that safe point; it does not roll back the
// already accepted numerical state.
using PeriodicFlowFrameSink = std::function<bool(
    const viewer::DiagnosticFrame& frame,
    std::string& error)>;

// Persistence remains outside the numerical adapter. The callback receives an
// immutable complete checkpoint while the worker is stopped at a safe point.
using PeriodicFlowCheckpointSink = std::function<bool(
    const PeriodicFlowCaseCheckpoint& checkpoint,
    std::string& error)>;

struct PeriodicFlowControlHooks {
    PeriodicFlowFrameSink publishFrame;
    PeriodicFlowCheckpointSink writeCheckpoint;
};

// Synchronous, transport-independent execution of decoded control commands.
// One owner thread calls execute() between solver steps. Advance commits and
// publishes each accepted step individually, checkpoint delegates persistence,
// and stop permanently prevents later mutation. Operational failures are valid
// Error responses; malformed commands fail transactionally without a response.
class PeriodicFlowControlSession final {
public:
    explicit PeriodicFlowControlSession(
        PeriodicFlowCase& simulation,
        PeriodicFlowControlHooks hooks = {});

    PeriodicFlowControlSession(const PeriodicFlowControlSession&) = delete;
    PeriodicFlowControlSession& operator=(
        const PeriodicFlowControlSession&) = delete;
    PeriodicFlowControlSession(PeriodicFlowControlSession&&) = delete;
    PeriodicFlowControlSession& operator=(
        PeriodicFlowControlSession&&) = delete;

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

    PeriodicFlowCase* simulation_ = nullptr;
    PeriodicFlowControlHooks hooks_;
    bool stopped_ = false;
};

} // namespace simwing::fsi
