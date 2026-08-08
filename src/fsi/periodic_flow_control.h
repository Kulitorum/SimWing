#pragma once

#include "periodic_flow_case.h"
#include "worker_control_session.h"

#include <functional>
#include <string>

namespace simwing::fsi {

// Called only after PeriodicFlowCase has committed one accepted step. Returning
// false reports an output failure at that safe point; it does not roll back the
// already accepted numerical state.
using PeriodicFlowFrameSink = WorkerControlFrameSink;

// Persistence remains outside the numerical adapter. The callback receives an
// immutable complete checkpoint while the worker is stopped at a safe point.
using PeriodicFlowCheckpointSink = std::function<bool(
    const PeriodicFlowCaseCheckpoint& checkpoint,
    std::string& error)>;

struct PeriodicFlowControlHooks {
    PeriodicFlowFrameSink publishFrame;
    PeriodicFlowCheckpointSink writeCheckpoint;
};

// Typed binding of PeriodicFlowCase to the shared synchronous safe-point
// control state machine.
class PeriodicFlowControlSession final : public WorkerControlSession {
public:
    explicit PeriodicFlowControlSession(
        PeriodicFlowCase& simulation,
        PeriodicFlowControlHooks hooks = {});
};

} // namespace simwing::fsi
