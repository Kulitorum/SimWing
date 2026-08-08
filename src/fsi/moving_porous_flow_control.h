#pragma once

#include "moving_porous_flow_case.h"
#include "worker_control_session.h"

#include <functional>
#include <string>

namespace simwing::fsi {

using MovingPorousFlowFrameSink = WorkerControlFrameSink;
using MovingPorousFlowCheckpointSink = std::function<bool(
    const MovingPorousFlowCaseCheckpoint& checkpoint,
    std::string& error)>;

struct MovingPorousFlowControlHooks {
    MovingPorousFlowFrameSink publishFrame;
    MovingPorousFlowCheckpointSink writeCheckpoint;
};

// Typed safe-point binding for the prescribed moving porous full-flow worker.
// Accepted frame publication and immutable checkpoint persistence remain host
// callbacks and cannot roll back already committed solver state.
class MovingPorousFlowControlSession final : public WorkerControlSession {
public:
    explicit MovingPorousFlowControlSession(
        MovingPorousFlowCase& simulation,
        MovingPorousFlowControlHooks hooks = {});
};

} // namespace simwing::fsi
