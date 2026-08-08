#pragma once

#include "open_piston_case.h"
#include "worker_control_session.h"

#include <functional>
#include <string>

namespace simwing::fsi {

using OpenPistonFrameSink = WorkerControlFrameSink;
using OpenPistonCheckpointSink = std::function<bool(
    const OpenPistonCaseCheckpoint& checkpoint,
    std::string& error)>;

struct OpenPistonControlHooks {
    OpenPistonFrameSink publishFrame;
    OpenPistonCheckpointSink writeCheckpoint;
};

// Typed binding of the complete open-piston structure/fluid worker to the
// shared synchronous safe-point control state machine.
class OpenPistonControlSession final : public WorkerControlSession {
public:
    explicit OpenPistonControlSession(
        OpenPistonCase& simulation,
        OpenPistonControlHooks hooks = {});
};

} // namespace simwing::fsi
