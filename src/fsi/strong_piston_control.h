#pragma once

#include "piston_case.h"
#include "worker_control_session.h"

#include <functional>
#include <string>

namespace simwing::fsi {

using StrongPistonFrameSink = WorkerControlFrameSink;
using StrongPistonCheckpointSink = std::function<bool(
    const StrongCoupledPistonCheckpoint& checkpoint,
    std::string& error)>;

struct StrongPistonControlHooks {
    StrongPistonFrameSink publishFrame;
    StrongPistonCheckpointSink writeCheckpoint;
};

// Typed safe-point binding for the accepted-only strong-coupling worker. Each
// Advance command completes whole macro-steps; rejected attempts and active
// fixed-point iterations are never externally observable or checkpointable.
class StrongPistonControlSession final : public WorkerControlSession {
public:
    explicit StrongPistonControlSession(
        StrongCoupledPistonWorkerCase& simulation,
        StrongPistonControlHooks hooks = {});
};

} // namespace simwing::fsi
