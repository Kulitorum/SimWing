#pragma once

#include "porous_sheet_case.h"
#include "worker_control_session.h"

#include <functional>
#include <string>

namespace simwing::fsi {

using PorousSheetFrameSink = WorkerControlFrameSink;
using PorousSheetCheckpointSink = std::function<bool(
    const CoupledPorousSheetCheckpoint& checkpoint,
    std::string& error)>;

struct PorousSheetControlHooks {
    PorousSheetFrameSink publishFrame;
    PorousSheetCheckpointSink writeCheckpoint;
};

// Typed binding of the coupled porous-sheet structure/fluid worker to the
// shared synchronous safe-point control state machine.
class PorousSheetControlSession final : public WorkerControlSession {
public:
    explicit PorousSheetControlSession(
        CoupledPorousSheetCase& simulation,
        PorousSheetControlHooks hooks = {});
};

} // namespace simwing::fsi
