#include "porous_sheet_control.h"

#include <utility>

namespace simwing::fsi {
namespace {

WorkerControlSessionHooks makeSessionHooks(
    CoupledPorousSheetCase& simulation,
    PorousSheetControlHooks hooks) {
    WorkerControlSessionHooks result;
    result.advance = [&simulation] { return simulation.advance(); };
    result.publishFrame = std::move(hooks.publishFrame);
    if (hooks.writeCheckpoint) {
        result.writeCheckpoint = [
            &simulation,
            sink = std::move(hooks.writeCheckpoint)](
                std::string& error) {
            return sink(simulation.checkpoint(), error);
        };
    }
    result.acceptedStepCount = [&simulation] {
        return simulation.acceptedStepCount();
    };
    result.simulationTimeSeconds = [&simulation] {
        return simulation.simulationTimeSeconds();
    };
    return result;
}

} // namespace

PorousSheetControlSession::PorousSheetControlSession(
    CoupledPorousSheetCase& simulation,
    PorousSheetControlHooks hooks)
    : WorkerControlSession(makeSessionHooks(
          simulation, std::move(hooks))) {}

} // namespace simwing::fsi
