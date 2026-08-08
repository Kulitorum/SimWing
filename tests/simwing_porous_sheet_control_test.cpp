#include "porous_sheet_control.h"

#include <cstdio>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

std::vector<std::uint8_t> serializedFrame(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "porous-sheet control frame serializes");
    return bytes;
}

void testCheckpointAndContinuation() {
    fsi::CoupledPorousSheetCase simulation;
    std::uint64_t publishedFrameCount = 0;
    std::optional<viewer::DiagnosticFrame> lastPublishedFrame;
    std::optional<fsi::CoupledPorousSheetCheckpoint> savedCheckpoint;
    fsi::PorousSheetControlHooks hooks;
    hooks.publishFrame = [&](const viewer::DiagnosticFrame& frame,
                             std::string&) {
        ++publishedFrameCount;
        lastPublishedFrame = frame;
        return true;
    };
    hooks.writeCheckpoint = [&savedCheckpoint](
                                const fsi::CoupledPorousSheetCheckpoint& value,
                                std::string&) {
        savedCheckpoint = value;
        return true;
    };
    fsi::PorousSheetControlSession session(simulation, std::move(hooks));

    const auto ready = session.readyResponse();
    check(ready.kind == fsi::WorkerControlResponseKind::Ready
              && ready.requestId == 0
              && ready.acceptedStepCount == 0
              && ready.simulationTimeSeconds == 0.0,
          "porous-sheet control reports its initial safe point");

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 401, 3},
              response, &protocolError)
              && !protocolError
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 401
              && response.acceptedStepCount == 3
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 3
              && publishedFrameCount == 3
              && lastPublishedFrame.has_value()
              && lastPublishedFrame->step == 3,
          "porous-sheet control publishes each accepted coupled step");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 402, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.requestId == 402
              && response.acceptedStepCount == 3
              && savedCheckpoint.has_value()
              && savedCheckpoint->acceptedStepCount == 3
              && savedCheckpoint->simulationTimeSeconds
                  == simulation.simulationTimeSeconds(),
          "porous-sheet control delegates the complete coupled checkpoint");

    fsi::CoupledPorousSheetCase rebuilt;
    rebuilt.restore(*savedCheckpoint);
    const auto expectedNext = rebuilt.advance();
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 403, 1},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.acceptedStepCount == 4
              && response.producedFrameCount == 1
              && publishedFrameCount == 4
              && serializedFrame(*lastPublishedFrame)
                  == serializedFrame(expectedNext),
          "delegated porous-sheet checkpoint reproduces the controlled next frame");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 404, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.requestId == 404
              && response.acceptedStepCount == 4
              && session.stopped(),
          "porous-sheet control stops at the final accepted safe point");
    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 405, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && response.acceptedStepCount == 4
              && simulation.acceptedStepCount() == 4,
          "stopped porous-sheet control rejects later checkpoint mutation");
}

} // namespace

int main() {
    testCheckpointAndContinuation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d porous-sheet control check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all porous-sheet control checks passed");
    return 0;
}
