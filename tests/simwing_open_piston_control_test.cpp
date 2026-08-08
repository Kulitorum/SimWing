#include "open_piston_control.h"

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
          "open-piston control frame serializes");
    return bytes;
}

void testRebaseCheckpointAndContinuation() {
    fsi::OpenPistonCase simulation;
    std::uint64_t publishedFrameCount = 0;
    std::optional<viewer::DiagnosticFrame> lastPublishedFrame;
    std::optional<fsi::OpenPistonCaseCheckpoint> savedCheckpoint;
    fsi::OpenPistonControlHooks hooks;
    hooks.publishFrame = [&](const viewer::DiagnosticFrame& frame,
                             std::string&) {
        ++publishedFrameCount;
        lastPublishedFrame = frame;
        return true;
    };
    hooks.writeCheckpoint = [&savedCheckpoint](
                                const fsi::OpenPistonCaseCheckpoint& value,
                                std::string&) {
        savedCheckpoint = value;
        return true;
    };
    fsi::OpenPistonControlSession session(simulation, std::move(hooks));

    const auto ready = session.readyResponse();
    check(ready.kind == fsi::WorkerControlResponseKind::Ready
              && ready.requestId == 0
              && ready.acceptedStepCount == 0
              && ready.simulationTimeSeconds == 0.0,
          "open-piston control reports its initial safe point");

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 301, 1200},
              response, &protocolError)
              && !protocolError
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 301
              && response.acceptedStepCount == 1200
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 1200
              && publishedFrameCount == 1200
              && simulation.topologyRebaseCount() == 1
              && simulation.movingPlaneCoordinate() == 7,
          "open-piston control crosses and publishes the first rebase");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 302, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.requestId == 302
              && response.acceptedStepCount == 1200
              && savedCheckpoint.has_value()
              && savedCheckpoint->acceptedStepCount == 1200
              && savedCheckpoint->topologyRebaseCount == 1
              && savedCheckpoint->movingPlaneCoordinate == 7,
          "open-piston control delegates the complete rebased checkpoint");

    fsi::OpenPistonCase rebuilt;
    rebuilt.restore(*savedCheckpoint);
    const auto expectedNext = rebuilt.advance();
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 303, 1},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.acceptedStepCount == 1201
              && response.producedFrameCount == 1
              && publishedFrameCount == 1201
              && lastPublishedFrame.has_value()
              && serializedFrame(*lastPublishedFrame)
                  == serializedFrame(expectedNext),
          "delegated rebased checkpoint reproduces the controlled next frame");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 304, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.requestId == 304
              && response.acceptedStepCount == 1201
              && session.stopped(),
          "open-piston control stops at the final accepted safe point");
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 305, 1},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && response.acceptedStepCount == 1201
              && simulation.acceptedStepCount() == 1201,
          "stopped open-piston control rejects later numerical mutation");
}

} // namespace

int main() {
    testRebaseCheckpointAndContinuation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d open-piston control check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all open-piston control checks passed");
    return 0;
}
