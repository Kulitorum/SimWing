#include "moving_porous_flow_control.h"

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
          "moving porous control frame serializes");
    return bytes;
}

void checkValidResponse(
    const fsi::WorkerControlResponse& response,
    const char* message) {
    fsi::WorkerControlProtocolError error;
    check(fsi::validateWorkerControlResponse(response, &error) && !error,
          message);
}

void testAdvanceCheckpointAndReplay() {
    fsi::MovingPorousFlowCase simulation;
    std::vector<std::vector<std::uint8_t>> frames;
    std::optional<fsi::MovingPorousFlowCaseCheckpoint> saved;
    fsi::MovingPorousFlowControlHooks hooks;
    hooks.publishFrame = [&frames](
                             const viewer::DiagnosticFrame& frame,
                             std::string&) {
        frames.push_back(serializedFrame(frame));
        return true;
    };
    hooks.writeCheckpoint = [&saved](
                                const auto& checkpoint,
                                std::string&) {
        saved = checkpoint;
        return true;
    };
    fsi::MovingPorousFlowControlSession session(
        simulation, std::move(hooks));
    const auto ready = session.readyResponse();
    check(ready.kind == fsi::WorkerControlResponseKind::Ready
              && ready.requestId == 0
              && ready.acceptedStepCount == 0
              && ready.simulationTimeSeconds == 0.0,
          "moving porous control reports its initial safe point");
    checkValidResponse(ready,
                       "moving porous initial ready response is valid");

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError error;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 11, 101},
              response, &error)
              && !error
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 11
              && response.acceptedStepCount == 101
              && response.producedFrameCount == 101
              && frames.size() == 101
              && simulation.topologyRebaseCount() == 5
              && simulation.porousTopology().faceCoordinate == 0
              && simulation.porousTopology().periodicImage == 2,
          "moving porous control advances through both periodic wraps");
    checkValidResponse(
        response, "moving porous advanced response is valid");
    check(frames.front() != frames.back(),
          "moving porous control publishes owning distinct safe-point frames");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 12, 0},
              response, &error)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.requestId == 12
              && response.acceptedStepCount == 101
              && saved.has_value()
              && saved->acceptedStepCount == 101
              && saved->topologyRebaseCount == 5
              && saved->porousTopology.periodicImage == 2,
          "moving porous control delegates the complete wrapped checkpoint");
    checkValidResponse(
        response, "moving porous checkpointed response is valid");

    fsi::MovingPorousFlowCase replay;
    replay.restore(*saved);
    check(serializedFrame(simulation.advance())
              == serializedFrame(replay.advance()),
          "moving porous control checkpoint replays the exact next frame");
}

void testFailuresAndStopRemainAtSafePoints() {
    fsi::MovingPorousFlowCase simulation;
    std::uint64_t publicationCount = 0;
    fsi::MovingPorousFlowControlHooks hooks;
    hooks.publishFrame = [&publicationCount](
                             const viewer::DiagnosticFrame&,
                             std::string& error) {
        ++publicationCount;
        if (publicationCount == 2) {
            error = "trace unavailable";
            return false;
        }
        return true;
    };
    hooks.writeCheckpoint = [](
                                const auto&,
                                std::string& error) {
        error = "storage unavailable";
        return false;
    };
    fsi::MovingPorousFlowControlSession session(
        simulation, std::move(hooks));
    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 21, 3},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InternalFailure
              && response.acceptedStepCount == 2
              && simulation.acceptedStepCount() == 2
              && publicationCount == 2,
          "moving porous publication failure exposes the accepted safe point");
    checkValidResponse(
        response, "moving porous publication error response is valid");

    const auto beforeCheckpoint = simulation.checkpoint();
    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 22, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::CheckpointFailure
              && response.acceptedStepCount == 2
              && simulation.acceptedStepCount()
                  == beforeCheckpoint.acceptedStepCount
              && simulation.porousTopology()
                  == beforeCheckpoint.porousTopology,
          "moving porous checkpoint output failure cannot advance state");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 23, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.acceptedStepCount == 2
              && session.stopped(),
          "moving porous stop is terminal at the accepted safe point");
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 24, 1},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && simulation.acceptedStepCount() == 2,
          "stopped moving porous control rejects later advance");
}

void testMalformedCommandIsTransactional() {
    fsi::MovingPorousFlowCase simulation;
    fsi::MovingPorousFlowControlSession session(simulation);
    fsi::WorkerControlResponse response{
        fsi::WorkerControlResponseKind::Ready, 0, 9, 0.9, 0,
        fsi::WorkerControlFailureCode::None, {}};
    const auto preserved = response;
    fsi::WorkerControlProtocolError error;
    check(!session.execute(
              {fsi::WorkerControlCommandKind::Advance, 0, 1},
              response, &error)
              && error.code
                  == fsi::WorkerControlProtocolErrorCode::InvalidData
              && response == preserved
              && simulation.acceptedStepCount() == 0
              && !session.stopped(),
          "malformed moving porous command fails before mutation");
}

} // namespace

int main() {
    testAdvanceCheckpointAndReplay();
    testFailuresAndStopRemainAtSafePoints();
    testMalformedCommandIsTransactional();
    if (failures != 0) {
        std::fprintf(
            stderr,
            "%d moving porous-flow control check(s) failed\n",
            failures);
        return 1;
    }
    std::puts("all moving porous-flow control checks passed");
    return 0;
}
