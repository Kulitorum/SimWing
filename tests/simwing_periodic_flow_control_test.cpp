#include "periodic_flow_control.h"

#include <cstdio>
#include <optional>
#include <stdexcept>
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

void checkValidResponse(
    const fsi::WorkerControlResponse& response,
    const char* message) {
    fsi::WorkerControlProtocolError error;
    check(fsi::validateWorkerControlResponse(response, &error) && !error,
          message);
}

std::vector<std::uint8_t> serializedFrame(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "control-session frame serializes");
    return bytes;
}

void testReadyAdvanceAndCheckpoint() {
    fsi::PeriodicFlowCase simulation;
    std::vector<std::vector<std::uint8_t>> publishedFrames;
    std::optional<fsi::PeriodicFlowCaseCheckpoint> savedCheckpoint;
    fsi::PeriodicFlowControlHooks hooks;
    hooks.publishFrame = [&](const viewer::DiagnosticFrame& frame,
                             std::string&) {
        publishedFrames.push_back(serializedFrame(frame));
        return true;
    };
    hooks.writeCheckpoint = [&savedCheckpoint](
                                const fsi::PeriodicFlowCaseCheckpoint& value,
                                std::string&) {
        savedCheckpoint = value;
        return true;
    };
    fsi::PeriodicFlowControlSession session(simulation, std::move(hooks));

    const auto ready = session.readyResponse();
    check(ready.kind == fsi::WorkerControlResponseKind::Ready
              && ready.requestId == 0
              && ready.acceptedStepCount == 0
              && ready.simulationTimeSeconds == 0.0,
          "ready response reports the initial absolute safe point");
    checkValidResponse(ready, "initial ready response is protocol-valid");

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError{
        fsi::WorkerControlProtocolErrorCode::InvalidData, "old"};
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 11, 3},
              response, &protocolError)
              && !protocolError
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 11
              && response.acceptedStepCount == 3
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 3
              && publishedFrames.size() == 3,
          "advance executes all requested accepted steps and correlates them");
    checkValidResponse(response, "advanced response is protocol-valid");
    check(publishedFrames[0] != publishedFrames[1]
              && publishedFrames[1] != publishedFrames[2],
          "advance publishes distinct owning frames at consecutive safe points");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 12, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.requestId == 12
              && response.acceptedStepCount == 3
              && savedCheckpoint.has_value()
              && savedCheckpoint->acceptedStepCount == 3
              && savedCheckpoint->simulationTimeSeconds
                  == simulation.simulationTimeSeconds(),
          "checkpoint delegates the complete current safe point");
    checkValidResponse(response, "checkpointed response is protocol-valid");

    fsi::PeriodicFlowCase rebuilt;
    rebuilt.restore(*savedCheckpoint);
    check(serializedFrame(simulation.advance())
              == serializedFrame(rebuilt.advance()),
          "delegated checkpoint resumes the exact next accepted frame");
}

void testPublicationFailureLeavesAcceptedSafePoint() {
    fsi::PeriodicFlowCase simulation;
    std::uint64_t callbackCount = 0;
    fsi::PeriodicFlowControlHooks hooks;
    hooks.publishFrame = [&](const viewer::DiagnosticFrame& frame,
                             std::string& error) {
        ++callbackCount;
        check(frame.step == callbackCount,
              "frame callback sees the newly accepted absolute step");
        if (callbackCount == 2) {
            error.assign(5000, 'x');
            return false;
        }
        return true;
    };
    fsi::PeriodicFlowControlSession session(simulation, std::move(hooks));
    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;

    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 21, 3},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.requestId == 21
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InternalFailure
              && response.acceptedStepCount == 2
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.errorMessage.size()
                  == fsi::WorkerControlProtocolLimits{}
                         .maximumErrorMessageBytes
              && callbackCount == 2
              && simulation.acceptedStepCount() == 2,
          "publication failure reports the last committed safe point boundedly");
    checkValidResponse(response,
                       "publication failure response is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 22, 1},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.acceptedStepCount == 3
              && response.producedFrameCount == 1
              && callbackCount == 3,
          "session continues from the accepted step after output failure");
}

void testCheckpointFailuresDoNotAdvance() {
    fsi::PeriodicFlowCase simulation;
    fsi::PeriodicFlowControlSession missingSink(simulation);
    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(missingSink.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 31, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::CheckpointFailure
              && response.acceptedStepCount == 0
              && simulation.acceptedStepCount() == 0,
          "missing checkpoint persistence reports an error without advancing");
    checkValidResponse(response,
                       "missing-checkpoint response is protocol-valid");

    fsi::PeriodicFlowControlHooks hooks;
    hooks.writeCheckpoint = [](
                                const fsi::PeriodicFlowCaseCheckpoint&,
                                std::string&) -> bool {
        throw std::runtime_error("storage unavailable");
    };
    fsi::PeriodicFlowControlSession throwingSink(
        simulation, std::move(hooks));
    check(throwingSink.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 32, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::CheckpointFailure
              && response.errorMessage.find("storage unavailable")
                  != std::string::npos
              && simulation.acceptedStepCount() == 0,
          "checkpoint callback exceptions become correlated bounded errors");
    checkValidResponse(response,
                       "thrown-checkpoint response is protocol-valid");
}

void testStopAndMalformedCommands() {
    fsi::PeriodicFlowCase simulation;
    std::uint64_t published = 0;
    fsi::PeriodicFlowControlHooks hooks;
    hooks.publishFrame = [&published](
                             const viewer::DiagnosticFrame&,
                             std::string&) {
        ++published;
        return true;
    };
    fsi::PeriodicFlowControlSession session(simulation, std::move(hooks));
    fsi::WorkerControlResponse response{
        fsi::WorkerControlResponseKind::Ready, 0, 9, 0.9, 0,
        fsi::WorkerControlFailureCode::None, {}};
    const auto preserved = response;
    fsi::WorkerControlProtocolError protocolError;

    check(!session.execute(
              {fsi::WorkerControlCommandKind::Advance, 0, 1},
              response, &protocolError)
              && protocolError.code
                  == fsi::WorkerControlProtocolErrorCode::InvalidData
              && response == preserved
              && simulation.acceptedStepCount() == 0
              && !session.stopped(),
          "malformed commands fail transactionally before worker mutation");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 41, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.requestId == 41
              && response.acceptedStepCount == 0
              && session.stopped(),
          "stop commits the terminal state at a safe point");
    checkValidResponse(response, "stopped response is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 42, 1},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && response.acceptedStepCount == 0
              && simulation.acceptedStepCount() == 0
              && published == 0,
          "a stopped session rejects later advance without mutation");
    checkValidResponse(response,
                       "post-stop rejection response is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 43, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && response.acceptedStepCount == 0
              && simulation.acceptedStepCount() == 0,
          "a stopped session rejects later checkpoint without mutation");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 44, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.requestId == 44
              && response.acceptedStepCount == 0,
          "repeated stop is idempotent and freshly correlated");
}

} // namespace

int main() {
    testReadyAdvanceAndCheckpoint();
    testPublicationFailureLeavesAcceptedSafePoint();
    testCheckpointFailuresDoNotAdvance();
    testStopAndMalformedCommands();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d periodic flow control check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all periodic flow control checks passed");
    return 0;
}
