#include "strong_piston_control.h"

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
          "strong control frame serializes");
    return bytes;
}

void testAcceptedSafePointsAndCheckpointReplay() {
    fsi::StrongCoupledPistonWorkerCase simulation;
    std::vector<std::vector<std::uint8_t>> publishedFrames;
    std::optional<fsi::StrongCoupledPistonCheckpoint> savedCheckpoint;
    fsi::StrongPistonControlHooks hooks;
    hooks.publishFrame = [&](const viewer::DiagnosticFrame& frame,
                             std::string&) {
        publishedFrames.push_back(serializedFrame(frame));
        return true;
    };
    hooks.writeCheckpoint = [&savedCheckpoint](
                                const fsi::StrongCoupledPistonCheckpoint& value,
                                std::string&) {
        savedCheckpoint = value;
        return true;
    };
    fsi::StrongPistonControlSession session(simulation, std::move(hooks));

    const auto ready = session.readyResponse();
    check(ready.kind == fsi::WorkerControlResponseKind::Ready
              && ready.acceptedStepCount == 0
              && ready.simulationTimeSeconds == 0.0,
          "strong control starts at the initial accepted safe point");
    checkValidResponse(ready, "strong control Ready is protocol-valid");

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 11, 2},
              response, &protocolError)
              && !protocolError
              && response.kind == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 11
              && response.acceptedStepCount == 2
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 2
              && publishedFrames.size() == 2,
          "strong control advances only complete accepted macro-steps");
    checkValidResponse(response, "strong control Advanced is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 12, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.acceptedStepCount == 2
              && savedCheckpoint.has_value()
              && savedCheckpoint->structure.acceptedStepCount == 2
              && savedCheckpoint->structure.simulationTimeSeconds
                  == simulation.simulationTimeSeconds(),
          "strong control checkpoints the complete accepted coupled epoch");
    checkValidResponse(response,
                       "strong control Checkpointed is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 13, 1},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Advanced
              && response.acceptedStepCount == 3
              && publishedFrames.size() == 3,
          "strong control continues after checkpointing");

    fsi::StrongCoupledPistonWorkerCase rebuilt;
    rebuilt.restore(*savedCheckpoint);
    check(serializedFrame(rebuilt.advance()) == publishedFrames.back(),
          "strong control checkpoint replays the exact next accepted frame");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 14, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.acceptedStepCount == 3
              && session.stopped(),
          "strong control stop is terminal at an accepted safe point");
    checkValidResponse(response, "strong control Stopped is protocol-valid");
}

void testOutputFailureKeepsCommittedState() {
    fsi::StrongCoupledPistonWorkerCase simulation;
    std::uint64_t publicationCount = 0;
    fsi::StrongPistonControlHooks hooks;
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
    fsi::StrongPistonControlSession session(simulation, std::move(hooks));
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
          "strong control reports publication failure at the committed epoch");
    checkValidResponse(response,
                       "strong publication failure is protocol-valid");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 22, 1},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Advanced
              && response.acceptedStepCount == 3
              && publicationCount == 3,
          "strong control resumes after an output-only failure");
}

void testMissingCheckpointSinkDoesNotAdvance() {
    fsi::StrongCoupledPistonWorkerCase simulation;
    fsi::StrongPistonControlSession session(simulation);
    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 31, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::CheckpointFailure
              && response.acceptedStepCount == 0
              && simulation.acceptedStepCount() == 0,
          "missing strong checkpoint sink cannot advance coupled state");
    checkValidResponse(response,
                       "missing strong checkpoint response is protocol-valid");
}

} // namespace

int main() {
    testAcceptedSafePointsAndCheckpointReplay();
    testOutputFailureKeepsCommittedState();
    testMissingCheckpointSinkDoesNotAdvance();
    if (failures != 0) {
        std::fprintf(stderr, "%d strong piston control test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("strong piston control tests passed\n");
    return 0;
}
