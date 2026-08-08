#include "porous_sheet_control.h"

#include <cstdio>
#include <exception>
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
              {fsi::WorkerControlCommandKind::Advance, 401, 330},
              response, &protocolError)
              && !protocolError
              && response.kind
                  == fsi::WorkerControlResponseKind::Advanced
              && response.requestId == 401
              && response.acceptedStepCount == 330
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 330
              && publishedFrameCount == 330
              && lastPublishedFrame.has_value()
              && lastPublishedFrame->step == 330
              && simulation.topologyRebaseCount() == 1
              && simulation.porousFaceCoordinate() == 4,
          "porous-sheet control publishes through the first topology rebase");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 402, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.requestId == 402
              && response.acceptedStepCount == 330
              && savedCheckpoint.has_value()
              && savedCheckpoint->acceptedStepCount == 330
              && savedCheckpoint->topologyRebaseCount == 1
              && savedCheckpoint->porousFaceCoordinate == 4
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
              && response.acceptedStepCount == 331
              && response.producedFrameCount == 1
              && publishedFrameCount == 331
              && serializedFrame(*lastPublishedFrame)
                  == serializedFrame(expectedNext),
          "delegated porous-sheet checkpoint reproduces the controlled next frame");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 404, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.requestId == 404
              && response.acceptedStepCount == 331
              && session.stopped(),
          "porous-sheet control stops at the final accepted safe point");
    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 405, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.failureCode
                  == fsi::WorkerControlFailureCode::InvalidCommand
              && response.acceptedStepCount == 331
              && simulation.acceptedStepCount() == 331,
          "stopped porous-sheet control rejects later checkpoint mutation");
}

void testNumericalFailureRemainsAtSafePoint() {
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

    fsi::WorkerControlResponse response;
    fsi::WorkerControlProtocolError protocolError;
    check(session.execute(
              {fsi::WorkerControlCommandKind::Advance, 501, 1000},
              response, &protocolError)
              && !protocolError
              && response.kind == fsi::WorkerControlResponseKind::Error
              && response.requestId == 501
              && response.failureCode
                  == fsi::WorkerControlFailureCode::NumericalFailure
              && response.errorMessage
                  == "worker advance failed: coupled porous sheet reached the pump-surface topology"
              && response.acceptedStepCount
                  == simulation.acceptedStepCount()
              && response.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && response.producedFrameCount == 0
              && publishedFrameCount == response.acceptedStepCount
              && lastPublishedFrame.has_value()
              && lastPublishedFrame->step == response.acceptedStepCount
              && response.acceptedStepCount > 330
              && simulation.topologyRebaseCount() == 1
              && simulation.porousFaceCoordinate() == 4,
          "porous-sheet control exposes the last accepted state after pump collision");

    const std::uint64_t safeStep = simulation.acceptedStepCount();
    const double safeTime = simulation.simulationTimeSeconds();
    check(session.execute(
              {fsi::WorkerControlCommandKind::Checkpoint, 502, 0},
              response, &protocolError)
              && response.kind
                  == fsi::WorkerControlResponseKind::Checkpointed
              && response.acceptedStepCount == safeStep
              && response.simulationTimeSeconds == safeTime
              && savedCheckpoint.has_value()
              && savedCheckpoint->acceptedStepCount == safeStep
              && savedCheckpoint->simulationTimeSeconds == safeTime
              && savedCheckpoint->topologyRebaseCount == 1
              && savedCheckpoint->porousFaceCoordinate == 4
              && publishedFrameCount == safeStep,
          "porous-sheet control checkpoints only the collision safe point");
    if (!savedCheckpoint) {
        return;
    }

    fsi::CoupledPorousSheetCase restored;
    restored.restore(*savedCheckpoint);
    bool repeatedFailure = false;
    try {
        static_cast<void>(restored.advance());
    } catch (const std::exception& exception) {
        repeatedFailure = std::string(exception.what())
            == "coupled porous sheet reached the pump-surface topology";
    }
    check(repeatedFailure
              && restored.acceptedStepCount() == safeStep
              && restored.simulationTimeSeconds() == safeTime
              && restored.topologyRebaseCount() == 1
              && restored.porousFaceCoordinate() == 4,
          "porous-sheet control collision checkpoint repeats without mutation");

    check(session.execute(
              {fsi::WorkerControlCommandKind::Stop, 503, 0},
              response, &protocolError)
              && response.kind == fsi::WorkerControlResponseKind::Stopped
              && response.acceptedStepCount == safeStep
              && response.simulationTimeSeconds == safeTime,
          "porous-sheet control can stop after a numerical rejection");
}

} // namespace

int main() {
    testCheckpointAndContinuation();
    testNumericalFailureRemainsAtSafePoint();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d porous-sheet control check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all porous-sheet control checks passed");
    return 0;
}
