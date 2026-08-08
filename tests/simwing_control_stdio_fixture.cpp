#include "moving_porous_flow_case.h"
#include "moving_porous_flow_checkpoint_persistence.h"
#include "open_piston_case.h"
#include "open_piston_checkpoint_persistence.h"
#include "periodic_flow_case.h"
#include "porous_sheet_case.h"
#include "porous_sheet_checkpoint_persistence.h"
#include "viewer_protocol.h"
#include "worker_control_stream.h"

#include <charconv>
#include <cstdio>
#include <exception>
#include <fstream>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace simwing;

bool readFile(
    const std::string& path,
    const std::uint64_t maximumBytes,
    std::vector<std::uint8_t>& bytes) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
        std::fprintf(stderr, "cannot open fixture file: %s\n", path.c_str());
        return false;
    }
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::uint64_t>(size) > maximumBytes
        || static_cast<std::uint64_t>(size)
            > std::numeric_limits<std::size_t>::max()) {
        std::fprintf(stderr, "fixture file has an invalid size: %s\n",
                     path.c_str());
        return false;
    }
    bytes.resize(static_cast<std::size_t>(size));
    input.seekg(0);
    if (!bytes.empty()) {
        input.read(reinterpret_cast<char*>(bytes.data()), size);
    }
    if (!input || input.peek() != std::char_traits<char>::eof()) {
        std::fprintf(stderr, "cannot read complete fixture file: %s\n",
                     path.c_str());
        return false;
    }
    return true;
}

bool readCheckpoint(
    const std::string& path,
    fsi::PeriodicFlowCaseCheckpoint& checkpoint) {
    const fsi::PeriodicFlowCaseCheckpointLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readFile(path, limits.maximumBytes, bytes)) {
        return false;
    }
    fsi::PeriodicFlowCaseCheckpointError error;
    if (!fsi::deserializePeriodicFlowCaseCheckpoint(
            bytes, checkpoint, &error, limits)) {
        std::fprintf(stderr, "control checkpoint is invalid: %s\n",
                     error.message.c_str());
        return false;
    }
    return true;
}

bool readMovingPorousFlowCheckpoint(
    const std::string& path,
    fsi::MovingPorousFlowCaseCheckpoint& checkpoint) {
    const fsi::MovingPorousFlowCaseCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readFile(path, limits.maximumEncodedBytes, bytes)) {
        return false;
    }
    fsi::MovingPorousFlowCaseCheckpointPersistenceError error;
    if (!fsi::deserializeMovingPorousFlowCaseCheckpoint(
            bytes, checkpoint, &error, limits)) {
        std::fprintf(
            stderr,
            "moving porous-flow control checkpoint is invalid: %s\n",
            error.message.c_str());
        return false;
    }
    return true;
}

bool readOpenPistonCheckpoint(
    const std::string& path,
    fsi::OpenPistonCaseCheckpoint& checkpoint) {
    const fsi::OpenPistonCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readFile(path, limits.maximumEncodedBytes, bytes)) {
        return false;
    }
    fsi::OpenPistonCheckpointPersistenceError error;
    if (!fsi::deserializeOpenPistonCheckpoint(
            bytes, checkpoint, &error, limits)) {
        std::fprintf(stderr, "open-piston control checkpoint is invalid: %s\n",
                     error.message.c_str());
        return false;
    }
    return true;
}

bool readPorousSheetCheckpoint(
    const std::string& path,
    fsi::CoupledPorousSheetCheckpoint& checkpoint) {
    const fsi::CoupledPorousSheetCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readFile(path, limits.maximumEncodedBytes, bytes)) {
        return false;
    }
    fsi::CoupledPorousSheetCase owner;
    fsi::CoupledPorousSheetCheckpointPersistenceError error;
    if (!fsi::deserializeCoupledPorousSheetCheckpoint(
            bytes, owner, checkpoint, &error, limits)) {
        std::fprintf(stderr,
                     "porous-sheet control checkpoint is invalid: %s\n",
                     error.message.c_str());
        return false;
    }
    return true;
}

bool verifyResponseFile(
    const std::string& path,
    const std::vector<fsi::WorkerControlResponse>& expected) {
    std::ifstream responses(path, std::ios::binary);
    if (!responses) {
        std::fprintf(stderr, "cannot open response fixture: %s\n",
                     path.c_str());
        return false;
    }
    fsi::WorkerControlStreamError streamError;
    for (const auto& value : expected) {
        fsi::WorkerControlResponse decoded;
        if (fsi::readWorkerControlResponse(
                responses, decoded, &streamError)
                != fsi::WorkerControlStreamResult::Message
            || decoded != value) {
            std::fprintf(stderr,
                         "control response sequence is not exact: %s\n",
                         streamError.message.c_str());
            return false;
        }
    }
    fsi::WorkerControlResponse trailing;
    if (fsi::readWorkerControlResponse(
            responses, trailing, &streamError)
        != fsi::WorkerControlStreamResult::EndOfStream) {
        std::fprintf(stderr,
                     "control stdout has a trailing or malformed message\n");
        return false;
    }
    return true;
}

std::vector<std::uint8_t> serializedFrame(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    if (!viewer::serializeFrame(frame, bytes, &error)) {
        std::fprintf(stderr, "cannot serialize fixture frame: %s\n",
                     error.message.c_str());
        return {};
    }
    return bytes;
}

int writeCommands(const std::string& path,
                  const std::uint64_t advanceStepCount,
                  const std::uint64_t failureAdvanceStepCount) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr, "cannot create command fixture: %s\n",
                     path.c_str());
        return 1;
    }
    std::vector<fsi::WorkerControlCommand> commands{
        {fsi::WorkerControlCommandKind::Advance, 101, advanceStepCount},
        {fsi::WorkerControlCommandKind::Checkpoint, 102, 0},
        {fsi::WorkerControlCommandKind::Advance, 103, 1},
    };
    if (failureAdvanceStepCount != 0) {
        commands.push_back({
            fsi::WorkerControlCommandKind::Advance, 104,
            failureAdvanceStepCount});
        commands.push_back({
            fsi::WorkerControlCommandKind::Stop, 105, 0});
    } else {
        commands.push_back({
            fsi::WorkerControlCommandKind::Stop, 104, 0});
    }
    fsi::WorkerControlStreamError error;
    for (const auto& command : commands) {
        if (!fsi::writeWorkerControlCommand(output, command, &error)) {
            std::fprintf(stderr, "cannot write command fixture: %s\n",
                         error.message.c_str());
            return 1;
        }
    }
    return 0;
}

int writeResumeCommands(const std::string& path) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr, "cannot create resume command fixture: %s\n",
                     path.c_str());
        return 1;
    }
    fsi::WorkerControlStreamError error;
    if (!fsi::writeWorkerControlCommand(
            output,
            {fsi::WorkerControlCommandKind::Advance, 201, 1}, &error)
        || !fsi::writeWorkerControlCommand(
            output,
            {fsi::WorkerControlCommandKind::Stop, 202, 0}, &error)) {
        std::fprintf(stderr, "cannot write resume command fixture: %s\n",
                     error.message.c_str());
        return 1;
    }
    return 0;
}

int writePorousCollisionCommands(const std::string& path,
                                 const bool resumed) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr,
                     "cannot create porous collision command fixture: %s\n",
                     path.c_str());
        return 1;
    }
    const std::vector<fsi::WorkerControlCommand> commands = resumed
        ? std::vector<fsi::WorkerControlCommand>{
              {fsi::WorkerControlCommandKind::Advance, 401, 1},
              {fsi::WorkerControlCommandKind::Stop, 402, 0},
          }
        : std::vector<fsi::WorkerControlCommand>{
              {fsi::WorkerControlCommandKind::Advance, 301, 5000},
              {fsi::WorkerControlCommandKind::Checkpoint, 302, 0},
              {fsi::WorkerControlCommandKind::Stop, 303, 0},
          };
    fsi::WorkerControlStreamError error;
    for (const auto& command : commands) {
        if (!fsi::writeWorkerControlCommand(output, command, &error)) {
            std::fprintf(stderr,
                         "cannot write porous collision command fixture: %s\n",
                         error.message.c_str());
            return 1;
        }
    }
    return 0;
}

int verifyResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::PeriodicFlowCase definition;
    const double stepSeconds =
        definition.stepSettings().flow.timeStepSeconds;
    const double timeAtStep2 = stepSeconds + stepSeconds;
    const double timeAtStep3 = timeAtStep2 + stepSeconds;
    const std::vector<fsi::WorkerControlResponse> expected{
        {fsi::WorkerControlResponseKind::Ready, 0, 0, 0.0, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 101, 2, timeAtStep2, 2,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Checkpointed, 102, 2, timeAtStep2, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 103, 3, timeAtStep3, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 104, 3, timeAtStep3, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expected)) {
        return 1;
    }

    fsi::PeriodicFlowCaseCheckpoint checkpoint;
    if (!readCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    if (checkpoint.acceptedStepCount != 2
        || checkpoint.simulationTimeSeconds != timeAtStep2) {
        std::fprintf(stderr,
                     "control checkpoint safe point is not step two\n");
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr, "cannot open control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::periodicFlowCaseChecksum
        || header.solverCommit != fsi::periodicFlowCaseSolverId) {
        std::fprintf(stderr, "control trace header is invalid\n");
        return 1;
    }
    std::vector<viewer::DiagnosticFrame> frames;
    for (;;) {
        viewer::DiagnosticFrame frame;
        const viewer::TraceReadStatus status = reader.readNext(frame);
        if (status == viewer::TraceReadStatus::Frame) {
            frames.push_back(std::move(frame));
        } else if (status == viewer::TraceReadStatus::End) {
            break;
        } else {
            std::fprintf(stderr,
                         "control trace is incomplete or invalid: %s\n",
                         reader.error().message.c_str());
            return 1;
        }
    }
    if (frames.size() != 3
        || frames[0].step != 1
        || frames[1].step != 2
        || frames[2].step != 3
        || frames[1].simulationTimeSeconds != timeAtStep2
        || frames[2].simulationTimeSeconds != timeAtStep3) {
        std::fprintf(stderr,
                     "control trace does not contain three accepted steps\n");
        return 1;
    }

    fsi::PeriodicFlowCase resumed;
    resumed.restore(checkpoint);
    if (serializedFrame(resumed.advance())
        != serializedFrame(frames[2])) {
        std::fprintf(stderr,
                     "control checkpoint does not replay trace step three\n");
        return 1;
    }
    return 0;
}

int verifyResume(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::PeriodicFlowCaseCheckpoint checkpoint;
    if (!readCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    fsi::PeriodicFlowCase simulation;
    simulation.restore(checkpoint);
    const viewer::DiagnosticFrame expectedFrame = simulation.advance();
    const std::vector<fsi::WorkerControlResponse> expectedResponses{
        {fsi::WorkerControlResponseKind::Ready, 0,
         checkpoint.acceptedStepCount, checkpoint.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 201,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 202,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expectedResponses)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr, "cannot open resumed control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::periodicFlowCaseChecksum
        || header.solverCommit != fsi::periodicFlowCaseSolverId
        || reader.readNext(frame) != viewer::TraceReadStatus::Frame
        || serializedFrame(frame) != serializedFrame(expectedFrame)
        || reader.readNext(frame) != viewer::TraceReadStatus::End) {
        std::fprintf(stderr,
                     "resumed control trace is not one exact continued frame\n");
        return 1;
    }
    return 0;
}

int verifyMovingPorousFlowResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    constexpr std::uint64_t checkpointStep = 101;
    constexpr std::uint64_t finalStep = checkpointStep + 1;
    fsi::MovingPorousFlowCase oracle;
    viewer::DiagnosticFrame checkpointFrame;
    for (std::uint64_t step = 0; step < checkpointStep; ++step) {
        checkpointFrame = oracle.advance();
    }
    const auto expectedCheckpoint = oracle.checkpoint();
    const auto finalFrame = oracle.advance();
    const std::vector<fsi::WorkerControlResponse> expectedResponses{
        {fsi::WorkerControlResponseKind::Ready, 0, 0, 0.0, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 101,
         checkpointStep, checkpointFrame.simulationTimeSeconds,
         checkpointStep, fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Checkpointed, 102,
         checkpointStep, checkpointFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 103,
         finalStep, finalFrame.simulationTimeSeconds, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 104,
         finalStep, finalFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expectedResponses)) {
        return 1;
    }

    fsi::MovingPorousFlowCaseCheckpoint checkpoint;
    if (!readMovingPorousFlowCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    if (checkpoint.acceptedStepCount != checkpointStep
        || checkpoint.simulationTimeSeconds
            != expectedCheckpoint.simulationTimeSeconds
        || checkpoint.sheetPositionMeters
            != expectedCheckpoint.sheetPositionMeters
        || checkpoint.topologyRebaseCount != 5
        || checkpoint.porousTopology.faceCoordinate != 0
        || checkpoint.porousTopology.periodicImage != 2) {
        std::fprintf(
            stderr,
            "moving porous-flow control checkpoint is not the second-wrap safe point\n");
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open moving porous-flow control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::movingPorousFlowCaseChecksum
        || header.solverCommit != fsi::movingPorousFlowCaseSolverId) {
        std::fprintf(
            stderr, "moving porous-flow control trace header is invalid\n");
        return 1;
    }
    fsi::MovingPorousFlowCase traceOracle;
    for (std::uint64_t step = 0; step < finalStep; ++step) {
        viewer::DiagnosticFrame actual;
        const auto expected = traceOracle.advance();
        if (reader.readNext(actual) != viewer::TraceReadStatus::Frame
            || serializedFrame(actual) != serializedFrame(expected)) {
            std::fprintf(
                stderr,
                "moving porous-flow control trace differs at step %llu\n",
                static_cast<unsigned long long>(step + 1));
            return 1;
        }
    }
    viewer::DiagnosticFrame trailing;
    if (reader.readNext(trailing) != viewer::TraceReadStatus::End) {
        std::fprintf(
            stderr,
            "moving porous-flow control trace is incomplete or has trailing frames\n");
        return 1;
    }

    fsi::MovingPorousFlowCase resumed;
    resumed.restore(checkpoint);
    if (serializedFrame(resumed.advance()) != serializedFrame(finalFrame)) {
        std::fprintf(
            stderr,
            "moving porous-flow checkpoint does not replay trace step 102\n");
        return 1;
    }
    return 0;
}

int verifyMovingPorousFlowResume(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::MovingPorousFlowCaseCheckpoint checkpoint;
    if (!readMovingPorousFlowCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    fsi::MovingPorousFlowCase simulation;
    simulation.restore(checkpoint);
    const auto expectedFrame = simulation.advance();
    const std::vector<fsi::WorkerControlResponse> expectedResponses{
        {fsi::WorkerControlResponseKind::Ready, 0,
         checkpoint.acceptedStepCount, checkpoint.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 201,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 202,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expectedResponses)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(
            stderr,
            "cannot open resumed moving porous-flow control trace: %s\n",
            tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::movingPorousFlowCaseChecksum
        || header.solverCommit != fsi::movingPorousFlowCaseSolverId
        || reader.readNext(frame) != viewer::TraceReadStatus::Frame
        || serializedFrame(frame) != serializedFrame(expectedFrame)
        || reader.readNext(frame) != viewer::TraceReadStatus::End) {
        std::fprintf(
            stderr,
            "resumed moving porous-flow trace is not one exact continued frame\n");
        return 1;
    }
    return 0;
}

int verifyOpenPistonResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::OpenPistonCase definition;
    const double stepSeconds = definition.stepSettings().timeStepSeconds;
    const double timeAtStep2 = stepSeconds + stepSeconds;
    const double timeAtStep3 = timeAtStep2 + stepSeconds;
    const std::vector<fsi::WorkerControlResponse> expected{
        {fsi::WorkerControlResponseKind::Ready, 0, 0, 0.0, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 101, 2, timeAtStep2, 2,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Checkpointed, 102, 2, timeAtStep2, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 103, 3, timeAtStep3, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 104, 3, timeAtStep3, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expected)) {
        return 1;
    }

    fsi::OpenPistonCaseCheckpoint checkpoint;
    if (!readOpenPistonCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    if (checkpoint.acceptedStepCount != 2
        || checkpoint.surfaceOffsetMeters <= 0.0) {
        std::fprintf(stderr,
                     "open-piston control checkpoint is not step two\n");
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr, "cannot open open-piston control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::openPistonCaseChecksum
        || header.solverCommit != fsi::openPistonCaseSolverId) {
        std::fprintf(stderr, "open-piston control trace header is invalid\n");
        return 1;
    }
    std::vector<viewer::DiagnosticFrame> frames;
    for (;;) {
        viewer::DiagnosticFrame frame;
        const viewer::TraceReadStatus status = reader.readNext(frame);
        if (status == viewer::TraceReadStatus::Frame) {
            frames.push_back(std::move(frame));
        } else if (status == viewer::TraceReadStatus::End) {
            break;
        } else {
            std::fprintf(stderr,
                         "open-piston control trace is incomplete: %s\n",
                         reader.error().message.c_str());
            return 1;
        }
    }
    if (frames.size() != 3
        || frames[0].step != 1
        || frames[1].step != 2
        || frames[2].step != 3
        || frames[1].simulationTimeSeconds != timeAtStep2
        || frames[2].simulationTimeSeconds != timeAtStep3) {
        std::fprintf(stderr,
                     "open-piston control trace lacks three accepted steps\n");
        return 1;
    }

    fsi::OpenPistonCase resumed;
    resumed.restore(checkpoint);
    if (serializedFrame(resumed.advance())
        != serializedFrame(frames[2])) {
        std::fprintf(stderr,
                     "open-piston checkpoint does not replay trace step three\n");
        return 1;
    }
    return 0;
}

int verifyOpenPistonResume(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::OpenPistonCaseCheckpoint checkpoint;
    if (!readOpenPistonCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    fsi::OpenPistonCase simulation;
    simulation.restore(checkpoint);
    const double checkpointTime = simulation.simulationTimeSeconds();
    const viewer::DiagnosticFrame expectedFrame = simulation.advance();
    const std::vector<fsi::WorkerControlResponse> expectedResponses{
        {fsi::WorkerControlResponseKind::Ready, 0,
         checkpoint.acceptedStepCount, checkpointTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 201,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 202,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expectedResponses)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open resumed open-piston control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::openPistonCaseChecksum
        || header.solverCommit != fsi::openPistonCaseSolverId
        || reader.readNext(frame) != viewer::TraceReadStatus::Frame
        || serializedFrame(frame) != serializedFrame(expectedFrame)
        || reader.readNext(frame) != viewer::TraceReadStatus::End) {
        std::fprintf(stderr,
                     "resumed open-piston trace is not one continued frame\n");
        return 1;
    }
    return 0;
}

int verifyPorousSheetResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::CoupledPorousSheetCase definition;
    const double stepSeconds = definition.stepSettings().timeStepSeconds;
    constexpr std::uint64_t checkpointStep = 330;
    constexpr std::uint64_t finalStep = checkpointStep + 1;
    double timeAtCheckpoint = 0.0;
    for (std::uint64_t step = 0; step < checkpointStep; ++step) {
        timeAtCheckpoint += stepSeconds;
    }
    const double timeAtContinuation = timeAtCheckpoint + stepSeconds;

    fsi::CoupledPorousSheetCheckpoint checkpoint;
    if (!readPorousSheetCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    if (checkpoint.acceptedStepCount != checkpointStep
        || checkpoint.simulationTimeSeconds != timeAtCheckpoint
        || checkpoint.topologyRebaseCount != 1
        || checkpoint.porousTopology.faceCoordinate != 4
        || checkpoint.porousTopology.axis
            != fsi::fluid::GridFaceAxis::X
        || checkpoint.porousTopology.periodicImage != 0) {
        std::fprintf(stderr,
                     "porous-sheet control checkpoint is not the rebased safe point\n");
        return 1;
    }

    fsi::CoupledPorousSheetCase collisionReplay;
    collisionReplay.restore(checkpoint);
    viewer::DiagnosticFrame lastAcceptedCollisionFrame;
    std::string collisionError;
    for (std::uint64_t attempt = 0;
         attempt < 5000 && collisionError.empty(); ++attempt) {
        try {
            lastAcceptedCollisionFrame = collisionReplay.advance();
        } catch (const std::exception& exception) {
            collisionError = exception.what();
        }
    }
    const std::uint64_t collisionSafeStep =
        collisionReplay.acceptedStepCount();
    const double collisionSafeTime =
        collisionReplay.simulationTimeSeconds();
    if (collisionSafeStep <= checkpointStep + 1
        || collisionReplay.topologyRebaseCount()
            != fsi::coupledPorousSheetMaximumOrdinaryRebaseCount
        || collisionReplay.porousTopology()
            != fsi::coupledPorousSheetTerminalSafeTopology
        || collisionError
            != "coupled porous sheet reached the pump-surface topology") {
        std::fprintf(stderr,
                     "porous-sheet collision oracle did not reach its safe point\n");
        return 1;
    }

    const std::vector<fsi::WorkerControlResponse> expected{
        {fsi::WorkerControlResponseKind::Ready, 0, 0, 0.0, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 101, checkpointStep,
         timeAtCheckpoint, checkpointStep,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Checkpointed, 102, checkpointStep,
         timeAtCheckpoint, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 103,
         checkpointStep + 1, timeAtContinuation, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Error, 104,
         collisionSafeStep, collisionSafeTime, 0,
         fsi::WorkerControlFailureCode::NumericalFailure,
         "worker advance failed: " + collisionError},
        {fsi::WorkerControlResponseKind::Stopped, 105,
         collisionSafeStep, collisionSafeTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expected)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open porous-sheet control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::coupledPorousSheetCaseChecksum
        || header.solverCommit != fsi::coupledPorousSheetCaseSolverId) {
        std::fprintf(stderr, "porous-sheet control trace header is invalid\n");
        return 1;
    }
    std::vector<viewer::DiagnosticFrame> frames;
    for (;;) {
        viewer::DiagnosticFrame frame;
        const viewer::TraceReadStatus status = reader.readNext(frame);
        if (status == viewer::TraceReadStatus::Frame) {
            frames.push_back(std::move(frame));
        } else if (status == viewer::TraceReadStatus::End) {
            break;
        } else {
            std::fprintf(stderr,
                         "porous-sheet control trace is incomplete: %s\n",
                         reader.error().message.c_str());
            return 1;
        }
    }
    if (frames.size() != collisionSafeStep
        || frames[0].step != 1
        || frames[checkpointStep - 1].step != checkpointStep
        || frames.back().step != collisionSafeStep
        || frames[checkpointStep - 1].simulationTimeSeconds
            != timeAtCheckpoint
        || frames.back().simulationTimeSeconds != collisionSafeTime
        || serializedFrame(frames.back())
            != serializedFrame(lastAcceptedCollisionFrame)) {
        std::fprintf(stderr,
                     "porous-sheet control trace crossed its collision safe point\n");
        return 1;
    }

    fsi::CoupledPorousSheetCase resumed;
    resumed.restore(checkpoint);
    if (serializedFrame(resumed.advance())
        != serializedFrame(frames[checkpointStep])) {
        std::fprintf(stderr,
                     "porous-sheet checkpoint does not replay the next trace step\n");
        return 1;
    }
    return 0;
}

int verifyPorousSheetResume(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::CoupledPorousSheetCheckpoint checkpoint;
    if (!readPorousSheetCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    fsi::CoupledPorousSheetCase simulation;
    simulation.restore(checkpoint);
    const double checkpointTime = simulation.simulationTimeSeconds();
    const viewer::DiagnosticFrame expectedFrame = simulation.advance();
    const std::vector<fsi::WorkerControlResponse> expectedResponses{
        {fsi::WorkerControlResponseKind::Ready, 0,
         checkpoint.acceptedStepCount, checkpointTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Advanced, 201,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 1,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 202,
         expectedFrame.step, expectedFrame.simulationTimeSeconds, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expectedResponses)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open resumed porous-sheet control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::coupledPorousSheetCaseChecksum
        || header.solverCommit != fsi::coupledPorousSheetCaseSolverId
        || reader.readNext(frame) != viewer::TraceReadStatus::Frame
        || serializedFrame(frame) != serializedFrame(expectedFrame)
        || reader.readNext(frame) != viewer::TraceReadStatus::End) {
        std::fprintf(stderr,
                     "resumed porous-sheet trace is not one continued frame\n");
        return 1;
    }
    return 0;
}

int verifyPorousCollisionResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::CoupledPorousSheetCase oracle;
    viewer::DiagnosticFrame lastAcceptedFrame;
    std::string collisionError;
    for (std::uint64_t attempt = 0;
         attempt < 5000 && collisionError.empty(); ++attempt) {
        try {
            lastAcceptedFrame = oracle.advance();
        } catch (const std::exception& exception) {
            collisionError = exception.what();
        }
    }
    const std::uint64_t safeStep = oracle.acceptedStepCount();
    const double safeTime = oracle.simulationTimeSeconds();
    if (collisionError
            != "coupled porous sheet reached the pump-surface topology"
        || safeStep == 0 || lastAcceptedFrame.step != safeStep) {
        std::fprintf(stderr,
                     "porous collision fixture did not reach its terminal safe point\n");
        return 1;
    }

    const std::vector<fsi::WorkerControlResponse> expected{
        {fsi::WorkerControlResponseKind::Ready, 0, 0, 0.0, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Error, 301,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::NumericalFailure,
         "worker advance failed: " + collisionError},
        {fsi::WorkerControlResponseKind::Checkpointed, 302,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Stopped, 303,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expected)) {
        return 1;
    }

    fsi::CoupledPorousSheetCheckpoint checkpoint;
    if (!readPorousSheetCheckpoint(checkpointPath, checkpoint)
        || checkpoint.acceptedStepCount != safeStep
        || checkpoint.simulationTimeSeconds != safeTime
        || checkpoint.topologyRebaseCount
            != fsi::coupledPorousSheetMaximumOrdinaryRebaseCount
        || checkpoint.porousTopology
            != fsi::coupledPorousSheetTerminalSafeTopology) {
        std::fprintf(stderr,
                     "porous collision checkpoint is not the terminal safe point\n");
        return 1;
    }
    fsi::CoupledPorousSheetCase restored;
    restored.restore(checkpoint);
    bool repeatedCollision = false;
    try {
        static_cast<void>(restored.advance());
    } catch (const std::exception& exception) {
        repeatedCollision = std::string(exception.what()) == collisionError;
    }
    if (!repeatedCollision
        || restored.acceptedStepCount() != safeStep
        || restored.simulationTimeSeconds() != safeTime) {
        std::fprintf(stderr,
                     "porous collision checkpoint did not repeat transactionally\n");
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open porous collision control trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    std::uint64_t frameCount = 0;
    std::uint64_t lastFrameStep = 0;
    double lastFrameTime = 0.0;
    std::vector<std::uint8_t> serializedLastFrame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::coupledPorousSheetCaseChecksum
        || header.solverCommit != fsi::coupledPorousSheetCaseSolverId) {
        std::fprintf(stderr, "porous collision trace header is invalid\n");
        return 1;
    }
    for (;;) {
        const viewer::TraceReadStatus status = reader.readNext(frame);
        if (status == viewer::TraceReadStatus::Frame) {
            ++frameCount;
            lastFrameStep = frame.step;
            lastFrameTime = frame.simulationTimeSeconds;
            serializedLastFrame = serializedFrame(frame);
        } else if (status == viewer::TraceReadStatus::End) {
            break;
        } else {
            std::fprintf(stderr,
                         "porous collision trace is incomplete: %s\n",
                         reader.error().message.c_str());
            return 1;
        }
    }
    if (frameCount != safeStep || lastFrameStep != safeStep
        || lastFrameTime != safeTime
        || serializedLastFrame != serializedFrame(lastAcceptedFrame)) {
        std::fprintf(stderr,
                     "porous collision trace published beyond its safe point\n");
        return 1;
    }
    return 0;
}

int verifyPorousCollisionResume(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    fsi::CoupledPorousSheetCheckpoint checkpoint;
    if (!readPorousSheetCheckpoint(checkpointPath, checkpoint)) {
        return 1;
    }
    fsi::CoupledPorousSheetCase restored;
    restored.restore(checkpoint);
    const std::uint64_t safeStep = restored.acceptedStepCount();
    const double safeTime = restored.simulationTimeSeconds();
    constexpr std::string_view collisionError =
        "coupled porous sheet reached the pump-surface topology";
    bool repeatedCollision = false;
    try {
        static_cast<void>(restored.advance());
    } catch (const std::exception& exception) {
        repeatedCollision = exception.what() == collisionError;
    }
    if (!repeatedCollision
        || restored.acceptedStepCount() != safeStep
        || restored.simulationTimeSeconds() != safeTime) {
        std::fprintf(stderr,
                     "restored porous collision oracle did not remain terminal\n");
        return 1;
    }

    const std::vector<fsi::WorkerControlResponse> expected{
        {fsi::WorkerControlResponseKind::Ready, 0,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
        {fsi::WorkerControlResponseKind::Error, 401,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::NumericalFailure,
         "worker advance failed: " + std::string(collisionError)},
        {fsi::WorkerControlResponseKind::Stopped, 402,
         safeStep, safeTime, 0,
         fsi::WorkerControlFailureCode::None, {}},
    };
    if (!verifyResponseFile(responsePath, expected)) {
        return 1;
    }

    std::ifstream traceInput(tracePath, std::ios::binary);
    if (!traceInput) {
        std::fprintf(stderr,
                     "cannot open resumed porous collision trace: %s\n",
                     tracePath.c_str());
        return 1;
    }
    viewer::TraceReader reader(traceInput, viewer::TraceReadMode::Follow);
    viewer::TraceHeader header;
    viewer::DiagnosticFrame frame;
    if (!reader.readHeader(header)
        || header.sceneChecksum != fsi::coupledPorousSheetCaseChecksum
        || header.solverCommit != fsi::coupledPorousSheetCaseSolverId
        || reader.readNext(frame) != viewer::TraceReadStatus::End) {
        std::fprintf(stderr,
                     "resumed porous collision trace contains a rejected frame\n");
        return 1;
    }
    return 0;
}

} // namespace

int main(int argc, char* argv[]) {
    if ((argc >= 3 && argc <= 5)
        && std::string_view(argv[1]) == "write") {
        std::uint64_t advanceStepCount = 2;
        std::uint64_t failureAdvanceStepCount = 0;
        const auto parsePositive = [](const std::string_view text,
                                      std::uint64_t& value) {
            const auto parsed = std::from_chars(
                text.data(), text.data() + text.size(), value);
            return parsed.ec == std::errc{}
                && parsed.ptr == text.data() + text.size()
                && value != 0;
        };
        if (argc >= 4) {
            const std::string_view text(argv[3]);
            if (!parsePositive(text, advanceStepCount)) {
                std::fprintf(stderr,
                             "write advance count must be a positive integer\n");
                return 2;
            }
        }
        if (argc == 5) {
            const std::string_view text(argv[4]);
            if (!parsePositive(text, failureAdvanceStepCount)) {
                std::fprintf(stderr,
                             "write failure advance count must be a positive integer\n");
                return 2;
            }
        }
        return writeCommands(
            argv[2], advanceStepCount, failureAdvanceStepCount);
    }
    if (argc == 3 && std::string_view(argv[1]) == "write-resume") {
        return writeResumeCommands(argv[2]);
    }
    if (argc == 3
        && std::string_view(argv[1]) == "write-porous-collision") {
        return writePorousCollisionCommands(argv[2], false);
    }
    if (argc == 3
        && std::string_view(argv[1])
            == "write-porous-collision-resume") {
        return writePorousCollisionCommands(argv[2], true);
    }
    if (argc == 5 && std::string_view(argv[1]) == "verify") {
        return verifyResponses(argv[2], argv[3], argv[4]);
    }
    if (argc == 5 && std::string_view(argv[1]) == "verify-resume") {
        return verifyResume(argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1]) == "verify-moving-porous") {
        return verifyMovingPorousFlowResponses(
            argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1])
            == "verify-moving-porous-resume") {
        return verifyMovingPorousFlowResume(
            argv[2], argv[3], argv[4]);
    }
    if (argc == 5 && std::string_view(argv[1]) == "verify-open") {
        return verifyOpenPistonResponses(argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1]) == "verify-open-resume") {
        return verifyOpenPistonResume(argv[2], argv[3], argv[4]);
    }
    if (argc == 5 && std::string_view(argv[1]) == "verify-porous") {
        return verifyPorousSheetResponses(argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1]) == "verify-porous-resume") {
        return verifyPorousSheetResume(argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1]) == "verify-porous-collision") {
        return verifyPorousCollisionResponses(argv[2], argv[3], argv[4]);
    }
    if (argc == 5
        && std::string_view(argv[1])
            == "verify-porous-collision-resume") {
        return verifyPorousCollisionResume(argv[2], argv[3], argv[4]);
    }
    std::fprintf(stderr,
                 "usage: simwing-control-stdio-fixture "
                 "write COMMANDS [ADVANCE_STEPS [FAILURE_STEPS]] | "
                 "write-resume COMMANDS | "
                 "write-porous-collision COMMANDS | "
                 "write-porous-collision-resume COMMANDS | "
                 "verify RESPONSES CHECKPOINT TRACE | "
                 "verify-resume RESPONSES CHECKPOINT TRACE | "
                 "verify-moving-porous RESPONSES CHECKPOINT TRACE | "
                 "verify-moving-porous-resume RESPONSES CHECKPOINT TRACE | "
                 "verify-open RESPONSES CHECKPOINT TRACE | "
                 "verify-open-resume RESPONSES CHECKPOINT TRACE | "
                 "verify-porous RESPONSES CHECKPOINT TRACE | "
                 "verify-porous-resume RESPONSES CHECKPOINT TRACE | "
                 "verify-porous-collision RESPONSES CHECKPOINT TRACE | "
                 "verify-porous-collision-resume RESPONSES CHECKPOINT TRACE\n");
    return 2;
}
