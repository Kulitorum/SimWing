#include "periodic_flow_case.h"
#include "viewer_protocol.h"
#include "worker_control_stream.h"

#include <cstdio>
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

int writeCommands(const std::string& path) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr, "cannot create command fixture: %s\n",
                     path.c_str());
        return 1;
    }
    const std::vector<fsi::WorkerControlCommand> commands{
        {fsi::WorkerControlCommandKind::Advance, 101, 2},
        {fsi::WorkerControlCommandKind::Checkpoint, 102, 0},
        {fsi::WorkerControlCommandKind::Advance, 103, 1},
        {fsi::WorkerControlCommandKind::Stop, 104, 0},
    };
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

int verifyResponses(
    const std::string& responsePath,
    const std::string& checkpointPath,
    const std::string& tracePath) {
    std::ifstream responses(responsePath, std::ios::binary);
    if (!responses) {
        std::fprintf(stderr, "cannot open response fixture: %s\n",
                     responsePath.c_str());
        return 1;
    }
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
            return 1;
        }
    }
    fsi::WorkerControlResponse trailing;
    if (fsi::readWorkerControlResponse(
            responses, trailing, &streamError)
        != fsi::WorkerControlStreamResult::EndOfStream) {
        std::fprintf(stderr,
                     "control stdout has a trailing or malformed message\n");
        return 1;
    }

    const fsi::PeriodicFlowCaseCheckpointLimits checkpointLimits;
    std::vector<std::uint8_t> checkpointBytes;
    if (!readFile(
            checkpointPath, checkpointLimits.maximumBytes,
            checkpointBytes)) {
        return 1;
    }
    fsi::PeriodicFlowCaseCheckpoint checkpoint;
    fsi::PeriodicFlowCaseCheckpointError checkpointError;
    if (!fsi::deserializePeriodicFlowCaseCheckpoint(
            checkpointBytes, checkpoint, &checkpointError,
            checkpointLimits)
        || checkpoint.acceptedStepCount != 2
        || checkpoint.simulationTimeSeconds != timeAtStep2) {
        std::fprintf(stderr, "control checkpoint is invalid: %s\n",
                     checkpointError.message.c_str());
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

} // namespace

int main(int argc, char* argv[]) {
    if (argc == 3 && std::string_view(argv[1]) == "write") {
        return writeCommands(argv[2]);
    }
    if (argc == 5 && std::string_view(argv[1]) == "verify") {
        return verifyResponses(argv[2], argv[3], argv[4]);
    }
    std::fprintf(stderr,
                 "usage: simwing-control-stdio-fixture "
                 "write COMMANDS | verify RESPONSES CHECKPOINT TRACE\n");
    return 2;
}
