#include "worker_control_stream.h"

#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

std::string encodedCommandBytes(const WorkerControlCommand& command) {
    std::ostringstream output(std::ios::out | std::ios::binary);
    WorkerControlStreamError error;
    check(writeWorkerControlCommand(output, command, &error) && !error,
          "stream helper writes a valid command");
    return output.str();
}

void testCommandSequenceAndCleanEof() {
    const std::vector<WorkerControlCommand> commands{
        {WorkerControlCommandKind::Advance, 1, 2},
        {WorkerControlCommandKind::Checkpoint, 2, 0},
        {WorkerControlCommandKind::Stop, 3, 0},
    };
    std::stringstream stream(
        std::ios::in | std::ios::out | std::ios::binary);
    WorkerControlStreamError error;
    for (const auto& command : commands) {
        check(writeWorkerControlCommand(stream, command, &error) && !error,
              "command stream writes one flushed self-framed message");
    }
    stream.seekg(0);
    for (const auto& expected : commands) {
        WorkerControlCommand decoded;
        check(readWorkerControlCommand(stream, decoded, &error)
                  == WorkerControlStreamResult::Message
                  && !error && decoded == expected,
              "command stream reads consecutive envelopes exactly");
    }
    WorkerControlCommand preserved{WorkerControlCommandKind::Advance, 9, 9};
    const auto before = preserved;
    check(readWorkerControlCommand(stream, preserved, &error)
              == WorkerControlStreamResult::EndOfStream
              && !error && preserved == before,
          "physical EOF between envelopes is clean and transactional");
}

void testResponseSequence() {
    const std::vector<WorkerControlResponse> responses{
        {WorkerControlResponseKind::Ready, 0, 4, 0.4, 0,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Advanced, 7, 6, 0.6, 2,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Stopped, 8, 6, 0.6, 0,
         WorkerControlFailureCode::None, {}},
    };
    std::stringstream stream(
        std::ios::in | std::ios::out | std::ios::binary);
    WorkerControlStreamError error;
    for (const auto& response : responses) {
        check(writeWorkerControlResponse(stream, response, &error) && !error,
              "response stream writes one flushed self-framed message");
    }
    stream.seekg(0);
    for (const auto& expected : responses) {
        WorkerControlResponse decoded;
        check(readWorkerControlResponse(stream, decoded, &error)
                  == WorkerControlStreamResult::Message
                  && !error && decoded == expected,
              "response stream reads consecutive envelopes exactly");
    }
}

void testReadFailuresAreBoundedAndTransactional() {
    const std::string valid = encodedCommandBytes(
        {WorkerControlCommandKind::Advance, 10, 3});
    const WorkerControlCommand preserved{
        WorkerControlCommandKind::Stop, 99, 0};
    WorkerControlCommand decoded = preserved;
    WorkerControlStreamError error;

    std::istringstream shortEnvelope(
        valid.substr(0, workerControlEnvelopeBytes - 1),
        std::ios::in | std::ios::binary);
    check(readWorkerControlCommand(shortEnvelope, decoded, &error)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::Truncated
              && decoded == preserved,
          "stream rejects an envelope cut short without output mutation");

    std::istringstream shortPayload(
        valid.substr(0, valid.size() - 1),
        std::ios::in | std::ios::binary);
    decoded = preserved;
    check(readWorkerControlCommand(shortPayload, decoded, &error)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::Truncated
              && decoded == preserved,
          "stream rejects a payload cut short without output mutation");

    std::string oversized = valid.substr(0, workerControlEnvelopeBytes);
    for (std::size_t byte = 0; byte < 4; ++byte) {
        oversized[workerControlPayloadSizeOffset + byte] =
            static_cast<char>(0xff);
    }
    std::istringstream oversizedStream(
        oversized, std::ios::in | std::ios::binary);
    decoded = preserved;
    check(readWorkerControlCommand(oversizedStream, decoded, &error)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::LimitExceeded
              && decoded == preserved,
          "declared oversized payload is rejected before allocation or read");

    std::string corrupt = valid;
    corrupt.back() ^= static_cast<char>(0x80);
    std::istringstream corruptStream(
        corrupt, std::ios::in | std::ios::binary);
    decoded = preserved;
    check(readWorkerControlCommand(corruptStream, decoded, &error)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::ProtocolFailure
              && error.protocolCode
                  == WorkerControlProtocolErrorCode::ChecksumMismatch
              && decoded == preserved,
          "stream forwards protocol corruption without output mutation");

    WorkerControlProtocolLimits smallLimit;
    smallLimit.maximumMessageBytes = valid.size() - 1;
    std::istringstream limitedStream(valid, std::ios::in | std::ios::binary);
    decoded = preserved;
    check(readWorkerControlCommand(
              limitedStream, decoded, &error, smallLimit)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::LimitExceeded
              && decoded == preserved,
          "stream enforces configured whole-message limits");
}

void testTypeAndIoFailures() {
    const std::string command = encodedCommandBytes(
        {WorkerControlCommandKind::Checkpoint, 20, 0});
    std::istringstream wrongType(command, std::ios::in | std::ios::binary);
    WorkerControlResponse response{
        WorkerControlResponseKind::Ready, 0, 1, 0.1, 0,
        WorkerControlFailureCode::None, {}};
    const auto before = response;
    WorkerControlStreamError error;
    check(readWorkerControlResponse(wrongType, response, &error)
              == WorkerControlStreamResult::Error
              && error.code == WorkerControlStreamErrorCode::ProtocolFailure
              && error.protocolCode
                  == WorkerControlProtocolErrorCode::InvalidMagic
              && response == before,
          "stream retains distinct command and response wire types");

    std::ostringstream failedOutput(std::ios::out | std::ios::binary);
    failedOutput.setstate(std::ios::badbit);
    check(!writeWorkerControlCommand(
              failedOutput,
              {WorkerControlCommandKind::Stop, 21, 0}, &error)
              && error.code == WorkerControlStreamErrorCode::IoFailure,
          "stream reports a flushed output failure");

    std::ostringstream protocolOutput(std::ios::out | std::ios::binary);
    check(!writeWorkerControlCommand(
              protocolOutput,
              {WorkerControlCommandKind::Advance, 0, 1}, &error)
              && error.code
                  == WorkerControlStreamErrorCode::ProtocolFailure
              && error.protocolCode
                  == WorkerControlProtocolErrorCode::InvalidData
              && protocolOutput.str().empty(),
          "invalid outgoing command writes no protocol bytes");
}

} // namespace

int main() {
    testCommandSequenceAndCleanEof();
    testResponseSequence();
    testReadFailuresAreBoundedAndTransactional();
    testTypeAndIoFailures();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d worker control stream check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all worker control stream checks passed");
    return 0;
}
