#include "worker_control_protocol.h"

#include <cstdio>
#include <limits>
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

std::vector<std::uint8_t> encoded(const WorkerControlCommand& command) {
    std::vector<std::uint8_t> bytes;
    WorkerControlProtocolError error;
    check(serializeWorkerControlCommand(command, bytes, &error),
          "valid worker control command serializes");
    check(!error, "command serialization clears its error");
    return bytes;
}

std::vector<std::uint8_t> encoded(const WorkerControlResponse& response) {
    std::vector<std::uint8_t> bytes;
    WorkerControlProtocolError error;
    check(serializeWorkerControlResponse(response, bytes, &error),
          "valid worker control response serializes");
    check(!error, "response serialization clears its error");
    return bytes;
}

void refreshChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 20;
    constexpr std::size_t checksumOffset = 12;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t checksum = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        checksum ^= bytes[index];
        checksum *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(checksum); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            checksum >> (8U * byte));
    }
}

void testCommandRoundTrips() {
    const std::vector<WorkerControlCommand> commands{
        {WorkerControlCommandKind::Advance, 1, 7},
        {WorkerControlCommandKind::Checkpoint, 2, 0},
        {WorkerControlCommandKind::Stop, 3, 0},
    };
    for (const auto& command : commands) {
        const auto first = encoded(command);
        const auto second = encoded(command);
        WorkerControlCommand decoded;
        WorkerControlProtocolError error;
        check(first == second && !first.empty(),
              "command encoding is byte deterministic");
        check(deserializeWorkerControlCommand(
                  first, decoded, &error)
                  && !error && decoded == command
                  && encoded(decoded) == first,
              "command round trip preserves every field and byte");
    }
}

void testResponseRoundTrips() {
    const std::vector<WorkerControlResponse> responses{
        {WorkerControlResponseKind::Ready, 0, 4, 0.25, 0,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Advanced, 11, 7, 0.5, 3,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Checkpointed, 12, 7, 0.5, 0,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Stopped, 13, 7, 0.5, 0,
         WorkerControlFailureCode::None, {}},
        {WorkerControlResponseKind::Error, 14, 7, 0.5, 0,
         WorkerControlFailureCode::CheckpointFailure,
         "checkpoint target is unavailable"},
    };
    for (const auto& response : responses) {
        const auto first = encoded(response);
        const auto second = encoded(response);
        WorkerControlResponse decoded;
        WorkerControlProtocolError error;
        check(first == second && !first.empty(),
              "response encoding is byte deterministic");
        check(deserializeWorkerControlResponse(
                  first, decoded, &error)
                  && !error && decoded == response
                  && encoded(decoded) == first,
              "response round trip preserves every field and byte");
    }
}

void testValidationAndLimits() {
    WorkerControlProtocolError error;
    std::vector<std::uint8_t> bytes{1, 2, 3};
    check(!serializeWorkerControlCommand(
              {WorkerControlCommandKind::Advance, 0, 1}, bytes, &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData
              && bytes.empty(),
          "command rejects a zero request ID transactionally");
    bytes = {1, 2, 3};
    check(!serializeWorkerControlCommand(
              {WorkerControlCommandKind::Advance, 1, 0}, bytes, &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData
              && bytes.empty(),
          "advance command rejects a zero step count");
    check(!validateWorkerControlCommand(
              {WorkerControlCommandKind::Stop, 1, 1}, &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData,
          "non-advance command rejects an advance count");
    WorkerControlProtocolLimits stepLimit;
    stepLimit.maximumAdvanceSteps = 2;
    check(!validateWorkerControlCommand(
              {WorkerControlCommandKind::Advance, 1, 3},
              &error, stepLimit)
              && error.code == WorkerControlProtocolErrorCode::LimitExceeded,
          "advance command honors the configured step limit");

    check(!validateWorkerControlResponse(
              {WorkerControlResponseKind::Ready, 1, 0, 0.0, 0,
               WorkerControlFailureCode::None, {}},
              &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData,
          "ready response reserves request ID zero");
    check(!validateWorkerControlResponse(
              {WorkerControlResponseKind::Advanced, 1, 1, 0.1, 0,
               WorkerControlFailureCode::None, {}},
              &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData,
          "advanced response requires produced frames");
    check(!validateWorkerControlResponse(
              {WorkerControlResponseKind::Error, 1, 1, 0.1, 0,
               WorkerControlFailureCode::None, "failed"},
              &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData,
          "error response requires a failure code");
    check(!validateWorkerControlResponse(
              {WorkerControlResponseKind::Stopped, 1, 1,
               std::numeric_limits<double>::infinity(), 0,
               WorkerControlFailureCode::None, {}},
              &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidData,
          "response rejects non-finite worker time");
    WorkerControlProtocolLimits stringLimit;
    stringLimit.maximumErrorMessageBytes = 3;
    bytes = {1, 2, 3};
    check(!serializeWorkerControlResponse(
              {WorkerControlResponseKind::Error, 1, 1, 0.1, 0,
               WorkerControlFailureCode::InternalFailure, "long"},
              bytes, &error, stringLimit)
              && error.code == WorkerControlProtocolErrorCode::LimitExceeded
              && bytes.empty(),
          "error response honors its configured text limit");
    WorkerControlProtocolLimits byteLimit;
    byteLimit.maximumMessageBytes = 19;
    bytes = {1, 2, 3};
    check(!serializeWorkerControlCommand(
              {WorkerControlCommandKind::Stop, 1, 0},
              bytes, &error, byteLimit)
              && error.code == WorkerControlProtocolErrorCode::LimitExceeded
              && bytes.empty(),
          "control encoding honors its configured byte limit");
}

void testCorruptionAndTransactionalDecode() {
    const WorkerControlCommand preservedCommand{
        WorkerControlCommandKind::Advance, 90, 5};
    WorkerControlCommand command = preservedCommand;
    WorkerControlProtocolError error;
    const auto valid = encoded(WorkerControlCommand{
        WorkerControlCommandKind::Checkpoint, 7, 0});
    const auto expectRejected = [&](
        const std::vector<std::uint8_t>& candidate,
        const WorkerControlProtocolErrorCode expected,
        const WorkerControlProtocolLimits& limits,
        const char* message) {
        error = {};
        check(!deserializeWorkerControlCommand(
                  candidate, command, &error, limits)
                  && error.code == expected
                  && command == preservedCommand,
              message);
    };

    auto corrupt = valid;
    corrupt[0] ^= 0xffU;
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::InvalidMagic, {},
        "control command rejects bad magic without changing output");
    corrupt = valid;
    ++corrupt[4];
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::UnsupportedVersion, {},
        "control command rejects unsupported versions");
    corrupt = valid;
    corrupt[7] = 1;
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::InvalidData, {},
        "control command rejects reserved bits");
    corrupt = valid;
    corrupt[6] = 0;
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::InvalidData, {},
        "control command rejects unknown command kinds");
    corrupt = valid;
    corrupt.pop_back();
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::Truncated, {},
        "control command rejects truncation transactionally");
    corrupt = valid;
    corrupt.push_back(0);
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::TrailingData, {},
        "control command rejects trailing data transactionally");
    corrupt = valid;
    corrupt.back() ^= 0x80U;
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::ChecksumMismatch, {},
        "control command detects payload corruption");
    WorkerControlProtocolLimits byteLimit;
    byteLimit.maximumMessageBytes = valid.size() - 1;
    expectRejected(
        valid, WorkerControlProtocolErrorCode::LimitExceeded, byteLimit,
        "control command decode enforces its byte limit");

    WorkerControlResponse wrongType;
    check(!deserializeWorkerControlResponse(valid, wrongType, &error)
              && error.code == WorkerControlProtocolErrorCode::InvalidMagic,
          "command and response envelopes cannot be confused");

    const WorkerControlResponse longError{
        WorkerControlResponseKind::Error, 8, 3, 0.3, 0,
        WorkerControlFailureCode::InternalFailure, "four"};
    const auto longErrorBytes = encoded(longError);
    WorkerControlProtocolLimits textLimit;
    textLimit.maximumErrorMessageBytes = 3;
    WorkerControlResponse preservedResponse{
        WorkerControlResponseKind::Stopped, 2, 2, 0.2, 0,
        WorkerControlFailureCode::None, {}};
    const auto before = preservedResponse;
    check(!deserializeWorkerControlResponse(
              longErrorBytes, preservedResponse, &error, textLimit)
              && error.code == WorkerControlProtocolErrorCode::LimitExceeded
              && preservedResponse == before,
          "response text limit rejection leaves output unchanged");

    corrupt = valid;
    corrupt[20] = 0;
    refreshChecksum(corrupt);
    expectRejected(
        corrupt, WorkerControlProtocolErrorCode::InvalidData, {},
        "control command validates decoded request IDs after checksum");
}

} // namespace

int main() {
    testCommandRoundTrips();
    testResponseRoundTrips();
    testValidationAndLimits();
    testCorruptionAndTransactionalDecode();
    if (failures != 0) {
        std::fprintf(stderr, "%d worker control protocol check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all worker control protocol checks passed");
    return 0;
}
