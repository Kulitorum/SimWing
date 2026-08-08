#include "worker_control_stream.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <istream>
#include <limits>
#include <new>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

void clearError(WorkerControlStreamError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

WorkerControlStreamResult failRead(
    WorkerControlStreamError* error,
    const WorkerControlStreamErrorCode code,
    std::string message,
    const WorkerControlProtocolErrorCode protocolCode =
        WorkerControlProtocolErrorCode::None) {
    if (error != nullptr) {
        error->code = code;
        error->protocolCode = protocolCode;
        error->message = std::move(message);
    }
    return WorkerControlStreamResult::Error;
}

bool failWrite(
    WorkerControlStreamError* error,
    const WorkerControlStreamErrorCode code,
    std::string message,
    const WorkerControlProtocolErrorCode protocolCode =
        WorkerControlProtocolErrorCode::None) {
    if (error != nullptr) {
        error->code = code;
        error->protocolCode = protocolCode;
        error->message = std::move(message);
    }
    return false;
}

std::uint32_t payloadSizeFromEnvelope(
    const std::span<const std::uint8_t> envelope) noexcept {
    std::uint32_t value = 0;
    for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
        value |= static_cast<std::uint32_t>(
                     envelope[workerControlPayloadSizeOffset + byte])
            << (8U * byte);
    }
    return value;
}

template<typename Message, typename Decode>
WorkerControlStreamResult readMessage(
    std::istream& input,
    Message& message,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits,
    Decode&& decode) {
    clearError(error);
    if (limits.maximumMessageBytes < workerControlEnvelopeBytes
        || limits.maximumMessageBytes
            > std::numeric_limits<std::size_t>::max()) {
        return failRead(
            error, WorkerControlStreamErrorCode::LimitExceeded,
            "configured control-stream byte limit is invalid");
    }

    std::array<std::uint8_t, workerControlEnvelopeBytes> envelope{};
    input.read(reinterpret_cast<char*>(envelope.data()),
               static_cast<std::streamsize>(envelope.size()));
    const std::streamsize envelopeRead = input.gcount();
    if (envelopeRead == 0 && input.eof() && !input.bad()) {
        clearError(error);
        return WorkerControlStreamResult::EndOfStream;
    }
    if (envelopeRead != static_cast<std::streamsize>(envelope.size())) {
        return failRead(
            error,
            input.bad() ? WorkerControlStreamErrorCode::IoFailure
                        : WorkerControlStreamErrorCode::Truncated,
            input.bad() ? "control stream failed while reading an envelope"
                        : "control stream ended inside an envelope");
    }

    const std::uint32_t payloadSize = payloadSizeFromEnvelope(envelope);
    const std::uint64_t totalSize =
        workerControlEnvelopeBytes + static_cast<std::uint64_t>(payloadSize);
    if (totalSize > limits.maximumMessageBytes
        || totalSize > std::numeric_limits<std::size_t>::max()) {
        return failRead(
            error, WorkerControlStreamErrorCode::LimitExceeded,
            "control stream message exceeds the byte limit");
    }

    try {
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(totalSize));
        std::ranges::copy(envelope, bytes.begin());
        if (payloadSize != 0) {
            input.read(
                reinterpret_cast<char*>(
                    bytes.data() + workerControlEnvelopeBytes),
                static_cast<std::streamsize>(payloadSize));
            if (input.gcount() != static_cast<std::streamsize>(payloadSize)) {
                return failRead(
                    error,
                    input.bad() ? WorkerControlStreamErrorCode::IoFailure
                                : WorkerControlStreamErrorCode::Truncated,
                    input.bad()
                        ? "control stream failed while reading a payload"
                        : "control stream ended inside a payload");
            }
        }
        WorkerControlProtocolError protocolError;
        Message candidate;
        if (!decode(bytes, candidate, &protocolError, limits)) {
            return failRead(
                error, WorkerControlStreamErrorCode::ProtocolFailure,
                protocolError.message, protocolError.code);
        }
        message = std::move(candidate);
        return WorkerControlStreamResult::Message;
    } catch (const std::bad_alloc&) {
        return failRead(
            error, WorkerControlStreamErrorCode::LimitExceeded,
            "unable to allocate the bounded control stream message");
    } catch (const std::length_error&) {
        return failRead(
            error, WorkerControlStreamErrorCode::LimitExceeded,
            "control stream message exceeds platform limits");
    }
}

template<typename Message, typename Encode>
bool writeMessage(
    std::ostream& output,
    const Message& message,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits,
    Encode&& encode) {
    clearError(error);
    std::vector<std::uint8_t> bytes;
    WorkerControlProtocolError protocolError;
    if (!encode(message, bytes, &protocolError, limits)) {
        return failWrite(
            error, WorkerControlStreamErrorCode::ProtocolFailure,
            protocolError.message, protocolError.code);
    }
    output.write(reinterpret_cast<const char*>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
    output.flush();
    if (!output) {
        return failWrite(
            error, WorkerControlStreamErrorCode::IoFailure,
            "control stream failed while writing a complete message");
    }
    return true;
}

} // namespace

WorkerControlStreamResult readWorkerControlCommand(
    std::istream& input,
    WorkerControlCommand& command,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits) {
    return readMessage(
        input, command, error, limits,
        deserializeWorkerControlCommand);
}

WorkerControlStreamResult readWorkerControlResponse(
    std::istream& input,
    WorkerControlResponse& response,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits) {
    return readMessage(
        input, response, error, limits,
        deserializeWorkerControlResponse);
}

bool writeWorkerControlCommand(
    std::ostream& output,
    const WorkerControlCommand& command,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits) {
    return writeMessage(
        output, command, error, limits,
        serializeWorkerControlCommand);
}

bool writeWorkerControlResponse(
    std::ostream& output,
    const WorkerControlResponse& response,
    WorkerControlStreamError* error,
    const WorkerControlProtocolLimits& limits) {
    return writeMessage(
        output, response, error, limits,
        serializeWorkerControlResponse);
}

} // namespace simwing::fsi
