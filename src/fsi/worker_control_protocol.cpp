#include "worker_control_protocol.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 4> commandMagic{'S', 'W', 'C', 'C'};
constexpr std::array<std::uint8_t, 4> responseMagic{'S', 'W', 'C', 'R'};
constexpr std::size_t envelopeBytes = 20;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(WorkerControlProtocolError* error,
          const WorkerControlProtocolErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(WorkerControlProtocolError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

std::uint64_t checksum(const std::span<const std::uint8_t> bytes) noexcept {
    std::uint64_t result = fnvOffsetBasis;
    for (const std::uint8_t value : bytes) {
        result ^= value;
        result *= fnvPrime;
    }
    return result;
}

class Writer final {
public:
    Writer(std::vector<std::uint8_t>& bytes, const std::size_t limit)
        : bytes_(bytes), limit_(limit) {}

    bool u8(const std::uint8_t value) {
        return raw(&value, sizeof(value));
    }

    bool u16(const std::uint16_t value) {
        return integer(value);
    }

    bool u32(const std::uint32_t value) {
        return integer(value);
    }

    bool u64(const std::uint64_t value) {
        return integer(value);
    }

    bool finiteDouble(const double value) {
        return std::isfinite(value)
            && u64(std::bit_cast<std::uint64_t>(value));
    }

    bool string(const std::string& value) {
        return value.size() <= std::numeric_limits<std::uint32_t>::max()
            && u32(static_cast<std::uint32_t>(value.size()))
            && raw(value.data(), value.size());
    }

    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }

private:
    template<typename Unsigned>
    bool integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        std::array<std::uint8_t, sizeof(Unsigned)> bytes{};
        for (std::size_t index = 0; index < bytes.size(); ++index) {
            bytes[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return raw(bytes.data(), bytes.size());
    }

    bool raw(const void* data, const std::size_t size) {
        if (size > limit_ || bytes_.size() > limit_ - size) {
            exceeded_ = true;
            return false;
        }
        const auto* begin = static_cast<const std::uint8_t*>(data);
        bytes_.insert(bytes_.end(), begin, begin + size);
        return true;
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t limit_ = 0;
    bool exceeded_ = false;
};

class Reader final {
public:
    explicit Reader(const std::span<const std::uint8_t> bytes)
        : bytes_(bytes) {}

    bool u8(std::uint8_t& value) { return integer(value); }
    bool u16(std::uint16_t& value) { return integer(value); }
    bool u32(std::uint32_t& value) { return integer(value); }
    bool u64(std::uint64_t& value) { return integer(value); }

    bool finiteDouble(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) {
            return false;
        }
        value = std::bit_cast<double>(bits);
        return std::isfinite(value);
    }

    bool string(std::string& value, const std::uint32_t maximumBytes) {
        std::uint32_t size = 0;
        if (!u32(size)) {
            return false;
        }
        if (size > maximumBytes) {
            limitExceeded_ = true;
            return false;
        }
        if (size > bytes_.size() - position_) {
            truncated_ = true;
            return false;
        }
        value.assign(
            reinterpret_cast<const char*>(bytes_.data() + position_), size);
        position_ += size;
        return true;
    }

    [[nodiscard]] bool atEnd() const noexcept {
        return position_ == bytes_.size();
    }
    [[nodiscard]] bool truncated() const noexcept { return truncated_; }
    [[nodiscard]] bool limitExceeded() const noexcept {
        return limitExceeded_;
    }

private:
    template<typename Unsigned>
    bool integer(Unsigned& value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        if (sizeof(Unsigned) > bytes_.size() - position_) {
            truncated_ = true;
            return false;
        }
        value = 0;
        for (std::size_t index = 0; index < sizeof(Unsigned); ++index) {
            value |= static_cast<Unsigned>(bytes_[position_ + index])
                << (8U * index);
        }
        position_ += sizeof(Unsigned);
        return true;
    }

    std::span<const std::uint8_t> bytes_;
    std::size_t position_ = 0;
    bool truncated_ = false;
    bool limitExceeded_ = false;
};

WorkerControlProtocolErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return WorkerControlProtocolErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return WorkerControlProtocolErrorCode::Truncated;
    }
    return WorkerControlProtocolErrorCode::InvalidData;
}

template<typename Kind>
bool validKind(const Kind kind, const Kind first, const Kind last) {
    const auto encoded = static_cast<std::uint8_t>(kind);
    return encoded >= static_cast<std::uint8_t>(first)
        && encoded <= static_cast<std::uint8_t>(last);
}

template<typename Message, typename WritePayload>
bool serializeMessage(
    const Message& message,
    const std::span<const std::uint8_t> magic,
    const std::uint8_t kind,
    std::vector<std::uint8_t>& bytes,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits,
    WritePayload&& writePayload) {
    bytes.clear();
    if (limits.maximumMessageBytes < envelopeBytes
        || limits.maximumMessageBytes
            > std::numeric_limits<std::size_t>::max()) {
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "configured control-message byte limit is invalid");
    }
    try {
        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload,
            static_cast<std::size_t>(limits.maximumMessageBytes)
                - envelopeBytes);
        if (!writePayload(message, payloadWriter)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? WorkerControlProtocolErrorCode::LimitExceeded
                    : WorkerControlProtocolErrorCode::InvalidData,
                "worker control payload cannot be serialized");
        }
        bytes.reserve(envelopeBytes + payload.size());
        bytes.insert(bytes.end(), magic.begin(), magic.end());
        Writer envelope(bytes,
                        static_cast<std::size_t>(limits.maximumMessageBytes));
        if (!envelope.u16(workerControlProtocolVersion)
            || !envelope.u8(kind)
            || !envelope.u8(0)
            || payload.size() > std::numeric_limits<std::uint32_t>::max()
            || !envelope.u32(static_cast<std::uint32_t>(payload.size()))
            || !envelope.u64(checksum(payload))) {
            bytes.clear();
            return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                        "worker control envelope exceeds the byte limit");
        }
        bytes.insert(bytes.end(), payload.begin(), payload.end());
        return true;
    } catch (const std::bad_alloc&) {
        bytes.clear();
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "unable to allocate worker control message");
    } catch (const std::length_error&) {
        bytes.clear();
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "worker control message exceeds platform limits");
    }
}

template<typename Kind, typename ReadPayload>
bool deserializeMessage(
    const std::span<const std::uint8_t> bytes,
    const std::span<const std::uint8_t> magic,
    Kind& kind,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits,
    ReadPayload&& readPayload) {
    if (limits.maximumMessageBytes < envelopeBytes
        || limits.maximumMessageBytes
            > std::numeric_limits<std::size_t>::max()
        || bytes.size() > limits.maximumMessageBytes) {
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "worker control message exceeds the byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error, WorkerControlProtocolErrorCode::Truncated,
                    "worker control envelope is truncated");
    }
    if (!std::ranges::equal(magic, bytes.first(magic.size()))) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidMagic,
                    "worker control message magic is invalid");
    }
    Reader envelope(bytes.subspan(magic.size()));
    std::uint16_t version = 0;
    std::uint8_t encodedKind = 0;
    std::uint8_t reserved = 0;
    std::uint32_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(version)
        || !envelope.u8(encodedKind)
        || !envelope.u8(reserved)
        || !envelope.u32(payloadSize)
        || !envelope.u64(expectedChecksum)) {
        return fail(error, WorkerControlProtocolErrorCode::Truncated,
                    "worker control envelope is truncated");
    }
    if (version != workerControlProtocolVersion) {
        return fail(error, WorkerControlProtocolErrorCode::UnsupportedVersion,
                    "worker control protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control reserved bits are nonzero");
    }
    const std::size_t availablePayload = bytes.size() - envelopeBytes;
    if (payloadSize > availablePayload) {
        return fail(error, WorkerControlProtocolErrorCode::Truncated,
                    "worker control payload is truncated");
    }
    if (payloadSize < availablePayload) {
        return fail(error, WorkerControlProtocolErrorCode::TrailingData,
                    "worker control message has trailing data");
    }
    const auto payload = bytes.subspan(envelopeBytes, payloadSize);
    if (checksum(payload) != expectedChecksum) {
        return fail(error, WorkerControlProtocolErrorCode::ChecksumMismatch,
                    "worker control checksum does not match");
    }
    kind = static_cast<Kind>(encodedKind);
    Reader reader(payload);
    if (!readPayload(reader)) {
        return fail(error, readerErrorCode(reader),
                    "worker control payload is invalid");
    }
    if (!reader.atEnd()) {
        return fail(error, WorkerControlProtocolErrorCode::TrailingData,
                    "worker control payload has trailing data");
    }
    return true;
}

} // namespace

bool validateWorkerControlCommand(
    const WorkerControlCommand& command,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    if (!validKind(
            command.kind,
            WorkerControlCommandKind::Advance,
            WorkerControlCommandKind::Stop)
        || command.requestId == 0) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control command identity is invalid");
    }
    if (command.kind == WorkerControlCommandKind::Advance) {
        if (command.advanceStepCount == 0) {
            return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                        "advance command requires at least one step");
        }
        if (command.advanceStepCount > limits.maximumAdvanceSteps) {
            return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                        "advance command exceeds the step limit");
        }
    } else if (command.advanceStepCount != 0) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "non-advance command carries an advance count");
    }
    return true;
}

bool validateWorkerControlResponse(
    const WorkerControlResponse& response,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    if (!validKind(
            response.kind,
            WorkerControlResponseKind::Ready,
            WorkerControlResponseKind::Error)
        || !std::isfinite(response.simulationTimeSeconds)
        || response.simulationTimeSeconds < 0.0) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control response state is invalid");
    }
    if ((response.kind == WorkerControlResponseKind::Ready)
            != (response.requestId == 0)) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control response request identity is invalid");
    }
    if (response.kind == WorkerControlResponseKind::Advanced) {
        if (response.producedFrameCount == 0
            || response.failureCode != WorkerControlFailureCode::None
            || !response.errorMessage.empty()) {
            return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                        "advanced response payload is invalid");
        }
    } else if (response.kind == WorkerControlResponseKind::Error) {
        if (response.producedFrameCount != 0
            || response.failureCode == WorkerControlFailureCode::None
            || response.errorMessage.empty()) {
            return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                        "error response payload is invalid");
        }
        if (response.errorMessage.size() > limits.maximumErrorMessageBytes) {
            return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                        "error response text exceeds the string limit");
        }
    } else if (response.producedFrameCount != 0
               || response.failureCode != WorkerControlFailureCode::None
               || !response.errorMessage.empty()) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control response carries fields for another kind");
    }
    if (static_cast<std::uint16_t>(response.failureCode)
        > static_cast<std::uint16_t>(
            WorkerControlFailureCode::InternalFailure)) {
        return fail(error, WorkerControlProtocolErrorCode::InvalidData,
                    "worker control failure code is invalid");
    }
    return true;
}

bool serializeWorkerControlCommand(
    const WorkerControlCommand& command,
    std::vector<std::uint8_t>& bytes,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    bytes.clear();
    if (!validateWorkerControlCommand(command, error, limits)) {
        return false;
    }
    return serializeMessage(
        command, commandMagic, static_cast<std::uint8_t>(command.kind),
        bytes, error, limits,
        [](const WorkerControlCommand& value, Writer& writer) {
            return writer.u64(value.requestId)
                && (value.kind != WorkerControlCommandKind::Advance
                    || writer.u64(value.advanceStepCount));
        });
}

bool deserializeWorkerControlCommand(
    const std::span<const std::uint8_t> bytes,
    WorkerControlCommand& command,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    try {
        WorkerControlCommand candidate;
        if (!deserializeMessage(
                bytes, commandMagic, candidate.kind, error, limits,
                [&](Reader& reader) {
                    return reader.u64(candidate.requestId)
                        && (candidate.kind
                                != WorkerControlCommandKind::Advance
                            || reader.u64(candidate.advanceStepCount));
                })
            || !validateWorkerControlCommand(candidate, error, limits)) {
            return false;
        }
        command = candidate;
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "unable to allocate worker control command");
    }
}

bool serializeWorkerControlResponse(
    const WorkerControlResponse& response,
    std::vector<std::uint8_t>& bytes,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    bytes.clear();
    if (!validateWorkerControlResponse(response, error, limits)) {
        return false;
    }
    return serializeMessage(
        response, responseMagic, static_cast<std::uint8_t>(response.kind),
        bytes, error, limits,
        [](const WorkerControlResponse& value, Writer& writer) {
            if (!writer.u64(value.requestId)
                || !writer.u64(value.acceptedStepCount)
                || !writer.finiteDouble(value.simulationTimeSeconds)) {
                return false;
            }
            if (value.kind == WorkerControlResponseKind::Advanced) {
                return writer.u64(value.producedFrameCount);
            }
            if (value.kind == WorkerControlResponseKind::Error) {
                return writer.u16(
                           static_cast<std::uint16_t>(value.failureCode))
                    && writer.string(value.errorMessage);
            }
            return true;
        });
}

bool deserializeWorkerControlResponse(
    const std::span<const std::uint8_t> bytes,
    WorkerControlResponse& response,
    WorkerControlProtocolError* error,
    const WorkerControlProtocolLimits& limits) {
    clearError(error);
    try {
        WorkerControlResponse candidate;
        if (!deserializeMessage(
                bytes, responseMagic, candidate.kind, error, limits,
                [&](Reader& reader) {
                    if (!reader.u64(candidate.requestId)
                        || !reader.u64(candidate.acceptedStepCount)
                        || !reader.finiteDouble(
                            candidate.simulationTimeSeconds)) {
                        return false;
                    }
                    if (candidate.kind
                        == WorkerControlResponseKind::Advanced) {
                        return reader.u64(candidate.producedFrameCount);
                    }
                    if (candidate.kind == WorkerControlResponseKind::Error) {
                        std::uint16_t failure = 0;
                        if (!reader.u16(failure)
                            || failure > static_cast<std::uint16_t>(
                                WorkerControlFailureCode::InternalFailure)
                            || !reader.string(
                                candidate.errorMessage,
                                limits.maximumErrorMessageBytes)) {
                            return false;
                        }
                        candidate.failureCode =
                            static_cast<WorkerControlFailureCode>(failure);
                    }
                    return true;
                })
            || !validateWorkerControlResponse(candidate, error, limits)) {
            return false;
        }
        response = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "unable to allocate worker control response");
    } catch (const std::length_error&) {
        return fail(error, WorkerControlProtocolErrorCode::LimitExceeded,
                    "worker control response exceeds platform limits");
    }
}

} // namespace simwing::fsi
