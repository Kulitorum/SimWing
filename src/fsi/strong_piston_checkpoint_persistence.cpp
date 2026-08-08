#include "strong_piston_checkpoint_persistence.h"

#include <algorithm>
#include <array>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 8> checkpointMagic{
    'S', 'W', 'S', 'P', 'C', 'K', 'P', '1'};
constexpr std::size_t checkpointEnvelopeBytes = 28;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(
    StrongPistonCheckpointPersistenceError* error,
    const StrongPistonCheckpointPersistenceErrorCode code,
    std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(StrongPistonCheckpointPersistenceError* error) {
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

class BufferWriter final {
public:
    BufferWriter(std::vector<std::uint8_t>& bytes, const std::size_t limit)
        : bytes_(bytes), limit_(limit) {}

    bool u16(const std::uint16_t value) { return integer(value); }
    bool u32(const std::uint32_t value) { return integer(value); }
    bool u64(const std::uint64_t value) { return integer(value); }

    bool raw(const std::span<const std::uint8_t> value) {
        if (value.size() > limit_ || bytes_.size() > limit_ - value.size()) {
            return false;
        }
        bytes_.insert(bytes_.end(), value.begin(), value.end());
        return true;
    }

private:
    template<typename Unsigned>
    bool integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        std::array<std::uint8_t, sizeof(Unsigned)> encoded{};
        for (std::size_t index = 0; index < encoded.size(); ++index) {
            encoded[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return raw(encoded);
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t limit_ = 0;
};

class BufferReader final {
public:
    explicit BufferReader(const std::span<const std::uint8_t> bytes)
        : bytes_(bytes) {}

    bool u16(std::uint16_t& value) { return integer(value); }
    bool u32(std::uint32_t& value) { return integer(value); }
    bool u64(std::uint64_t& value) { return integer(value); }

    bool raw(const std::size_t size, std::span<const std::uint8_t>& value) {
        if (size > remaining()) {
            truncated_ = true;
            return false;
        }
        value = bytes_.subspan(position_, size);
        position_ += size;
        return true;
    }

    [[nodiscard]] std::size_t remaining() const noexcept {
        return bytes_.size() - position_;
    }

    [[nodiscard]] bool atEnd() const noexcept {
        return position_ == bytes_.size();
    }

    [[nodiscard]] bool truncated() const noexcept { return truncated_; }

private:
    template<typename Unsigned>
    bool integer(Unsigned& value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        if (sizeof(Unsigned) > remaining()) {
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
};

StrongPistonCheckpointPersistenceErrorCode mapStructureError(
    const StructureCheckpointPersistenceErrorCode code) {
    if (code == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
        return StrongPistonCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code == StructureCheckpointPersistenceErrorCode::TopologyMismatch) {
        return StrongPistonCheckpointPersistenceErrorCode::TopologyMismatch;
    }
    return StrongPistonCheckpointPersistenceErrorCode::InvalidData;
}

StrongPistonCheckpointPersistenceErrorCode mapFluidError(
    const fluid::MovingInterfaceFluidCheckpointErrorCode code) {
    if (code
        == fluid::MovingInterfaceFluidCheckpointErrorCode::LimitExceeded) {
        return StrongPistonCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code
        == fluid::MovingInterfaceFluidCheckpointErrorCode::ChecksumMismatch) {
        return StrongPistonCheckpointPersistenceErrorCode::ChecksumMismatch;
    }
    return StrongPistonCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const StrongPistonCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= checkpointEnvelopeBytes
        && limits.structure.maximumEncodedBytes > 0
        && limits.fluid.maximumBytes > 0;
}

} // namespace

bool serializeStrongPistonCheckpoint(
    const StrongCoupledPistonCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    StrongPistonCheckpointPersistenceError* error,
    const StrongPistonCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
            "strong piston checkpoint limits are invalid");
    }
    try {
        StrongCoupledPistonCase validator;
        validator.restore(checkpoint);

        std::vector<std::uint8_t> structureBytes;
        StructureCheckpointPersistenceError structureError;
        if (!serializeStructureCheckpoint(
                validator.structure(), checkpoint.structure,
                structureBytes, &structureError, limits.structure)) {
            return fail(
                error, mapStructureError(structureError.code),
                "strong piston Structure payload: "
                    + structureError.message);
        }
        std::vector<std::uint8_t> fluidBytes;
        fluid::MovingInterfaceFluidCheckpointError fluidError;
        if (!fluid::serializeMovingInterfaceFluidCheckpoint(
                checkpoint.fluid, fluidBytes, &fluidError, limits.fluid)) {
            return fail(
                error, mapFluidError(fluidError.code),
                "strong piston fluid payload: " + fluidError.message);
        }

        std::vector<std::uint8_t> payload;
        BufferWriter payloadWriter(
            payload, limits.maximumEncodedBytes - checkpointEnvelopeBytes);
        if (!payloadWriter.u32(checkpoint.version)
            || !payloadWriter.u64(
                checkpoint.interfaceDefinitionFingerprint)
            || !payloadWriter.u64(
                static_cast<std::uint64_t>(structureBytes.size()))
            || !payloadWriter.raw(structureBytes)
            || !payloadWriter.u64(
                static_cast<std::uint64_t>(fluidBytes.size()))
            || !payloadWriter.raw(fluidBytes)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
                "strong piston checkpoint exceeds its byte limit");
        }

        std::vector<std::uint8_t> encoded;
        BufferWriter writer(encoded, limits.maximumEncodedBytes);
        if (!writer.raw(checkpointMagic)
            || !writer.u16(strongPistonCheckpointProtocolVersion)
            || !writer.u16(0)
            || !writer.u64(static_cast<std::uint64_t>(payload.size()))
            || !writer.u64(checksum(payload))
            || !writer.raw(payload)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
                "strong piston checkpoint envelope exceeds its byte limit");
        }
        bytes.swap(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
            "strong piston checkpoint allocation exceeded available memory");
    } catch (const std::exception& exception) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::InvalidData,
            exception.what());
    }
}

bool deserializeStrongPistonCheckpoint(
    const std::span<const std::uint8_t> bytes,
    StrongCoupledPistonCheckpoint& checkpoint,
    StrongPistonCheckpointPersistenceError* error,
    const StrongPistonCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits) || bytes.size() > limits.maximumEncodedBytes) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
            "strong piston checkpoint exceeds its byte limit");
    }
    if (bytes.size() < checkpointEnvelopeBytes) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::Truncated,
            "strong piston checkpoint envelope is truncated");
    }
    try {
        BufferReader reader(bytes);
        std::span<const std::uint8_t> magic;
        std::uint16_t protocolVersion = 0;
        std::uint16_t reserved = 0;
        std::uint64_t payloadSize = 0;
        std::uint64_t storedChecksum = 0;
        if (!reader.raw(checkpointMagic.size(), magic)
            || !reader.u16(protocolVersion)
            || !reader.u16(reserved)
            || !reader.u64(payloadSize)
            || !reader.u64(storedChecksum)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston checkpoint envelope is truncated");
        }
        if (!std::equal(magic.begin(), magic.end(), checkpointMagic.begin())) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::InvalidMagic,
                "strong piston checkpoint magic is invalid");
        }
        if (protocolVersion != strongPistonCheckpointProtocolVersion
            || reserved != 0) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::UnsupportedVersion,
                "strong piston checkpoint protocol is unsupported");
        }
        if (payloadSize > std::numeric_limits<std::size_t>::max()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
                "strong piston checkpoint payload size is not representable");
        }
        const std::size_t expectedPayloadSize =
            static_cast<std::size_t>(payloadSize);
        if (expectedPayloadSize > reader.remaining()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston checkpoint payload is truncated");
        }
        if (expectedPayloadSize < reader.remaining()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::TrailingData,
                "strong piston checkpoint has trailing data");
        }
        std::span<const std::uint8_t> payload;
        if (!reader.raw(expectedPayloadSize, payload)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston checkpoint payload is truncated");
        }
        if (checksum(payload) != storedChecksum) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::ChecksumMismatch,
                "strong piston checkpoint checksum does not match");
        }

        BufferReader payloadReader(payload);
        StrongCoupledPistonCheckpoint candidate;
        std::uint64_t structureSize = 0;
        std::uint64_t fluidSize = 0;
        if (!payloadReader.u32(candidate.version)
            || !payloadReader.u64(
                candidate.interfaceDefinitionFingerprint)
            || !payloadReader.u64(structureSize)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston checkpoint metadata is truncated");
        }
        if (structureSize > limits.structure.maximumEncodedBytes
            || structureSize > std::numeric_limits<std::size_t>::max()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
                "strong piston Structure payload exceeds its limit");
        }
        std::span<const std::uint8_t> structureBytes;
        if (!payloadReader.raw(
                static_cast<std::size_t>(structureSize), structureBytes)
            || !payloadReader.u64(fluidSize)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston nested payload is truncated");
        }
        if (fluidSize > limits.fluid.maximumBytes
            || fluidSize > std::numeric_limits<std::size_t>::max()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
                "strong piston fluid payload exceeds its limit");
        }
        std::span<const std::uint8_t> fluidBytes;
        if (!payloadReader.raw(
                static_cast<std::size_t>(fluidSize), fluidBytes)) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::Truncated,
                "strong piston fluid payload is truncated");
        }
        if (!payloadReader.atEnd()) {
            return fail(
                error,
                StrongPistonCheckpointPersistenceErrorCode::TrailingData,
                "strong piston payload has trailing data");
        }

        StrongCoupledPistonCase validator;
        StructureCheckpointPersistenceError structureError;
        if (!deserializeStructureCheckpoint(
                structureBytes, validator.structure(), candidate.structure,
                &structureError, limits.structure)) {
            return fail(
                error, mapStructureError(structureError.code),
                "strong piston Structure payload: "
                    + structureError.message);
        }
        fluid::MovingInterfaceFluidCheckpointError fluidError;
        if (!fluid::deserializeMovingInterfaceFluidCheckpoint(
                fluidBytes, candidate.fluid, &fluidError, limits.fluid)) {
            return fail(
                error, mapFluidError(fluidError.code),
                "strong piston fluid payload: " + fluidError.message);
        }
        validator.restore(candidate);
        checkpoint = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
            "strong piston checkpoint allocation exceeded available memory");
    } catch (const std::exception& exception) {
        return fail(
            error,
            StrongPistonCheckpointPersistenceErrorCode::InvalidData,
            exception.what());
    }
}

} // namespace simwing::fsi
