#include "structure_checkpoint_persistence.h"

#include "structure_checkpoint_detail.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {

struct StructureCheckpointPersistenceAccess {
    [[nodiscard]] static const StructureCheckpoint::Detail* detail(
        const StructureCheckpoint& checkpoint) noexcept {
        return checkpoint.detail.get();
    }

    static void setDetail(
        StructureCheckpoint& checkpoint,
        std::shared_ptr<const StructureCheckpoint::Detail> detail) {
        checkpoint.detail = std::move(detail);
    }
};

namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'S', 'T'};
constexpr std::uint32_t stateVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t nodeRecordBytes = 72;
constexpr std::size_t forceRecordBytes = 24;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(StructureCheckpointPersistenceError* error,
          const StructureCheckpointPersistenceErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(StructureCheckpointPersistenceError* error) {
    if (error != nullptr) *error = {};
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
    bool u8(const std::uint8_t value) { return raw(&value, sizeof(value)); }
    bool u16(const std::uint16_t value) { return integer(value); }
    bool u32(const std::uint32_t value) { return integer(value); }
    bool u64(const std::uint64_t value) { return integer(value); }
    bool count(const std::size_t value) {
        return u64(static_cast<std::uint64_t>(value));
    }
    bool finiteDouble(const double value) {
        if (!std::isfinite(value)) return false;
        return u64(std::bit_cast<std::uint64_t>(value));
    }
    bool vector3(const StructureVector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool bytes(const std::span<const std::uint8_t> value) {
        return raw(value.data(), value.size());
    }
    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }

private:
    template<typename Unsigned>
    bool integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        std::array<std::uint8_t, sizeof(Unsigned)> encoded{};
        for (std::size_t index = 0; index < encoded.size(); ++index) {
            encoded[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return raw(encoded.data(), encoded.size());
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
    bool count(std::size_t& value, const std::size_t maximum) {
        std::uint64_t encoded = 0;
        if (!u64(encoded)) return false;
        if (encoded > maximum
            || encoded > std::numeric_limits<std::size_t>::max()) {
            limitExceeded_ = true;
            return false;
        }
        value = static_cast<std::size_t>(encoded);
        return true;
    }
    bool finiteDouble(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) return false;
        value = std::bit_cast<double>(bits);
        return std::isfinite(value);
    }
    bool vector3(StructureVector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded) || encoded > 1) return false;
        value = encoded != 0;
        return true;
    }
    bool fixedRecords(const std::size_t count,
                      const std::size_t bytesPerRecord) {
        if (count > (bytes_.size() - position_) / bytesPerRecord) {
            truncated_ = true;
            return false;
        }
        return true;
    }
    [[nodiscard]] std::span<const std::uint8_t> remainingBytes(
        const std::size_t count) {
        if (count > bytes_.size() - position_) {
            truncated_ = true;
            return {};
        }
        const auto result = bytes_.subspan(position_, count);
        position_ += count;
        return result;
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

StructureCheckpointPersistenceErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return StructureCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return StructureCheckpointPersistenceErrorCode::Truncated;
    }
    return StructureCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const StructureCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.maximumNodes > 0;
}

StructureCheckpointPersistenceErrorCode bodyErrorCode(
    const softwing::SoftBodyCheckpointPersistenceErrorCode code) {
    using Source = softwing::SoftBodyCheckpointPersistenceErrorCode;
    if (code == Source::LimitExceeded) {
        return StructureCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code == Source::TopologyMismatch) {
        return StructureCheckpointPersistenceErrorCode::TopologyMismatch;
    }
    return StructureCheckpointPersistenceErrorCode::InvalidData;
}

StructureCheckpointPersistenceErrorCode suspensionErrorCode(
    const softwing::SuspensionCheckpointPersistenceErrorCode code) {
    using Source = softwing::SuspensionCheckpointPersistenceErrorCode;
    if (code == Source::LimitExceeded) {
        return StructureCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code == Source::TopologyMismatch) {
        return StructureCheckpointPersistenceErrorCode::TopologyMismatch;
    }
    return StructureCheckpointPersistenceErrorCode::InvalidData;
}

bool writePayload(Writer& writer,
                  const StructureCheckpoint& checkpoint,
                  const std::vector<std::uint8_t>& body,
                  const std::vector<std::uint8_t>& suspension) {
    if (!writer.u32(stateVersion) || !writer.u32(0)
        || !writer.u32(checkpoint.version) || !writer.u32(0)
        || !writer.u64(checkpoint.definitionFingerprint)
        || !writer.u64(checkpoint.acceptedStepCount)
        || !writer.finiteDouble(checkpoint.simulationTimeSeconds)
        || !writer.count(checkpoint.nodes.size())) {
        return false;
    }
    for (const StructureNodeState& node : checkpoint.nodes) {
        if (!writer.vector3(node.positionMeters)
            || !writer.vector3(node.previousPositionMeters)
            || !writer.vector3(node.velocityMetersPerSecond)) {
            return false;
        }
    }
    if (!writer.count(checkpoint.pendingExternalForcesNewtons.size())) {
        return false;
    }
    for (const StructureVector3& force :
         checkpoint.pendingExternalForcesNewtons) {
        if (!writer.vector3(force)) return false;
    }
    return writer.vector3(checkpoint.lastAppliedExternalForceNewtons)
        && writer.count(body.size()) && writer.bytes(body)
        && writer.u8(suspension.empty() ? 0 : 1)
        && (suspension.empty()
                || (writer.count(suspension.size())
                    && writer.bytes(suspension)));
}

} // namespace

bool serializeStructureCheckpoint(
    const Structure& owner,
    const StructureCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    StructureCheckpointPersistenceError* error,
    const StructureCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)
        || checkpoint.definitionFingerprint
            != owner.definitionFingerprint()) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::InvalidData,
                    "invalid Structure checkpoint or persistence limits");
    }
    try {
        Structure validator(owner.definition());
        validator.restore(checkpoint);
        const StructureCheckpoint canonical = validator.checkpoint();
        const StructureCheckpoint::Detail* detail =
            StructureCheckpointPersistenceAccess::detail(canonical);
        if (detail == nullptr || canonical.nodes.size() > limits.maximumNodes) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::LimitExceeded,
                        "Structure checkpoint exceeds persistence limits");
        }

        std::vector<std::uint8_t> bodyBytes;
        softwing::SoftBodyCheckpointPersistenceError bodyError;
        if (!softwing::serializeSoftBodyCheckpoint(
                detail->body, bodyBytes, &bodyError, limits.body)) {
            return fail(error, bodyErrorCode(bodyError.code),
                        "cannot encode Structure body checkpoint: "
                            + bodyError.message);
        }
        std::vector<std::uint8_t> suspensionBytes;
        if (detail->suspension) {
            softwing::SuspensionCheckpointPersistenceError suspensionError;
            if (!softwing::serializeSuspensionCheckpoint(
                    *detail->suspension, suspensionBytes,
                    &suspensionError, limits.suspension)) {
                return fail(error,
                            suspensionErrorCode(suspensionError.code),
                            "cannot encode Structure suspension checkpoint: "
                                + suspensionError.message);
            }
        }

        std::vector<std::uint8_t> payload;
        Writer payloadWriter(payload,
                             limits.maximumEncodedBytes - envelopeBytes);
        if (!writePayload(payloadWriter, canonical, bodyBytes,
                          suspensionBytes)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? StructureCheckpointPersistenceErrorCode::LimitExceeded
                    : StructureCheckpointPersistenceErrorCode::InvalidData,
                "Structure checkpoint payload is invalid or too large");
        }
        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t byte : magic) {
            if (!writer.u8(byte)) return false;
        }
        if (!writer.u16(structureCheckpointProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::LimitExceeded,
                        "Structure checkpoint exceeds the encoded byte limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate Structure checkpoint encoding");
    } catch (const std::length_error&) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::LimitExceeded,
                    "Structure checkpoint encoding is too large");
    } catch (const std::exception& exception) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::InvalidData,
                    std::string("Structure checkpoint validation failed: ")
                        + exception.what());
    }
}

bool deserializeStructureCheckpoint(
    const std::span<const std::uint8_t> bytes,
    const Structure& owner,
    StructureCheckpoint& checkpoint,
    StructureCheckpointPersistenceError* error,
    const StructureCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::InvalidData,
                    "invalid Structure checkpoint persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::LimitExceeded,
                    "Structure checkpoint exceeds the encoded byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::Truncated,
                    "Structure checkpoint envelope is truncated");
    }

    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) return false;
        if (actual != expected) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::InvalidMagic,
                        "Structure checkpoint magic is invalid");
        }
    }
    std::uint16_t protocol = 0;
    std::uint16_t reserved = 0;
    std::size_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(protocol) || !envelope.u16(reserved)
        || !envelope.count(payloadSize,
                           limits.maximumEncodedBytes - envelopeBytes)
        || !envelope.u64(expectedChecksum)) {
        return fail(error, readerErrorCode(envelope),
                    "Structure checkpoint envelope is invalid");
    }
    if (protocol != structureCheckpointProtocolVersion) {
        return fail(error,
            StructureCheckpointPersistenceErrorCode::UnsupportedVersion,
            "Structure checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::InvalidData,
                    "Structure checkpoint reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::Truncated,
                    "Structure checkpoint payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::TrailingData,
                    "Structure checkpoint has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::ChecksumMismatch,
                    "Structure checkpoint checksum does not match");
    }

    try {
        Structure validator(owner.definition());
        const StructureCheckpoint topologyTemplate = validator.checkpoint();
        const StructureCheckpoint::Detail* templateDetail =
            StructureCheckpointPersistenceAccess::detail(topologyTemplate);
        if (templateDetail == nullptr) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::InvalidData,
                        "Structure topology template is incomplete");
        }

        Reader reader(payload);
        StructureCheckpoint candidate;
        std::uint32_t payloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t checkpointReserved = 0;
        if (!reader.u32(payloadVersion) || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(checkpointReserved)
            || !reader.u64(candidate.definitionFingerprint)
            || !reader.u64(candidate.acceptedStepCount)
            || !reader.finiteDouble(candidate.simulationTimeSeconds)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint payload header is invalid");
        }
        if (payloadVersion != stateVersion
            || candidate.version != structureCheckpointVersion) {
            return fail(error,
                StructureCheckpointPersistenceErrorCode::UnsupportedVersion,
                "Structure checkpoint state version is unsupported");
        }
        if (payloadReserved != 0 || checkpointReserved != 0) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::InvalidData,
                        "Structure checkpoint payload reserved bits are nonzero");
        }
        if (candidate.definitionFingerprint
            != owner.definitionFingerprint()) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::TopologyMismatch,
                        "Structure checkpoint belongs to a different definition");
        }

        const std::size_t expectedNodes = owner.definition().nodes.size();
        std::size_t nodeCount = 0;
        if (!reader.count(nodeCount, limits.maximumNodes)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint node count is invalid");
        }
        if (nodeCount != expectedNodes) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::TopologyMismatch,
                        "Structure checkpoint node count does not match topology");
        }
        if (!reader.fixedRecords(nodeCount, nodeRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint node state is truncated");
        }
        candidate.nodes.resize(nodeCount);
        for (StructureNodeState& node : candidate.nodes) {
            if (!reader.vector3(node.positionMeters)
                || !reader.vector3(node.previousPositionMeters)
                || !reader.vector3(node.velocityMetersPerSecond)) {
                return fail(error, readerErrorCode(reader),
                            "Structure checkpoint node state is invalid");
            }
        }
        std::size_t forceCount = 0;
        if (!reader.count(forceCount, limits.maximumNodes)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint force count is invalid");
        }
        if (forceCount != expectedNodes) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::TopologyMismatch,
                        "Structure checkpoint force count does not match topology");
        }
        if (!reader.fixedRecords(forceCount, forceRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint force state is truncated");
        }
        candidate.pendingExternalForcesNewtons.resize(forceCount);
        for (StructureVector3& force :
             candidate.pendingExternalForcesNewtons) {
            if (!reader.vector3(force)) {
                return fail(error, readerErrorCode(reader),
                            "Structure checkpoint force state is invalid");
            }
        }
        if (!reader.vector3(candidate.lastAppliedExternalForceNewtons)) {
            return fail(error, readerErrorCode(reader),
                        "Structure checkpoint applied force is invalid");
        }

        std::size_t bodySize = 0;
        if (!reader.count(bodySize, limits.body.maximumEncodedBytes)) {
            return fail(error, readerErrorCode(reader),
                        "Structure body checkpoint length is invalid");
        }
        const auto bodyBytes = reader.remainingBytes(bodySize);
        if (reader.truncated()) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::Truncated,
                        "Structure body checkpoint is truncated");
        }
        softwing::SoftBodyCheckpoint body;
        softwing::SoftBodyCheckpointPersistenceError bodyError;
        if (!softwing::deserializeSoftBodyCheckpoint(
                bodyBytes, templateDetail->body, body, &bodyError,
                limits.body)) {
            return fail(error, bodyErrorCode(bodyError.code),
                        "cannot decode Structure body checkpoint: "
                            + bodyError.message);
        }

        bool hasSuspension = false;
        if (!reader.boolean(hasSuspension)) {
            return fail(error, readerErrorCode(reader),
                        "Structure suspension marker is invalid");
        }
        if (hasSuspension != templateDetail->suspension.has_value()) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::TopologyMismatch,
                        "Structure suspension topology does not match");
        }
        std::optional<softwing::SuspensionCheckpoint> suspension;
        if (hasSuspension) {
            std::size_t suspensionSize = 0;
            if (!reader.count(
                    suspensionSize,
                    limits.suspension.maximumEncodedBytes)) {
                return fail(error, readerErrorCode(reader),
                            "Structure suspension checkpoint length is invalid");
            }
            const auto suspensionBytes =
                reader.remainingBytes(suspensionSize);
            if (reader.truncated()) {
                return fail(error,
                            StructureCheckpointPersistenceErrorCode::Truncated,
                            "Structure suspension checkpoint is truncated");
            }
            softwing::SuspensionCheckpoint decoded;
            softwing::SuspensionCheckpointPersistenceError suspensionError;
            if (!softwing::deserializeSuspensionCheckpoint(
                    suspensionBytes, *templateDetail->suspension, decoded,
                    &suspensionError, limits.suspension)) {
                return fail(error,
                            suspensionErrorCode(suspensionError.code),
                            "cannot decode Structure suspension checkpoint: "
                                + suspensionError.message);
            }
            suspension = std::move(decoded);
        }
        if (!reader.atEnd()) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::TrailingData,
                        "Structure checkpoint payload has trailing data");
        }

        auto detail = std::make_shared<StructureCheckpoint::Detail>();
        detail->body = std::move(body);
        detail->suspension = std::move(suspension);
        detail->publicNodes = candidate.nodes;
        StructureCheckpointPersistenceAccess::setDetail(
            candidate, std::move(detail));
        validator.restore(candidate);
        const StructureCheckpoint canonical = validator.checkpoint();
        if (canonical.nodes != candidate.nodes
            || canonical.pendingExternalForcesNewtons
                != candidate.pendingExternalForcesNewtons
            || canonical.lastAppliedExternalForceNewtons
                != candidate.lastAppliedExternalForceNewtons) {
            return fail(error,
                        StructureCheckpointPersistenceErrorCode::InvalidData,
                        "Structure checkpoint public and solver state disagree");
        }
        checkpoint = canonical;
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate decoded Structure checkpoint");
    } catch (const std::length_error&) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::LimitExceeded,
                    "decoded Structure checkpoint is too large");
    } catch (const std::exception& exception) {
        return fail(error,
                    StructureCheckpointPersistenceErrorCode::InvalidData,
                    std::string("Structure checkpoint validation failed: ")
                        + exception.what());
    }
}

} // namespace simwing::fsi
