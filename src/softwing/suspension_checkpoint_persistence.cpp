#include "softwing/suspension_checkpoint_persistence.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace softwing {
namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'S', 'C'};
constexpr std::uint32_t stateVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(SuspensionCheckpointPersistenceError* error,
          const SuspensionCheckpointPersistenceErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(SuspensionCheckpointPersistenceError* error) {
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

bool validUtf8(const std::string_view text) {
    std::size_t index = 0;
    while (index < text.size()) {
        const auto first = static_cast<unsigned char>(text[index]);
        std::size_t trailing = 0;
        std::uint32_t codePoint = 0;
        if (first <= 0x7fU) {
            ++index;
            continue;
        }
        if ((first & 0xe0U) == 0xc0U) {
            trailing = 1;
            codePoint = first & 0x1fU;
        } else if ((first & 0xf0U) == 0xe0U) {
            trailing = 2;
            codePoint = first & 0x0fU;
        } else if ((first & 0xf8U) == 0xf0U) {
            trailing = 3;
            codePoint = first & 0x07U;
        } else {
            return false;
        }
        if (trailing > text.size() - index - 1) return false;
        for (std::size_t offset = 1; offset <= trailing; ++offset) {
            const auto next = static_cast<unsigned char>(text[index + offset]);
            if ((next & 0xc0U) != 0x80U) return false;
            codePoint = (codePoint << 6U) | (next & 0x3fU);
        }
        if ((trailing == 1 && codePoint < 0x80U)
            || (trailing == 2 && codePoint < 0x800U)
            || (trailing == 3 && codePoint < 0x10000U)
            || codePoint > 0x10ffffU
            || (codePoint >= 0xd800U && codePoint <= 0xdfffU)) {
            return false;
        }
        index += trailing + 1;
    }
    return true;
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
        if (!std::isfinite(value)) {
            invalid_ = true;
            return false;
        }
        return u64(std::bit_cast<std::uint64_t>(value));
    }
    bool vector3(const Vec3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool bytes(const std::span<const std::uint8_t> values) {
        return raw(values.data(), values.size());
    }
    void markLimitExceeded() noexcept { exceeded_ = true; }
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
    bool invalid_ = false;
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
        if (!std::isfinite(value)) {
            invalid_ = true;
            return false;
        }
        return true;
    }
    bool vector3(Vec3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded)) return false;
        if (encoded > 1) {
            invalid_ = true;
            return false;
        }
        value = encoded != 0;
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
    void markLimitExceeded() noexcept { limitExceeded_ = true; }

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
    bool invalid_ = false;
};

SuspensionCheckpointPersistenceErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return SuspensionCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return SuspensionCheckpointPersistenceErrorCode::Truncated;
    }
    return SuspensionCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const SuspensionCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.maximumStringBytes <= limits.maximumTotalTextBytes;
}

bool addText(const std::size_t amount,
             std::size_t& total,
             const SuspensionCheckpointPersistenceLimits& limits) {
    if (amount > limits.maximumStringBytes
        || amount > limits.maximumTotalTextBytes - total) {
        return false;
    }
    total += amount;
    return true;
}

bool writeString(Writer& writer,
                 const std::string_view text,
                 std::size_t& totalText,
                 const SuspensionCheckpointPersistenceLimits& limits) {
    if (!validUtf8(text)) {
        return false;
    }
    if (!addText(text.size(), totalText, limits)) {
        writer.markLimitExceeded();
        return false;
    }
    if (!writer.count(text.size())) return false;
    return writer.bytes(std::span<const std::uint8_t>{
        reinterpret_cast<const std::uint8_t*>(text.data()), text.size()});
}

bool readString(Reader& reader,
                std::string& text,
                std::size_t& totalText,
                const SuspensionCheckpointPersistenceLimits& limits) {
    std::size_t size = 0;
    if (!reader.count(size, limits.maximumStringBytes)) {
        return false;
    }
    if (!addText(size, totalText, limits)) {
        reader.markLimitExceeded();
        return false;
    }
    const auto bytes = reader.remainingBytes(size);
    if (reader.truncated()) return false;
    text.assign(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    return validUtf8(text);
}

bool writeQuaternion(Writer& writer, const Quaternion& value) {
    return writer.finiteDouble(value.w) && writer.finiteDouble(value.x)
        && writer.finiteDouble(value.y) && writer.finiteDouble(value.z);
}

bool readQuaternion(Reader& reader, Quaternion& value) {
    return reader.finiteDouble(value.w) && reader.finiteDouble(value.x)
        && reader.finiteDouble(value.y) && reader.finiteDouble(value.z);
}

bool writePayloadState(Writer& writer, const RigidPayloadState& value) {
    return writer.vector3(value.centreOfMassWorld)
        && writeQuaternion(writer, value.orientation)
        && writer.vector3(value.linearVelocity)
        && writer.vector3(value.angularVelocity);
}

bool readPayloadState(Reader& reader, RigidPayloadState& value) {
    return reader.vector3(value.centreOfMassWorld)
        && readQuaternion(reader, value.orientation)
        && reader.vector3(value.linearVelocity)
        && reader.vector3(value.angularVelocity);
}

bool writeWrench(Writer& writer, const Wrench& value) {
    return writer.vector3(value.force) && writer.vector3(value.moment);
}

bool readWrench(Reader& reader, Wrench& value) {
    return reader.vector3(value.force) && reader.vector3(value.moment);
}

bool writeEndpoint(Writer& writer,
                   const SuspensionEndpoint& value,
                   std::size_t& totalText,
                   const SuspensionCheckpointPersistenceLimits& limits) {
    return writer.u8(static_cast<std::uint8_t>(value.kind))
        && writeString(writer, value.id, totalText, limits);
}

bool readEndpoint(Reader& reader,
                  SuspensionEndpoint& value,
                  std::size_t& totalText,
                  const SuspensionCheckpointPersistenceLimits& limits) {
    std::uint8_t kind = 0;
    if (!reader.u8(kind)
        || kind > static_cast<std::uint8_t>(
            SuspensionEndpointKind::HangPoint)
        || !readString(reader, value.id, totalText, limits)) {
        return false;
    }
    value.kind = static_cast<SuspensionEndpointKind>(kind);
    return true;
}

bool writeStringVector(
    Writer& writer,
    const std::vector<std::string>& values,
    std::size_t& totalText,
    const SuspensionCheckpointPersistenceLimits& limits) {
    if (values.size() > limits.maximumRecords
        || !writer.count(values.size())) {
        return false;
    }
    for (const std::string& value : values) {
        if (!writeString(writer, value, totalText, limits)) return false;
    }
    return true;
}

bool readStringVector(
    Reader& reader,
    std::vector<std::string>& values,
    std::size_t& totalText,
    const SuspensionCheckpointPersistenceLimits& limits) {
    std::size_t count = 0;
    if (!reader.count(count, limits.maximumRecords)) return false;
    values.clear();
    values.resize(count);
    for (std::string& value : values) {
        if (!readString(reader, value, totalText, limits)) return false;
    }
    return true;
}

bool writeSegmentDiagnostics(
    Writer& writer,
    const SuspensionSegmentDiagnostics& value,
    std::size_t& totalText,
    const SuspensionCheckpointPersistenceLimits& limits) {
    if (!writeString(writer, value.id, totalText, limits)
        || !writeEndpoint(writer, value.from, totalText, limits)
        || !writeEndpoint(writer, value.to, totalText, limits)
        || !writeStringVector(writer, value.paths, totalText, limits)
        || !writeStringVector(writer, value.groups, totalText, limits)
        || !writer.finiteDouble(value.length)
        || !writer.finiteDouble(value.commandedRestLength)
        || !writer.finiteDouble(value.stretch)
        || !writer.finiteDouble(value.strain)
        || !writer.u8(value.taut ? 1 : 0)
        || !writer.finiteDouble(value.multiplier)
        || !writer.finiteDouble(value.tension)
        || !writer.finiteDouble(value.residual)
        || !writer.finiteDouble(value.elasticEnergy)
        || !writer.vector3(value.dampingImpulse)
        || !writer.finiteDouble(value.dampingWork)
        || !writer.finiteDouble(value.controlWork)
        || !writer.vector3(value.fromImpulse)
        || !writer.vector3(value.toImpulse)
        || !writer.vector3(value.fromMoment)
        || !writer.vector3(value.toMoment)) {
        return false;
    }
    return true;
}

bool readSegmentDiagnostics(
    Reader& reader,
    SuspensionSegmentDiagnostics& value,
    std::size_t& totalText,
    const SuspensionCheckpointPersistenceLimits& limits) {
    return readString(reader, value.id, totalText, limits)
        && readEndpoint(reader, value.from, totalText, limits)
        && readEndpoint(reader, value.to, totalText, limits)
        && readStringVector(reader, value.paths, totalText, limits)
        && readStringVector(reader, value.groups, totalText, limits)
        && reader.finiteDouble(value.length)
        && reader.finiteDouble(value.commandedRestLength)
        && reader.finiteDouble(value.stretch)
        && reader.finiteDouble(value.strain)
        && reader.boolean(value.taut)
        && reader.finiteDouble(value.multiplier)
        && reader.finiteDouble(value.tension)
        && reader.finiteDouble(value.residual)
        && reader.finiteDouble(value.elasticEnergy)
        && reader.vector3(value.dampingImpulse)
        && reader.finiteDouble(value.dampingWork)
        && reader.finiteDouble(value.controlWork)
        && reader.vector3(value.fromImpulse)
        && reader.vector3(value.toImpulse)
        && reader.vector3(value.fromMoment)
        && reader.vector3(value.toMoment);
}

bool writePayloadDiagnostics(Writer& writer,
                             const PayloadDiagnostics& value) {
    return writePayloadState(writer, value.state)
        && writer.vector3(value.linearMomentum)
        && writer.vector3(value.angularMomentum)
        && writer.finiteDouble(value.translationalKineticEnergy)
        && writer.finiteDouble(value.rotationalKineticEnergy)
        && writer.finiteDouble(value.gravitationalEnergy)
        && writeWrench(writer, value.appliedWrench)
        && writeWrench(writer, value.lineWrench)
        && writeWrench(writer, value.groundWrench)
        && writeWrench(writer, value.anchorWrench);
}

bool readPayloadDiagnostics(Reader& reader, PayloadDiagnostics& value) {
    return readPayloadState(reader, value.state)
        && reader.vector3(value.linearMomentum)
        && reader.vector3(value.angularMomentum)
        && reader.finiteDouble(value.translationalKineticEnergy)
        && reader.finiteDouble(value.rotationalKineticEnergy)
        && reader.finiteDouble(value.gravitationalEnergy)
        && readWrench(reader, value.appliedWrench)
        && readWrench(reader, value.lineWrench)
        && readWrench(reader, value.groundWrench)
        && readWrench(reader, value.anchorWrench);
}

bool writeLoads(Writer& writer,
                const std::vector<std::pair<std::string, Vec3>>& values,
                std::size_t& totalText,
                const SuspensionCheckpointPersistenceLimits& limits) {
    if (values.size() > limits.maximumRecords
        || !writer.count(values.size())) {
        return false;
    }
    for (const auto& [id, load] : values) {
        if (!writeString(writer, id, totalText, limits)
            || !writer.vector3(load)) {
            return false;
        }
    }
    return true;
}

bool readLoads(Reader& reader,
               std::vector<std::pair<std::string, Vec3>>& values,
               std::size_t& totalText,
               const SuspensionCheckpointPersistenceLimits& limits) {
    std::size_t count = 0;
    if (!reader.count(count, limits.maximumRecords)) return false;
    values.clear();
    values.resize(count);
    for (auto& [id, load] : values) {
        if (!readString(reader, id, totalText, limits)
            || !reader.vector3(load)) {
            return false;
        }
    }
    return true;
}

bool writeDiagnostics(Writer& writer,
                      const SuspensionDiagnostics& value,
                      std::size_t& totalText,
                      const SuspensionCheckpointPersistenceLimits& limits) {
    return writer.u8(value.registered ? 1 : 0)
        && writer.u8(value.converged ? 1 : 0)
        && writer.u8(value.allSlack ? 1 : 0)
        && writer.u8(value.anchored ? 1 : 0)
        && writer.u8(value.grounded ? 1 : 0)
        && writer.u8(value.failedTrial ? 1 : 0)
        && writer.count(value.attachmentCount)
        && writer.count(value.junctionCount)
        && writer.count(value.segmentCount)
        && writer.count(value.tautCount)
        && writer.count(value.slackCount)
        && writer.count(value.solverIterations)
        && writer.finiteDouble(value.maximumResidual)
        && writer.finiteDouble(value.maximumGroundPenetration)
        && writer.finiteDouble(value.elasticEnergy)
        && writer.finiteDouble(value.dampingWork)
        && writer.finiteDouble(value.controlWork)
        && writer.vector3(value.netInternalImpulse)
        && writer.vector3(value.netInternalMoment)
        && writeWrench(writer, value.fixedSupportReaction)
        && writeWrench(writer, value.canopySupportReaction)
        && writeWrench(writer, value.groundReaction)
        && writeLoads(writer, value.attachmentLoads, totalText, limits)
        && writeLoads(writer, value.groupLoads, totalText, limits)
        && writeString(writer, value.provenance, totalText, limits)
        && writer.u8(static_cast<std::uint8_t>(value.failurePhase))
        && writeString(writer, value.failureEntity, totalText, limits);
}

bool readDiagnostics(Reader& reader,
                     SuspensionDiagnostics& value,
                     std::size_t& totalText,
                     const SuspensionCheckpointPersistenceLimits& limits) {
    std::uint8_t phase = 0;
    if (!reader.boolean(value.registered)
        || !reader.boolean(value.converged)
        || !reader.boolean(value.allSlack)
        || !reader.boolean(value.anchored)
        || !reader.boolean(value.grounded)
        || !reader.boolean(value.failedTrial)
        || !reader.count(value.attachmentCount, limits.maximumRecords)
        || !reader.count(value.junctionCount, limits.maximumRecords)
        || !reader.count(value.segmentCount, limits.maximumRecords)
        || !reader.count(value.tautCount, limits.maximumRecords)
        || !reader.count(value.slackCount, limits.maximumRecords)
        || !reader.count(value.solverIterations, limits.maximumRecords)
        || !reader.finiteDouble(value.maximumResidual)
        || !reader.finiteDouble(value.maximumGroundPenetration)
        || !reader.finiteDouble(value.elasticEnergy)
        || !reader.finiteDouble(value.dampingWork)
        || !reader.finiteDouble(value.controlWork)
        || !reader.vector3(value.netInternalImpulse)
        || !reader.vector3(value.netInternalMoment)
        || !readWrench(reader, value.fixedSupportReaction)
        || !readWrench(reader, value.canopySupportReaction)
        || !readWrench(reader, value.groundReaction)
        || !readLoads(reader, value.attachmentLoads, totalText, limits)
        || !readLoads(reader, value.groupLoads, totalText, limits)
        || !readString(reader, value.provenance, totalText, limits)
        || !reader.u8(phase)
        || phase > static_cast<std::uint8_t>(SuspensionPhase::Diagnostics)
        || !readString(reader, value.failureEntity, totalText, limits)) {
        return false;
    }
    value.failurePhase = static_cast<SuspensionPhase>(phase);
    return true;
}

bool exactCount(Reader& reader,
                const std::size_t expected,
                const SuspensionCheckpointPersistenceLimits& limits,
                bool& topologyMismatch) {
    std::size_t count = 0;
    if (!reader.count(count, limits.maximumRecords)) return false;
    if (count != expected) {
        topologyMismatch = true;
        return false;
    }
    return true;
}

bool writeState(Writer& writer,
                const SuspensionCheckpoint& value,
                std::size_t& totalText,
                const SuspensionCheckpointPersistenceLimits& limits) {
    if (!writer.u32(stateVersion) || !writer.u32(0)
        || !writer.count(value.schemaMajor)
        || !writer.u64(value.topologyFingerprint)
        || !writer.u64(value.stateFingerprint)
        || !writePayloadState(writer, value.payloadState)
        || !writePayloadState(writer, value.previousPayloadState)
        || value.commandedHangPointPositions.size() > limits.maximumRecords
        || !writer.count(value.commandedHangPointPositions.size())) {
        return false;
    }
    for (const Vec3& point : value.commandedHangPointPositions) {
        if (!writer.vector3(point)) return false;
    }
    if (value.segments.size() > limits.maximumRecords
        || !writer.count(value.segments.size())) {
        return false;
    }
    for (const SuspensionSegmentCheckpoint& segment : value.segments) {
        if (!writeString(writer, segment.id, totalText, limits)
            || !writer.finiteDouble(segment.commandedRestLength)
            || !writer.finiteDouble(segment.accumulatedLambda)) {
            return false;
        }
    }
    if (value.controls.size() > limits.maximumRecords
        || !writer.count(value.controls.size())) {
        return false;
    }
    for (const SuspensionControlState& control : value.controls) {
        if (!writeString(writer, control.id, totalText, limits)
            || !writer.finiteDouble(control.targetCommand)
            || !writer.finiteDouble(control.actualCommand)
            || !writer.finiteDouble(control.travel)
            || !writer.finiteDouble(control.signedWork)) {
            return false;
        }
    }
    const auto doubles = [&](const std::vector<double>& values) {
        if (values.size() > limits.maximumRecords
            || !writer.count(values.size())) {
            return false;
        }
        for (const double value : values) {
            if (!writer.finiteDouble(value)) return false;
        }
        return true;
    };
    if (!doubles(value.groundMultipliers)
        || !writeWrench(writer, value.appliedWrench)
        || !writer.vector3(value.currentGravity)
        || !doubles(value.pendingSegmentControlWork)
        || value.segmentDiagnostics.size() > limits.maximumRecords
        || !writer.count(value.segmentDiagnostics.size())) {
        return false;
    }
    for (const SuspensionSegmentDiagnostics& diagnostic :
         value.segmentDiagnostics) {
        if (!writeSegmentDiagnostics(writer, diagnostic, totalText, limits)) {
            return false;
        }
    }
    return writePayloadDiagnostics(writer, value.payloadDiagnostics)
        && writeDiagnostics(writer, value.diagnostics, totalText, limits)
        && writeDiagnostics(
            writer, value.committedDiagnostics, totalText, limits);
}

bool readState(Reader& reader,
               const SuspensionCheckpoint& topologyTemplate,
               SuspensionCheckpoint& value,
               std::size_t& totalText,
               const SuspensionCheckpointPersistenceLimits& limits,
               bool& topologyMismatch) {
    std::uint32_t version = 0;
    std::uint32_t reserved = 0;
    std::size_t schema = 0;
    if (!reader.u32(version) || !reader.u32(reserved)
        || !reader.count(schema, std::numeric_limits<std::size_t>::max())
        || !reader.u64(value.topologyFingerprint)
        || !reader.u64(value.stateFingerprint)) {
        return false;
    }
    if (version != stateVersion || reserved != 0) return false;
    if (schema != topologyTemplate.schemaMajor
        || value.topologyFingerprint
            != topologyTemplate.topologyFingerprint) {
        topologyMismatch = true;
        return false;
    }
    value.schemaMajor = schema;
    if (!readPayloadState(reader, value.payloadState)
        || !readPayloadState(reader, value.previousPayloadState)
        || !exactCount(reader,
                       topologyTemplate.commandedHangPointPositions.size(),
                       limits, topologyMismatch)) {
        return false;
    }
    value.commandedHangPointPositions.resize(
        topologyTemplate.commandedHangPointPositions.size());
    for (Vec3& point : value.commandedHangPointPositions) {
        if (!reader.vector3(point)) return false;
    }
    if (!exactCount(reader, topologyTemplate.segments.size(), limits,
                    topologyMismatch)) {
        return false;
    }
    value.segments.resize(topologyTemplate.segments.size());
    for (std::size_t index = 0; index < value.segments.size(); ++index) {
        auto& segment = value.segments[index];
        if (!readString(reader, segment.id, totalText, limits)
            || segment.id != topologyTemplate.segments[index].id
            || !reader.finiteDouble(segment.commandedRestLength)
            || !reader.finiteDouble(segment.accumulatedLambda)) {
            topologyMismatch = segment.id
                != topologyTemplate.segments[index].id;
            return false;
        }
    }
    if (!exactCount(reader, topologyTemplate.controls.size(), limits,
                    topologyMismatch)) {
        return false;
    }
    value.controls.resize(topologyTemplate.controls.size());
    for (std::size_t index = 0; index < value.controls.size(); ++index) {
        auto& control = value.controls[index];
        if (!readString(reader, control.id, totalText, limits)
            || control.id != topologyTemplate.controls[index].id
            || !reader.finiteDouble(control.targetCommand)
            || !reader.finiteDouble(control.actualCommand)
            || !reader.finiteDouble(control.travel)
            || !reader.finiteDouble(control.signedWork)) {
            topologyMismatch = control.id
                != topologyTemplate.controls[index].id;
            return false;
        }
    }
    const auto doubles = [&](std::vector<double>& values,
                             const std::size_t expected) {
        if (!exactCount(reader, expected, limits, topologyMismatch)) {
            return false;
        }
        values.resize(expected);
        for (double& item : values) {
            if (!reader.finiteDouble(item)) return false;
        }
        return true;
    };
    if (!doubles(value.groundMultipliers,
                 topologyTemplate.groundMultipliers.size())
        || !readWrench(reader, value.appliedWrench)
        || !reader.vector3(value.currentGravity)
        || !doubles(value.pendingSegmentControlWork,
                    topologyTemplate.pendingSegmentControlWork.size())
        || !exactCount(reader, topologyTemplate.segmentDiagnostics.size(),
                       limits, topologyMismatch)) {
        return false;
    }
    value.segmentDiagnostics.resize(
        topologyTemplate.segmentDiagnostics.size());
    for (SuspensionSegmentDiagnostics& diagnostic :
         value.segmentDiagnostics) {
        if (!readSegmentDiagnostics(reader, diagnostic, totalText, limits)) {
            return false;
        }
    }
    return readPayloadDiagnostics(reader, value.payloadDiagnostics)
        && readDiagnostics(reader, value.diagnostics, totalText, limits)
        && readDiagnostics(
            reader, value.committedDiagnostics, totalText, limits);
}

} // namespace

bool serializeSuspensionCheckpoint(
    const SuspensionCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    SuspensionCheckpointPersistenceError* error,
    const SuspensionCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)
        || checkpoint.schemaMajor != suspensionCheckpointSchemaMajor) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::InvalidData,
                    "invalid suspension checkpoint or persistence limits");
    }
    if (checkpoint.stateFingerprint
        != suspensionCheckpointStateFingerprint(checkpoint)) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::InvalidData,
                    "suspension checkpoint state fingerprint is invalid");
    }
    try {
        std::vector<std::uint8_t> payload;
        Writer payloadWriter(payload,
                             limits.maximumEncodedBytes - envelopeBytes);
        std::size_t totalText = 0;
        if (!writeState(payloadWriter, checkpoint, totalText, limits)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? SuspensionCheckpointPersistenceErrorCode::LimitExceeded
                    : SuspensionCheckpointPersistenceErrorCode::InvalidData,
                payloadWriter.exceeded()
                    ? "suspension checkpoint exceeds the encoded byte limit"
                    : "suspension checkpoint contains invalid state");
        }
        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t byte : magic) {
            if (!writer.u8(byte)) return false;
        }
        if (!writer.u16(suspensionCheckpointProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error,
                SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                "suspension checkpoint exceeds the encoded byte limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate suspension checkpoint encoding");
    } catch (const std::length_error&) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                    "suspension checkpoint encoding is too large");
    }
}

bool deserializeSuspensionCheckpoint(
    const std::span<const std::uint8_t> bytes,
    const SuspensionCheckpoint& topologyTemplate,
    SuspensionCheckpoint& checkpoint,
    SuspensionCheckpointPersistenceError* error,
    const SuspensionCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)
        || topologyTemplate.schemaMajor != suspensionCheckpointSchemaMajor) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::InvalidData,
                    "invalid topology template or persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                    "suspension checkpoint exceeds the encoded byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::Truncated,
                    "suspension checkpoint envelope is truncated");
    }
    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) return false;
        if (actual != expected) {
            return fail(error,
                        SuspensionCheckpointPersistenceErrorCode::InvalidMagic,
                        "suspension checkpoint magic is invalid");
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
                    "suspension checkpoint envelope is invalid");
    }
    if (protocol != suspensionCheckpointProtocolVersion) {
        return fail(error,
            SuspensionCheckpointPersistenceErrorCode::UnsupportedVersion,
            "suspension checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::InvalidData,
                    "suspension checkpoint reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::Truncated,
                    "suspension checkpoint payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::TrailingData,
                    "suspension checkpoint has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::ChecksumMismatch,
                    "suspension checkpoint checksum does not match");
    }
    try {
        SuspensionCheckpoint candidate;
        Reader reader(payload);
        std::size_t totalText = 0;
        bool topologyMismatch = false;
        if (!readState(reader, topologyTemplate, candidate, totalText, limits,
                       topologyMismatch)) {
            return fail(
                error,
                topologyMismatch
                    ? SuspensionCheckpointPersistenceErrorCode::TopologyMismatch
                    : readerErrorCode(reader),
                topologyMismatch
                    ? "suspension checkpoint belongs to a different topology"
                    : "suspension checkpoint state is invalid");
        }
        if (!reader.atEnd()) {
            return fail(error,
                        SuspensionCheckpointPersistenceErrorCode::TrailingData,
                        "suspension checkpoint payload has trailing data");
        }
        if (candidate.stateFingerprint
            != suspensionCheckpointStateFingerprint(candidate)) {
            return fail(error,
                        SuspensionCheckpointPersistenceErrorCode::InvalidData,
                        "suspension checkpoint state fingerprint is invalid");
        }
        checkpoint = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate decoded suspension checkpoint");
    } catch (const std::length_error&) {
        return fail(error,
                    SuspensionCheckpointPersistenceErrorCode::LimitExceeded,
                    "decoded suspension checkpoint is too large");
    }
}

} // namespace softwing
