#include "softwing/checkpoint_persistence.h"

#include "soft_body_checkpoint_detail.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace softwing {

struct SoftBodyCheckpointPersistenceAccess {
    [[nodiscard]] static const SoftBodyCheckpoint::State* state(
        const SoftBodyCheckpoint& checkpoint) noexcept {
        return checkpoint.state_.get();
    }

    [[nodiscard]] static SoftBodyCheckpoint make(
        std::shared_ptr<const SoftBodyCheckpoint::State> state) {
        return SoftBodyCheckpoint(std::move(state));
    }
};

namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'B', 'C'};
constexpr std::uint32_t stateVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t nodeRecordBytes = 96;
constexpr std::size_t triangleRecordBytes = 8;
constexpr std::size_t constraintRecordBytes = 8;
constexpr std::size_t membraneRecordBytes = 56;
constexpr std::size_t dihedralRecordBytes = 8;
constexpr std::size_t contactKeyBytes = 64;
constexpr std::size_t contactMultiplierRecordBytes = 72;
constexpr std::size_t contactRecordBytes = 634;
constexpr std::size_t contactDiagnosticsBytes = 307;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(SoftBodyCheckpointPersistenceError* error,
          const SoftBodyCheckpointPersistenceErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(SoftBodyCheckpointPersistenceError* error) {
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

    bool u8(const std::uint8_t value) { return raw(&value, sizeof(value)); }
    bool u16(const std::uint16_t value) { return integer(value); }
    bool u32(const std::uint32_t value) { return integer(value); }
    bool u64(const std::uint64_t value) { return integer(value); }

    bool i32(const int value) {
        if (value < std::numeric_limits<std::int32_t>::min()
            || value > std::numeric_limits<std::int32_t>::max()) {
            invalid_ = true;
            return false;
        }
        return u32(std::bit_cast<std::uint32_t>(
            static_cast<std::int32_t>(value)));
    }

    bool count(const std::size_t value) {
        if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
            if (value > std::numeric_limits<std::uint64_t>::max()) {
                invalid_ = true;
                return false;
            }
        }
        return u64(static_cast<std::uint64_t>(value));
    }

    bool finiteDouble(const double value) {
        if (!std::isfinite(value)) {
            invalid_ = true;
            return false;
        }
        return rawDouble(value);
    }

    bool diagnosticGap(const double value) {
        if (!std::isfinite(value)
            && !(std::isinf(value) && value > 0.0)) {
            invalid_ = true;
            return false;
        }
        return rawDouble(value);
    }

    bool vector3(const Vec3& value) {
        return finiteDouble(value.x)
            && finiteDouble(value.y)
            && finiteDouble(value.z);
    }

    bool bytes(const std::span<const std::uint8_t> values) {
        return raw(values.data(), values.size());
    }

    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }
    [[nodiscard]] bool invalid() const noexcept { return invalid_; }

private:
    bool rawDouble(const double value) {
        return u64(std::bit_cast<std::uint64_t>(value));
    }

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

    bool i32(int& value) {
        std::uint32_t bits = 0;
        if (!u32(bits)) {
            return false;
        }
        value = static_cast<int>(std::bit_cast<std::int32_t>(bits));
        return true;
    }

    bool count(std::size_t& value, const std::size_t maximum) {
        std::uint64_t encoded = 0;
        if (!u64(encoded)) {
            return false;
        }
        if (encoded > maximum
            || encoded > std::numeric_limits<std::size_t>::max()) {
            limitExceeded_ = true;
            return false;
        }
        value = static_cast<std::size_t>(encoded);
        return true;
    }

    bool finiteDouble(double& value) {
        if (!rawDouble(value)) {
            return false;
        }
        if (!std::isfinite(value)) {
            invalid_ = true;
            return false;
        }
        return true;
    }

    bool diagnosticGap(double& value) {
        if (!rawDouble(value)) {
            return false;
        }
        if (!std::isfinite(value)
            && !(std::isinf(value) && value > 0.0)) {
            invalid_ = true;
            return false;
        }
        return true;
    }

    bool vector3(Vec3& value) {
        return finiteDouble(value.x)
            && finiteDouble(value.y)
            && finiteDouble(value.z);
    }

    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded)) {
            return false;
        }
        if (encoded > 1) {
            invalid_ = true;
            return false;
        }
        value = encoded != 0;
        return true;
    }

    bool zeroBytes(const std::size_t count) {
        for (std::size_t index = 0; index < count; ++index) {
            std::uint8_t value = 0;
            if (!u8(value)) {
                return false;
            }
            if (value != 0) {
                invalid_ = true;
                return false;
            }
        }
        return true;
    }

    bool fixedRecords(const std::size_t count,
                      const std::size_t bytesPerRecord) {
        if (bytesPerRecord == 0
            || count > (bytes_.size() - position_) / bytesPerRecord) {
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
    [[nodiscard]] bool invalid() const noexcept { return invalid_; }

private:
    bool rawDouble(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) {
            return false;
        }
        value = std::bit_cast<double>(bits);
        return true;
    }

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

SoftBodyCheckpointPersistenceErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return SoftBodyCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return SoftBodyCheckpointPersistenceErrorCode::Truncated;
    }
    return SoftBodyCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const SoftBodyCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.maximumNodes > 0
        && limits.maximumTriangles > 0;
}

bool withinStaticLimits(
    const SoftBodyCheckpoint::State& state,
    const SoftBodyCheckpointPersistenceLimits& limits) {
    return state.nodes.size() <= limits.maximumNodes
        && state.triangles.size() <= limits.maximumTriangles
        && state.constraints.size() <= limits.maximumConstraints
        && state.membranes.size() <= limits.maximumMembranes
        && state.dihedrals.size() <= limits.maximumDihedrals
        && state.contactPairDiagnostics.size()
            <= limits.maximumContactPairs;
}

bool withinDynamicLimits(
    const SoftBodyCheckpoint::State& state,
    const SoftBodyCheckpointPersistenceLimits& limits) {
    if (state.contactMultipliers.size()
            > limits.maximumContactMultipliers
        || state.contactRecords.size() > limits.maximumContactRecords) {
        return false;
    }
    std::size_t auditCount = 0;
    for (const std::size_t count : {
             state.iterationCandidateKeys.size(),
             state.iterationQueryKeys.size(),
             state.certificationCandidateKeys.size(),
             state.certificationQueryKeys.size()}) {
        if (count > limits.maximumAuditKeys - auditCount) {
            return false;
        }
        auditCount += count;
    }
    return true;
}

bool writeKey(Writer& writer, const ContactFeatureKey& value) {
    return writer.count(value.pair)
        && writer.u8(static_cast<std::uint8_t>(value.kind))
        && writer.u8(0) && writer.u8(0) && writer.u8(0)
        && writer.u8(0) && writer.u8(0) && writer.u8(0) && writer.u8(0)
        && writer.count(value.firstPrimitive[0])
        && writer.count(value.firstPrimitive[1])
        && writer.count(value.firstPrimitive[2])
        && writer.count(value.secondPrimitive[0])
        && writer.count(value.secondPrimitive[1])
        && writer.count(value.secondPrimitive[2]);
}

bool readKey(Reader& reader, ContactFeatureKey& value) {
    std::uint8_t kind = 0;
    return reader.count(value.pair,
                        std::numeric_limits<std::size_t>::max())
        && reader.u8(kind)
        && kind <= static_cast<std::uint8_t>(
            ContactFeatureKind::SegmentTriangle)
        && (value.kind = static_cast<ContactFeatureKind>(kind), true)
        && reader.zeroBytes(7)
        && reader.count(value.firstPrimitive[0],
                        std::numeric_limits<std::size_t>::max())
        && reader.count(value.firstPrimitive[1],
                        std::numeric_limits<std::size_t>::max())
        && reader.count(value.firstPrimitive[2],
                        std::numeric_limits<std::size_t>::max())
        && reader.count(value.secondPrimitive[0],
                        std::numeric_limits<std::size_t>::max())
        && reader.count(value.secondPrimitive[1],
                        std::numeric_limits<std::size_t>::max())
        && reader.count(value.secondPrimitive[2],
                        std::numeric_limits<std::size_t>::max());
}

bool validKey(const ContactFeatureKey& key,
              const SoftBodyCheckpoint::State& state) {
    if (key.pair >= state.contactPairDiagnostics.size()) {
        return false;
    }
    const auto node = [&](const std::size_t index) {
        return index < state.nodes.size();
    };
    const auto triangle = [&](const std::size_t index) {
        return index < state.triangles.size();
    };
    switch (key.kind) {
    case ContactFeatureKind::VertexTriangle:
        return node(key.firstPrimitive[0])
            && triangle(key.secondPrimitive[0]);
    case ContactFeatureKind::EdgeEdge:
        return node(key.firstPrimitive[0])
            && node(key.firstPrimitive[1])
            && node(key.secondPrimitive[0])
            && node(key.secondPrimitive[1]);
    case ContactFeatureKind::SegmentTriangle:
        return node(key.firstPrimitive[0])
            && node(key.firstPrimitive[1])
            && triangle(key.secondPrimitive[0]);
    }
    return false;
}

bool writeDiagnostics(Writer& writer, const ContactDiagnostics& value) {
    return writer.u8(value.registered ? 1 : 0)
        && writer.u8(value.solveSucceeded ? 1 : 0)
        && writer.count(value.possibleCount)
        && writer.count(value.candidateCount)
        && writer.count(value.excludedCount)
        && writer.count(value.queryCount)
        && writer.count(value.activeCount)
        && writer.count(value.indeterminateCount)
        && writer.diagnosticGap(value.minimumGap)
        && writer.finiteDouble(value.maximumPenetration)
        && writer.finiteDouble(value.maximumNormalResidual)
        && writer.finiteDouble(value.maximumFrictionResidual)
        && writer.vector3(value.firstImpulse)
        && writer.vector3(value.secondImpulse)
        && writer.vector3(value.netInternalImpulse)
        && writer.vector3(value.firstMoment)
        && writer.vector3(value.secondMoment)
        && writer.vector3(value.netInternalMoment)
        && writer.finiteDouble(value.frictionWork)
        && writer.i32(value.maximumCcdIterations)
        && writer.i32(value.maximumIntervalSubdivisions)
        && writer.u8(value.hasFailure ? 1 : 0)
        && writeKey(writer, value.failureKey);
}

bool readDiagnostics(Reader& reader,
                     ContactDiagnostics& value,
                     const SoftBodyCheckpointPersistenceLimits& limits,
                     const SoftBodyCheckpoint::State& state) {
    if (!reader.boolean(value.registered)
        || !reader.boolean(value.solveSucceeded)
        || !reader.count(value.possibleCount, limits.maximumAuditKeys)
        || !reader.count(value.candidateCount, limits.maximumAuditKeys)
        || !reader.count(value.excludedCount, limits.maximumAuditKeys)
        || !reader.count(value.queryCount, limits.maximumAuditKeys)
        || !reader.count(value.activeCount, limits.maximumContactRecords)
        || !reader.count(value.indeterminateCount,
                         limits.maximumContactRecords)
        || !reader.diagnosticGap(value.minimumGap)
        || !reader.finiteDouble(value.maximumPenetration)
        || !reader.finiteDouble(value.maximumNormalResidual)
        || !reader.finiteDouble(value.maximumFrictionResidual)
        || !reader.vector3(value.firstImpulse)
        || !reader.vector3(value.secondImpulse)
        || !reader.vector3(value.netInternalImpulse)
        || !reader.vector3(value.firstMoment)
        || !reader.vector3(value.secondMoment)
        || !reader.vector3(value.netInternalMoment)
        || !reader.finiteDouble(value.frictionWork)
        || !reader.i32(value.maximumCcdIterations)
        || !reader.i32(value.maximumIntervalSubdivisions)
        || !reader.boolean(value.hasFailure)
        || !readKey(reader, value.failureKey)) {
        return false;
    }
    return value.maximumCcdIterations >= 0
        && value.maximumIntervalSubdivisions >= 0
        && (!value.hasFailure || validKey(value.failureKey, state));
}

bool writeRecord(Writer& writer, const ContactRecord& value) {
    if (!writeKey(writer, value.key)
        || !writer.u8(static_cast<std::uint8_t>(value.kind))
        || !writer.u8(static_cast<std::uint8_t>(value.ccdState))
        || !writer.u8(value.usedIntervalFallback ? 1 : 0)
        || !writer.u8(static_cast<std::uint8_t>(value.frictionState))
        || !writer.finiteDouble(value.timeOfImpact)
        || !writer.finiteDouble(value.bracketLower)
        || !writer.finiteDouble(value.bracketUpper)
        || !writer.i32(value.ccdIterations)
        || !writer.i32(value.intervalSubdivisions)) {
        return false;
    }
    for (const std::size_t node : value.firstNodes) {
        if (!writer.count(node)) return false;
    }
    for (const std::size_t node : value.secondNodes) {
        if (!writer.count(node)) return false;
    }
    if (value.firstNodeCount > 3 || value.secondNodeCount > 3
        || !writer.u8(static_cast<std::uint8_t>(value.firstNodeCount))
        || !writer.u8(static_cast<std::uint8_t>(value.secondNodeCount))) {
        return false;
    }
    for (const double weight : value.firstWeights) {
        if (!writer.finiteDouble(weight)) return false;
    }
    for (const double weight : value.secondWeights) {
        if (!writer.finiteDouble(weight)) return false;
    }
    if (!writer.vector3(value.firstPoint)
        || !writer.vector3(value.secondPoint)
        || !writer.vector3(value.normal)
        || !writer.finiteDouble(value.pairSeparation)
        || !writer.finiteDouble(value.gap)
        || !writer.finiteDouble(value.penetration)
        || !writer.finiteDouble(value.normalMultiplier)
        || !writer.finiteDouble(value.normalForceEstimate)
        || !writer.finiteDouble(value.normalImpulseMagnitude)
        || !writer.vector3(value.tangentialImpulse)
        || !writer.finiteDouble(value.frictionConeRatio)
        || !writer.finiteDouble(value.normalResidual)
        || !writer.finiteDouble(value.frictionResidual)
        || !writer.finiteDouble(value.frictionWork)
        || !writer.finiteDouble(value.tangentSpeedBefore)
        || !writer.finiteDouble(value.tangentSpeedAfter)
        || !writer.vector3(value.firstImpulse)
        || !writer.vector3(value.secondImpulse)
        || !writer.vector3(value.firstMoment)
        || !writer.vector3(value.secondMoment)) {
        return false;
    }
    for (const Vec3& impulse : value.firstNodeImpulses) {
        if (!writer.vector3(impulse)) return false;
    }
    for (const Vec3& impulse : value.secondNodeImpulses) {
        if (!writer.vector3(impulse)) return false;
    }
    return writer.i32(value.solverVisits);
}

bool readRecord(Reader& reader,
                ContactRecord& value,
                const SoftBodyCheckpoint::State& state) {
    std::uint8_t kind = 0;
    std::uint8_t ccdState = 0;
    std::uint8_t frictionState = 0;
    std::uint8_t firstCount = 0;
    std::uint8_t secondCount = 0;
    if (!readKey(reader, value.key)
        || !reader.u8(kind)
        || kind > static_cast<std::uint8_t>(
            ContactFeatureKind::SegmentTriangle)
        || !reader.u8(ccdState)
        || ccdState > static_cast<std::uint8_t>(CcdState::Indeterminate)
        || !reader.boolean(value.usedIntervalFallback)
        || !reader.u8(frictionState)
        || frictionState > static_cast<std::uint8_t>(
            ContactFrictionState::Sliding)
        || !reader.finiteDouble(value.timeOfImpact)
        || !reader.finiteDouble(value.bracketLower)
        || !reader.finiteDouble(value.bracketUpper)
        || !reader.i32(value.ccdIterations)
        || !reader.i32(value.intervalSubdivisions)) {
        return false;
    }
    value.kind = static_cast<ContactFeatureKind>(kind);
    value.ccdState = static_cast<CcdState>(ccdState);
    value.frictionState = static_cast<ContactFrictionState>(frictionState);
    for (std::size_t& node : value.firstNodes) {
        if (!reader.count(node, std::numeric_limits<std::size_t>::max())) {
            return false;
        }
    }
    for (std::size_t& node : value.secondNodes) {
        if (!reader.count(node, std::numeric_limits<std::size_t>::max())) {
            return false;
        }
    }
    if (!reader.u8(firstCount) || !reader.u8(secondCount)
        || firstCount > 3 || secondCount > 3) {
        return false;
    }
    value.firstNodeCount = firstCount;
    value.secondNodeCount = secondCount;
    for (double& weight : value.firstWeights) {
        if (!reader.finiteDouble(weight)) return false;
    }
    for (double& weight : value.secondWeights) {
        if (!reader.finiteDouble(weight)) return false;
    }
    if (!reader.vector3(value.firstPoint)
        || !reader.vector3(value.secondPoint)
        || !reader.vector3(value.normal)
        || !reader.finiteDouble(value.pairSeparation)
        || !reader.finiteDouble(value.gap)
        || !reader.finiteDouble(value.penetration)
        || !reader.finiteDouble(value.normalMultiplier)
        || !reader.finiteDouble(value.normalForceEstimate)
        || !reader.finiteDouble(value.normalImpulseMagnitude)
        || !reader.vector3(value.tangentialImpulse)
        || !reader.finiteDouble(value.frictionConeRatio)
        || !reader.finiteDouble(value.normalResidual)
        || !reader.finiteDouble(value.frictionResidual)
        || !reader.finiteDouble(value.frictionWork)
        || !reader.finiteDouble(value.tangentSpeedBefore)
        || !reader.finiteDouble(value.tangentSpeedAfter)
        || !reader.vector3(value.firstImpulse)
        || !reader.vector3(value.secondImpulse)
        || !reader.vector3(value.firstMoment)
        || !reader.vector3(value.secondMoment)) {
        return false;
    }
    for (Vec3& impulse : value.firstNodeImpulses) {
        if (!reader.vector3(impulse)) return false;
    }
    for (Vec3& impulse : value.secondNodeImpulses) {
        if (!reader.vector3(impulse)) return false;
    }
    if (!reader.i32(value.solverVisits)
        || value.kind != value.key.kind
        || value.ccdIterations < 0 || value.intervalSubdivisions < 0
        || value.solverVisits < 0 || !validKey(value.key, state)) {
        return false;
    }
    for (std::size_t index = 0; index < value.firstNodeCount; ++index) {
        if (value.firstNodes[index] >= state.nodes.size()) return false;
    }
    for (std::size_t index = 0; index < value.secondNodeCount; ++index) {
        if (value.secondNodes[index] >= state.nodes.size()) return false;
    }
    return true;
}

bool writeMutableState(Writer& writer,
                       const SoftBodyCheckpoint::State& state) {
    if (!writer.u32(stateVersion) || !writer.u32(0)
        || !writer.u64(state.topologyFingerprint)
        || !writer.count(state.nodes.size())
        || !writer.count(state.triangles.size())
        || !writer.count(state.constraints.size())
        || !writer.count(state.membranes.size())
        || !writer.count(state.dihedrals.size())
        || !writer.count(state.contactPairDiagnostics.size())) {
        return false;
    }
    for (const Node& value : state.nodes) {
        if (!writer.vector3(value.position)
            || !writer.vector3(value.previousPosition)
            || !writer.vector3(value.velocity)
            || !writer.vector3(value.force)) {
            return false;
        }
    }
    for (const Triangle& value : state.triangles) {
        if (!writer.finiteDouble(value.pressureDifference)) return false;
    }
    for (const DistanceConstraint& value : state.constraints) {
        if (!writer.finiteDouble(value.accumulatedLambda)) return false;
    }
    for (const MembraneElement& value : state.membranes) {
        if (!writer.vector3(value.multiplier)
            || !writer.vector3(value.solverResultantEstimate)
            || !writer.finiteDouble(value.normalizedResidual)) {
            return false;
        }
    }
    for (const DihedralBendingConstraint& value : state.dihedrals) {
        if (!writer.finiteDouble(value.accumulatedLambda)) return false;
    }
    if (!writer.count(state.contactMultipliers.size())) return false;
    for (const auto& [key, multiplier] : state.contactMultipliers) {
        if (!writeKey(writer, key) || !writer.finiteDouble(multiplier)) {
            return false;
        }
    }
    if (!writer.count(state.contactRecords.size())) return false;
    for (const ContactRecord& value : state.contactRecords) {
        if (!writeRecord(writer, value)) return false;
    }
    if (!writeDiagnostics(writer, state.contactDiagnostics)) return false;
    for (const ContactDiagnostics& value : state.contactPairDiagnostics) {
        if (!writeDiagnostics(writer, value)) return false;
    }
    const auto keys = [&](const auto& values) {
        if (!writer.count(values.size())) return false;
        for (const ContactFeatureKey& value : values) {
            if (!writeKey(writer, value)) return false;
        }
        return true;
    };
    return keys(state.iterationCandidateKeys)
        && keys(state.iterationQueryKeys)
        && keys(state.certificationCandidateKeys)
        && keys(state.certificationQueryKeys);
}

bool readExactCount(Reader& reader,
                    const std::size_t expected,
                    const std::size_t maximum,
                    bool& topologyMismatch) {
    std::size_t value = 0;
    if (!reader.count(value, maximum)) {
        return false;
    }
    if (value != expected) {
        topologyMismatch = true;
        return false;
    }
    return true;
}

bool readMutableState(
    Reader& reader,
    SoftBodyCheckpoint::State& state,
    const SoftBodyCheckpointPersistenceLimits& limits,
    bool& topologyMismatch) {
    std::uint32_t version = 0;
    std::uint32_t reserved = 0;
    std::uint64_t fingerprint = 0;
    if (!reader.u32(version) || !reader.u32(reserved)
        || !reader.u64(fingerprint)) {
        return false;
    }
    if (version != stateVersion || reserved != 0) {
        return false;
    }
    if (fingerprint != state.topologyFingerprint) {
        topologyMismatch = true;
        return false;
    }
    if (!readExactCount(reader, state.nodes.size(), limits.maximumNodes,
                        topologyMismatch)
        || !readExactCount(reader, state.triangles.size(),
                           limits.maximumTriangles, topologyMismatch)
        || !readExactCount(reader, state.constraints.size(),
                           limits.maximumConstraints, topologyMismatch)
        || !readExactCount(reader, state.membranes.size(),
                           limits.maximumMembranes, topologyMismatch)
        || !readExactCount(reader, state.dihedrals.size(),
                           limits.maximumDihedrals, topologyMismatch)
        || !readExactCount(reader, state.contactPairDiagnostics.size(),
                           limits.maximumContactPairs,
                           topologyMismatch)
        || !reader.fixedRecords(state.nodes.size(), nodeRecordBytes)) {
        return false;
    }
    for (Node& value : state.nodes) {
        if (!reader.vector3(value.position)
            || !reader.vector3(value.previousPosition)
            || !reader.vector3(value.velocity)
            || !reader.vector3(value.force)) {
            return false;
        }
    }
    if (!reader.fixedRecords(state.triangles.size(), triangleRecordBytes)) {
        return false;
    }
    for (Triangle& value : state.triangles) {
        if (!reader.finiteDouble(value.pressureDifference)) return false;
    }
    if (!reader.fixedRecords(state.constraints.size(), constraintRecordBytes)) {
        return false;
    }
    for (DistanceConstraint& value : state.constraints) {
        if (!reader.finiteDouble(value.accumulatedLambda)) return false;
    }
    if (!reader.fixedRecords(state.membranes.size(), membraneRecordBytes)) {
        return false;
    }
    for (MembraneElement& value : state.membranes) {
        if (!reader.vector3(value.multiplier)
            || !reader.vector3(value.solverResultantEstimate)
            || !reader.finiteDouble(value.normalizedResidual)) {
            return false;
        }
    }
    if (!reader.fixedRecords(state.dihedrals.size(), dihedralRecordBytes)) {
        return false;
    }
    for (DihedralBendingConstraint& value : state.dihedrals) {
        if (!reader.finiteDouble(value.accumulatedLambda)) return false;
    }

    std::size_t multiplierCount = 0;
    if (!reader.count(multiplierCount, limits.maximumContactMultipliers)
        || !reader.fixedRecords(multiplierCount,
                                contactMultiplierRecordBytes)) {
        return false;
    }
    state.contactMultipliers.clear();
    state.contactMultipliers.reserve(multiplierCount);
    for (std::size_t index = 0; index < multiplierCount; ++index) {
        ContactFeatureKey key;
        double multiplier = 0.0;
        if (!readKey(reader, key) || !reader.finiteDouble(multiplier)
            || !validKey(key, state)) {
            return false;
        }
        state.contactMultipliers.emplace_back(key, multiplier);
    }

    std::size_t recordCount = 0;
    if (!reader.count(recordCount, limits.maximumContactRecords)
        || !reader.fixedRecords(recordCount, contactRecordBytes)) {
        return false;
    }
    state.contactRecords.clear();
    state.contactRecords.resize(recordCount);
    for (ContactRecord& value : state.contactRecords) {
        if (!readRecord(reader, value, state)) return false;
    }
    if (!reader.fixedRecords(1 + state.contactPairDiagnostics.size(),
                             contactDiagnosticsBytes)
        || !readDiagnostics(reader, state.contactDiagnostics, limits, state)) {
        return false;
    }
    for (ContactDiagnostics& value : state.contactPairDiagnostics) {
        if (!readDiagnostics(reader, value, limits, state)) return false;
    }

    std::size_t auditTotal = 0;
    const auto keys = [&](auto& values) {
        std::size_t count = 0;
        if (!reader.count(count, limits.maximumAuditKeys)
            || count > limits.maximumAuditKeys - auditTotal
            || !reader.fixedRecords(count, contactKeyBytes)) {
            return false;
        }
        auditTotal += count;
        values.clear();
        values.resize(count);
        for (ContactFeatureKey& value : values) {
            if (!readKey(reader, value) || !validKey(value, state)) {
                return false;
            }
        }
        return true;
    };
    return keys(state.iterationCandidateKeys)
        && keys(state.iterationQueryKeys)
        && keys(state.certificationCandidateKeys)
        && keys(state.certificationQueryKeys);
}

} // namespace

bool serializeSoftBodyCheckpoint(
    const SoftBodyCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    SoftBodyCheckpointPersistenceError* error,
    const SoftBodyCheckpointPersistenceLimits& limits) {
    clearError(error);
    const SoftBodyCheckpoint::State* state =
        SoftBodyCheckpointPersistenceAccess::state(checkpoint);
    if (!validLimits(limits) || state == nullptr) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::InvalidData,
                    "invalid SoftBody checkpoint or persistence limits");
    }
    if (!withinStaticLimits(*state, limits)
        || !withinDynamicLimits(*state, limits)) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "SoftBody checkpoint exceeds persistence limits");
    }
    try {
        std::vector<std::uint8_t> payload;
        payload.reserve(std::min<std::size_t>(
            limits.maximumEncodedBytes - envelopeBytes, 4096));
        Writer payloadWriter(payload,
                             limits.maximumEncodedBytes - envelopeBytes);
        if (!writeMutableState(payloadWriter, *state)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? SoftBodyCheckpointPersistenceErrorCode::LimitExceeded
                    : SoftBodyCheckpointPersistenceErrorCode::InvalidData,
                payloadWriter.exceeded()
                    ? "SoftBody checkpoint exceeds the encoded byte limit"
                    : "SoftBody checkpoint contains invalid mutable state");
        }

        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) {
                return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "SoftBody checkpoint envelope exceeds byte limit");
            }
        }
        if (!writer.u16(softBodyCheckpointProtocolVersion)
            || !writer.u16(0)
            || !writer.count(payload.size())
            || !writer.u64(checksum(payload))
            || !writer.bytes(payload)) {
            return fail(error,
                SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                "SoftBody checkpoint exceeds the encoded byte limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate SoftBody checkpoint encoding");
    } catch (const std::length_error&) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "SoftBody checkpoint encoding is too large");
    }
}

bool deserializeSoftBodyCheckpoint(
    const std::span<const std::uint8_t> bytes,
    const SoftBodyCheckpoint& topologyTemplate,
    SoftBodyCheckpoint& checkpoint,
    SoftBodyCheckpointPersistenceError* error,
    const SoftBodyCheckpointPersistenceLimits& limits) {
    clearError(error);
    const SoftBodyCheckpoint::State* templateState =
        SoftBodyCheckpointPersistenceAccess::state(topologyTemplate);
    if (!validLimits(limits) || templateState == nullptr) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::InvalidData,
                    "invalid topology template or persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "SoftBody checkpoint exceeds the encoded byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::Truncated,
                    "SoftBody checkpoint envelope is truncated");
    }
    if (!withinStaticLimits(*templateState, limits)) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "topology template exceeds persistence limits");
    }

    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) {
            return fail(error,
                        SoftBodyCheckpointPersistenceErrorCode::Truncated,
                        "SoftBody checkpoint envelope is truncated");
        }
        if (actual != expected) {
            return fail(error,
                        SoftBodyCheckpointPersistenceErrorCode::InvalidMagic,
                        "SoftBody checkpoint magic is invalid");
        }
    }
    std::uint16_t version = 0;
    std::uint16_t reserved = 0;
    std::size_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(version) || !envelope.u16(reserved)
        || !envelope.count(payloadSize,
                           limits.maximumEncodedBytes - envelopeBytes)
        || !envelope.u64(expectedChecksum)) {
        return fail(error, readerErrorCode(envelope),
                    "SoftBody checkpoint envelope is invalid");
    }
    if (version != softBodyCheckpointProtocolVersion) {
        return fail(error,
            SoftBodyCheckpointPersistenceErrorCode::UnsupportedVersion,
            "SoftBody checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::InvalidData,
                    "SoftBody checkpoint reserved bits are nonzero");
    }
    const std::span<const std::uint8_t> payload =
        envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::Truncated,
                    "SoftBody checkpoint payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::TrailingData,
                    "SoftBody checkpoint has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::ChecksumMismatch,
                    "SoftBody checkpoint checksum does not match");
    }

    try {
        auto candidate =
            std::make_shared<SoftBodyCheckpoint::State>(*templateState);
        Reader reader(payload);
        bool topologyMismatch = false;
        if (!readMutableState(reader, *candidate, limits,
                              topologyMismatch)) {
            return fail(
                error,
                topologyMismatch
                    ? SoftBodyCheckpointPersistenceErrorCode::TopologyMismatch
                    : readerErrorCode(reader),
                topologyMismatch
                    ? "SoftBody checkpoint belongs to a different topology"
                    : "SoftBody checkpoint mutable state is invalid");
        }
        if (!reader.atEnd()) {
            return fail(error,
                        SoftBodyCheckpointPersistenceErrorCode::TrailingData,
                        "SoftBody checkpoint payload has trailing data");
        }
        checkpoint = SoftBodyCheckpointPersistenceAccess::make(
            std::move(candidate));
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate decoded SoftBody checkpoint");
    } catch (const std::length_error&) {
        return fail(error,
                    SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
                    "decoded SoftBody checkpoint is too large");
    }
}

} // namespace softwing
