#include "scene_fluid_mimetic_pressure_state_persistence.h"

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

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'M', 'P'};
constexpr std::uint32_t payloadVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t controlRecordBytes = 56;
constexpr std::size_t traceRecordBytes = 48;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

using ErrorCode = SceneFluidMimeticPressureStatePersistenceErrorCode;
using Error = SceneFluidMimeticPressureStatePersistenceError;
using Limits = SceneFluidMimeticPressureStatePersistenceLimits;

bool fail(Error* error, const ErrorCode code, std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(Error* error) {
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
        return std::isfinite(value)
            && u64(std::bit_cast<std::uint64_t>(value));
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
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded) || encoded > 1) return false;
        value = encoded != 0;
        return true;
    }
    bool fixedRecords(const std::size_t count,
                      const std::size_t recordBytes) {
        if (count > (bytes_.size() - position_) / recordBytes) {
            truncated_ = true;
            return false;
        }
        return true;
    }
    std::span<const std::uint8_t> remainingBytes(const std::size_t count) {
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

ErrorCode readerErrorCode(const Reader& reader) noexcept {
    if (reader.limitExceeded()) return ErrorCode::LimitExceeded;
    if (reader.truncated()) return ErrorCode::Truncated;
    return ErrorCode::InvalidData;
}

bool validLimits(const Limits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.maximumControlCells > 0
        && limits.maximumReducedTraces > 0;
}

bool writePayload(Writer& writer,
                  const SceneFluidMimeticPressureState& state) {
    if (!writer.u32(payloadVersion) || !writer.u32(0)
        || !writer.u32(state.version) || !writer.u32(0)
        || !writer.u64(state.fingerprint)
        || !writer.u64(state.mimeticControlCellFingerprint)
        || !writer.u64(state.fullTraceSystemFingerprint)
        || !writer.u64(state.condensedTraceSystemFingerprint)
        || !writer.u64(state.pressureSourceFingerprint)
        || !writer.u64(state.structureDefinitionFingerprint)
        || !writer.u64(state.acceptedStepCount)
        || !writer.finiteDouble(state.simulationTimeSeconds)
        || !writer.count(state.ownedStorageBytes)
        || !writer.count(state.componentCount)
        || !writer.finiteDouble(
            state.maximumAbsoluteControlPressurePascals)
        || !writer.finiteDouble(
            state.maximumAbsoluteTracePressurePascals)
        || !writer.count(state.controls.size())) {
        return false;
    }
    for (const auto& control : state.controls) {
        if (!writer.count(control.controlCellIndex)
            || !writer.count(control.controlVolumeIndex)
            || !writer.u64(control.stableId)
            || !writer.count(control.cellIndex)
            || !writer.u64(control.regionId)
            || !writer.count(control.componentIndex)
            || !writer.finiteDouble(control.pressurePascals)) {
            return false;
        }
    }
    if (!writer.count(state.traces.size())) return false;
    for (const auto& trace : state.traces) {
        if (!writer.count(trace.reducedTraceIndex)
            || !writer.count(trace.fullTraceIndex)
            || !writer.u64(trace.stableId)
            || !writer.u8(static_cast<std::uint8_t>(trace.kind))
            || !writer.u8(trace.isGauge ? 1 : 0)
            || !writer.u16(0) || !writer.u32(0)
            || !writer.count(trace.componentIndex)
            || !writer.finiteDouble(trace.pressurePascals)) {
            return false;
        }
    }
    return true;
}

bool topologyMatches(
    const SceneFluidMimeticPressureState& state,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem) {
    return state.mimeticControlCellFingerprint == controlCells.fingerprint
        && state.fullTraceSystemFingerprint == fullSystem.fingerprint
        && state.condensedTraceSystemFingerprint
            == condensedSystem.fingerprint
        && state.structureDefinitionFingerprint
            == controlCells.structureDefinitionFingerprint
        && state.acceptedStepCount == controlCells.acceptedStepCount
        && state.simulationTimeSeconds
            == controlCells.simulationTimeSeconds
        && state.componentCount == condensedSystem.componentCount
        && state.controls.size() == controlCells.controlCells.size()
        && state.traces.size() == condensedSystem.traces.size();
}

} // namespace

bool serializeSceneFluidMimeticPressureState(
    const SceneFluidMimeticPressureState& state,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    std::vector<std::uint8_t>& bytes,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "invalid mimetic pressure-state persistence limits");
    }
    try {
        validateSceneFluidMimeticPressureState(
            state, controlCells, fullSystem, condensedSystem);
        if (state.controls.size() > limits.maximumControlCells
            || state.traces.size() > limits.maximumReducedTraces) {
            return fail(error, ErrorCode::LimitExceeded,
                        "mimetic pressure state exceeds record limits");
        }
        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload, limits.maximumEncodedBytes - envelopeBytes);
        if (!writePayload(payloadWriter, state)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? ErrorCode::LimitExceeded : ErrorCode::InvalidData,
                "mimetic pressure-state payload is invalid or too large");
        }
        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) {
                return fail(error, ErrorCode::LimitExceeded,
                            "mimetic pressure-state envelope is too large");
            }
        }
        if (!writer.u16(sceneFluidMimeticPressureStateProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "mimetic pressure state exceeds its encoded limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate mimetic pressure-state encoding");
    } catch (const std::length_error&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "mimetic pressure-state encoding is too large");
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData,
                    std::string("mimetic pressure-state validation failed: ")
                        + exception.what());
    }
}

bool deserializeSceneFluidMimeticPressureState(
    const std::span<const std::uint8_t> bytes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    SceneFluidMimeticPressureState& state,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "invalid mimetic pressure-state persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error, ErrorCode::LimitExceeded,
                    "mimetic pressure state exceeds its encoded limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error, ErrorCode::Truncated,
                    "mimetic pressure-state envelope is truncated");
    }
    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual) || actual != expected) {
            return fail(error, ErrorCode::InvalidMagic,
                        "mimetic pressure-state magic is invalid");
        }
    }
    std::uint16_t protocol = 0;
    std::uint16_t reserved = 0;
    std::size_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(protocol) || !envelope.u16(reserved)
        || !envelope.count(
            payloadSize, limits.maximumEncodedBytes - envelopeBytes)
        || !envelope.u64(expectedChecksum)) {
        return fail(error, readerErrorCode(envelope),
                    "mimetic pressure-state envelope is invalid");
    }
    if (protocol != sceneFluidMimeticPressureStateProtocolVersion) {
        return fail(error, ErrorCode::UnsupportedVersion,
                    "mimetic pressure-state protocol is unsupported");
    }
    if (reserved != 0) {
        return fail(error, ErrorCode::InvalidData,
                    "mimetic pressure-state reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error, ErrorCode::Truncated,
                    "mimetic pressure-state payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error, ErrorCode::TrailingData,
                    "mimetic pressure state has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error, ErrorCode::ChecksumMismatch,
                    "mimetic pressure-state checksum does not match");
    }

    try {
        validateSceneFluidMimeticTraceSystem(fullSystem, controlCells);
        validateSceneFluidMimeticCondensedTraceSystem(
            condensedSystem, fullSystem);
        Reader reader(payload);
        SceneFluidMimeticPressureState candidate;
        std::uint32_t decodedPayloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t stateReserved = 0;
        if (!reader.u32(decodedPayloadVersion)
            || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(stateReserved)
            || !reader.u64(candidate.fingerprint)
            || !reader.u64(candidate.mimeticControlCellFingerprint)
            || !reader.u64(candidate.fullTraceSystemFingerprint)
            || !reader.u64(candidate.condensedTraceSystemFingerprint)
            || !reader.u64(candidate.pressureSourceFingerprint)
            || !reader.u64(candidate.structureDefinitionFingerprint)
            || !reader.u64(candidate.acceptedStepCount)
            || !reader.finiteDouble(candidate.simulationTimeSeconds)
            || !reader.count(
                candidate.ownedStorageBytes,
                std::numeric_limits<std::size_t>::max())
            || !reader.count(
                candidate.componentCount,
                std::numeric_limits<std::size_t>::max())
            || !reader.finiteDouble(
                candidate.maximumAbsoluteControlPressurePascals)
            || !reader.finiteDouble(
                candidate.maximumAbsoluteTracePressurePascals)) {
            return fail(error, readerErrorCode(reader),
                        "mimetic pressure-state payload header is invalid");
        }
        if (decodedPayloadVersion != payloadVersion
            || candidate.version
                != sceneFluidMimeticPressureStateVersion) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "mimetic pressure-state payload version is unsupported");
        }
        if (payloadReserved != 0 || stateReserved != 0) {
            return fail(error, ErrorCode::InvalidData,
                        "mimetic pressure-state payload reserved bits are nonzero");
        }
        std::size_t controlCount = 0;
        if (!reader.count(controlCount, limits.maximumControlCells)) {
            return fail(error, readerErrorCode(reader),
                        "mimetic pressure-state control count is invalid");
        }
        if (controlCount != controlCells.controlCells.size()) {
            return fail(error, ErrorCode::TopologyMismatch,
                        "mimetic pressure-state control count does not match topology");
        }
        if (!reader.fixedRecords(controlCount, controlRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "mimetic pressure-state controls are truncated");
        }
        candidate.controls.resize(controlCount);
        for (auto& control : candidate.controls) {
            if (!reader.count(
                    control.controlCellIndex,
                    limits.maximumControlCells)
                || !reader.count(
                    control.controlVolumeIndex,
                    limits.maximumControlCells)
                || !reader.u64(control.stableId)
                || !reader.count(
                    control.cellIndex,
                    std::numeric_limits<std::size_t>::max())
                || !reader.u64(control.regionId)
                || !reader.count(
                    control.componentIndex,
                    std::numeric_limits<std::size_t>::max())
                || !reader.finiteDouble(control.pressurePascals)) {
                return fail(error, readerErrorCode(reader),
                            "mimetic pressure-state control is invalid");
            }
        }
        std::size_t traceCount = 0;
        if (!reader.count(traceCount, limits.maximumReducedTraces)) {
            return fail(error, readerErrorCode(reader),
                        "mimetic pressure-state trace count is invalid");
        }
        if (traceCount != condensedSystem.traces.size()) {
            return fail(error, ErrorCode::TopologyMismatch,
                        "mimetic pressure-state trace count does not match topology");
        }
        if (!reader.fixedRecords(traceCount, traceRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "mimetic pressure-state traces are truncated");
        }
        candidate.traces.resize(traceCount);
        for (auto& trace : candidate.traces) {
            std::uint8_t kind = 0;
            std::uint16_t rowReserved16 = 0;
            std::uint32_t rowReserved32 = 0;
            if (!reader.count(
                    trace.reducedTraceIndex,
                    limits.maximumReducedTraces)
                || !reader.count(
                    trace.fullTraceIndex,
                    std::numeric_limits<std::size_t>::max())
                || !reader.u64(trace.stableId) || !reader.u8(kind)
                || !reader.boolean(trace.isGauge)
                || !reader.u16(rowReserved16)
                || !reader.u32(rowReserved32)
                || !reader.count(
                    trace.componentIndex,
                    std::numeric_limits<std::size_t>::max())
                || !reader.finiteDouble(trace.pressurePascals)) {
                return fail(error, readerErrorCode(reader),
                            "mimetic pressure-state trace is invalid");
            }
            if (rowReserved16 != 0 || rowReserved32 != 0) {
                return fail(error, ErrorCode::InvalidData,
                            "mimetic pressure-state trace reserved bits are nonzero");
            }
            trace.kind = static_cast<SceneFluidMimeticHalfFaceKind>(kind);
        }
        if (!reader.atEnd()) {
            return fail(error, ErrorCode::TrailingData,
                        "mimetic pressure-state payload has trailing data");
        }
        if (!topologyMatches(
                candidate, controlCells, fullSystem, condensedSystem)) {
            return fail(error, ErrorCode::TopologyMismatch,
                        "mimetic pressure state belongs to a different topology");
        }
        validateSceneFluidMimeticPressureState(
            candidate, controlCells, fullSystem, condensedSystem);
        state = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate decoded mimetic pressure state");
    } catch (const std::length_error&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "decoded mimetic pressure state is too large");
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData,
                    std::string("decoded mimetic pressure state is invalid: ")
                        + exception.what());
    }
}

} // namespace simwing::fsi
