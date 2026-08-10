#include "scene_fluid_regional_opening_momentum_wall_cycle_state_persistence.h"

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

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'R', 'W'};
constexpr std::uint32_t payloadVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t adjustmentControlRecordBytes = 88;
constexpr std::size_t wallTractionRecordBytes = 32;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

using ErrorCode =
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode;
using Error =
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError;
using Limits =
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits;

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
    bool bytes(const std::span<const std::uint8_t> values) {
        return raw(values.data(), values.size());
    }
    bool zeroes(const std::size_t count) {
        for (std::size_t index = 0; index < count; ++index) {
            if (!u8(0)) return false;
        }
        return true;
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
    bool zeroes(const std::size_t count) {
        for (std::size_t index = 0; index < count; ++index) {
            std::uint8_t value = 0;
            if (!u8(value) || value != 0) return false;
        }
        return true;
    }
    bool fixedRecords(const std::size_t count,
                      const std::size_t recordBytes) {
        if (recordBytes == 0
            || count > (bytes_.size() - position_) / recordBytes) {
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
        && limits.maximumAdjustmentControls > 0
        && limits.maximumWallTractions > 0
        && limits.state.maximumWallTractions > 0
        && limits.state.maximumOwnedBytes > 0
        && limits.state.adjustmentLimits.maximumFragments > 0
        && limits.state.adjustmentLimits.maximumOwnedBytes > 0
        && limits.state.acceptedStateLimits.maximumTopologyLinkVelocities > 0
        && limits.state.acceptedStateLimits.maximumOpeningSamples > 0
        && limits.state.acceptedStateLimits.maximumPressureCorrections > 0
        && limits.state.acceptedStateLimits.maximumOwnedBytes > 0
        && limits.state.acceptedStateLimits.maximumWorkingBytes > 0;
}

fluid::PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits
acceptedPersistenceLimits(const Limits& limits) {
    fluid::PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits
        result;
    result.stateLimits = limits.state.acceptedStateLimits;
    result.maximumEncodedBytes = limits.maximumEncodedBytes;
    result.maximumTopologyLinkVelocities =
        limits.state.acceptedStateLimits.maximumTopologyLinkVelocities;
    result.maximumOpeningSamples =
        limits.state.acceptedStateLimits.maximumOpeningSamples;
    result.maximumPressureCorrections =
        limits.state.acceptedStateLimits.maximumPressureCorrections;
    return result;
}

template<typename Vector>
bool writeVector(Writer& writer, const Vector& value) {
    return writer.finiteDouble(value.x)
        && writer.finiteDouble(value.y)
        && writer.finiteDouble(value.z);
}

template<typename Vector>
bool readVector(Reader& reader, Vector& value) {
    return reader.finiteDouble(value.x)
        && reader.finiteDouble(value.y)
        && reader.finiteDouble(value.z);
}

bool writeAdjustmentState(
    Writer& writer,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState&
        state) {
    if (!writer.u32(state.version) || !writer.u32(0)
        || !writer.u64(state.fingerprint)
        || !writer.u64(state.sourceTransportFingerprint)
        || !writer.u64(state.sourceAdjustmentFingerprint)
        || !writer.u64(state.sourceMetricFingerprint)
        || !writer.finiteDouble(state.densityKgPerCubicMeter)
        || !writer.finiteDouble(state.timeStepSeconds)
        || !writer.finiteDouble(
            state.settings
                .absoluteMomentumToleranceKilogramMetersPerSecond)
        || !writer.finiteDouble(state.settings.relativeMomentumTolerance)
        || !writeVector(
            writer, state.sourceMomentumKilogramMetersPerSecond)
        || !writeVector(
            writer, state.adjustedMomentumKilogramMetersPerSecond)
        || !writeVector(
            writer, state.adjustmentImpulseKilogramMetersPerSecond)
        || !writer.finiteDouble(state.sourceKineticEnergyJoules)
        || !writer.finiteDouble(state.adjustedKineticEnergyJoules)
        || !writer.finiteDouble(state.kineticEnergyChangeJoules)
        || !writer.finiteDouble(
            state.maximumMomentumVelocityResidualKilogramMetersPerSecond)
        || !writer.count(state.ownedStorageBytes)
        || !writer.count(state.controls.size())) {
        return false;
    }
    for (const auto& control : state.controls) {
        if (!writer.count(control.fragmentIndex)
            || !writer.u64(control.stableId)
            || !writer.u64(control.regionStableId)
            || !writer.count(control.connectedComponentIndex)
            || !writer.finiteDouble(control.volumeCubicMeters)
            || !writeVector(writer, control.velocityMetersPerSecond)
            || !writeVector(
                writer, control.momentumKilogramMetersPerSecond)) {
            return false;
        }
    }
    return true;
}

bool readAdjustmentState(
    Reader& reader,
    fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState& state,
    const Limits& limits) {
    std::uint32_t reserved = 0;
    std::size_t controlCount = 0;
    if (!reader.u32(state.version) || !reader.u32(reserved)
        || reserved != 0 || !reader.u64(state.fingerprint)
        || !reader.u64(state.sourceTransportFingerprint)
        || !reader.u64(state.sourceAdjustmentFingerprint)
        || !reader.u64(state.sourceMetricFingerprint)
        || !reader.finiteDouble(state.densityKgPerCubicMeter)
        || !reader.finiteDouble(state.timeStepSeconds)
        || !reader.finiteDouble(
            state.settings
                .absoluteMomentumToleranceKilogramMetersPerSecond)
        || !reader.finiteDouble(state.settings.relativeMomentumTolerance)
        || !readVector(
            reader, state.sourceMomentumKilogramMetersPerSecond)
        || !readVector(
            reader, state.adjustedMomentumKilogramMetersPerSecond)
        || !readVector(
            reader, state.adjustmentImpulseKilogramMetersPerSecond)
        || !reader.finiteDouble(state.sourceKineticEnergyJoules)
        || !reader.finiteDouble(state.adjustedKineticEnergyJoules)
        || !reader.finiteDouble(state.kineticEnergyChangeJoules)
        || !reader.finiteDouble(
            state.maximumMomentumVelocityResidualKilogramMetersPerSecond)
        || !reader.count(
            state.ownedStorageBytes,
            limits.state.adjustmentLimits.maximumOwnedBytes)
        || !reader.count(
            controlCount, limits.maximumAdjustmentControls)) {
        return false;
    }
    if (!reader.fixedRecords(
            controlCount, adjustmentControlRecordBytes)) {
        return false;
    }
    state.controls.resize(controlCount);
    for (auto& control : state.controls) {
        if (!reader.count(
                control.fragmentIndex,
                limits.maximumAdjustmentControls)
            || !reader.u64(control.stableId)
            || !reader.u64(control.regionStableId)
            || !reader.count(
                control.connectedComponentIndex,
                limits.state.adjustmentLimits.maximumFragments)
            || !reader.finiteDouble(control.volumeCubicMeters)
            || !readVector(reader, control.velocityMetersPerSecond)
            || !readVector(
                reader, control.momentumKilogramMetersPerSecond)) {
            return false;
        }
    }
    return true;
}

bool writeWallTractions(
    Writer& writer,
    const SceneFluidAcceptedWallTractionSet& tractions) {
    if (!writer.u32(tractions.version) || !writer.u32(0)
        || !writer.u64(tractions.fingerprint)
        || !writer.u64(tractions.wallExchangeFingerprint)
        || !writer.u64(tractions.quadratureFingerprint)
        || !writer.u64(tractions.surfaceDefinitionFingerprint)
        || !writer.u64(tractions.surfaceStateFingerprint)
        || !writer.u64(tractions.structureDefinitionFingerprint)
        || !writer.u64(tractions.acceptedStepCount)
        || !writer.finiteDouble(tractions.simulationTimeSeconds)
        || !writer.count(tractions.ownedStorageBytes)
        || !writer.count(tractions.tractions.size())) {
        return false;
    }
    for (const auto& traction : tractions.tractions) {
        if (!writer.u64(traction.stableId)
            || !writeVector(writer, traction.tractionPascals)) {
            return false;
        }
    }
    return true;
}

bool readWallTractions(
    Reader& reader,
    SceneFluidAcceptedWallTractionSet& tractions,
    const Limits& limits) {
    std::uint32_t reserved = 0;
    std::size_t tractionCount = 0;
    if (!reader.u32(tractions.version) || !reader.u32(reserved)
        || reserved != 0 || !reader.u64(tractions.fingerprint)
        || !reader.u64(tractions.wallExchangeFingerprint)
        || !reader.u64(tractions.quadratureFingerprint)
        || !reader.u64(tractions.surfaceDefinitionFingerprint)
        || !reader.u64(tractions.surfaceStateFingerprint)
        || !reader.u64(tractions.structureDefinitionFingerprint)
        || !reader.u64(tractions.acceptedStepCount)
        || !reader.finiteDouble(tractions.simulationTimeSeconds)
        || !reader.count(
            tractions.ownedStorageBytes, limits.state.maximumOwnedBytes)
        || !reader.count(
            tractionCount, limits.maximumWallTractions)) {
        return false;
    }
    if (!reader.fixedRecords(tractionCount, wallTractionRecordBytes)) {
        return false;
    }
    tractions.tractions.resize(tractionCount);
    for (auto& traction : tractions.tractions) {
        if (!reader.u64(traction.stableId)
            || !readVector(reader, traction.tractionPascals)) {
            return false;
        }
    }
    return true;
}

ErrorCode nestedErrorCode(
    const fluid::
        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode
            code) {
    using Nested = fluid::
        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode;
    switch (code) {
    case Nested::None:
        return ErrorCode::None;
    case Nested::LimitExceeded:
        return ErrorCode::LimitExceeded;
    case Nested::SourceMismatch:
        return ErrorCode::SourceMismatch;
    case Nested::InvalidData:
    case Nested::InvalidMagic:
    case Nested::UnsupportedVersion:
    case Nested::Truncated:
    case Nested::TrailingData:
    case Nested::ChecksumMismatch:
        return ErrorCode::InvalidData;
    }
    return ErrorCode::InvalidData;
}

bool hasForeignLiveSources(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature) {
    return state.transportMetricFingerprint != transportMetric.fingerprint
        || state.adjustedMomentum.sourceMetricFingerprint
            != transportMetric.fingerprint
        || state.acceptedMetricFingerprint != acceptedMetric.fingerprint
        || state.wallTractions.quadratureFingerprint
            != quadrature.fingerprint
        || state.wallTractions.surfaceDefinitionFingerprint
            != quadrature.surfaceDefinitionFingerprint
        || state.wallTractions.surfaceStateFingerprint
            != quadrature.surfaceStateFingerprint
        || state.wallTractions.structureDefinitionFingerprint
            != quadrature.structureDefinitionFingerprint
        || state.wallTractions.acceptedStepCount
            != quadrature.acceptedStepCount
        || state.wallTractions.simulationTimeSeconds
            != quadrature.simulationTimeSeconds;
}

} // namespace

bool serializeSceneFluidRegionalOpeningMomentumWallCycleState(
    const SceneFluidRegionalOpeningMomentumWallCycleState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    std::vector<std::uint8_t>& bytes,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "regional opening wall-cycle persistence limits are invalid");
    }
    try {
        validateSceneFluidRegionalOpeningMomentumWallCycleState(
            state, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, quadrature, limits.state);
        if (state.adjustedMomentum.controls.size()
                > limits.maximumAdjustmentControls
            || state.wallTractions.tractions.size()
                > limits.maximumWallTractions) {
            return fail(error, ErrorCode::LimitExceeded,
                        "regional opening wall-cycle record limit exceeded");
        }

        std::vector<std::uint8_t> acceptedBytes;
        fluid::
            PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
                acceptedError;
        if (!fluid::
                serializePlanarPressureRegionFragmentOpeningAcceptedState(
                    state.acceptedPressure, acceptedPressureOperator,
                    acceptedBasePressureOperator, grid, acceptedSweep,
                    acceptedFragments, acceptedTopology,
                    acceptedVolumeRates, acceptedOpeningDefinitions,
                    acceptedOpenings, acceptedResistanceDefinitions,
                    acceptedBytes, &acceptedError,
                    acceptedPersistenceLimits(limits))) {
            return fail(
                error, nestedErrorCode(acceptedError.code),
                std::string("cannot encode nested opening accepted state: ")
                    + acceptedError.message);
        }

        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload, limits.maximumEncodedBytes - envelopeBytes);
        if (!payloadWriter.u32(payloadVersion)
            || !payloadWriter.u32(0)
            || !payloadWriter.u32(state.version)
            || !payloadWriter.u32(0)
            || !payloadWriter.u64(state.fingerprint)
            || !payloadWriter.u64(
                state.sourceWallPressureEpochFingerprint)
            || !payloadWriter.u64(state.sourceWallExchangeFingerprint)
            || !payloadWriter.u64(state.transportMetricFingerprint)
            || !payloadWriter.u64(state.acceptedMetricFingerprint)
            || !payloadWriter.u64(state.predictionFingerprint)
            || !payloadWriter.u64(state.pressureWarmStartFingerprint)
            || !payloadWriter.u64(state.predictedOpeningFluxFingerprint)
            || !payloadWriter.count(state.ownedStorageBytes)
            || !writeAdjustmentState(
                payloadWriter, state.adjustedMomentum)
            || !writeWallTractions(payloadWriter, state.wallTractions)
            || !payloadWriter.count(acceptedBytes.size())
            || !payloadWriter.bytes(acceptedBytes)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? ErrorCode::LimitExceeded : ErrorCode::InvalidData,
                "regional opening wall-cycle payload is invalid or too large");
        }

        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) {
                return fail(error, ErrorCode::LimitExceeded,
                            "regional opening wall-cycle envelope is too large");
            }
        }
        if (!writer.u16(
                sceneFluidRegionalOpeningMomentumWallCycleStateProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "regional opening wall-cycle envelope is too large");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate encoded regional opening wall-cycle state");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(
            error, ErrorCode::InvalidData,
            std::string("regional opening wall-cycle state is invalid: ")
                + exception.what());
    }
}

bool deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
    const std::span<const std::uint8_t> bytes,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidQuadratureDefinition& quadrature,
    SceneFluidRegionalOpeningMomentumWallCycleState& state,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "regional opening wall-cycle persistence limits are invalid");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error, ErrorCode::LimitExceeded,
                    "regional opening wall-cycle envelope exceeds byte limit");
    }

    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) {
            return fail(error, ErrorCode::Truncated,
                        "regional opening wall-cycle envelope is truncated");
        }
        if (actual != expected) {
            return fail(error, ErrorCode::InvalidMagic,
                        "regional opening wall-cycle magic is invalid");
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
                    "regional opening wall-cycle envelope header is invalid");
    }
    if (protocol
        != sceneFluidRegionalOpeningMomentumWallCycleStateProtocolVersion) {
        return fail(error, ErrorCode::UnsupportedVersion,
                    "regional opening wall-cycle protocol is unsupported");
    }
    if (reserved != 0) {
        return fail(error, ErrorCode::InvalidData,
                    "regional opening wall-cycle envelope reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error, ErrorCode::Truncated,
                    "regional opening wall-cycle payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error, ErrorCode::TrailingData,
                    "regional opening wall-cycle envelope has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error, ErrorCode::ChecksumMismatch,
                    "regional opening wall-cycle payload checksum differs");
    }

    try {
        Reader reader(payload);
        SceneFluidRegionalOpeningMomentumWallCycleState candidate;
        std::uint32_t decodedPayloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t stateReserved = 0;
        if (!reader.u32(decodedPayloadVersion)
            || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(stateReserved)) {
            return fail(error, readerErrorCode(reader),
                        "regional opening wall-cycle payload header is invalid");
        }
        if (decodedPayloadVersion != payloadVersion) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "regional opening wall-cycle payload version is unsupported");
        }
        if (payloadReserved != 0 || stateReserved != 0) {
            return fail(error, ErrorCode::InvalidData,
                        "regional opening wall-cycle payload reserved bits are nonzero");
        }
        if (!reader.u64(candidate.fingerprint)
            || !reader.u64(
                candidate.sourceWallPressureEpochFingerprint)
            || !reader.u64(candidate.sourceWallExchangeFingerprint)
            || !reader.u64(candidate.transportMetricFingerprint)
            || !reader.u64(candidate.acceptedMetricFingerprint)
            || !reader.u64(candidate.predictionFingerprint)
            || !reader.u64(candidate.pressureWarmStartFingerprint)
            || !reader.u64(candidate.predictedOpeningFluxFingerprint)
            || !reader.count(
                candidate.ownedStorageBytes,
                limits.state.maximumOwnedBytes)
            || !readAdjustmentState(
                reader, candidate.adjustedMomentum, limits)
            || !readWallTractions(
                reader, candidate.wallTractions, limits)) {
            return fail(error, readerErrorCode(reader),
                        "regional opening wall-cycle payload is invalid");
        }
        if (hasForeignLiveSources(
                candidate, transportMetric, acceptedMetric, quadrature)) {
            return fail(error, ErrorCode::SourceMismatch,
                        "regional opening wall-cycle state belongs to different restart sources");
        }

        std::size_t acceptedByteCount = 0;
        if (!reader.count(
                acceptedByteCount,
                limits.maximumEncodedBytes - envelopeBytes)) {
            return fail(error, readerErrorCode(reader),
                        "regional opening wall-cycle nested state size is invalid");
        }
        const auto acceptedBytes = reader.remainingBytes(acceptedByteCount);
        if (reader.truncated()) {
            return fail(error, ErrorCode::Truncated,
                        "regional opening wall-cycle nested state is truncated");
        }
        if (!reader.atEnd()) {
            return fail(error, ErrorCode::TrailingData,
                        "regional opening wall-cycle payload has trailing data");
        }

        fluid::
            PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
                acceptedError;
        if (!fluid::
                deserializePlanarPressureRegionFragmentOpeningAcceptedState(
                    acceptedBytes, acceptedPressureOperator,
                    acceptedBasePressureOperator, grid, acceptedSweep,
                    acceptedFragments, acceptedTopology,
                    acceptedVolumeRates, acceptedOpeningDefinitions,
                    acceptedOpenings, acceptedResistanceDefinitions,
                    candidate.acceptedPressure, &acceptedError,
                    acceptedPersistenceLimits(limits))) {
            return fail(
                error, nestedErrorCode(acceptedError.code),
                std::string("cannot decode nested opening accepted state: ")
                    + acceptedError.message);
        }
        validateSceneFluidRegionalOpeningMomentumWallCycleState(
            candidate, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, quadrature, limits.state);
        state = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate decoded regional opening wall-cycle state");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(
            error, ErrorCode::InvalidData,
            std::string("decoded regional opening wall-cycle state is invalid: ")
                + exception.what());
    }
}

} // namespace simwing::fsi
