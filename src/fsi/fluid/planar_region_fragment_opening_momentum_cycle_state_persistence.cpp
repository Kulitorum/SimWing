#include "fluid/planar_region_fragment_opening_momentum_cycle_state_persistence.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'R', 'M'};
constexpr std::uint32_t payloadVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t transportControlRecordBytes = 88;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

using ErrorCode =
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode;
using Error =
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError;
using Limits =
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits;

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
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded) || encoded > 1) return false;
        value = encoded != 0;
        return true;
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
        && limits.maximumTransportControls > 0
        && limits.state.maximumOwnedBytes > 0
        && limits.state.transport.maximumFragments > 0
        && limits.state.transport.maximumDofs > 0
        && limits.state.transport.maximumOwnedBytes > 0
        && limits.state.transport.maximumWorkingBytes > 0
        && limits.state.acceptedState.maximumTopologyLinkVelocities > 0
        && limits.state.acceptedState.maximumOpeningSamples > 0
        && limits.state.acceptedState.maximumPressureCorrections > 0
        && limits.state.acceptedState.maximumOwnedBytes > 0
        && limits.state.acceptedState.maximumWorkingBytes > 0;
}

PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits
acceptedPersistenceLimits(const Limits& limits) {
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits result;
    result.stateLimits = limits.state.acceptedState;
    result.maximumEncodedBytes = limits.maximumEncodedBytes;
    result.maximumTopologyLinkVelocities =
        limits.state.acceptedState.maximumTopologyLinkVelocities;
    result.maximumOpeningSamples =
        limits.state.acceptedState.maximumOpeningSamples;
    result.maximumPressureCorrections =
        limits.state.acceptedState.maximumPressureCorrections;
    return result;
}

bool writeVector(Writer& writer, const Vector3& value) {
    return writer.finiteDouble(value.x)
        && writer.finiteDouble(value.y)
        && writer.finiteDouble(value.z);
}

bool readVector(Reader& reader, Vector3& value) {
    return reader.finiteDouble(value.x)
        && reader.finiteDouble(value.y)
        && reader.finiteDouble(value.z);
}

bool writeTransportSettings(
    Writer& writer,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings) {
    return writer.finiteDouble(settings.maximumOutgoingCourantNumber)
        && writer.count(settings.maximumSubsteps)
        && writer.finiteDouble(
            settings.absoluteContinuityToleranceCubicMetersPerSecond)
        && writer.finiteDouble(settings.relativeContinuityTolerance)
        && writer.finiteDouble(
            settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && writer.finiteDouble(settings.relativeMomentumTolerance)
        && writer.finiteDouble(settings.absoluteEnergyToleranceJoules)
        && writer.finiteDouble(settings.relativeEnergyTolerance);
}

bool readTransportSettings(
    Reader& reader,
    PlanarPressureRegionFragmentOpeningMomentumTransportSettings& settings) {
    return reader.finiteDouble(settings.maximumOutgoingCourantNumber)
        && reader.count(
            settings.maximumSubsteps,
            std::numeric_limits<std::size_t>::max())
        && reader.finiteDouble(
            settings.absoluteContinuityToleranceCubicMetersPerSecond)
        && reader.finiteDouble(settings.relativeContinuityTolerance)
        && reader.finiteDouble(
            settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && reader.finiteDouble(settings.relativeMomentumTolerance)
        && reader.finiteDouble(settings.absoluteEnergyToleranceJoules)
        && reader.finiteDouble(settings.relativeEnergyTolerance);
}

bool writeTransportDiagnostics(
    Writer& writer,
    const PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics&
        diagnostics) {
    return writer.count(diagnostics.fragmentCount)
        && writer.count(diagnostics.transportDofCount)
        && writer.count(diagnostics.openingDofCount)
        && writer.count(diagnostics.substepCount)
        && writer.finiteDouble(
            diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.maximumContinuityResidualCubicMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.continuityToleranceCubicMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.maximumFullStepOutgoingCourantNumber)
        && writer.finiteDouble(
            diagnostics.maximumAcceptedSubstepOutgoingCourantNumber)
        && writeVector(
            writer, diagnostics.momentumBeforeKilogramMetersPerSecond)
        && writeVector(
            writer, diagnostics.momentumAfterKilogramMetersPerSecond)
        && writeVector(
            writer, diagnostics.momentumResidualKilogramMetersPerSecond)
        && writer.finiteDouble(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        && writer.finiteDouble(diagnostics.kineticEnergyBeforeJoules)
        && writer.finiteDouble(diagnostics.kineticEnergyAfterJoules)
        && writer.finiteDouble(
            diagnostics.advectiveKineticEnergyLossJoules)
        && writer.finiteDouble(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        && writer.u8(static_cast<std::uint8_t>(diagnostics.failureStage))
        && writer.u8(diagnostics.finite ? 1U : 0U)
        && writer.u8(diagnostics.accepted ? 1U : 0U)
        && writer.zeroes(5);
}

bool readTransportDiagnostics(
    Reader& reader,
    PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics&
        diagnostics,
    const Limits& limits) {
    std::uint8_t failureStage = 0;
    if (!reader.count(
            diagnostics.fragmentCount,
            limits.maximumTransportControls)
        || !reader.count(
            diagnostics.transportDofCount,
            limits.state.transport.maximumDofs)
        || !reader.count(
            diagnostics.openingDofCount,
            limits.state.transport.maximumDofs)
        || !reader.count(
            diagnostics.substepCount,
            std::numeric_limits<std::size_t>::max())
        || !reader.finiteDouble(
            diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.maximumContinuityResidualCubicMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.continuityToleranceCubicMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.maximumFullStepOutgoingCourantNumber)
        || !reader.finiteDouble(
            diagnostics.maximumAcceptedSubstepOutgoingCourantNumber)
        || !readVector(
            reader, diagnostics.momentumBeforeKilogramMetersPerSecond)
        || !readVector(
            reader, diagnostics.momentumAfterKilogramMetersPerSecond)
        || !readVector(
            reader, diagnostics.momentumResidualKilogramMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        || !reader.finiteDouble(diagnostics.kineticEnergyBeforeJoules)
        || !reader.finiteDouble(diagnostics.kineticEnergyAfterJoules)
        || !reader.finiteDouble(
            diagnostics.advectiveKineticEnergyLossJoules)
        || !reader.finiteDouble(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        || !reader.u8(failureStage)
        || !reader.boolean(diagnostics.finite)
        || !reader.boolean(diagnostics.accepted)
        || !reader.zeroes(5)) {
        return false;
    }
    if (failureStage
        > static_cast<std::uint8_t>(
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                NonFinite)) {
        return false;
    }
    diagnostics.failureStage = static_cast<
        PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage>(
        failureStage);
    return true;
}

bool writeTransport(
    Writer& writer,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport) {
    if (!writer.u32(transport.version) || !writer.u32(0)
        || !writer.u64(transport.fingerprint)
        || !writer.u64(transport.sourceStateFingerprint)
        || !writer.u64(transport.sourceTransportFingerprint)
        || !writer.u64(transport.sourceMetricFingerprint)
        || !writer.u64(transport.targetFlowStateFingerprint)
        || !writer.u64(transport.targetMetricFingerprint)
        || !writer.u64(transport.targetVolumeRateFingerprint)
        || !writer.finiteDouble(transport.densityKgPerCubicMeter)
        || !writer.finiteDouble(transport.timeStepSeconds)
        || !writeTransportSettings(writer, transport.settings)
        || !writeTransportDiagnostics(writer, transport.diagnostics)
        || !writer.count(transport.ownedStorageBytes)
        || !writer.count(transport.workingStorageBytes)
        || !writer.count(transport.controls.size())) {
        return false;
    }
    for (const auto& control : transport.controls) {
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

bool readTransport(
    Reader& reader,
    PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const Limits& limits) {
    std::uint32_t reserved = 0;
    std::size_t controlCount = 0;
    if (!reader.u32(transport.version) || !reader.u32(reserved)
        || reserved != 0 || !reader.u64(transport.fingerprint)
        || !reader.u64(transport.sourceStateFingerprint)
        || !reader.u64(transport.sourceTransportFingerprint)
        || !reader.u64(transport.sourceMetricFingerprint)
        || !reader.u64(transport.targetFlowStateFingerprint)
        || !reader.u64(transport.targetMetricFingerprint)
        || !reader.u64(transport.targetVolumeRateFingerprint)
        || !reader.finiteDouble(transport.densityKgPerCubicMeter)
        || !reader.finiteDouble(transport.timeStepSeconds)
        || !readTransportSettings(reader, transport.settings)
        || !readTransportDiagnostics(reader, transport.diagnostics, limits)
        || !reader.count(
            transport.ownedStorageBytes,
            limits.state.transport.maximumOwnedBytes)
        || !reader.count(
            transport.workingStorageBytes,
            limits.state.transport.maximumWorkingBytes)
        || !reader.count(controlCount, limits.maximumTransportControls)) {
        return false;
    }
    if (!reader.fixedRecords(controlCount, transportControlRecordBytes)) {
        return false;
    }
    transport.controls.resize(controlCount);
    for (auto& control : transport.controls) {
        if (!reader.count(
                control.fragmentIndex,
                limits.maximumTransportControls)
            || !reader.u64(control.stableId)
            || !reader.u64(control.regionStableId)
            || !reader.count(
                control.connectedComponentIndex,
                limits.state.transport.maximumFragments)
            || !reader.finiteDouble(control.volumeCubicMeters)
            || !readVector(reader, control.velocityMetersPerSecond)
            || !readVector(
                reader, control.momentumKilogramMetersPerSecond)) {
            return false;
        }
    }
    return true;
}

ErrorCode nestedErrorCode(
    const PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode
        code) {
    using Nested =
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

} // namespace

bool serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    const PlanarPressureRegionFragmentVolumeRateSet& transportVolumeRates,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& acceptedSweep,
    const PlanarPressureRegionFragmentSet& acceptedFragments,
    const PlanarPressureRegionFragmentTopology& acceptedTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    std::vector<std::uint8_t>& bytes,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "opening momentum-cycle persistence limits are invalid");
    }
    try {
        validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
            state, transportVolumeRates, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedMetric, limits.state);
        if (state.transport.controls.size()
            > limits.maximumTransportControls) {
            return fail(error, ErrorCode::LimitExceeded,
                        "opening momentum-cycle transport count exceeds persistence limit");
        }

        std::vector<std::uint8_t> acceptedBytes;
        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
            acceptedError;
        if (!serializePlanarPressureRegionFragmentOpeningAcceptedState(
                state.acceptedState, acceptedPressureOperator,
                acceptedBasePressureOperator, grid, acceptedSweep,
                acceptedFragments, acceptedTopology, acceptedVolumeRates,
                acceptedOpeningDefinitions, acceptedOpenings,
                acceptedResistanceDefinitions, acceptedBytes,
                &acceptedError, acceptedPersistenceLimits(limits))) {
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
            || !payloadWriter.u64(state.transportMetricFingerprint)
            || !payloadWriter.u64(state.acceptedMetricFingerprint)
            || !payloadWriter.u64(
                state.currentAcceptedFlowStateFingerprint)
            || !payloadWriter.u64(state.predictionFingerprint)
            || !payloadWriter.u64(state.pressureWarmStartFingerprint)
            || !payloadWriter.u64(state.predictedOpeningFluxFingerprint)
            || !payloadWriter.count(state.ownedStorageBytes)
            || !writeTransport(payloadWriter, state.transport)
            || !payloadWriter.count(acceptedBytes.size())
            || !payloadWriter.bytes(acceptedBytes)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? ErrorCode::LimitExceeded : ErrorCode::InvalidData,
                "opening momentum-cycle payload is invalid or too large");
        }

        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) {
                return fail(error, ErrorCode::LimitExceeded,
                            "opening momentum-cycle envelope is too large");
            }
        }
        if (!writer.u16(
                planarPressureRegionFragmentOpeningMomentumCycleStateProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "opening momentum-cycle envelope is too large");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate encoded opening momentum-cycle state");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(
            error, ErrorCode::InvalidData,
            std::string("opening momentum-cycle state is invalid: ")
                + exception.what());
    }
}

bool deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
    const std::span<const std::uint8_t> bytes,
    const PlanarPressureRegionFragmentVolumeRateSet& transportVolumeRates,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& transportMetric,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& acceptedSweep,
    const PlanarPressureRegionFragmentSet& acceptedFragments,
    const PlanarPressureRegionFragmentTopology& acceptedTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& acceptedVolumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& acceptedMetric,
    PlanarPressureRegionFragmentOpeningMomentumCycleState& state,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "opening momentum-cycle persistence limits are invalid");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error, ErrorCode::LimitExceeded,
                    "opening momentum-cycle envelope exceeds byte limit");
    }

    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) {
            return fail(error, ErrorCode::Truncated,
                        "opening momentum-cycle envelope is truncated");
        }
        if (actual != expected) {
            return fail(error, ErrorCode::InvalidMagic,
                        "opening momentum-cycle magic is invalid");
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
                    "opening momentum-cycle envelope header is invalid");
    }
    if (protocol
        != planarPressureRegionFragmentOpeningMomentumCycleStateProtocolVersion) {
        return fail(error, ErrorCode::UnsupportedVersion,
                    "opening momentum-cycle protocol is unsupported");
    }
    if (reserved != 0) {
        return fail(error, ErrorCode::InvalidData,
                    "opening momentum-cycle envelope reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error, ErrorCode::Truncated,
                    "opening momentum-cycle payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error, ErrorCode::TrailingData,
                    "opening momentum-cycle envelope has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error, ErrorCode::ChecksumMismatch,
                    "opening momentum-cycle payload checksum differs");
    }

    try {
        Reader reader(payload);
        PlanarPressureRegionFragmentOpeningMomentumCycleState candidate;
        std::uint32_t decodedPayloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t stateReserved = 0;
        if (!reader.u32(decodedPayloadVersion)
            || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(stateReserved)) {
            return fail(error, readerErrorCode(reader),
                        "opening momentum-cycle payload header is invalid");
        }
        if (decodedPayloadVersion != payloadVersion) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "opening momentum-cycle payload version is unsupported");
        }
        if (payloadReserved != 0 || stateReserved != 0) {
            return fail(error, ErrorCode::InvalidData,
                        "opening momentum-cycle payload reserved bits are nonzero");
        }
        if (!reader.u64(candidate.fingerprint)
            || !reader.u64(candidate.transportMetricFingerprint)
            || !reader.u64(candidate.acceptedMetricFingerprint)
            || !reader.u64(
                candidate.currentAcceptedFlowStateFingerprint)
            || !reader.u64(candidate.predictionFingerprint)
            || !reader.u64(candidate.pressureWarmStartFingerprint)
            || !reader.u64(candidate.predictedOpeningFluxFingerprint)
            || !reader.count(
                candidate.ownedStorageBytes,
                limits.state.maximumOwnedBytes)
            || !readTransport(reader, candidate.transport, limits)) {
            return fail(error, readerErrorCode(reader),
                        "opening momentum-cycle payload is invalid");
        }
        if (candidate.transportMetricFingerprint
                != transportMetric.fingerprint
            || candidate.transport.targetMetricFingerprint
                != transportMetric.fingerprint
            || candidate.transport.targetVolumeRateFingerprint
                != transportVolumeRates.fingerprint
            || candidate.acceptedMetricFingerprint
                != acceptedMetric.fingerprint) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening momentum-cycle state belongs to different restart geometry");
        }

        std::size_t acceptedByteCount = 0;
        if (!reader.count(
                acceptedByteCount,
                limits.maximumEncodedBytes - envelopeBytes)) {
            return fail(error, readerErrorCode(reader),
                        "opening momentum-cycle nested state size is invalid");
        }
        const auto acceptedBytes = reader.remainingBytes(acceptedByteCount);
        if (reader.truncated()) {
            return fail(error, ErrorCode::Truncated,
                        "opening momentum-cycle nested state is truncated");
        }
        if (!reader.atEnd()) {
            return fail(error, ErrorCode::TrailingData,
                        "opening momentum-cycle payload has trailing data");
        }

        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
            acceptedError;
        if (!deserializePlanarPressureRegionFragmentOpeningAcceptedState(
                acceptedBytes, acceptedPressureOperator,
                acceptedBasePressureOperator, grid, acceptedSweep,
                acceptedFragments, acceptedTopology, acceptedVolumeRates,
                acceptedOpeningDefinitions, acceptedOpenings,
                acceptedResistanceDefinitions, candidate.acceptedState,
                &acceptedError, acceptedPersistenceLimits(limits))) {
            return fail(
                error, nestedErrorCode(acceptedError.code),
                std::string("cannot decode nested opening accepted state: ")
                    + acceptedError.message);
        }
        validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
            candidate, transportVolumeRates, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedMetric, limits.state);
        state = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate decoded opening momentum-cycle state");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(
            error, ErrorCode::InvalidData,
            std::string("decoded opening momentum-cycle state is invalid: ")
                + exception.what());
    }
}

} // namespace simwing::fsi::fluid
