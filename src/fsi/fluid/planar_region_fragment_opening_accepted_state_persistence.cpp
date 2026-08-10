#include "fluid/planar_region_fragment_opening_accepted_state_persistence.h"

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

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'R', 'O'};
constexpr std::uint32_t payloadVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t topologyVelocityRecordBytes = 8;
constexpr std::size_t openingSampleRecordBytes = 16;
constexpr std::size_t pressureCorrectionRecordBytes = 8;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

using ErrorCode =
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode;
using Error =
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError;
using Limits =
    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits;

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
        && limits.maximumTopologyLinkVelocities > 0
        && limits.maximumOpeningSamples > 0
        && limits.maximumPressureCorrections > 0;
}

bool writeSettings(
    Writer& writer,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings) {
    const auto& projection = settings.projection;
    const auto& solve = projection.pressureSolve;
    return writer.finiteDouble(projection.densityKgPerCubicMeter)
        && writer.finiteDouble(projection.timeStepSeconds)
        && writer.finiteDouble(
            projection.absoluteContinuityToleranceCubicMetersPerSecond)
        && writer.finiteDouble(projection.relativeContinuityTolerance)
        && writer.finiteDouble(
            projection
                .absoluteMomentumResidualToleranceKilogramMetersPerSecond)
        && writer.finiteDouble(
            projection.relativeMomentumResidualTolerance)
        && writer.finiteDouble(
            projection.absoluteEnergyResidualToleranceJoules)
        && writer.finiteDouble(projection.relativeEnergyResidualTolerance)
        && writer.finiteDouble(
            solve.absoluteResidualTolerancePascalsMeters)
        && writer.finiteDouble(solve.relativeResidualTolerance)
        && writer.finiteDouble(
            solve.absoluteComponentCompatibilityTolerancePascalsMeters)
        && writer.count(solve.maximumIterations)
        && writer.u8(settings.useAuthoredPressureDrive ? 1U : 0U)
        && writer.zeroes(7);
}

bool readSettings(
    Reader& reader,
    PlanarPressureRegionFragmentOpeningPressureStepSettings& settings) {
    auto& projection = settings.projection;
    auto& solve = projection.pressureSolve;
    return reader.finiteDouble(projection.densityKgPerCubicMeter)
        && reader.finiteDouble(projection.timeStepSeconds)
        && reader.finiteDouble(
            projection.absoluteContinuityToleranceCubicMetersPerSecond)
        && reader.finiteDouble(projection.relativeContinuityTolerance)
        && reader.finiteDouble(
            projection
                .absoluteMomentumResidualToleranceKilogramMetersPerSecond)
        && reader.finiteDouble(
            projection.relativeMomentumResidualTolerance)
        && reader.finiteDouble(
            projection.absoluteEnergyResidualToleranceJoules)
        && reader.finiteDouble(projection.relativeEnergyResidualTolerance)
        && reader.finiteDouble(
            solve.absoluteResidualTolerancePascalsMeters)
        && reader.finiteDouble(solve.relativeResidualTolerance)
        && reader.finiteDouble(
            solve.absoluteComponentCompatibilityTolerancePascalsMeters)
        && reader.count(
            solve.maximumIterations,
            std::numeric_limits<std::size_t>::max())
        && reader.boolean(settings.useAuthoredPressureDrive)
        && reader.zeroes(7);
}

bool writePayload(
    Writer& writer,
    const PlanarPressureRegionFragmentOpeningAcceptedState& state) {
    if (!writer.u32(payloadVersion) || !writer.u32(0)
        || !writer.u32(state.version) || !writer.u32(0)
        || !writer.u64(state.fingerprint)
        || !writer.u64(state.sourcePressureOperatorFingerprint)
        || !writer.u64(state.sourceBasePressureOperatorFingerprint)
        || !writer.u64(state.sourceOpeningFingerprint)
        || !writer.u64(state.sourceFragmentFingerprint)
        || !writer.u64(state.sourceTopologyFingerprint)
        || !writer.u64(state.sourceVolumeRateFingerprint)
        || !writer.u64(state.sourceOpeningFluxFingerprint)
        || !writer.u64(state.resultOpeningFluxFingerprint)
        || !writer.u64(state.resistanceDefinitionFingerprint)
        || !writeSettings(writer, state.settings)
        || !writer.count(state.ownedStorageBytes)
        || !writer.count(state.workingStorageBytes)
        || !writer.finiteDouble(
            state.correctedContinuityResidualL2CubicMetersPerSecond)
        || !writer.finiteDouble(
            state.correctedContinuityResidualMaximumCubicMetersPerSecond)
        || !writer.finiteDouble(
            state.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond)
        || !writer.finiteDouble(
            state.continuityToleranceCubicMetersPerSecond)
        || !writer.finiteDouble(
            state.maximumAbsoluteCorrectionVolumeMeanPascals)
        || !writer.finiteDouble(
            state.maximumAbsolutePressureCorrectionPascals)
        || !writer.finiteDouble(
            state.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond)
        || !writer.finiteDouble(
            state.momentumResidualToleranceKilogramMetersPerSecond)
        || !writer.finiteDouble(state.kineticEnergyBeforeJoules)
        || !writer.finiteDouble(state.kineticEnergyAfterJoules)
        || !writer.finiteDouble(state.kineticEnergyChangeJoules)
        || !writer.finiteDouble(state.authoredPressureWorkJoules)
        || !writer.finiteDouble(state.geometryPressureWorkJoules)
        || !writer.finiteDouble(state.correctionKineticEnergyJoules)
        || !writer.finiteDouble(state.dissipatedEnergyJoules)
        || !writer.finiteDouble(state.energyResidualJoules)
        || !writer.finiteDouble(state.energyToleranceJoules)
        || !writer.count(state.pressureSolveIterationCount)
        || !writer.finiteDouble(
            state.pressureSolveFinalResidualL2PascalsMeters)
        || !writer.finiteDouble(
            state.pressureSolveFinalResidualMaximumPascalsMeters)
        || !writer.u8(state.accepted ? 1U : 0U) || !writer.zeroes(7)
        || !writer.count(
            state.orientedTopologyLinkVelocityMetersPerSecond.size())) {
        return false;
    }
    for (const double value
         : state.orientedTopologyLinkVelocityMetersPerSecond) {
        if (!writer.finiteDouble(value)) return false;
    }
    if (!writer.count(state.openingVelocitySamples.size())) return false;
    for (const auto& sample : state.openingVelocitySamples) {
        if (!writer.u64(sample.patchStableId)
            || !writer.finiteDouble(
                sample.relativeNormalVelocityMetersPerSecond)) {
            return false;
        }
    }
    if (!writer.count(state.pressureCorrectionPascals.size())) return false;
    for (const double value : state.pressureCorrectionPascals) {
        if (!writer.finiteDouble(value)) return false;
    }
    return true;
}

bool sourceIdentityMatches(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    return state.sourcePressureOperatorFingerprint
            == pressureOperator.fingerprint
        && state.sourceBasePressureOperatorFingerprint
            == basePressureOperator.fingerprint
        && state.sourceOpeningFingerprint == openings.fingerprint
        && state.sourceFragmentFingerprint == fragments.fingerprint
        && state.sourceTopologyFingerprint == topology.fingerprint
        && state.sourceVolumeRateFingerprint == volumeRates.fingerprint
        && state.resultOpeningFluxFingerprint != 0
        && state.sourceOpeningFluxFingerprint != 0
        && state.orientedTopologyLinkVelocityMetersPerSecond.size()
            == topology.links.size()
        && state.openingVelocitySamples.size() == openings.patches.size()
        && state.pressureCorrectionPascals.size()
            == pressureOperator.rows.size();
}

} // namespace

bool serializePlanarPressureRegionFragmentOpeningAcceptedState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    std::vector<std::uint8_t>& bytes,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "invalid opening accepted-state persistence limits");
    }
    try {
        validatePlanarPressureRegionFragmentOpeningAcceptedState(
            state, pressureOperator, basePressureOperator, grid, sweep,
            fragments, topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, limits.stateLimits);
        if (state.orientedTopologyLinkVelocityMetersPerSecond.size()
                > limits.maximumTopologyLinkVelocities
            || state.openingVelocitySamples.size()
                > limits.maximumOpeningSamples
            || state.pressureCorrectionPascals.size()
                > limits.maximumPressureCorrections) {
            return fail(error, ErrorCode::LimitExceeded,
                        "opening accepted state exceeds record limits");
        }
        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload, limits.maximumEncodedBytes - envelopeBytes);
        if (!writePayload(payloadWriter, state)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? ErrorCode::LimitExceeded : ErrorCode::InvalidData,
                "opening accepted-state payload is invalid or too large");
        }
        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) {
                return fail(error, ErrorCode::LimitExceeded,
                            "opening accepted-state envelope is too large");
            }
        }
        if (!writer.u16(
                planarPressureRegionFragmentOpeningAcceptedStateProtocolVersion)
            || !writer.u16(0) || !writer.count(payload.size())
            || !writer.u64(checksum(payload)) || !writer.bytes(payload)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "opening accepted state exceeds its encoded limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate opening accepted-state encoding");
    } catch (const std::length_error&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "opening accepted-state encoding is too large");
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData,
                    std::string("opening accepted-state validation failed: ")
                        + exception.what());
    }
}

bool deserializePlanarPressureRegionFragmentOpeningAcceptedState(
    const std::span<const std::uint8_t> bytes,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    PlanarPressureRegionFragmentOpeningAcceptedState& state,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::InvalidData,
                    "invalid opening accepted-state persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error, ErrorCode::LimitExceeded,
                    "opening accepted state exceeds its encoded limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error, ErrorCode::Truncated,
                    "opening accepted-state envelope is truncated");
    }
    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual) || actual != expected) {
            return fail(error, ErrorCode::InvalidMagic,
                        "opening accepted-state magic is invalid");
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
                    "opening accepted-state envelope is invalid");
    }
    if (protocol
        != planarPressureRegionFragmentOpeningAcceptedStateProtocolVersion) {
        return fail(error, ErrorCode::UnsupportedVersion,
                    "opening accepted-state protocol is unsupported");
    }
    if (reserved != 0) {
        return fail(error, ErrorCode::InvalidData,
                    "opening accepted-state reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error, ErrorCode::Truncated,
                    "opening accepted-state payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error, ErrorCode::TrailingData,
                    "opening accepted state has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error, ErrorCode::ChecksumMismatch,
                    "opening accepted-state checksum does not match");
    }

    try {
        validatePlanarPressureRegionFragmentOpeningPressureOperator(
            pressureOperator, basePressureOperator, grid, sweep, fragments,
            topology, openingDefinitions, openings,
            limits.stateLimits.pressureOperatorLimits);
        validatePlanarPressureRegionFragmentVolumeRates(
            volumeRates, grid, sweep, fragments, topology,
            limits.stateLimits.volumeRateLimits);
        Reader reader(payload);
        PlanarPressureRegionFragmentOpeningAcceptedState candidate;
        std::uint32_t decodedPayloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t stateReserved = 0;
        if (!reader.u32(decodedPayloadVersion)
            || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(stateReserved)
            || !reader.u64(candidate.fingerprint)
            || !reader.u64(candidate.sourcePressureOperatorFingerprint)
            || !reader.u64(candidate.sourceBasePressureOperatorFingerprint)
            || !reader.u64(candidate.sourceOpeningFingerprint)
            || !reader.u64(candidate.sourceFragmentFingerprint)
            || !reader.u64(candidate.sourceTopologyFingerprint)
            || !reader.u64(candidate.sourceVolumeRateFingerprint)
            || !reader.u64(candidate.sourceOpeningFluxFingerprint)
            || !reader.u64(candidate.resultOpeningFluxFingerprint)
            || !reader.u64(candidate.resistanceDefinitionFingerprint)
            || !readSettings(reader, candidate.settings)
            || !reader.count(
                candidate.ownedStorageBytes,
                std::numeric_limits<std::size_t>::max())
            || !reader.count(
                candidate.workingStorageBytes,
                std::numeric_limits<std::size_t>::max())
            || !reader.finiteDouble(
                candidate.correctedContinuityResidualL2CubicMetersPerSecond)
            || !reader.finiteDouble(
                candidate.correctedContinuityResidualMaximumCubicMetersPerSecond)
            || !reader.finiteDouble(
                candidate.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond)
            || !reader.finiteDouble(
                candidate.continuityToleranceCubicMetersPerSecond)
            || !reader.finiteDouble(
                candidate.maximumAbsoluteCorrectionVolumeMeanPascals)
            || !reader.finiteDouble(
                candidate.maximumAbsolutePressureCorrectionPascals)
            || !reader.finiteDouble(
                candidate.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond)
            || !reader.finiteDouble(
                candidate.momentumResidualToleranceKilogramMetersPerSecond)
            || !reader.finiteDouble(candidate.kineticEnergyBeforeJoules)
            || !reader.finiteDouble(candidate.kineticEnergyAfterJoules)
            || !reader.finiteDouble(candidate.kineticEnergyChangeJoules)
            || !reader.finiteDouble(candidate.authoredPressureWorkJoules)
            || !reader.finiteDouble(candidate.geometryPressureWorkJoules)
            || !reader.finiteDouble(
                candidate.correctionKineticEnergyJoules)
            || !reader.finiteDouble(candidate.dissipatedEnergyJoules)
            || !reader.finiteDouble(candidate.energyResidualJoules)
            || !reader.finiteDouble(candidate.energyToleranceJoules)
            || !reader.count(
                candidate.pressureSolveIterationCount,
                std::numeric_limits<std::size_t>::max())
            || !reader.finiteDouble(
                candidate.pressureSolveFinalResidualL2PascalsMeters)
            || !reader.finiteDouble(
                candidate.pressureSolveFinalResidualMaximumPascalsMeters)
            || !reader.boolean(candidate.accepted)
            || !reader.zeroes(7)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state payload header is invalid");
        }
        if (decodedPayloadVersion != payloadVersion
            || candidate.version
                != planarPressureRegionFragmentOpeningAcceptedStateVersion) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "opening accepted-state payload version is unsupported");
        }
        if (payloadReserved != 0 || stateReserved != 0) {
            return fail(error, ErrorCode::InvalidData,
                        "opening accepted-state payload reserved bits are nonzero");
        }

        std::size_t topologyVelocityCount = 0;
        if (!reader.count(
                topologyVelocityCount,
                limits.maximumTopologyLinkVelocities)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state topology velocity count is invalid");
        }
        if (topologyVelocityCount != topology.links.size()) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening accepted-state topology velocity count does not match sources");
        }
        if (!reader.fixedRecords(
                topologyVelocityCount, topologyVelocityRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state topology velocities are truncated");
        }
        candidate.orientedTopologyLinkVelocityMetersPerSecond.resize(
            topologyVelocityCount);
        for (double& value
             : candidate.orientedTopologyLinkVelocityMetersPerSecond) {
            if (!reader.finiteDouble(value)) {
                return fail(error, readerErrorCode(reader),
                            "opening accepted-state topology velocity is invalid");
            }
        }

        std::size_t openingSampleCount = 0;
        if (!reader.count(
                openingSampleCount, limits.maximumOpeningSamples)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state aperture sample count is invalid");
        }
        if (openingSampleCount != openings.patches.size()) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening accepted-state aperture count does not match sources");
        }
        if (!reader.fixedRecords(
                openingSampleCount, openingSampleRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state aperture samples are truncated");
        }
        candidate.openingVelocitySamples.resize(openingSampleCount);
        for (auto& sample : candidate.openingVelocitySamples) {
            if (!reader.u64(sample.patchStableId)
                || !reader.finiteDouble(
                    sample.relativeNormalVelocityMetersPerSecond)) {
                return fail(error, readerErrorCode(reader),
                            "opening accepted-state aperture sample is invalid");
            }
        }

        std::size_t pressureCorrectionCount = 0;
        if (!reader.count(
                pressureCorrectionCount,
                limits.maximumPressureCorrections)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state pressure count is invalid");
        }
        if (pressureCorrectionCount != pressureOperator.rows.size()) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening accepted-state pressure count does not match sources");
        }
        if (!reader.fixedRecords(
                pressureCorrectionCount, pressureCorrectionRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "opening accepted-state pressure corrections are truncated");
        }
        candidate.pressureCorrectionPascals.resize(
            pressureCorrectionCount);
        for (double& value : candidate.pressureCorrectionPascals) {
            if (!reader.finiteDouble(value)) {
                return fail(error, readerErrorCode(reader),
                            "opening accepted-state pressure correction is invalid");
            }
        }
        if (!reader.atEnd()) {
            return fail(error, ErrorCode::TrailingData,
                        "opening accepted-state payload has trailing data");
        }
        if (!sourceIdentityMatches(
                candidate, pressureOperator, basePressureOperator,
                fragments, topology, volumeRates, openings)) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening accepted state belongs to different sources");
        }
        candidate.openingFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                grid, sweep, fragments, topology, openingDefinitions,
                openings, candidate.openingVelocitySamples,
                limits.stateLimits.openingFluxLimits);
        if (candidate.openingFlux.fingerprint
            != candidate.resultOpeningFluxFingerprint) {
            return fail(error, ErrorCode::SourceMismatch,
                        "opening accepted-state rebuilt flux does not match sources");
        }
        validatePlanarPressureRegionFragmentOpeningAcceptedState(
            candidate, pressureOperator, basePressureOperator, grid, sweep,
            fragments, topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions, limits.stateLimits);
        state = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate decoded opening accepted state");
    } catch (const std::length_error&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "decoded opening accepted state is too large");
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData,
                    std::string("decoded opening accepted state is invalid: ")
                        + exception.what());
    }
}

} // namespace simwing::fsi::fluid
