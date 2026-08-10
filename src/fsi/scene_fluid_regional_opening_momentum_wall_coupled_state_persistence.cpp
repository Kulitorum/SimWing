#include "scene_fluid_regional_opening_momentum_wall_coupled_state_persistence.h"

#include <algorithm>
#include <array>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'R', 'C'};
constexpr std::uint32_t stateVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

using ErrorCode =
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode;
using Error =
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError;
using Limits =
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits;

void clearError(Error* error) {
    if (error != nullptr) *error = {};
}

bool fail(Error* error, const ErrorCode code, std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

std::uint64_t checksum(const std::span<const std::uint8_t> bytes) {
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
        return raw(std::span<const std::uint8_t>(&value, 1));
    }

    bool u16(std::uint16_t value) { return integer(value); }
    bool u32(std::uint32_t value) { return integer(value); }
    bool u64(std::uint64_t value) { return integer(value); }

    bool count(const std::size_t value) {
        return u64(static_cast<std::uint64_t>(value));
    }

    bool raw(const std::span<const std::uint8_t> values) {
        if (values.size() > limit_ - bytes_.size()) return false;
        bytes_.insert(bytes_.end(), values.begin(), values.end());
        return true;
    }

private:
    template<typename Unsigned>
    bool integer(Unsigned value) {
        std::array<std::uint8_t, sizeof(value)> encoded{};
        for (std::size_t index = 0; index < encoded.size(); ++index) {
            encoded[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return raw(encoded);
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t limit_ = 0;
};

class Reader final {
public:
    explicit Reader(const std::span<const std::uint8_t> bytes)
        : bytes_(bytes) {}

    bool u8(std::uint8_t& value) {
        if (remaining() < 1) {
            truncated_ = true;
            return false;
        }
        value = bytes_[offset_++];
        return true;
    }

    bool u16(std::uint16_t& value) { return integer(value); }
    bool u32(std::uint32_t& value) { return integer(value); }
    bool u64(std::uint64_t& value) { return integer(value); }

    bool count(std::size_t& value, const std::size_t limit) {
        std::uint64_t encoded = 0;
        if (!u64(encoded)
            || encoded > std::numeric_limits<std::size_t>::max()) {
            return false;
        }
        if (encoded > limit) {
            limitExceeded_ = true;
            return false;
        }
        value = static_cast<std::size_t>(encoded);
        return true;
    }

    std::span<const std::uint8_t> raw(const std::size_t count) {
        if (count > remaining()) {
            truncated_ = true;
            return {};
        }
        const auto result = bytes_.subspan(offset_, count);
        offset_ += count;
        return result;
    }

    [[nodiscard]] std::size_t remaining() const noexcept {
        return bytes_.size() - offset_;
    }

    [[nodiscard]] bool truncated() const noexcept { return truncated_; }
    [[nodiscard]] bool limitExceeded() const noexcept {
        return limitExceeded_;
    }

private:
    template<typename Unsigned>
    bool integer(Unsigned& value) {
        if (remaining() < sizeof(value)) {
            truncated_ = true;
            return false;
        }
        value = 0;
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            value |= static_cast<Unsigned>(bytes_[offset_++])
                << (8U * index);
        }
        return true;
    }

    std::span<const std::uint8_t> bytes_;
    std::size_t offset_ = 0;
    bool truncated_ = false;
    bool limitExceeded_ = false;
};

ErrorCode readerErrorCode(const Reader& reader) {
    return reader.limitExceeded()
        ? ErrorCode::LimitExceeded
        : ErrorCode::Truncated;
}

bool validLimits(const Limits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.cycleState.maximumEncodedBytes >= envelopeBytes
        && limits.state.maximumOwnedBytes > 0
        && limits.state.structureStep.checkpointPersistence
               .maximumEncodedBytes >= envelopeBytes;
}

ErrorCode cycleErrorCode(
    const SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode
        code) {
    using CycleCode =
        SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode;
    switch (code) {
    case CycleCode::LimitExceeded:
        return ErrorCode::LimitExceeded;
    case CycleCode::SourceMismatch:
        return ErrorCode::SourceMismatch;
    case CycleCode::InvalidMagic:
    case CycleCode::UnsupportedVersion:
    case CycleCode::Truncated:
    case CycleCode::TrailingData:
    case CycleCode::ChecksumMismatch:
    case CycleCode::InvalidData:
    case CycleCode::None:
        return ErrorCode::InvalidData;
    }
    return ErrorCode::InvalidData;
}

ErrorCode structureErrorCode(
    const StructureCheckpointPersistenceErrorCode code) {
    if (code == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
        return ErrorCode::LimitExceeded;
    }
    if (code == StructureCheckpointPersistenceErrorCode::TopologyMismatch) {
        return ErrorCode::TopologyMismatch;
    }
    return ErrorCode::InvalidData;
}

bool encodeCycleState(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
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
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError nested;
    if (!serializeSceneFluidRegionalOpeningMomentumWallCycleState(
            state.cycleState, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, quadrature, bytes, &nested,
            limits.cycleState)) {
        return fail(error, cycleErrorCode(nested.code),
                    "nested SWRW state: " + nested.message);
    }
    return true;
}

bool decodeCycleState(
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
    SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError nested;
    if (!deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
            bytes, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, quadrature, state, &nested,
            limits.cycleState)) {
        return fail(error, cycleErrorCode(nested.code),
                    "nested SWRW state: " + nested.message);
    }
    return true;
}

} // namespace

bool serializeSceneFluidRegionalOpeningMomentumWallCoupledState(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& structureDefinitionOwner,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    std::vector<std::uint8_t>& bytes,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error, ErrorCode::LimitExceeded,
                    "regional opening wall coupled-state persistence limits are invalid");
    }
    try {
        validateSceneFluidRegionalOpeningMomentumWallCoupledState(
            state, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, surface, surfaceState, transfer, quadrature,
            structureDefinitionOwner, settings, limits.state);

        std::vector<std::uint8_t> cycleBytes;
        if (!encodeCycleState(
                state, transportMetric, acceptedPressureOperator,
                acceptedBasePressureOperator, grid, acceptedSweep,
                acceptedFragments, acceptedTopology, acceptedVolumeRates,
                acceptedOpeningDefinitions, acceptedOpenings,
                acceptedResistanceDefinitions, acceptedBaseMetric,
                acceptedMetric, quadrature, cycleBytes, error, limits)) {
            return false;
        }
        const auto& before = state.structureStep.beforeStructureCheckpoint;
        const auto& after = state.structureStep.afterStructureCheckpoint;
        const std::size_t checkpointLimit = limits.state.structureStep
            .checkpointPersistence.maximumEncodedBytes;
        if (before.size() > checkpointLimit || after.size() > checkpointLimit) {
            return fail(error, ErrorCode::LimitExceeded,
                        "nested Structure checkpoint exceeds its byte limit");
        }

        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload, limits.maximumEncodedBytes - envelopeBytes);
        if (!payloadWriter.u32(stateVersion)
            || !payloadWriter.u32(state.version)
            || !payloadWriter.u64(state.fingerprint)
            || !payloadWriter.count(cycleBytes.size())
            || !payloadWriter.count(before.size())
            || !payloadWriter.count(after.size())
            || !payloadWriter.raw(cycleBytes)
            || !payloadWriter.raw(before)
            || !payloadWriter.raw(after)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "regional opening wall coupled-state payload exceeds its byte limit");
        }
        std::vector<std::uint8_t> encoded;
        Writer writer(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t value : magic) {
            if (!writer.u8(value)) return false;
        }
        if (!writer.u16(
                sceneFluidRegionalOpeningMomentumWallCoupledStateProtocolVersion)
            || !writer.u16(0)
            || !writer.count(payload.size())
            || !writer.u64(checksum(payload))
            || !writer.raw(payload)) {
            return fail(error, ErrorCode::LimitExceeded,
                        "regional opening wall coupled-state envelope exceeds its byte limit");
        }
        bytes.swap(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate regional opening wall coupled-state encoding");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData, exception.what());
    }
}

bool deserializeSceneFluidRegionalOpeningMomentumWallCoupledState(
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& structureDefinitionOwner,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    Error* error,
    const Limits& limits) {
    clearError(error);
    if (!validLimits(limits) || bytes.size() > limits.maximumEncodedBytes) {
        return fail(error, ErrorCode::LimitExceeded,
                    "regional opening wall coupled-state input exceeds its byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error, ErrorCode::Truncated,
                    "regional opening wall coupled-state envelope is truncated");
    }
    try {
        Reader envelope(bytes);
        for (const std::uint8_t expected : magic) {
            std::uint8_t actual = 0;
            if (!envelope.u8(actual)) {
                return fail(error, ErrorCode::Truncated,
                            "regional opening wall coupled-state magic is truncated");
            }
            if (actual != expected) {
                return fail(error, ErrorCode::InvalidMagic,
                            "regional opening wall coupled-state magic is invalid");
            }
        }
        std::uint16_t protocol = 0;
        std::uint16_t reserved = 0;
        std::size_t payloadSize = 0;
        std::uint64_t storedChecksum = 0;
        if (!envelope.u16(protocol) || !envelope.u16(reserved)
            || !envelope.count(
                payloadSize, limits.maximumEncodedBytes - envelopeBytes)
            || !envelope.u64(storedChecksum)) {
            return fail(error, readerErrorCode(envelope),
                        "regional opening wall coupled-state envelope is invalid");
        }
        if (protocol
                != sceneFluidRegionalOpeningMomentumWallCoupledStateProtocolVersion
            || reserved != 0) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "regional opening wall coupled-state protocol is unsupported");
        }
        if (payloadSize > envelope.remaining()) {
            return fail(error, ErrorCode::Truncated,
                        "regional opening wall coupled-state payload is truncated");
        }
        if (payloadSize < envelope.remaining()) {
            return fail(error, ErrorCode::TrailingData,
                        "regional opening wall coupled-state has trailing data");
        }
        const auto payload = envelope.raw(payloadSize);
        if (checksum(payload) != storedChecksum) {
            return fail(error, ErrorCode::ChecksumMismatch,
                        "regional opening wall coupled-state checksum does not match");
        }

        Reader reader(payload);
        std::uint32_t payloadStateVersion = 0;
        std::uint32_t coupledStateVersion = 0;
        std::uint64_t retainedFingerprint = 0;
        std::size_t cycleSize = 0;
        std::size_t beforeSize = 0;
        std::size_t afterSize = 0;
        const std::size_t checkpointLimit = limits.state.structureStep
            .checkpointPersistence.maximumEncodedBytes;
        if (!reader.u32(payloadStateVersion)
            || !reader.u32(coupledStateVersion)
            || !reader.u64(retainedFingerprint)
            || !reader.count(
                cycleSize, limits.cycleState.maximumEncodedBytes)
            || !reader.count(beforeSize, checkpointLimit)
            || !reader.count(afterSize, checkpointLimit)) {
            return fail(error, readerErrorCode(reader),
                        "regional opening wall coupled-state metadata is invalid");
        }
        if (payloadStateVersion != stateVersion
            || coupledStateVersion
                != sceneFluidRegionalOpeningMomentumWallCoupledStateVersion) {
            return fail(error, ErrorCode::UnsupportedVersion,
                        "regional opening wall coupled-state payload is unsupported");
        }
        const auto cycleBytes = reader.raw(cycleSize);
        const auto beforeBytes = reader.raw(beforeSize);
        const auto afterBytes = reader.raw(afterSize);
        if (reader.truncated()) {
            return fail(error, ErrorCode::Truncated,
                        "regional opening wall coupled-state nested payload is truncated");
        }
        if (reader.remaining() != 0) {
            return fail(error, ErrorCode::TrailingData,
                        "regional opening wall coupled-state payload has trailing data");
        }

        SceneFluidRegionalOpeningMomentumWallCycleState cycleState;
        if (!decodeCycleState(
                cycleBytes, transportMetric, acceptedPressureOperator,
                acceptedBasePressureOperator, grid, acceptedSweep,
                acceptedFragments, acceptedTopology, acceptedVolumeRates,
                acceptedOpeningDefinitions, acceptedOpenings,
                acceptedResistanceDefinitions, acceptedBaseMetric,
                acceptedMetric, quadrature, cycleState, error, limits)) {
            return false;
        }
        StructureCheckpoint beforeCheckpoint;
        StructureCheckpointPersistenceError structureError;
        if (!deserializeStructureCheckpoint(
                beforeBytes, structureDefinitionOwner, beforeCheckpoint,
                &structureError,
                limits.state.structureStep.checkpointPersistence)) {
            return fail(error, structureErrorCode(structureError.code),
                        "nested pre-step Structure state: "
                            + structureError.message);
        }
        Structure replay(structureDefinitionOwner.definition());
        replay.restore(beforeCheckpoint);
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner owner(
            std::move(replay));
        owner.advance(
            cycleState, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, surface, surfaceState, transfer, quadrature,
            settings, limits.state);
        auto candidate = owner.checkpoint();
        if (candidate.fingerprint != retainedFingerprint
            || !std::ranges::equal(
                candidate.structureStep.beforeStructureCheckpoint,
                beforeBytes)
            || !std::ranges::equal(
                candidate.structureStep.afterStructureCheckpoint,
                afterBytes)) {
            return fail(error, ErrorCode::ReplayMismatch,
                        "regional opening wall coupled-state replay differs");
        }
        static_assert(std::is_nothrow_move_assignable_v<
            SceneFluidRegionalOpeningMomentumWallCoupledState>);
        state = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error, ErrorCode::LimitExceeded,
                    "cannot allocate decoded regional opening wall coupled state");
    } catch (const std::length_error& exception) {
        return fail(error, ErrorCode::LimitExceeded, exception.what());
    } catch (const std::exception& exception) {
        return fail(error, ErrorCode::InvalidData, exception.what());
    }
}

} // namespace simwing::fsi
