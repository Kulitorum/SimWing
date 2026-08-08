#include "moving_porous_flow_checkpoint_persistence.h"

#include "moving_porous_flow_checkpoint_detail.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'M', 'F'};
constexpr std::size_t envelopePrefixBytes = 16;
constexpr std::size_t envelopeChecksumBytes = 8;
constexpr std::size_t payloadMetadataBytes = 164;
constexpr std::size_t scalarBytesPerCell = 4 * sizeof(double);
constexpr std::size_t pressureJumpRecordBytes = 72;
constexpr std::uint64_t fnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

void clearError(MovingPorousFlowCaseCheckpointPersistenceError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

bool fail(
    MovingPorousFlowCaseCheckpointPersistenceError* error,
    const MovingPorousFlowCaseCheckpointPersistenceErrorCode code,
    std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

std::uint64_t checksum(const std::span<const std::uint8_t> bytes) {
    std::uint64_t result = fnvOffset;
    for (const std::uint8_t byte : bytes) {
        result ^= byte;
        result *= fnvPrime;
    }
    return result;
}

class Writer final {
public:
    void u16(const std::uint16_t value) {
        for (unsigned shift = 0; shift < 16; shift += 8) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> shift));
        }
    }

    void u32(const std::uint32_t value) {
        for (unsigned shift = 0; shift < 32; shift += 8) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> shift));
        }
    }

    void u64(const std::uint64_t value) {
        for (unsigned shift = 0; shift < 64; shift += 8) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> shift));
        }
    }

    void real(const double value) {
        u64(std::bit_cast<std::uint64_t>(value));
    }

    void raw(const std::span<const std::uint8_t> value) {
        bytes_.insert(bytes_.end(), value.begin(), value.end());
    }

    [[nodiscard]] const std::vector<std::uint8_t>& bytes() const noexcept {
        return bytes_;
    }

    [[nodiscard]] std::vector<std::uint8_t> take() {
        return std::move(bytes_);
    }

private:
    std::vector<std::uint8_t> bytes_;
};

class Reader final {
public:
    explicit Reader(const std::span<const std::uint8_t> bytes)
        : bytes_(bytes) {}

    [[nodiscard]] std::uint16_t u16() {
        require(2);
        std::uint16_t result = 0;
        for (unsigned shift = 0; shift < 16; shift += 8) {
            result |= static_cast<std::uint16_t>(bytes_[offset_++]) << shift;
        }
        return result;
    }

    [[nodiscard]] std::uint32_t u32() {
        require(4);
        std::uint32_t result = 0;
        for (unsigned shift = 0; shift < 32; shift += 8) {
            result |= static_cast<std::uint32_t>(bytes_[offset_++]) << shift;
        }
        return result;
    }

    [[nodiscard]] std::uint64_t u64() {
        require(8);
        std::uint64_t result = 0;
        for (unsigned shift = 0; shift < 64; shift += 8) {
            result |= static_cast<std::uint64_t>(bytes_[offset_++]) << shift;
        }
        return result;
    }

    [[nodiscard]] double real() {
        return std::bit_cast<double>(u64());
    }

    [[nodiscard]] std::span<const std::uint8_t> raw(
        const std::size_t count) {
        require(count);
        const auto result = bytes_.subspan(offset_, count);
        offset_ += count;
        return result;
    }

    [[nodiscard]] bool empty() const noexcept {
        return offset_ == bytes_.size();
    }

private:
    void require(const std::size_t count) const {
        if (count > bytes_.size() - offset_) {
            throw std::out_of_range("checkpoint input is truncated");
        }
    }

    std::span<const std::uint8_t> bytes_;
    std::size_t offset_ = 0;
};

bool validLimits(
    const MovingPorousFlowCaseCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes
            >= envelopePrefixBytes + envelopeChecksumBytes
        && limits.maximumScalarSamples > 0
        && limits.maximumPressureJumpFaces > 0
        && limits.maximumReplaySteps > 0;
}

void writeField(Writer& writer, const std::span<const double> values) {
    for (const double value : values) {
        writer.real(value);
    }
}

void readField(Reader& reader, const std::span<double> values) {
    for (double& value : values) {
        value = reader.real();
    }
}

void writeJump(Writer& writer, const fluid::GridFacePressureJump& jump) {
    writer.u64(jump.surfaceStableId);
    writer.u64(jump.minusRegionStableId);
    writer.u64(jump.plusRegionStableId);
    writer.u32(static_cast<std::uint32_t>(jump.axis));
    writer.u32(0);
    writer.u64(jump.i);
    writer.u64(jump.j);
    writer.u64(jump.k);
    writer.real(jump.pressureJumpPascals);
    writer.real(jump.crossingFraction);
}

bool samePublicCheckpoint(
    const MovingPorousFlowCaseCheckpoint& first,
    const MovingPorousFlowCaseCheckpoint& second) {
    return first.version == second.version
        && first.caseDefinitionFingerprint
            == second.caseDefinitionFingerprint
        && first.cellCounts == second.cellCounts
        && first.lowerMeters == second.lowerMeters
        && first.upperMeters == second.upperMeters
        && first.scalarSampleCount == second.scalarSampleCount
        && first.pressureJumpCount == second.pressureJumpCount
        && first.acceptedStepCount == second.acceptedStepCount
        && first.simulationTimeSeconds == second.simulationTimeSeconds
        && first.sheetPositionMeters == second.sheetPositionMeters
        && first.sheetVelocityMetersPerSecond
            == second.sheetVelocityMetersPerSecond
        && first.topologyRebaseCount == second.topologyRebaseCount
        && first.porousTopology == second.porousTopology;
}

} // namespace

bool serializeMovingPorousFlowCaseCheckpoint(
    const MovingPorousFlowCaseCheckpoint& checkpointValue,
    std::vector<std::uint8_t>& bytes,
    MovingPorousFlowCaseCheckpointPersistenceError* error,
    const MovingPorousFlowCaseCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                    "moving porous-flow checkpoint limits are invalid");
    }
    try {
        MovingPorousFlowCase validator;
        validator.restore(checkpointValue);
        const auto& detail =
            MovingPorousFlowCaseCheckpointCodecAccess::detail(
                checkpointValue);
        if (checkpointValue.scalarSampleCount
                > limits.maximumScalarSamples / 4
            || checkpointValue.pressureJumpCount
                > limits.maximumPressureJumpFaces
            || checkpointValue.acceptedStepCount
                > limits.maximumReplaySteps) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint exceeds a configured record limit");
        }

        Writer payload;
        payload.u32(checkpointValue.version);
        payload.u64(checkpointValue.caseDefinitionFingerprint);
        payload.u64(checkpointValue.cellCounts.x);
        payload.u64(checkpointValue.cellCounts.y);
        payload.u64(checkpointValue.cellCounts.z);
        payload.real(checkpointValue.lowerMeters.x);
        payload.real(checkpointValue.lowerMeters.y);
        payload.real(checkpointValue.lowerMeters.z);
        payload.real(checkpointValue.upperMeters.x);
        payload.real(checkpointValue.upperMeters.y);
        payload.real(checkpointValue.upperMeters.z);
        payload.u64(checkpointValue.scalarSampleCount);
        payload.u64(checkpointValue.pressureJumpCount);
        payload.u64(checkpointValue.acceptedStepCount);
        payload.real(checkpointValue.simulationTimeSeconds);
        payload.real(checkpointValue.sheetPositionMeters);
        payload.real(checkpointValue.sheetVelocityMetersPerSecond);
        payload.u64(checkpointValue.topologyRebaseCount);
        payload.u32(checkpointValue.porousTopology.version);
        payload.u32(static_cast<std::uint32_t>(
            checkpointValue.porousTopology.axis));
        payload.u64(checkpointValue.porousTopology.faceCoordinate);
        payload.u64(std::bit_cast<std::uint64_t>(
            checkpointValue.porousTopology.periodicImage));
        writeField(payload, detail.velocityMetersPerSecond.xFaces());
        writeField(payload, detail.velocityMetersPerSecond.yFaces());
        writeField(payload, detail.velocityMetersPerSecond.zFaces());
        writeField(payload, detail.pressurePascals.values());
        for (const auto& jump : detail.pressureJumps.faces()) {
            writeJump(payload, jump);
        }
        if (payload.bytes().size()
            > limits.maximumEncodedBytes
                - envelopePrefixBytes - envelopeChecksumBytes) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint exceeds the byte limit");
        }

        Writer envelope;
        envelope.raw(magic);
        envelope.u16(movingPorousFlowCaseCheckpointProtocolVersion);
        envelope.u16(0);
        envelope.u64(payload.bytes().size());
        envelope.raw(payload.bytes());
        envelope.u64(checksum(payload.bytes()));
        auto encoded = envelope.take();
        if (encoded.size() > limits.maximumEncodedBytes) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint exceeds the byte limit");
        }
        bytes.swap(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                    "unable to allocate moving porous-flow checkpoint bytes");
    } catch (const std::exception& exception) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                    exception.what());
    }
}

bool deserializeMovingPorousFlowCaseCheckpoint(
    const std::span<const std::uint8_t> bytes,
    MovingPorousFlowCaseCheckpoint& checkpointValue,
    MovingPorousFlowCaseCheckpointPersistenceError* error,
    const MovingPorousFlowCaseCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                    "moving porous-flow checkpoint limits are invalid");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                    "moving porous-flow checkpoint exceeds the byte limit");
    }
    if (bytes.size() < envelopePrefixBytes + envelopeChecksumBytes) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::Truncated,
                    "moving porous-flow checkpoint envelope is truncated");
    }
    try {
        Reader envelope(bytes);
        const auto actualMagic = envelope.raw(magic.size());
        if (!std::equal(actualMagic.begin(), actualMagic.end(),
                        magic.begin())) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidMagic,
                        "moving porous-flow checkpoint magic is invalid");
        }
        if (envelope.u16()
            != movingPorousFlowCaseCheckpointProtocolVersion) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::UnsupportedVersion,
                        "moving porous-flow checkpoint version is unsupported");
        }
        if (envelope.u16() != 0) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                        "moving porous-flow checkpoint reserved bits are nonzero");
        }
        const std::uint64_t payloadSize64 = envelope.u64();
        if (payloadSize64 > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint payload exceeds platform limits");
        }
        const std::size_t payloadSize =
            static_cast<std::size_t>(payloadSize64);
        if (payloadSize > bytes.size() - envelopePrefixBytes
                - envelopeChecksumBytes) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::Truncated,
                        "moving porous-flow checkpoint payload is truncated");
        }
        if (payloadSize != bytes.size() - envelopePrefixBytes
                - envelopeChecksumBytes) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::TrailingData,
                        "moving porous-flow checkpoint has trailing data");
        }
        const auto payloadBytes = envelope.raw(payloadSize);
        const std::uint64_t expectedChecksum = envelope.u64();
        if (!envelope.empty()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::TrailingData,
                        "moving porous-flow checkpoint has trailing data");
        }
        if (checksum(payloadBytes) != expectedChecksum) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::ChecksumMismatch,
                        "moving porous-flow checkpoint checksum does not match");
        }

        MovingPorousFlowCase owner;
        MovingPorousFlowCaseCheckpoint metadata;
        Reader payload(payloadBytes);
        metadata.version = payload.u32();
        metadata.caseDefinitionFingerprint = payload.u64();
        const std::uint64_t countX = payload.u64();
        const std::uint64_t countY = payload.u64();
        const std::uint64_t countZ = payload.u64();
        if (countX > std::numeric_limits<std::size_t>::max()
            || countY > std::numeric_limits<std::size_t>::max()
            || countZ > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint grid exceeds platform limits");
        }
        metadata.cellCounts = {
            static_cast<std::size_t>(countX),
            static_cast<std::size_t>(countY),
            static_cast<std::size_t>(countZ),
        };
        metadata.lowerMeters = {
            payload.real(), payload.real(), payload.real()};
        metadata.upperMeters = {
            payload.real(), payload.real(), payload.real()};
        const std::uint64_t scalarSampleCount64 = payload.u64();
        const std::uint64_t pressureJumpCount64 = payload.u64();
        metadata.acceptedStepCount = payload.u64();
        metadata.simulationTimeSeconds = payload.real();
        metadata.sheetPositionMeters = payload.real();
        metadata.sheetVelocityMetersPerSecond = payload.real();
        metadata.topologyRebaseCount = payload.u64();
        const std::uint32_t topologyVersion = payload.u32();
        const std::uint32_t topologyAxis = payload.u32();
        const std::uint64_t topologyFace64 = payload.u64();
        const std::int64_t topologyImage =
            std::bit_cast<std::int64_t>(payload.u64());
        if (scalarSampleCount64 > std::numeric_limits<std::size_t>::max()
            || pressureJumpCount64
                > std::numeric_limits<std::size_t>::max()
            || topologyFace64 > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint metadata exceeds platform limits");
        }
        metadata.scalarSampleCount =
            static_cast<std::size_t>(scalarSampleCount64);
        metadata.pressureJumpCount =
            static_cast<std::size_t>(pressureJumpCount64);
        if (topologyAxis
            > static_cast<std::uint32_t>(fluid::GridFaceAxis::Z)) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                        "moving porous-flow checkpoint topology axis is invalid");
        }
        metadata.porousTopology = {
            topologyVersion,
            static_cast<fluid::GridFaceAxis>(topologyAxis),
            static_cast<std::size_t>(topologyFace64),
            topologyImage,
        };
        if (metadata.version != movingPorousFlowCaseCheckpointVersion) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::UnsupportedVersion,
                        "moving porous-flow checkpoint state version is unsupported");
        }
        if (metadata.caseDefinitionFingerprint
                != movingPorousFlowCaseDefinitionFingerprint
            || metadata.cellCounts != owner.grid().cellCounts()
            || metadata.lowerMeters != owner.grid().lowerMeters()
            || metadata.upperMeters != owner.grid().upperMeters()
            || metadata.scalarSampleCount != owner.grid().cellCount()
            || !std::isfinite(metadata.simulationTimeSeconds)
            || !std::isfinite(metadata.sheetPositionMeters)
            || !std::isfinite(metadata.sheetVelocityMetersPerSecond)) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                        "moving porous-flow checkpoint metadata is foreign or invalid");
        }
        if (metadata.scalarSampleCount
                > limits.maximumScalarSamples / 4
            || metadata.pressureJumpCount
                > limits.maximumPressureJumpFaces
            || metadata.acceptedStepCount
                > limits.maximumReplaySteps) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint exceeds a configured record limit");
        }
        if (metadata.scalarSampleCount
                > (std::numeric_limits<std::size_t>::max()
                   - payloadMetadataBytes) / scalarBytesPerCell) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint field size exceeds platform limits");
        }
        std::size_t expectedPayloadBytes = payloadMetadataBytes
            + metadata.scalarSampleCount * scalarBytesPerCell;
        if (metadata.pressureJumpCount
            > (std::numeric_limits<std::size_t>::max()
               - expectedPayloadBytes) / pressureJumpRecordBytes) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                        "moving porous-flow checkpoint crossing size exceeds platform limits");
        }
        expectedPayloadBytes += metadata.pressureJumpCount
            * pressureJumpRecordBytes;
        if (expectedPayloadBytes > payloadBytes.size()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::Truncated,
                        "moving porous-flow checkpoint field payload is truncated");
        }
        if (expectedPayloadBytes < payloadBytes.size()) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::TrailingData,
                        "moving porous-flow checkpoint field payload has trailing data");
        }

        fluid::MacVelocityField decodedVelocity(owner.grid());
        fluid::CellScalarField decodedPressure(owner.grid());
        readField(payload, decodedVelocity.xFaces());
        readField(payload, decodedVelocity.yFaces());
        readField(payload, decodedVelocity.zFaces());
        readField(payload, decodedPressure.values());
        std::vector<fluid::GridFacePressureJump> decodedFaces;
        decodedFaces.reserve(metadata.pressureJumpCount);
        for (std::size_t index = 0;
             index < metadata.pressureJumpCount; ++index) {
            fluid::GridFacePressureJump jump;
            jump.surfaceStableId = payload.u64();
            jump.minusRegionStableId = payload.u64();
            jump.plusRegionStableId = payload.u64();
            const std::uint32_t axis = payload.u32();
            const std::uint32_t reserved = payload.u32();
            const std::uint64_t i = payload.u64();
            const std::uint64_t j = payload.u64();
            const std::uint64_t k = payload.u64();
            jump.pressureJumpPascals = payload.real();
            jump.crossingFraction = payload.real();
            if (axis > static_cast<std::uint32_t>(fluid::GridFaceAxis::Z)
                || reserved != 0
                || i > std::numeric_limits<std::size_t>::max()
                || j > std::numeric_limits<std::size_t>::max()
                || k > std::numeric_limits<std::size_t>::max()) {
                return fail(error,
                            MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                            "moving porous-flow checkpoint crossing is invalid");
            }
            jump.axis = static_cast<fluid::GridFaceAxis>(axis);
            jump.i = static_cast<std::size_t>(i);
            jump.j = static_cast<std::size_t>(j);
            jump.k = static_cast<std::size_t>(k);
            decodedFaces.push_back(jump);
        }
        if (!payload.empty() || !fluid::isFinite(decodedVelocity)
            || !fluid::isFinite(decodedPressure)) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                        "moving porous-flow checkpoint fields are invalid");
        }
        const fluid::SharpPressureJumpField decodedJumps(
            owner.grid(), std::move(decodedFaces));

        for (std::uint64_t step = 0;
             step < metadata.acceptedStepCount; ++step) {
            static_cast<void>(owner.advance());
        }
        const auto replayCheckpoint = owner.checkpoint();
        const auto& replayDetail =
            MovingPorousFlowCaseCheckpointCodecAccess::detail(
                replayCheckpoint);
        if (!samePublicCheckpoint(metadata, replayCheckpoint)
            || replayDetail.velocityMetersPerSecond != decodedVelocity
            || replayDetail.pressurePascals != decodedPressure
            || replayDetail.pressureJumps != decodedJumps) {
            return fail(error,
                        MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                        "moving porous-flow checkpoint does not match deterministic replay");
        }
        checkpointValue = replayCheckpoint;
        return true;
    } catch (const std::out_of_range&) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::Truncated,
                    "moving porous-flow checkpoint payload is truncated");
    } catch (const std::bad_alloc&) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::LimitExceeded,
                    "unable to allocate moving porous-flow checkpoint state");
    } catch (const std::exception& exception) {
        return fail(error,
                    MovingPorousFlowCaseCheckpointPersistenceErrorCode::InvalidData,
                    exception.what());
    }
}

} // namespace simwing::fsi
