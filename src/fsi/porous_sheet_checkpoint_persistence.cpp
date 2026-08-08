#include "porous_sheet_checkpoint_persistence.h"

#include "porous_sheet_checkpoint_detail.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'P', 'S'};
constexpr std::size_t envelopePrefixBytes = 16;
constexpr std::size_t envelopeChecksumBytes = 8;
constexpr std::uint64_t fnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

void clearError(CoupledPorousSheetCheckpointPersistenceError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

bool fail(CoupledPorousSheetCheckpointPersistenceError* error,
          const CoupledPorousSheetCheckpointPersistenceErrorCode code,
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
    const CoupledPorousSheetCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes
            >= envelopePrefixBytes + envelopeChecksumBytes
        && limits.maximumScalarSamples > 0
        && limits.maximumReplaySteps > 0;
}

bool samePublicStructure(const StructureCheckpoint& first,
                         const StructureCheckpoint& second) {
    return first.version == second.version
        && first.definitionFingerprint == second.definitionFingerprint
        && first.acceptedStepCount == second.acceptedStepCount
        && first.simulationTimeSeconds == second.simulationTimeSeconds
        && first.nodes == second.nodes
        && first.pendingExternalForcesNewtons
            == second.pendingExternalForcesNewtons
        && first.lastAppliedExternalForceNewtons
            == second.lastAppliedExternalForceNewtons;
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

} // namespace

bool serializeCoupledPorousSheetCheckpoint(
    const CoupledPorousSheetCase& owner,
    const CoupledPorousSheetCheckpoint& checkpointValue,
    std::vector<std::uint8_t>& bytes,
    CoupledPorousSheetCheckpointPersistenceError* error,
    const CoupledPorousSheetCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                    "coupled porous sheet checkpoint limits are invalid");
    }
    try {
        CoupledPorousSheetCase validator;
        validator.restore(checkpointValue);
        const auto& detail =
            CoupledPorousSheetCheckpointCodecAccess::detail(checkpointValue);
        const std::size_t sampleCount = owner.grid().cellCount();
        if (sampleCount > limits.maximumScalarSamples / 4) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet field samples exceed the configured limit");
        }
        std::vector<std::uint8_t> structureBytes;
        StructureCheckpointPersistenceError structureError;
        if (!serializeStructureCheckpoint(
                owner.structure(), detail.structure, structureBytes,
                &structureError, limits.structure)) {
            const auto code = structureError.code
                    == StructureCheckpointPersistenceErrorCode::LimitExceeded
                ? CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded
                : CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData;
            return fail(error, code,
                        "nested Structure checkpoint is invalid: "
                            + structureError.message);
        }

        Writer payload;
        payload.u32(checkpointValue.version);
        payload.u64(checkpointValue.caseFingerprint);
        payload.u64(checkpointValue.acceptedStepCount);
        payload.real(checkpointValue.simulationTimeSeconds);
        const auto counts = owner.grid().cellCounts();
        payload.u64(counts.x);
        payload.u64(counts.y);
        payload.u64(counts.z);
        const auto lower = owner.grid().lowerMeters();
        const auto upper = owner.grid().upperMeters();
        payload.real(lower.x);
        payload.real(lower.y);
        payload.real(lower.z);
        payload.real(upper.x);
        payload.real(upper.y);
        payload.real(upper.z);
        payload.u64(sampleCount);
        payload.u64(structureBytes.size());
        payload.raw(structureBytes);
        writeField(payload, detail.velocity.xFaces());
        writeField(payload, detail.velocity.yFaces());
        writeField(payload, detail.velocity.zFaces());
        writeField(payload, detail.pressure.values());
        if (payload.bytes().size()
            > limits.maximumEncodedBytes
                - envelopePrefixBytes - envelopeChecksumBytes) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet checkpoint exceeds the byte limit");
        }

        Writer envelope;
        envelope.raw(magic);
        envelope.u16(coupledPorousSheetCheckpointProtocolVersion);
        envelope.u16(0);
        envelope.u64(payload.bytes().size());
        envelope.raw(payload.bytes());
        envelope.u64(checksum(payload.bytes()));
        auto encoded = envelope.take();
        if (encoded.size() > limits.maximumEncodedBytes) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet checkpoint exceeds the byte limit");
        }
        bytes.swap(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                    "unable to allocate coupled porous sheet checkpoint bytes");
    } catch (const std::exception& exception) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                    exception.what());
    }
}

bool deserializeCoupledPorousSheetCheckpoint(
    const std::span<const std::uint8_t> bytes,
    const CoupledPorousSheetCase& owner,
    CoupledPorousSheetCheckpoint& checkpointValue,
    CoupledPorousSheetCheckpointPersistenceError* error,
    const CoupledPorousSheetCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                    "coupled porous sheet checkpoint limits are invalid");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                    "coupled porous sheet checkpoint exceeds the byte limit");
    }
    if (bytes.size() < envelopePrefixBytes + envelopeChecksumBytes) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::Truncated,
                    "coupled porous sheet checkpoint envelope is truncated");
    }
    try {
        Reader envelope(bytes);
        const auto actualMagic = envelope.raw(magic.size());
        if (!std::equal(actualMagic.begin(), actualMagic.end(), magic.begin())) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidMagic,
                        "coupled porous sheet checkpoint magic is invalid");
        }
        if (envelope.u16()
            != coupledPorousSheetCheckpointProtocolVersion) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::UnsupportedVersion,
                        "coupled porous sheet checkpoint version is unsupported");
        }
        if (envelope.u16() != 0) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                        "coupled porous sheet checkpoint reserved bits are nonzero");
        }
        const std::uint64_t payloadSize64 = envelope.u64();
        if (payloadSize64 > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet checkpoint payload exceeds platform limits");
        }
        const std::size_t payloadSize =
            static_cast<std::size_t>(payloadSize64);
        if (payloadSize > bytes.size() - envelopePrefixBytes
                - envelopeChecksumBytes) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::Truncated,
                        "coupled porous sheet checkpoint payload is truncated");
        }
        if (payloadSize != bytes.size() - envelopePrefixBytes
                - envelopeChecksumBytes) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::TrailingData,
                        "coupled porous sheet checkpoint has trailing data");
        }
        const auto payloadBytes = envelope.raw(payloadSize);
        const std::uint64_t expectedChecksum = envelope.u64();
        if (!envelope.empty()) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::TrailingData,
                        "coupled porous sheet checkpoint has trailing data");
        }
        if (checksum(payloadBytes) != expectedChecksum) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::ChecksumMismatch,
                        "coupled porous sheet checkpoint checksum does not match");
        }

        Reader payload(payloadBytes);
        const std::uint32_t stateVersion = payload.u32();
        const std::uint64_t caseFingerprint = payload.u64();
        const std::uint64_t acceptedStepCount = payload.u64();
        const double simulationTimeSeconds = payload.real();
        if (stateVersion != coupledPorousSheetCheckpointVersion
            || caseFingerprint != coupledPorousSheetCaseFingerprint
            || !std::isfinite(simulationTimeSeconds)
            || simulationTimeSeconds < 0.0) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                        "coupled porous sheet checkpoint metadata is invalid");
        }
        if (acceptedStepCount > limits.maximumReplaySteps) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet checkpoint replay exceeds the configured limit");
        }
        const std::uint64_t countX = payload.u64();
        const std::uint64_t countY = payload.u64();
        const std::uint64_t countZ = payload.u64();
        if (countX > std::numeric_limits<std::size_t>::max()
            || countY > std::numeric_limits<std::size_t>::max()
            || countZ > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet grid exceeds platform limits");
        }
        const auto expectedCounts = owner.grid().cellCounts();
        const fluid::GridCellCounts counts{
            static_cast<std::size_t>(countX),
            static_cast<std::size_t>(countY),
            static_cast<std::size_t>(countZ),
        };
        const fluid::Vector3 lower{
            payload.real(), payload.real(), payload.real()};
        const fluid::Vector3 upper{
            payload.real(), payload.real(), payload.real()};
        const std::uint64_t sampleCount64 = payload.u64();
        const std::uint64_t structureSize64 = payload.u64();
        if (counts != expectedCounts
            || lower != owner.grid().lowerMeters()
            || upper != owner.grid().upperMeters()
            || sampleCount64 != owner.grid().cellCount()) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                        "coupled porous sheet checkpoint grid is foreign");
        }
        if (sampleCount64 > limits.maximumScalarSamples / 4
            || structureSize64 > limits.maximumEncodedBytes
            || structureSize64 > std::numeric_limits<std::size_t>::max()) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                        "coupled porous sheet checkpoint allocation exceeds its limit");
        }
        const auto structureBytes = payload.raw(
            static_cast<std::size_t>(structureSize64));
        StructureCheckpoint decodedStructure;
        StructureCheckpointPersistenceError structureError;
        if (!deserializeStructureCheckpoint(
                structureBytes, owner.structure(), decodedStructure,
                &structureError, limits.structure)) {
            const auto code = structureError.code
                    == StructureCheckpointPersistenceErrorCode::LimitExceeded
                ? CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded
                : CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData;
            return fail(error, code,
                        "nested Structure checkpoint is invalid: "
                            + structureError.message);
        }
        fluid::MacVelocityField decodedVelocity(owner.grid());
        fluid::CellScalarField decodedPressure(owner.grid());
        readField(payload, decodedVelocity.xFaces());
        readField(payload, decodedVelocity.yFaces());
        readField(payload, decodedVelocity.zFaces());
        readField(payload, decodedPressure.values());
        if (!payload.empty() || !fluid::isFinite(decodedVelocity)
            || !fluid::isFinite(decodedPressure)) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                        "coupled porous sheet checkpoint fields are invalid");
        }

        CoupledPorousSheetCase replay;
        for (std::uint64_t step = 0; step < acceptedStepCount; ++step) {
            static_cast<void>(replay.advance());
        }
        const CoupledPorousSheetCheckpoint replayCheckpoint =
            replay.checkpoint();
        const auto& replayDetail =
            CoupledPorousSheetCheckpointCodecAccess::detail(
                replayCheckpoint);
        if (replayCheckpoint.version != stateVersion
            || replayCheckpoint.caseFingerprint != caseFingerprint
            || replayCheckpoint.acceptedStepCount != acceptedStepCount
            || replayCheckpoint.simulationTimeSeconds
                != simulationTimeSeconds
            || !samePublicStructure(
                replayDetail.structure, decodedStructure)
            || replayDetail.velocity != decodedVelocity
            || replayDetail.pressure != decodedPressure) {
            return fail(error,
                        CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                        "coupled porous sheet checkpoint does not match deterministic replay");
        }
        checkpointValue = replayCheckpoint;
        return true;
    } catch (const std::out_of_range&) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::Truncated,
                    "coupled porous sheet checkpoint payload is truncated");
    } catch (const std::bad_alloc&) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::LimitExceeded,
                    "unable to allocate coupled porous sheet checkpoint state");
    } catch (const std::exception& exception) {
        return fail(error,
                    CoupledPorousSheetCheckpointPersistenceErrorCode::InvalidData,
                    exception.what());
    }
}

} // namespace simwing::fsi
