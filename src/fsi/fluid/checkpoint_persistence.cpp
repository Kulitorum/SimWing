#include "fluid/checkpoint.h"

#include <algorithm>
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

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'M', 'F'};
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t interfaceFaceRecordBytes = 57;
constexpr std::size_t regionDiagnosticRecordBytes = 40;
constexpr std::size_t faceDiagnosticRecordBytes = 249;
constexpr std::size_t surfaceDiagnosticRecordBytes = 192;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(
    MovingInterfaceFluidCheckpointError* error,
    const MovingInterfaceFluidCheckpointErrorCode code,
    std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(MovingInterfaceFluidCheckpointError* error) {
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
        return u64(std::bit_cast<std::uint64_t>(value));
    }

    bool vector3(const Vector3& value) {
        return finiteDouble(value.x)
            && finiteDouble(value.y)
            && finiteDouble(value.z);
    }

    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }
    [[nodiscard]] bool invalid() const noexcept { return invalid_; }

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

    bool count(std::size_t& value, const std::uint64_t maximum) {
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
        std::uint64_t bits = 0;
        if (!u64(bits)) {
            return false;
        }
        value = std::bit_cast<double>(bits);
        if (!std::isfinite(value)) {
            invalid_ = true;
            return false;
        }
        return true;
    }

    bool vector3(Vector3& value) {
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

    bool fixedRecords(
        const std::size_t count,
        const std::size_t bytesPerRecord) {
        if (bytesPerRecord == 0
            || count > (bytes_.size() - position_) / bytesPerRecord) {
            truncated_ = true;
            return false;
        }
        return true;
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

MovingInterfaceFluidCheckpointErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return MovingInterfaceFluidCheckpointErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return MovingInterfaceFluidCheckpointErrorCode::Truncated;
    }
    return MovingInterfaceFluidCheckpointErrorCode::InvalidData;
}

bool validLimits(const MovingInterfaceFluidCheckpointLimits& limits) {
    return limits.maximumBytes >= envelopeBytes
        && limits.maximumBytes <= std::numeric_limits<std::size_t>::max()
        && limits.maximumScalarSamples > 0;
}

bool writeProjectionDiagnostics(
    Writer& writer,
    const ProjectionDiagnostics& value) {
    return writer.u8(value.converged ? 1 : 0)
        && writer.count(value.iterationCount)
        && writer.finiteDouble(value.compatibilityDivergencePerSecond)
        && writer.finiteDouble(
            value.initialResidualPascalsPerSquareMeter)
        && writer.finiteDouble(
            value.finalResidualPascalsPerSquareMeter)
        && writer.finiteDouble(value.divergenceL2BeforePerSecond)
        && writer.finiteDouble(value.divergenceL2AfterPerSecond)
        && writer.finiteDouble(value.divergenceMaximumBeforePerSecond)
        && writer.finiteDouble(value.divergenceMaximumAfterPerSecond)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.pressureMeanPascals)
        && writer.count(value.pressureJumpFaceCount)
        && writer.finiteDouble(
            value.pressureJumpSourceCompatibilityPascalsPerSquareMeter);
}

bool readProjectionDiagnostics(
    Reader& reader,
    ProjectionDiagnostics& value,
    const MovingInterfaceFluidCheckpointLimits& limits) {
    return reader.boolean(value.converged)
        && reader.count(value.iterationCount, limits.maximumScalarSamples)
        && reader.finiteDouble(value.compatibilityDivergencePerSecond)
        && reader.finiteDouble(
            value.initialResidualPascalsPerSquareMeter)
        && reader.finiteDouble(
            value.finalResidualPascalsPerSquareMeter)
        && reader.finiteDouble(value.divergenceL2BeforePerSecond)
        && reader.finiteDouble(value.divergenceL2AfterPerSecond)
        && reader.finiteDouble(value.divergenceMaximumBeforePerSecond)
        && reader.finiteDouble(value.divergenceMaximumAfterPerSecond)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.pressureMeanPascals)
        && reader.count(
            value.pressureJumpFaceCount, limits.maximumInterfaceFaces)
        && reader.finiteDouble(
            value.pressureJumpSourceCompatibilityPascalsPerSquareMeter);
}

bool writeMovingDiagnostics(
    Writer& writer,
    const MovingInterfaceProjectionDiagnostics& value) {
    if (!writeProjectionDiagnostics(writer, value.projection)
        || !writer.u32(value.interfaceVersion)
        || !writer.count(value.interfaceFaceCount)
        || !writer.count(value.fluidRegionCount)
        || !writer.finiteDouble(
            value.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond)
        || !writer.finiteDouble(
            value.maximumNormalVelocityErrorMetersPerSecond)
        || !writer.vector3(value.totalPressureForceNewtons)
        || !writer.vector3(value.totalPressureImpulseNewtonSeconds)
        || !writer.finiteDouble(value.totalPressurePowerWatts)
        || !writer.finiteDouble(value.totalPressureWorkJoules)
        || !writer.vector3(value.totalConstraintReactionForceNewtons)
        || !writer.vector3(
            value.totalConstraintReactionImpulseNewtonSeconds)
        || !writer.finiteDouble(value.totalConstraintReactionPowerWatts)
        || !writer.finiteDouble(value.totalConstraintReactionWorkJoules)
        || !writer.count(value.regions.size())) {
        return false;
    }
    for (const auto& region : value.regions) {
        if (!writer.u64(region.stableId)
            || !writer.count(region.cellCount)
            || !writer.finiteDouble(
                region.compatibilityVolumeRateCubicMetersPerSecond)
            || !writer.finiteDouble(region.pressureMeanBeforePascals)
            || !writer.finiteDouble(region.pressureMeanAfterPascals)) {
            return false;
        }
    }
    if (!writer.count(value.faces.size())) {
        return false;
    }
    for (const auto& face : value.faces) {
        if (!writer.u64(face.surfaceStableId)
            || !writer.u64(face.minusRegionStableId)
            || !writer.u64(face.plusRegionStableId)
            || !writer.u8(static_cast<std::uint8_t>(face.axis))
            || !writer.count(face.i)
            || !writer.count(face.j)
            || !writer.count(face.k)
            || !writer.vector3(face.lowerCornerMeters)
            || !writer.vector3(face.upperCornerMeters)
            || !writer.finiteDouble(face.areaSquareMeters)
            || !writer.finiteDouble(
                face.normalVelocityMetersPerSecond)
            || !writer.vector3(face.pressureTractionPascals)
            || !writer.vector3(face.pressureForceNewtons)
            || !writer.finiteDouble(face.pressurePowerWatts)
            || !writer.vector3(face.directConstraintForceNewtons)
            || !writer.vector3(
                face.constraintReactionTractionPascals)
            || !writer.vector3(face.constraintReactionForceNewtons)
            || !writer.finiteDouble(
                face.constraintReactionPowerWatts)) {
            return false;
        }
    }
    if (!writer.count(value.surfaces.size())) {
        return false;
    }
    for (const auto& surface : value.surfaces) {
        if (!writer.u64(surface.stableId)
            || !writer.count(surface.faceCount)
            || !writer.finiteDouble(surface.areaSquareMeters)
            || !writer.vector3(surface.pressureForceNewtons)
            || !writer.vector3(surface.pressureImpulseNewtonSeconds)
            || !writer.finiteDouble(
                surface.maximumPressureTractionDeviationPascals)
            || !writer.finiteDouble(surface.pressurePowerWatts)
            || !writer.finiteDouble(surface.pressureWorkJoules)
            || !writer.vector3(surface.directConstraintForceNewtons)
            || !writer.vector3(
                surface.constraintReactionForceNewtons)
            || !writer.vector3(
                surface.constraintReactionImpulseNewtonSeconds)
            || !writer.finiteDouble(
                surface.maximumConstraintReactionTractionDeviationPascals)
            || !writer.finiteDouble(
                surface.constraintReactionPowerWatts)
            || !writer.finiteDouble(
                surface.constraintReactionWorkJoules)) {
            return false;
        }
    }
    return writer.u8(value.finite ? 1 : 0);
}

bool readAxis(Reader& reader, GridFaceAxis& axis) {
    std::uint8_t encoded = 0;
    if (!reader.u8(encoded)
        || encoded > static_cast<std::uint8_t>(GridFaceAxis::Z)) {
        return false;
    }
    axis = static_cast<GridFaceAxis>(encoded);
    return true;
}

bool readMovingDiagnostics(
    Reader& reader,
    MovingInterfaceProjectionDiagnostics& value,
    const MovingInterfaceFluidCheckpointLimits& limits) {
    if (!readProjectionDiagnostics(reader, value.projection, limits)
        || !reader.u32(value.interfaceVersion)
        || !reader.count(
            value.interfaceFaceCount, limits.maximumInterfaceFaces)
        || !reader.count(
            value.fluidRegionCount, limits.maximumFluidRegions)
        || !reader.finiteDouble(
            value.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond)
        || !reader.finiteDouble(
            value.maximumNormalVelocityErrorMetersPerSecond)
        || !reader.vector3(value.totalPressureForceNewtons)
        || !reader.vector3(value.totalPressureImpulseNewtonSeconds)
        || !reader.finiteDouble(value.totalPressurePowerWatts)
        || !reader.finiteDouble(value.totalPressureWorkJoules)
        || !reader.vector3(value.totalConstraintReactionForceNewtons)
        || !reader.vector3(
            value.totalConstraintReactionImpulseNewtonSeconds)
        || !reader.finiteDouble(value.totalConstraintReactionPowerWatts)
        || !reader.finiteDouble(value.totalConstraintReactionWorkJoules)) {
        return false;
    }
    std::size_t regionCount = 0;
    if (!reader.count(regionCount, limits.maximumFluidRegions)) {
        return false;
    }
    if (!reader.fixedRecords(regionCount, regionDiagnosticRecordBytes)) {
        return false;
    }
    value.regions.resize(regionCount);
    for (auto& region : value.regions) {
        if (!reader.u64(region.stableId)
            || !reader.count(
                region.cellCount, limits.maximumScalarSamples)
            || !reader.finiteDouble(
                region.compatibilityVolumeRateCubicMetersPerSecond)
            || !reader.finiteDouble(region.pressureMeanBeforePascals)
            || !reader.finiteDouble(region.pressureMeanAfterPascals)) {
            return false;
        }
    }
    std::size_t faceCount = 0;
    if (!reader.count(faceCount, limits.maximumInterfaceFaces)) {
        return false;
    }
    if (!reader.fixedRecords(faceCount, faceDiagnosticRecordBytes)) {
        return false;
    }
    value.faces.resize(faceCount);
    for (auto& face : value.faces) {
        if (!reader.u64(face.surfaceStableId)
            || !reader.u64(face.minusRegionStableId)
            || !reader.u64(face.plusRegionStableId)
            || !readAxis(reader, face.axis)
            || !reader.count(face.i, limits.maximumScalarSamples)
            || !reader.count(face.j, limits.maximumScalarSamples)
            || !reader.count(face.k, limits.maximumScalarSamples)
            || !reader.vector3(face.lowerCornerMeters)
            || !reader.vector3(face.upperCornerMeters)
            || !reader.finiteDouble(face.areaSquareMeters)
            || !reader.finiteDouble(
                face.normalVelocityMetersPerSecond)
            || !reader.vector3(face.pressureTractionPascals)
            || !reader.vector3(face.pressureForceNewtons)
            || !reader.finiteDouble(face.pressurePowerWatts)
            || !reader.vector3(face.directConstraintForceNewtons)
            || !reader.vector3(
                face.constraintReactionTractionPascals)
            || !reader.vector3(face.constraintReactionForceNewtons)
            || !reader.finiteDouble(
                face.constraintReactionPowerWatts)) {
            return false;
        }
    }
    std::size_t surfaceCount = 0;
    if (!reader.count(
            surfaceCount, limits.maximumDiagnosticSurfaces)) {
        return false;
    }
    if (!reader.fixedRecords(
            surfaceCount, surfaceDiagnosticRecordBytes)) {
        return false;
    }
    value.surfaces.resize(surfaceCount);
    for (auto& surface : value.surfaces) {
        if (!reader.u64(surface.stableId)
            || !reader.count(
                surface.faceCount, limits.maximumInterfaceFaces)
            || !reader.finiteDouble(surface.areaSquareMeters)
            || !reader.vector3(surface.pressureForceNewtons)
            || !reader.vector3(surface.pressureImpulseNewtonSeconds)
            || !reader.finiteDouble(
                surface.maximumPressureTractionDeviationPascals)
            || !reader.finiteDouble(surface.pressurePowerWatts)
            || !reader.finiteDouble(surface.pressureWorkJoules)
            || !reader.vector3(surface.directConstraintForceNewtons)
            || !reader.vector3(
                surface.constraintReactionForceNewtons)
            || !reader.vector3(
                surface.constraintReactionImpulseNewtonSeconds)
            || !reader.finiteDouble(
                surface.maximumConstraintReactionTractionDeviationPascals)
            || !reader.finiteDouble(
                surface.constraintReactionPowerWatts)
            || !reader.finiteDouble(
                surface.constraintReactionWorkJoules)) {
            return false;
        }
    }
    return reader.boolean(value.finite);
}

template<typename Values>
bool writeField(Writer& writer, const Values& values) {
    if (!writer.count(values.size())) {
        return false;
    }
    for (const double value : values) {
        if (!writer.finiteDouble(value)) {
            return false;
        }
    }
    return true;
}

template<typename Values>
bool readField(
    Reader& reader,
    Values values,
    const std::size_t expectedCount,
    const MovingInterfaceFluidCheckpointLimits& limits) {
    std::size_t count = 0;
    if (!reader.count(count, limits.maximumScalarSamples)
        || count != expectedCount || values.size() != expectedCount) {
        return false;
    }
    for (double& value : values) {
        if (!reader.finiteDouble(value)) {
            return false;
        }
    }
    return true;
}

bool writePayload(
    Writer& writer,
    const MovingInterfaceFluidCheckpoint& checkpoint,
    const MovingInterfaceFluidState& state) {
    const auto counts = checkpoint.cellCounts;
    if (!writer.u32(checkpoint.version)
        || !writer.count(counts.x)
        || !writer.count(counts.y)
        || !writer.count(counts.z)
        || !writer.vector3(checkpoint.lowerMeters)
        || !writer.vector3(checkpoint.upperMeters)
        || !writer.count(checkpoint.scalarSampleCount)
        || !writer.u64(checkpoint.topologyFingerprint)
        || !writer.u32(state.interfaces.version())
        || !writer.count(state.interfaces.faceCount())) {
        return false;
    }
    for (const auto& face : state.interfaces.faces()) {
        if (!writer.u64(face.surfaceStableId)
            || !writer.u64(face.minusRegionStableId)
            || !writer.u64(face.plusRegionStableId)
            || !writer.u8(static_cast<std::uint8_t>(face.axis))
            || !writer.count(face.i)
            || !writer.count(face.j)
            || !writer.count(face.k)
            || !writer.finiteDouble(
                face.normalVelocityMetersPerSecond)) {
            return false;
        }
    }
    return writeField(writer, state.velocityMetersPerSecond.xFaces())
        && writeField(writer, state.velocityMetersPerSecond.yFaces())
        && writeField(writer, state.velocityMetersPerSecond.zFaces())
        && writeField(writer, state.pressurePascals.values())
        && writeMovingDiagnostics(writer, state.diagnostics);
}

} // namespace

bool serializeMovingInterfaceFluidCheckpoint(
    const MovingInterfaceFluidCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    MovingInterfaceFluidCheckpointError* error,
    const MovingInterfaceFluidCheckpointLimits& limits) {
    clearError(error);
    bytes.clear();
    if (!validLimits(limits)) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "configured fluid checkpoint limits are invalid");
    }
    if (checkpoint.scalarSampleCount > limits.maximumScalarSamples) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "fluid checkpoint exceeds the scalar-sample limit");
    }
    try {
        const PeriodicCartesianGrid grid(
            checkpoint.cellCounts,
            checkpoint.lowerMeters, checkpoint.upperMeters);
        const MovingInterfaceFluidState state =
            restoreMovingInterfaceFluidState(grid, checkpoint);
        if (state.interfaces.faceCount() > limits.maximumInterfaceFaces
            || state.interfaces.regionCount() > limits.maximumFluidRegions
            || state.diagnostics.surfaces.size()
                > limits.maximumDiagnosticSurfaces) {
            return fail(
                error,
                MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                "fluid checkpoint topology exceeds configured limits");
        }
        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload,
            static_cast<std::size_t>(limits.maximumBytes) - envelopeBytes);
        if (!writePayload(payloadWriter, checkpoint, state)) {
            return fail(
                error,
                payloadWriter.exceeded()
                    ? MovingInterfaceFluidCheckpointErrorCode::LimitExceeded
                    : MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                "fluid checkpoint payload cannot be serialized");
        }
        bytes.reserve(envelopeBytes + payload.size());
        bytes.insert(bytes.end(), magic.begin(), magic.end());
        Writer envelope(
            bytes, static_cast<std::size_t>(limits.maximumBytes));
        if (!envelope.u16(movingInterfaceFluidCheckpointProtocolVersion)
            || !envelope.u16(0)
            || !envelope.u64(static_cast<std::uint64_t>(payload.size()))
            || !envelope.u64(checksum(payload))) {
            bytes.clear();
            return fail(
                error,
                MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                "fluid checkpoint envelope exceeds the byte limit");
        }
        bytes.insert(bytes.end(), payload.begin(), payload.end());
        return true;
    } catch (const std::bad_alloc&) {
        bytes.clear();
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "unable to allocate fluid checkpoint encoding");
    } catch (const std::length_error&) {
        bytes.clear();
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "fluid checkpoint exceeds platform limits");
    } catch (const std::exception& exception) {
        bytes.clear();
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                    exception.what());
    }
}

bool deserializeMovingInterfaceFluidCheckpoint(
    const std::span<const std::uint8_t> bytes,
    MovingInterfaceFluidCheckpoint& checkpoint,
    MovingInterfaceFluidCheckpointError* error,
    const MovingInterfaceFluidCheckpointLimits& limits) {
    clearError(error);
    if (!validLimits(limits) || bytes.size() > limits.maximumBytes) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "fluid checkpoint exceeds the byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::Truncated,
                    "fluid checkpoint envelope is truncated");
    }
    if (!std::ranges::equal(magic, bytes.first(magic.size()))) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::InvalidMagic,
                    "fluid checkpoint magic is invalid");
    }
    Reader envelope(bytes.subspan(magic.size()));
    std::uint16_t protocolVersion = 0;
    std::uint16_t reserved = 0;
    std::uint64_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(protocolVersion)
        || !envelope.u16(reserved)
        || !envelope.u64(payloadSize)
        || !envelope.u64(expectedChecksum)) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::Truncated,
                    "fluid checkpoint envelope is truncated");
    }
    if (protocolVersion
        != movingInterfaceFluidCheckpointProtocolVersion) {
        return fail(
            error,
            MovingInterfaceFluidCheckpointErrorCode::UnsupportedVersion,
            "fluid checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                    "fluid checkpoint reserved bits are nonzero");
    }
    const std::uint64_t availablePayload = bytes.size() - envelopeBytes;
    if (payloadSize > availablePayload) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::Truncated,
                    "fluid checkpoint payload is truncated");
    }
    if (payloadSize < availablePayload) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::TrailingData,
                    "fluid checkpoint has trailing data");
    }
    const auto payload = bytes.subspan(
        envelopeBytes, static_cast<std::size_t>(payloadSize));
    if (checksum(payload) != expectedChecksum) {
        return fail(
            error,
            MovingInterfaceFluidCheckpointErrorCode::ChecksumMismatch,
            "fluid checkpoint checksum does not match");
    }

    try {
        Reader reader(payload);
        std::uint32_t stateVersion = 0;
        GridCellCounts counts;
        Vector3 lowerMeters;
        Vector3 upperMeters;
        std::size_t scalarSampleCount = 0;
        std::uint64_t topologyFingerprint = 0;
        std::uint32_t interfaceVersion = 0;
        std::size_t interfaceFaceCount = 0;
        if (!reader.u32(stateVersion)
            || !reader.count(counts.x, limits.maximumScalarSamples)
            || !reader.count(counts.y, limits.maximumScalarSamples)
            || !reader.count(counts.z, limits.maximumScalarSamples)
            || !reader.vector3(lowerMeters)
            || !reader.vector3(upperMeters)
            || !reader.count(
                scalarSampleCount, limits.maximumScalarSamples)
            || !reader.u64(topologyFingerprint)
            || !reader.u32(interfaceVersion)
            || !reader.count(
                interfaceFaceCount, limits.maximumInterfaceFaces)) {
            return fail(error, readerErrorCode(reader),
                        "fluid checkpoint metadata is invalid");
        }
        if (stateVersion != movingInterfaceFluidCheckpointVersion
            || interfaceVersion != faceAlignedMovingInterfaceVersion) {
            return fail(
                error,
                MovingInterfaceFluidCheckpointErrorCode::UnsupportedVersion,
                "fluid checkpoint state version is unsupported");
        }
        const PeriodicCartesianGrid grid(
            counts, lowerMeters, upperMeters);
        if (grid.cellCount() != scalarSampleCount) {
            return fail(error,
                        MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                        "fluid checkpoint grid sample count is inconsistent");
        }
        if (!reader.fixedRecords(
                interfaceFaceCount, interfaceFaceRecordBytes)) {
            return fail(error, readerErrorCode(reader),
                        "fluid checkpoint interface is truncated");
        }
        std::vector<GridFaceMovingInterface> faces(interfaceFaceCount);
        for (auto& face : faces) {
            if (!reader.u64(face.surfaceStableId)
                || !reader.u64(face.minusRegionStableId)
                || !reader.u64(face.plusRegionStableId)
                || !readAxis(reader, face.axis)
                || !reader.count(face.i, counts.x)
                || !reader.count(face.j, counts.y)
                || !reader.count(face.k, counts.z)
                || !reader.finiteDouble(
                    face.normalVelocityMetersPerSecond)) {
                return fail(error, readerErrorCode(reader),
                            "fluid checkpoint interface is invalid");
            }
        }
        FaceAlignedMovingInterface interfaces(grid, std::move(faces));
        if (interfaces.regionCount() > limits.maximumFluidRegions) {
            return fail(
                error,
                MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                "fluid checkpoint exceeds the region limit");
        }
        MacVelocityField velocity(grid);
        CellScalarField pressure(grid);
        if (!readField(
                reader, velocity.xFaces(), scalarSampleCount, limits)
            || !readField(
                reader, velocity.yFaces(), scalarSampleCount, limits)
            || !readField(
                reader, velocity.zFaces(), scalarSampleCount, limits)
            || !readField(
                reader, pressure.values(), scalarSampleCount, limits)) {
            return fail(error, readerErrorCode(reader),
                        "fluid checkpoint field payload is invalid");
        }
        MovingInterfaceProjectionDiagnostics diagnostics;
        if (!readMovingDiagnostics(reader, diagnostics, limits)) {
            return fail(error, readerErrorCode(reader),
                        "fluid checkpoint diagnostics are invalid");
        }
        if (!reader.atEnd()) {
            return fail(error,
                        MovingInterfaceFluidCheckpointErrorCode::TrailingData,
                        "fluid checkpoint payload has trailing data");
        }

        MovingInterfaceFluidCheckpoint candidate =
            checkpointMovingInterfaceFluidState(
                grid, velocity, pressure, interfaces, diagnostics);
        if (candidate.version != stateVersion
            || candidate.cellCounts != counts
            || candidate.lowerMeters != lowerMeters
            || candidate.upperMeters != upperMeters
            || candidate.scalarSampleCount != scalarSampleCount
            || candidate.topologyFingerprint != topologyFingerprint) {
            return fail(
                error,
                MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                "fluid checkpoint public metadata does not match its state");
        }
        checkpoint = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "unable to allocate fluid checkpoint state");
    } catch (const std::length_error&) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::LimitExceeded,
                    "fluid checkpoint state exceeds platform limits");
    } catch (const std::exception& exception) {
        return fail(error,
                    MovingInterfaceFluidCheckpointErrorCode::InvalidData,
                    exception.what());
    }
}

} // namespace simwing::fsi::fluid
