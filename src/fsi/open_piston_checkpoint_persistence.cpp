#include "open_piston_checkpoint_persistence.h"

#include "open_piston_checkpoint_detail.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {

struct OpenPistonCheckpointPersistenceAccess {
    [[nodiscard]] static const OpenPistonCaseCheckpoint::Detail* detail(
        const OpenPistonCaseCheckpoint& checkpoint) noexcept {
        return checkpoint.detail.get();
    }
    static void setDetail(
        OpenPistonCaseCheckpoint& checkpoint,
        std::shared_ptr<const OpenPistonCaseCheckpoint::Detail> detail) {
        checkpoint.detail = std::move(detail);
    }
};

namespace {

constexpr std::array<std::uint8_t, 4> magic{'S', 'W', 'O', 'P'};
constexpr std::uint32_t stateVersion = 1;
constexpr std::size_t envelopeBytes = 24;
constexpr std::size_t cutFaceRecordBytes = 217;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(OpenPistonCheckpointPersistenceError* error,
          const OpenPistonCheckpointPersistenceErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(OpenPistonCheckpointPersistenceError* error) {
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
        if (!std::isfinite(value)) return false;
        return u64(std::bit_cast<std::uint64_t>(value));
    }
    bool vector3(const StructureVector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool vector3(const fluid::Vector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
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
    bool vector3(StructureVector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool vector3(fluid::Vector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded) || encoded > 1) return false;
        value = encoded != 0;
        return true;
    }
    bool fixedRecords(const std::size_t count,
                      const std::size_t bytesPerRecord) {
        if (count > (bytes_.size() - position_) / bytesPerRecord) {
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

OpenPistonCheckpointPersistenceErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return OpenPistonCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return OpenPistonCheckpointPersistenceErrorCode::Truncated;
    }
    return OpenPistonCheckpointPersistenceErrorCode::InvalidData;
}

bool writeAxis(Writer& writer, const fluid::GridFaceAxis value) {
    return writer.u8(static_cast<std::uint8_t>(value));
}

bool readAxis(Reader& reader, fluid::GridFaceAxis& value) {
    std::uint8_t encoded = 0;
    if (!reader.u8(encoded)
        || encoded > static_cast<std::uint8_t>(fluid::GridFaceAxis::Z)) {
        return false;
    }
    value = static_cast<fluid::GridFaceAxis>(encoded);
    return true;
}

bool writeControlVolume(
    Writer& writer,
    const fluid::PlanarControlVolumeDiagnostics& value) {
    return writer.u32(value.version)
        && writer.u64(value.movingSurfaceStableId)
        && writer.u64(value.fluidRegionStableId)
        && writeAxis(writer, value.axis)
        && writer.count(value.movingPlaneCoordinate)
        && writer.count(value.openingPlaneCoordinate)
        && writer.count(value.movingSurfaceFaceCount)
        && writer.count(value.openingFaceCount)
        && writer.finiteDouble(value.crossSectionAreaSquareMeters)
        && writer.finiteDouble(value.referenceVolumeCubicMeters)
        && writer.finiteDouble(value.startVolumeCubicMeters)
        && writer.finiteDouble(value.endVolumeCubicMeters)
        && writer.finiteDouble(value.startCutCellVolumeCubicMeters)
        && writer.finiteDouble(value.endCutCellVolumeCubicMeters)
        && writer.finiteDouble(value.startCutCellVolumeFraction)
        && writer.finiteDouble(value.endCutCellVolumeFraction)
        && writer.finiteDouble(value.geometryVolumeChangeCubicMeters)
        && writer.finiteDouble(value.surfaceSweptVolumeCubicMeters)
        && writer.finiteDouble(value.openingTransportVolumeCubicMeters)
        && writer.finiteDouble(value.surfaceGeometryResidualCubicMeters)
        && writer.finiteDouble(value.continuityResidualCubicMeters)
        && writer.finiteDouble(
            value.maximumSurfaceVelocityErrorMetersPerSecond)
        && writer.finiteDouble(value.surfacePressurePowerWatts)
        && writer.finiteDouble(value.rectangularSurfacePressureWorkJoules)
        && writer.u8(value.finite ? 1 : 0)
        && writer.u8(value.accepted ? 1 : 0);
}

bool readControlVolume(Reader& reader,
                       fluid::PlanarControlVolumeDiagnostics& value,
                       const std::size_t maximum) {
    return reader.u32(value.version)
        && value.version == fluid::planarMovingControlVolumeVersion
        && reader.u64(value.movingSurfaceStableId)
        && reader.u64(value.fluidRegionStableId)
        && readAxis(reader, value.axis)
        && reader.count(value.movingPlaneCoordinate, maximum)
        && reader.count(value.openingPlaneCoordinate, maximum)
        && reader.count(value.movingSurfaceFaceCount, maximum)
        && reader.count(value.openingFaceCount, maximum)
        && reader.finiteDouble(value.crossSectionAreaSquareMeters)
        && reader.finiteDouble(value.referenceVolumeCubicMeters)
        && reader.finiteDouble(value.startVolumeCubicMeters)
        && reader.finiteDouble(value.endVolumeCubicMeters)
        && reader.finiteDouble(value.startCutCellVolumeCubicMeters)
        && reader.finiteDouble(value.endCutCellVolumeCubicMeters)
        && reader.finiteDouble(value.startCutCellVolumeFraction)
        && reader.finiteDouble(value.endCutCellVolumeFraction)
        && reader.finiteDouble(value.geometryVolumeChangeCubicMeters)
        && reader.finiteDouble(value.surfaceSweptVolumeCubicMeters)
        && reader.finiteDouble(value.openingTransportVolumeCubicMeters)
        && reader.finiteDouble(value.surfaceGeometryResidualCubicMeters)
        && reader.finiteDouble(value.continuityResidualCubicMeters)
        && reader.finiteDouble(
            value.maximumSurfaceVelocityErrorMetersPerSecond)
        && reader.finiteDouble(value.surfacePressurePowerWatts)
        && reader.finiteDouble(value.rectangularSurfacePressureWorkJoules)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeRebase(
    Writer& writer,
    const fluid::PlanarControlVolumeRebaseDiagnostics& value) {
    return writer.u32(value.version)
        && writer.u64(value.movingSurfaceStableId)
        && writer.u64(value.fluidRegionStableId)
        && writeAxis(writer, value.axis)
        && writer.count(value.previousMovingPlaneCoordinate)
        && writer.count(value.rebasedMovingPlaneCoordinate)
        && writer.count(value.openingPlaneCoordinate)
        && writer.finiteDouble(value.completedCellOffsetMeters)
        && writer.finiteDouble(value.previousTerminalVolumeCubicMeters)
        && writer.finiteDouble(value.rebasedReferenceVolumeCubicMeters)
        && writer.finiteDouble(value.volumeContinuityResidualCubicMeters)
        && writer.u8(value.finite ? 1 : 0)
        && writer.u8(value.accepted ? 1 : 0);
}

bool readRebase(Reader& reader,
                fluid::PlanarControlVolumeRebaseDiagnostics& value,
                const std::size_t maximum) {
    return reader.u32(value.version)
        && value.version == fluid::planarControlVolumeRebaseVersion
        && reader.u64(value.movingSurfaceStableId)
        && reader.u64(value.fluidRegionStableId)
        && readAxis(reader, value.axis)
        && reader.count(value.previousMovingPlaneCoordinate, maximum)
        && reader.count(value.rebasedMovingPlaneCoordinate, maximum)
        && reader.count(value.openingPlaneCoordinate, maximum)
        && reader.finiteDouble(value.completedCellOffsetMeters)
        && reader.finiteDouble(value.previousTerminalVolumeCubicMeters)
        && reader.finiteDouble(value.rebasedReferenceVolumeCubicMeters)
        && reader.finiteDouble(value.volumeContinuityResidualCubicMeters)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeCutFace(
    Writer& writer,
    const fluid::PlanarCutSurfacePressureFaceDiagnostics& value) {
    return writer.u64(value.surfaceStableId)
        && writer.u64(value.minusRegionStableId)
        && writer.u64(value.plusRegionStableId)
        && writeAxis(writer, value.axis)
        && writer.count(value.i) && writer.count(value.j)
        && writer.count(value.k)
        && writer.vector3(value.gridLowerCornerMeters)
        && writer.vector3(value.gridUpperCornerMeters)
        && writer.vector3(value.physicalLowerCornerMeters)
        && writer.vector3(value.physicalUpperCornerMeters)
        && writer.finiteDouble(value.areaSquareMeters)
        && writer.finiteDouble(value.normalVelocityMetersPerSecond)
        && writer.vector3(value.pressureTractionPascals)
        && writer.vector3(value.pressureForceNewtons)
        && writer.finiteDouble(value.pressurePowerWatts);
}

bool readCutFace(
    Reader& reader,
    fluid::PlanarCutSurfacePressureFaceDiagnostics& value,
    const std::size_t maximum) {
    return reader.u64(value.surfaceStableId)
        && reader.u64(value.minusRegionStableId)
        && reader.u64(value.plusRegionStableId)
        && readAxis(reader, value.axis)
        && reader.count(value.i, maximum) && reader.count(value.j, maximum)
        && reader.count(value.k, maximum)
        && reader.vector3(value.gridLowerCornerMeters)
        && reader.vector3(value.gridUpperCornerMeters)
        && reader.vector3(value.physicalLowerCornerMeters)
        && reader.vector3(value.physicalUpperCornerMeters)
        && reader.finiteDouble(value.areaSquareMeters)
        && reader.finiteDouble(value.normalVelocityMetersPerSecond)
        && reader.vector3(value.pressureTractionPascals)
        && reader.vector3(value.pressureForceNewtons)
        && reader.finiteDouble(value.pressurePowerWatts);
}

bool writeCutSurface(
    Writer& writer,
    const fluid::PlanarCutSurfacePressureDiagnostics& value,
    const OpenPistonCheckpointPersistenceLimits& limits) {
    if (value.faces.size() > limits.maximumCutSurfaceFaces) return false;
    if (!writer.u32(value.version)
        || !writer.u32(value.sourceInterfaceVersion)
        || !writer.u64(value.surfaceStableId)
        || !writer.u64(value.fluidRegionStableId)
        || !writeAxis(writer, value.axis)
        || !writer.count(value.movingPlaneCoordinate)
        || !writer.count(value.faceCount)
        || !writer.vector3(value.momentReferenceMeters)
        || !writer.finiteDouble(value.surfaceOffsetMeters)
        || !writer.finiteDouble(value.gridPlaneCoordinateMeters)
        || !writer.finiteDouble(value.physicalPlaneCoordinateMeters)
        || !writer.finiteDouble(value.periodicPositionResidualMeters)
        || !writer.finiteDouble(value.normalVelocityMetersPerSecond)
        || !writer.finiteDouble(
            value.maximumNormalVelocitySpreadMetersPerSecond)
        || !writer.u8(value.kinematicsResampled ? 1 : 0)
        || !writer.finiteDouble(
            value.reactionSourcePhysicalPlaneCoordinateMeters)
        || !writer.finiteDouble(
            value.reactionSourceNormalVelocityMetersPerSecond)
        || !writer.finiteDouble(value.areaSquareMeters)
        || !writer.finiteDouble(value.sourceAreaSquareMeters)
        || !writer.finiteDouble(value.areaResidualSquareMeters)
        || !writer.vector3(value.pressureForceNewtons)
        || !writer.vector3(value.sourcePressureForceNewtons)
        || !writer.vector3(value.forceResidualNewtons)
        || !writer.finiteDouble(value.forceResidualNormNewtons)
        || !writer.vector3(value.pressureMomentNewtonMeters)
        || !writer.finiteDouble(value.pressurePowerWatts)
        || !writer.finiteDouble(value.sourcePressurePowerWatts)
        || !writer.finiteDouble(value.powerResidualWatts)
        || !writer.count(value.faces.size())) {
        return false;
    }
    for (const auto& face : value.faces) {
        if (!writeCutFace(writer, face)) return false;
    }
    return writer.u8(value.finite ? 1 : 0)
        && writer.u8(value.accepted ? 1 : 0);
}

bool readCutSurface(
    Reader& reader,
    fluid::PlanarCutSurfacePressureDiagnostics& value,
    const OpenPistonCheckpointPersistenceLimits& limits) {
    std::size_t count = 0;
    if (!reader.u32(value.version)
        || value.version != fluid::planarCutSurfacePressureVersion
        || !reader.u32(value.sourceInterfaceVersion)
        || !reader.u64(value.surfaceStableId)
        || !reader.u64(value.fluidRegionStableId)
        || !readAxis(reader, value.axis)
        || !reader.count(value.movingPlaneCoordinate,
                         limits.maximumCutSurfaceFaces)
        || !reader.count(value.faceCount,
                         limits.maximumCutSurfaceFaces)
        || !reader.vector3(value.momentReferenceMeters)
        || !reader.finiteDouble(value.surfaceOffsetMeters)
        || !reader.finiteDouble(value.gridPlaneCoordinateMeters)
        || !reader.finiteDouble(value.physicalPlaneCoordinateMeters)
        || !reader.finiteDouble(value.periodicPositionResidualMeters)
        || !reader.finiteDouble(value.normalVelocityMetersPerSecond)
        || !reader.finiteDouble(
            value.maximumNormalVelocitySpreadMetersPerSecond)
        || !reader.boolean(value.kinematicsResampled)
        || !reader.finiteDouble(
            value.reactionSourcePhysicalPlaneCoordinateMeters)
        || !reader.finiteDouble(
            value.reactionSourceNormalVelocityMetersPerSecond)
        || !reader.finiteDouble(value.areaSquareMeters)
        || !reader.finiteDouble(value.sourceAreaSquareMeters)
        || !reader.finiteDouble(value.areaResidualSquareMeters)
        || !reader.vector3(value.pressureForceNewtons)
        || !reader.vector3(value.sourcePressureForceNewtons)
        || !reader.vector3(value.forceResidualNewtons)
        || !reader.finiteDouble(value.forceResidualNormNewtons)
        || !reader.vector3(value.pressureMomentNewtonMeters)
        || !reader.finiteDouble(value.pressurePowerWatts)
        || !reader.finiteDouble(value.sourcePressurePowerWatts)
        || !reader.finiteDouble(value.powerResidualWatts)
        || !reader.count(count, limits.maximumCutSurfaceFaces)
        || !reader.fixedRecords(count, cutFaceRecordBytes)) {
        return false;
    }
    value.faces.resize(count);
    for (auto& face : value.faces) {
        if (!readCutFace(reader, face, limits.maximumCutSurfaceFaces)) {
            return false;
        }
    }
    return reader.boolean(value.finite) && reader.boolean(value.accepted)
        && value.faceCount == value.faces.size();
}

bool writeBridge(Writer& writer,
                 const PlanarFaceResolvedBridgeDiagnostics& value) {
    return writer.u32(value.version)
        && writer.u64(value.fluidSurfaceStableId)
        && writer.count(value.fluidFaceCount)
        && writer.count(value.structureTriangleCount)
        && writer.count(value.overlapPatchCount)
        && writer.finiteDouble(value.fluidAreaSquareMeters)
        && writer.finiteDouble(value.referenceStructureAreaSquareMeters)
        && writer.finiteDouble(value.areaResidualSquareMeters)
        && writer.vector3(value.fluidPressureForceNewtons)
        && writer.vector3(value.structureSurfaceForceNewtons)
        && writer.vector3(value.forceResidualNewtons)
        && writer.finiteDouble(value.forceResidualNormNewtons)
        && writer.vector3(value.fluidPressureMomentNewtonMeters)
        && writer.vector3(value.structureSurfaceMomentNewtonMeters)
        && writer.vector3(value.momentResidualNewtonMeters)
        && writer.finiteDouble(value.momentResidualNormNewtonMeters)
        && writer.finiteDouble(value.fluidPressurePowerWatts)
        && writer.finiteDouble(value.structureSurfacePowerWatts)
        && writer.finiteDouble(value.powerResidualWatts)
        && writer.finiteDouble(value.maximumFacePowerResidualWatts)
        && writer.u8(static_cast<std::uint8_t>(value.correspondenceMode))
        && writer.finiteDouble(value.gridPlaneCoordinateMeters)
        && writer.finiteDouble(value.physicalPlaneCoordinateMeters)
        && writer.finiteDouble(value.normalTranslationFromReferenceMeters)
        && writer.finiteDouble(value.maximumRigidPositionResidualMeters)
        && writer.finiteDouble(
            value.maximumRigidVelocityResidualMetersPerSecond)
        && writer.u8(value.finite ? 1 : 0);
}

bool readBridge(Reader& reader,
                PlanarFaceResolvedBridgeDiagnostics& value,
                const std::size_t maximum) {
    std::uint8_t mode = 0;
    if (!reader.u32(value.version)
        || value.version != planarFaceResolvedBridgeVersion
        || !reader.u64(value.fluidSurfaceStableId)
        || !reader.count(value.fluidFaceCount, maximum)
        || !reader.count(value.structureTriangleCount, maximum)
        || !reader.count(value.overlapPatchCount, maximum)
        || !reader.finiteDouble(value.fluidAreaSquareMeters)
        || !reader.finiteDouble(value.referenceStructureAreaSquareMeters)
        || !reader.finiteDouble(value.areaResidualSquareMeters)
        || !reader.vector3(value.fluidPressureForceNewtons)
        || !reader.vector3(value.structureSurfaceForceNewtons)
        || !reader.vector3(value.forceResidualNewtons)
        || !reader.finiteDouble(value.forceResidualNormNewtons)
        || !reader.vector3(value.fluidPressureMomentNewtonMeters)
        || !reader.vector3(value.structureSurfaceMomentNewtonMeters)
        || !reader.vector3(value.momentResidualNewtonMeters)
        || !reader.finiteDouble(value.momentResidualNormNewtonMeters)
        || !reader.finiteDouble(value.fluidPressurePowerWatts)
        || !reader.finiteDouble(value.structureSurfacePowerWatts)
        || !reader.finiteDouble(value.powerResidualWatts)
        || !reader.finiteDouble(value.maximumFacePowerResidualWatts)
        || !reader.u8(mode)
        || mode > static_cast<std::uint8_t>(
            PlanarFaceCorrespondenceMode::RigidNormalTranslation)
        || !reader.finiteDouble(value.gridPlaneCoordinateMeters)
        || !reader.finiteDouble(value.physicalPlaneCoordinateMeters)
        || !reader.finiteDouble(value.normalTranslationFromReferenceMeters)
        || !reader.finiteDouble(value.maximumRigidPositionResidualMeters)
        || !reader.finiteDouble(
            value.maximumRigidVelocityResidualMetersPerSecond)
        || !reader.boolean(value.finite)) {
        return false;
    }
    value.correspondenceMode =
        static_cast<PlanarFaceCorrespondenceMode>(mode);
    return true;
}

bool writeConservation(Writer& writer,
                       const OpenPistonConservationDiagnostics& value) {
    return writer.u32(value.version)
        && writer.vector3(value.structureMomentumChangeNewtonSeconds)
        && writer.vector3(value.fluidMomentumChangeNewtonSeconds)
        && writer.vector3(value.pressureImpulseToStructureNewtonSeconds)
        && writer.vector3(value.actuatorImpulseNewtonSeconds)
        && writer.vector3(value.structureMomentumResidualNewtonSeconds)
        && writer.finiteDouble(
            value.structureMomentumResidualNormNewtonSeconds)
        && writer.vector3(value.fluidMomentumResidualNewtonSeconds)
        && writer.finiteDouble(
            value.fluidMomentumResidualNormNewtonSeconds)
        && writer.vector3(value.systemMomentumResidualNewtonSeconds)
        && writer.finiteDouble(
            value.systemMomentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.structureKineticEnergyChangeJoules)
        && writer.finiteDouble(value.fluidKineticEnergyChangeJoules)
        && writer.finiteDouble(value.pressureWorkToStructureJoules)
        && writer.finiteDouble(value.actuatorWorkJoules)
        && writer.finiteDouble(value.structureEnergyResidualJoules)
        && writer.finiteDouble(value.fluidEnergyResidualJoules)
        && writer.finiteDouble(value.systemEnergyResidualJoules)
        && writer.u8(value.finite ? 1 : 0)
        && writer.u8(value.accepted ? 1 : 0);
}

bool readConservation(Reader& reader,
                      OpenPistonConservationDiagnostics& value) {
    return reader.u32(value.version)
        && value.version == openPistonConservationVersion
        && reader.vector3(value.structureMomentumChangeNewtonSeconds)
        && reader.vector3(value.fluidMomentumChangeNewtonSeconds)
        && reader.vector3(value.pressureImpulseToStructureNewtonSeconds)
        && reader.vector3(value.actuatorImpulseNewtonSeconds)
        && reader.vector3(value.structureMomentumResidualNewtonSeconds)
        && reader.finiteDouble(
            value.structureMomentumResidualNormNewtonSeconds)
        && reader.vector3(value.fluidMomentumResidualNewtonSeconds)
        && reader.finiteDouble(
            value.fluidMomentumResidualNormNewtonSeconds)
        && reader.vector3(value.systemMomentumResidualNewtonSeconds)
        && reader.finiteDouble(
            value.systemMomentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.structureKineticEnergyChangeJoules)
        && reader.finiteDouble(value.fluidKineticEnergyChangeJoules)
        && reader.finiteDouble(value.pressureWorkToStructureJoules)
        && reader.finiteDouble(value.actuatorWorkJoules)
        && reader.finiteDouble(value.structureEnergyResidualJoules)
        && reader.finiteDouble(value.fluidEnergyResidualJoules)
        && reader.finiteDouble(value.systemEnergyResidualJoules)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

OpenPistonCheckpointPersistenceErrorCode structureErrorCode(
    const StructureCheckpointPersistenceErrorCode code) {
    if (code == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
        return OpenPistonCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code == StructureCheckpointPersistenceErrorCode::TopologyMismatch) {
        return OpenPistonCheckpointPersistenceErrorCode::TopologyMismatch;
    }
    return OpenPistonCheckpointPersistenceErrorCode::InvalidData;
}

OpenPistonCheckpointPersistenceErrorCode fluidErrorCode(
    const fluid::MovingInterfaceFluidCheckpointErrorCode code) {
    using Source = fluid::MovingInterfaceFluidCheckpointErrorCode;
    if (code == Source::LimitExceeded) {
        return OpenPistonCheckpointPersistenceErrorCode::LimitExceeded;
    }
    return OpenPistonCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const OpenPistonCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= envelopeBytes
        && limits.maximumCutSurfaceFaces > 0;
}

} // namespace

bool serializeOpenPistonCheckpoint(
    const OpenPistonCaseCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    OpenPistonCheckpointPersistenceError* error,
    const OpenPistonCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                    "invalid open-piston persistence limits");
    }
    try {
        OpenPistonCase validator;
        validator.restore(checkpoint);
        const OpenPistonCaseCheckpoint canonical = validator.checkpoint();
        const auto* detail =
            OpenPistonCheckpointPersistenceAccess::detail(canonical);
        if (detail == nullptr) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                        "open-piston checkpoint payload is missing");
        }
        std::vector<std::uint8_t> structureBytes;
        StructureCheckpointPersistenceError structureError;
        if (!serializeStructureCheckpoint(
                validator.structure(), detail->structure, structureBytes,
                &structureError, limits.structure)) {
            return fail(error, structureErrorCode(structureError.code),
                        "cannot encode open-piston Structure: "
                            + structureError.message);
        }
        std::vector<std::uint8_t> fluidBytes;
        fluid::MovingInterfaceFluidCheckpointError fluidError;
        if (!fluid::serializeMovingInterfaceFluidCheckpoint(
                detail->fluid, fluidBytes, &fluidError, limits.fluid)) {
            return fail(error, fluidErrorCode(fluidError.code),
                        "cannot encode open-piston fluid: "
                            + fluidError.message);
        }

        std::vector<std::uint8_t> payload;
        Writer writer(payload, limits.maximumEncodedBytes - envelopeBytes);
        if (!writer.u32(stateVersion) || !writer.u32(0)
            || !writer.u32(canonical.version) || !writer.u32(0)
            || !writer.u64(canonical.caseDefinitionFingerprint)
            || !writer.u64(canonical.acceptedStepCount)
            || !writer.u64(canonical.topologyRebaseCount)
            || !writer.count(canonical.movingPlaneCoordinate)
            || !writer.finiteDouble(canonical.surfaceOffsetMeters)
            || !writer.count(structureBytes.size())
            || !writer.bytes(structureBytes)
            || !writer.count(fluidBytes.size()) || !writer.bytes(fluidBytes)
            || !writeControlVolume(writer, detail->controlVolumeDiagnostics)
            || !writeRebase(writer, detail->lastRebaseDiagnostics)
            || !writeCutSurface(writer, detail->cutSurfaceDiagnostics, limits)
            || !writeBridge(writer, detail->bridgeDiagnostics)
            || !writeConservation(writer, detail->conservationDiagnostics)
            || !writer.finiteDouble(
                detail->lastRebaseVelocityResidualMetersPerSecond)
            || !writer.u64(detail->topologyRebaseCount)
            || !writer.count(detail->movingPlaneCoordinate)
            || !writer.finiteDouble(detail->surfaceOffsetMeters)) {
            return fail(
                error,
                writer.exceeded()
                    ? OpenPistonCheckpointPersistenceErrorCode::LimitExceeded
                    : OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                "open-piston checkpoint payload is invalid or too large");
        }
        std::vector<std::uint8_t> encoded;
        encoded.reserve(envelopeBytes + payload.size());
        Writer envelope(encoded, limits.maximumEncodedBytes);
        for (const std::uint8_t byte : magic) {
            if (!envelope.u8(byte)) return false;
        }
        if (!envelope.u16(openPistonCheckpointProtocolVersion)
            || !envelope.u16(0) || !envelope.count(payload.size())
            || !envelope.u64(checksum(payload)) || !envelope.bytes(payload)) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                        "open-piston checkpoint exceeds encoded byte limit");
        }
        bytes = std::move(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate open-piston checkpoint encoding");
    } catch (const std::length_error&) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                    "open-piston checkpoint encoding is too large");
    } catch (const std::exception& exception) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                    std::string("open-piston checkpoint validation failed: ")
                        + exception.what());
    }
}

bool deserializeOpenPistonCheckpoint(
    const std::span<const std::uint8_t> bytes,
    OpenPistonCaseCheckpoint& checkpoint,
    OpenPistonCheckpointPersistenceError* error,
    const OpenPistonCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                    "invalid open-piston persistence limits");
    }
    if (bytes.size() > limits.maximumEncodedBytes) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                    "open-piston checkpoint exceeds encoded byte limit");
    }
    if (bytes.size() < envelopeBytes) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::Truncated,
                    "open-piston checkpoint envelope is truncated");
    }
    Reader envelope(bytes);
    for (const std::uint8_t expected : magic) {
        std::uint8_t actual = 0;
        if (!envelope.u8(actual)) return false;
        if (actual != expected) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::InvalidMagic,
                        "open-piston checkpoint magic is invalid");
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
                    "open-piston checkpoint envelope is invalid");
    }
    if (protocol != openPistonCheckpointProtocolVersion) {
        return fail(error,
            OpenPistonCheckpointPersistenceErrorCode::UnsupportedVersion,
            "open-piston checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                    "open-piston checkpoint reserved bits are nonzero");
    }
    const auto payload = envelope.remainingBytes(payloadSize);
    if (envelope.truncated()) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::Truncated,
                    "open-piston checkpoint payload is truncated");
    }
    if (!envelope.atEnd()) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::TrailingData,
                    "open-piston checkpoint has trailing data");
    }
    if (checksum(payload) != expectedChecksum) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::ChecksumMismatch,
                    "open-piston checkpoint checksum does not match");
    }

    try {
        Reader reader(payload);
        OpenPistonCaseCheckpoint candidate;
        std::uint32_t payloadVersion = 0;
        std::uint32_t payloadReserved = 0;
        std::uint32_t checkpointReserved = 0;
        if (!reader.u32(payloadVersion) || !reader.u32(payloadReserved)
            || !reader.u32(candidate.version)
            || !reader.u32(checkpointReserved)
            || !reader.u64(candidate.caseDefinitionFingerprint)
            || !reader.u64(candidate.acceptedStepCount)
            || !reader.u64(candidate.topologyRebaseCount)
            || !reader.count(candidate.movingPlaneCoordinate,
                             limits.maximumCutSurfaceFaces)
            || !reader.finiteDouble(candidate.surfaceOffsetMeters)) {
            return fail(error, readerErrorCode(reader),
                        "open-piston checkpoint payload header is invalid");
        }
        if (payloadVersion != stateVersion
            || candidate.version != openPistonCaseCheckpointVersion) {
            return fail(error,
                OpenPistonCheckpointPersistenceErrorCode::UnsupportedVersion,
                "open-piston checkpoint state version is unsupported");
        }
        if (payloadReserved != 0 || checkpointReserved != 0) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                        "open-piston payload reserved bits are nonzero");
        }
        if (candidate.caseDefinitionFingerprint
            != openPistonCaseDefinitionFingerprint) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::TopologyMismatch,
                        "open-piston checkpoint belongs to a foreign case");
        }

        OpenPistonCase validator;
        const auto templateCheckpoint = validator.checkpoint();
        const auto* templateDetail =
            OpenPistonCheckpointPersistenceAccess::detail(templateCheckpoint);
        if (templateDetail == nullptr) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                        "open-piston topology template is incomplete");
        }
        std::size_t structureSize = 0;
        if (!reader.count(
                structureSize, limits.structure.maximumEncodedBytes)) {
            return fail(error, readerErrorCode(reader),
                        "open-piston Structure length is invalid");
        }
        const auto structureBytes = reader.remainingBytes(structureSize);
        if (reader.truncated()) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::Truncated,
                        "open-piston Structure checkpoint is truncated");
        }
        StructureCheckpoint structure;
        StructureCheckpointPersistenceError structureError;
        if (!deserializeStructureCheckpoint(
                structureBytes, validator.structure(), structure,
                &structureError, limits.structure)) {
            return fail(error, structureErrorCode(structureError.code),
                        "cannot decode open-piston Structure: "
                            + structureError.message);
        }
        std::size_t fluidSize = 0;
        if (!reader.count(fluidSize, limits.fluid.maximumBytes)) {
            return fail(error, readerErrorCode(reader),
                        "open-piston fluid length is invalid");
        }
        const auto fluidBytes = reader.remainingBytes(fluidSize);
        if (reader.truncated()) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::Truncated,
                        "open-piston fluid checkpoint is truncated");
        }
        fluid::MovingInterfaceFluidCheckpoint fluidState;
        fluid::MovingInterfaceFluidCheckpointError fluidError;
        if (!fluid::deserializeMovingInterfaceFluidCheckpoint(
                fluidBytes, fluidState, &fluidError, limits.fluid)) {
            return fail(error, fluidErrorCode(fluidError.code),
                        "cannot decode open-piston fluid: "
                            + fluidError.message);
        }

        auto detail = std::make_shared<OpenPistonCaseCheckpoint::Detail>();
        detail->structure = std::move(structure);
        detail->fluid = std::move(fluidState);
        if (!readControlVolume(reader, detail->controlVolumeDiagnostics,
                               limits.maximumCutSurfaceFaces)
            || !readRebase(reader, detail->lastRebaseDiagnostics,
                           limits.maximumCutSurfaceFaces)
            || !readCutSurface(reader, detail->cutSurfaceDiagnostics, limits)
            || !readBridge(reader, detail->bridgeDiagnostics,
                           limits.maximumCutSurfaceFaces)
            || !readConservation(reader, detail->conservationDiagnostics)
            || !reader.finiteDouble(
                detail->lastRebaseVelocityResidualMetersPerSecond)
            || !reader.u64(detail->topologyRebaseCount)
            || !reader.count(detail->movingPlaneCoordinate,
                             limits.maximumCutSurfaceFaces)
            || !reader.finiteDouble(detail->surfaceOffsetMeters)) {
            return fail(error, readerErrorCode(reader),
                        "open-piston diagnostic payload is invalid");
        }
        if (!reader.atEnd()) {
            return fail(error,
                        OpenPistonCheckpointPersistenceErrorCode::TrailingData,
                        "open-piston payload has trailing data");
        }
        OpenPistonCheckpointPersistenceAccess::setDetail(candidate, detail);
        validator.restore(candidate);
        checkpoint = validator.checkpoint();
        return true;
    } catch (const std::bad_alloc&) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                    "cannot allocate decoded open-piston checkpoint");
    } catch (const std::length_error&) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::LimitExceeded,
                    "decoded open-piston checkpoint is too large");
    } catch (const std::exception& exception) {
        return fail(error,
                    OpenPistonCheckpointPersistenceErrorCode::InvalidData,
                    std::string("open-piston checkpoint validation failed: ")
                        + exception.what());
    }
}

} // namespace simwing::fsi
