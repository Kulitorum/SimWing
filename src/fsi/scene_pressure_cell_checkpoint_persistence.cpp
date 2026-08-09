#include "scene_pressure_cell_checkpoint_persistence.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 9> checkpointMagic{
    'S', 'W', 'P', 'C', 'E', 'L', 'L', '1', '0'};
constexpr std::uint32_t checkpointStateVersion = 10;
constexpr std::size_t checkpointEnvelopeBytes = 29;
constexpr std::size_t solveComponentRecordBytes = 56;
constexpr std::size_t controlVolumeRecordBytes = 72;
constexpr std::size_t linkRecordBytes = 73;
constexpr std::size_t momentumControlVolumeRecordBytes = 152;
constexpr std::size_t wallTractionRecordBytes = 32;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(
    ScenePressureCellCheckpointPersistenceError* error,
    const ScenePressureCellCheckpointPersistenceErrorCode code,
    std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(ScenePressureCellCheckpointPersistenceError* error) {
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

    bool u8(const std::uint8_t value) { return integer(value); }
    bool u16(const std::uint16_t value) { return integer(value); }
    bool u32(const std::uint32_t value) { return integer(value); }
    bool u64(const std::uint64_t value) { return integer(value); }
    bool count(const std::size_t value) {
        return u64(static_cast<std::uint64_t>(value));
    }
    bool boolean(const bool value) { return u8(value ? 1U : 0U); }
    bool finiteDouble(const double value) {
        return std::isfinite(value)
            && u64(std::bit_cast<std::uint64_t>(value));
    }
    bool vector3(const fluid::Vector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool bytes(const std::span<const std::uint8_t> value) {
        if (value.size() > limit_ || bytes_.size() > limit_ - value.size()) {
            return false;
        }
        bytes_.insert(bytes_.end(), value.begin(), value.end());
        return true;
    }

private:
    template<typename Unsigned>
    bool integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        std::array<std::uint8_t, sizeof(Unsigned)> encoded{};
        for (std::size_t index = 0; index < encoded.size(); ++index) {
            encoded[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return bytes(encoded);
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t limit_ = 0;
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
    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded) || encoded > 1) {
            return false;
        }
        value = encoded != 0;
        return true;
    }
    bool finiteDouble(double& value) {
        std::uint64_t encoded = 0;
        if (!u64(encoded)) {
            return false;
        }
        value = std::bit_cast<double>(encoded);
        return std::isfinite(value);
    }
    bool vector3(fluid::Vector3& value) {
        return finiteDouble(value.x) && finiteDouble(value.y)
            && finiteDouble(value.z);
    }
    bool bytes(const std::size_t size, std::span<const std::uint8_t>& value) {
        if (size > remaining()) {
            truncated_ = true;
            return false;
        }
        value = bytes_.subspan(position_, size);
        position_ += size;
        return true;
    }
    bool fixedRecords(const std::size_t count,
                      const std::size_t bytesPerRecord) {
        if (count > remaining() / bytesPerRecord) {
            truncated_ = true;
            return false;
        }
        return true;
    }
    [[nodiscard]] std::size_t remaining() const noexcept {
        return bytes_.size() - position_;
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
        if (sizeof(Unsigned) > remaining()) {
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

ScenePressureCellCheckpointPersistenceErrorCode readerErrorCode(
    const Reader& reader) noexcept {
    if (reader.limitExceeded()) {
        return ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return ScenePressureCellCheckpointPersistenceErrorCode::Truncated;
    }
    return ScenePressureCellCheckpointPersistenceErrorCode::InvalidData;
}

ScenePressureCellCheckpointPersistenceErrorCode structureErrorCode(
    const StructureCheckpointPersistenceErrorCode code) noexcept {
    if (code == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
        return ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code == StructureCheckpointPersistenceErrorCode::TopologyMismatch) {
        return ScenePressureCellCheckpointPersistenceErrorCode::TopologyMismatch;
    }
    return ScenePressureCellCheckpointPersistenceErrorCode::InvalidData;
}

ScenePressureCellCheckpointPersistenceErrorCode mimeticStateErrorCode(
    const SceneFluidMimeticPressureStatePersistenceErrorCode code) noexcept {
    if (code
        == SceneFluidMimeticPressureStatePersistenceErrorCode::
            LimitExceeded) {
        return ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded;
    }
    if (code
        == SceneFluidMimeticPressureStatePersistenceErrorCode::
            TopologyMismatch) {
        return ScenePressureCellCheckpointPersistenceErrorCode::
            TopologyMismatch;
    }
    return ScenePressureCellCheckpointPersistenceErrorCode::InvalidData;
}

bool validLimits(const ScenePressureCellCheckpointPersistenceLimits& limits) {
    return limits.maximumEncodedBytes >= checkpointEnvelopeBytes
        && limits.maximumControlVolumes > 0
        && limits.maximumControlVolumes
            <= std::numeric_limits<std::size_t>::max() / 3
        && limits.maximumLinks > 0
        && limits.maximumSolveComponents > 0
        && limits.maximumProjectionStorageBytes > 0
        && limits.maximumMomentumStorageBytes > 0
        && limits.maximumWallTractions > 0
        && limits.maximumWallTractionStorageBytes > 0
        && limits.mimeticPressureState.maximumEncodedBytes > 0
        && limits.mimeticPressureState.maximumControlCells > 0
        && limits.mimeticPressureState.maximumReducedTraces > 0
        && limits.structure.maximumEncodedBytes > 0;
}

bool writeProjectionSettings(
    Writer& writer,
    const SceneFluidPressureProjectionSettings& settings) {
    return writer.finiteDouble(settings.densityKgPerCubicMeter)
        && writer.finiteDouble(settings.timeStepSeconds)
        && writer.finiteDouble(settings
            .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond)
        && writer.finiteDouble(
            settings.relativeCorrectedVolumeRateTolerance)
        && writer.finiteDouble(settings.pressureSolve
            .absoluteResidualTolerancePascalsMeters)
        && writer.finiteDouble(
            settings.pressureSolve.relativeResidualTolerance)
        && writer.finiteDouble(settings.pressureSolve
            .absoluteComponentCompatibilityTolerancePascalsMeters)
        && writer.count(settings.pressureSolve.maximumIterations);
}

bool readProjectionSettings(
    Reader& reader,
    SceneFluidPressureProjectionSettings& settings) {
    return reader.finiteDouble(settings.densityKgPerCubicMeter)
        && reader.finiteDouble(settings.timeStepSeconds)
        && reader.finiteDouble(settings
            .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond)
        && reader.finiteDouble(
            settings.relativeCorrectedVolumeRateTolerance)
        && reader.finiteDouble(settings.pressureSolve
            .absoluteResidualTolerancePascalsMeters)
        && reader.finiteDouble(
            settings.pressureSolve.relativeResidualTolerance)
        && reader.finiteDouble(settings.pressureSolve
            .absoluteComponentCompatibilityTolerancePascalsMeters)
        && reader.count(
            settings.pressureSolve.maximumIterations,
            std::numeric_limits<std::size_t>::max());
}

bool writeSolveDiagnostics(
    Writer& writer,
    const SceneFluidPressureSolveDiagnostics& diagnostics) {
    if (!writer.boolean(diagnostics.compatible)
        || !writer.boolean(diagnostics.converged)
        || !writer.boolean(diagnostics.finite)
        || !writer.u64(diagnostics.pressureOperatorFingerprint)
        || !writer.count(diagnostics.rowCount)
        || !writer.count(diagnostics.componentCount)
        || !writer.count(diagnostics.iterationCount)
        || !writer.finiteDouble(
            diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters)
        || !writer.finiteDouble(
            diagnostics.initialResidualL2PascalsMeters)
        || !writer.finiteDouble(diagnostics.finalResidualL2PascalsMeters)
        || !writer.finiteDouble(
            diagnostics.finalResidualMaximumPascalsMeters)
        || !writer.count(diagnostics.components.size())) {
        return false;
    }
    for (const auto& component : diagnostics.components) {
        if (!writer.count(component.componentIndex)
            || !writer.count(component.controlVolumeCount)
            || !writer.count(component.gaugeControlVolumeIndex)
            || !writer.finiteDouble(
                component.rightHandSideSumPascalsMeters)
            || !writer.finiteDouble(
                component.compatibilityCorrectionPascalsMeters)
            || !writer.finiteDouble(component.pressureGaugeBeforePascals)
            || !writer.finiteDouble(component.pressureGaugeAfterPascals)) {
            return false;
        }
    }
    return true;
}

bool readSolveDiagnostics(
    Reader& reader,
    SceneFluidPressureSolveDiagnostics& diagnostics,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    std::size_t componentRecords = 0;
    if (!reader.boolean(diagnostics.compatible)
        || !reader.boolean(diagnostics.converged)
        || !reader.boolean(diagnostics.finite)
        || !reader.u64(diagnostics.pressureOperatorFingerprint)
        || !reader.count(
            diagnostics.rowCount, limits.maximumControlVolumes)
        || !reader.count(
            diagnostics.componentCount, limits.maximumSolveComponents)
        || !reader.count(
            diagnostics.iterationCount,
            std::numeric_limits<std::size_t>::max())
        || !reader.finiteDouble(
            diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters)
        || !reader.finiteDouble(
            diagnostics.initialResidualL2PascalsMeters)
        || !reader.finiteDouble(diagnostics.finalResidualL2PascalsMeters)
        || !reader.finiteDouble(
            diagnostics.finalResidualMaximumPascalsMeters)
        || !reader.count(
            componentRecords, limits.maximumSolveComponents)
        || !reader.fixedRecords(
            componentRecords, solveComponentRecordBytes)) {
        return false;
    }
    diagnostics.components.resize(componentRecords);
    for (auto& component : diagnostics.components) {
        if (!reader.count(
                component.componentIndex,
                limits.maximumSolveComponents)
            || !reader.count(
                component.controlVolumeCount,
                limits.maximumControlVolumes)
            || !reader.count(
                component.gaugeControlVolumeIndex,
                limits.maximumControlVolumes)
            || !reader.finiteDouble(
                component.rightHandSideSumPascalsMeters)
            || !reader.finiteDouble(
                component.compatibilityCorrectionPascalsMeters)
            || !reader.finiteDouble(
                component.pressureGaugeBeforePascals)
            || !reader.finiteDouble(
                component.pressureGaugeAfterPascals)) {
            return false;
        }
    }
    return true;
}

bool writeProjectionDiagnostics(
    Writer& writer,
    const SceneFluidPressureProjectionDiagnostics& diagnostics) {
    return writer.boolean(diagnostics.accepted)
        && writer.boolean(diagnostics.finite)
        && writer.count(diagnostics.controlVolumeCount)
        && writer.count(diagnostics.linkCount)
        && writer.count(diagnostics.authoredOpeningLinkCount)
        && writer.boolean(diagnostics.usesMovingVolumeRates)
        && writer.finiteDouble(diagnostics
            .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .predictedNetOutwardVolumeRateL2CubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .maximumPredictedComponentBalanceResidualCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .predictedContinuityResidualL2CubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .predictedContinuityResidualMaximumCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .maximumPredictedComponentContinuityResidualCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .correctedNetOutwardVolumeRateL2CubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .correctedContinuityResidualL2CubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .correctedContinuityResidualMaximumCubicMetersPerSecond)
        && writer.finiteDouble(diagnostics
            .maximumCorrectedComponentContinuityResidualCubicMetersPerSecond)
        && writeSolveDiagnostics(writer, diagnostics.pressureSolve);
}

bool readProjectionDiagnostics(
    Reader& reader,
    SceneFluidPressureProjectionDiagnostics& diagnostics,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    return reader.boolean(diagnostics.accepted)
        && reader.boolean(diagnostics.finite)
        && reader.count(
            diagnostics.controlVolumeCount,
            limits.maximumControlVolumes)
        && reader.count(diagnostics.linkCount, limits.maximumLinks)
        && reader.count(
            diagnostics.authoredOpeningLinkCount, limits.maximumLinks)
        && reader.boolean(diagnostics.usesMovingVolumeRates)
        && reader.finiteDouble(diagnostics
            .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .predictedNetOutwardVolumeRateL2CubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .maximumPredictedComponentBalanceResidualCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .predictedContinuityResidualL2CubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .predictedContinuityResidualMaximumCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .maximumPredictedComponentContinuityResidualCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .correctedNetOutwardVolumeRateL2CubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .correctedContinuityResidualL2CubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .correctedContinuityResidualMaximumCubicMetersPerSecond)
        && reader.finiteDouble(diagnostics
            .maximumCorrectedComponentContinuityResidualCubicMetersPerSecond)
        && readSolveDiagnostics(reader, diagnostics.pressureSolve, limits);
}

bool writeProjection(
    Writer& writer,
    const SceneFluidPressureProjection& projection) {
    if (!writer.u32(projection.version)
        || !writer.u64(projection.fingerprint)
        || !writer.u64(projection.pressureOperatorFingerprint)
        || !writer.u64(projection.pressureFaceLinkFingerprint)
        || !writer.u64(projection.pressureControlVolumeFingerprint)
        || !writer.u64(projection.pressureVolumeRateFingerprint)
        || !writer.u64(projection.openingFluxFingerprint)
        || !writer.u64(projection.velocityFingerprint)
        || !writer.u64(projection.linkFlowContinuationFingerprint)
        || !writer.u64(projection.regionLinkFlowPredictionFingerprint)
        || !writer.u64(projection.regionWallExchangeFingerprint)
        || !writer.u64(projection.acceptedStepCount)
        || !writer.finiteDouble(projection.simulationTimeSeconds)
        || !writer.count(projection.cellCounts.x)
        || !writer.count(projection.cellCounts.y)
        || !writer.count(projection.cellCounts.z)
        || !writer.vector3(projection.lowerMeters)
        || !writer.vector3(projection.upperMeters)
        || !writeProjectionSettings(writer, projection.settings)
        || !writer.count(projection.ownedStorageBytes)
        || !writeProjectionDiagnostics(writer, projection.diagnostics)
        || !writer.count(projection.pressurePascals.size())) {
        return false;
    }
    for (const double pressure : projection.pressurePascals) {
        if (!writer.finiteDouble(pressure)) {
            return false;
        }
    }
    if (!writer.count(projection.controlVolumes.size())) {
        return false;
    }
    for (const auto& control : projection.controlVolumes) {
        if (!writer.count(control.controlVolumeIndex)
            || !writer.u64(control.stableId)
            || !writer.count(control.componentIndex)
            || !writer.finiteDouble(control
                .geometryVolumeChangeRateCubicMetersPerSecond)
            || !writer.finiteDouble(control
                .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || !writer.finiteDouble(control
                .predictedContinuityResidualCubicMetersPerSecond)
            || !writer.finiteDouble(
                control.integratedRightHandSidePascalsMeters)
            || !writer.finiteDouble(control
                .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || !writer.finiteDouble(control
                .correctedContinuityResidualCubicMetersPerSecond)) {
            return false;
        }
    }
    if (!writer.count(projection.links.size())) {
        return false;
    }
    for (const auto& link : projection.links) {
        if (!writer.count(link.linkIndex)
            || !writer.u64(link.stableId)
            || !writer.count(link.faceIndex)
            || !writer.u8(static_cast<std::uint8_t>(link.kind))
            || !writer.count(link.minusControlVolumeIndex)
            || !writer.count(link.plusControlVolumeIndex)
            || !writer.u64(link.openingPatchStableId)
            || !writer.finiteDouble(link
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond)
            || !writer.finiteDouble(link
                .pressureCorrectionVolumeFlowRateCubicMetersPerSecond)
            || !writer.finiteDouble(link
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond)) {
            return false;
        }
    }
    return true;
}

bool readProjection(
    Reader& reader,
    SceneFluidPressureProjection& projection,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    if (!reader.u32(projection.version)
        || !reader.u64(projection.fingerprint)
        || !reader.u64(projection.pressureOperatorFingerprint)
        || !reader.u64(projection.pressureFaceLinkFingerprint)
        || !reader.u64(projection.pressureControlVolumeFingerprint)
        || !reader.u64(projection.pressureVolumeRateFingerprint)
        || !reader.u64(projection.openingFluxFingerprint)
        || !reader.u64(projection.velocityFingerprint)
        || !reader.u64(projection.linkFlowContinuationFingerprint)
        || !reader.u64(projection.regionLinkFlowPredictionFingerprint)
        || !reader.u64(projection.regionWallExchangeFingerprint)
        || !reader.u64(projection.acceptedStepCount)
        || !reader.finiteDouble(projection.simulationTimeSeconds)
        || !reader.count(
            projection.cellCounts.x,
            std::numeric_limits<std::size_t>::max())
        || !reader.count(
            projection.cellCounts.y,
            std::numeric_limits<std::size_t>::max())
        || !reader.count(
            projection.cellCounts.z,
            std::numeric_limits<std::size_t>::max())
        || !reader.vector3(projection.lowerMeters)
        || !reader.vector3(projection.upperMeters)
        || !readProjectionSettings(reader, projection.settings)
        || !reader.count(
            projection.ownedStorageBytes,
            limits.maximumProjectionStorageBytes)
        || !readProjectionDiagnostics(
            reader, projection.diagnostics, limits)) {
        return false;
    }

    std::size_t count = 0;
    if (!reader.count(count, limits.maximumControlVolumes)
        || !reader.fixedRecords(count, sizeof(double))) {
        return false;
    }
    projection.pressurePascals.resize(count);
    for (double& pressure : projection.pressurePascals) {
        if (!reader.finiteDouble(pressure)) {
            return false;
        }
    }

    if (!reader.count(count, limits.maximumControlVolumes)
        || !reader.fixedRecords(count, controlVolumeRecordBytes)) {
        return false;
    }
    projection.controlVolumes.resize(count);
    for (auto& control : projection.controlVolumes) {
        if (!reader.count(
                control.controlVolumeIndex, limits.maximumControlVolumes)
            || !reader.u64(control.stableId)
            || !reader.count(
                control.componentIndex, limits.maximumSolveComponents)
            || !reader.finiteDouble(control
                .geometryVolumeChangeRateCubicMetersPerSecond)
            || !reader.finiteDouble(control
                .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || !reader.finiteDouble(control
                .predictedContinuityResidualCubicMetersPerSecond)
            || !reader.finiteDouble(
                control.integratedRightHandSidePascalsMeters)
            || !reader.finiteDouble(control
                .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || !reader.finiteDouble(control
                .correctedContinuityResidualCubicMetersPerSecond)) {
            return false;
        }
    }

    if (!reader.count(count, limits.maximumLinks)
        || !reader.fixedRecords(count, linkRecordBytes)) {
        return false;
    }
    projection.links.resize(count);
    for (auto& link : projection.links) {
        std::uint8_t kind = 0;
        if (!reader.count(link.linkIndex, limits.maximumLinks)
            || !reader.u64(link.stableId)
            || !reader.count(
                link.faceIndex, std::numeric_limits<std::size_t>::max())
            || !reader.u8(kind)
            || (kind != static_cast<std::uint8_t>(
                    SceneFluidPressureFaceLinkKind::SameRegion)
                && kind != static_cast<std::uint8_t>(
                    SceneFluidPressureFaceLinkKind::AuthoredOpening))
            || !reader.count(
                link.minusControlVolumeIndex,
                limits.maximumControlVolumes)
            || !reader.count(
                link.plusControlVolumeIndex,
                limits.maximumControlVolumes)
            || !reader.u64(link.openingPatchStableId)
            || !reader.finiteDouble(link
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond)
            || !reader.finiteDouble(link
                .pressureCorrectionVolumeFlowRateCubicMetersPerSecond)
            || !reader.finiteDouble(link
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond)) {
            return false;
        }
        link.kind = static_cast<SceneFluidPressureFaceLinkKind>(kind);
    }
    validateSceneFluidPressureProjectionIntegrity(projection);
    return true;
}

bool writeWallTractions(
    Writer& writer,
    const SceneFluidAcceptedWallTractionSet& wall) {
    validateSceneFluidAcceptedWallTractionSetIntegrity(wall);
    if (!writer.u32(wall.version)
        || !writer.u64(wall.fingerprint)
        || !writer.u64(wall.wallExchangeFingerprint)
        || !writer.u64(wall.quadratureFingerprint)
        || !writer.u64(wall.surfaceDefinitionFingerprint)
        || !writer.u64(wall.surfaceStateFingerprint)
        || !writer.u64(wall.structureDefinitionFingerprint)
        || !writer.u64(wall.acceptedStepCount)
        || !writer.finiteDouble(wall.simulationTimeSeconds)
        || !writer.count(wall.ownedStorageBytes)
        || !writer.count(wall.tractions.size())) {
        return false;
    }
    for (const auto& traction : wall.tractions) {
        if (!writer.u64(traction.stableId)
            || !writer.finiteDouble(traction.tractionPascals.x)
            || !writer.finiteDouble(traction.tractionPascals.y)
            || !writer.finiteDouble(traction.tractionPascals.z)) {
            return false;
        }
    }
    return true;
}

bool readWallTractions(
    Reader& reader,
    SceneFluidAcceptedWallTractionSet& wall,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    std::size_t count = 0;
    if (!reader.u32(wall.version)
        || !reader.u64(wall.fingerprint)
        || !reader.u64(wall.wallExchangeFingerprint)
        || !reader.u64(wall.quadratureFingerprint)
        || !reader.u64(wall.surfaceDefinitionFingerprint)
        || !reader.u64(wall.surfaceStateFingerprint)
        || !reader.u64(wall.structureDefinitionFingerprint)
        || !reader.u64(wall.acceptedStepCount)
        || !reader.finiteDouble(wall.simulationTimeSeconds)
        || !reader.count(
            wall.ownedStorageBytes,
            limits.maximumWallTractionStorageBytes)
        || !reader.count(count, limits.maximumWallTractions)
        || !reader.fixedRecords(count, wallTractionRecordBytes)) {
        return false;
    }
    wall.tractions.resize(count);
    for (auto& traction : wall.tractions) {
        if (!reader.u64(traction.stableId)
            || !reader.finiteDouble(traction.tractionPascals.x)
            || !reader.finiteDouble(traction.tractionPascals.y)
            || !reader.finiteDouble(traction.tractionPascals.z)) {
            return false;
        }
    }
    validateSceneFluidAcceptedWallTractionSetIntegrity(wall);
    return true;
}

bool writeRegionMomentum(
    Writer& writer,
    const SceneFluidRegionMomentumState& momentum) {
    validateSceneFluidRegionMomentumStateIntegrity(momentum);
    const auto& diagnostics = momentum.diagnostics;
    if (!writer.u32(momentum.version)
        || !writer.u64(momentum.fingerprint)
        || !writer.u64(momentum.pressureProjectionFingerprint)
        || !writer.u64(momentum.pressureControlVolumeFingerprint)
        || !writer.u64(momentum.pressureFaceLinkFingerprint)
        || !writer.u64(momentum.openingPatchFingerprint)
        || !writer.u64(momentum.fallbackVelocityFingerprint)
        || !writer.u64(momentum.acceptedStepCount)
        || !writer.finiteDouble(momentum.simulationTimeSeconds)
        || !writer.finiteDouble(momentum.densityKgPerCubicMeter)
        || !writer.count(momentum.cellCounts.x)
        || !writer.count(momentum.cellCounts.y)
        || !writer.count(momentum.cellCounts.z)
        || !writer.vector3(momentum.lowerMeters)
        || !writer.vector3(momentum.upperMeters)
        || !writer.count(momentum.ownedStorageBytes)
        || !writer.count(diagnostics.controlVolumeCount)
        || !writer.count(diagnostics.linkCount)
        || !writer.count(diagnostics.openingLinkCount)
        || !writer.count(diagnostics.embeddedOpeningLinkCount)
        || !writer.count(diagnostics.normalEquationControlCount)
        || !writer.count(diagnostics.sampledComponentCount)
        || !writer.count(diagnostics.fallbackComponentCount)
        || !writer.vector3(
            diagnostics.totalMomentumKilogramMetersPerSecond)
        || !writer.finiteDouble(diagnostics.kineticEnergyJoules)
        || !writer.finiteDouble(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond)
        || !writer.finiteDouble(
            diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond)
        || !writer.boolean(diagnostics.finite)
        || !writer.count(momentum.controlVolumes.size())) {
        return false;
    }
    for (const auto& control : momentum.controlVolumes) {
        if (!writer.count(control.controlVolumeIndex)
            || !writer.u64(control.stableId)
            || !writer.count(control.cellIndex)
            || !writer.count(control.regionIndex)
            || !writer.u64(control.regionId)
            || !writer.count(control.componentIndex)
            || !writer.finiteDouble(control.volumeCubicMeters)
            || !writer.vector3(control.velocityMetersPerSecond)
            || !writer.vector3(control.momentumKilogramMetersPerSecond)) {
            return false;
        }
        for (const double area : control.sampledFaceAreaSquareMeters) {
            if (!writer.finiteDouble(area)) {
                return false;
            }
        }
        for (const std::size_t count : control.sampledLinkCounts) {
            if (!writer.count(count)) {
                return false;
            }
        }
    }
    return true;
}

bool readRegionMomentum(
    Reader& reader,
    SceneFluidRegionMomentumState& momentum,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    auto& diagnostics = momentum.diagnostics;
    if (!reader.u32(momentum.version)
        || !reader.u64(momentum.fingerprint)
        || !reader.u64(momentum.pressureProjectionFingerprint)
        || !reader.u64(momentum.pressureControlVolumeFingerprint)
        || !reader.u64(momentum.pressureFaceLinkFingerprint)
        || !reader.u64(momentum.openingPatchFingerprint)
        || !reader.u64(momentum.fallbackVelocityFingerprint)
        || !reader.u64(momentum.acceptedStepCount)
        || !reader.finiteDouble(momentum.simulationTimeSeconds)
        || !reader.finiteDouble(momentum.densityKgPerCubicMeter)
        || !reader.count(
            momentum.cellCounts.x, std::numeric_limits<std::size_t>::max())
        || !reader.count(
            momentum.cellCounts.y, std::numeric_limits<std::size_t>::max())
        || !reader.count(
            momentum.cellCounts.z, std::numeric_limits<std::size_t>::max())
        || !reader.vector3(momentum.lowerMeters)
        || !reader.vector3(momentum.upperMeters)
        || !reader.count(
            momentum.ownedStorageBytes,
            limits.maximumMomentumStorageBytes)
        || !reader.count(
            diagnostics.controlVolumeCount, limits.maximumControlVolumes)
        || !reader.count(diagnostics.linkCount, limits.maximumLinks)
        || !reader.count(diagnostics.openingLinkCount, limits.maximumLinks)
        || !reader.count(
            diagnostics.embeddedOpeningLinkCount, limits.maximumLinks)
        || !reader.count(
            diagnostics.normalEquationControlCount,
            limits.maximumControlVolumes)
        || !reader.count(
            diagnostics.sampledComponentCount,
            3 * limits.maximumControlVolumes)
        || !reader.count(
            diagnostics.fallbackComponentCount,
            3 * limits.maximumControlVolumes)
        || !reader.vector3(
            diagnostics.totalMomentumKilogramMetersPerSecond)
        || !reader.finiteDouble(diagnostics.kineticEnergyJoules)
        || !reader.finiteDouble(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond)
        || !reader.finiteDouble(
            diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond)
        || !reader.boolean(diagnostics.finite)) {
        return false;
    }
    std::size_t count = 0;
    if (!reader.count(count, limits.maximumControlVolumes)
        || !reader.fixedRecords(
            count, momentumControlVolumeRecordBytes)) {
        return false;
    }
    momentum.controlVolumes.resize(count);
    for (auto& control : momentum.controlVolumes) {
        if (!reader.count(
                control.controlVolumeIndex, limits.maximumControlVolumes)
            || !reader.u64(control.stableId)
            || !reader.count(
                control.cellIndex, std::numeric_limits<std::size_t>::max())
            || !reader.count(
                control.regionIndex, limits.maximumControlVolumes)
            || !reader.u64(control.regionId)
            || !reader.count(
                control.componentIndex, limits.maximumControlVolumes)
            || !reader.finiteDouble(control.volumeCubicMeters)
            || !reader.vector3(control.velocityMetersPerSecond)
            || !reader.vector3(control.momentumKilogramMetersPerSecond)) {
            return false;
        }
        for (double& area : control.sampledFaceAreaSquareMeters) {
            if (!reader.finiteDouble(area)) {
                return false;
            }
        }
        for (std::size_t& linkCount : control.sampledLinkCounts) {
            if (!reader.count(linkCount, limits.maximumLinks)) {
                return false;
            }
        }
    }
    validateSceneFluidRegionMomentumStateIntegrity(momentum);
    return true;
}

} // namespace

bool serializeScenePressureCellCheckpoint(
    const ScenePressureCellCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    ScenePressureCellCheckpointPersistenceError* error,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits)) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell checkpoint limits are invalid");
    }
    try {
        ScenePressureCellCase validator(
            checkpoint.coupling
                    .mimeticPressureAuditSettingsFingerprint
                != 0);
        validator.restore(checkpoint);

        std::vector<std::uint8_t> structureBytes;
        StructureCheckpointPersistenceError structureError;
        if (!serializeStructureCheckpoint(
                validator.structure(), checkpoint.coupling.structure,
                structureBytes, &structureError, limits.structure)) {
            return fail(
                error, structureErrorCode(structureError.code),
                "scene pressure cell Structure payload: "
                    + structureError.message);
        }
        if (checkpoint.regionMomentum
            && (checkpoint.regionMomentum->controlVolumes.size()
                    > limits.maximumControlVolumes
                || checkpoint.regionMomentum->ownedStorageBytes
                    > limits.maximumMomentumStorageBytes)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell region momentum exceeds its limit");
        }
        if (checkpoint.coupling.wallTractions
            && (checkpoint.coupling.wallTractions->tractions.size()
                    > limits.maximumWallTractions
                || checkpoint.coupling.wallTractions->ownedStorageBytes
                    > limits.maximumWallTractionStorageBytes)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell wall traction exceeds its limit");
        }
        std::vector<std::uint8_t> mimeticPressureStateBytes;
        if (checkpoint.coupling.mimeticPressureState) {
            const auto topology =
                validator.rebuildMimeticPressureAuditTopology(
                    checkpoint.coupling.structure);
            SceneFluidMimeticPressureStatePersistenceError mimeticError;
            if (!serializeSceneFluidMimeticPressureState(
                    *checkpoint.coupling.mimeticPressureState,
                    topology.controlCells, topology.fullTraceSystem,
                    topology.condensedTraceSystem,
                    mimeticPressureStateBytes, &mimeticError,
                    limits.mimeticPressureState)) {
                return fail(
                    error, mimeticStateErrorCode(mimeticError.code),
                    "scene pressure cell mimetic pressure state: "
                        + mimeticError.message);
            }
        }

        std::vector<std::uint8_t> payload;
        Writer payloadWriter(
            payload, limits.maximumEncodedBytes - checkpointEnvelopeBytes);
        if (!payloadWriter.u32(checkpointStateVersion)
            || !payloadWriter.u32(checkpoint.version)
            || !payloadWriter.u64(checkpoint.coupling
                .mimeticPressureAuditSettingsFingerprint)
            || !payloadWriter.count(structureBytes.size())
            || !payloadWriter.bytes(structureBytes)
            || !payloadWriter.boolean(
                checkpoint.coupling.pressureProjection.has_value())
            || (checkpoint.coupling.pressureProjection
                && !writeProjection(
                    payloadWriter,
                    *checkpoint.coupling.pressureProjection))
            || !payloadWriter.boolean(
                checkpoint.coupling.wallTractions.has_value())
            || (checkpoint.coupling.wallTractions
                && !writeWallTractions(
                    payloadWriter, *checkpoint.coupling.wallTractions))
            || !payloadWriter.boolean(
                checkpoint.regionMomentum.has_value())
            || (checkpoint.regionMomentum
                && !writeRegionMomentum(
                    payloadWriter, *checkpoint.regionMomentum))
            || !payloadWriter.boolean(
                checkpoint.coupling.mimeticPressureState.has_value())
            || (checkpoint.coupling.mimeticPressureState
                && (!payloadWriter.count(
                        mimeticPressureStateBytes.size())
                    || !payloadWriter.bytes(
                        mimeticPressureStateBytes)))) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell checkpoint exceeds its byte limit");
        }

        std::vector<std::uint8_t> encoded;
        Writer writer(encoded, limits.maximumEncodedBytes);
        if (!writer.bytes(checkpointMagic)
            || !writer.u16(scenePressureCellCheckpointProtocolVersion)
            || !writer.u16(0)
            || !writer.u64(static_cast<std::uint64_t>(payload.size()))
            || !writer.u64(checksum(payload))
            || !writer.bytes(payload)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell checkpoint envelope exceeds its byte limit");
        }
        bytes.swap(encoded);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell checkpoint allocation exceeded available memory");
    } catch (const std::exception& exception) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::InvalidData,
            exception.what());
    }
}

bool deserializeScenePressureCellCheckpoint(
    const std::span<const std::uint8_t> bytes,
    ScenePressureCellCheckpoint& checkpoint,
    ScenePressureCellCheckpointPersistenceError* error,
    const ScenePressureCellCheckpointPersistenceLimits& limits) {
    clearError(error);
    if (!validLimits(limits) || bytes.size() > limits.maximumEncodedBytes) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell checkpoint exceeds its byte limit");
    }
    if (bytes.size() < checkpointEnvelopeBytes) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
            "scene pressure cell checkpoint envelope is truncated");
    }
    try {
        Reader reader(bytes);
        std::span<const std::uint8_t> magic;
        std::uint16_t protocolVersion = 0;
        std::uint16_t reserved = 0;
        std::uint64_t payloadSize = 0;
        std::uint64_t storedChecksum = 0;
        if (!reader.bytes(checkpointMagic.size(), magic)
            || !reader.u16(protocolVersion)
            || !reader.u16(reserved)
            || !reader.u64(payloadSize)
            || !reader.u64(storedChecksum)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
                "scene pressure cell checkpoint envelope is truncated");
        }
        if (!std::equal(magic.begin(), magic.end(), checkpointMagic.begin())) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::InvalidMagic,
                "scene pressure cell checkpoint magic is invalid");
        }
        if (protocolVersion != scenePressureCellCheckpointProtocolVersion
            || reserved != 0) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::UnsupportedVersion,
                "scene pressure cell checkpoint protocol is unsupported");
        }
        if (payloadSize > std::numeric_limits<std::size_t>::max()) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell checkpoint payload is not representable");
        }
        const std::size_t expectedSize =
            static_cast<std::size_t>(payloadSize);
        if (expectedSize > reader.remaining()) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
                "scene pressure cell checkpoint payload is truncated");
        }
        if (expectedSize < reader.remaining()) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::TrailingData,
                "scene pressure cell checkpoint has trailing data");
        }
        std::span<const std::uint8_t> payload;
        if (!reader.bytes(expectedSize, payload)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
                "scene pressure cell checkpoint payload is truncated");
        }
        if (checksum(payload) != storedChecksum) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::ChecksumMismatch,
                "scene pressure cell checkpoint checksum does not match");
        }

        Reader payloadReader(payload);
        std::uint32_t stateVersion = 0;
        std::uint32_t checkpointVersion = 0;
        std::uint64_t mimeticSettingsFingerprint = 0;
        std::uint64_t structureSize = 0;
        if (!payloadReader.u32(stateVersion)
            || !payloadReader.u32(checkpointVersion)
            || !payloadReader.u64(mimeticSettingsFingerprint)
            || !payloadReader.u64(structureSize)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
                "scene pressure cell checkpoint metadata is truncated");
        }
        if (stateVersion != checkpointStateVersion) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::UnsupportedVersion,
                "scene pressure cell checkpoint state is unsupported");
        }
        auto validator = std::make_unique<ScenePressureCellCase>(
            mimeticSettingsFingerprint != 0);
        ScenePressureCellCheckpoint candidate = validator->checkpoint();
        candidate.version = checkpointVersion;
        candidate.coupling.mimeticPressureAuditSettingsFingerprint =
            mimeticSettingsFingerprint;
        if (structureSize > limits.structure.maximumEncodedBytes
            || structureSize > std::numeric_limits<std::size_t>::max()) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
                "scene pressure cell Structure payload exceeds its limit");
        }
        std::span<const std::uint8_t> structureBytes;
        if (!payloadReader.bytes(
                static_cast<std::size_t>(structureSize), structureBytes)) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
                "scene pressure cell Structure payload is truncated");
        }
        StructureCheckpointPersistenceError structureError;
        if (!deserializeStructureCheckpoint(
                structureBytes, validator->structure(),
                candidate.coupling.structure, &structureError,
                limits.structure)) {
            return fail(
                error, structureErrorCode(structureError.code),
                "scene pressure cell Structure payload: "
                    + structureError.message);
        }

        bool hasProjection = false;
        if (!payloadReader.boolean(hasProjection)) {
            return fail(
                error, readerErrorCode(payloadReader),
                "scene pressure cell pressure marker is invalid");
        }
        if (hasProjection) {
            SceneFluidPressureProjection projection;
            if (!readProjection(payloadReader, projection, limits)) {
                return fail(
                    error, readerErrorCode(payloadReader),
                    "scene pressure cell pressure payload is invalid");
            }
            candidate.coupling.pressureProjection = std::move(projection);
        } else {
            candidate.coupling.pressureProjection.reset();
        }
        bool hasWallTractions = false;
        if (!payloadReader.boolean(hasWallTractions)) {
            return fail(
                error, readerErrorCode(payloadReader),
                "scene pressure cell wall-traction marker is invalid");
        }
        if (hasWallTractions) {
            SceneFluidAcceptedWallTractionSet wall;
            if (!readWallTractions(payloadReader, wall, limits)) {
                return fail(
                    error, readerErrorCode(payloadReader),
                    "scene pressure cell wall-traction payload is invalid");
            }
            candidate.coupling.wallTractions = std::move(wall);
        } else {
            candidate.coupling.wallTractions.reset();
        }
        bool hasRegionMomentum = false;
        if (!payloadReader.boolean(hasRegionMomentum)) {
            return fail(
                error, readerErrorCode(payloadReader),
                "scene pressure cell region-momentum marker is invalid");
        }
        if (hasRegionMomentum) {
            SceneFluidRegionMomentumState momentum;
            if (!readRegionMomentum(payloadReader, momentum, limits)) {
                return fail(
                    error, readerErrorCode(payloadReader),
                    "scene pressure cell region-momentum payload is invalid");
            }
            candidate.regionMomentum = std::move(momentum);
        } else {
            candidate.regionMomentum.reset();
        }
        bool hasMimeticPressureState = false;
        if (!payloadReader.boolean(hasMimeticPressureState)) {
            return fail(
                error, readerErrorCode(payloadReader),
                "scene pressure cell mimetic pressure marker is invalid");
        }
        if (hasMimeticPressureState) {
            std::size_t mimeticStateSize = 0;
            if (!payloadReader.count(
                    mimeticStateSize,
                    limits.mimeticPressureState.maximumEncodedBytes)) {
                return fail(
                    error, readerErrorCode(payloadReader),
                    "scene pressure cell mimetic pressure size is invalid");
            }
            std::span<const std::uint8_t> mimeticStateBytes;
            if (!payloadReader.bytes(
                    mimeticStateSize, mimeticStateBytes)) {
                return fail(
                    error, readerErrorCode(payloadReader),
                    "scene pressure cell mimetic pressure payload is truncated");
            }
            const auto topology =
                validator->rebuildMimeticPressureAuditTopology(
                    candidate.coupling.structure);
            SceneFluidMimeticPressureState state;
            SceneFluidMimeticPressureStatePersistenceError mimeticError;
            if (!deserializeSceneFluidMimeticPressureState(
                    mimeticStateBytes, topology.controlCells,
                    topology.fullTraceSystem,
                    topology.condensedTraceSystem, state,
                    &mimeticError, limits.mimeticPressureState)) {
                return fail(
                    error, mimeticStateErrorCode(mimeticError.code),
                    "scene pressure cell mimetic pressure state: "
                        + mimeticError.message);
            }
            candidate.coupling.mimeticPressureState = std::move(state);
        } else {
            candidate.coupling.mimeticPressureState.reset();
        }
        if (!payloadReader.atEnd()) {
            return fail(
                error,
                ScenePressureCellCheckpointPersistenceErrorCode::TrailingData,
                "scene pressure cell payload has trailing data");
        }

        validator->restore(candidate);
        checkpoint = std::move(candidate);
        return true;
    } catch (const std::bad_alloc&) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell checkpoint allocation exceeded available memory");
    } catch (const std::exception& exception) {
        return fail(
            error,
            ScenePressureCellCheckpointPersistenceErrorCode::InvalidData,
            exception.what());
    }
}

} // namespace simwing::fsi
