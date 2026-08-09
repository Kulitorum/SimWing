#include "scene_fluid_pressure_operator_response_audit.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double corrected = value - correction_;
        const double updated = sum_ + corrected;
        correction_ = (updated - sum_) - corrected;
        sum_ = updated;
    }

    [[nodiscard]] double value() const noexcept { return sum_; }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

bool finiteSolveSettings(const SceneFluidPressureSolveSettings& settings) {
    return std::isfinite(settings.absoluteResidualTolerancePascalsMeters)
        && settings.absoluteResidualTolerancePascalsMeters >= 0.0
        && std::isfinite(settings.relativeResidualTolerance)
        && settings.relativeResidualTolerance >= 0.0
        && (settings.absoluteResidualTolerancePascalsMeters > 0.0
            || settings.relativeResidualTolerance > 0.0)
        && std::isfinite(
            settings.absoluteComponentCompatibilityTolerancePascalsMeters)
        && settings
                .absoluteComponentCompatibilityTolerancePascalsMeters
            >= 0.0
        && settings.maximumIterations > 0;
}

bool finiteSolveSettings(
    const SceneFluidMimeticTraceSolveSettings& settings) {
    return std::isfinite(settings.absoluteResidualTolerancePascalsMeters)
        && settings.absoluteResidualTolerancePascalsMeters >= 0.0
        && std::isfinite(settings.relativeResidualTolerance)
        && settings.relativeResidualTolerance >= 0.0
        && (settings.absoluteResidualTolerancePascalsMeters > 0.0
            || settings.relativeResidualTolerance > 0.0)
        && std::isfinite(
            settings.absoluteComponentCompatibilityTolerancePascalsMeters)
        && settings
                .absoluteComponentCompatibilityTolerancePascalsMeters
            >= 0.0
        && settings.maximumIterations > 0;
}

void validateSettings(
    const SceneFluidPressureOperatorResponseAuditSettings& settings) {
    if (!std::isfinite(settings.manufacturedPressureL2Pascals)
        || !(settings.manufacturedPressureL2Pascals > 0.0)
        || !finiteSolveSettings(settings.graphSolve)
        || !finiteSolveSettings(settings.shadowSolve)) {
        throw std::invalid_argument(
            "invalid scene fluid pressure operator-response settings");
    }
}

std::size_t storageBytesForCounts(
    const std::size_t modeCount,
    const std::size_t responseCount) {
    if (modeCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidPressureOperatorResponseModeDiagnostics)
        || responseCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidPressureOperatorResponseRecord)) {
        throw std::length_error(
            "scene fluid pressure operator-response storage overflows");
    }
    const std::size_t modeBytes = modeCount
        * sizeof(SceneFluidPressureOperatorResponseModeDiagnostics);
    const std::size_t responseBytes = responseCount
        * sizeof(SceneFluidPressureOperatorResponseRecord);
    if (responseBytes > std::numeric_limits<std::size_t>::max() - modeBytes) {
        throw std::length_error(
            "scene fluid pressure operator-response storage overflows");
    }
    return modeBytes + responseBytes;
}

std::size_t storageBytes(
    const SceneFluidPressureOperatorResponseAudit& audit) {
    return storageBytesForCounts(audit.modes.size(), audit.responses.size());
}

void fingerprintGraphSettings(
    Fingerprint& fingerprint,
    const SceneFluidPressureSolveSettings& settings) {
    fingerprint.real(settings.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(settings.relativeResidualTolerance);
    fingerprint.real(
        settings.absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        settings.maximumIterations));
}

void fingerprintShadowSettings(
    Fingerprint& fingerprint,
    const SceneFluidMimeticTraceSolveSettings& settings) {
    fingerprint.real(settings.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(settings.relativeResidualTolerance);
    fingerprint.real(
        settings.absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        settings.maximumIterations));
}

std::uint64_t productFingerprint(
    const SceneFluidPressureOperatorResponseAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.graphPressureOperatorFingerprint);
    fingerprint.integer(audit.mimeticControlCellFingerprint);
    fingerprint.integer(audit.fullTraceSystemFingerprint);
    fingerprint.integer(audit.condensedTraceSystemFingerprint);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprint.integer(audit.acceptedStepCount);
    fingerprint.real(audit.simulationTimeSeconds);
    fingerprint.real(audit.settings.manufacturedPressureL2Pascals);
    fingerprintGraphSettings(fingerprint, audit.settings.graphSolve);
    fingerprintShadowSettings(fingerprint, audit.settings.shadowSolve);
    fingerprint.integer(static_cast<std::uint64_t>(audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(audit.componentCount));
    fingerprint.integer(static_cast<std::uint8_t>(
        audit.includesAcceptedSource));
    fingerprint.integer(static_cast<std::uint64_t>(audit.modes.size()));
    for (const auto& mode : audit.modes) {
        fingerprint.integer(static_cast<std::uint64_t>(mode.modeIndex));
        fingerprint.integer(static_cast<std::uint8_t>(mode.kind));
        fingerprint.integer(static_cast<std::uint64_t>(mode.firstResponse));
        fingerprint.integer(static_cast<std::uint64_t>(mode.responseCount));
        fingerprint.real(mode.sourceL2PascalsMeters);
        fingerprint.real(mode.maximumAbsoluteSourcePascalsMeters);
        fingerprint.real(mode.graphPressureL2Pascals);
        fingerprint.real(mode.shadowPressureL2Pascals);
        fingerprint.real(mode.pressureDotProductPascalsSquared);
        fingerprint.real(mode.bestFitShadowPressureScale);
        fingerprint.real(mode.pressureCosineSimilarity);
        fingerprint.real(mode.bestFitShapeResidualL2Pascals);
        fingerprint.real(mode.relativeBestFitShapeResidualL2);
        fingerprint.real(mode.maximumAbsoluteBestFitShapeResidualPascals);
        fingerprint.integer(static_cast<std::uint64_t>(
            mode.graphIterationCount));
        fingerprint.real(mode.graphFinalResidualL2PascalsMeters);
        fingerprint.real(mode.graphFinalResidualMaximumPascalsMeters);
        fingerprint.integer(static_cast<std::uint64_t>(
            mode.shadowIterationCount));
        fingerprint.real(mode.shadowFinalResidualL2PascalsMeters);
        fingerprint.real(mode.shadowFinalResidualMaximumPascalsMeters);
        fingerprint.real(mode.shadowMaximumCellConservationResidual);
        fingerprint.integer(static_cast<std::uint8_t>(mode.finite));
    }
    fingerprint.integer(static_cast<std::uint64_t>(audit.responses.size()));
    for (const auto& response : audit.responses) {
        fingerprint.integer(static_cast<std::uint64_t>(response.modeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            response.controlVolumeIndex));
        fingerprint.integer(response.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            response.componentIndex));
        fingerprint.real(response.integratedSourcePascalsMeters);
        fingerprint.real(response.graphGaugeAlignedPressurePascals);
        fingerprint.real(response.shadowGaugeAlignedPressurePascals);
        fingerprint.real(
            response.shadowMinusBestFitGraphPressurePascals);
    }
    return fingerprint.value();
}

void centerByComponent(
    const SceneFluidPressureOperator& graphOperator,
    std::vector<double>& values) {
    for (const auto& component : graphOperator.components) {
        CompensatedSum sum;
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            sum.add(values[graphOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]]);
        }
        const double mean = sum.value()
            / static_cast<double>(component.controlVolumeCount);
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            values[graphOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]] -= mean;
        }
    }
}

double l2(const std::span<const double> values) {
    CompensatedSum squared;
    for (const double value : values) squared.add(value * value);
    return std::sqrt(std::max(0.0, squared.value()));
}

double manufacturedValue(
    const SceneFluidMimeticControlCell& control,
    const SceneFluidPressureOperatorResponseModeKind kind) {
    const auto& centroid = control.centroidMeters;
    switch (kind) {
    case SceneFluidPressureOperatorResponseModeKind::CoordinateX:
        return centroid.x;
    case SceneFluidPressureOperatorResponseModeKind::CoordinateY:
        return centroid.y;
    case SceneFluidPressureOperatorResponseModeKind::CoordinateZ:
        return centroid.z;
    case SceneFluidPressureOperatorResponseModeKind::MixedCoordinate:
        return centroid.x * centroid.y + 0.37 * centroid.z;
    case SceneFluidPressureOperatorResponseModeKind::StableIdPattern:
        return static_cast<double>(control.stableId % 31U) - 15.0;
    case SceneFluidPressureOperatorResponseModeKind::AcceptedSource:
        break;
    }
    throw std::logic_error(
        "accepted pressure source is not a manufactured mode");
}

std::vector<double> manufacturedSource(
    const SceneFluidPressureOperator& graphOperator,
    const SceneFluidMimeticControlCellSet& controls,
    const SceneFluidPressureOperatorResponseModeKind kind,
    const double targetL2Pascals) {
    std::vector<double> pressure(controls.controlCells.size(), 0.0);
    for (const auto& control : controls.controlCells) {
        pressure[control.controlVolumeIndex] =
            manufacturedValue(control, kind);
    }
    centerByComponent(graphOperator, pressure);
    const double norm = l2(pressure);
    if (!std::isfinite(norm) || !(norm > 0.0)) {
        throw std::invalid_argument(
            "scene fluid pressure operator-response mode is constant");
    }
    for (double& value : pressure) value *= targetL2Pascals / norm;
    return applySceneFluidPressureOperator(graphOperator, pressure);
}

SceneFluidPressureOperatorResponseModeDiagnostics summarizeMode(
    const std::span<const SceneFluidPressureOperatorResponseRecord> responses,
    SceneFluidPressureOperatorResponseModeDiagnostics result) {
    CompensatedSum sourceSquared;
    CompensatedSum graphSquared;
    CompensatedSum shadowSquared;
    CompensatedSum cross;
    for (const auto& response : responses) {
        sourceSquared.add(response.integratedSourcePascalsMeters
            * response.integratedSourcePascalsMeters);
        graphSquared.add(response.graphGaugeAlignedPressurePascals
            * response.graphGaugeAlignedPressurePascals);
        shadowSquared.add(response.shadowGaugeAlignedPressurePascals
            * response.shadowGaugeAlignedPressurePascals);
        cross.add(response.graphGaugeAlignedPressurePascals
            * response.shadowGaugeAlignedPressurePascals);
        result.maximumAbsoluteSourcePascalsMeters = std::max(
            result.maximumAbsoluteSourcePascalsMeters,
            std::abs(response.integratedSourcePascalsMeters));
    }
    result.sourceL2PascalsMeters = std::sqrt(
        std::max(0.0, sourceSquared.value()));
    result.graphPressureL2Pascals = std::sqrt(
        std::max(0.0, graphSquared.value()));
    result.shadowPressureL2Pascals = std::sqrt(
        std::max(0.0, shadowSquared.value()));
    result.pressureDotProductPascalsSquared = cross.value();
    result.bestFitShadowPressureScale = graphSquared.value() > 0.0
        ? cross.value() / graphSquared.value() : 0.0;
    const double normProduct = result.graphPressureL2Pascals
        * result.shadowPressureL2Pascals;
    result.pressureCosineSimilarity = normProduct > 0.0
        ? cross.value() / normProduct : 0.0;
    CompensatedSum residualSquared;
    for (const auto& response : responses) {
        const double residual = response.shadowGaugeAlignedPressurePascals
            - result.bestFitShadowPressureScale
                * response.graphGaugeAlignedPressurePascals;
        residualSquared.add(residual * residual);
        result.maximumAbsoluteBestFitShapeResidualPascals = std::max(
            result.maximumAbsoluteBestFitShapeResidualPascals,
            std::abs(residual));
    }
    result.bestFitShapeResidualL2Pascals = std::sqrt(
        std::max(0.0, residualSquared.value()));
    result.relativeBestFitShapeResidualL2 = result.shadowPressureL2Pascals > 0.0
        ? result.bestFitShapeResidualL2Pascals
            / result.shadowPressureL2Pascals
        : 0.0;
    result.finite = result.sourceL2PascalsMeters > 0.0
        && result.graphPressureL2Pascals > 0.0
        && result.shadowPressureL2Pascals > 0.0
        && std::isfinite(result.maximumAbsoluteSourcePascalsMeters)
        && std::isfinite(result.pressureDotProductPascalsSquared)
        && std::isfinite(result.bestFitShadowPressureScale)
        && std::isfinite(result.pressureCosineSimilarity)
        && std::isfinite(result.bestFitShapeResidualL2Pascals)
        && std::isfinite(result.relativeBestFitShapeResidualL2)
        && std::isfinite(
            result.maximumAbsoluteBestFitShapeResidualPascals)
        && std::isfinite(result.graphFinalResidualL2PascalsMeters)
        && std::isfinite(result.graphFinalResidualMaximumPascalsMeters)
        && std::isfinite(result.shadowFinalResidualL2PascalsMeters)
        && std::isfinite(result.shadowFinalResidualMaximumPascalsMeters)
        && std::isfinite(result.shadowMaximumCellConservationResidual);
    return result;
}

} // namespace

SceneFluidPressureOperatorResponseAudit
auditSceneFluidPressureOperatorResponses(
    const SceneFluidPressureOperator& graphOperator,
    const SceneFluidMimeticControlCellSet& mimeticControlCells,
    const SceneFluidMimeticTraceSystem& fullTraceSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedTraceSystem,
    const std::span<const double> acceptedIntegratedSourcePascalsMeters,
    const SceneFluidPressureOperatorResponseAuditSettings& settings,
    const SceneFluidPressureOperatorResponseAuditLimits& limits) {
    validateSettings(settings);
    validateSceneFluidMimeticTraceSystem(
        fullTraceSystem, mimeticControlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedTraceSystem, fullTraceSystem);
    const std::size_t controlCount = mimeticControlCells.controlCells.size();
    const std::size_t componentCount = graphOperator.components.size();
    const bool includesAcceptedSource =
        !acceptedIntegratedSourcePascalsMeters.empty();
    const std::size_t modeCount = includesAcceptedSource ? 6 : 5;
    if (graphOperator.fingerprint == 0
        || controlCount == 0
        || graphOperator.rows.size() != controlCount
        || graphOperator.pressureControlVolumeFingerprint
            != mimeticControlCells.pressureControlVolumeFingerprint
        || graphOperator.structureDefinitionFingerprint
            != mimeticControlCells.structureDefinitionFingerprint
        || graphOperator.acceptedStepCount
            != mimeticControlCells.acceptedStepCount
        || graphOperator.simulationTimeSeconds
            != mimeticControlCells.simulationTimeSeconds
        || fullTraceSystem.mimeticControlCellFingerprint
            != mimeticControlCells.fingerprint
        || condensedTraceSystem.fullTraceSystemFingerprint
            != fullTraceSystem.fingerprint
        || componentCount == 0
        || fullTraceSystem.componentCount != componentCount
        || condensedTraceSystem.componentCount != componentCount
        || (includesAcceptedSource
            && acceptedIntegratedSourcePascalsMeters.size() != controlCount)) {
        throw std::invalid_argument(
            "scene fluid pressure operator-response topology is invalid");
    }
    for (std::size_t index = 0; index < controlCount; ++index) {
        const auto& row = graphOperator.rows[index];
        const auto& control = mimeticControlCells.controlCells[index];
        if (row.rowIndex != index || row.controlVolumeIndex != index
            || control.controlCellIndex != index
            || control.controlVolumeIndex != index
            || row.stableId != control.stableId
            || row.componentIndex != control.componentIndex) {
            throw std::invalid_argument(
                "scene fluid pressure operator-response control is misbound");
        }
    }
    if (controlCount > limits.maximumControlVolumes
        || modeCount > limits.maximumModes
        || (controlCount != 0
            && modeCount > limits.maximumResponseRecords / controlCount)) {
        throw std::length_error(
            "scene fluid pressure operator-response count limit exceeded");
    }
    const std::size_t responseCount = modeCount * controlCount;
    if (storageBytesForCounts(modeCount, responseCount)
        > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure operator-response byte limit exceeded");
    }
    if (!std::ranges::all_of(
            acceptedIntegratedSourcePascalsMeters,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "scene fluid pressure operator-response source is non-finite");
    }

    SceneFluidPressureOperatorResponseAudit result;
    result.graphPressureOperatorFingerprint = graphOperator.fingerprint;
    result.mimeticControlCellFingerprint = mimeticControlCells.fingerprint;
    result.fullTraceSystemFingerprint = fullTraceSystem.fingerprint;
    result.condensedTraceSystemFingerprint = condensedTraceSystem.fingerprint;
    result.structureDefinitionFingerprint =
        graphOperator.structureDefinitionFingerprint;
    result.acceptedStepCount = graphOperator.acceptedStepCount;
    result.simulationTimeSeconds = graphOperator.simulationTimeSeconds;
    result.settings = settings;
    result.controlVolumeCount = controlCount;
    result.componentCount = componentCount;
    result.includesAcceptedSource = includesAcceptedSource;
    result.modes.reserve(modeCount);
    result.responses.reserve(responseCount);

    std::vector<SceneFluidPressureOperatorResponseModeKind> kinds;
    kinds.reserve(modeCount);
    if (includesAcceptedSource) {
        kinds.push_back(
            SceneFluidPressureOperatorResponseModeKind::AcceptedSource);
    }
    kinds.insert(kinds.end(), {
        SceneFluidPressureOperatorResponseModeKind::CoordinateX,
        SceneFluidPressureOperatorResponseModeKind::CoordinateY,
        SceneFluidPressureOperatorResponseModeKind::CoordinateZ,
        SceneFluidPressureOperatorResponseModeKind::MixedCoordinate,
        SceneFluidPressureOperatorResponseModeKind::StableIdPattern,
    });

    for (std::size_t modeIndex = 0; modeIndex < kinds.size(); ++modeIndex) {
        const auto kind = kinds[modeIndex];
        std::vector<double> source = kind
                == SceneFluidPressureOperatorResponseModeKind::AcceptedSource
            ? std::vector<double>(
                acceptedIntegratedSourcePascalsMeters.begin(),
                acceptedIntegratedSourcePascalsMeters.end())
            : manufacturedSource(
                graphOperator, mimeticControlCells, kind,
                settings.manufacturedPressureL2Pascals);
        std::vector<double> graphPressure(controlCount, 0.0);
        const auto graphDiagnostics = solveSceneFluidPressureSystem(
            graphOperator, source, graphPressure, settings.graphSolve);
        const std::vector<double> zeroWarmStart(
            condensedTraceSystem.traces.size(), 0.0);
        const auto shadow = solveSceneFluidMimeticPressureSystem(
            condensedTraceSystem, fullTraceSystem, source, zeroWarmStart,
            settings.shadowSolve);
        if (!graphDiagnostics.compatible || !graphDiagnostics.converged
            || !graphDiagnostics.finite || !shadow.diagnostics.accepted
            || shadow.evaluation.cellScalars.size() != controlCount) {
            throw std::runtime_error(
                "scene fluid pressure operator-response solve did not accept");
        }
        auto shadowPressure = shadow.evaluation.cellScalars;
        centerByComponent(graphOperator, graphPressure);
        centerByComponent(graphOperator, shadowPressure);

        SceneFluidPressureOperatorResponseModeDiagnostics mode;
        mode.modeIndex = modeIndex;
        mode.kind = kind;
        mode.firstResponse = result.responses.size();
        mode.responseCount = controlCount;
        mode.graphIterationCount = graphDiagnostics.iterationCount;
        mode.graphFinalResidualL2PascalsMeters =
            graphDiagnostics.finalResidualL2PascalsMeters;
        mode.graphFinalResidualMaximumPascalsMeters =
            graphDiagnostics.finalResidualMaximumPascalsMeters;
        mode.shadowIterationCount =
            shadow.diagnostics.reducedTraceSolve.iterationCount;
        mode.shadowFinalResidualL2PascalsMeters =
            shadow.diagnostics.reducedTraceSolve
                .finalResidualL2PascalsMeters;
        mode.shadowFinalResidualMaximumPascalsMeters =
            shadow.diagnostics.reducedTraceSolve
                .finalResidualMaximumPascalsMeters;
        mode.shadowMaximumCellConservationResidual =
            shadow.diagnostics.maximumCellConservationResidual;
        for (std::size_t index = 0; index < controlCount; ++index) {
            result.responses.push_back({
                modeIndex,
                index,
                graphOperator.rows[index].stableId,
                graphOperator.rows[index].componentIndex,
                source[index],
                graphPressure[index],
                shadowPressure[index],
                0.0,
            });
        }
        const auto responseSpan = std::span(
            result.responses).subspan(mode.firstResponse, controlCount);
        mode = summarizeMode(responseSpan, mode);
        for (auto& response : responseSpan) {
            response.shadowMinusBestFitGraphPressurePascals =
                response.shadowGaugeAlignedPressurePascals
                - mode.bestFitShadowPressureScale
                    * response.graphGaugeAlignedPressurePascals;
        }
        if (!mode.finite) {
            throw std::overflow_error(
                "scene fluid pressure operator-response mode is non-finite");
        }
        result.modes.push_back(mode);
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure operator-response byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidPressureOperatorResponseAuditIntegrity(result);
    return result;
}

void validateSceneFluidPressureOperatorResponseAuditIntegrity(
    const SceneFluidPressureOperatorResponseAudit& audit) {
    validateSettings(audit.settings);
    if (audit.version != sceneFluidPressureOperatorResponseAuditVersion
        || audit.fingerprint == 0
        || audit.graphPressureOperatorFingerprint == 0
        || audit.mimeticControlCellFingerprint == 0
        || audit.fullTraceSystemFingerprint == 0
        || audit.condensedTraceSystemFingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || !std::isfinite(audit.simulationTimeSeconds)
        || audit.controlVolumeCount == 0 || audit.componentCount == 0
        || audit.modes.size() != (audit.includesAcceptedSource ? 6U : 5U)
        || audit.controlVolumeCount
            > std::numeric_limits<std::size_t>::max() / audit.modes.size()
        || audit.responses.size()
            != audit.controlVolumeCount * audit.modes.size()
        || audit.ownedStorageBytes != storageBytes(audit)) {
        throw std::invalid_argument(
            "scene fluid pressure operator-response integrity is invalid");
    }
    for (std::size_t modeIndex = 0;
         modeIndex < audit.modes.size(); ++modeIndex) {
        const auto& mode = audit.modes[modeIndex];
        const auto expectedKind = audit.includesAcceptedSource
            ? static_cast<SceneFluidPressureOperatorResponseModeKind>(
                static_cast<std::uint8_t>(
                    SceneFluidPressureOperatorResponseModeKind::AcceptedSource)
                + static_cast<std::uint8_t>(modeIndex))
            : static_cast<SceneFluidPressureOperatorResponseModeKind>(
                static_cast<std::uint8_t>(
                    SceneFluidPressureOperatorResponseModeKind::CoordinateX)
                + static_cast<std::uint8_t>(modeIndex));
        if (mode.modeIndex != modeIndex || mode.kind != expectedKind
            || mode.firstResponse != modeIndex * audit.controlVolumeCount
            || mode.responseCount != audit.controlVolumeCount
            || !mode.finite) {
            throw std::invalid_argument(
                "scene fluid pressure operator-response mode is invalid");
        }
        const auto responses = std::span(audit.responses).subspan(
            mode.firstResponse, mode.responseCount);
        for (std::size_t index = 0; index < responses.size(); ++index) {
            const auto& response = responses[index];
            if (response.modeIndex != modeIndex
                || response.controlVolumeIndex != index
                || response.stableId == 0
                || response.componentIndex >= audit.componentCount
                || !std::isfinite(response.integratedSourcePascalsMeters)
                || !std::isfinite(
                    response.graphGaugeAlignedPressurePascals)
                || !std::isfinite(
                    response.shadowGaugeAlignedPressurePascals)) {
                throw std::invalid_argument(
                    "scene fluid pressure operator-response record is invalid");
            }
        }
        auto expected = mode;
        expected.sourceL2PascalsMeters = 0.0;
        expected.maximumAbsoluteSourcePascalsMeters = 0.0;
        expected.graphPressureL2Pascals = 0.0;
        expected.shadowPressureL2Pascals = 0.0;
        expected.pressureDotProductPascalsSquared = 0.0;
        expected.bestFitShadowPressureScale = 0.0;
        expected.pressureCosineSimilarity = 0.0;
        expected.bestFitShapeResidualL2Pascals = 0.0;
        expected.relativeBestFitShapeResidualL2 = 0.0;
        expected.maximumAbsoluteBestFitShapeResidualPascals = 0.0;
        expected.finite = false;
        expected = summarizeMode(responses, expected);
        if (expected != mode) {
            throw std::invalid_argument(
                "scene fluid pressure operator-response summary is invalid");
        }
        for (const auto& response : responses) {
            if (response.shadowMinusBestFitGraphPressurePascals
                != response.shadowGaugeAlignedPressurePascals
                    - mode.bestFitShadowPressureScale
                        * response.graphGaugeAlignedPressurePascals) {
                throw std::invalid_argument(
                    "scene fluid pressure operator-response residual is invalid");
            }
        }
    }
    if (audit.fingerprint != productFingerprint(audit)) {
        throw std::invalid_argument(
            "scene fluid pressure operator-response fingerprint is invalid");
    }
}

} // namespace simwing::fsi
