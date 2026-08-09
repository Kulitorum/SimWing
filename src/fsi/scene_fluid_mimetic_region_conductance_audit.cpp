#include "scene_fluid_mimetic_region_conductance_audit.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

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
        const double next = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) {
            correction_ += (sum_ - next) + value;
        } else {
            correction_ += (value - next) + sum_;
        }
        sum_ = next;
    }

    [[nodiscard]] double value() const noexcept {
        return sum_ + correction_;
    }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

struct OpeningHalfFaceReference {
    std::uint64_t traceStableId = 0;
    SceneFluidMimeticHalfFaceKind traceKind =
        SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace;
    std::size_t halfFaceIndex = 0;
    std::size_t controlCellIndex = 0;
    std::size_t otherControlCellIndex = 0;
    StableId regionId = invalidStableId;
    double areaSquareMeters = 0.0;
};

bool checkedMultiply(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

void validateSettings(
    const SceneFluidMimeticRegionConductanceAuditSettings& settings) {
    if (!std::isfinite(settings.terminalIntegratedTransferPascalsMeters)
        || !(settings.terminalIntegratedTransferPascalsMeters > 0.0)
        || !std::isfinite(
            settings.absoluteOpeningPairAreaToleranceSquareMeters)
        || settings.absoluteOpeningPairAreaToleranceSquareMeters < 0.0
        || !std::isfinite(settings.relativeOpeningPairAreaTolerance)
        || settings.relativeOpeningPairAreaTolerance < 0.0
        || !std::isfinite(
            settings.solve.absoluteResidualTolerancePascalsMeters)
        || settings.solve.absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(settings.solve.relativeResidualTolerance)
        || settings.solve.relativeResidualTolerance < 0.0
        || !std::isfinite(settings.solve
            .absoluteComponentCompatibilityTolerancePascalsMeters)
        || settings.solve
            .absoluteComponentCompatibilityTolerancePascalsMeters < 0.0
        || settings.solve.maximumIterations == 0) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance settings are invalid");
    }
}

std::size_t storageBytes(
    const SceneFluidMimeticRegionConductanceAudit& audit) {
    std::size_t openingBytes = 0;
    std::size_t responseBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t result = 0;
    if (!checkedMultiply(audit.openings.size(),
                         sizeof(SceneFluidMimeticRegionConductanceOpening),
                         openingBytes)
        || !checkedMultiply(audit.responses.size(),
                            sizeof(SceneFluidMimeticRegionConductanceResponse),
                            responseBytes)
        || !checkedMultiply(
            audit.solveDiagnostics.reducedTraceSolve.components.size(),
            sizeof(SceneFluidMimeticTraceSolveComponentDiagnostics),
            componentBytes)
        || !checkedAdd(openingBytes, responseBytes, result)
        || !checkedAdd(result, componentBytes, result)) {
        throw std::length_error(
            "scene fluid mimetic region-conductance storage overflows");
    }
    return result;
}

void fingerprintSolveDiagnostics(
    Fingerprint& fingerprint,
    const SceneFluidMimeticPressureSolveDiagnostics& diagnostics) {
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.accepted));
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.reconstructedFullResidualConverged));
    fingerprint.real(
        diagnostics.reconstructedFullResidualTolerancePascalsMeters);
    fingerprint.real(
        diagnostics.reconstructedFullResidualL2PascalsMeters);
    fingerprint.real(
        diagnostics.reconstructedFullResidualMaximumPascalsMeters);
    fingerprint.real(diagnostics.maximumCellConservationResidual);
    const auto& reduced = diagnostics.reducedTraceSolve;
    fingerprint.integer(static_cast<std::uint8_t>(reduced.compatible));
    fingerprint.integer(static_cast<std::uint8_t>(reduced.converged));
    fingerprint.integer(static_cast<std::uint8_t>(reduced.finite));
    fingerprint.integer(reduced.traceSystemFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(reduced.traceCount));
    fingerprint.integer(static_cast<std::uint64_t>(reduced.componentCount));
    fingerprint.integer(static_cast<std::uint64_t>(reduced.iterationCount));
    fingerprint.real(
        reduced.maximumAbsoluteComponentCompatibilityPascalsMeters);
    fingerprint.real(reduced.initialResidualL2PascalsMeters);
    fingerprint.real(reduced.finalResidualL2PascalsMeters);
    fingerprint.real(reduced.finalResidualMaximumPascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(reduced.components.size()));
    for (const auto& component : reduced.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.traceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeTraceIndex));
        fingerprint.real(component.rightHandSideSumPascalsMeters);
        fingerprint.real(component.compatibilityCorrectionPascalsMeters);
        fingerprint.real(component.traceGaugeBeforePascals);
        fingerprint.real(component.traceGaugeAfterPascals);
    }
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticRegionConductanceAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.mimeticControlCellFingerprint);
    fingerprint.integer(audit.fullTraceSystemFingerprint);
    fingerprint.integer(audit.condensedTraceSystemFingerprint);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprint.integer(audit.acceptedStepCount);
    fingerprint.real(audit.simulationTimeSeconds);
    fingerprint.real(
        audit.settings.terminalIntegratedTransferPascalsMeters);
    fingerprint.real(
        audit.settings.absoluteOpeningPairAreaToleranceSquareMeters);
    fingerprint.real(audit.settings.relativeOpeningPairAreaTolerance);
    fingerprint.real(
        audit.settings.solve.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(audit.settings.solve.relativeResidualTolerance);
    fingerprint.real(audit.settings.solve
        .absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.settings.solve.maximumIterations));
    fingerprint.integer(static_cast<std::uint64_t>(audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(audit.componentCount));
    fingerprint.integer(audit.lowerTerminalRegionId);
    fingerprint.integer(audit.upperTerminalRegionId);
    fingerprint.real(audit.openingAreaSquareMeters);
    fingerprint.real(audit.maximumOpeningPairAreaMismatchSquareMeters);
    fingerprint.real(audit.achievedIntegratedTransferPascalsMeters);
    fingerprint.real(audit.lowerTerminalIntegratedSourcePascalsMeters);
    fingerprint.real(audit.upperTerminalIntegratedSourcePascalsMeters);
    fingerprint.real(audit.componentIntegratedSourcePascalsMeters);
    fingerprint.real(audit.sourceL2PascalsMeters);
    fingerprint.real(audit.maximumAbsoluteSourcePascalsMeters);
    fingerprint.real(audit.pressureL2Pascals);
    fingerprint.real(audit.sourcePressureWorkPascalsSquaredMeters);
    fingerprint.real(audit.effectiveTerminalPressureDifferencePascals);
    fingerprint.real(audit.conductanceMeters);
    fingerprintSolveDiagnostics(fingerprint, audit.solveDiagnostics);
    fingerprint.integer(static_cast<std::uint64_t>(audit.openings.size()));
    for (const auto& opening : audit.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(opening.openingIndex));
        fingerprint.integer(opening.traceStableId);
        fingerprint.integer(static_cast<std::uint8_t>(opening.traceKind));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.lowerHalfFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.upperHalfFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.lowerControlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.upperControlCellIndex));
        fingerprint.real(opening.areaSquareMeters);
        fingerprint.real(opening.areaMismatchSquareMeters);
        fingerprint.real(opening.integratedTransferPascalsMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(audit.responses.size()));
    for (const auto& response : audit.responses) {
        fingerprint.integer(static_cast<std::uint64_t>(
            response.controlCellIndex));
        fingerprint.integer(response.stableId);
        fingerprint.integer(response.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            response.componentIndex));
        fingerprint.real(response.integratedSourcePascalsMeters);
        fingerprint.real(response.gaugeAlignedPressurePascals);
    }
    return fingerprint.value();
}

} // namespace

SceneFluidMimeticRegionConductanceAudit
auditSceneFluidMimeticRegionConductance(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullTraceSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedTraceSystem,
    const SceneFluidMimeticRegionConductanceAuditSettings& settings,
    const SceneFluidMimeticRegionConductanceAuditLimits& limits) {
    validateSettings(settings);
    validateSceneFluidMimeticTraceSystem(fullTraceSystem, controlCells);
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedTraceSystem, fullTraceSystem);
    const std::size_t controlCount = controlCells.controlCells.size();
    const auto isOpeningHalfFace = [&](const auto& halfFace) {
        return halfFace.kind
                != SceneFluidMimeticHalfFaceKind::MaterialWall
            && halfFace.otherControlVolumeIndex < controlCount
            && controlCells.controlCells[halfFace.controlVolumeIndex].regionId
                != controlCells.controlCells[
                    halfFace.otherControlVolumeIndex].regionId;
    };
    const std::size_t openingHalfFaceCount = static_cast<std::size_t>(
        std::ranges::count_if(controlCells.halfFaces, isOpeningHalfFace));
    if (controlCount == 0 || fullTraceSystem.componentCount != 1
        || openingHalfFaceCount == 0 || openingHalfFaceCount % 2 != 0) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance topology is unsupported");
    }
    const std::size_t openingCount = openingHalfFaceCount / 2;
    if (controlCount > limits.maximumControlCells
        || openingCount > limits.maximumOpeningTraces) {
        throw std::length_error(
            "scene fluid mimetic region-conductance count limit exceeded");
    }
    std::size_t staticBytes = 0;
    std::size_t openingBytes = 0;
    if (!checkedMultiply(
            controlCount,
            sizeof(SceneFluidMimeticRegionConductanceResponse), staticBytes)
        || !checkedMultiply(
            openingCount,
            sizeof(SceneFluidMimeticRegionConductanceOpening), openingBytes)
        || !checkedAdd(staticBytes, openingBytes, staticBytes)
        || staticBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic region-conductance byte limit exceeded");
    }

    std::vector<StableId> regionIds;
    regionIds.reserve(controlCount);
    for (const auto& control : controlCells.controlCells) {
        if (control.componentIndex != 0) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance component is invalid");
        }
        regionIds.push_back(control.regionId);
    }
    std::ranges::sort(regionIds);
    regionIds.erase(
        std::unique(regionIds.begin(), regionIds.end()), regionIds.end());
    if (regionIds.size() != 2 || regionIds[0] == invalidStableId) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance requires two regions");
    }

    std::vector<OpeningHalfFaceReference> halfFaces;
    halfFaces.reserve(openingHalfFaceCount);
    for (const auto& halfFace : controlCells.halfFaces) {
        if (!isOpeningHalfFace(halfFace)) {
            continue;
        }
        if (halfFace.controlVolumeIndex >= controlCount
            || halfFace.otherControlVolumeIndex >= controlCount) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance opening is unbound");
        }
        halfFaces.push_back({
            halfFace.traceStableId,
            halfFace.kind,
            halfFace.halfFaceIndex,
            halfFace.controlVolumeIndex,
            halfFace.otherControlVolumeIndex,
            controlCells.controlCells[
                halfFace.controlVolumeIndex].regionId,
            halfFace.areaSquareMeters,
        });
    }
    std::ranges::sort(
        halfFaces, [](const auto& first, const auto& second) {
            return std::tie(first.traceStableId, first.traceKind,
                            first.regionId,
                            first.halfFaceIndex)
                < std::tie(second.traceStableId, second.traceKind,
                           second.regionId,
                           second.halfFaceIndex);
        });

    SceneFluidMimeticRegionConductanceAudit result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.fullTraceSystemFingerprint = fullTraceSystem.fingerprint;
    result.condensedTraceSystemFingerprint =
        condensedTraceSystem.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.settings = settings;
    result.componentCount = 1;
    result.lowerTerminalRegionId = regionIds[0];
    result.upperTerminalRegionId = regionIds[1];
    result.openings.reserve(openingCount);
    CompensatedSum openingArea;
    for (std::size_t index = 0; index < halfFaces.size(); index += 2) {
        const auto& first = halfFaces[index];
        const auto& second = halfFaces[index + 1];
        if (first.traceStableId == 0
            || first.traceStableId != second.traceStableId
            || first.traceKind != second.traceKind
            || (index + 2 < halfFaces.size()
                && halfFaces[index + 2].traceStableId
                    == first.traceStableId)
            || first.regionId != regionIds[0]
            || second.regionId != regionIds[1]
            || first.otherControlCellIndex != second.controlCellIndex
            || second.otherControlCellIndex != first.controlCellIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance opening pair is invalid");
        }
        const double mismatch = std::abs(
            first.areaSquareMeters - second.areaSquareMeters);
        const double area = 0.5
            * (first.areaSquareMeters + second.areaSquareMeters);
        const double tolerance =
            settings.absoluteOpeningPairAreaToleranceSquareMeters
            + settings.relativeOpeningPairAreaTolerance * area;
        if (!std::isfinite(area) || !(area > 0.0)
            || !std::isfinite(mismatch) || mismatch > tolerance) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance opening area is invalid");
        }
        result.openings.push_back({
            result.openings.size(), first.traceStableId,
            first.traceKind,
            first.halfFaceIndex, second.halfFaceIndex,
            first.controlCellIndex, second.controlCellIndex,
            area, mismatch, 0.0,
        });
        openingArea.add(area);
        result.maximumOpeningPairAreaMismatchSquareMeters = std::max(
            result.maximumOpeningPairAreaMismatchSquareMeters, mismatch);
    }
    result.openingAreaSquareMeters = openingArea.value();
    if (!std::isfinite(result.openingAreaSquareMeters)
        || !(result.openingAreaSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance opening area is empty");
    }

    std::vector<double> sources(controlCount, 0.0);
    double assignedTransfer = 0.0;
    for (std::size_t index = 0; index < result.openings.size(); ++index) {
        auto& opening = result.openings[index];
        opening.integratedTransferPascalsMeters =
            index + 1 == result.openings.size()
            ? settings.terminalIntegratedTransferPascalsMeters
                - assignedTransfer
            : settings.terminalIntegratedTransferPascalsMeters
                * opening.areaSquareMeters
                / result.openingAreaSquareMeters;
        if (!std::isfinite(opening.integratedTransferPascalsMeters)
            || !(opening.integratedTransferPascalsMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance transfer is invalid");
        }
        assignedTransfer += opening.integratedTransferPascalsMeters;
        sources[opening.lowerControlCellIndex] +=
            opening.integratedTransferPascalsMeters;
        sources[opening.upperControlCellIndex] -=
            opening.integratedTransferPascalsMeters;
    }

    const std::vector<double> zeroWarmStart(
        condensedTraceSystem.traces.size(), 0.0);
    const auto solve = solveSceneFluidMimeticPressureSystem(
        condensedTraceSystem, fullTraceSystem, sources, zeroWarmStart,
        settings.solve);
    if (!solve.diagnostics.accepted
        || solve.evaluation.cellScalars.size() != controlCount) {
        throw std::runtime_error(
            "scene fluid mimetic region-conductance solve did not accept");
    }
    auto pressure = solve.evaluation.cellScalars;
    CompensatedSum pressureMean;
    for (const double value : pressure) pressureMean.add(value);
    const double mean = pressureMean.value()
        / static_cast<double>(pressure.size());
    for (double& value : pressure) value -= mean;

    result.solveDiagnostics = solve.diagnostics;
    result.responses.reserve(controlCount);
    CompensatedSum lowerSource;
    CompensatedSum upperSource;
    CompensatedSum componentSource;
    CompensatedSum sourceSquared;
    CompensatedSum pressureSquared;
    CompensatedSum sourceWork;
    for (std::size_t index = 0; index < controlCount; ++index) {
        const auto& control = controlCells.controlCells[index];
        result.responses.push_back({
            index, control.stableId, control.regionId,
            control.componentIndex, sources[index], pressure[index],
        });
        if (control.regionId == result.lowerTerminalRegionId) {
            lowerSource.add(sources[index]);
        } else {
            upperSource.add(sources[index]);
        }
        componentSource.add(sources[index]);
        sourceSquared.add(sources[index] * sources[index]);
        pressureSquared.add(pressure[index] * pressure[index]);
        sourceWork.add(sources[index] * pressure[index]);
        result.maximumAbsoluteSourcePascalsMeters = std::max(
            result.maximumAbsoluteSourcePascalsMeters,
            std::abs(sources[index]));
    }
    result.lowerTerminalIntegratedSourcePascalsMeters = lowerSource.value();
    result.upperTerminalIntegratedSourcePascalsMeters = upperSource.value();
    result.componentIntegratedSourcePascalsMeters = componentSource.value();
    result.achievedIntegratedTransferPascalsMeters = 0.5
        * (std::abs(result.lowerTerminalIntegratedSourcePascalsMeters)
           + std::abs(result.upperTerminalIntegratedSourcePascalsMeters));
    result.sourceL2PascalsMeters = std::sqrt(
        std::max(0.0, sourceSquared.value()));
    result.pressureL2Pascals = std::sqrt(
        std::max(0.0, pressureSquared.value()));
    result.sourcePressureWorkPascalsSquaredMeters = sourceWork.value();
    result.effectiveTerminalPressureDifferencePascals =
        result.sourcePressureWorkPascalsSquaredMeters
        / result.achievedIntegratedTransferPascalsMeters;
    result.conductanceMeters =
        result.achievedIntegratedTransferPascalsMeters
        * result.achievedIntegratedTransferPascalsMeters
        / result.sourcePressureWorkPascalsSquaredMeters;
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic region-conductance byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticRegionConductanceAuditIntegrity(result);
    return result;
}

void validateSceneFluidMimeticRegionConductanceAuditIntegrity(
    const SceneFluidMimeticRegionConductanceAudit& audit) {
    validateSettings(audit.settings);
    if (audit.version != sceneFluidMimeticRegionConductanceAuditVersion
        || audit.fingerprint == 0
        || audit.mimeticControlCellFingerprint == 0
        || audit.fullTraceSystemFingerprint == 0
        || audit.condensedTraceSystemFingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || !std::isfinite(audit.simulationTimeSeconds)
        || audit.componentCount != 1
        || audit.lowerTerminalRegionId == invalidStableId
        || audit.upperTerminalRegionId == invalidStableId
        || audit.lowerTerminalRegionId >= audit.upperTerminalRegionId
        || audit.openings.empty() || audit.responses.empty()
        || !std::isfinite(audit.openingAreaSquareMeters)
        || !(audit.openingAreaSquareMeters > 0.0)
        || !std::isfinite(
            audit.maximumOpeningPairAreaMismatchSquareMeters)
        || audit.maximumOpeningPairAreaMismatchSquareMeters < 0.0
        || !std::isfinite(
            audit.achievedIntegratedTransferPascalsMeters)
        || !(audit.achievedIntegratedTransferPascalsMeters > 0.0)
        || !std::isfinite(
            audit.lowerTerminalIntegratedSourcePascalsMeters)
        || !std::isfinite(
            audit.upperTerminalIntegratedSourcePascalsMeters)
        || !std::isfinite(audit.componentIntegratedSourcePascalsMeters)
        || !std::isfinite(audit.sourceL2PascalsMeters)
        || audit.sourceL2PascalsMeters < 0.0
        || !std::isfinite(audit.maximumAbsoluteSourcePascalsMeters)
        || audit.maximumAbsoluteSourcePascalsMeters < 0.0
        || !std::isfinite(audit.pressureL2Pascals)
        || audit.pressureL2Pascals < 0.0
        || !std::isfinite(
            audit.sourcePressureWorkPascalsSquaredMeters)
        || !(audit.sourcePressureWorkPascalsSquaredMeters > 0.0)
        || !std::isfinite(
            audit.effectiveTerminalPressureDifferencePascals)
        || !(audit.effectiveTerminalPressureDifferencePascals > 0.0)
        || !std::isfinite(audit.conductanceMeters)
        || !(audit.conductanceMeters > 0.0)
        || audit.ownedStorageBytes != storageBytes(audit)
        || !audit.solveDiagnostics.accepted
        || !audit.solveDiagnostics.reconstructedFullResidualConverged
        || !audit.solveDiagnostics.reducedTraceSolve.compatible
        || !audit.solveDiagnostics.reducedTraceSolve.converged
        || !audit.solveDiagnostics.reducedTraceSolve.finite
        || audit.solveDiagnostics.reducedTraceSolve.componentCount != 1
        || audit.solveDiagnostics.reducedTraceSolve.components.size() != 1) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance integrity is invalid");
    }

    std::vector<double> expectedSources(audit.responses.size(), 0.0);
    CompensatedSum openingArea;
    double assignedTransfer = 0.0;
    double maximumAreaMismatch = 0.0;
    std::uint64_t previousTraceStableId = 0;
    for (std::size_t index = 0; index < audit.openings.size(); ++index) {
        const auto& opening = audit.openings[index];
        const double expectedTransfer = index + 1 == audit.openings.size()
            ? audit.settings.terminalIntegratedTransferPascalsMeters
                - assignedTransfer
            : audit.settings.terminalIntegratedTransferPascalsMeters
                * opening.areaSquareMeters / audit.openingAreaSquareMeters;
        const double areaTolerance =
            audit.settings.absoluteOpeningPairAreaToleranceSquareMeters
            + audit.settings.relativeOpeningPairAreaTolerance
                * opening.areaSquareMeters;
        if (opening.openingIndex != index
            || opening.traceStableId <= previousTraceStableId
            || (opening.traceKind
                    != SceneFluidMimeticHalfFaceKind::CartesianTrace
                && opening.traceKind
                    != SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace)
            || opening.lowerControlCellIndex >= audit.responses.size()
            || opening.upperControlCellIndex >= audit.responses.size()
            || opening.lowerControlCellIndex == opening.upperControlCellIndex
            || !std::isfinite(opening.areaSquareMeters)
            || !(opening.areaSquareMeters > 0.0)
            || !std::isfinite(opening.areaMismatchSquareMeters)
            || opening.areaMismatchSquareMeters < 0.0
            || opening.areaMismatchSquareMeters > areaTolerance
            || opening.integratedTransferPascalsMeters != expectedTransfer
            || !(opening.integratedTransferPascalsMeters > 0.0)
            || audit.responses[opening.lowerControlCellIndex].regionId
                != audit.lowerTerminalRegionId
            || audit.responses[opening.upperControlCellIndex].regionId
                != audit.upperTerminalRegionId) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance opening is invalid");
        }
        openingArea.add(opening.areaSquareMeters);
        maximumAreaMismatch = std::max(
            maximumAreaMismatch, opening.areaMismatchSquareMeters);
        assignedTransfer += opening.integratedTransferPascalsMeters;
        expectedSources[opening.lowerControlCellIndex] +=
            opening.integratedTransferPascalsMeters;
        expectedSources[opening.upperControlCellIndex] -=
            opening.integratedTransferPascalsMeters;
        previousTraceStableId = opening.traceStableId;
    }

    CompensatedSum lowerSource;
    CompensatedSum upperSource;
    CompensatedSum componentSource;
    CompensatedSum sourceSquared;
    CompensatedSum pressureSquared;
    CompensatedSum sourceWork;
    double maximumSource = 0.0;
    for (std::size_t index = 0; index < audit.responses.size(); ++index) {
        const auto& response = audit.responses[index];
        if (response.controlCellIndex != index || response.stableId == 0
            || response.componentIndex != 0
            || (response.regionId != audit.lowerTerminalRegionId
                && response.regionId != audit.upperTerminalRegionId)
            || response.integratedSourcePascalsMeters
                != expectedSources[index]
            || !std::isfinite(response.gaugeAlignedPressurePascals)) {
            throw std::invalid_argument(
                "scene fluid mimetic region-conductance response is invalid");
        }
        if (response.regionId == audit.lowerTerminalRegionId) {
            lowerSource.add(response.integratedSourcePascalsMeters);
        } else {
            upperSource.add(response.integratedSourcePascalsMeters);
        }
        componentSource.add(response.integratedSourcePascalsMeters);
        sourceSquared.add(response.integratedSourcePascalsMeters
            * response.integratedSourcePascalsMeters);
        pressureSquared.add(response.gaugeAlignedPressurePascals
            * response.gaugeAlignedPressurePascals);
        sourceWork.add(response.integratedSourcePascalsMeters
            * response.gaugeAlignedPressurePascals);
        maximumSource = std::max(
            maximumSource,
            std::abs(response.integratedSourcePascalsMeters));
    }
    const double achievedTransfer = 0.5
        * (std::abs(lowerSource.value()) + std::abs(upperSource.value()));
    const double work = sourceWork.value();
    if (audit.openingAreaSquareMeters != openingArea.value()
        || audit.maximumOpeningPairAreaMismatchSquareMeters
            != maximumAreaMismatch
        || audit.achievedIntegratedTransferPascalsMeters
            != achievedTransfer
        || audit.lowerTerminalIntegratedSourcePascalsMeters
            != lowerSource.value()
        || audit.upperTerminalIntegratedSourcePascalsMeters
            != upperSource.value()
        || audit.componentIntegratedSourcePascalsMeters
            != componentSource.value()
        || audit.sourceL2PascalsMeters
            != std::sqrt(std::max(0.0, sourceSquared.value()))
        || audit.maximumAbsoluteSourcePascalsMeters != maximumSource
        || audit.pressureL2Pascals
            != std::sqrt(std::max(0.0, pressureSquared.value()))
        || audit.sourcePressureWorkPascalsSquaredMeters != work
        || !std::isfinite(work) || !(work > 0.0)
        || audit.effectiveTerminalPressureDifferencePascals
            != work / achievedTransfer
        || audit.conductanceMeters
            != achievedTransfer * achievedTransfer / work
        || !std::isfinite(audit.conductanceMeters)
        || !(audit.conductanceMeters > 0.0)) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance summary is invalid");
    }
    if (audit.fingerprint != productFingerprint(audit)) {
        throw std::invalid_argument(
            "scene fluid mimetic region-conductance fingerprint is invalid");
    }
}

} // namespace simwing::fsi
