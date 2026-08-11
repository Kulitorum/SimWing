#include "scene_fluid_mimetic_pressure_flow.h"

#include "scene_fluid_mimetic_condensed_trace_system.h"
#include "scene_fluid_mimetic_pressure_source.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
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

bool finite(const double value) noexcept {
    return std::isfinite(value);
}

std::size_t storageBytes(
    const SceneFluidMimeticCorrectedTraceFlow& flow) {
    if (flow.traces.size()
            > std::numeric_limits<std::size_t>::max()
                / sizeof(SceneFluidMimeticCorrectedTrace)
        || flow.correctedContinuityResidualsCubicMetersPerSecond.size()
            > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        throw std::length_error(
            "scene fluid mimetic corrected-flow storage overflows");
    }
    const std::size_t traceBytes = flow.traces.size()
        * sizeof(SceneFluidMimeticCorrectedTrace);
    const std::size_t residualBytes =
        flow.correctedContinuityResidualsCubicMetersPerSecond.size()
        * sizeof(double);
    if (traceBytes
        > std::numeric_limits<std::size_t>::max() - residualBytes) {
        throw std::length_error(
            "scene fluid mimetic corrected-flow storage overflows");
    }
    return traceBytes + residualBytes;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticCorrectedTraceFlow& flow) {
    Fingerprint fingerprint;
    fingerprint.integer(flow.version);
    fingerprint.integer(flow.pressureAuditFingerprint);
    fingerprint.integer(flow.mimeticControlCellFingerprint);
    fingerprint.integer(flow.fullTraceSystemFingerprint);
    fingerprint.integer(flow.condensedTraceSystemFingerprint);
    fingerprint.integer(flow.predictedTraceFlowFingerprint);
    fingerprint.integer(flow.pressureSourceFingerprint);
    fingerprint.integer(flow.pressureStateFingerprint);
    fingerprint.integer(flow.pressureFaceLinkFingerprint);
    fingerprint.integer(flow.openingPatchFingerprint);
    fingerprint.integer(flow.structureDefinitionFingerprint);
    fingerprint.integer(flow.acceptedStepCount);
    fingerprint.real(flow.simulationTimeSeconds);
    fingerprint.real(flow.densityKgPerCubicMeter);
    fingerprint.real(flow.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(flow.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        flow.cartesianTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        flow.authoredOpeningTraceCount));
    fingerprint.real(
        flow.maximumAbsolutePressureCorrectionCubicMetersPerSecond);
    fingerprint.real(flow.maximumAbsoluteCorrectedFlowCubicMetersPerSecond);
    fingerprint.real(
        flow.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond);
    fingerprint.real(flow.correctedContinuityToleranceCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint8_t>(flow.finite));
    fingerprint.integer(static_cast<std::uint8_t>(flow.accepted));
    fingerprint.integer(static_cast<std::uint64_t>(flow.traces.size()));
    for (const auto& trace : flow.traces) {
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.sharedTraceOrdinal));
        fingerprint.integer(static_cast<std::uint64_t>(trace.traceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(trace.sourceIndex));
        fingerprint.integer(trace.stableId);
        fingerprint.integer(static_cast<std::uint8_t>(trace.kind));
        fingerprint.integer(trace.sourceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.minusControlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.plusControlCellIndex));
        fingerprint.real(
            trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            trace.pressureCorrectionVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        flow.correctedContinuityResidualsCubicMetersPerSecond.size()));
    for (const double residual :
         flow.correctedContinuityResidualsCubicMetersPerSecond) {
        fingerprint.real(residual);
    }
    return fingerprint.value();
}

const SceneFluidMimeticHalfFace& orientedMinusHalfFace(
    const SceneFluidMimeticPredictedTraceFlow& predicted,
    const SceneFluidMimeticTraceSystem& system,
    const SceneFluidMimeticControlCellSet& cells) {
    if (predicted.traceIndex >= system.traces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow trace is unavailable");
    }
    const auto& trace = system.traces[predicted.traceIndex];
    const SceneFluidMimeticHalfFace* found = nullptr;
    for (std::size_t offset = 0; offset < trace.incidenceCount; ++offset) {
        const auto& incidence = system.incidences[
            trace.firstIncidence + offset];
        if (incidence.controlCellIndex != predicted.minusControlCellIndex) {
            continue;
        }
        if (incidence.halfFaceIndex >= cells.halfFaces.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow incidence is invalid");
        }
        const auto& halfFace = cells.halfFaces[incidence.halfFaceIndex];
        if (found != nullptr
            || halfFace.side
                != SceneFluidMimeticHalfFaceSide::MinusOrNegative) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow orientation is ambiguous");
        }
        found = &halfFace;
    }
    if (found == nullptr) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow minus incidence is missing");
    }
    return *found;
}

double axisComponent(const Vec3 value,
                     const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return value.x;
    case fluid::GridFaceAxis::Y:
        return value.y;
    case fluid::GridFaceAxis::Z:
        return value.z;
    }
    throw std::invalid_argument(
        "scene fluid mimetic MAC-collapse axis is invalid");
}

} // namespace

SceneFluidMimeticCorrectedTraceFlow
correctSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint) {
    validateSceneFluidMimeticPressureAuditEndpointIntegrity(endpoint);
    const auto& epoch = endpoint.pressureEpoch;
    const auto& state = epoch.acceptedPressureState;
    if (!epoch.diagnostics.accepted || state.fingerprint == 0) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow pressure is not accepted");
    }

    const auto integratedSources =
        sceneFluidMimeticIntegratedCellSources(endpoint.pressureSources);
    const auto fullRightHandSide =
        buildSceneFluidMimeticTraceRightHandSide(
            endpoint.fullTraceSystem, integratedSources);
    std::vector<double> reducedTracePascals(
        endpoint.condensedTraceSystem.traces.size(), 0.0);
    for (const auto& trace : state.traces) {
        if (trace.reducedTraceIndex >= reducedTracePascals.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow pressure trace is invalid");
        }
        reducedTracePascals[trace.reducedTraceIndex] = trace.pressurePascals;
    }
    const auto fullTracePascals = reconstructSceneFluidMimeticFullTraces(
        endpoint.condensedTraceSystem, endpoint.fullTraceSystem,
        fullRightHandSide, reducedTracePascals);
    const auto evaluation = evaluateSceneFluidMimeticTraceSystem(
        endpoint.fullTraceSystem, fullTracePascals, integratedSources);
    if (evaluation.cellScalars.size() != state.controls.size()
        || evaluation.halfFaceIntegratedOutwardFluxes.size()
            != endpoint.controlCells.halfFaces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow reconstruction is incomplete");
    }
    for (std::size_t index = 0; index < state.controls.size(); ++index) {
        if (state.controls[index].controlCellIndex != index
            || state.controls[index].pressurePascals
                != evaluation.cellScalars[index]) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow pressure state changed");
        }
    }

    SceneFluidMimeticCorrectedTraceFlow result;
    result.pressureAuditFingerprint = endpoint.fingerprint;
    result.mimeticControlCellFingerprint = endpoint.controlCells.fingerprint;
    result.fullTraceSystemFingerprint = endpoint.fullTraceSystem.fingerprint;
    result.condensedTraceSystemFingerprint =
        endpoint.condensedTraceSystem.fingerprint;
    result.predictedTraceFlowFingerprint =
        endpoint.predictedTraceFlows.fingerprint;
    result.pressureSourceFingerprint = endpoint.pressureSources.fingerprint;
    result.pressureStateFingerprint = state.fingerprint;
    result.pressureFaceLinkFingerprint =
        endpoint.predictedTraceFlows.pressureFaceLinkFingerprint;
    result.openingPatchFingerprint =
        endpoint.controlCells.openingPatchFingerprint;
    result.structureDefinitionFingerprint =
        endpoint.structureDefinitionFingerprint;
    result.acceptedStepCount = endpoint.acceptedStepCount;
    result.simulationTimeSeconds = endpoint.simulationTimeSeconds;
    result.densityKgPerCubicMeter =
        endpoint.pressureSources.settings.densityKgPerCubicMeter;
    result.timeStepSeconds =
        endpoint.pressureSources.settings.timeStepSeconds;
    const double correctionScale = result.timeStepSeconds
        / result.densityKgPerCubicMeter;
    std::vector<CompensatedSum> correctedNetOutward(
        endpoint.controlCells.controlCells.size());
    result.traces.reserve(endpoint.predictedTraceFlows.traces.size());
    for (const auto& predicted : endpoint.predictedTraceFlows.traces) {
        const auto& halfFace = orientedMinusHalfFace(
            predicted, endpoint.fullTraceSystem, endpoint.controlCells);
        if (halfFace.halfFaceIndex
            >= evaluation.halfFaceIntegratedOutwardFluxes.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow half-face is invalid");
        }
        const double correction = correctionScale
            * evaluation.halfFaceIntegratedOutwardFluxes[
                halfFace.halfFaceIndex];
        const double corrected =
            predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond
            + correction;
        if (!finite(correction) || !finite(corrected)) {
            throw std::overflow_error(
                "scene fluid mimetic corrected-flow value is non-finite");
        }
        result.traces.push_back({
            predicted.sharedTraceOrdinal,
            predicted.traceIndex,
            halfFace.sourceIndex,
            predicted.stableId,
            predicted.kind,
            predicted.sourceStableId,
            predicted.componentIndex,
            predicted.minusControlCellIndex,
            predicted.plusControlCellIndex,
            predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            correction,
            corrected,
        });
        correctedNetOutward[predicted.minusControlCellIndex].add(corrected);
        correctedNetOutward[predicted.plusControlCellIndex].add(-corrected);
        result.maximumAbsolutePressureCorrectionCubicMetersPerSecond =
            std::max(
                result.maximumAbsolutePressureCorrectionCubicMetersPerSecond,
                std::abs(correction));
        result.maximumAbsoluteCorrectedFlowCubicMetersPerSecond = std::max(
            result.maximumAbsoluteCorrectedFlowCubicMetersPerSecond,
            std::abs(corrected));
        if (predicted.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            ++result.cartesianTraceCount;
        } else if (predicted.kind
                   == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            ++result.authoredOpeningTraceCount;
        } else {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow shared kind is invalid");
        }
    }
    result.correctedContinuityResidualsCubicMetersPerSecond.resize(
        endpoint.controlCells.controlCells.size(), 0.0);
    double continuityReference = 1.0;
    for (std::size_t index = 0;
         index < endpoint.controlCells.controlCells.size(); ++index) {
        const auto& source = endpoint.pressureSources.controls[index];
        if (source.controlCellIndex != index) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow source order changed");
        }
        const double residual =
            source.geometryVolumeChangeRateCubicMetersPerSecond
            + correctedNetOutward[index].value();
        if (!finite(residual)) {
            throw std::overflow_error(
                "scene fluid mimetic corrected-flow residual is non-finite");
        }
        result.correctedContinuityResidualsCubicMetersPerSecond[index] =
            residual;
        result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond,
                std::abs(residual));
        continuityReference = std::max({
            continuityReference,
            std::abs(source.geometryVolumeChangeRateCubicMetersPerSecond),
            std::abs(source.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond),
        });
    }
    result.correctedContinuityToleranceCubicMetersPerSecond =
        correctionScale
            * (epoch.diagnostics.pressureSolve.maximumCellConservationResidual
               + static_cast<double>(
                     endpoint.controlCells.maximumHalfFaceCountPerControl)
                   * epoch.diagnostics.pressureSolve
                         .reconstructedFullResidualMaximumPascalsMeters)
        + 64.0 * std::numeric_limits<double>::epsilon()
            * continuityReference;
    result.ownedStorageBytes = storageBytes(result);
    result.finite = finite(
        result.maximumAbsolutePressureCorrectionCubicMetersPerSecond)
        && finite(result.maximumAbsoluteCorrectedFlowCubicMetersPerSecond)
        && finite(
            result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond)
        && finite(result.correctedContinuityToleranceCubicMetersPerSecond);
    result.accepted = result.finite
        && result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond
            <= result.correctedContinuityToleranceCubicMetersPerSecond;
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(result);
    return result;
}

void validateSceneFluidMimeticCorrectedTraceFlowIntegrity(
    const SceneFluidMimeticCorrectedTraceFlow& flow) {
    if (flow.version != sceneFluidMimeticCorrectedTraceFlowVersion
        || flow.fingerprint == 0 || flow.pressureAuditFingerprint == 0
        || flow.mimeticControlCellFingerprint == 0
        || flow.fullTraceSystemFingerprint == 0
        || flow.condensedTraceSystemFingerprint == 0
        || flow.predictedTraceFlowFingerprint == 0
        || flow.pressureSourceFingerprint == 0
        || flow.pressureStateFingerprint == 0
        || flow.pressureFaceLinkFingerprint == 0
        || flow.openingPatchFingerprint == 0
        || flow.structureDefinitionFingerprint == 0
        || !finite(flow.simulationTimeSeconds)
        || flow.simulationTimeSeconds < 0.0
        || !finite(flow.densityKgPerCubicMeter)
        || !(flow.densityKgPerCubicMeter > 0.0)
        || !finite(flow.timeStepSeconds) || !(flow.timeStepSeconds > 0.0)
        || flow.traces.empty()
        || flow.correctedContinuityResidualsCubicMetersPerSecond.empty()
        || flow.cartesianTraceCount + flow.authoredOpeningTraceCount
            != flow.traces.size()
        || flow.ownedStorageBytes != storageBytes(flow)) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow metadata is invalid");
    }
    std::size_t cartesianCount = 0;
    std::size_t openingCount = 0;
    double maximumCorrection = 0.0;
    double maximumFlow = 0.0;
    for (std::size_t index = 0; index < flow.traces.size(); ++index) {
        const auto& trace = flow.traces[index];
        if (trace.sharedTraceOrdinal != index || trace.stableId == 0
            || trace.sourceStableId == 0
            || trace.minusControlCellIndex
                >= flow.correctedContinuityResidualsCubicMetersPerSecond.size()
            || trace.plusControlCellIndex
                >= flow.correctedContinuityResidualsCubicMetersPerSecond.size()
            || trace.minusControlCellIndex == trace.plusControlCellIndex
            || !finite(
                trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond)
            || !finite(
                trace.pressureCorrectionVolumeFlowRateCubicMetersPerSecond)
            || !finite(
                trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond)
            || trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond
                != trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond
                    + trace.pressureCorrectionVolumeFlowRateCubicMetersPerSecond) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow trace is invalid");
        }
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            ++cartesianCount;
        } else if (trace.kind
                   == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            ++openingCount;
        } else {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow trace kind is invalid");
        }
        maximumCorrection = std::max(
            maximumCorrection,
            std::abs(
                trace.pressureCorrectionVolumeFlowRateCubicMetersPerSecond));
        maximumFlow = std::max(
            maximumFlow,
            std::abs(trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond));
    }
    double maximumResidual = 0.0;
    for (const double residual :
         flow.correctedContinuityResidualsCubicMetersPerSecond) {
        if (!finite(residual)) {
            throw std::invalid_argument(
                "scene fluid mimetic corrected-flow residual is invalid");
        }
        maximumResidual = std::max(maximumResidual, std::abs(residual));
    }
    const bool finiteExpected = finite(maximumCorrection)
        && finite(maximumFlow) && finite(maximumResidual)
        && finite(flow.correctedContinuityToleranceCubicMetersPerSecond);
    const bool acceptedExpected = finiteExpected
        && maximumResidual
            <= flow.correctedContinuityToleranceCubicMetersPerSecond;
    if (cartesianCount != flow.cartesianTraceCount
        || openingCount != flow.authoredOpeningTraceCount
        || maximumCorrection
            != flow.maximumAbsolutePressureCorrectionCubicMetersPerSecond
        || maximumFlow != flow.maximumAbsoluteCorrectedFlowCubicMetersPerSecond
        || maximumResidual
            != flow.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond
        || flow.correctedContinuityToleranceCubicMetersPerSecond < 0.0
        || flow.finite != finiteExpected || flow.accepted != acceptedExpected
        || flow.fingerprint != productFingerprint(flow)) {
        throw std::invalid_argument(
            "scene fluid mimetic corrected-flow diagnostics are invalid");
    }
}

SceneFluidMimeticMacVelocityCollapse
collapseSceneFluidMimeticCorrectedMacVelocity(
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const fluid::PeriodicCartesianGrid& grid) {
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(correctedFlow);
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningGridPatchIntegrity(openingPatches);
    if (!correctedFlow.accepted
        || correctedFlow.mimeticControlCellFingerprint
            != controlCells.fingerprint
        || correctedFlow.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || correctedFlow.openingPatchFingerprint != openingPatches.fingerprint
        || faceLinks.openingPatchFingerprint != openingPatches.fingerprint
        || controlCells.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || controlCells.openingPatchFingerprint != openingPatches.fingerprint
        || correctedFlow.structureDefinitionFingerprint
            != faceLinks.structureDefinitionFingerprint
        || correctedFlow.acceptedStepCount != faceLinks.acceptedStepCount
        || correctedFlow.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || faceLinks.cellCounts != grid.cellCounts()
        || faceLinks.lowerMeters != grid.lowerMeters()
        || faceLinks.upperMeters != grid.upperMeters()) {
        throw std::invalid_argument(
            "scene fluid mimetic MAC-collapse source is foreign");
    }

    std::vector<const SceneFluidMimeticCorrectedTrace*> linkTraces(
        faceLinks.links.size(), nullptr);
    std::map<std::uint64_t, const SceneFluidMimeticCorrectedTrace*>
        embeddedTracesByPatchId;
    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*> patchesById;
    for (const auto& patch : openingPatches.patches) {
        if (patch.stableId == 0
            || !patchesById.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene fluid mimetic MAC-collapse opening identity is invalid");
        }
    }

    SceneFluidMimeticMacVelocityCollapse result(grid);
    auto& diagnostics = result.diagnostics;
    diagnostics.correctedTraceFlowFingerprint = correctedFlow.fingerprint;
    diagnostics.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    diagnostics.openingPatchFingerprint = openingPatches.fingerprint;
    diagnostics.acceptedStepCount = correctedFlow.acceptedStepCount;
    diagnostics.simulationTimeSeconds = correctedFlow.simulationTimeSeconds;
    diagnostics.faceCount = faceLinks.faces.size();
    diagnostics.cartesianTraceCount = correctedFlow.cartesianTraceCount;
    diagnostics.authoredOpeningTraceCount =
        correctedFlow.authoredOpeningTraceCount;

    for (const auto& trace : correctedFlow.traces) {
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            if (trace.sourceIndex >= faceLinks.links.size()
                || linkTraces[trace.sourceIndex] != nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse Cartesian trace is invalid");
            }
            const auto& link = faceLinks.links[trace.sourceIndex];
            if (link.linkIndex != trace.sourceIndex
                || link.stableId != trace.sourceStableId
                || link.geometryKind
                    != SceneFluidPressureLinkGeometryKind::CartesianFace) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse link binding is invalid");
            }
            linkTraces[trace.sourceIndex] = &trace;
        } else {
            if (trace.sourceIndex
                >= controlCells.openingTraceGroups.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse embedded trace is invalid");
            }
            const auto& group = controlCells.openingTraceGroups[
                trace.sourceIndex];
            if (group.stableId != trace.sourceStableId
                || !finite(
                    trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond
                    + group.surfaceSweepRateCubicMetersPerSecond)) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse embedded opening is foreign");
            }
            for (std::size_t offset = 0;
                 offset < group.patchCount; ++offset) {
                const std::size_t patchIndex =
                    controlCells.openingTracePatchIndices[
                        group.firstPatchIndex + offset];
                if (patchIndex >= openingPatches.patches.size()) {
                    throw std::invalid_argument(
                        "scene fluid mimetic MAC-collapse embedded opening member is invalid");
                }
                const auto& patch = openingPatches.patches[patchIndex];
                if (patch.ownerKind
                        != SceneFluidOpeningPatchOwnerKind::Cell
                    || !embeddedTracesByPatchId
                            .emplace(patch.stableId, &trace).second) {
                    throw std::invalid_argument(
                        "scene fluid mimetic MAC-collapse embedded opening member is foreign");
                }
            }
            ++diagnostics.embeddedOpeningTraceCount;
        }
    }

    const auto absoluteFlow = [&](const SceneFluidPressureFaceLink& link,
                                  const SceneFluidPressureFace& face,
                                  const SceneFluidMimeticCorrectedTrace& trace) {
        double flow =
            trace.correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        if (link.kind == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            const auto found = patchesById.find(link.openingPatchStableId);
            if (found == patchesById.end()) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse face opening is missing");
            }
            const auto& patch = *found->second;
            const double normal = axisComponent(
                patch.unitNormalNegativeToPositive, face.axis);
            if (patch.openingId != link.openingId
                || patch.areaSquareMeters != link.areaSquareMeters
                || std::abs(std::abs(normal) - 1.0) > 1.0e-10) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse face opening is foreign");
            }
            flow += (normal > 0.0 ? 1.0 : -1.0)
                * patch.surfaceSweepRateCubicMetersPerSecond;
        }
        return flow;
    };

    for (const auto& face : faceLinks.faces) {
        if (!(face.faceAreaSquareMeters > 0.0) || face.linkCount == 0
            || face.firstLink > faceLinks.links.size()
            || face.linkCount > faceLinks.links.size() - face.firstLink) {
            throw std::invalid_argument(
                "scene fluid mimetic MAC-collapse face is unresolved");
        }
        if (face.linkCount > 1) {
            ++diagnostics.multiTraceFaceCount;
        }
        double totalFlow = 0.0;
        for (std::size_t offset = 0; offset < face.linkCount; ++offset) {
            const std::size_t linkIndex = face.firstLink + offset;
            const auto& link = faceLinks.links[linkIndex];
            const auto* trace = linkTraces[linkIndex];
            if (link.linkIndex != linkIndex || link.faceIndex != face.faceIndex
                || trace == nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse face binding is invalid");
            }
            totalFlow += absoluteFlow(link, face, *trace);
        }
        const double velocity = totalFlow / face.faceAreaSquareMeters;
        if (!finite(velocity)) {
            throw std::overflow_error(
                "scene fluid mimetic MAC-collapse velocity is non-finite");
        }
        diagnostics.maximumAbsoluteVelocityMetersPerSecond = std::max(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond,
            std::abs(velocity));
        diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond = std::max(
            diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond,
            std::abs(velocity * face.faceAreaSquareMeters - totalFlow));
        for (std::size_t offset = 0; offset < face.linkCount; ++offset) {
            const std::size_t linkIndex = face.firstLink + offset;
            const auto& link = faceLinks.links[linkIndex];
            const double subfaceVelocity = absoluteFlow(
                link, face, *linkTraces[linkIndex])
                / link.areaSquareMeters;
            diagnostics.maximumSubfaceVelocityDeviationMetersPerSecond =
                std::max(
                    diagnostics.maximumSubfaceVelocityDeviationMetersPerSecond,
                    std::abs(subfaceVelocity - velocity));
        }
        const std::size_t cellIndex = grid.cellIndex(face.i, face.j, face.k);
        switch (face.axis) {
        case fluid::GridFaceAxis::X:
            result.velocityMetersPerSecond.xFaces()[cellIndex] = velocity;
            break;
        case fluid::GridFaceAxis::Y:
            result.velocityMetersPerSecond.yFaces()[cellIndex] = velocity;
            break;
        case fluid::GridFaceAxis::Z:
            result.velocityMetersPerSecond.zFaces()[cellIndex] = velocity;
            break;
        }
    }
    for (const auto& link : faceLinks.links) {
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::CartesianFace) {
            if (linkTraces[link.linkIndex] == nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic MAC-collapse missed a Cartesian link");
            }
            continue;
        }
        const auto found = embeddedTracesByPatchId.find(
            link.openingPatchStableId);
        if (link.geometryKind
                != SceneFluidPressureLinkGeometryKind::EmbeddedOpening
            || link.kind != SceneFluidPressureFaceLinkKind::AuthoredOpening
            || link.faceIndex != invalidSceneFluidPressureFaceIndex
            || found == embeddedTracesByPatchId.end()) {
            throw std::invalid_argument(
                "scene fluid mimetic MAC-collapse embedded link is invalid");
        }
    }
    diagnostics.finite = fluid::isFinite(result.velocityMetersPerSecond)
        && finite(diagnostics.maximumAbsoluteVelocityMetersPerSecond)
        && finite(diagnostics.maximumSubfaceVelocityDeviationMetersPerSecond)
        && finite(diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene fluid mimetic MAC-collapse diagnostics are non-finite");
    }
    return result;
}

} // namespace simwing::fsi
