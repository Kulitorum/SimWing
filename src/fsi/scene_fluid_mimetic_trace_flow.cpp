#include "scene_fluid_mimetic_trace_flow.h"

#include "scene_fluid_mimetic_pressure_flow.h"
#include "scene_fluid_region_transport.h"

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

    template<typename Enum>
    void enumeration(const Enum value) {
        integer(static_cast<std::underlying_type_t<Enum>>(value));
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

std::size_t storageBytesForCounts(const std::size_t traceCount,
                                  const std::size_t componentCount) {
    if (traceCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticPredictedTraceFlow)
        || componentCount > std::numeric_limits<std::size_t>::max()
            / sizeof(double)) {
        throw std::length_error(
            "scene fluid mimetic trace-flow storage overflows");
    }
    const std::size_t traceBytes = traceCount
        * sizeof(SceneFluidMimeticPredictedTraceFlow);
    const std::size_t componentBytes = componentCount * sizeof(double);
    if (componentBytes
        > std::numeric_limits<std::size_t>::max() - traceBytes) {
        throw std::length_error(
            "scene fluid mimetic trace-flow storage overflows");
    }
    return traceBytes + componentBytes;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticTraceFlowPrediction& prediction) {
    Fingerprint fingerprint;
    fingerprint.integer(prediction.version);
    fingerprint.integer(prediction.mimeticControlCellFingerprint);
    fingerprint.integer(prediction.mimeticTraceSystemFingerprint);
    fingerprint.integer(prediction.pressureFaceLinkFingerprint);
    fingerprint.integer(prediction.openingFluxFingerprint);
    fingerprint.integer(prediction.velocityFingerprint);
    fingerprint.integer(prediction.regionWallExchangeFingerprint);
    fingerprint.integer(prediction.structureDefinitionFingerprint);
    fingerprint.integer(prediction.acceptedStepCount);
    fingerprint.real(prediction.simulationTimeSeconds);
    fingerprint.real(prediction.sourceDensityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.componentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.cartesianTraceCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.authoredOpeningTraceCount));
    fingerprint.real(prediction
        .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond);
    fingerprint.real(prediction
        .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint64_t>(prediction.traces.size()));
    for (const auto& trace : prediction.traces) {
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.sharedTraceOrdinal));
        fingerprint.integer(static_cast<std::uint64_t>(trace.traceIndex));
        fingerprint.integer(trace.stableId);
        fingerprint.enumeration(trace.kind);
        fingerprint.integer(trace.sourceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.minusControlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.plusControlCellIndex));
        fingerprint.real(
            trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.componentBalanceResidualsCubicMetersPerSecond.size()));
    for (const double value :
         prediction.componentBalanceResidualsCubicMetersPerSecond) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

std::uint64_t continuationFingerprint(
    const SceneFluidMimeticTraceFlowContinuation& continuation) {
    Fingerprint fingerprint;
    fingerprint.integer(continuation.version);
    fingerprint.integer(
        continuation.previousCorrectedTraceFlowFingerprint);
    fingerprint.integer(
        continuation.previousBaselinePredictionFingerprint);
    fingerprint.integer(
        continuation.currentBaselinePredictionFingerprint);
    fingerprint.integer(continuation.mimeticControlCellFingerprint);
    fingerprint.integer(continuation.mimeticTraceSystemFingerprint);
    fingerprint.integer(continuation.pressureFaceLinkFingerprint);
    fingerprint.integer(continuation.structureDefinitionFingerprint);
    fingerprint.integer(continuation.acceptedStepCount);
    fingerprint.real(continuation.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.ownedStorageBytes));
    fingerprint.real(
        continuation.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond);
    fingerprint.real(
        continuation.maximumAbsoluteBulkIncrementCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint8_t>(continuation.finite));
    fingerprint.integer(continuation.prediction.fingerprint);
    return fingerprint.value();
}

std::uint64_t regionTransportPredictionFingerprint(
    const SceneFluidMimeticRegionTransportFlowPrediction& prediction) {
    Fingerprint fingerprint;
    fingerprint.integer(prediction.version);
    fingerprint.integer(prediction.regionTransportFingerprint);
    fingerprint.integer(
        prediction.sourceCorrectedTraceFlowFingerprint);
    fingerprint.integer(
        prediction.currentBulkBaselinePredictionFingerprint);
    fingerprint.integer(prediction.mimeticControlCellFingerprint);
    fingerprint.integer(prediction.mimeticTraceSystemFingerprint);
    fingerprint.integer(prediction.pressureFaceLinkFingerprint);
    fingerprint.integer(prediction.structureDefinitionFingerprint);
    fingerprint.integer(prediction.acceptedStepCount);
    fingerprint.real(prediction.simulationTimeSeconds);
    fingerprint.real(prediction.transportTargetSimulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.ownedStorageBytes));
    fingerprint.real(
        prediction
            .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint8_t>(prediction.finite));
    fingerprint.integer(prediction.prediction.fingerprint);
    return fingerprint.value();
}

double faceVelocity(const SceneFluidPressureFace& face,
                    const fluid::PeriodicCartesianGrid& grid,
                    const fluid::MacVelocityField& velocity) {
    const std::size_t index = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case fluid::GridFaceAxis::X: return velocity.xFaces()[index];
    case fluid::GridFaceAxis::Y: return velocity.yFaces()[index];
    case fluid::GridFaceAxis::Z: return velocity.zFaces()[index];
    }
    throw std::invalid_argument(
        "scene fluid mimetic trace-flow face axis is invalid");
}

double orientedOpeningFlow(
    const SceneFluidPressureFaceLink& link,
    const SceneFluidOpeningFluxSample& sample) {
    if (sample.openingId != link.openingId
        || sample.areaSquareMeters != link.areaSquareMeters) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow opening identity is inconsistent");
    }
    if (sample.negativeSideRegionId == link.minusRegionId
        && sample.positiveSideRegionId == link.plusRegionId) {
        return sample.relativeVolumeFlowRateCubicMetersPerSecond;
    }
    if (sample.negativeSideRegionId == link.plusRegionId
        && sample.positiveSideRegionId == link.minusRegionId) {
        return -sample.relativeVolumeFlowRateCubicMetersPerSecond;
    }
    throw std::invalid_argument(
        "scene fluid mimetic trace-flow opening orientation is inconsistent");
}

double orientedOpeningSweep(
    const SceneFluidPressureFaceLink& link,
    const SceneFluidOpeningFluxSample& sample) {
    if (sample.openingId != link.openingId
        || sample.areaSquareMeters != link.areaSquareMeters) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow opening identity is inconsistent");
    }
    if (sample.negativeSideRegionId == link.minusRegionId
        && sample.positiveSideRegionId == link.plusRegionId) {
        return sample.surfaceSweepRateCubicMetersPerSecond;
    }
    if (sample.negativeSideRegionId == link.plusRegionId
        && sample.positiveSideRegionId == link.minusRegionId) {
        return -sample.surfaceSweepRateCubicMetersPerSecond;
    }
    throw std::invalid_argument(
        "scene fluid mimetic trace-flow opening orientation is inconsistent");
}

double dot(const fluid::Vector3& first, const fluid::Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

fluid::Vector3 averageVelocity(
    const SceneFluidRegionWallExchange& wallExchange,
    const std::size_t minusControlCellIndex,
    const std::size_t plusControlCellIndex) {
    const auto& minus = wallExchange.controlVolumes[
        minusControlCellIndex].velocityMetersPerSecond;
    const auto& plus = wallExchange.controlVolumes[
        plusControlCellIndex].velocityMetersPerSecond;
    return {
        0.5 * (minus.x + plus.x),
        0.5 * (minus.y + plus.y),
        0.5 * (minus.z + plus.z),
    };
}

fluid::Vector3 averageVelocity(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidRegionTransport& regionTransport,
    const std::size_t minusControlCellIndex,
    const std::size_t plusControlCellIndex) {
    const auto& minusCell = controlCells.controlCells[
        minusControlCellIndex];
    const auto& plusCell = controlCells.controlCells[
        plusControlCellIndex];
    const auto& minus = regionTransport.controlVolumes[
        minusCell.controlVolumeIndex].velocityMetersPerSecond;
    const auto& plus = regionTransport.controlVolumes[
        plusCell.controlVolumeIndex].velocityMetersPerSecond;
    return {
        0.5 * (minus.x + plus.x),
        0.5 * (minus.y + plus.y),
        0.5 * (minus.z + plus.z),
    };
}

struct OrientedTrace {
    const SceneFluidMimeticHalfFace* minus = nullptr;
    const SceneFluidMimeticHalfFace* plus = nullptr;
    std::size_t minusControlCellIndex = 0;
    std::size_t plusControlCellIndex = 0;
};

OrientedTrace orientTrace(
    const SceneFluidMimeticTrace& trace,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidMimeticControlCellSet& controlCells) {
    if (trace.incidenceCount != 2
        || trace.firstIncidence > traceSystem.incidences.size()
        || trace.incidenceCount
            > traceSystem.incidences.size() - trace.firstIncidence) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow shared incidence is invalid");
    }
    OrientedTrace result;
    for (std::size_t offset = 0; offset < trace.incidenceCount; ++offset) {
        const auto& incidence = traceSystem.incidences[
            trace.firstIncidence + offset];
        if (incidence.traceIndex != trace.traceIndex
            || incidence.halfFaceIndex >= controlCells.halfFaces.size()
            || incidence.controlCellIndex
                >= controlCells.controlCells.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow incidence binding is invalid");
        }
        const auto& halfFace = controlCells.halfFaces[
            incidence.halfFaceIndex];
        if (halfFace.kind != trace.kind
            || halfFace.sourceStableId != trace.sourceStableId
            || halfFace.controlVolumeIndex
                != controlCells.controlCells[incidence.controlCellIndex]
                    .controlVolumeIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow half-face binding is invalid");
        }
        if (halfFace.side
            == SceneFluidMimeticHalfFaceSide::MinusOrNegative) {
            if (result.minus != nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow has duplicate minus side");
            }
            result.minus = &halfFace;
            result.minusControlCellIndex = incidence.controlCellIndex;
        } else if (halfFace.side
                   == SceneFluidMimeticHalfFaceSide::PlusOrPositive) {
            if (result.plus != nullptr) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow has duplicate plus side");
            }
            result.plus = &halfFace;
            result.plusControlCellIndex = incidence.controlCellIndex;
        } else {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow side is invalid");
        }
    }
    if (result.minus == nullptr || result.plus == nullptr
        || result.minus->sourceIndex != result.plus->sourceIndex
        || result.minus->sourceStableId != result.plus->sourceStableId
        || result.minus->areaSquareMeters != result.plus->areaSquareMeters) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow pair is inconsistent");
    }
    return result;
}

void validateInputIdentity(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond) {
    validateSceneFluidMimeticTraceSystem(traceSystem, controlCells);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningFluxIntegrity(openingFlux);
    if (controlCells.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || controlCells.openingPatchFingerprint
            != openingFlux.openingPatchFingerprint
        || controlCells.surfaceDefinitionFingerprint
            != openingFlux.surfaceDefinitionFingerprint
        || controlCells.surfaceStateFingerprint
            != openingFlux.surfaceStateFingerprint
        || controlCells.structureDefinitionFingerprint
            != openingFlux.structureDefinitionFingerprint
        || controlCells.acceptedStepCount != faceLinks.acceptedStepCount
        || controlCells.acceptedStepCount != openingFlux.acceptedStepCount
        || controlCells.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || controlCells.simulationTimeSeconds
            != openingFlux.simulationTimeSeconds
        || controlCells.cellCounts != grid.cellCounts()
        || controlCells.lowerMeters != grid.lowerMeters()
        || controlCells.upperMeters != grid.upperMeters()
        || openingFlux.cellCounts != grid.cellCounts()
        || openingFlux.lowerMeters != grid.lowerMeters()
        || openingFlux.upperMeters != grid.upperMeters()
        || openingFlux.velocityFingerprint
            != sceneFluidOpeningFluxVelocityFingerprint(
                grid, predictedVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow input identity is invalid");
    }
}

void validateWallInputIdentity(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange) {
    validateSceneFluidMimeticTraceSystem(traceSystem, controlCells);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningFluxIntegrity(openingFlux);
    validateSceneFluidRegionWallExchangeIntegrity(wallExchange);
    if (!wallExchange.diagnostics.accepted
        || controlCells.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || controlCells.pressureControlVolumeFingerprint
            != wallExchange.currentPressureControlVolumeFingerprint
        || controlCells.openingPatchFingerprint
            != openingFlux.openingPatchFingerprint
        || controlCells.surfaceDefinitionFingerprint
            != openingFlux.surfaceDefinitionFingerprint
        || controlCells.surfaceDefinitionFingerprint
            != wallExchange.surfaceDefinitionFingerprint
        || controlCells.surfaceStateFingerprint
            != openingFlux.surfaceStateFingerprint
        || controlCells.surfaceStateFingerprint
            != wallExchange.surfaceStateFingerprint
        || controlCells.structureDefinitionFingerprint
            != openingFlux.structureDefinitionFingerprint
        || controlCells.structureDefinitionFingerprint
            != wallExchange.structureDefinitionFingerprint
        || controlCells.acceptedStepCount != faceLinks.acceptedStepCount
        || controlCells.acceptedStepCount != openingFlux.acceptedStepCount
        || controlCells.acceptedStepCount != wallExchange.acceptedStepCount
        || controlCells.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || controlCells.simulationTimeSeconds
            != openingFlux.simulationTimeSeconds
        || controlCells.simulationTimeSeconds
            != wallExchange.simulationTimeSeconds
        || controlCells.cellCounts != wallExchange.cellCounts
        || controlCells.lowerMeters != wallExchange.lowerMeters
        || controlCells.upperMeters != wallExchange.upperMeters
        || openingFlux.cellCounts != wallExchange.cellCounts
        || openingFlux.lowerMeters != wallExchange.lowerMeters
        || openingFlux.upperMeters != wallExchange.upperMeters
        || wallExchange.controlVolumes.size()
            != controlCells.controlCells.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow wall input identity is invalid");
    }
    for (const auto& cell : controlCells.controlCells) {
        if (cell.controlCellIndex >= wallExchange.controlVolumes.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow wall control index is invalid");
        }
        const auto& wall = wallExchange.controlVolumes[
            cell.controlCellIndex];
        if (wall.controlVolumeIndex != cell.controlVolumeIndex
            || wall.stableId != cell.stableId
            || wall.volumeCubicMeters != cell.volumeCubicMeters) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow wall control topology is invalid");
        }
    }
}

void validateRegionTransportInputIdentity(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidMimeticCorrectedTraceFlow& sourceCorrectedFlow,
    const SceneFluidRegionTransport& regionTransport,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline) {
    validateSceneFluidMimeticTraceSystem(traceSystem, controlCells);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningFluxIntegrity(openingFlux);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(
        sourceCorrectedFlow);
    validateSceneFluidRegionTransportIntegrity(regionTransport);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        currentBulkBaseline);
    if (!regionTransport.diagnostics.accepted
        || !sourceCorrectedFlow.accepted
        || !regionTransport.diagnostics.usesBulkVelocityIncrement
        || regionTransport.pressureProjectionFingerprint
            != sourceCorrectedFlow.fingerprint
        || sourceCorrectedFlow.mimeticControlCellFingerprint
            != controlCells.fingerprint
        || sourceCorrectedFlow.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || sourceCorrectedFlow.openingPatchFingerprint
            != openingFlux.openingPatchFingerprint
        || sourceCorrectedFlow.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || sourceCorrectedFlow.acceptedStepCount
            != controlCells.acceptedStepCount
        || sourceCorrectedFlow.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || sourceCorrectedFlow.densityKgPerCubicMeter
            != regionTransport.densityKgPerCubicMeter
        || regionTransport.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || regionTransport.currentBulkVelocityFingerprint
            != openingFlux.velocityFingerprint
        || regionTransport.currentBulkVelocityFingerprint
            != currentBulkBaseline.velocityFingerprint
        || regionTransport.acceptedStepCount
            != controlCells.acceptedStepCount
        || regionTransport.sourceSimulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || regionTransport.controlVolumes.size()
            != controlCells.controlCells.size()
        || currentBulkBaseline.mimeticControlCellFingerprint
            != controlCells.fingerprint
        || currentBulkBaseline.mimeticTraceSystemFingerprint
            != traceSystem.fingerprint
        || currentBulkBaseline.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || currentBulkBaseline.openingFluxFingerprint
            != openingFlux.fingerprint
        || currentBulkBaseline.regionWallExchangeFingerprint != 0
        || currentBulkBaseline.sourceDensityKgPerCubicMeter != 0.0
        || currentBulkBaseline.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || currentBulkBaseline.acceptedStepCount
            != controlCells.acceptedStepCount
        || currentBulkBaseline.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || currentBulkBaseline.componentCount
            != traceSystem.componentCount
        || currentBulkBaseline.traces.size()
            != traceSystem.sharedTraceCount) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace-flow input identity is invalid");
    }
    for (const auto& cell : controlCells.controlCells) {
        if (cell.controlVolumeIndex
            >= regionTransport.controlVolumes.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport control index is invalid");
        }
        const auto& transported = regionTransport.controlVolumes[
            cell.controlVolumeIndex];
        if (transported.controlVolumeIndex != cell.controlVolumeIndex
            || transported.stableId != cell.stableId
            || transported.volumeCubicMeters != cell.volumeCubicMeters) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport control topology is invalid");
        }
    }
}

double regionTransportPredictedFlow(
    const SceneFluidMimeticTrace& trace,
    const OrientedTrace& oriented,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidRegionTransport& regionTransport,
    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>&
        openingSamples) {
    const fluid::Vector3 velocity = averageVelocity(
        controlCells, regionTransport,
        oriented.minusControlCellIndex,
        oriented.plusControlCellIndex);
    if (trace.kind == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
        if (oriented.minus->sourceIndex >= faceLinks.links.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport face-link index is invalid");
        }
        const auto& link = faceLinks.links[oriented.minus->sourceIndex];
        if (link.linkIndex != oriented.minus->sourceIndex
            || link.stableId != trace.sourceStableId
            || link.geometryKind
                != SceneFluidPressureLinkGeometryKind::CartesianFace
            || link.minusControlVolumeIndex
                != controlCells.controlCells[
                    oriented.minusControlCellIndex].controlVolumeIndex
            || link.plusControlVolumeIndex
                != controlCells.controlCells[
                    oriented.plusControlCellIndex].controlVolumeIndex
            || link.areaSquareMeters != oriented.minus->areaSquareMeters
            || link.faceIndex >= faceLinks.faces.size()) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport face-link binding is invalid");
        }
        double result = link.areaSquareMeters
            * dot(velocity, link.unitNormalMinusToPlus);
        if (link.kind == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            const auto found = openingSamples.find(
                link.openingPatchStableId);
            if (found == openingSamples.end()) {
                throw std::invalid_argument(
                    "scene fluid mimetic regional-transport opening sample is missing");
            }
            result -= orientedOpeningSweep(link, *found->second);
            openingSamples.erase(found);
        } else if (link.kind
                   != SceneFluidPressureFaceLinkKind::SameRegion) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport face-link kind is invalid");
        }
        return result;
    }
    if (trace.kind
        != SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace kind is invalid");
    }
    const auto found = openingSamples.find(trace.sourceStableId);
    if (found == openingSamples.end()) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport embedded opening sample is missing");
    }
    const auto& sample = *found->second;
    const auto& minusCell = controlCells.controlCells[
        oriented.minusControlCellIndex];
    const auto& plusCell = controlCells.controlCells[
        oriented.plusControlCellIndex];
    if (sample.patchStableId != trace.sourceStableId
        || sample.areaSquareMeters != oriented.minus->areaSquareMeters
        || sample.negativeSideRegionId != minusCell.regionId
        || sample.positiveSideRegionId != plusCell.regionId) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport embedded opening identity is inconsistent");
    }
    const double result = sample.areaSquareMeters
            * dot(velocity, oriented.minus->outwardUnitNormal)
        - sample.surfaceSweepRateCubicMetersPerSecond;
    openingSamples.erase(found);
    return result;
}

} // namespace

SceneFluidMimeticTraceFlowPrediction sampleSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidMimeticTraceFlowLimits& limits) {
    validateInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux, grid,
        predictedVelocityMetersPerSecond);
    if (traceSystem.sharedTraceCount > limits.maximumSharedTraces
        || traceSystem.componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid mimetic trace-flow count limit exceeded");
    }
    const std::size_t expectedBytes = storageBytesForCounts(
        traceSystem.sharedTraceCount, traceSystem.componentCount);
    if (expectedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic trace-flow byte limit exceeded");
    }

    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        openingSamples;
    for (const auto& sample : openingFlux.samples) {
        if (sample.patchStableId == 0
            || !openingSamples.emplace(sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow has duplicate opening identity");
        }
    }

    SceneFluidMimeticTraceFlowPrediction result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.mimeticTraceSystemFingerprint = traceSystem.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.openingFluxFingerprint = openingFlux.fingerprint;
    result.velocityFingerprint = openingFlux.velocityFingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.ownedStorageBytes = expectedBytes;
    result.componentCount = traceSystem.componentCount;
    result.traces.reserve(traceSystem.sharedTraceCount);
    std::vector<CompensatedSum> componentBalances(result.componentCount);

    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        if (trace.kind != SceneFluidMimeticHalfFaceKind::CartesianTrace
            && trace.kind
                != SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow kind is invalid");
        }
        const OrientedTrace oriented = orientTrace(
            trace, traceSystem, controlCells);
        double predictedFlow = 0.0;
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            if (oriented.minus->sourceIndex >= faceLinks.links.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow face-link index is invalid");
            }
            const auto& link = faceLinks.links[oriented.minus->sourceIndex];
            if (link.linkIndex != oriented.minus->sourceIndex
                || link.stableId != trace.sourceStableId
                || link.geometryKind
                    != SceneFluidPressureLinkGeometryKind::CartesianFace
                || link.minusControlVolumeIndex
                    != controlCells.controlCells[
                        oriented.minusControlCellIndex].controlVolumeIndex
                || link.plusControlVolumeIndex
                    != controlCells.controlCells[
                        oriented.plusControlCellIndex].controlVolumeIndex
                || link.areaSquareMeters
                    != oriented.minus->areaSquareMeters
                || link.faceIndex >= faceLinks.faces.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow face-link binding is invalid");
            }
            if (link.kind == SceneFluidPressureFaceLinkKind::SameRegion) {
                predictedFlow = link.areaSquareMeters * faceVelocity(
                    faceLinks.faces[link.faceIndex], grid,
                    predictedVelocityMetersPerSecond);
            } else if (link.kind
                       == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
                const auto found = openingSamples.find(
                    link.openingPatchStableId);
                if (found == openingSamples.end()) {
                    throw std::invalid_argument(
                        "scene fluid mimetic trace-flow opening sample is missing");
                }
                predictedFlow = orientedOpeningFlow(link, *found->second);
                openingSamples.erase(found);
            } else {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow face-link kind is invalid");
            }
            ++result.cartesianTraceCount;
        } else {
            const auto found = openingSamples.find(trace.sourceStableId);
            if (found == openingSamples.end()) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow embedded opening sample is missing");
            }
            const auto& sample = *found->second;
            const auto& minusCell = controlCells.controlCells[
                oriented.minusControlCellIndex];
            const auto& plusCell = controlCells.controlCells[
                oriented.plusControlCellIndex];
            if (sample.patchStableId != trace.sourceStableId
                || sample.areaSquareMeters
                    != oriented.minus->areaSquareMeters
                || sample.negativeSideRegionId != minusCell.regionId
                || sample.positiveSideRegionId != plusCell.regionId) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace-flow embedded opening identity is inconsistent");
            }
            predictedFlow =
                sample.relativeVolumeFlowRateCubicMetersPerSecond;
            openingSamples.erase(found);
            ++result.authoredOpeningTraceCount;
        }
        if (!std::isfinite(predictedFlow)) {
            throw std::overflow_error(
                "scene fluid mimetic trace-flow prediction is non-finite");
        }
        SceneFluidMimeticPredictedTraceFlow predicted;
        predicted.sharedTraceOrdinal = result.traces.size();
        predicted.traceIndex = trace.traceIndex;
        predicted.stableId = trace.stableId;
        predicted.kind = trace.kind;
        predicted.sourceStableId = trace.sourceStableId;
        predicted.componentIndex = trace.componentIndex;
        predicted.minusControlCellIndex = oriented.minusControlCellIndex;
        predicted.plusControlCellIndex = oriented.plusControlCellIndex;
        predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
            predictedFlow;
        result.traces.push_back(predicted);
        componentBalances[trace.componentIndex].add(predictedFlow);
        componentBalances[trace.componentIndex].add(-predictedFlow);
        result
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(predictedFlow));
    }
    if (result.traces.size() != traceSystem.sharedTraceCount
        || !openingSamples.empty()) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow did not consume its complete topology");
    }
    result.componentBalanceResidualsCubicMetersPerSecond.resize(
        result.componentCount);
    for (std::size_t component = 0;
         component < result.componentCount; ++component) {
        result.componentBalanceResidualsCubicMetersPerSecond[component] =
            componentBalances[component].value();
        result.maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond,
                std::abs(result
                    .componentBalanceResidualsCubicMetersPerSecond[
                        component]));
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticTraceFlowPrediction(
        result, controlCells, traceSystem, faceLinks, openingFlux, grid,
        predictedVelocityMetersPerSecond);
    return result;
}

SceneFluidMimeticTraceFlowPrediction sampleSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange,
    const SceneFluidMimeticTraceFlowLimits& limits) {
    validateWallInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux, wallExchange);
    if (traceSystem.sharedTraceCount > limits.maximumSharedTraces
        || traceSystem.componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid mimetic wall trace-flow count limit exceeded");
    }
    const std::size_t expectedBytes = storageBytesForCounts(
        traceSystem.sharedTraceCount, traceSystem.componentCount);
    if (expectedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic wall trace-flow byte limit exceeded");
    }

    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        openingSamples;
    for (const auto& sample : openingFlux.samples) {
        if (sample.patchStableId == 0
            || !openingSamples.emplace(sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid mimetic wall trace-flow has duplicate opening identity");
        }
    }

    SceneFluidMimeticTraceFlowPrediction result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.mimeticTraceSystemFingerprint = traceSystem.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.openingFluxFingerprint = openingFlux.fingerprint;
    result.velocityFingerprint = openingFlux.velocityFingerprint;
    result.regionWallExchangeFingerprint = wallExchange.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.sourceDensityKgPerCubicMeter =
        wallExchange.densityKgPerCubicMeter;
    result.ownedStorageBytes = expectedBytes;
    result.componentCount = traceSystem.componentCount;
    result.traces.reserve(traceSystem.sharedTraceCount);
    std::vector<CompensatedSum> componentBalances(result.componentCount);

    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        if (trace.kind != SceneFluidMimeticHalfFaceKind::CartesianTrace
            && trace.kind
                != SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            throw std::invalid_argument(
                "scene fluid mimetic wall trace-flow kind is invalid");
        }
        const OrientedTrace oriented = orientTrace(
            trace, traceSystem, controlCells);
        const fluid::Vector3 velocity = averageVelocity(
            wallExchange, oriented.minusControlCellIndex,
            oriented.plusControlCellIndex);
        double predictedFlow = 0.0;
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            if (oriented.minus->sourceIndex >= faceLinks.links.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic wall trace-flow face-link index is invalid");
            }
            const auto& link = faceLinks.links[oriented.minus->sourceIndex];
            if (link.linkIndex != oriented.minus->sourceIndex
                || link.stableId != trace.sourceStableId
                || link.geometryKind
                    != SceneFluidPressureLinkGeometryKind::CartesianFace
                || link.minusControlVolumeIndex
                    != controlCells.controlCells[
                        oriented.minusControlCellIndex].controlVolumeIndex
                || link.plusControlVolumeIndex
                    != controlCells.controlCells[
                        oriented.plusControlCellIndex].controlVolumeIndex
                || link.areaSquareMeters
                    != oriented.minus->areaSquareMeters
                || link.faceIndex >= faceLinks.faces.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic wall trace-flow face-link binding is invalid");
            }
            predictedFlow = link.areaSquareMeters
                * dot(velocity, link.unitNormalMinusToPlus);
            if (link.kind
                == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
                const auto found = openingSamples.find(
                    link.openingPatchStableId);
                if (found == openingSamples.end()) {
                    throw std::invalid_argument(
                        "scene fluid mimetic wall trace-flow opening sample is missing");
                }
                predictedFlow -= orientedOpeningSweep(link, *found->second);
                openingSamples.erase(found);
            } else if (link.kind
                       != SceneFluidPressureFaceLinkKind::SameRegion) {
                throw std::invalid_argument(
                    "scene fluid mimetic wall trace-flow face-link kind is invalid");
            }
            ++result.cartesianTraceCount;
        } else {
            const auto found = openingSamples.find(trace.sourceStableId);
            if (found == openingSamples.end()) {
                throw std::invalid_argument(
                    "scene fluid mimetic wall trace-flow embedded opening sample is missing");
            }
            const auto& sample = *found->second;
            const auto& minusCell = controlCells.controlCells[
                oriented.minusControlCellIndex];
            const auto& plusCell = controlCells.controlCells[
                oriented.plusControlCellIndex];
            if (sample.patchStableId != trace.sourceStableId
                || sample.areaSquareMeters
                    != oriented.minus->areaSquareMeters
                || sample.negativeSideRegionId != minusCell.regionId
                || sample.positiveSideRegionId != plusCell.regionId) {
                throw std::invalid_argument(
                    "scene fluid mimetic wall trace-flow embedded opening identity is inconsistent");
            }
            predictedFlow = sample.areaSquareMeters
                    * dot(velocity, oriented.minus->outwardUnitNormal)
                - sample.surfaceSweepRateCubicMetersPerSecond;
            openingSamples.erase(found);
            ++result.authoredOpeningTraceCount;
        }
        if (!std::isfinite(predictedFlow)) {
            throw std::overflow_error(
                "scene fluid mimetic wall trace-flow prediction is non-finite");
        }
        SceneFluidMimeticPredictedTraceFlow predicted;
        predicted.sharedTraceOrdinal = result.traces.size();
        predicted.traceIndex = trace.traceIndex;
        predicted.stableId = trace.stableId;
        predicted.kind = trace.kind;
        predicted.sourceStableId = trace.sourceStableId;
        predicted.componentIndex = trace.componentIndex;
        predicted.minusControlCellIndex = oriented.minusControlCellIndex;
        predicted.plusControlCellIndex = oriented.plusControlCellIndex;
        predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
            predictedFlow;
        result.traces.push_back(predicted);
        componentBalances[trace.componentIndex].add(predictedFlow);
        componentBalances[trace.componentIndex].add(-predictedFlow);
        result
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(predictedFlow));
    }
    if (result.traces.size() != traceSystem.sharedTraceCount
        || !openingSamples.empty()) {
        throw std::invalid_argument(
            "scene fluid mimetic wall trace-flow did not consume its complete topology");
    }
    result.componentBalanceResidualsCubicMetersPerSecond.resize(
        result.componentCount);
    for (std::size_t component = 0;
         component < result.componentCount; ++component) {
        result.componentBalanceResidualsCubicMetersPerSecond[component] =
            componentBalances[component].value();
        result.maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond,
                std::abs(result
                    .componentBalanceResidualsCubicMetersPerSecond[
                        component]));
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticTraceFlowPrediction(
        result, controlCells, traceSystem, faceLinks, openingFlux,
        wallExchange);
    return result;
}

SceneFluidMimeticTraceFlowContinuation
continueSceneFluidMimeticTraceFlowsFixedTopology(
    const SceneFluidMimeticCorrectedTraceFlow& previousCorrectedFlow,
    const SceneFluidMimeticTraceFlowPrediction& previousBulkBaseline,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline,
    const SceneFluidMimeticTraceFlowLimits& limits) {
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(
        previousCorrectedFlow);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        previousBulkBaseline);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        currentBulkBaseline);
    const std::size_t traceCount = previousCorrectedFlow.traces.size();
    if (!previousCorrectedFlow.accepted
        || previousBulkBaseline.regionWallExchangeFingerprint != 0
        || currentBulkBaseline.regionWallExchangeFingerprint != 0
        || previousBulkBaseline.sourceDensityKgPerCubicMeter != 0.0
        || currentBulkBaseline.sourceDensityKgPerCubicMeter != 0.0
        || previousCorrectedFlow.mimeticControlCellFingerprint
            != previousBulkBaseline.mimeticControlCellFingerprint
        || previousCorrectedFlow.mimeticControlCellFingerprint
            != currentBulkBaseline.mimeticControlCellFingerprint
        || previousCorrectedFlow.fullTraceSystemFingerprint
            != previousBulkBaseline.mimeticTraceSystemFingerprint
        || previousCorrectedFlow.fullTraceSystemFingerprint
            != currentBulkBaseline.mimeticTraceSystemFingerprint
        || previousCorrectedFlow.pressureFaceLinkFingerprint
            != previousBulkBaseline.pressureFaceLinkFingerprint
        || previousCorrectedFlow.pressureFaceLinkFingerprint
            != currentBulkBaseline.pressureFaceLinkFingerprint
        || previousCorrectedFlow.structureDefinitionFingerprint
            != previousBulkBaseline.structureDefinitionFingerprint
        || previousCorrectedFlow.structureDefinitionFingerprint
            != currentBulkBaseline.structureDefinitionFingerprint
        || previousCorrectedFlow.acceptedStepCount
            != previousBulkBaseline.acceptedStepCount
        || previousCorrectedFlow.acceptedStepCount
            != currentBulkBaseline.acceptedStepCount
        || previousCorrectedFlow.simulationTimeSeconds
            != previousBulkBaseline.simulationTimeSeconds
        || previousCorrectedFlow.simulationTimeSeconds
            != currentBulkBaseline.simulationTimeSeconds
        || traceCount != previousBulkBaseline.traces.size()
        || traceCount != currentBulkBaseline.traces.size()
        || previousBulkBaseline.componentCount
            != currentBulkBaseline.componentCount) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow continuation is foreign");
    }
    if (traceCount > limits.maximumSharedTraces
        || currentBulkBaseline.componentCount > limits.maximumComponents
        || currentBulkBaseline.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic trace-flow continuation exceeds its limits");
    }

    SceneFluidMimeticTraceFlowContinuation result;
    result.previousCorrectedTraceFlowFingerprint =
        previousCorrectedFlow.fingerprint;
    result.previousBaselinePredictionFingerprint =
        previousBulkBaseline.fingerprint;
    result.currentBaselinePredictionFingerprint =
        currentBulkBaseline.fingerprint;
    result.mimeticControlCellFingerprint =
        currentBulkBaseline.mimeticControlCellFingerprint;
    result.mimeticTraceSystemFingerprint =
        currentBulkBaseline.mimeticTraceSystemFingerprint;
    result.pressureFaceLinkFingerprint =
        currentBulkBaseline.pressureFaceLinkFingerprint;
    result.structureDefinitionFingerprint =
        currentBulkBaseline.structureDefinitionFingerprint;
    result.acceptedStepCount = currentBulkBaseline.acceptedStepCount;
    result.simulationTimeSeconds = currentBulkBaseline.simulationTimeSeconds;
    result.prediction = currentBulkBaseline;
    result.prediction.maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    std::vector<CompensatedSum> componentBalances(
        result.prediction.componentCount);
    for (std::size_t index = 0; index < traceCount; ++index) {
        const auto& corrected = previousCorrectedFlow.traces[index];
        const auto& previous = previousBulkBaseline.traces[index];
        const auto& current = currentBulkBaseline.traces[index];
        auto& continued = result.prediction.traces[index];
        if (corrected.sharedTraceOrdinal != index
            || previous.sharedTraceOrdinal != index
            || current.sharedTraceOrdinal != index
            || corrected.traceIndex != previous.traceIndex
            || corrected.traceIndex != current.traceIndex
            || corrected.stableId != previous.stableId
            || corrected.stableId != current.stableId
            || corrected.kind != previous.kind
            || corrected.kind != current.kind
            || corrected.sourceStableId != previous.sourceStableId
            || corrected.sourceStableId != current.sourceStableId
            || corrected.componentIndex != previous.componentIndex
            || corrected.componentIndex != current.componentIndex
            || corrected.minusControlCellIndex
                != previous.minusControlCellIndex
            || corrected.minusControlCellIndex
                != current.minusControlCellIndex
            || corrected.plusControlCellIndex
                != previous.plusControlCellIndex
            || corrected.plusControlCellIndex
                != current.plusControlCellIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow continuation topology changed");
        }
        const double carriedCorrection =
            corrected.correctedRelativeVolumeFlowRateCubicMetersPerSecond
            - previous.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        const double bulkIncrement =
            current.predictedRelativeVolumeFlowRateCubicMetersPerSecond
            - previous.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        const double predicted =
            corrected.correctedRelativeVolumeFlowRateCubicMetersPerSecond
            + bulkIncrement;
        if (!std::isfinite(carriedCorrection)
            || !std::isfinite(bulkIncrement) || !std::isfinite(predicted)) {
            throw std::overflow_error(
                "scene fluid mimetic trace-flow continuation is non-finite");
        }
        continued.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
            predicted;
        result.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond,
                std::abs(carriedCorrection));
        result.maximumAbsoluteBulkIncrementCubicMetersPerSecond = std::max(
            result.maximumAbsoluteBulkIncrementCubicMetersPerSecond,
            std::abs(bulkIncrement));
        result.prediction
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                result.prediction
                    .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(predicted));
        componentBalances[continued.componentIndex].add(predicted);
        componentBalances[continued.componentIndex].add(-predicted);
    }
    result.prediction.maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond =
        0.0;
    for (std::size_t component = 0;
         component < result.prediction.componentCount; ++component) {
        const double balance = componentBalances[component].value();
        result.prediction
            .componentBalanceResidualsCubicMetersPerSecond[component] =
            balance;
        result.prediction
            .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond =
            std::max(
                result.prediction
                    .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond,
                std::abs(balance));
    }
    result.prediction.fingerprint = productFingerprint(result.prediction);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(result.prediction);
    result.ownedStorageBytes = result.prediction.ownedStorageBytes;
    result.finite = std::isfinite(
            result.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond)
        && std::isfinite(
            result.maximumAbsoluteBulkIncrementCubicMetersPerSecond);
    result.fingerprint = continuationFingerprint(result);
    validateSceneFluidMimeticTraceFlowContinuation(
        result, previousCorrectedFlow, previousBulkBaseline,
        currentBulkBaseline);
    return result;
}

SceneFluidMimeticRegionTransportFlowPrediction
predictSceneFluidMimeticTraceFlowsFromRegionTransportFixedTopology(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidMimeticCorrectedTraceFlow& sourceCorrectedFlow,
    const SceneFluidRegionTransport& regionTransport,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline,
    const SceneFluidMimeticTraceFlowLimits& limits) {
    validateRegionTransportInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux,
        sourceCorrectedFlow, regionTransport, currentBulkBaseline);
    if (traceSystem.sharedTraceCount > limits.maximumSharedTraces
        || traceSystem.componentCount > limits.maximumComponents
        || currentBulkBaseline.ownedStorageBytes
            > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic regional-transport trace-flow limit exceeded");
    }

    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        openingSamples;
    for (const auto& sample : openingFlux.samples) {
        if (sample.patchStableId == 0
            || !openingSamples.emplace(
                sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport trace-flow has duplicate opening identity");
        }
    }

    SceneFluidMimeticRegionTransportFlowPrediction result;
    result.regionTransportFingerprint = regionTransport.fingerprint;
    result.sourceCorrectedTraceFlowFingerprint =
        sourceCorrectedFlow.fingerprint;
    result.currentBulkBaselinePredictionFingerprint =
        currentBulkBaseline.fingerprint;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.mimeticTraceSystemFingerprint = traceSystem.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.transportTargetSimulationTimeSeconds =
        regionTransport.targetSimulationTimeSeconds;
    result.prediction = currentBulkBaseline;
    result.prediction
        .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    std::vector<CompensatedSum> componentBalances(
        result.prediction.componentCount);
    std::size_t sharedOrdinal = 0;
    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        const auto oriented = orientTrace(
            trace, traceSystem, controlCells);
        const double predictedFlow = regionTransportPredictedFlow(
            trace, oriented, controlCells, faceLinks,
            regionTransport, openingSamples);
        if (!std::isfinite(predictedFlow)) {
            throw std::overflow_error(
                "scene fluid mimetic regional-transport trace-flow prediction is non-finite");
        }
        auto& predicted = result.prediction.traces[sharedOrdinal];
        const double baselineFlow =
            predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        predicted.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
            predictedFlow;
        result
            .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond,
                std::abs(predictedFlow - baselineFlow));
        result.prediction
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                result.prediction
                    .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(predictedFlow));
        componentBalances[trace.componentIndex].add(predictedFlow);
        componentBalances[trace.componentIndex].add(-predictedFlow);
        ++sharedOrdinal;
    }
    if (sharedOrdinal != result.prediction.traces.size()
        || !openingSamples.empty()) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace-flow did not consume its complete topology");
    }
    result.prediction
        .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond = 0.0;
    for (std::size_t component = 0;
         component < result.prediction.componentCount; ++component) {
        const double balance = componentBalances[component].value();
        result.prediction
            .componentBalanceResidualsCubicMetersPerSecond[component] =
            balance;
        result.prediction
            .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond =
            std::max(
                result.prediction
                    .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond,
                std::abs(balance));
    }
    result.prediction.fingerprint = productFingerprint(result.prediction);
    result.ownedStorageBytes = result.prediction.ownedStorageBytes;
    result.finite = std::isfinite(
        result
            .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond);
    result.fingerprint = regionTransportPredictionFingerprint(result);
    validateSceneFluidMimeticRegionTransportFlowPrediction(
        result, controlCells, traceSystem, faceLinks, openingFlux,
        sourceCorrectedFlow, regionTransport, currentBulkBaseline);
    return result;
}

void validateSceneFluidMimeticTraceFlowPredictionIntegrity(
    const SceneFluidMimeticTraceFlowPrediction& prediction) {
    const bool fixedEpochSource = prediction.regionWallExchangeFingerprint == 0
        && prediction.sourceDensityKgPerCubicMeter == 0.0;
    const bool wallAdjustedSource =
        prediction.regionWallExchangeFingerprint != 0
        && std::isfinite(prediction.sourceDensityKgPerCubicMeter)
        && prediction.sourceDensityKgPerCubicMeter > 0.0;
    if (prediction.version != sceneFluidMimeticTraceFlowVersion
        || prediction.fingerprint == 0
        || prediction.mimeticControlCellFingerprint == 0
        || prediction.mimeticTraceSystemFingerprint == 0
        || prediction.pressureFaceLinkFingerprint == 0
        || prediction.openingFluxFingerprint == 0
        || prediction.velocityFingerprint == 0
        || prediction.structureDefinitionFingerprint == 0
        || !std::isfinite(prediction.simulationTimeSeconds)
        || (!fixedEpochSource && !wallAdjustedSource)
        || prediction.componentCount == 0
        || prediction.traces.empty()
        || prediction.componentBalanceResidualsCubicMetersPerSecond.size()
            != prediction.componentCount
        || prediction.cartesianTraceCount
                + prediction.authoredOpeningTraceCount
            != prediction.traces.size()
        || prediction.ownedStorageBytes != storageBytesForCounts(
            prediction.traces.size(), prediction.componentCount)) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow integrity is invalid");
    }
    std::vector<CompensatedSum> componentBalances(
        prediction.componentCount);
    std::size_t cartesianCount = 0;
    std::size_t openingCount = 0;
    double maximumFlow = 0.0;
    std::size_t previousTraceIndex = 0;
    bool hasPreviousTrace = false;
    for (const auto& trace : prediction.traces) {
        if (trace.sharedTraceOrdinal
                != &trace - prediction.traces.data()
            || trace.stableId == 0 || trace.sourceStableId == 0
            || trace.componentIndex >= prediction.componentCount
            || (hasPreviousTrace && trace.traceIndex <= previousTraceIndex)
            || !std::isfinite(
                trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow trace is invalid");
        }
        hasPreviousTrace = true;
        previousTraceIndex = trace.traceIndex;
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            ++cartesianCount;
        } else if (trace.kind
                   == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            ++openingCount;
        } else {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow contains a wall trace");
        }
        const double flow =
            trace.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        componentBalances[trace.componentIndex].add(flow);
        componentBalances[trace.componentIndex].add(-flow);
        maximumFlow = std::max(maximumFlow, std::abs(flow));
    }
    double maximumComponentBalance = 0.0;
    for (std::size_t component = 0;
         component < prediction.componentCount; ++component) {
        const double balance = componentBalances[component].value();
        if (prediction
                .componentBalanceResidualsCubicMetersPerSecond[component]
            != balance) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow component balance is invalid");
        }
        maximumComponentBalance = std::max(
            maximumComponentBalance, std::abs(balance));
    }
    if (prediction.cartesianTraceCount != cartesianCount
        || prediction.authoredOpeningTraceCount != openingCount
        || prediction
                .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond
            != maximumFlow
        || prediction
                .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond
            != maximumComponentBalance
        || productFingerprint(prediction) != prediction.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow summary is invalid");
    }
}

void validateSceneFluidMimeticTraceFlowContinuationIntegrity(
    const SceneFluidMimeticTraceFlowContinuation& continuation) {
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        continuation.prediction);
    if (continuation.version
            != sceneFluidMimeticTraceFlowContinuationVersion
        || continuation.fingerprint == 0
        || continuation.previousCorrectedTraceFlowFingerprint == 0
        || continuation.previousBaselinePredictionFingerprint == 0
        || continuation.currentBaselinePredictionFingerprint == 0
        || continuation.mimeticControlCellFingerprint == 0
        || continuation.mimeticTraceSystemFingerprint == 0
        || continuation.pressureFaceLinkFingerprint == 0
        || continuation.structureDefinitionFingerprint == 0
        || !std::isfinite(continuation.simulationTimeSeconds)
        || continuation.ownedStorageBytes
            != continuation.prediction.ownedStorageBytes
        || !std::isfinite(
            continuation.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond)
        || continuation.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond
            < 0.0
        || !std::isfinite(
            continuation.maximumAbsoluteBulkIncrementCubicMetersPerSecond)
        || continuation.maximumAbsoluteBulkIncrementCubicMetersPerSecond
            < 0.0
        || !continuation.finite
        || continuation.mimeticControlCellFingerprint
            != continuation.prediction.mimeticControlCellFingerprint
        || continuation.mimeticTraceSystemFingerprint
            != continuation.prediction.mimeticTraceSystemFingerprint
        || continuation.pressureFaceLinkFingerprint
            != continuation.prediction.pressureFaceLinkFingerprint
        || continuation.structureDefinitionFingerprint
            != continuation.prediction.structureDefinitionFingerprint
        || continuation.acceptedStepCount
            != continuation.prediction.acceptedStepCount
        || continuation.simulationTimeSeconds
            != continuation.prediction.simulationTimeSeconds
        || continuation.fingerprint
            != continuationFingerprint(continuation)) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow continuation integrity is invalid");
    }
}

void validateSceneFluidMimeticTraceFlowContinuation(
    const SceneFluidMimeticTraceFlowContinuation& continuation,
    const SceneFluidMimeticCorrectedTraceFlow& previousCorrectedFlow,
    const SceneFluidMimeticTraceFlowPrediction& previousBulkBaseline,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline) {
    validateSceneFluidMimeticTraceFlowContinuationIntegrity(continuation);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(
        previousCorrectedFlow);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        previousBulkBaseline);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        currentBulkBaseline);
    if (continuation.previousCorrectedTraceFlowFingerprint
            != previousCorrectedFlow.fingerprint
        || continuation.previousBaselinePredictionFingerprint
            != previousBulkBaseline.fingerprint
        || continuation.currentBaselinePredictionFingerprint
            != currentBulkBaseline.fingerprint
        || continuation.prediction.traces.size()
            != previousCorrectedFlow.traces.size()
        || continuation.prediction.traces.size()
            != previousBulkBaseline.traces.size()
        || continuation.prediction.traces.size()
            != currentBulkBaseline.traces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow continuation binding is invalid");
    }
    double maximumCorrection = 0.0;
    double maximumIncrement = 0.0;
    for (std::size_t index = 0;
         index < continuation.prediction.traces.size(); ++index) {
        const auto& corrected = previousCorrectedFlow.traces[index];
        const auto& previous = previousBulkBaseline.traces[index];
        const auto& current = currentBulkBaseline.traces[index];
        const auto& continued = continuation.prediction.traces[index];
        const double carriedCorrection =
            corrected.correctedRelativeVolumeFlowRateCubicMetersPerSecond
            - previous.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        const double bulkIncrement =
            current.predictedRelativeVolumeFlowRateCubicMetersPerSecond
            - previous.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        if (continued.sharedTraceOrdinal != index
            || continued.traceIndex != corrected.traceIndex
            || continued.stableId != corrected.stableId
            || continued.kind != corrected.kind
            || continued.sourceStableId != corrected.sourceStableId
            || continued.componentIndex != corrected.componentIndex
            || continued.minusControlCellIndex
                != corrected.minusControlCellIndex
            || continued.plusControlCellIndex
                != corrected.plusControlCellIndex
            || continued.predictedRelativeVolumeFlowRateCubicMetersPerSecond
                != corrected
                        .correctedRelativeVolumeFlowRateCubicMetersPerSecond
                    + bulkIncrement) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow continuation value changed");
        }
        maximumCorrection = std::max(
            maximumCorrection, std::abs(carriedCorrection));
        maximumIncrement = std::max(
            maximumIncrement, std::abs(bulkIncrement));
    }
    if (continuation.maximumAbsoluteCarriedCorrectionCubicMetersPerSecond
            != maximumCorrection
        || continuation.maximumAbsoluteBulkIncrementCubicMetersPerSecond
            != maximumIncrement) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow continuation summary changed");
    }
}

void validateSceneFluidMimeticRegionTransportFlowPredictionIntegrity(
    const SceneFluidMimeticRegionTransportFlowPrediction& prediction) {
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(
        prediction.prediction);
    if (prediction.version
            != sceneFluidMimeticRegionTransportFlowPredictionVersion
        || prediction.fingerprint == 0
        || prediction.regionTransportFingerprint == 0
        || prediction.sourceCorrectedTraceFlowFingerprint == 0
        || prediction.currentBulkBaselinePredictionFingerprint == 0
        || prediction.mimeticControlCellFingerprint == 0
        || prediction.mimeticTraceSystemFingerprint == 0
        || prediction.pressureFaceLinkFingerprint == 0
        || prediction.structureDefinitionFingerprint == 0
        || !std::isfinite(prediction.simulationTimeSeconds)
        || !std::isfinite(
            prediction.transportTargetSimulationTimeSeconds)
        || !(prediction.transportTargetSimulationTimeSeconds
            > prediction.simulationTimeSeconds)
        || prediction.ownedStorageBytes
            != prediction.prediction.ownedStorageBytes
        || !std::isfinite(
            prediction
                .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond)
        || prediction
               .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond
            < 0.0
        || !prediction.finite
        || prediction.mimeticControlCellFingerprint
            != prediction.prediction.mimeticControlCellFingerprint
        || prediction.mimeticTraceSystemFingerprint
            != prediction.prediction.mimeticTraceSystemFingerprint
        || prediction.pressureFaceLinkFingerprint
            != prediction.prediction.pressureFaceLinkFingerprint
        || prediction.structureDefinitionFingerprint
            != prediction.prediction.structureDefinitionFingerprint
        || prediction.acceptedStepCount
            != prediction.prediction.acceptedStepCount
        || prediction.simulationTimeSeconds
            != prediction.prediction.simulationTimeSeconds
        || prediction.prediction.regionWallExchangeFingerprint != 0
        || prediction.prediction.sourceDensityKgPerCubicMeter != 0.0
        || prediction.fingerprint
            != regionTransportPredictionFingerprint(prediction)) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace-flow integrity is invalid");
    }
}

void validateSceneFluidMimeticRegionTransportFlowPrediction(
    const SceneFluidMimeticRegionTransportFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidMimeticCorrectedTraceFlow& sourceCorrectedFlow,
    const SceneFluidRegionTransport& regionTransport,
    const SceneFluidMimeticTraceFlowPrediction& currentBulkBaseline) {
    validateRegionTransportInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux,
        sourceCorrectedFlow, regionTransport, currentBulkBaseline);
    validateSceneFluidMimeticRegionTransportFlowPredictionIntegrity(
        prediction);
    if (prediction.regionTransportFingerprint
            != regionTransport.fingerprint
        || prediction.sourceCorrectedTraceFlowFingerprint
            != sourceCorrectedFlow.fingerprint
        || prediction.currentBulkBaselinePredictionFingerprint
            != currentBulkBaseline.fingerprint
        || prediction.mimeticControlCellFingerprint
            != controlCells.fingerprint
        || prediction.mimeticTraceSystemFingerprint
            != traceSystem.fingerprint
        || prediction.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || prediction.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || prediction.acceptedStepCount
            != controlCells.acceptedStepCount
        || prediction.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || prediction.transportTargetSimulationTimeSeconds
            != regionTransport.targetSimulationTimeSeconds
        || prediction.prediction.traces.size()
            != currentBulkBaseline.traces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace-flow binding is invalid");
    }
    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        openingSamples;
    for (const auto& sample : openingFlux.samples) {
        if (sample.patchStableId == 0
            || !openingSamples.emplace(
                sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport trace-flow has duplicate opening identity");
        }
    }
    double maximumDifference = 0.0;
    std::size_t sharedOrdinal = 0;
    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        const auto oriented = orientTrace(
            trace, traceSystem, controlCells);
        const double expectedFlow = regionTransportPredictedFlow(
            trace, oriented, controlCells, faceLinks,
            regionTransport, openingSamples);
        const auto& actual = prediction.prediction.traces[sharedOrdinal];
        const auto& baseline = currentBulkBaseline.traces[sharedOrdinal];
        if (actual.sharedTraceOrdinal != baseline.sharedTraceOrdinal
            || actual.traceIndex != baseline.traceIndex
            || actual.stableId != baseline.stableId
            || actual.kind != baseline.kind
            || actual.sourceStableId != baseline.sourceStableId
            || actual.componentIndex != baseline.componentIndex
            || actual.minusControlCellIndex
                != baseline.minusControlCellIndex
            || actual.plusControlCellIndex
                != baseline.plusControlCellIndex
            || actual.predictedRelativeVolumeFlowRateCubicMetersPerSecond
                != expectedFlow) {
            throw std::invalid_argument(
                "scene fluid mimetic regional-transport trace-flow value changed");
        }
        maximumDifference = std::max(
            maximumDifference,
            std::abs(
                expectedFlow
                - baseline
                    .predictedRelativeVolumeFlowRateCubicMetersPerSecond));
        ++sharedOrdinal;
    }
    if (sharedOrdinal != prediction.prediction.traces.size()
        || !openingSamples.empty()
        || prediction
               .maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond
            != maximumDifference) {
        throw std::invalid_argument(
            "scene fluid mimetic regional-transport trace-flow summary changed");
    }
}

void validateSceneFluidMimeticTraceFlowPrediction(
    const SceneFluidMimeticTraceFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond) {
    validateInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux, grid,
        predictedVelocityMetersPerSecond);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(prediction);
    if (prediction.mimeticControlCellFingerprint != controlCells.fingerprint
        || prediction.mimeticTraceSystemFingerprint != traceSystem.fingerprint
        || prediction.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || prediction.openingFluxFingerprint != openingFlux.fingerprint
        || prediction.velocityFingerprint != openingFlux.velocityFingerprint
        || prediction.regionWallExchangeFingerprint != 0
        || prediction.sourceDensityKgPerCubicMeter != 0.0
        || prediction.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || prediction.acceptedStepCount != controlCells.acceptedStepCount
        || prediction.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || prediction.componentCount != traceSystem.componentCount
        || prediction.traces.size() != traceSystem.sharedTraceCount) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-flow prediction is foreign");
    }
    std::size_t sharedOrdinal = 0;
    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        const auto& predicted = prediction.traces[sharedOrdinal];
        const auto oriented = orientTrace(trace, traceSystem, controlCells);
        if (predicted.sharedTraceOrdinal != sharedOrdinal
            || predicted.traceIndex != trace.traceIndex
            || predicted.stableId != trace.stableId
            || predicted.kind != trace.kind
            || predicted.sourceStableId != trace.sourceStableId
            || predicted.componentIndex != trace.componentIndex
            || predicted.minusControlCellIndex
                != oriented.minusControlCellIndex
            || predicted.plusControlCellIndex
                != oriented.plusControlCellIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-flow topology is invalid");
        }
        ++sharedOrdinal;
    }
}

void validateSceneFluidMimeticTraceFlowPrediction(
    const SceneFluidMimeticTraceFlowPrediction& prediction,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& traceSystem,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningFluxSet& openingFlux,
    const SceneFluidRegionWallExchange& wallExchange) {
    validateWallInputIdentity(
        controlCells, traceSystem, faceLinks, openingFlux, wallExchange);
    validateSceneFluidMimeticTraceFlowPredictionIntegrity(prediction);
    if (prediction.mimeticControlCellFingerprint != controlCells.fingerprint
        || prediction.mimeticTraceSystemFingerprint != traceSystem.fingerprint
        || prediction.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || prediction.openingFluxFingerprint != openingFlux.fingerprint
        || prediction.velocityFingerprint != openingFlux.velocityFingerprint
        || prediction.regionWallExchangeFingerprint
            != wallExchange.fingerprint
        || prediction.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || prediction.acceptedStepCount != controlCells.acceptedStepCount
        || prediction.simulationTimeSeconds
            != controlCells.simulationTimeSeconds
        || prediction.sourceDensityKgPerCubicMeter
            != wallExchange.densityKgPerCubicMeter
        || prediction.componentCount != traceSystem.componentCount
        || prediction.traces.size() != traceSystem.sharedTraceCount) {
        throw std::invalid_argument(
            "scene fluid mimetic wall trace-flow prediction is foreign");
    }
    std::size_t sharedOrdinal = 0;
    for (const auto& trace : traceSystem.traces) {
        if (trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        const auto& predicted = prediction.traces[sharedOrdinal];
        const auto oriented = orientTrace(trace, traceSystem, controlCells);
        if (predicted.sharedTraceOrdinal != sharedOrdinal
            || predicted.traceIndex != trace.traceIndex
            || predicted.stableId != trace.stableId
            || predicted.kind != trace.kind
            || predicted.sourceStableId != trace.sourceStableId
            || predicted.componentIndex != trace.componentIndex
            || predicted.minusControlCellIndex
                != oriented.minusControlCellIndex
            || predicted.plusControlCellIndex
                != oriented.plusControlCellIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic wall trace-flow topology is invalid");
        }
        ++sharedOrdinal;
    }
}

} // namespace simwing::fsi
