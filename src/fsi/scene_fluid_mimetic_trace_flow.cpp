#include "scene_fluid_mimetic_trace_flow.h"

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
    fingerprint.integer(prediction.structureDefinitionFingerprint);
    fingerprint.integer(prediction.acceptedStepCount);
    fingerprint.real(prediction.simulationTimeSeconds);
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

void validateSceneFluidMimeticTraceFlowPredictionIntegrity(
    const SceneFluidMimeticTraceFlowPrediction& prediction) {
    if (prediction.version != sceneFluidMimeticTraceFlowVersion
        || prediction.fingerprint == 0
        || prediction.mimeticControlCellFingerprint == 0
        || prediction.mimeticTraceSystemFingerprint == 0
        || prediction.pressureFaceLinkFingerprint == 0
        || prediction.openingFluxFingerprint == 0
        || prediction.velocityFingerprint == 0
        || prediction.structureDefinitionFingerprint == 0
        || !std::isfinite(prediction.simulationTimeSeconds)
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

} // namespace simwing::fsi
