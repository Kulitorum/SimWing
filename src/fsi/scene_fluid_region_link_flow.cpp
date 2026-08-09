#include "scene_fluid_region_link_flow.h"

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
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

fluid::Vector3 add(const fluid::Vector3& first,
                   const fluid::Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

double component(const fluid::Vector3& value,
                 const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return value.x;
    case fluid::GridFaceAxis::Y: return value.y;
    case fluid::GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "scene fluid region link-flow axis is invalid");
}

std::size_t storageBytesForLinks(const std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max()
                    / sizeof(SceneFluidRegionPredictedLinkFlow)) {
        throw std::length_error(
            "scene fluid region link-flow storage size overflows");
    }
    return count * sizeof(SceneFluidRegionPredictedLinkFlow);
}

std::uint64_t predictionFingerprint(
    const SceneFluidRegionLinkFlowPrediction& prediction) {
    Fingerprint fingerprint;
    fingerprint.integer(prediction.version);
    fingerprint.integer(prediction.sourceTransportFingerprint);
    fingerprint.integer(prediction.currentPressureControlVolumeFingerprint);
    fingerprint.integer(prediction.currentPressureFaceLinkFingerprint);
    fingerprint.integer(prediction.currentOpeningFluxFingerprint);
    fingerprint.integer(prediction.currentVelocityFingerprint);
    fingerprint.integer(prediction.previousAcceptedStepCount);
    fingerprint.integer(prediction.currentAcceptedStepCount);
    fingerprint.real(prediction.previousSimulationTimeSeconds);
    fingerprint.real(prediction.currentSimulationTimeSeconds);
    fingerprint.real(prediction.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(prediction.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(prediction.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(prediction.cellCounts.z));
    for (const double value : {
             prediction.lowerMeters.x,
             prediction.lowerMeters.y,
             prediction.lowerMeters.z,
             prediction.upperMeters.x,
             prediction.upperMeters.y,
             prediction.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        prediction.ownedStorageBytes));
    const auto& diagnostics = prediction.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.faceCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.linkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.openingLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.multiLinkFaceCount));
    for (const double value : {
             diagnostics.sourceMomentumKilogramMetersPerSecond.x,
             diagnostics.sourceMomentumKilogramMetersPerSecond.y,
             diagnostics.sourceMomentumKilogramMetersPerSecond.z,
             diagnostics.remappedMomentumKilogramMetersPerSecond.x,
             diagnostics.remappedMomentumKilogramMetersPerSecond.y,
             diagnostics.remappedMomentumKilogramMetersPerSecond.z,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.x,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.y,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.z,
             diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
             diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond,
             diagnostics
                 .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.finite ? 1 : 0));
    fingerprint.integer(static_cast<std::uint64_t>(prediction.links.size()));
    for (const auto& link : prediction.links) {
        fingerprint.integer(static_cast<std::uint64_t>(link.linkIndex));
        fingerprint.integer(link.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(link.faceIndex));
        fingerprint.enumeration(link.kind);
        fingerprint.integer(link.openingPatchStableId);
        fingerprint.real(link.predictedAbsoluteVelocityMetersPerSecond);
        fingerprint.real(
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    return fingerprint.value();
}

void validateGridIdentity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::GridCellCounts counts,
    const fluid::Vector3 lower,
    const fluid::Vector3 upper,
    const char* message) {
    if (counts != grid.cellCounts()
        || lower != grid.lowerMeters()
        || upper != grid.upperMeters()) {
        throw std::invalid_argument(message);
    }
}

double orientedOpeningSweep(
    const SceneFluidPressureFaceLink& link,
    const SceneFluidOpeningFluxSample& sample) {
    if (sample.patchStableId != link.openingPatchStableId
        || sample.openingId != link.openingId
        || sample.areaSquareMeters != link.areaSquareMeters) {
        throw std::invalid_argument(
            "scene fluid region link-flow opening identity is inconsistent");
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
        "scene fluid region link-flow opening orientation is inconsistent");
}

} // namespace

SceneFluidRegionLinkFlowPrediction predictSceneFluidRegionLinkFlows(
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux,
    const SceneFluidRegionLinkFlowLimits& limits) {
    validateSceneFluidRegionTransportIntegrity(transport);
    validateSceneFluidPressureControlVolumeIntegrity(currentPressureVolumes);
    validateSceneFluidOpeningFluxIntegrity(currentOpeningFlux);
    if (!transport.diagnostics.accepted) {
        throw std::invalid_argument(
            "scene fluid region link-flow source transport was not accepted");
    }
    validateGridIdentity(
        grid, currentPressureVolumes.cellCounts,
        currentPressureVolumes.lowerMeters,
        currentPressureVolumes.upperMeters,
        "scene fluid region link-flow pressure grid is foreign");
    validateGridIdentity(
        grid, currentFaceLinks.cellCounts,
        currentFaceLinks.lowerMeters, currentFaceLinks.upperMeters,
        "scene fluid region link-flow face grid is foreign");
    validateGridIdentity(
        grid, currentOpeningFlux.cellCounts,
        currentOpeningFlux.lowerMeters, currentOpeningFlux.upperMeters,
        "scene fluid region link-flow opening grid is foreign");
    if (currentFaceLinks.version != sceneFluidPressureFaceLinkVersion
        || currentFaceLinks.fingerprint == 0
        || currentFaceLinks.pressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || currentFaceLinks.openingPatchFingerprint
            != currentOpeningFlux.openingPatchFingerprint
        || currentFaceLinks.acceptedStepCount
            != transport.acceptedStepCount + 1
        || currentPressureVolumes.acceptedStepCount
            != currentFaceLinks.acceptedStepCount
        || currentOpeningFlux.acceptedStepCount
            != currentFaceLinks.acceptedStepCount
        || currentFaceLinks.simulationTimeSeconds
            != transport.targetSimulationTimeSeconds
        || currentPressureVolumes.simulationTimeSeconds
            != currentFaceLinks.simulationTimeSeconds
        || currentOpeningFlux.simulationTimeSeconds
            != currentFaceLinks.simulationTimeSeconds
        || transport.controlVolumes.size()
            != currentPressureVolumes.controlVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid region link-flow identity is invalid");
    }
    const std::size_t storageBytes = storageBytesForLinks(
        currentFaceLinks.links.size());
    if (currentPressureVolumes.controlVolumes.size()
            > limits.maximumControlVolumes
        || currentFaceLinks.links.size() > limits.maximumLinks
        || storageBytes > limits.maximumLinkFlowBytes) {
        throw std::length_error(
            "scene fluid region link-flow exceeds its limits");
    }

    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*> sampleById;
    for (const auto& sample : currentOpeningFlux.samples) {
        if (sample.patchStableId == 0
            || !sampleById.emplace(sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid region link-flow opening identity is invalid");
        }
    }

    SceneFluidRegionLinkFlowPrediction result;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.currentPressureControlVolumeFingerprint =
        currentPressureVolumes.fingerprint;
    result.currentPressureFaceLinkFingerprint = currentFaceLinks.fingerprint;
    result.currentOpeningFluxFingerprint = currentOpeningFlux.fingerprint;
    result.currentVelocityFingerprint = currentOpeningFlux.velocityFingerprint;
    result.previousAcceptedStepCount = transport.acceptedStepCount;
    result.currentAcceptedStepCount = currentFaceLinks.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        transport.sourceSimulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentFaceLinks.simulationTimeSeconds;
    result.densityKgPerCubicMeter = transport.densityKgPerCubicMeter;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.ownedStorageBytes = storageBytes;
    result.links.reserve(currentFaceLinks.links.size());
    auto& diagnostics = result.diagnostics;
    diagnostics.controlVolumeCount = transport.controlVolumes.size();
    diagnostics.faceCount = currentFaceLinks.faces.size();
    diagnostics.linkCount = currentFaceLinks.links.size();
    diagnostics.sourceMomentumKilogramMetersPerSecond =
        transport.diagnostics.momentumAfterKilogramMetersPerSecond;

    for (std::size_t index = 0;
         index < currentPressureVolumes.controlVolumes.size(); ++index) {
        const auto& current = currentPressureVolumes.controlVolumes[index];
        const auto& source = transport.controlVolumes[index];
        if (current.controlVolumeIndex != index
            || source.controlVolumeIndex != index
            || current.stableId != source.stableId
            || !(current.volumeCubicMeters > 0.0)
            || !(source.volumeCubicMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region link-flow control topology changed");
        }
        diagnostics.maximumAbsoluteVolumeChangeCubicMeters = std::max(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
            std::abs(current.volumeCubicMeters - source.volumeCubicMeters));
        const double mass = result.densityKgPerCubicMeter
            * current.volumeCubicMeters;
        diagnostics.remappedMomentumKilogramMetersPerSecond = add(
            diagnostics.remappedMomentumKilogramMetersPerSecond,
            {
                mass * source.velocityMetersPerSecond.x,
                mass * source.velocityMetersPerSecond.y,
                mass * source.velocityMetersPerSecond.z,
            });
    }
    diagnostics.geometricMomentumChangeKilogramMetersPerSecond = subtract(
        diagnostics.remappedMomentumKilogramMetersPerSecond,
        diagnostics.sourceMomentumKilogramMetersPerSecond);

    std::size_t consumedOpeningSamples = 0;
    for (const auto& face : currentFaceLinks.faces) {
        if (face.linkCount > 1) {
            ++diagnostics.multiLinkFaceCount;
        }
    }
    for (std::size_t index = 0;
         index < currentFaceLinks.links.size(); ++index) {
        const auto& source = currentFaceLinks.links[index];
        if (source.linkIndex != index
            || source.faceIndex >= currentFaceLinks.faces.size()
            || source.minusControlVolumeIndex >= transport.controlVolumes.size()
            || source.plusControlVolumeIndex >= transport.controlVolumes.size()
            || !(source.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region link-flow binding is invalid");
        }
        const auto axis = currentFaceLinks.faces[source.faceIndex].axis;
        const double velocity = 0.5
            * (component(
                   transport.controlVolumes[source.minusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   axis)
               + component(
                   transport.controlVolumes[source.plusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   axis));
        double relativeFlow = source.areaSquareMeters * velocity;
        if (source.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            const auto found = sampleById.find(
                source.openingPatchStableId);
            if (found == sampleById.end()) {
                throw std::invalid_argument(
                    "scene fluid region link-flow opening sample is missing");
            }
            relativeFlow -= orientedOpeningSweep(source, *found->second);
            sampleById.erase(found);
            ++consumedOpeningSamples;
            ++diagnostics.openingLinkCount;
        }
        if (!std::isfinite(velocity) || !std::isfinite(relativeFlow)) {
            throw std::overflow_error(
                "scene fluid region link-flow prediction is non-finite");
        }
        result.links.push_back({
            index,
            source.stableId,
            source.faceIndex,
            source.kind,
            source.openingPatchStableId,
            velocity,
            relativeFlow,
        });
        diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond =
            std::max(
                diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond,
                std::abs(velocity));
        diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                std::abs(relativeFlow));
    }
    if (consumedOpeningSamples != currentOpeningFlux.samples.size()
        || !sampleById.empty()) {
        throw std::invalid_argument(
            "scene fluid region link-flow did not consume every opening sample");
    }
    diagnostics.finite = finite(
            diagnostics.sourceMomentumKilogramMetersPerSecond)
        && finite(diagnostics.remappedMomentumKilogramMetersPerSecond)
        && finite(diagnostics.geometricMomentumChangeKilogramMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters)
        && std::isfinite(
            diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond)
        && std::isfinite(diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene fluid region link-flow diagnostics are non-finite");
    }
    result.fingerprint = predictionFingerprint(result);
    validateSceneFluidRegionLinkFlowPrediction(
        result, transport, grid, currentPressureVolumes,
        currentFaceLinks, currentOpeningFlux);
    return result;
}

void validateSceneFluidRegionLinkFlowPredictionIntegrity(
    const SceneFluidRegionLinkFlowPrediction& prediction) {
    const auto& diagnostics = prediction.diagnostics;
    bool linksValid = true;
    for (std::size_t index = 0; index < prediction.links.size(); ++index) {
        const auto& link = prediction.links[index];
        linksValid = linksValid
            && link.linkIndex == index
            && link.stableId != 0
            && link.faceIndex < diagnostics.faceCount
            && (link.kind == SceneFluidPressureFaceLinkKind::SameRegion
                || (link.kind
                        == SceneFluidPressureFaceLinkKind::AuthoredOpening
                    && link.openingPatchStableId != 0))
            && std::isfinite(
                link.predictedAbsoluteVelocityMetersPerSecond)
            && std::isfinite(
                link.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    if (prediction.version != sceneFluidRegionLinkFlowVersion
        || prediction.fingerprint == 0
        || prediction.sourceTransportFingerprint == 0
        || prediction.currentPressureControlVolumeFingerprint == 0
        || prediction.currentPressureFaceLinkFingerprint == 0
        || prediction.currentOpeningFluxFingerprint == 0
        || prediction.currentVelocityFingerprint == 0
        || prediction.currentAcceptedStepCount
            != prediction.previousAcceptedStepCount + 1
        || !std::isfinite(prediction.previousSimulationTimeSeconds)
        || !std::isfinite(prediction.currentSimulationTimeSeconds)
        || !(prediction.currentSimulationTimeSeconds
             > prediction.previousSimulationTimeSeconds)
        || !std::isfinite(prediction.densityKgPerCubicMeter)
        || !(prediction.densityKgPerCubicMeter > 0.0)
        || prediction.cellCounts.x == 0
        || prediction.cellCounts.y == 0
        || prediction.cellCounts.z == 0
        || !finite(prediction.lowerMeters)
        || !finite(prediction.upperMeters)
        || !(prediction.upperMeters.x > prediction.lowerMeters.x)
        || !(prediction.upperMeters.y > prediction.lowerMeters.y)
        || !(prediction.upperMeters.z > prediction.lowerMeters.z)
        || prediction.ownedStorageBytes
            != storageBytesForLinks(prediction.links.size())
        || diagnostics.controlVolumeCount == 0
        || diagnostics.faceCount == 0
        || diagnostics.linkCount != prediction.links.size()
        || diagnostics.openingLinkCount > diagnostics.linkCount
        || diagnostics.multiLinkFaceCount > diagnostics.faceCount
        || !finite(diagnostics.sourceMomentumKilogramMetersPerSecond)
        || !finite(diagnostics.remappedMomentumKilogramMetersPerSecond)
        || !finite(diagnostics.geometricMomentumChangeKilogramMetersPerSecond)
        || diagnostics.geometricMomentumChangeKilogramMetersPerSecond
            != subtract(
                diagnostics.remappedMomentumKilogramMetersPerSecond,
                diagnostics.sourceMomentumKilogramMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters)
        || diagnostics.maximumAbsoluteVolumeChangeCubicMeters < 0.0
        || !std::isfinite(
            diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond)
        || diagnostics.maximumAbsolutePredictedVelocityMetersPerSecond < 0.0
        || !std::isfinite(diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond)
        || diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond
            < 0.0
        || !diagnostics.finite
        || !linksValid
        || prediction.fingerprint != predictionFingerprint(prediction)) {
        throw std::invalid_argument(
            "scene fluid region link-flow integrity is invalid");
    }
}

void validateSceneFluidRegionLinkFlowPrediction(
    const SceneFluidRegionLinkFlowPrediction& prediction,
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux) {
    validateSceneFluidRegionLinkFlowPredictionIntegrity(prediction);
    validateSceneFluidRegionTransportIntegrity(transport);
    validateSceneFluidOpeningFluxIntegrity(currentOpeningFlux);
    validateGridIdentity(
        grid, prediction.cellCounts,
        prediction.lowerMeters, prediction.upperMeters,
        "scene fluid region link-flow prediction grid is foreign");
    if (prediction.sourceTransportFingerprint != transport.fingerprint
        || prediction.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || prediction.currentPressureFaceLinkFingerprint
            != currentFaceLinks.fingerprint
        || prediction.currentOpeningFluxFingerprint
            != currentOpeningFlux.fingerprint
        || prediction.currentVelocityFingerprint
            != currentOpeningFlux.velocityFingerprint
        || prediction.previousAcceptedStepCount
            != transport.acceptedStepCount
        || prediction.currentAcceptedStepCount
            != currentFaceLinks.acceptedStepCount
        || prediction.currentSimulationTimeSeconds
            != currentFaceLinks.simulationTimeSeconds
        || prediction.densityKgPerCubicMeter
            != transport.densityKgPerCubicMeter
        || prediction.links.size() != currentFaceLinks.links.size()) {
        throw std::invalid_argument(
            "scene fluid region link-flow prediction binding is invalid");
    }
    for (std::size_t index = 0; index < prediction.links.size(); ++index) {
        const auto& predicted = prediction.links[index];
        const auto& source = currentFaceLinks.links[index];
        if (predicted.linkIndex != source.linkIndex
            || predicted.stableId != source.stableId
            || predicted.faceIndex != source.faceIndex
            || predicted.kind != source.kind
            || predicted.openingPatchStableId
                != source.openingPatchStableId) {
            throw std::invalid_argument(
                "scene fluid region link-flow prediction link is foreign");
        }
    }
}

} // namespace simwing::fsi
