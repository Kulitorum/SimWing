#include "scene_fluid_pressure_link_flow.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

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

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1 : 0));
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

std::size_t storageBytesForLinks(const std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max()
                    / sizeof(SceneFluidPressureContinuedLinkFlow)) {
        throw std::length_error(
            "scene pressure link-flow continuation storage size overflows");
    }
    return count * sizeof(SceneFluidPressureContinuedLinkFlow);
}

std::uint64_t continuationFingerprint(
    const SceneFluidPressureLinkFlowContinuation& continuation) {
    Fingerprint fingerprint;
    fingerprint.integer(continuation.version);
    fingerprint.integer(
        continuation.previousPressureProjectionFingerprint);
    fingerprint.integer(continuation.previousPressureFaceLinkFingerprint);
    fingerprint.integer(continuation.previousOpeningPatchFingerprint);
    fingerprint.integer(continuation.currentPressureFaceLinkFingerprint);
    fingerprint.integer(continuation.currentOpeningFluxFingerprint);
    fingerprint.integer(continuation.currentVelocityFingerprint);
    fingerprint.integer(continuation.previousAcceptedStepCount);
    fingerprint.integer(continuation.currentAcceptedStepCount);
    fingerprint.real(continuation.previousSimulationTimeSeconds);
    fingerprint.real(continuation.currentSimulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.cellCounts.z));
    for (const double value : {
             continuation.lowerMeters.x,
             continuation.lowerMeters.y,
             continuation.lowerMeters.z,
             continuation.upperMeters.x,
             continuation.upperMeters.y,
             continuation.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.ownedStorageBytes));
    const auto& diagnostics = continuation.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.faceCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.linkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.openingLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.multiLinkFaceCount));
    fingerprint.real(
        diagnostics.maximumCarriedAbsoluteVelocityDeviationMetersPerSecond);
    fingerprint.real(
        diagnostics.maximumAreaRenormalizationMetersPerSecond);
    fingerprint.real(diagnostics
        .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond);
    fingerprint.real(
        diagnostics.maximumAbsoluteFaceFlowClosureCubicMetersPerSecond);
    fingerprint.boolean(diagnostics.finite);
    fingerprint.integer(static_cast<std::uint64_t>(
        continuation.links.size()));
    for (const auto& link : continuation.links) {
        fingerprint.integer(static_cast<std::uint64_t>(link.linkIndex));
        fingerprint.integer(link.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(link.faceIndex));
        fingerprint.enumeration(link.kind);
        fingerprint.integer(link.openingPatchStableId);
        fingerprint.real(
            link.carriedAbsoluteVelocityDeviationMetersPerSecond);
        fingerprint.real(
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    return fingerprint.value();
}

double faceVelocity(
    const SceneFluidPressureFace& face,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    const std::size_t index = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case fluid::GridFaceAxis::X: return velocity.xFaces()[index];
    case fluid::GridFaceAxis::Y: return velocity.yFaces()[index];
    case fluid::GridFaceAxis::Z: return velocity.zFaces()[index];
    }
    throw std::invalid_argument(
        "scene pressure link-flow continuation face axis is invalid");
}

double axisComponent(const Vec3& value, const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return value.x;
    case fluid::GridFaceAxis::Y: return value.y;
    case fluid::GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "scene pressure link-flow continuation patch axis is invalid");
}

double orientedOpeningFlow(
    const SceneFluidPressureFaceLink& link,
    const SceneFluidOpeningFluxSample& sample,
    const double SceneFluidOpeningFluxSample::* member) {
    if (sample.patchStableId != link.openingPatchStableId
        || sample.openingId != link.openingId
        || sample.areaSquareMeters != link.areaSquareMeters) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation opening identity is inconsistent");
    }
    const double value = sample.*member;
    if (sample.negativeSideRegionId == link.minusRegionId
        && sample.positiveSideRegionId == link.plusRegionId) {
        return value;
    }
    if (sample.negativeSideRegionId == link.plusRegionId
        && sample.positiveSideRegionId == link.minusRegionId) {
        return -value;
    }
    throw std::invalid_argument(
        "scene pressure link-flow continuation opening orientation is inconsistent");
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

} // namespace

SceneFluidPressureLinkFlowContinuation
continueSceneFluidPressureLinkFlows(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureFaceLinkSet& previousFaceLinks,
    const SceneFluidOpeningGridPatchSet& previousOpeningPatches,
    const SceneFluidPressureProjection& previousProjection,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux,
    const fluid::MacVelocityField& currentBulkVelocityMetersPerSecond,
    const SceneFluidPressureLinkFlowContinuationLimits& limits) {
    validateSceneFluidPressureProjectionIntegrity(previousProjection);
    validateSceneFluidOpeningFluxIntegrity(currentOpeningFlux);
    if (!currentBulkVelocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(currentBulkVelocityMetersPerSecond)
        || currentOpeningFlux.velocityFingerprint
            != sceneFluidOpeningFluxVelocityFingerprint(
                grid, currentBulkVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation bulk velocity is invalid");
    }
    validateGridIdentity(
        grid, previousFaceLinks.cellCounts,
        previousFaceLinks.lowerMeters, previousFaceLinks.upperMeters,
        "scene pressure link-flow continuation previous grid is foreign");
    validateGridIdentity(
        grid, previousOpeningPatches.cellCounts,
        previousOpeningPatches.lowerMeters,
        previousOpeningPatches.upperMeters,
        "scene pressure link-flow continuation previous opening grid is foreign");
    validateGridIdentity(
        grid, currentFaceLinks.cellCounts,
        currentFaceLinks.lowerMeters, currentFaceLinks.upperMeters,
        "scene pressure link-flow continuation current grid is foreign");
    validateGridIdentity(
        grid, currentOpeningFlux.cellCounts,
        currentOpeningFlux.lowerMeters, currentOpeningFlux.upperMeters,
        "scene pressure link-flow continuation current opening grid is foreign");
    if (previousFaceLinks.version != sceneFluidPressureFaceLinkVersion
        || previousFaceLinks.fingerprint == 0
        || currentFaceLinks.version != sceneFluidPressureFaceLinkVersion
        || currentFaceLinks.fingerprint == 0
        || previousFaceLinks.surfaceDefinitionFingerprint
            != currentFaceLinks.surfaceDefinitionFingerprint
        || previousFaceLinks.structureDefinitionFingerprint
            != currentFaceLinks.structureDefinitionFingerprint
        || previousFaceLinks.settings != currentFaceLinks.settings
        || !previousProjection.diagnostics.accepted
        || previousProjection.pressureFaceLinkFingerprint
            != previousFaceLinks.fingerprint
        || previousProjection.acceptedStepCount
            != previousFaceLinks.acceptedStepCount
        || previousProjection.simulationTimeSeconds
            != previousFaceLinks.simulationTimeSeconds
        || previousProjection.links.size() != previousFaceLinks.links.size()
        || previousOpeningPatches.version
            != sceneFluidOpeningGridPatchVersion
        || previousOpeningPatches.fingerprint == 0
        || previousOpeningPatches.surfaceStateFingerprint
            != previousFaceLinks.surfaceStateFingerprint
        || previousOpeningPatches.fingerprint
            != previousFaceLinks.openingPatchFingerprint
        || previousOpeningPatches.acceptedStepCount
            != previousFaceLinks.acceptedStepCount
        || previousOpeningPatches.simulationTimeSeconds
            != previousFaceLinks.simulationTimeSeconds
        || currentOpeningFlux.openingPatchFingerprint
            != currentFaceLinks.openingPatchFingerprint
        || currentOpeningFlux.acceptedStepCount
            != currentFaceLinks.acceptedStepCount
        || currentOpeningFlux.simulationTimeSeconds
            != currentFaceLinks.simulationTimeSeconds
        || currentOpeningFlux.surfaceStateFingerprint
            != currentFaceLinks.surfaceStateFingerprint
        || currentFaceLinks.acceptedStepCount
            != previousFaceLinks.acceptedStepCount + 1
        || !(currentFaceLinks.simulationTimeSeconds
             > previousFaceLinks.simulationTimeSeconds)
        || previousFaceLinks.faces.size() != currentFaceLinks.faces.size()
        || previousFaceLinks.links.size() != currentFaceLinks.links.size()) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation identity is invalid");
    }
    const std::size_t continuationBytes =
        storageBytesForLinks(currentFaceLinks.links.size());
    if (currentFaceLinks.links.size() > limits.maximumLinks
        || continuationBytes > limits.maximumContinuationBytes) {
        throw std::length_error(
            "scene pressure link-flow continuation exceeds its limits");
    }

    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*>
        previousPatchById;
    for (const auto& patch : previousOpeningPatches.patches) {
        if (patch.stableId == 0
            || !previousPatchById.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene pressure link-flow continuation previous opening identity is invalid");
        }
    }
    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        currentSampleById;
    for (const auto& sample : currentOpeningFlux.samples) {
        if (sample.patchStableId == 0
            || !currentSampleById.emplace(
                sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene pressure link-flow continuation current opening identity is invalid");
        }
    }

    SceneFluidPressureLinkFlowContinuation result;
    result.previousPressureProjectionFingerprint =
        previousProjection.fingerprint;
    result.previousPressureFaceLinkFingerprint = previousFaceLinks.fingerprint;
    result.previousOpeningPatchFingerprint = previousOpeningPatches.fingerprint;
    result.currentPressureFaceLinkFingerprint = currentFaceLinks.fingerprint;
    result.currentOpeningFluxFingerprint = currentOpeningFlux.fingerprint;
    result.currentVelocityFingerprint = currentOpeningFlux.velocityFingerprint;
    result.previousAcceptedStepCount = previousFaceLinks.acceptedStepCount;
    result.currentAcceptedStepCount = currentFaceLinks.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousFaceLinks.simulationTimeSeconds;
    result.currentSimulationTimeSeconds = currentFaceLinks.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.links.reserve(currentFaceLinks.links.size());
    auto& diagnostics = result.diagnostics;
    diagnostics.faceCount = currentFaceLinks.faces.size();
    diagnostics.linkCount = currentFaceLinks.links.size();

    for (std::size_t faceIndex = 0;
         faceIndex < currentFaceLinks.faces.size(); ++faceIndex) {
        const auto& previousFace = previousFaceLinks.faces[faceIndex];
        const auto& currentFace = currentFaceLinks.faces[faceIndex];
        if (previousFace.faceIndex != faceIndex
            || currentFace.faceIndex != faceIndex
            || previousFace.stableId != currentFace.stableId
            || previousFace.axis != currentFace.axis
            || previousFace.i != currentFace.i
            || previousFace.j != currentFace.j
            || previousFace.k != currentFace.k
            || previousFace.minusCellIndex != currentFace.minusCellIndex
            || previousFace.plusCellIndex != currentFace.plusCellIndex
            || previousFace.status != currentFace.status
            || previousFace.linkCount == 0
            || previousFace.firstLink != currentFace.firstLink
            || previousFace.linkCount != currentFace.linkCount
            || previousFace.firstLink > previousFaceLinks.links.size()
            || previousFace.linkCount
                > previousFaceLinks.links.size() - previousFace.firstLink
            || currentFace.firstLink > currentFaceLinks.links.size()
            || currentFace.linkCount
                > currentFaceLinks.links.size() - currentFace.firstLink
            || !(previousFace.faceAreaSquareMeters > 0.0)
            || currentFace.faceAreaSquareMeters
                != previousFace.faceAreaSquareMeters) {
            throw std::invalid_argument(
                "scene pressure link-flow continuation face topology changed");
        }
        if (currentFace.linkCount > 1) {
            ++diagnostics.multiLinkFaceCount;
        }

        std::vector<double> previousDeviations;
        previousDeviations.reserve(currentFace.linkCount);
        double previousAbsoluteFlow = 0.0;
        for (std::size_t offset = 0;
             offset < currentFace.linkCount; ++offset) {
            const std::size_t previousIndex = previousFace.firstLink + offset;
            const std::size_t currentIndex = currentFace.firstLink + offset;
            const auto& previousLink = previousFaceLinks.links[previousIndex];
            const auto& currentLink = currentFaceLinks.links[currentIndex];
            const auto& projected = previousProjection.links[previousIndex];
            if (previousLink.linkIndex != previousIndex
                || currentLink.linkIndex != currentIndex
                || projected.linkIndex != previousIndex
                || previousLink.stableId != currentLink.stableId
                || projected.stableId != previousLink.stableId
                || previousLink.kind != currentLink.kind
                || projected.kind != previousLink.kind
                || previousLink.minusRegionId != currentLink.minusRegionId
                || previousLink.plusRegionId != currentLink.plusRegionId
                || previousLink.openingId != currentLink.openingId
                || previousLink.openingPatchStableId
                    != currentLink.openingPatchStableId
                || previousLink.minusRegionIndex
                    != currentLink.minusRegionIndex
                || previousLink.plusRegionIndex
                    != currentLink.plusRegionIndex
                || previousLink.componentIndex != currentLink.componentIndex
                || !(previousLink.areaSquareMeters > 0.0)
                || !(currentLink.areaSquareMeters > 0.0)) {
                throw std::invalid_argument(
                    "scene pressure link-flow continuation link topology changed");
            }
            double absoluteFlow = projected
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond;
            if (previousLink.kind
                == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
                const auto found = previousPatchById.find(
                    previousLink.openingPatchStableId);
                if (found == previousPatchById.end()) {
                    throw std::invalid_argument(
                        "scene pressure link-flow continuation previous opening patch is missing");
                }
                const auto& patch = *found->second;
                const double normal = axisComponent(
                    patch.unitNormalNegativeToPositive,
                    previousFace.axis);
                const bool forwardRegions =
                    patch.negativeSideRegionId == previousLink.minusRegionId
                    && patch.positiveSideRegionId
                        == previousLink.plusRegionId;
                const bool reverseRegions =
                    patch.negativeSideRegionId == previousLink.plusRegionId
                    && patch.positiveSideRegionId
                        == previousLink.minusRegionId;
                if (patch.openingId != previousLink.openingId
                    || patch.areaSquareMeters
                        != previousLink.areaSquareMeters
                    || (!forwardRegions && !reverseRegions)
                    || std::abs(std::abs(normal) - 1.0) > 1.0e-10) {
                    throw std::invalid_argument(
                        "scene pressure link-flow continuation previous opening patch is foreign");
                }
                const double orientation = forwardRegions ? 1.0 : -1.0;
                if (orientation * normal < 1.0 - 1.0e-10) {
                    throw std::invalid_argument(
                        "scene pressure link-flow continuation previous opening orientation is inconsistent");
                }
                absoluteFlow += orientation
                    * patch.surfaceSweepRateCubicMetersPerSecond;
                ++diagnostics.openingLinkCount;
            }
            previousAbsoluteFlow += absoluteFlow;
            previousDeviations.push_back(
                absoluteFlow / previousLink.areaSquareMeters);
        }
        const double previousMeanVelocity =
            previousAbsoluteFlow / previousFace.faceAreaSquareMeters;
        double currentWeightedDeviation = 0.0;
        for (std::size_t offset = 0;
             offset < currentFace.linkCount; ++offset) {
            previousDeviations[offset] -= previousMeanVelocity;
            currentWeightedDeviation +=
                currentFaceLinks.links[currentFace.firstLink + offset]
                    .areaSquareMeters
                * previousDeviations[offset];
        }
        const double renormalization =
            currentWeightedDeviation / currentFace.faceAreaSquareMeters;
        diagnostics.maximumAreaRenormalizationMetersPerSecond = std::max(
            diagnostics.maximumAreaRenormalizationMetersPerSecond,
            std::abs(renormalization));

        const double bulkFaceVelocity = faceVelocity(
            currentFace, grid, currentBulkVelocityMetersPerSecond);
        double predictedAbsoluteFaceFlow = 0.0;
        for (std::size_t offset = 0;
             offset < currentFace.linkCount; ++offset) {
            const std::size_t currentIndex = currentFace.firstLink + offset;
            const auto& currentLink = currentFaceLinks.links[currentIndex];
            const double carriedDeviation =
                previousDeviations[offset] - renormalization;
            double predictedRelativeFlow =
                currentLink.areaSquareMeters
                * (bulkFaceVelocity + carriedDeviation);
            double currentOpeningSweep = 0.0;
            if (currentLink.kind
                == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
                const auto found = currentSampleById.find(
                    currentLink.openingPatchStableId);
                if (found == currentSampleById.end()) {
                    throw std::invalid_argument(
                        "scene pressure link-flow continuation current opening sample is missing");
                }
                predictedRelativeFlow = orientedOpeningFlow(
                    currentLink, *found->second,
                    &SceneFluidOpeningFluxSample::
                        relativeVolumeFlowRateCubicMetersPerSecond)
                    + currentLink.areaSquareMeters * carriedDeviation;
                currentOpeningSweep = orientedOpeningFlow(
                    currentLink, *found->second,
                    &SceneFluidOpeningFluxSample::
                        surfaceSweepRateCubicMetersPerSecond);
            }
            SceneFluidPressureContinuedLinkFlow link;
            link.linkIndex = currentIndex;
            link.stableId = currentLink.stableId;
            link.faceIndex = currentLink.faceIndex;
            link.kind = currentLink.kind;
            link.openingPatchStableId = currentLink.openingPatchStableId;
            link.carriedAbsoluteVelocityDeviationMetersPerSecond =
                carriedDeviation;
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
                predictedRelativeFlow;
            result.links.push_back(link);
            predictedAbsoluteFaceFlow +=
                predictedRelativeFlow + currentOpeningSweep;
            diagnostics
                .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond =
                std::max(
                    diagnostics
                        .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond,
                    std::abs(carriedDeviation));
            diagnostics
                .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond =
                std::max(
                    diagnostics
                        .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond,
                    std::abs(predictedRelativeFlow));
        }
        diagnostics.maximumAbsoluteFaceFlowClosureCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteFaceFlowClosureCubicMetersPerSecond,
                std::abs(
                    predictedAbsoluteFaceFlow
                    - bulkFaceVelocity * currentFace.faceAreaSquareMeters));
    }
    if (result.links.size() != currentFaceLinks.links.size()
        || diagnostics.openingLinkCount != currentOpeningFlux.samples.size()) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation did not consume every link");
    }
    result.ownedStorageBytes = continuationBytes;
    diagnostics.finite = std::isfinite(
            diagnostics
                .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAreaRenormalizationMetersPerSecond)
        && std::isfinite(diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsoluteFaceFlowClosureCubicMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene pressure link-flow continuation is non-finite");
    }
    result.fingerprint = continuationFingerprint(result);
    validateSceneFluidPressureLinkFlowContinuation(
        result, grid, currentFaceLinks, currentOpeningFlux);
    return result;
}

void validateSceneFluidPressureLinkFlowContinuationIntegrity(
    const SceneFluidPressureLinkFlowContinuation& continuation) {
    const auto& diagnostics = continuation.diagnostics;
    std::unordered_set<std::uint64_t> stableIds;
    stableIds.reserve(continuation.links.size());
    bool linksFinite = true;
    for (std::size_t index = 0; index < continuation.links.size(); ++index) {
        const auto& link = continuation.links[index];
        linksFinite = linksFinite
            && link.linkIndex == index
            && link.stableId != 0
            && link.faceIndex < diagnostics.faceCount
            && (link.kind == SceneFluidPressureFaceLinkKind::SameRegion
                || (link.kind
                        == SceneFluidPressureFaceLinkKind::AuthoredOpening
                    && link.openingPatchStableId != 0))
            && std::isfinite(
                link.carriedAbsoluteVelocityDeviationMetersPerSecond)
            && std::isfinite(
                link.predictedRelativeVolumeFlowRateCubicMetersPerSecond)
            && stableIds.insert(link.stableId).second;
    }
    if (continuation.version
            != sceneFluidPressureLinkFlowContinuationVersion
        || continuation.fingerprint == 0
        || continuation.previousPressureProjectionFingerprint == 0
        || continuation.previousPressureFaceLinkFingerprint == 0
        || continuation.previousOpeningPatchFingerprint == 0
        || continuation.currentPressureFaceLinkFingerprint == 0
        || continuation.currentOpeningFluxFingerprint == 0
        || continuation.currentVelocityFingerprint == 0
        || continuation.currentAcceptedStepCount
            != continuation.previousAcceptedStepCount + 1
        || !std::isfinite(continuation.previousSimulationTimeSeconds)
        || !std::isfinite(continuation.currentSimulationTimeSeconds)
        || !(continuation.currentSimulationTimeSeconds
             > continuation.previousSimulationTimeSeconds)
        || continuation.cellCounts.x == 0
        || continuation.cellCounts.y == 0
        || continuation.cellCounts.z == 0
        || !std::isfinite(continuation.lowerMeters.x)
        || !std::isfinite(continuation.lowerMeters.y)
        || !std::isfinite(continuation.lowerMeters.z)
        || !std::isfinite(continuation.upperMeters.x)
        || !std::isfinite(continuation.upperMeters.y)
        || !std::isfinite(continuation.upperMeters.z)
        || !(continuation.upperMeters.x > continuation.lowerMeters.x)
        || !(continuation.upperMeters.y > continuation.lowerMeters.y)
        || !(continuation.upperMeters.z > continuation.lowerMeters.z)
        || continuation.ownedStorageBytes
            != storageBytesForLinks(continuation.links.size())
        || diagnostics.faceCount == 0
        || diagnostics.linkCount != continuation.links.size()
        || diagnostics.openingLinkCount > diagnostics.linkCount
        || diagnostics.multiLinkFaceCount > diagnostics.faceCount
        || !diagnostics.finite
        || !std::isfinite(diagnostics
            .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond)
        || diagnostics
            .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond < 0.0
        || !std::isfinite(
            diagnostics.maximumAreaRenormalizationMetersPerSecond)
        || diagnostics.maximumAreaRenormalizationMetersPerSecond < 0.0
        || !std::isfinite(diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond)
        || diagnostics
            .maximumAbsolutePredictedRelativeVolumeFlowRateCubicMetersPerSecond
            < 0.0
        || !std::isfinite(
            diagnostics.maximumAbsoluteFaceFlowClosureCubicMetersPerSecond)
        || diagnostics.maximumAbsoluteFaceFlowClosureCubicMetersPerSecond
            < 0.0
        || !linksFinite
        || continuation.fingerprint
            != continuationFingerprint(continuation)) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation integrity is invalid");
    }
}

void validateSceneFluidPressureLinkFlowContinuation(
    const SceneFluidPressureLinkFlowContinuation& continuation,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureFaceLinkSet& currentFaceLinks,
    const SceneFluidOpeningFluxSet& currentOpeningFlux) {
    validateSceneFluidPressureLinkFlowContinuationIntegrity(continuation);
    validateSceneFluidOpeningFluxIntegrity(currentOpeningFlux);
    validateGridIdentity(
        grid, continuation.cellCounts,
        continuation.lowerMeters, continuation.upperMeters,
        "scene pressure link-flow continuation grid is foreign");
    if (continuation.currentPressureFaceLinkFingerprint
            != currentFaceLinks.fingerprint
        || continuation.currentOpeningFluxFingerprint
            != currentOpeningFlux.fingerprint
        || continuation.currentVelocityFingerprint
            != currentOpeningFlux.velocityFingerprint
        || continuation.currentAcceptedStepCount
            != currentFaceLinks.acceptedStepCount
        || continuation.currentAcceptedStepCount
            != currentOpeningFlux.acceptedStepCount
        || continuation.currentSimulationTimeSeconds
            != currentFaceLinks.simulationTimeSeconds
        || continuation.currentSimulationTimeSeconds
            != currentOpeningFlux.simulationTimeSeconds
        || continuation.links.size() != currentFaceLinks.links.size()) {
        throw std::invalid_argument(
            "scene pressure link-flow continuation current identity is invalid");
    }
    for (std::size_t index = 0;
         index < continuation.links.size(); ++index) {
        const auto& continued = continuation.links[index];
        const auto& source = currentFaceLinks.links[index];
        if (continued.linkIndex != source.linkIndex
            || continued.stableId != source.stableId
            || continued.faceIndex != source.faceIndex
            || continued.kind != source.kind
            || continued.openingPatchStableId
                != source.openingPatchStableId) {
            throw std::invalid_argument(
                "scene pressure link-flow continuation binding is invalid");
        }
    }
}

} // namespace simwing::fsi
