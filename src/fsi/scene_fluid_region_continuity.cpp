#include "scene_fluid_region_continuity.h"

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

void validateSettings(const SceneFluidRegionContinuitySettings& settings) {
    if (!std::isfinite(settings.absoluteVolumeToleranceCubicMeters)
        || settings.absoluteVolumeToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeVolumeTolerance)
        || settings.relativeVolumeTolerance < 0.0
        || (settings.absoluteVolumeToleranceCubicMeters == 0.0
            && settings.relativeVolumeTolerance == 0.0)) {
        throw std::invalid_argument(
            "scene fluid region-continuity tolerances are invalid");
    }
}

bool sameGrid(const SceneFluidCellVolumeSet& volumes,
              const SceneFluidOpeningFluxSet& flux) {
    return volumes.cellCounts == flux.cellCounts
        && volumes.lowerMeters == flux.lowerMeters
        && volumes.upperMeters == flux.upperMeters;
}

bool sameGrid(const SceneFluidCellVolumeSet& first,
              const SceneFluidCellVolumeSet& second) {
    return first.cellCounts == second.cellCounts
        && first.lowerMeters == second.lowerMeters
        && first.upperMeters == second.upperMeters;
}

bool sameOpeningTopology(const SceneFluidOpeningFluxSet& first,
                         const SceneFluidOpeningFluxSet& second) {
    if (first.openings.size() != second.openings.size()) return false;
    for (std::size_t index = 0; index < first.openings.size(); ++index) {
        const auto& a = first.openings[index];
        const auto& b = second.openings[index];
        if (a.openingIndex != b.openingIndex
            || a.openingId != b.openingId
            || a.negativeSideRegionId != b.negativeSideRegionId
            || a.positiveSideRegionId != b.positiveSideRegionId
            || a.role != b.role) {
            return false;
        }
    }
    return true;
}

void validateSources(const SceneFluidCellVolumeSet& previousVolumes,
                     const SceneFluidCellVolumeSet& currentVolumes,
                     const SceneFluidOpeningFluxSet& previousFlux,
                     const SceneFluidOpeningFluxSet& currentFlux) {
    validateSceneFluidCellVolumeIntegrity(previousVolumes);
    validateSceneFluidCellVolumeIntegrity(currentVolumes);
    validateSceneFluidOpeningFluxIntegrity(previousFlux);
    validateSceneFluidOpeningFluxIntegrity(currentFlux);
    if (previousVolumes.version != sceneFluidCellVolumeVersion
        || currentVolumes.version != sceneFluidCellVolumeVersion
        || previousFlux.version != sceneFluidOpeningFluxVersion
        || currentFlux.version != sceneFluidOpeningFluxVersion
        || previousVolumes.fingerprint == 0
        || currentVolumes.fingerprint == 0
        || previousFlux.fingerprint == 0
        || currentFlux.fingerprint == 0
        || previousVolumes.surfaceDefinitionFingerprint == 0
        || previousVolumes.structureDefinitionFingerprint == 0
        || previousVolumes.surfaceStateFingerprint == 0
        || currentVolumes.surfaceStateFingerprint == 0
        || previousVolumes.surfaceDefinitionFingerprint
            != currentVolumes.surfaceDefinitionFingerprint
        || previousVolumes.surfaceDefinitionFingerprint
            != previousFlux.surfaceDefinitionFingerprint
        || previousVolumes.surfaceDefinitionFingerprint
            != currentFlux.surfaceDefinitionFingerprint
        || previousVolumes.structureDefinitionFingerprint
            != currentVolumes.structureDefinitionFingerprint
        || previousVolumes.structureDefinitionFingerprint
            != previousFlux.structureDefinitionFingerprint
        || previousVolumes.structureDefinitionFingerprint
            != currentFlux.structureDefinitionFingerprint
        || previousVolumes.surfaceStateFingerprint
            != previousFlux.surfaceStateFingerprint
        || currentVolumes.surfaceStateFingerprint
            != currentFlux.surfaceStateFingerprint
        || previousVolumes.acceptedStepCount
            != previousFlux.acceptedStepCount
        || currentVolumes.acceptedStepCount
            != currentFlux.acceptedStepCount
        || previousVolumes.simulationTimeSeconds
            != previousFlux.simulationTimeSeconds
        || currentVolumes.simulationTimeSeconds
            != currentFlux.simulationTimeSeconds
        || previousVolumes.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentVolumes.acceptedStepCount
            != previousVolumes.acceptedStepCount + 1
        || !std::isfinite(previousVolumes.simulationTimeSeconds)
        || !std::isfinite(currentVolumes.simulationTimeSeconds)
        || !(currentVolumes.simulationTimeSeconds
             > previousVolumes.simulationTimeSeconds)
        || previousVolumes.settings != currentVolumes.settings
        || previousVolumes.outsideRegionId != currentVolumes.outsideRegionId
        || !sameGrid(previousVolumes, currentVolumes)
        || !sameGrid(previousVolumes, previousFlux)
        || !sameGrid(currentVolumes, currentFlux)
        || !sameOpeningTopology(previousFlux, currentFlux)
        || previousVolumes.regionVolumes.size()
            != currentVolumes.regionVolumes.size()
        || previousVolumes.regionVolumes.size()
            != previousFlux.regions.size()
        || previousVolumes.regionVolumes.size()
            != currentFlux.regions.size()) {
        throw std::invalid_argument(
            "scene fluid region-continuity source identity is invalid");
    }
    for (std::size_t index = 0;
         index < previousVolumes.regionVolumes.size(); ++index) {
        const auto& previousVolume = previousVolumes.regionVolumes[index];
        const auto& currentVolume = currentVolumes.regionVolumes[index];
        const auto& previousRegionFlux = previousFlux.regions[index];
        const auto& currentRegionFlux = currentFlux.regions[index];
        if (previousRegionFlux.regionIndex != index
            || currentRegionFlux.regionIndex != index
            || previousVolume.regionId != currentVolume.regionId
            || previousVolume.regionId != previousRegionFlux.regionId
            || previousVolume.regionId != currentRegionFlux.regionId
            || previousRegionFlux.kind != currentRegionFlux.kind
            || !std::isfinite(previousVolume.summedCellVolumeCubicMeters)
            || !std::isfinite(currentVolume.summedCellVolumeCubicMeters)
            || previousVolume.summedCellVolumeCubicMeters < 0.0
            || currentVolume.summedCellVolumeCubicMeters < 0.0
            || !std::isfinite(
                previousRegionFlux
                    .outwardRelativeVolumeFlowRateCubicMetersPerSecond)
            || !std::isfinite(
                currentRegionFlux
                    .outwardRelativeVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid region-continuity region sources are invalid");
        }
    }
}

std::size_t storageBytesForCount(const std::size_t regionCount) {
    if (regionCount != 0
        && regionCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidRegionContinuity)) {
        throw std::length_error(
            "scene fluid region-continuity storage size overflows");
    }
    return regionCount * sizeof(SceneFluidRegionContinuity);
}

std::uint64_t continuityFingerprint(
    const SceneFluidRegionContinuitySet& continuity) {
    Fingerprint fingerprint;
    fingerprint.integer(continuity.version);
    for (const std::uint64_t value : {
             continuity.surfaceDefinitionFingerprint,
             continuity.structureDefinitionFingerprint,
             continuity.previousSurfaceStateFingerprint,
             continuity.currentSurfaceStateFingerprint,
             continuity.previousCellVolumeFingerprint,
             continuity.currentCellVolumeFingerprint,
             continuity.previousOpeningFluxFingerprint,
             continuity.currentOpeningFluxFingerprint,
             continuity.previousAcceptedStepCount,
             continuity.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    for (const double value : {
             continuity.previousSimulationTimeSeconds,
             continuity.currentSimulationTimeSeconds,
             continuity.durationSeconds,
             continuity.lowerMeters.x,
             continuity.lowerMeters.y,
             continuity.lowerMeters.z,
             continuity.upperMeters.x,
             continuity.upperMeters.y,
             continuity.upperMeters.z,
             continuity.settings.absoluteVolumeToleranceCubicMeters,
             continuity.settings.relativeVolumeTolerance,
             continuity.maximumAbsoluteContinuityResidualCubicMeters,
             continuity.globalGeometryVolumeChangeCubicMeters,
             continuity.globalIntegratedOutwardRelativeVolumeCubicMeters,
             continuity.globalContinuityResidualCubicMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(continuity.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(continuity.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(continuity.cellCounts.z));
    fingerprint.integer(static_cast<std::uint64_t>(
        continuity.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        continuity.failedRegionCount));
    fingerprint.integer(static_cast<std::uint8_t>(
        continuity.allRegionsWithinTolerance));
    fingerprint.integer(static_cast<std::uint64_t>(continuity.regions.size()));
    for (const auto& region : continuity.regions) {
        fingerprint.integer(static_cast<std::uint64_t>(region.regionIndex));
        fingerprint.integer(region.regionId);
        fingerprint.enumeration(region.kind);
        for (const double value : {
                 region.previousVolumeCubicMeters,
                 region.currentVolumeCubicMeters,
                 region.geometryVolumeChangeCubicMeters,
                 region.previousOutwardRelativeFlowRateCubicMetersPerSecond,
                 region.currentOutwardRelativeFlowRateCubicMetersPerSecond,
                 region.integratedOutwardRelativeVolumeCubicMeters,
                 region.continuityResidualCubicMeters,
                 region.toleranceCubicMeters}) {
            fingerprint.real(value);
        }
        fingerprint.integer(static_cast<std::uint8_t>(
            region.withinTolerance));
    }
    return fingerprint.value();
}

SceneFluidRegionContinuitySet buildContinuity(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionContinuitySettings& settings,
    const SceneFluidRegionContinuityLimits& limits) {
    validateSettings(settings);
    validateSources(
        previousVolumes, currentVolumes, previousFlux, currentFlux);
    const std::size_t regionCount = previousVolumes.regionVolumes.size();
    if (regionCount > limits.maximumRegions) {
        throw std::length_error(
            "scene fluid region-continuity exceeds its region limit");
    }
    const std::size_t storageBytes = storageBytesForCount(regionCount);
    if (storageBytes > limits.maximumContinuityBytes) {
        throw std::length_error(
            "scene fluid region-continuity exceeds its byte limit");
    }

    SceneFluidRegionContinuitySet result;
    result.surfaceDefinitionFingerprint =
        previousVolumes.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        previousVolumes.structureDefinitionFingerprint;
    result.previousSurfaceStateFingerprint =
        previousVolumes.surfaceStateFingerprint;
    result.currentSurfaceStateFingerprint =
        currentVolumes.surfaceStateFingerprint;
    result.previousCellVolumeFingerprint = previousVolumes.fingerprint;
    result.currentCellVolumeFingerprint = currentVolumes.fingerprint;
    result.previousOpeningFluxFingerprint = previousFlux.fingerprint;
    result.currentOpeningFluxFingerprint = currentFlux.fingerprint;
    result.previousAcceptedStepCount = previousVolumes.acceptedStepCount;
    result.currentAcceptedStepCount = currentVolumes.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousVolumes.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentVolumes.simulationTimeSeconds;
    result.durationSeconds = result.currentSimulationTimeSeconds
        - result.previousSimulationTimeSeconds;
    result.cellCounts = previousVolumes.cellCounts;
    result.lowerMeters = previousVolumes.lowerMeters;
    result.upperMeters = previousVolumes.upperMeters;
    result.settings = settings;
    result.ownedStorageBytes = storageBytes;
    result.regions.reserve(regionCount);

    for (std::size_t index = 0; index < regionCount; ++index) {
        const auto& previousVolume = previousVolumes.regionVolumes[index];
        const auto& currentVolume = currentVolumes.regionVolumes[index];
        const auto& previousRegionFlux = previousFlux.regions[index];
        const auto& currentRegionFlux = currentFlux.regions[index];
        SceneFluidRegionContinuity region;
        region.regionIndex = index;
        region.regionId = previousVolume.regionId;
        region.kind = previousRegionFlux.kind;
        region.previousVolumeCubicMeters =
            previousVolume.summedCellVolumeCubicMeters;
        region.currentVolumeCubicMeters =
            currentVolume.summedCellVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters =
            region.currentVolumeCubicMeters
            - region.previousVolumeCubicMeters;
        region.previousOutwardRelativeFlowRateCubicMetersPerSecond =
            previousRegionFlux
                .outwardRelativeVolumeFlowRateCubicMetersPerSecond;
        region.currentOutwardRelativeFlowRateCubicMetersPerSecond =
            currentRegionFlux
                .outwardRelativeVolumeFlowRateCubicMetersPerSecond;
        region.integratedOutwardRelativeVolumeCubicMeters =
            0.5 * result.durationSeconds
            * (region.previousOutwardRelativeFlowRateCubicMetersPerSecond
               + region.currentOutwardRelativeFlowRateCubicMetersPerSecond);
        region.continuityResidualCubicMeters =
            region.geometryVolumeChangeCubicMeters
            + region.integratedOutwardRelativeVolumeCubicMeters;
        const double reference = std::max(
            std::abs(region.geometryVolumeChangeCubicMeters),
            std::abs(region.integratedOutwardRelativeVolumeCubicMeters));
        region.toleranceCubicMeters = std::max(
            settings.absoluteVolumeToleranceCubicMeters,
            settings.relativeVolumeTolerance * reference);
        region.withinTolerance =
            std::isfinite(region.continuityResidualCubicMeters)
            && std::abs(region.continuityResidualCubicMeters)
                <= region.toleranceCubicMeters;
        result.failedRegionCount += !region.withinTolerance;
        result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
            result.maximumAbsoluteContinuityResidualCubicMeters,
            std::abs(region.continuityResidualCubicMeters));
        result.globalGeometryVolumeChangeCubicMeters +=
            region.geometryVolumeChangeCubicMeters;
        result.globalIntegratedOutwardRelativeVolumeCubicMeters +=
            region.integratedOutwardRelativeVolumeCubicMeters;
        result.globalContinuityResidualCubicMeters +=
            region.continuityResidualCubicMeters;
        result.regions.push_back(region);
    }
    result.allRegionsWithinTolerance = result.failedRegionCount == 0;
    if (!std::isfinite(result.maximumAbsoluteContinuityResidualCubicMeters)
        || !std::isfinite(result.globalGeometryVolumeChangeCubicMeters)
        || !std::isfinite(
            result.globalIntegratedOutwardRelativeVolumeCubicMeters)
        || !std::isfinite(result.globalContinuityResidualCubicMeters)) {
        throw std::invalid_argument(
            "scene fluid region-continuity ledger is non-finite");
    }
    result.fingerprint = continuityFingerprint(result);
    return result;
}

} // namespace

SceneFluidRegionContinuitySet auditSceneFluidRegionContinuity(
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionContinuitySettings& settings,
    const SceneFluidRegionContinuityLimits& limits) {
    auto result = buildContinuity(
        previousVolumes, currentVolumes, previousFlux, currentFlux,
        settings, limits);
    validateSceneFluidRegionContinuity(
        result, previousVolumes, currentVolumes, previousFlux, currentFlux);
    return result;
}

void validateSceneFluidRegionContinuity(
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux) {
    validateSources(
        previousVolumes, currentVolumes, previousFlux, currentFlux);
    if (continuity.version != sceneFluidRegionContinuityVersion
        || continuity.fingerprint == 0
        || continuity.previousCellVolumeFingerprint
            != previousVolumes.fingerprint
        || continuity.currentCellVolumeFingerprint
            != currentVolumes.fingerprint
        || continuity.previousOpeningFluxFingerprint
            != previousFlux.fingerprint
        || continuity.currentOpeningFluxFingerprint
            != currentFlux.fingerprint) {
        throw std::invalid_argument(
            "scene fluid region-continuity identity is invalid");
    }
    const SceneFluidRegionContinuityLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildContinuity(
        previousVolumes, currentVolumes, previousFlux, currentFlux,
        continuity.settings, unlimited);
    if (continuity != expected
        || continuity.ownedStorageBytes
            != storageBytesForCount(continuity.regions.size())
        || continuity.fingerprint != continuityFingerprint(continuity)) {
        throw std::invalid_argument(
            "scene fluid region-continuity payload is invalid");
    }
}

} // namespace simwing::fsi
