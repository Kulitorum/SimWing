#include "fluid/planar_region_sweep.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

std::size_t checkedProduct(const std::size_t count,
                           const std::size_t itemBytes) {
    if (count != 0
        && itemBytes > std::numeric_limits<std::size_t>::max() / count) {
        throw std::length_error(
            "planar pressure region sweep storage size overflows");
    }
    return count * itemBytes;
}

void checkedAdd(std::size_t& total, const std::size_t bytes) {
    if (bytes > std::numeric_limits<std::size_t>::max() - total) {
        throw std::length_error(
            "planar pressure region sweep storage size overflows");
    }
    total += bytes;
}

std::size_t ownedStorageBytes(const std::size_t layerCount,
                              const std::size_t regionCount) {
    std::size_t result = 0;
    checkedAdd(result, checkedProduct(
        layerCount, 2 * sizeof(PlanarPressureRegionInterval)));
    checkedAdd(result, checkedProduct(
        regionCount, 2 * sizeof(PlanarPressureRegionSummary)));
    checkedAdd(result, checkedProduct(
        layerCount, sizeof(PlanarPressureRegionIntervalSweep)));
    checkedAdd(result, checkedProduct(
        regionCount, sizeof(PlanarPressureRegionSweepSummary)));
    return result;
}

double transverseArea(const PeriodicCartesianGrid& grid,
                      const GridFaceAxis axis) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    switch (axis) {
    case GridFaceAxis::X:
        return (upper.y - lower.y) * (upper.z - lower.z);
    case GridFaceAxis::Y:
        return (upper.x - lower.x) * (upper.z - lower.z);
    case GridFaceAxis::Z:
        return (upper.x - lower.x) * (upper.y - lower.y);
    }
    throw std::invalid_argument(
        "planar pressure region sweep axis is invalid");
}

using LayerBySurface = std::map<
    std::uint64_t, const PlanarPressureJumpLayerDefinition*>;

LayerBySurface indexLayers(
    const std::span<const PlanarPressureJumpLayerDefinition> layers) {
    LayerBySurface result;
    for (const auto& layer : layers) {
        if (!result.emplace(layer.surfaceStableId, &layer).second) {
            throw std::invalid_argument(
                "planar pressure region sweep surface identity is duplicated");
        }
    }
    return result;
}

const PlanarPressureJumpLayerDefinition& layer(
    const LayerBySurface& layers,
    const std::uint64_t surfaceStableId) {
    const auto found = layers.find(surfaceStableId);
    if (found == layers.end()) {
        throw std::invalid_argument(
            "planar pressure region sweep surface identity is missing");
    }
    return *found->second;
}

bool finiteInterval(const PlanarPressureRegionIntervalSweep& value) {
    return std::ranges::all_of(
        std::array{
            value.previousVolumeCubicMeters,
            value.currentVolumeCubicMeters,
            value.geometryVolumeChangeCubicMeters,
            value.lowerSurfaceDisplacementMeters,
            value.upperSurfaceDisplacementMeters,
            value.lowerSurfaceVelocityMetersPerSecond,
            value.upperSurfaceVelocityMetersPerSecond,
            value.boundarySweptVolumeCubicMeters,
            value.surfaceGeometryResidualCubicMeters,
        },
        [](const double sample) { return std::isfinite(sample); });
}

bool finiteSummary(const PlanarPressureRegionSweepSummary& value) {
    return std::ranges::all_of(
        std::array{
            value.previousVolumeCubicMeters,
            value.currentVolumeCubicMeters,
            value.geometryVolumeChangeCubicMeters,
            value.boundarySweptVolumeCubicMeters,
            value.surfaceGeometryResidualCubicMeters,
        },
        [](const double sample) { return std::isfinite(sample); });
}

bool validAxis(const GridFaceAxis axis) noexcept {
    return axis == GridFaceAxis::X
        || axis == GridFaceAxis::Y
        || axis == GridFaceAxis::Z;
}

void validateProfile(
    const StaticPlanarPressureRegionProfile& profile,
    const PlanarPressureRegionSweepLimits& limits) {
    if (!validAxis(profile.axis)
        || profile.intervals.empty()
        || profile.intervals.size() > limits.maximumLayers
        || profile.regions.empty()
        || profile.regions.size() > limits.maximumRegions
        || !std::isfinite(profile.windowLowerCoordinateMeters)
        || !std::isfinite(profile.windowUpperCoordinateMeters)
        || !(profile.windowUpperCoordinateMeters
             > profile.windowLowerCoordinateMeters)
        || !std::isfinite(profile.geometricDomainVolumeCubicMeters)
        || !(profile.geometricDomainVolumeCubicMeters > 0.0)
        || !std::isfinite(profile.intervalVolumeCubicMeters)
        || !(profile.intervalVolumeCubicMeters > 0.0)
        || !std::isfinite(profile.volumeClosureResidualCubicMeters)
        || !std::isfinite(profile.requestedVolumeMeanPressurePascals)
        || !std::isfinite(profile.achievedVolumeMeanPressurePascals)) {
        throw std::invalid_argument(
            "planar pressure region sweep profile metadata is invalid");
    }
    std::map<std::uint64_t, PlanarPressureRegionSummary> regions;
    double totalVolume = 0.0;
    double weightedPressure = 0.0;
    for (std::size_t index = 0; index < profile.intervals.size(); ++index) {
        const auto& interval = profile.intervals[index];
        const auto& next = profile.intervals[
            (index + 1) % profile.intervals.size()];
        if (interval.lowerSurfaceStableId == 0
            || interval.upperSurfaceStableId == 0
            || interval.regionStableId == 0
            || interval.upperSurfaceStableId
                != next.lowerSurfaceStableId
            || !std::isfinite(interval.lowerCoordinateMeters)
            || !std::isfinite(interval.upperCoordinateMeters)
            || !(interval.upperCoordinateMeters
                 > interval.lowerCoordinateMeters)
            || !std::isfinite(interval.volumeCubicMeters)
            || !(interval.volumeCubicMeters > 0.0)
            || !std::isfinite(interval.pressurePascals)) {
            throw std::invalid_argument(
                "planar pressure region sweep profile interval is invalid");
        }
        auto [found, inserted] = regions.emplace(
            interval.regionStableId,
            PlanarPressureRegionSummary{
                interval.regionStableId, 0.0,
                interval.pressurePascals});
        if (!inserted
            && found->second.pressurePascals
                != interval.pressurePascals) {
            throw std::invalid_argument(
                "planar pressure region sweep profile region pressure is inconsistent");
        }
        found->second.volumeCubicMeters += interval.volumeCubicMeters;
        totalVolume += interval.volumeCubicMeters;
        weightedPressure += interval.volumeCubicMeters
            * interval.pressurePascals;
        if (!std::isfinite(found->second.volumeCubicMeters)
            || !std::isfinite(totalVolume)
            || !std::isfinite(weightedPressure)) {
            throw std::invalid_argument(
                "planar pressure region sweep profile accumulation is non-finite");
        }
    }
    std::vector<PlanarPressureRegionSummary> expectedRegions;
    expectedRegions.reserve(regions.size());
    for (const auto& entry : regions) {
        expectedRegions.push_back(entry.second);
    }
    if (profile.windowLowerCoordinateMeters
            != profile.intervals.front().lowerCoordinateMeters
        || profile.windowUpperCoordinateMeters
            != profile.intervals.back().upperCoordinateMeters
        || profile.intervalVolumeCubicMeters != totalVolume
        || profile.volumeClosureResidualCubicMeters
            != totalVolume - profile.geometricDomainVolumeCubicMeters
        || profile.achievedVolumeMeanPressurePascals
            != weightedPressure / totalVolume
        || profile.regions != expectedRegions) {
        throw std::invalid_argument(
            "planar pressure region sweep profile ledger is invalid");
    }
}

} // namespace

PlanarPressureRegionSweepLedger makePlanarPressureRegionSweepLedger(
    const PeriodicCartesianGrid& grid,
    const std::span<const PlanarPressureJumpLayerDefinition> previousLayers,
    const std::span<const PlanarPressureJumpLayerDefinition> currentLayers,
    const double durationSeconds,
    const PlanarPressureRegionSweepLimits& limits) {
    if (!std::isfinite(durationSeconds) || !(durationSeconds > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region sweep duration must be finite and positive");
    }
    if (limits.maximumLayers == 0 || limits.maximumRegions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region sweep limits are invalid");
    }
    if (previousLayers.empty()
        || previousLayers.size() != currentLayers.size()) {
        throw std::invalid_argument(
            "planar pressure region sweep endpoint layer counts do not match");
    }
    if (previousLayers.size() > limits.maximumLayers) {
        throw std::length_error(
            "planar pressure region sweep exceeds its layer limit");
    }
    const std::size_t maximumStorage = ownedStorageBytes(
        previousLayers.size(), previousLayers.size());
    if (maximumStorage > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region sweep exceeds its byte limit");
    }

    PlanarPressureRegionSweepLedger result;
    result.durationSeconds = durationSeconds;
    result.previousProfile = makeStaticPlanarPressureRegionProfile(
        grid, previousLayers);
    result.currentProfile = makeStaticPlanarPressureRegionProfile(
        grid, currentLayers);
    if (result.previousProfile.axis != result.currentProfile.axis
        || result.previousProfile.intervals.size()
            != result.currentProfile.intervals.size()) {
        throw std::invalid_argument(
            "planar pressure region sweep endpoint topology does not match");
    }
    result.axis = result.previousProfile.axis;
    result.crossSectionAreaSquareMeters = transverseArea(grid, result.axis);
    if (!std::isfinite(result.crossSectionAreaSquareMeters)
        || !(result.crossSectionAreaSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region sweep cross-section area is invalid");
    }

    const LayerBySurface previousBySurface = indexLayers(previousLayers);
    const LayerBySurface currentBySurface = indexLayers(currentLayers);
    for (const auto& [surfaceStableId, previous] : previousBySurface) {
        const auto& current = layer(currentBySurface, surfaceStableId);
        if (previous->surfaceStableId != current.surfaceStableId
            || previous->minusRegionStableId
                != current.minusRegionStableId
            || previous->plusRegionStableId
                != current.plusRegionStableId
            || previous->pressureJumpPascals
                != current.pressureJumpPascals
            || previous->topology.axis != current.topology.axis) {
            throw std::invalid_argument(
                "planar pressure region sweep layer identity changed");
        }
        const auto selection = selectMovingPlanarTopology(
            grid, previous->topology,
            current.physicalPlaneCoordinateMeters);
        if (selection.topology != current.topology) {
            throw std::invalid_argument(
                "planar pressure region sweep endpoint topology is inconsistent");
        }
    }

    result.intervals.reserve(
        result.previousProfile.intervals.size());
    std::map<std::uint64_t, PlanarPressureRegionSweepSummary> regions;
    for (std::size_t index = 0;
         index < result.previousProfile.intervals.size(); ++index) {
        const auto& previous = result.previousProfile.intervals[index];
        const auto& current = result.currentProfile.intervals[index];
        if (previous.lowerSurfaceStableId
                != current.lowerSurfaceStableId
            || previous.upperSurfaceStableId
                != current.upperSurfaceStableId
            || previous.regionStableId != current.regionStableId) {
            throw std::invalid_argument(
                "planar pressure region sweep layer order changed");
        }
        const auto& previousLower = layer(
            previousBySurface, previous.lowerSurfaceStableId);
        const auto& currentLower = layer(
            currentBySurface, current.lowerSurfaceStableId);
        const auto& previousUpper = layer(
            previousBySurface, previous.upperSurfaceStableId);
        const auto& currentUpper = layer(
            currentBySurface, current.upperSurfaceStableId);
        PlanarPressureRegionIntervalSweep interval;
        interval.lowerSurfaceStableId = previous.lowerSurfaceStableId;
        interval.upperSurfaceStableId = previous.upperSurfaceStableId;
        interval.regionStableId = previous.regionStableId;
        interval.previousVolumeCubicMeters = previous.volumeCubicMeters;
        interval.currentVolumeCubicMeters = current.volumeCubicMeters;
        interval.geometryVolumeChangeCubicMeters =
            interval.currentVolumeCubicMeters
            - interval.previousVolumeCubicMeters;
        interval.lowerSurfaceDisplacementMeters =
            currentLower.physicalPlaneCoordinateMeters
            - previousLower.physicalPlaneCoordinateMeters;
        interval.upperSurfaceDisplacementMeters =
            currentUpper.physicalPlaneCoordinateMeters
            - previousUpper.physicalPlaneCoordinateMeters;
        interval.lowerSurfaceVelocityMetersPerSecond =
            interval.lowerSurfaceDisplacementMeters / durationSeconds;
        interval.upperSurfaceVelocityMetersPerSecond =
            interval.upperSurfaceDisplacementMeters / durationSeconds;
        interval.boundarySweptVolumeCubicMeters =
            result.crossSectionAreaSquareMeters
            * (interval.upperSurfaceDisplacementMeters
               - interval.lowerSurfaceDisplacementMeters);
        interval.surfaceGeometryResidualCubicMeters =
            interval.boundarySweptVolumeCubicMeters
            - interval.geometryVolumeChangeCubicMeters;
        if (!finiteInterval(interval)) {
            throw std::invalid_argument(
                "planar pressure region sweep interval ledger is non-finite");
        }
        result.maximumAbsoluteSurfaceGeometryResidualCubicMeters = std::max(
            result.maximumAbsoluteSurfaceGeometryResidualCubicMeters,
            std::abs(interval.surfaceGeometryResidualCubicMeters));
        auto& region = regions[interval.regionStableId];
        region.regionStableId = interval.regionStableId;
        region.previousVolumeCubicMeters +=
            interval.previousVolumeCubicMeters;
        region.currentVolumeCubicMeters +=
            interval.currentVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters +=
            interval.geometryVolumeChangeCubicMeters;
        region.boundarySweptVolumeCubicMeters +=
            interval.boundarySweptVolumeCubicMeters;
        region.surfaceGeometryResidualCubicMeters +=
            interval.surfaceGeometryResidualCubicMeters;
        result.intervals.push_back(interval);
    }

    if (regions.size() > limits.maximumRegions) {
        throw std::length_error(
            "planar pressure region sweep exceeds its region limit");
    }
    result.regions.reserve(regions.size());
    for (const auto& entry : regions) {
        const auto& region = entry.second;
        if (!finiteSummary(region)) {
            throw std::invalid_argument(
                "planar pressure region sweep region ledger is non-finite");
        }
        result.globalGeometryVolumeChangeCubicMeters +=
            region.geometryVolumeChangeCubicMeters;
        result.globalBoundarySweptVolumeCubicMeters +=
            region.boundarySweptVolumeCubicMeters;
        result.maximumAbsoluteSurfaceGeometryResidualCubicMeters = std::max(
            result.maximumAbsoluteSurfaceGeometryResidualCubicMeters,
            std::abs(region.surfaceGeometryResidualCubicMeters));
        result.regions.push_back(region);
    }
    result.globalSurfaceGeometryResidualCubicMeters =
        result.globalBoundarySweptVolumeCubicMeters
        - result.globalGeometryVolumeChangeCubicMeters;
    result.maximumAbsoluteSurfaceGeometryResidualCubicMeters = std::max(
        result.maximumAbsoluteSurfaceGeometryResidualCubicMeters,
        std::abs(result.globalSurfaceGeometryResidualCubicMeters));
    result.ownedStorageBytes = ownedStorageBytes(
        result.intervals.size(), result.regions.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region sweep exceeds its byte limit");
    }
    if (!std::isfinite(result.globalGeometryVolumeChangeCubicMeters)
        || !std::isfinite(result.globalBoundarySweptVolumeCubicMeters)
        || !std::isfinite(result.globalSurfaceGeometryResidualCubicMeters)) {
        throw std::invalid_argument(
            "planar pressure region sweep global ledger is non-finite");
    }
    validatePlanarPressureRegionSweepLedger(result, limits);
    return result;
}

void validatePlanarPressureRegionSweepLedger(
    const PlanarPressureRegionSweepLedger& ledger,
    const PlanarPressureRegionSweepLimits& limits) {
    if (limits.maximumLayers == 0 || limits.maximumRegions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region sweep limits are invalid");
    }
    if (ledger.version != planarPressureRegionSweepVersion
        || !validAxis(ledger.axis)
        || !std::isfinite(ledger.durationSeconds)
        || !(ledger.durationSeconds > 0.0)
        || !std::isfinite(ledger.crossSectionAreaSquareMeters)
        || !(ledger.crossSectionAreaSquareMeters > 0.0)
        || ledger.intervals.empty()
        || ledger.intervals.size() > limits.maximumLayers
        || ledger.regions.empty()
        || ledger.regions.size() > limits.maximumRegions) {
        throw std::invalid_argument(
            "planar pressure region sweep ledger metadata is invalid");
    }
    validateProfile(ledger.previousProfile, limits);
    validateProfile(ledger.currentProfile, limits);
    if (ledger.axis != ledger.previousProfile.axis
        || ledger.axis != ledger.currentProfile.axis
        || ledger.intervals.size()
            != ledger.previousProfile.intervals.size()
        || ledger.intervals.size()
            != ledger.currentProfile.intervals.size()) {
        throw std::invalid_argument(
            "planar pressure region sweep profile identity is invalid");
    }

    std::map<std::uint64_t, PlanarPressureRegionSweepSummary> regions;
    double maximumResidual = 0.0;
    for (std::size_t index = 0; index < ledger.intervals.size(); ++index) {
        const auto& interval = ledger.intervals[index];
        const auto& previous = ledger.previousProfile.intervals[index];
        const auto& current = ledger.currentProfile.intervals[index];
        const double geometryChange = current.volumeCubicMeters
            - previous.volumeCubicMeters;
        const double lowerVelocity =
            interval.lowerSurfaceDisplacementMeters
            / ledger.durationSeconds;
        const double upperVelocity =
            interval.upperSurfaceDisplacementMeters
            / ledger.durationSeconds;
        const double boundarySweep = ledger.crossSectionAreaSquareMeters
            * (interval.upperSurfaceDisplacementMeters
               - interval.lowerSurfaceDisplacementMeters);
        const double residual = boundarySweep - geometryChange;
        if (!finiteInterval(interval)
            || interval.lowerSurfaceStableId
                != previous.lowerSurfaceStableId
            || interval.upperSurfaceStableId
                != previous.upperSurfaceStableId
            || interval.regionStableId != previous.regionStableId
            || interval.lowerSurfaceStableId
                != current.lowerSurfaceStableId
            || interval.upperSurfaceStableId
                != current.upperSurfaceStableId
            || interval.regionStableId != current.regionStableId
            || interval.previousVolumeCubicMeters
                != previous.volumeCubicMeters
            || interval.currentVolumeCubicMeters
                != current.volumeCubicMeters
            || interval.geometryVolumeChangeCubicMeters != geometryChange
            || interval.lowerSurfaceVelocityMetersPerSecond
                != lowerVelocity
            || interval.upperSurfaceVelocityMetersPerSecond
                != upperVelocity
            || interval.boundarySweptVolumeCubicMeters != boundarySweep
            || interval.surfaceGeometryResidualCubicMeters != residual) {
            throw std::invalid_argument(
                "planar pressure region sweep interval ledger is invalid");
        }
        auto& region = regions[interval.regionStableId];
        region.regionStableId = interval.regionStableId;
        region.previousVolumeCubicMeters +=
            interval.previousVolumeCubicMeters;
        region.currentVolumeCubicMeters +=
            interval.currentVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters += geometryChange;
        region.boundarySweptVolumeCubicMeters += boundarySweep;
        region.surfaceGeometryResidualCubicMeters += residual;
        maximumResidual = std::max(maximumResidual, std::abs(residual));
    }

    std::vector<PlanarPressureRegionSweepSummary> expectedRegions;
    expectedRegions.reserve(regions.size());
    double globalGeometryChange = 0.0;
    double globalBoundarySweep = 0.0;
    for (const auto& entry : regions) {
        const auto& region = entry.second;
        if (!finiteSummary(region)) {
            throw std::invalid_argument(
                "planar pressure region sweep region ledger is non-finite");
        }
        expectedRegions.push_back(region);
        globalGeometryChange += region.geometryVolumeChangeCubicMeters;
        globalBoundarySweep += region.boundarySweptVolumeCubicMeters;
        maximumResidual = std::max(
            maximumResidual,
            std::abs(region.surfaceGeometryResidualCubicMeters));
    }
    const double globalResidual = globalBoundarySweep
        - globalGeometryChange;
    maximumResidual = std::max(maximumResidual, std::abs(globalResidual));
    const std::size_t storageBytes = ownedStorageBytes(
        ledger.intervals.size(), ledger.regions.size());
    if (storageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region sweep exceeds its byte limit");
    }
    if (!std::isfinite(globalGeometryChange)
        || !std::isfinite(globalBoundarySweep)
        || !std::isfinite(globalResidual)
        || ledger.regions != expectedRegions
        || ledger.maximumAbsoluteSurfaceGeometryResidualCubicMeters
            != maximumResidual
        || ledger.globalGeometryVolumeChangeCubicMeters
            != globalGeometryChange
        || ledger.globalBoundarySweptVolumeCubicMeters
            != globalBoundarySweep
        || ledger.globalSurfaceGeometryResidualCubicMeters
            != globalResidual
        || ledger.ownedStorageBytes != storageBytes) {
        throw std::invalid_argument(
            "planar pressure region sweep aggregate ledger is invalid");
    }
}

} // namespace simwing::fsi::fluid
