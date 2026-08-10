#include "fluid/planar_region_flux.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

void validateSettings(const PlanarPressureRegionFluxSettings& settings) {
    if (!std::isfinite(settings.absoluteVelocityToleranceMetersPerSecond)
        || settings.absoluteVelocityToleranceMetersPerSecond < 0.0
        || !std::isfinite(settings.relativeVelocityTolerance)
        || settings.relativeVelocityTolerance < 0.0
        || (settings.absoluteVelocityToleranceMetersPerSecond == 0.0
            && settings.relativeVelocityTolerance == 0.0)
        || !std::isfinite(settings.absoluteVolumeToleranceCubicMeters)
        || settings.absoluteVolumeToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeVolumeTolerance)
        || settings.relativeVolumeTolerance < 0.0
        || (settings.absoluteVolumeToleranceCubicMeters == 0.0
            && settings.relativeVolumeTolerance == 0.0)) {
        throw std::invalid_argument(
            "planar pressure region flux tolerances are invalid");
    }
}

std::size_t checkedProduct(const std::size_t count,
                           const std::size_t itemBytes) {
    if (count != 0
        && itemBytes > std::numeric_limits<std::size_t>::max() / count) {
        throw std::length_error(
            "planar pressure region flux storage size overflows");
    }
    return count * itemBytes;
}

std::size_t ownedStorageBytes(const std::size_t intervalCount,
                              const std::size_t regionCount) {
    const std::size_t intervals = checkedProduct(
        intervalCount, sizeof(PlanarPressureRegionIntervalFlux));
    const std::size_t regions = checkedProduct(
        regionCount, sizeof(PlanarPressureRegionFluxSummary));
    if (regions > std::numeric_limits<std::size_t>::max() - intervals) {
        throw std::length_error(
            "planar pressure region flux storage size overflows");
    }
    return intervals + regions;
}

double velocityTolerance(
    const PlanarPressureRegionIntervalSweep& source,
    const double fluidVelocity,
    const PlanarPressureRegionFluxSettings& settings) {
    const double reference = std::max({
        std::abs(source.lowerSurfaceVelocityMetersPerSecond),
        std::abs(source.upperSurfaceVelocityMetersPerSecond),
        std::abs(fluidVelocity),
    });
    const double tolerance = std::max(
        settings.absoluteVelocityToleranceMetersPerSecond,
        settings.relativeVelocityTolerance * reference);
    if (!std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "planar pressure region flux velocity tolerance is non-finite");
    }
    return tolerance;
}

double volumeTolerance(
    const double previousVolume,
    const double currentVolume,
    const double geometryChange,
    const double integratedRelativeVolume,
    const PlanarPressureRegionFluxSettings& settings) {
    const double reference = std::max({
        std::abs(previousVolume),
        std::abs(currentVolume),
        std::abs(geometryChange),
        std::abs(integratedRelativeVolume),
    });
    const double tolerance = std::max(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance * reference);
    if (!std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "planar pressure region flux volume tolerance is non-finite");
    }
    return tolerance;
}

bool finiteInterval(const PlanarPressureRegionIntervalFlux& value) {
    return std::ranges::all_of(
        std::array{
            value.previousVolumeCubicMeters,
            value.currentVolumeCubicMeters,
            value.geometryVolumeChangeCubicMeters,
            value.lowerSurfaceVelocityMetersPerSecond,
            value.upperSurfaceVelocityMetersPerSecond,
            value.leastSquaresFluidVelocityMetersPerSecond,
            value.lowerOutwardRelativeVelocityMetersPerSecond,
            value.upperOutwardRelativeVelocityMetersPerSecond,
            value.lowerOutwardRelativeFlowRateCubicMetersPerSecond,
            value.upperOutwardRelativeFlowRateCubicMetersPerSecond,
            value.outwardRelativeFlowRateCubicMetersPerSecond,
            value.integratedOutwardRelativeVolumeCubicMeters,
            value.continuityResidualCubicMeters,
            value.maximumAbsoluteInterfaceSlipMetersPerSecond,
            value.rmsInterfaceSlipMetersPerSecond,
            value.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond,
            value.velocityToleranceMetersPerSecond,
            value.continuityToleranceCubicMeters,
        },
        [](const double sample) { return std::isfinite(sample); });
}

bool finiteSummary(const PlanarPressureRegionFluxSummary& value) {
    return std::ranges::all_of(
        std::array{
            value.previousVolumeCubicMeters,
            value.currentVolumeCubicMeters,
            value.geometryVolumeChangeCubicMeters,
            value.outwardRelativeFlowRateCubicMetersPerSecond,
            value.integratedOutwardRelativeVolumeCubicMeters,
            value.continuityResidualCubicMeters,
            value.maximumAbsoluteInterfaceSlipMetersPerSecond,
            value.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond,
            value.continuityToleranceCubicMeters,
        },
        [](const double sample) { return std::isfinite(sample); });
}

} // namespace

PlanarPressureRegionFluxCompatibility
assessPlanarPressureRegionFluxCompatibility(
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFluxSettings& settings,
    const PlanarPressureRegionFluxLimits& limits) {
    validateSettings(settings);
    if (limits.maximumIntervals == 0 || limits.maximumRegions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region flux limits are invalid");
    }
    validatePlanarPressureRegionSweepLedger(
        sweep,
        {
            limits.maximumIntervals,
            limits.maximumRegions,
            std::numeric_limits<std::size_t>::max(),
        });
    if (sweep.intervals.size() > limits.maximumIntervals) {
        throw std::length_error(
            "planar pressure region flux exceeds its interval limit");
    }
    const std::size_t maximumStorage = ownedStorageBytes(
        sweep.intervals.size(), sweep.intervals.size());
    if (maximumStorage > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region flux exceeds its byte limit");
    }

    PlanarPressureRegionFluxCompatibility result;
    result.sourceSweepVersion = sweep.version;
    result.axis = sweep.axis;
    result.durationSeconds = sweep.durationSeconds;
    result.crossSectionAreaSquareMeters =
        sweep.crossSectionAreaSquareMeters;
    result.settings = settings;
    result.intervals.reserve(sweep.intervals.size());
    std::map<std::uint64_t, PlanarPressureRegionFluxSummary> regions;
    for (const auto& source : sweep.intervals) {
        PlanarPressureRegionIntervalFlux interval;
        interval.lowerSurfaceStableId = source.lowerSurfaceStableId;
        interval.upperSurfaceStableId = source.upperSurfaceStableId;
        interval.regionStableId = source.regionStableId;
        interval.previousVolumeCubicMeters =
            source.previousVolumeCubicMeters;
        interval.currentVolumeCubicMeters =
            source.currentVolumeCubicMeters;
        interval.geometryVolumeChangeCubicMeters =
            source.geometryVolumeChangeCubicMeters;
        interval.lowerSurfaceVelocityMetersPerSecond =
            source.lowerSurfaceVelocityMetersPerSecond;
        interval.upperSurfaceVelocityMetersPerSecond =
            source.upperSurfaceVelocityMetersPerSecond;
        interval.leastSquaresFluidVelocityMetersPerSecond =
            0.5 * interval.lowerSurfaceVelocityMetersPerSecond
            + 0.5 * interval.upperSurfaceVelocityMetersPerSecond;
        interval.lowerOutwardRelativeVelocityMetersPerSecond =
            interval.lowerSurfaceVelocityMetersPerSecond
            - interval.leastSquaresFluidVelocityMetersPerSecond;
        interval.upperOutwardRelativeVelocityMetersPerSecond =
            interval.leastSquaresFluidVelocityMetersPerSecond
            - interval.upperSurfaceVelocityMetersPerSecond;
        interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond =
            result.crossSectionAreaSquareMeters
            * interval.lowerOutwardRelativeVelocityMetersPerSecond;
        interval.upperOutwardRelativeFlowRateCubicMetersPerSecond =
            result.crossSectionAreaSquareMeters
            * interval.upperOutwardRelativeVelocityMetersPerSecond;
        interval.outwardRelativeFlowRateCubicMetersPerSecond =
            interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond
            + interval.upperOutwardRelativeFlowRateCubicMetersPerSecond;
        interval.integratedOutwardRelativeVolumeCubicMeters =
            result.durationSeconds
            * interval.outwardRelativeFlowRateCubicMetersPerSecond;
        interval.continuityResidualCubicMeters =
            interval.geometryVolumeChangeCubicMeters
            + interval.integratedOutwardRelativeVolumeCubicMeters;
        const double lowerSlip =
            interval.leastSquaresFluidVelocityMetersPerSecond
            - interval.lowerSurfaceVelocityMetersPerSecond;
        const double upperSlip =
            interval.leastSquaresFluidVelocityMetersPerSecond
            - interval.upperSurfaceVelocityMetersPerSecond;
        interval.maximumAbsoluteInterfaceSlipMetersPerSecond = std::max(
            std::abs(lowerSlip), std::abs(upperSlip));
        interval.rmsInterfaceSlipMetersPerSecond =
            std::hypot(lowerSlip, upperSlip) / std::sqrt(2.0);
        interval.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond =
            std::abs(
                interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond)
            + std::abs(
                interval.upperOutwardRelativeFlowRateCubicMetersPerSecond);
        interval.velocityToleranceMetersPerSecond = velocityTolerance(
            source, interval.leastSquaresFluidVelocityMetersPerSecond,
            settings);
        interval.continuityToleranceCubicMeters = volumeTolerance(
            interval.previousVolumeCubicMeters,
            interval.currentVolumeCubicMeters,
            interval.geometryVolumeChangeCubicMeters,
            interval.integratedOutwardRelativeVolumeCubicMeters,
            settings);
        interval.impermeableWithinTolerance =
            interval.maximumAbsoluteInterfaceSlipMetersPerSecond
            <= interval.velocityToleranceMetersPerSecond;
        interval.continuityWithinTolerance =
            std::abs(interval.continuityResidualCubicMeters)
            <= interval.continuityToleranceCubicMeters;
        if (!finiteInterval(interval)) {
            throw std::invalid_argument(
                "planar pressure region flux interval is non-finite");
        }
        result.failedImpermeableIntervalCount +=
            !interval.impermeableWithinTolerance;
        result.failedContinuityIntervalCount +=
            !interval.continuityWithinTolerance;
        result.maximumAbsoluteInterfaceSlipMetersPerSecond = std::max(
            result.maximumAbsoluteInterfaceSlipMetersPerSecond,
            interval.maximumAbsoluteInterfaceSlipMetersPerSecond);
        result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
            result.maximumAbsoluteContinuityResidualCubicMeters,
            std::abs(interval.continuityResidualCubicMeters));

        auto [regionEntry, inserted] = regions.try_emplace(
            interval.regionStableId);
        auto& region = regionEntry->second;
        if (inserted) {
            region.regionStableId = interval.regionStableId;
            region.impermeableWithinTolerance = true;
            region.continuityWithinTolerance = true;
        }
        region.previousVolumeCubicMeters +=
            interval.previousVolumeCubicMeters;
        region.currentVolumeCubicMeters +=
            interval.currentVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters +=
            interval.geometryVolumeChangeCubicMeters;
        region.outwardRelativeFlowRateCubicMetersPerSecond +=
            interval.outwardRelativeFlowRateCubicMetersPerSecond;
        region.integratedOutwardRelativeVolumeCubicMeters +=
            interval.integratedOutwardRelativeVolumeCubicMeters;
        region.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond +=
            interval.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond;
        region.maximumAbsoluteInterfaceSlipMetersPerSecond = std::max(
            region.maximumAbsoluteInterfaceSlipMetersPerSecond,
            interval.maximumAbsoluteInterfaceSlipMetersPerSecond);
        region.impermeableWithinTolerance =
            region.impermeableWithinTolerance
            && interval.impermeableWithinTolerance;
        result.intervals.push_back(interval);
    }

    if (regions.size() > limits.maximumRegions) {
        throw std::length_error(
            "planar pressure region flux exceeds its region limit");
    }
    result.regions.reserve(regions.size());
    for (auto& entry : regions) {
        auto& region = entry.second;
        region.continuityResidualCubicMeters =
            region.geometryVolumeChangeCubicMeters
            + region.integratedOutwardRelativeVolumeCubicMeters;
        region.continuityToleranceCubicMeters = volumeTolerance(
            region.previousVolumeCubicMeters,
            region.currentVolumeCubicMeters,
            region.geometryVolumeChangeCubicMeters,
            region.integratedOutwardRelativeVolumeCubicMeters,
            settings);
        region.continuityWithinTolerance =
            std::abs(region.continuityResidualCubicMeters)
            <= region.continuityToleranceCubicMeters;
        if (!finiteSummary(region)) {
            throw std::invalid_argument(
                "planar pressure region flux summary is non-finite");
        }
        result.failedImpermeableRegionCount +=
            !region.impermeableWithinTolerance;
        result.failedContinuityRegionCount +=
            !region.continuityWithinTolerance;
        result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
            result.maximumAbsoluteContinuityResidualCubicMeters,
            std::abs(region.continuityResidualCubicMeters));
        result.globalGeometryVolumeChangeCubicMeters +=
            region.geometryVolumeChangeCubicMeters;
        result.globalIntegratedOutwardRelativeVolumeCubicMeters +=
            region.integratedOutwardRelativeVolumeCubicMeters;
        result.regions.push_back(region);
    }
    result.globalContinuityResidualCubicMeters =
        result.globalGeometryVolumeChangeCubicMeters
        + result.globalIntegratedOutwardRelativeVolumeCubicMeters;
    result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
        result.maximumAbsoluteContinuityResidualCubicMeters,
        std::abs(result.globalContinuityResidualCubicMeters));
    result.allIntervalsImpermeableWithinTolerance =
        result.failedImpermeableIntervalCount == 0;
    result.allIntervalsContinuousWithinTolerance =
        result.failedContinuityIntervalCount == 0;
    result.allRegionsImpermeableWithinTolerance =
        result.failedImpermeableRegionCount == 0;
    result.allRegionsContinuousWithinTolerance =
        result.failedContinuityRegionCount == 0;
    result.ownedStorageBytes = ownedStorageBytes(
        result.intervals.size(), result.regions.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region flux exceeds its byte limit");
    }
    if (!std::isfinite(result.maximumAbsoluteInterfaceSlipMetersPerSecond)
        || !std::isfinite(
            result.maximumAbsoluteContinuityResidualCubicMeters)
        || !std::isfinite(result.globalGeometryVolumeChangeCubicMeters)
        || !std::isfinite(
            result.globalIntegratedOutwardRelativeVolumeCubicMeters)
        || !std::isfinite(result.globalContinuityResidualCubicMeters)) {
        throw std::invalid_argument(
            "planar pressure region flux aggregate is non-finite");
    }
    return result;
}

} // namespace simwing::fsi::fluid
