#include "fluid/planar_region_flux.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    void integer(const std::uint64_t value) {
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value >> (8 * byte));
            value_ *= fnvPrime;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

bool validAxis(const GridFaceAxis axis) noexcept {
    return axis == GridFaceAxis::X
        || axis == GridFaceAxis::Y
        || axis == GridFaceAxis::Z;
}

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
    const double lowerVelocity,
    const double upperVelocity,
    const double fluidVelocity,
    const PlanarPressureRegionFluxSettings& settings) {
    const double reference = std::max({
        std::abs(lowerVelocity),
        std::abs(upperVelocity),
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

PlanarPressureRegionIntervalFlux deriveInterval(
    const std::uint64_t lowerSurfaceStableId,
    const std::uint64_t upperSurfaceStableId,
    const std::uint64_t regionStableId,
    const double previousVolumeCubicMeters,
    const double currentVolumeCubicMeters,
    const double lowerSurfaceVelocityMetersPerSecond,
    const double upperSurfaceVelocityMetersPerSecond,
    const double durationSeconds,
    const double crossSectionAreaSquareMeters,
    const PlanarPressureRegionFluxSettings& settings) {
    if (lowerSurfaceStableId == 0 || upperSurfaceStableId == 0
        || regionStableId == 0
        || !std::isfinite(previousVolumeCubicMeters)
        || !(previousVolumeCubicMeters > 0.0)
        || !std::isfinite(currentVolumeCubicMeters)
        || !(currentVolumeCubicMeters > 0.0)
        || !std::isfinite(lowerSurfaceVelocityMetersPerSecond)
        || !std::isfinite(upperSurfaceVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "planar pressure region flux primitive interval is invalid");
    }

    PlanarPressureRegionIntervalFlux interval;
    interval.lowerSurfaceStableId = lowerSurfaceStableId;
    interval.upperSurfaceStableId = upperSurfaceStableId;
    interval.regionStableId = regionStableId;
    interval.previousVolumeCubicMeters = previousVolumeCubicMeters;
    interval.currentVolumeCubicMeters = currentVolumeCubicMeters;
    interval.geometryVolumeChangeCubicMeters =
        currentVolumeCubicMeters - previousVolumeCubicMeters;
    interval.lowerSurfaceVelocityMetersPerSecond =
        lowerSurfaceVelocityMetersPerSecond;
    interval.upperSurfaceVelocityMetersPerSecond =
        upperSurfaceVelocityMetersPerSecond;
    interval.leastSquaresFluidVelocityMetersPerSecond =
        0.5 * lowerSurfaceVelocityMetersPerSecond
        + 0.5 * upperSurfaceVelocityMetersPerSecond;
    interval.lowerOutwardRelativeVelocityMetersPerSecond =
        lowerSurfaceVelocityMetersPerSecond
        - interval.leastSquaresFluidVelocityMetersPerSecond;
    interval.upperOutwardRelativeVelocityMetersPerSecond =
        interval.leastSquaresFluidVelocityMetersPerSecond
        - upperSurfaceVelocityMetersPerSecond;
    interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond =
        crossSectionAreaSquareMeters
        * interval.lowerOutwardRelativeVelocityMetersPerSecond;
    interval.upperOutwardRelativeFlowRateCubicMetersPerSecond =
        crossSectionAreaSquareMeters
        * interval.upperOutwardRelativeVelocityMetersPerSecond;
    interval.outwardRelativeFlowRateCubicMetersPerSecond =
        interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond
        + interval.upperOutwardRelativeFlowRateCubicMetersPerSecond;
    interval.integratedOutwardRelativeVolumeCubicMeters =
        durationSeconds
        * interval.outwardRelativeFlowRateCubicMetersPerSecond;
    interval.continuityResidualCubicMeters =
        interval.geometryVolumeChangeCubicMeters
        + interval.integratedOutwardRelativeVolumeCubicMeters;
    const double lowerSlip =
        interval.leastSquaresFluidVelocityMetersPerSecond
        - lowerSurfaceVelocityMetersPerSecond;
    const double upperSlip =
        interval.leastSquaresFluidVelocityMetersPerSecond
        - upperSurfaceVelocityMetersPerSecond;
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
        lowerSurfaceVelocityMetersPerSecond,
        upperSurfaceVelocityMetersPerSecond,
        interval.leastSquaresFluidVelocityMetersPerSecond,
        settings);
    interval.continuityToleranceCubicMeters = volumeTolerance(
        previousVolumeCubicMeters,
        currentVolumeCubicMeters,
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
    return interval;
}

std::uint64_t compatibilityFingerprint(
    const PlanarPressureRegionFluxCompatibility& compatibility) {
    Fingerprint fingerprint;
    fingerprint.integer(compatibility.version);
    fingerprint.integer(compatibility.sourceSweepVersion);
    fingerprint.integer(static_cast<std::uint64_t>(compatibility.axis));
    fingerprint.real(compatibility.durationSeconds);
    fingerprint.real(compatibility.crossSectionAreaSquareMeters);
    fingerprint.real(
        compatibility.settings.absoluteVelocityToleranceMetersPerSecond);
    fingerprint.real(compatibility.settings.relativeVelocityTolerance);
    fingerprint.real(
        compatibility.settings.absoluteVolumeToleranceCubicMeters);
    fingerprint.real(compatibility.settings.relativeVolumeTolerance);
    fingerprint.integer(compatibility.intervals.size());
    for (const auto& interval : compatibility.intervals) {
        fingerprint.integer(interval.lowerSurfaceStableId);
        fingerprint.integer(interval.upperSurfaceStableId);
        fingerprint.integer(interval.regionStableId);
        fingerprint.real(interval.previousVolumeCubicMeters);
        fingerprint.real(interval.currentVolumeCubicMeters);
        fingerprint.real(interval.geometryVolumeChangeCubicMeters);
        fingerprint.real(interval.lowerSurfaceVelocityMetersPerSecond);
        fingerprint.real(interval.upperSurfaceVelocityMetersPerSecond);
        fingerprint.real(
            interval.leastSquaresFluidVelocityMetersPerSecond);
        fingerprint.real(
            interval.lowerOutwardRelativeVelocityMetersPerSecond);
        fingerprint.real(
            interval.upperOutwardRelativeVelocityMetersPerSecond);
        fingerprint.real(
            interval.lowerOutwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            interval.upperOutwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            interval.outwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            interval.integratedOutwardRelativeVolumeCubicMeters);
        fingerprint.real(interval.continuityResidualCubicMeters);
        fingerprint.real(
            interval.maximumAbsoluteInterfaceSlipMetersPerSecond);
        fingerprint.real(interval.rmsInterfaceSlipMetersPerSecond);
        fingerprint.real(
            interval.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(interval.velocityToleranceMetersPerSecond);
        fingerprint.real(interval.continuityToleranceCubicMeters);
        fingerprint.integer(interval.impermeableWithinTolerance);
        fingerprint.integer(interval.continuityWithinTolerance);
    }
    fingerprint.integer(compatibility.regions.size());
    for (const auto& region : compatibility.regions) {
        fingerprint.integer(region.regionStableId);
        fingerprint.real(region.previousVolumeCubicMeters);
        fingerprint.real(region.currentVolumeCubicMeters);
        fingerprint.real(region.geometryVolumeChangeCubicMeters);
        fingerprint.real(region.outwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(region.integratedOutwardRelativeVolumeCubicMeters);
        fingerprint.real(region.continuityResidualCubicMeters);
        fingerprint.real(
            region.maximumAbsoluteInterfaceSlipMetersPerSecond);
        fingerprint.real(
            region.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(region.continuityToleranceCubicMeters);
        fingerprint.integer(region.impermeableWithinTolerance);
        fingerprint.integer(region.continuityWithinTolerance);
    }
    fingerprint.integer(compatibility.failedImpermeableIntervalCount);
    fingerprint.integer(compatibility.failedContinuityIntervalCount);
    fingerprint.integer(compatibility.failedImpermeableRegionCount);
    fingerprint.integer(compatibility.failedContinuityRegionCount);
    fingerprint.real(
        compatibility.maximumAbsoluteInterfaceSlipMetersPerSecond);
    fingerprint.real(
        compatibility.maximumAbsoluteContinuityResidualCubicMeters);
    fingerprint.real(compatibility.globalGeometryVolumeChangeCubicMeters);
    fingerprint.real(
        compatibility.globalIntegratedOutwardRelativeVolumeCubicMeters);
    fingerprint.real(compatibility.globalContinuityResidualCubicMeters);
    fingerprint.integer(
        compatibility.allIntervalsImpermeableWithinTolerance);
    fingerprint.integer(compatibility.allIntervalsContinuousWithinTolerance);
    fingerprint.integer(
        compatibility.allRegionsImpermeableWithinTolerance);
    fingerprint.integer(compatibility.allRegionsContinuousWithinTolerance);
    fingerprint.integer(compatibility.ownedStorageBytes);
    return fingerprint.value();
}

void finalizeCompatibility(
    PlanarPressureRegionFluxCompatibility& result,
    const PlanarPressureRegionFluxLimits& limits) {
    if (result.intervals.size() > limits.maximumIntervals) {
        throw std::length_error(
            "planar pressure region flux exceeds its interval limit");
    }
    const std::size_t maximumStorage = ownedStorageBytes(
        result.intervals.size(), result.intervals.size());
    if (maximumStorage > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region flux exceeds its byte limit");
    }

    result.regions.clear();
    result.failedImpermeableIntervalCount = 0;
    result.failedContinuityIntervalCount = 0;
    result.failedImpermeableRegionCount = 0;
    result.failedContinuityRegionCount = 0;
    result.maximumAbsoluteInterfaceSlipMetersPerSecond = 0.0;
    result.maximumAbsoluteContinuityResidualCubicMeters = 0.0;
    result.globalGeometryVolumeChangeCubicMeters = 0.0;
    result.globalIntegratedOutwardRelativeVolumeCubicMeters = 0.0;
    result.globalContinuityResidualCubicMeters = 0.0;
    result.allIntervalsImpermeableWithinTolerance = false;
    result.allIntervalsContinuousWithinTolerance = false;
    result.allRegionsImpermeableWithinTolerance = false;
    result.allRegionsContinuousWithinTolerance = false;
    result.ownedStorageBytes = 0;
    result.fingerprint = 0;

    std::map<std::uint64_t, PlanarPressureRegionFluxSummary> regions;
    for (const auto& interval : result.intervals) {
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
            result.settings);
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
    result.fingerprint = compatibilityFingerprint(result);
    if (result.fingerprint == 0) {
        throw std::invalid_argument(
            "planar pressure region flux fingerprint is invalid");
    }
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

    PlanarPressureRegionFluxCompatibility result;
    result.sourceSweepVersion = sweep.version;
    result.axis = sweep.axis;
    result.durationSeconds = sweep.durationSeconds;
    result.crossSectionAreaSquareMeters =
        sweep.crossSectionAreaSquareMeters;
    result.settings = settings;
    result.intervals.reserve(sweep.intervals.size());
    for (const auto& source : sweep.intervals) {
        result.intervals.push_back(deriveInterval(
            source.lowerSurfaceStableId,
            source.upperSurfaceStableId,
            source.regionStableId,
            source.previousVolumeCubicMeters,
            source.currentVolumeCubicMeters,
            source.lowerSurfaceVelocityMetersPerSecond,
            source.upperSurfaceVelocityMetersPerSecond,
            result.durationSeconds,
            result.crossSectionAreaSquareMeters,
            settings));
    }
    finalizeCompatibility(result, limits);
    validatePlanarPressureRegionFluxCompatibility(result, limits);
    return result;
}

void validatePlanarPressureRegionFluxCompatibility(
    const PlanarPressureRegionFluxCompatibility& compatibility,
    const PlanarPressureRegionFluxLimits& limits) {
    validateSettings(compatibility.settings);
    if (limits.maximumIntervals == 0 || limits.maximumRegions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region flux limits are invalid");
    }
    if (compatibility.version != planarPressureRegionFluxVersion
        || compatibility.fingerprint == 0
        || compatibility.sourceSweepVersion
            != planarPressureRegionSweepVersion
        || !validAxis(compatibility.axis)
        || !std::isfinite(compatibility.durationSeconds)
        || !(compatibility.durationSeconds > 0.0)
        || !std::isfinite(compatibility.crossSectionAreaSquareMeters)
        || !(compatibility.crossSectionAreaSquareMeters > 0.0)
        || compatibility.intervals.empty()
        || compatibility.intervals.size() > limits.maximumIntervals
        || compatibility.regions.empty()
        || compatibility.regions.size() > limits.maximumRegions) {
        throw std::invalid_argument(
            "planar pressure region flux compatibility metadata is invalid");
    }

    PlanarPressureRegionFluxCompatibility expected;
    expected.sourceSweepVersion = compatibility.sourceSweepVersion;
    expected.axis = compatibility.axis;
    expected.durationSeconds = compatibility.durationSeconds;
    expected.crossSectionAreaSquareMeters =
        compatibility.crossSectionAreaSquareMeters;
    expected.settings = compatibility.settings;
    expected.intervals.reserve(compatibility.intervals.size());
    for (std::size_t index = 0;
         index < compatibility.intervals.size(); ++index) {
        const auto& source = compatibility.intervals[index];
        const auto& next = compatibility.intervals[
            (index + 1) % compatibility.intervals.size()];
        if (source.upperSurfaceStableId
            != next.lowerSurfaceStableId) {
            throw std::invalid_argument(
                "planar pressure region flux interval chain is invalid");
        }
        expected.intervals.push_back(deriveInterval(
            source.lowerSurfaceStableId,
            source.upperSurfaceStableId,
            source.regionStableId,
            source.previousVolumeCubicMeters,
            source.currentVolumeCubicMeters,
            source.lowerSurfaceVelocityMetersPerSecond,
            source.upperSurfaceVelocityMetersPerSecond,
            compatibility.durationSeconds,
            compatibility.crossSectionAreaSquareMeters,
            compatibility.settings));
    }
    finalizeCompatibility(expected, limits);
    if (expected != compatibility) {
        throw std::invalid_argument(
            "planar pressure region flux compatibility ledger is invalid");
    }
}

} // namespace simwing::fsi::fluid
