#include "fluid/planar_region_opening_power.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
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
    std::uint64_t value_ = fnvOffsetBasis;
};

struct EndpointPressure {
    double previousPascals = 0.0;
    double currentPascals = 0.0;
};

void validateSettings(
    const PlanarPressureRegionOpeningPowerSettings& settings) {
    if (!std::isfinite(settings.absolutePowerToleranceWatts)
        || settings.absolutePowerToleranceWatts < 0.0
        || !std::isfinite(settings.relativePowerTolerance)
        || settings.relativePowerTolerance < 0.0
        || (settings.absolutePowerToleranceWatts == 0.0
            && settings.relativePowerTolerance == 0.0)) {
        throw std::invalid_argument(
            "planar pressure region opening-power settings are invalid");
    }
}

void validateLimits(
    const PlanarPressureRegionOpeningPowerLimits& limits) {
    if (limits.maximumOpenings == 0 || limits.maximumRegions == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar pressure region opening-power limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t count,
                           const std::size_t itemBytes) {
    if (count != 0
        && itemBytes > std::numeric_limits<std::size_t>::max() / count) {
        throw std::length_error(
            "planar pressure region opening-power storage size overflows");
    }
    return count * itemBytes;
}

std::size_t ownedStorageBytes(const std::size_t openingCount,
                              const std::size_t regionCount) {
    const std::size_t openings = checkedProduct(
        openingCount, sizeof(PlanarPressureRegionOpeningPower));
    const std::size_t regions = checkedProduct(
        regionCount, sizeof(PlanarPressureRegionPressureVolumePower));
    if (regions > std::numeric_limits<std::size_t>::max() - openings) {
        throw std::length_error(
            "planar pressure region opening-power storage size overflows");
    }
    return openings + regions;
}

double powerTolerance(
    const double reference,
    const PlanarPressureRegionOpeningPowerSettings& settings) {
    const double tolerance = std::max(
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance * std::abs(reference));
    if (!std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "planar pressure region opening-power tolerance is non-finite");
    }
    return tolerance;
}

std::map<std::uint64_t, EndpointPressure> endpointPressures(
    const PlanarPressureRegionSweepLedger& sweep) {
    std::map<std::uint64_t, EndpointPressure> pressures;
    for (const auto& region : sweep.previousProfile.regions) {
        pressures.emplace(
            region.regionStableId,
            EndpointPressure{region.pressurePascals, 0.0});
    }
    if (pressures.size() != sweep.previousProfile.regions.size()) {
        throw std::invalid_argument(
            "planar pressure region opening-power pressure identity is invalid");
    }
    for (const auto& region : sweep.currentProfile.regions) {
        const auto found = pressures.find(region.regionStableId);
        if (found == pressures.end()) {
            throw std::invalid_argument(
                "planar pressure region opening-power pressure identity changed");
        }
        found->second.currentPascals = region.pressurePascals;
    }
    if (pressures.size() != sweep.currentProfile.regions.size()) {
        throw std::invalid_argument(
            "planar pressure region opening-power pressure identity changed");
    }
    return pressures;
}

std::uint64_t auditFingerprint(
    const PlanarPressureRegionOpeningPowerAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.sourceSweepVersion);
    fingerprint.integer(audit.sourceOpeningFlowVersion);
    fingerprint.integer(audit.sourceOpeningFlowFingerprint);
    fingerprint.enumeration(audit.axis);
    fingerprint.real(audit.durationSeconds);
    fingerprint.real(audit.settings.absolutePowerToleranceWatts);
    fingerprint.real(audit.settings.relativePowerTolerance);
    fingerprint.integer(static_cast<std::uint64_t>(audit.openings.size()));
    for (const auto& opening : audit.openings) {
        fingerprint.integer(opening.openingStableId);
        fingerprint.integer(opening.negativeSideRegionStableId);
        fingerprint.integer(opening.positiveSideRegionStableId);
        fingerprint.real(opening.midpointNegativePressurePascals);
        fingerprint.real(opening.midpointPositivePressurePascals);
        fingerprint.real(opening.negativeToPositivePressureDropPascals);
        fingerprint.real(
            opening.relativeVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(opening.pressurePowerWatts);
        fingerprint.real(opening.minimumExternalPowerWatts);
        fingerprint.real(opening.passiveToleranceWatts);
        fingerprint.integer(static_cast<std::uint8_t>(
            opening.passiveWithinTolerance));
    }
    fingerprint.integer(static_cast<std::uint64_t>(audit.regions.size()));
    for (const auto& region : audit.regions) {
        fingerprint.integer(region.regionStableId);
        fingerprint.real(region.previousPressurePascals);
        fingerprint.real(region.currentPressurePascals);
        fingerprint.real(region.midpointPressurePascals);
        fingerprint.real(region.geometryVolumeRateCubicMetersPerSecond);
        fingerprint.real(region.pressureVolumePowerWatts);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(audit.failedPassiveOpeningCount));
    fingerprint.real(audit.maximumOpeningExternalPowerWatts);
    fingerprint.real(audit.summedOpeningExternalPowerDeficitWatts);
    fingerprint.real(audit.totalOpeningPressurePowerWatts);
    fingerprint.real(audit.minimumNetExternalPowerWatts);
    fingerprint.real(audit.totalRegionPressureVolumePowerWatts);
    fingerprint.real(audit.pressurePowerClosureResidualWatts);
    fingerprint.real(audit.pressurePowerClosureToleranceWatts);
    fingerprint.integer(static_cast<std::uint8_t>(
        audit.allOpeningsPassiveWithinTolerance));
    fingerprint.integer(static_cast<std::uint8_t>(
        audit.pressurePowerClosesWithinTolerance));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionOpeningPowerAudit buildAudit(
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition> definitions,
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionOpeningPowerSettings& settings,
    const PlanarPressureRegionOpeningPowerLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    validatePlanarPressureRegionOpeningFlow(
        allocation, sweep, definitions, limits.sourceLimits);
    if (!allocation.allComponentsFeasible
        || !allocation.allRegionsWithinTolerance) {
        throw std::invalid_argument(
            "planar pressure region opening-power source is not kinematically feasible");
    }
    if (allocation.openings.size() > limits.maximumOpenings
        || allocation.regions.size() > limits.maximumRegions) {
        throw std::length_error(
            "planar pressure region opening-power count limit exceeded");
    }
    const std::size_t storageBytes = ownedStorageBytes(
        allocation.openings.size(), allocation.regions.size());
    if (storageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region opening-power byte limit exceeded");
    }
    const auto pressures = endpointPressures(sweep);

    PlanarPressureRegionOpeningPowerAudit result;
    result.sourceSweepVersion = sweep.version;
    result.sourceOpeningFlowVersion = allocation.version;
    result.sourceOpeningFlowFingerprint = allocation.fingerprint;
    result.axis = sweep.axis;
    result.durationSeconds = sweep.durationSeconds;
    result.settings = settings;
    result.ownedStorageBytes = storageBytes;
    result.regions.reserve(allocation.regions.size());
    double regionPowerMagnitudeWatts = 0.0;
    for (const auto& source : allocation.regions) {
        const auto pressure = pressures.find(source.regionStableId);
        if (pressure == pressures.end()) {
            throw std::invalid_argument(
                "planar pressure region opening-power region is foreign");
        }
        PlanarPressureRegionPressureVolumePower region;
        region.regionStableId = source.regionStableId;
        region.previousPressurePascals =
            pressure->second.previousPascals;
        region.currentPressurePascals =
            pressure->second.currentPascals;
        region.midpointPressurePascals =
            0.5 * region.previousPressurePascals
            + 0.5 * region.currentPressurePascals;
        region.geometryVolumeRateCubicMetersPerSecond =
            source.geometryVolumeChangeCubicMeters
            / result.durationSeconds;
        region.pressureVolumePowerWatts =
            region.midpointPressurePascals
            * region.geometryVolumeRateCubicMetersPerSecond;
        if (!std::isfinite(region.midpointPressurePascals)
            || !std::isfinite(
                region.geometryVolumeRateCubicMetersPerSecond)
            || !std::isfinite(region.pressureVolumePowerWatts)) {
            throw std::invalid_argument(
                "planar pressure region opening-power region is non-finite");
        }
        result.totalRegionPressureVolumePowerWatts +=
            region.pressureVolumePowerWatts;
        regionPowerMagnitudeWatts +=
            std::abs(region.pressureVolumePowerWatts);
        result.regions.push_back(region);
    }

    result.openings.reserve(allocation.openings.size());
    double openingPowerMagnitudeWatts = 0.0;
    for (const auto& source : allocation.openings) {
        const auto negative = pressures.find(
            source.negativeSideRegionStableId);
        const auto positive = pressures.find(
            source.positiveSideRegionStableId);
        if (negative == pressures.end() || positive == pressures.end()) {
            throw std::invalid_argument(
                "planar pressure region opening-power opening is foreign");
        }
        PlanarPressureRegionOpeningPower opening;
        opening.openingStableId = source.openingStableId;
        opening.negativeSideRegionStableId =
            source.negativeSideRegionStableId;
        opening.positiveSideRegionStableId =
            source.positiveSideRegionStableId;
        opening.midpointNegativePressurePascals = 0.5
            * (negative->second.previousPascals
               + negative->second.currentPascals);
        opening.midpointPositivePressurePascals = 0.5
            * (positive->second.previousPascals
               + positive->second.currentPascals);
        opening.negativeToPositivePressureDropPascals =
            opening.midpointNegativePressurePascals
            - opening.midpointPositivePressurePascals;
        opening.relativeVolumeFlowRateCubicMetersPerSecond =
            source.relativeVolumeFlowRateCubicMetersPerSecond;
        opening.pressurePowerWatts =
            opening.negativeToPositivePressureDropPascals
            * opening.relativeVolumeFlowRateCubicMetersPerSecond;
        opening.minimumExternalPowerWatts = std::max(
            0.0, -opening.pressurePowerWatts);
        opening.passiveToleranceWatts = powerTolerance(
            opening.pressurePowerWatts, settings);
        opening.passiveWithinTolerance = opening.pressurePowerWatts
            >= -opening.passiveToleranceWatts;
        if (!std::isfinite(opening.pressurePowerWatts)
            || !std::isfinite(opening.minimumExternalPowerWatts)) {
            throw std::invalid_argument(
                "planar pressure region opening power is non-finite");
        }
        result.failedPassiveOpeningCount +=
            !opening.passiveWithinTolerance;
        result.maximumOpeningExternalPowerWatts = std::max(
            result.maximumOpeningExternalPowerWatts,
            opening.minimumExternalPowerWatts);
        result.summedOpeningExternalPowerDeficitWatts +=
            opening.minimumExternalPowerWatts;
        result.totalOpeningPressurePowerWatts +=
            opening.pressurePowerWatts;
        openingPowerMagnitudeWatts +=
            std::abs(opening.pressurePowerWatts);
        result.openings.push_back(opening);
    }

    result.minimumNetExternalPowerWatts = std::max(
        0.0, -result.totalOpeningPressurePowerWatts);
    result.pressurePowerClosureResidualWatts =
        result.totalOpeningPressurePowerWatts
        + result.totalRegionPressureVolumePowerWatts;
    const double closureReference = std::max(
        openingPowerMagnitudeWatts, regionPowerMagnitudeWatts);
    result.pressurePowerClosureToleranceWatts = powerTolerance(
        closureReference, settings);
    result.allOpeningsPassiveWithinTolerance =
        result.failedPassiveOpeningCount == 0;
    result.pressurePowerClosesWithinTolerance =
        std::isfinite(result.pressurePowerClosureResidualWatts)
        && std::abs(result.pressurePowerClosureResidualWatts)
            <= result.pressurePowerClosureToleranceWatts;
    if (!std::isfinite(result.maximumOpeningExternalPowerWatts)
        || !std::isfinite(result.summedOpeningExternalPowerDeficitWatts)
        || !std::isfinite(result.totalOpeningPressurePowerWatts)
        || !std::isfinite(result.minimumNetExternalPowerWatts)
        || !std::isfinite(result.totalRegionPressureVolumePowerWatts)) {
        throw std::invalid_argument(
            "planar pressure region opening-power aggregate is non-finite");
    }
    result.fingerprint = auditFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionOpeningPowerAudit
auditPlanarPressureRegionOpeningPower(
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition>
        openingDefinitions,
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionOpeningPowerSettings& settings,
    const PlanarPressureRegionOpeningPowerLimits& limits) {
    auto result = buildAudit(
        sweep, openingDefinitions, allocation, settings, limits);
    validatePlanarPressureRegionOpeningPower(
        result, sweep, openingDefinitions, allocation, limits);
    return result;
}

void validatePlanarPressureRegionOpeningPower(
    const PlanarPressureRegionOpeningPowerAudit& audit,
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition>
        openingDefinitions,
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionOpeningPowerLimits& limits) {
    validateLimits(limits);
    if (audit.openings.size() > limits.maximumOpenings
        || audit.regions.size() > limits.maximumRegions) {
        throw std::length_error(
            "planar pressure region opening-power validation limit exceeded");
    }
    const auto expected = buildAudit(
        sweep, openingDefinitions, allocation, audit.settings, limits);
    if (expected != audit) {
        throw std::invalid_argument(
            "planar pressure region opening-power audit is invalid");
    }
}

} // namespace simwing::fsi::fluid
