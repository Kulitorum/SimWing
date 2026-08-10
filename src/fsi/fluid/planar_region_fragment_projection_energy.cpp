#include "fluid/planar_region_fragment_projection_energy.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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

void validateTolerancePair(const double absolute,
                           const double relative,
                           const char* message) {
    if (!std::isfinite(absolute) || absolute < 0.0
        || !std::isfinite(relative) || relative < 0.0
        || (absolute == 0.0 && relative == 0.0)) {
        throw std::invalid_argument(message);
    }
}

void validateSettings(
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "planar regional projection-energy physical settings are "
            "invalid");
    }
    validateTolerancePair(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance,
        "planar regional projection-energy continuity tolerance is invalid");
    validateTolerancePair(
        settings.absoluteVelocityResidualToleranceMetersPerSecond,
        settings.relativeVelocityResidualTolerance,
        "planar regional projection-energy velocity tolerance is invalid");
    validateTolerancePair(
        settings
            .absoluteMomentumResidualToleranceKilogramMetersPerSecond,
        settings.relativeMomentumResidualTolerance,
        "planar regional projection-energy momentum tolerance is invalid");
    validateTolerancePair(
        settings.absoluteEnergyResidualToleranceJoules,
        settings.relativeEnergyResidualTolerance,
        "planar regional projection-energy energy tolerance is invalid");
    validateTolerancePair(
        settings.absolutePressureGaugeTolerancePascals,
        settings.relativePressureGaugeTolerance,
        "planar regional projection-energy gauge tolerance is invalid");
}

void validateLimits(
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    if (limits.maximumCorrections == 0
        || limits.maximumPressureSamples == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional projection-energy limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional projection-energy storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional projection-energy storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t pressureCount,
                              const std::size_t correctionCount,
                              const std::size_t componentCount) {
    return checkedSum(
        checkedSum(
            checkedProduct(pressureCount, sizeof(double)),
            checkedProduct(
                correctionCount,
                sizeof(
                    PlanarPressureRegionFragmentProjectionEnergyCorrection))),
        checkedProduct(
            componentCount,
            sizeof(
                PlanarPressureRegionFragmentProjectionEnergyComponent)));
}

bool isStaticGeometry(const PlanarPressureRegionSweepLedger& sweep) {
    if (sweep.previousProfile != sweep.currentProfile
        || sweep.globalGeometryVolumeChangeCubicMeters != 0.0
        || sweep.globalBoundarySweptVolumeCubicMeters != 0.0
        || sweep.globalSurfaceGeometryResidualCubicMeters != 0.0) {
        return false;
    }
    return std::ranges::all_of(
        sweep.intervals,
        [](const PlanarPressureRegionIntervalSweep& interval) {
            return interval.previousVolumeCubicMeters
                    == interval.currentVolumeCubicMeters
                && interval.geometryVolumeChangeCubicMeters == 0.0
                && interval.lowerSurfaceDisplacementMeters == 0.0
                && interval.upperSurfaceDisplacementMeters == 0.0
                && interval.lowerSurfaceVelocityMetersPerSecond == 0.0
                && interval.upperSurfaceVelocityMetersPerSecond == 0.0
                && interval.boundarySweptVolumeCubicMeters == 0.0
                && interval.surfaceGeometryResidualCubicMeters == 0.0;
        });
}

double& vectorCoordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "planar regional projection-energy axis is invalid");
}

Vector3 vectorDifference(const Vector3& after, const Vector3& before) {
    return {
        after.x - before.x,
        after.y - before.y,
        after.z - before.z,
    };
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double scaledTolerance(const double absolute,
                       const double relative,
                       const std::initializer_list<double> scales) {
    double scale = 0.0;
    for (const double value : scales) {
        scale = std::max(scale, std::abs(value));
    }
    return std::max(absolute, relative * scale);
}

double rootMeanSquare(const std::vector<double>& values) {
    double squaredSum = 0.0;
    for (const double value : values) squaredSum += value * value;
    return std::sqrt(squaredSum / static_cast<double>(values.size()));
}

double maximumAbsolute(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings) {
    for (const double value : {
             settings.densityKgPerCubicMeter,
             settings.timeStepSeconds,
             settings.absoluteContinuityToleranceCubicMetersPerSecond,
             settings.relativeContinuityTolerance,
             settings.absoluteVelocityResidualToleranceMetersPerSecond,
             settings.relativeVelocityResidualTolerance,
             settings
                 .absoluteMomentumResidualToleranceKilogramMetersPerSecond,
             settings.relativeMomentumResidualTolerance,
             settings.absoluteEnergyResidualToleranceJoules,
             settings.relativeEnergyResidualTolerance,
             settings.absolutePressureGaugeTolerancePascals,
             settings.relativePressureGaugeTolerance}) {
        fingerprint.real(value);
    }
}

std::uint64_t auditFingerprint(
    const PlanarPressureRegionFragmentProjectionEnergyAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.sourceMetricFingerprint);
    fingerprint.integer(audit.sourceTopologyFingerprint);
    fingerprint.integer(audit.sourceFragmentFingerprint);
    fingerprint.integer(audit.volumeRateFingerprint);
    fingerprint.integer(audit.beforeVelocityStateFingerprint);
    fingerprint.integer(audit.afterVelocityStateFingerprint);
    fingerprint.boolean(audit.staticGeometry);
    fingerprint.boolean(audit.usesMovingVolumeRates);
    fingerprintSettings(fingerprint, audit.settings);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.pressureCorrectionPascals.size()));
    for (const double pressure : audit.pressureCorrectionPascals) {
        fingerprint.real(pressure);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.corrections.size()));
    for (const auto& correction : audit.corrections) {
        fingerprint.integer(static_cast<std::uint64_t>(
            correction.correctionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            correction.dofIndex));
        fingerprint.integer(correction.dofStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            correction.sourceFaceLinkIndex));
        fingerprint.integer(correction.sourceFaceLinkStableId);
        fingerprint.enumeration(correction.axis);
        fingerprint.integer(static_cast<std::uint64_t>(
            correction.componentIndex));
        fingerprint.integer(correction.regionStableId);
        for (const double value : {
                 correction.areaSquareMeters,
                 correction.centerDistanceMeters,
                 correction.dualVolumeCubicMeters,
                 correction.diagonalMassKilograms,
                 correction.pressureDifferencePascals,
                 correction.velocityBeforeMetersPerSecond,
                 correction.velocityAfterMetersPerSecond,
                 correction.expectedVelocityChangeMetersPerSecond,
                 correction.velocityChangeResidualMetersPerSecond,
                 correction.momentumChangeKilogramMetersPerSecond,
                 correction.pressureImpulseKilogramMetersPerSecond,
                 correction
                     .momentumImpulseResidualKilogramMetersPerSecond,
                 correction.kineticEnergyChangeJoules,
                 correction.midpointPressureWorkJoules,
                 correction.workEnergyResidualJoules,
                 correction.correctionKineticEnergyJoules,
                 correction.finalPressureWorkJoules,
                 correction.affineWorkResidualJoules}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.components.size()));
    for (const auto& component : audit.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.correctionCount));
        fingerprint.real(component.pressureCorrectionVolumeMeanPascals);
        fingerprint.real(
            component.predictedContinuityResidualCubicMetersPerSecond);
        fingerprint.real(
            component.correctedContinuityResidualCubicMetersPerSecond);
        fingerprint.real(
            component.geometryVolumeRateCubicMetersPerSecond);
        fingerprintVector(
            fingerprint, component.momentumChangeKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint, component.pressureImpulseKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint,
            component
                .momentumImpulseResidualKilogramMetersPerSecond);
        fingerprint.real(component.kineticEnergyChangeJoules);
        fingerprint.real(component.midpointPressureWorkJoules);
        fingerprint.real(component.workEnergyResidualJoules);
        fingerprint.real(component.correctionKineticEnergyJoules);
        fingerprint.real(component.finalPressureWorkJoules);
        fingerprint.real(component.geometryPressureWorkJoules);
        fingerprint.real(component.finalGeometryWorkResidualJoules);
        fingerprint.real(component.affineEnergyResidualJoules);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.pressureLayerTraceCount));
    for (const double value : {
             audit.predictedContinuityResidualL2CubicMetersPerSecond,
             audit.predictedContinuityResidualMaximumCubicMetersPerSecond,
             audit.correctedContinuityResidualL2CubicMetersPerSecond,
             audit.correctedContinuityResidualMaximumCubicMetersPerSecond,
             audit.continuityToleranceCubicMetersPerSecond,
             audit.maximumAbsoluteVelocityChangeResidualMetersPerSecond,
             audit
                 .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
             audit.maximumAbsolutePressureGaugePascals,
             audit.maximumAbsoluteWallTraceVelocityResidualMetersPerSecond,
             audit
                 .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprintVector(
        fingerprint, audit.momentumBeforeKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, audit.momentumAfterKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, audit.momentumChangeKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, audit.pressureImpulseKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint,
        audit.momentumImpulseResidualKilogramMetersPerSecond);
    for (const double value : {
             audit.kineticEnergyBeforeJoules,
             audit.kineticEnergyAfterJoules,
             audit.kineticEnergyChangeJoules,
             audit.midpointPressureWorkJoules,
             audit.workEnergyResidualJoules,
             audit.correctionKineticEnergyJoules,
             audit.finalPressureWorkJoules,
             audit.geometryPressureWorkJoules,
             audit.finalGeometryWorkResidualJoules,
             audit.affineEnergyResidualJoules,
             audit.kineticEnergyRemovedJoules}) {
        fingerprint.real(value);
    }
    fingerprint.boolean(audit.nonIncreasingKineticEnergy);
    fingerprint.boolean(audit.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentProjectionEnergyAudit buildAudit(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet* volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    const bool staticGeometry = isStaticGeometry(sweep);
    if (volumeRates == nullptr && !staticGeometry) {
        throw std::invalid_argument(
            "planar regional projection-energy audit requires static "
            "geometry");
    }
    if (volumeRates != nullptr) {
        validatePlanarPressureRegionFragmentVolumeRates(
            *volumeRates, grid, sweep, fragments, topology,
            limits.volumeRateLimits);
        if (volumeRates->durationSeconds != settings.timeStepSeconds) {
            throw std::invalid_argument(
                "planar regional projection-energy duration does not "
                "match its volume rates");
        }
    }
    validatePlanarPressureRegionFragmentVelocityState(
        before, grid, sweep, fragments, topology, metric,
        limits.velocityStateLimits);
    validatePlanarPressureRegionFragmentVelocityState(
        after, grid, sweep, fragments, topology, metric,
        limits.velocityStateLimits);
    if (before.densityKgPerCubicMeter != settings.densityKgPerCubicMeter
        || after.densityKgPerCubicMeter
            != settings.densityKgPerCubicMeter
        || before.samples.size() != after.samples.size()
        || pressureCorrectionPascals.size() != fragments.fragments.size()
        || pressureCorrectionPascals.size() > limits.maximumPressureSamples
        || metric.sharedRegionGridDofCount > limits.maximumCorrections
        || metric.components.size() > limits.maximumComponents) {
        throw std::invalid_argument(
            "planar regional projection-energy sources are incompatible");
    }
    if (!std::ranges::all_of(
            pressureCorrectionPascals,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "planar regional projection-energy pressure is not finite");
    }

    PlanarPressureRegionFragmentProjectionEnergyAudit result;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.volumeRateFingerprint =
        volumeRates == nullptr ? 0 : volumeRates->fingerprint;
    result.beforeVelocityStateFingerprint = before.fingerprint;
    result.afterVelocityStateFingerprint = after.fingerprint;
    result.staticGeometry = staticGeometry;
    result.usesMovingVolumeRates = volumeRates != nullptr;
    result.settings = settings;
    result.ownedStorageBytes = ownedStorageBytes(
        pressureCorrectionPascals.size(), metric.sharedRegionGridDofCount,
        metric.components.size());
    result.workingStorageBytes = checkedSum(
        checkedProduct(
            checkedProduct(fragments.fragments.size(), sizeof(double)), 2),
        checkedProduct(metric.components.size(), sizeof(double)));
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional projection-energy storage limit exceeded");
    }
    result.pressureCorrectionPascals = pressureCorrectionPascals;
    result.corrections.reserve(metric.sharedRegionGridDofCount);
    result.components.resize(metric.components.size());
    for (std::size_t index = 0; index < metric.components.size(); ++index) {
        const auto& source = metric.components[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.regionStableId = source.regionStableId;
    }

    double maximumPressure = 0.0;
    std::vector<double> pressureVolumeMoment(metric.components.size(), 0.0);
    for (std::size_t index = 0;
         index < pressureCorrectionPascals.size(); ++index) {
        const double pressure = pressureCorrectionPascals[index];
        maximumPressure = std::max(maximumPressure, std::abs(pressure));
        pressureVolumeMoment[metric.fragments[index].componentIndex] +=
            pressure * metric.fragments[index].sourceVolumeCubicMeters;
    }
    const double gaugeTolerance = scaledTolerance(
        settings.absolutePressureGaugeTolerancePascals,
        settings.relativePressureGaugeTolerance, {maximumPressure});
    for (std::size_t index = 0; index < result.components.size(); ++index) {
        auto& component = result.components[index];
        component.pressureCorrectionVolumeMeanPascals =
            pressureVolumeMoment[index]
            / metric.components[index].sourceVolumeCubicMeters;
        result.maximumAbsolutePressureGaugePascals = std::max(
            result.maximumAbsolutePressureGaugePascals,
            std::abs(component.pressureCorrectionVolumeMeanPascals));
        if (!std::isfinite(component.pressureCorrectionVolumeMeanPascals)
            || std::abs(component.pressureCorrectionVolumeMeanPascals)
                > gaugeTolerance) {
            throw std::invalid_argument(
                "planar regional projection-energy correction gauge is "
                "invalid");
        }
    }

    std::vector<double> predictedFlow(fragments.fragments.size(), 0.0);
    std::vector<double> correctedFlow(fragments.fragments.size(), 0.0);
    for (std::size_t index = 0; index < metric.dofs.size(); ++index) {
        const auto& dof = metric.dofs[index];
        const auto& beforeSample = before.samples[index];
        const auto& afterSample = after.samples[index];
        if (dof.kind
            != PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            ++result.pressureLayerTraceCount;
            double expectedWallVelocity = 0.0;
            if (volumeRates != nullptr) {
                const auto& rate = volumeRates->fragments[
                    dof.ownerFragmentIndex];
                expectedWallVelocity = dof.kind
                        == PlanarPressureRegionFragmentVelocityDofKind::
                            PressureLayerMinusTrace
                    ? rate.upperBoundaryVelocityMetersPerSecond
                    : rate.lowerBoundaryVelocityMetersPerSecond;
            }
            const double beforeResidual =
                beforeSample.normalVelocityMetersPerSecond
                - expectedWallVelocity;
            const double afterResidual =
                afterSample.normalVelocityMetersPerSecond
                - expectedWallVelocity;
            result.maximumAbsoluteWallTraceVelocityResidualMetersPerSecond =
                std::max({
                    result
                        .maximumAbsoluteWallTraceVelocityResidualMetersPerSecond,
                    std::abs(beforeResidual), std::abs(afterResidual)});
            if (!std::isfinite(expectedWallVelocity)
                || beforeResidual != 0.0 || afterResidual != 0.0) {
                throw std::invalid_argument(
                    "planar regional projection-energy wall trace is "
                    "inconsistent with material motion");
            }
            continue;
        }

        const auto& link = topology.links[dof.sourceFaceLinkIndex];
        const double pressureDifference =
            pressureCorrectionPascals[link.minusFragmentIndex]
            - pressureCorrectionPascals[link.plusFragmentIndex];
        const double expectedVelocityChange = settings.timeStepSeconds
            / settings.densityKgPerCubicMeter * pressureDifference
            / link.centerDistanceMeters;
        const double actualVelocityChange =
            afterSample.normalVelocityMetersPerSecond
            - beforeSample.normalVelocityMetersPerSecond;
        const double velocityResidual =
            actualVelocityChange - expectedVelocityChange;
        const double velocityTolerance = scaledTolerance(
            settings.absoluteVelocityResidualToleranceMetersPerSecond,
            settings.relativeVelocityResidualTolerance,
            {beforeSample.normalVelocityMetersPerSecond,
             afterSample.normalVelocityMetersPerSecond,
             expectedVelocityChange});
        const double momentumChange =
            afterSample.normalMomentumKilogramMetersPerSecond
            - beforeSample.normalMomentumKilogramMetersPerSecond;
        const double pressureImpulse = settings.timeStepSeconds
            * link.areaSquareMeters * pressureDifference;
        const double momentumResidual = momentumChange - pressureImpulse;
        const double momentumTolerance = scaledTolerance(
            settings
                .absoluteMomentumResidualToleranceKilogramMetersPerSecond,
            settings.relativeMomentumResidualTolerance,
            {momentumChange, pressureImpulse});
        const double kineticEnergyChange =
            afterSample.kineticEnergyJoules
            - beforeSample.kineticEnergyJoules;
        const double midpointPressureWork = pressureImpulse * 0.5
            * (beforeSample.normalVelocityMetersPerSecond
               + afterSample.normalVelocityMetersPerSecond);
        const double workEnergyResidual =
            kineticEnergyChange - midpointPressureWork;
        const double correctionKineticEnergy = 0.5
            * settings.densityKgPerCubicMeter
            * dof.dualVolumeCubicMeters
            * actualVelocityChange * actualVelocityChange;
        const double finalPressureWork = pressureImpulse
            * afterSample.normalVelocityMetersPerSecond;
        const double affineWorkResidual = midpointPressureWork
            - (finalPressureWork - correctionKineticEnergy);
        const double energyTolerance = scaledTolerance(
            settings.absoluteEnergyResidualToleranceJoules,
            settings.relativeEnergyResidualTolerance,
            {beforeSample.kineticEnergyJoules,
             afterSample.kineticEnergyJoules,
             kineticEnergyChange,
             midpointPressureWork,
             correctionKineticEnergy,
             finalPressureWork});
        if (!std::isfinite(expectedVelocityChange)
            || !std::isfinite(velocityResidual)
            || !std::isfinite(momentumResidual)
            || !std::isfinite(workEnergyResidual)
            || !std::isfinite(correctionKineticEnergy)
            || !std::isfinite(finalPressureWork)
            || !std::isfinite(affineWorkResidual)
            || std::abs(velocityResidual) > velocityTolerance
            || std::abs(momentumResidual) > momentumTolerance
            || std::abs(workEnergyResidual) > energyTolerance
            || std::abs(affineWorkResidual) > energyTolerance) {
            throw std::invalid_argument(
                "planar regional projection-energy correction closure "
                "failed");
        }

        result.corrections.push_back({
            result.corrections.size(),
            dof.dofIndex,
            dof.stableId,
            dof.sourceFaceLinkIndex,
            dof.sourceFaceLinkStableId,
            dof.axis,
            dof.componentIndex,
            dof.regionStableId,
            dof.areaSquareMeters,
            link.centerDistanceMeters,
            dof.dualVolumeCubicMeters,
            settings.densityKgPerCubicMeter
                * dof.dualVolumeCubicMeters,
            pressureDifference,
            beforeSample.normalVelocityMetersPerSecond,
            afterSample.normalVelocityMetersPerSecond,
            expectedVelocityChange,
            velocityResidual,
            momentumChange,
            pressureImpulse,
            momentumResidual,
            kineticEnergyChange,
            midpointPressureWork,
            workEnergyResidual,
            correctionKineticEnergy,
            finalPressureWork,
            affineWorkResidual,
        });
        result.maximumAbsoluteVelocityChangeResidualMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteVelocityChangeResidualMetersPerSecond,
                std::abs(velocityResidual));
        result
            .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
                std::abs(momentumResidual));

        auto& component = result.components[dof.componentIndex];
        ++component.correctionCount;
        vectorCoordinate(
            component.pressureImpulseKilogramMetersPerSecond,
            dof.axis) += pressureImpulse;
        component.midpointPressureWorkJoules += midpointPressureWork;
        component.correctionKineticEnergyJoules +=
            correctionKineticEnergy;
        component.finalPressureWorkJoules += finalPressureWork;
        vectorCoordinate(
            result.pressureImpulseKilogramMetersPerSecond,
            dof.axis) += pressureImpulse;
        result.midpointPressureWorkJoules += midpointPressureWork;
        result.correctionKineticEnergyJoules += correctionKineticEnergy;
        result.finalPressureWorkJoules += finalPressureWork;

        const double predictedVolumeFlow = link.areaSquareMeters
            * beforeSample.normalVelocityMetersPerSecond;
        const double correctedVolumeFlow = link.areaSquareMeters
            * afterSample.normalVelocityMetersPerSecond;
        predictedFlow[link.minusFragmentIndex] += predictedVolumeFlow;
        predictedFlow[link.plusFragmentIndex] -= predictedVolumeFlow;
        correctedFlow[link.minusFragmentIndex] += correctedVolumeFlow;
        correctedFlow[link.plusFragmentIndex] -= correctedVolumeFlow;
    }
    if (result.corrections.size() != metric.sharedRegionGridDofCount
        || result.pressureLayerTraceCount
            != metric.pressureLayerTraceDofCount) {
        throw std::invalid_argument(
            "planar regional projection-energy DOF coverage is incomplete");
    }

    for (std::size_t index = 0; index < fragments.fragments.size(); ++index) {
        const double geometryVolumeRate = volumeRates == nullptr
            ? 0.0
            : volumeRates->fragments[index]
                .geometryVolumeChangeRateCubicMetersPerSecond;
        result.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
                std::abs(geometryVolumeRate));
        predictedFlow[index] += geometryVolumeRate;
        correctedFlow[index] += geometryVolumeRate;
        const std::size_t componentIndex =
            metric.fragments[index].componentIndex;
        result.components[componentIndex]
            .geometryVolumeRateCubicMetersPerSecond += geometryVolumeRate;
        const double geometryPressureWork = -settings.timeStepSeconds
            * pressureCorrectionPascals[index] * geometryVolumeRate;
        result.components[componentIndex].geometryPressureWorkJoules +=
            geometryPressureWork;
        result.geometryPressureWorkJoules += geometryPressureWork;
    }

    result.predictedContinuityResidualL2CubicMetersPerSecond =
        rootMeanSquare(predictedFlow);
    result.predictedContinuityResidualMaximumCubicMetersPerSecond =
        maximumAbsolute(predictedFlow);
    result.correctedContinuityResidualL2CubicMetersPerSecond =
        rootMeanSquare(correctedFlow);
    result.correctedContinuityResidualMaximumCubicMetersPerSecond =
        maximumAbsolute(correctedFlow);
    result.continuityToleranceCubicMetersPerSecond = std::max(
        settings.absoluteContinuityToleranceCubicMetersPerSecond,
        settings.relativeContinuityTolerance
            * result
                .predictedContinuityResidualMaximumCubicMetersPerSecond);
    if (!std::isfinite(
            result.predictedContinuityResidualL2CubicMetersPerSecond)
        || !std::isfinite(
            result.correctedContinuityResidualL2CubicMetersPerSecond)
        || result.correctedContinuityResidualMaximumCubicMetersPerSecond
            > result.continuityToleranceCubicMetersPerSecond) {
        throw std::invalid_argument(
            "planar regional projection-energy continuity closure failed");
    }

    result.momentumBeforeKilogramMetersPerSecond =
        before.momentumKilogramMetersPerSecond;
    result.momentumAfterKilogramMetersPerSecond =
        after.momentumKilogramMetersPerSecond;
    result.momentumChangeKilogramMetersPerSecond = vectorDifference(
        result.momentumAfterKilogramMetersPerSecond,
        result.momentumBeforeKilogramMetersPerSecond);
    result.momentumImpulseResidualKilogramMetersPerSecond = vectorDifference(
        result.momentumChangeKilogramMetersPerSecond,
        result.pressureImpulseKilogramMetersPerSecond);
    result.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        std::max(
            result
                .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
            maximumAbsoluteComponent(
                result.momentumImpulseResidualKilogramMetersPerSecond));
    const double globalMomentumTolerance = scaledTolerance(
        settings.absoluteMomentumResidualToleranceKilogramMetersPerSecond,
        settings.relativeMomentumResidualTolerance,
        {maximumAbsoluteComponent(
             result.momentumChangeKilogramMetersPerSecond),
         maximumAbsoluteComponent(
             result.pressureImpulseKilogramMetersPerSecond)});
    if (!finiteVector(result.momentumImpulseResidualKilogramMetersPerSecond)
        || maximumAbsoluteComponent(
               result.momentumImpulseResidualKilogramMetersPerSecond)
            > globalMomentumTolerance) {
        throw std::invalid_argument(
            "planar regional projection-energy global momentum closure "
            "failed");
    }

    result.kineticEnergyBeforeJoules = before.kineticEnergyJoules;
    result.kineticEnergyAfterJoules = after.kineticEnergyJoules;
    result.kineticEnergyChangeJoules =
        result.kineticEnergyAfterJoules
        - result.kineticEnergyBeforeJoules;
    result.workEnergyResidualJoules = result.kineticEnergyChangeJoules
        - result.midpointPressureWorkJoules;
    result.finalGeometryWorkResidualJoules =
        result.finalPressureWorkJoules
        - result.geometryPressureWorkJoules;
    result.affineEnergyResidualJoules = result.kineticEnergyChangeJoules
        - (result.geometryPressureWorkJoules
           - result.correctionKineticEnergyJoules);
    result.kineticEnergyRemovedJoules =
        -result.kineticEnergyChangeJoules;
    const double globalEnergyTolerance = scaledTolerance(
        settings.absoluteEnergyResidualToleranceJoules,
        settings.relativeEnergyResidualTolerance,
        {result.kineticEnergyBeforeJoules,
         result.kineticEnergyAfterJoules,
         result.kineticEnergyChangeJoules,
         result.midpointPressureWorkJoules,
         result.correctionKineticEnergyJoules,
         result.finalPressureWorkJoules,
         result.geometryPressureWorkJoules});
    result.nonIncreasingKineticEnergy =
        result.kineticEnergyAfterJoules
        <= result.kineticEnergyBeforeJoules + globalEnergyTolerance;
    if (!std::isfinite(result.workEnergyResidualJoules)
        || !std::isfinite(result.finalGeometryWorkResidualJoules)
        || !std::isfinite(result.affineEnergyResidualJoules)
        || std::abs(result.workEnergyResidualJoules)
            > globalEnergyTolerance
        || std::abs(result.finalGeometryWorkResidualJoules)
            > globalEnergyTolerance
        || std::abs(result.affineEnergyResidualJoules)
            > globalEnergyTolerance
        || (volumeRates == nullptr
            && !result.nonIncreasingKineticEnergy)) {
        throw std::invalid_argument(
            "planar regional projection-energy global energy closure "
            "failed");
    }

    for (std::size_t fragmentIndex = 0;
         fragmentIndex < fragments.fragments.size(); ++fragmentIndex) {
        auto& component = result.components[
            metric.fragments[fragmentIndex].componentIndex];
        component.predictedContinuityResidualCubicMetersPerSecond +=
            predictedFlow[fragmentIndex];
        component.correctedContinuityResidualCubicMetersPerSecond +=
            correctedFlow[fragmentIndex];
    }
    for (std::size_t index = 0; index < result.components.size(); ++index) {
        auto& component = result.components[index];
        component.momentumChangeKilogramMetersPerSecond = vectorDifference(
            after.components[index].momentumKilogramMetersPerSecond,
            before.components[index].momentumKilogramMetersPerSecond);
        component.momentumImpulseResidualKilogramMetersPerSecond =
            vectorDifference(
                component.momentumChangeKilogramMetersPerSecond,
                component.pressureImpulseKilogramMetersPerSecond);
        result
            .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
                maximumAbsoluteComponent(
                    component
                        .momentumImpulseResidualKilogramMetersPerSecond));
        component.kineticEnergyChangeJoules =
            after.components[index].kineticEnergyJoules
            - before.components[index].kineticEnergyJoules;
        component.workEnergyResidualJoules =
            component.kineticEnergyChangeJoules
            - component.midpointPressureWorkJoules;
        component.finalGeometryWorkResidualJoules =
            component.finalPressureWorkJoules
            - component.geometryPressureWorkJoules;
        component.affineEnergyResidualJoules =
            component.kineticEnergyChangeJoules
            - (component.geometryPressureWorkJoules
               - component.correctionKineticEnergyJoules);
        const double componentMomentumTolerance = scaledTolerance(
            settings
                .absoluteMomentumResidualToleranceKilogramMetersPerSecond,
            settings.relativeMomentumResidualTolerance,
            {maximumAbsoluteComponent(
                 component.momentumChangeKilogramMetersPerSecond),
             maximumAbsoluteComponent(
                 component.pressureImpulseKilogramMetersPerSecond)});
        const double componentEnergyTolerance = scaledTolerance(
            settings.absoluteEnergyResidualToleranceJoules,
            settings.relativeEnergyResidualTolerance,
            {before.components[index].kineticEnergyJoules,
             after.components[index].kineticEnergyJoules,
             component.kineticEnergyChangeJoules,
             component.midpointPressureWorkJoules,
             component.correctionKineticEnergyJoules,
             component.finalPressureWorkJoules,
             component.geometryPressureWorkJoules});
        if (!finiteVector(
                component
                    .momentumImpulseResidualKilogramMetersPerSecond)
            || maximumAbsoluteComponent(
                   component
                       .momentumImpulseResidualKilogramMetersPerSecond)
                > componentMomentumTolerance
            || !std::isfinite(component.workEnergyResidualJoules)
            || !std::isfinite(component.finalGeometryWorkResidualJoules)
            || !std::isfinite(component.affineEnergyResidualJoules)
            || std::abs(component.workEnergyResidualJoules)
                > componentEnergyTolerance
            || std::abs(component.finalGeometryWorkResidualJoules)
                > componentEnergyTolerance
            || std::abs(component.affineEnergyResidualJoules)
                > componentEnergyTolerance
            || std::abs(
                   component
                       .correctedContinuityResidualCubicMetersPerSecond)
                > result.continuityToleranceCubicMetersPerSecond) {
            throw std::invalid_argument(
                "planar regional projection-energy component closure "
                "failed");
        }
    }

    result.accepted = true;
    result.fingerprint = auditFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentProjectionEnergyAudit
auditStaticPlanarPressureRegionFragmentProjectionEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    return buildAudit(
        grid, sweep, fragments, topology, nullptr, metric, before, after,
        pressureCorrectionPascals, settings, limits);
}

void validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
    const PlanarPressureRegionFragmentProjectionEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    validateLimits(limits);
    if (audit.pressureCorrectionPascals.size()
            > limits.maximumPressureSamples
        || audit.corrections.size() > limits.maximumCorrections
        || audit.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional projection-energy validation limit exceeded");
    }
    if (audit != buildAudit(
                     grid, sweep, fragments, topology, nullptr, metric,
                     before, after, audit.pressureCorrectionPascals,
                     audit.settings, limits)) {
        throw std::invalid_argument(
            "planar regional projection-energy audit is corrupted");
    }
}

PlanarPressureRegionFragmentProjectionEnergyAudit
auditMovingPlanarPressureRegionFragmentProjectionEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    return buildAudit(
        grid, sweep, fragments, topology, &volumeRates, metric, before, after,
        pressureCorrectionPascals, settings, limits);
}

void validateMovingPlanarPressureRegionFragmentProjectionEnergyAudit(
    const PlanarPressureRegionFragmentProjectionEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits) {
    validateLimits(limits);
    if (audit.pressureCorrectionPascals.size()
            > limits.maximumPressureSamples
        || audit.corrections.size() > limits.maximumCorrections
        || audit.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "moving planar regional projection-energy validation limit "
            "exceeded");
    }
    if (audit != buildAudit(
                     grid, sweep, fragments, topology, &volumeRates, metric,
                     before, after, audit.pressureCorrectionPascals,
                     audit.settings, limits)) {
        throw std::invalid_argument(
            "moving planar regional projection-energy audit is corrupted");
    }
}

} // namespace simwing::fsi::fluid
