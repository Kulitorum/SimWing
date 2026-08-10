#include "fluid/planar_region_fragment_pressure_jump_energy.h"

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
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings) {
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "planar regional pressure-jump energy time step is invalid");
    }
    validateTolerancePair(
        settings.absolutePressureResidualTolerancePascals,
        settings.relativePressureResidualTolerance,
        "planar regional pressure-jump pressure tolerance is invalid");
    validateTolerancePair(
        settings.absoluteVelocityResidualToleranceMetersPerSecond,
        settings.relativeVelocityResidualTolerance,
        "planar regional pressure-jump velocity tolerance is invalid");
    validateTolerancePair(
        settings.absoluteForceResidualToleranceNewtons,
        settings.relativeForceResidualTolerance,
        "planar regional pressure-jump force tolerance is invalid");
    validateTolerancePair(
        settings.absoluteWorkResidualToleranceJoules,
        settings.relativeWorkResidualTolerance,
        "planar regional pressure-jump work tolerance is invalid");
}

void validateLimits(
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    if (limits.maximumLayers == 0 || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional pressure-jump energy limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional pressure-jump energy storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional pressure-jump energy storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t layerCount,
                              const std::size_t componentCount) {
    return checkedSum(
        checkedProduct(
            layerCount,
            sizeof(
                PlanarPressureRegionFragmentPressureJumpEnergyLayer)),
        checkedProduct(
            componentCount,
            sizeof(
                PlanarPressureRegionFragmentPressureJumpEnergyComponent)));
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

double scaledTolerance(const double absolute,
                       const double relative,
                       const std::initializer_list<double> scales) {
    double scale = 0.0;
    for (const double value : scales) {
        scale = std::max(scale, std::abs(value));
    }
    return std::max(absolute, relative * scale);
}

Vector3 scaledVector(const Vector3& value, const double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vector3 vectorSum(const Vector3& first, const Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

Vector3 vectorDifference(const Vector3& first, const Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

void addVector(Vector3& target, const Vector3& value) {
    target.x += value.x;
    target.y += value.y;
    target.z += value.z;
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings) {
    for (const double value : {
             settings.timeStepSeconds,
             settings.absolutePressureResidualTolerancePascals,
             settings.relativePressureResidualTolerance,
             settings.absoluteVelocityResidualToleranceMetersPerSecond,
             settings.relativeVelocityResidualTolerance,
             settings.absoluteForceResidualToleranceNewtons,
             settings.relativeForceResidualTolerance,
             settings.absoluteWorkResidualToleranceJoules,
             settings.relativeWorkResidualTolerance}) {
        fingerprint.real(value);
    }
}

std::uint64_t auditFingerprint(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.sourceMetricFingerprint);
    fingerprint.integer(audit.sourceVelocityStateFingerprint);
    fingerprint.integer(audit.sourceTopologyFingerprint);
    fingerprint.integer(audit.sourceFragmentFingerprint);
    fingerprint.integer(audit.volumeRateFingerprint);
    fingerprint.boolean(audit.staticGeometry);
    fingerprint.boolean(audit.usesMovingVolumeRates);
    fingerprintSettings(fingerprint, audit.settings);
    fingerprint.integer(static_cast<std::uint64_t>(audit.layers.size()));
    for (const auto& layer : audit.layers) {
        for (const std::size_t value : {
                 layer.layerIndex,
                 layer.sourceFaceLinkIndex,
                 layer.minusFragmentIndex,
                 layer.plusFragmentIndex,
                 layer.minusComponentIndex,
                 layer.plusComponentIndex,
                 layer.minusTraceDofIndex,
                 layer.plusTraceDofIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(layer.sourceFaceLinkStableId);
        fingerprint.integer(layer.surfaceStableId);
        fingerprint.enumeration(layer.axis);
        fingerprint.integer(layer.minusRegionStableId);
        fingerprint.integer(layer.plusRegionStableId);
        fingerprint.real(layer.areaSquareMeters);
        fingerprintVector(fingerprint, layer.unitNormalMinusToPlus);
        for (const double value : {
                 layer.minusPressurePascals,
                 layer.plusPressurePascals,
                 layer.authoredPressureJumpPascals,
                 layer.reconstructedPressureJumpPascals,
                 layer.pressureJumpResidualPascals,
                 layer.materialWallVelocityMetersPerSecond,
                 layer.minusTraceVelocityMetersPerSecond,
                 layer.plusTraceVelocityMetersPerSecond,
                 layer.maximumAbsoluteWallVelocityResidualMetersPerSecond}) {
            fingerprint.real(value);
        }
        fingerprintVector(
            fingerprint, layer.resolvedPressureForceOnMinusFluidNewtons);
        fingerprintVector(
            fingerprint, layer.resolvedPressureForceOnPlusFluidNewtons);
        fingerprintVector(
            fingerprint, layer.resolvedPressureForceOnFluidNewtons);
        fingerprintVector(
            fingerprint, layer.authoredPressureJumpForceOnFluidNewtons);
        fingerprintVector(
            fingerprint, layer.pressureForceClosureResidualNewtons);
        fingerprintVector(
            fingerprint, layer.pressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, layer.actionReactionForceResidualNewtons);
        fingerprintVector(
            fingerprint, layer.pressureJumpImpulseOnFluidNewtonSeconds);
        fingerprintVector(
            fingerprint, layer.pressureImpulseOnSheetNewtonSeconds);
        fingerprintVector(
            fingerprint, layer.actionReactionImpulseResidualNewtonSeconds);
        for (const double value : {
                 layer.resolvedPressurePowerToFluidWatts,
                 layer.authoredPressureJumpPowerToFluidWatts,
                 layer.pressurePowerClosureResidualWatts,
                 layer.pressureJumpWorkToFluidJoules,
                 layer.pressureWorkToSheetJoules,
                 layer.actionReactionWorkResidualJoules}) {
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
            component.pressureLayerSideCount));
        fingerprintVector(
            fingerprint, component.resolvedPressureForceOnFluidNewtons);
        fingerprintVector(
            fingerprint, component.closedBoundaryForceResidualNewtons);
        fingerprintVector(
            fingerprint, component.pressureImpulseOnFluidNewtonSeconds);
        for (const double value : {
                 component.resolvedPressurePowerToFluidWatts,
                 component.pressureWorkToFluidJoules,
                 component.geometryPressureWorkToFluidJoules,
                 component.workGeometryResidualJoules}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.pressureLayerTraceCount));
    for (const double value : {
             audit.maximumAbsolutePressureJumpResidualPascals,
             audit.maximumAbsoluteWallVelocityResidualMetersPerSecond,
             audit.maximumAbsoluteForceClosureResidualNewtons,
             audit
                 .maximumAbsoluteComponentClosedBoundaryForceResidualNewtons,
             audit.maximumAbsoluteWorkClosureResidualJoules}) {
        fingerprint.real(value);
    }
    fingerprintVector(
        fingerprint, audit.resolvedPressureForceOnFluidNewtons);
    fingerprintVector(
        fingerprint, audit.authoredPressureJumpForceOnFluidNewtons);
    fingerprintVector(
        fingerprint, audit.pressureForceClosureResidualNewtons);
    fingerprintVector(
        fingerprint, audit.pressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, audit.actionReactionForceResidualNewtons);
    fingerprintVector(
        fingerprint, audit.pressureJumpImpulseOnFluidNewtonSeconds);
    fingerprintVector(
        fingerprint, audit.pressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint, audit.actionReactionImpulseResidualNewtonSeconds);
    for (const double value : {
             audit.resolvedPressurePowerToFluidWatts,
             audit.authoredPressureJumpPowerToFluidWatts,
             audit.pressurePowerClosureResidualWatts,
             audit.pressureJumpWorkToFluidJoules,
             audit.pressureWorkToSheetJoules,
             audit.actionReactionWorkResidualJoules,
             audit.geometryPressureWorkToFluidJoules,
             audit.workGeometryResidualJoules}) {
        fingerprint.real(value);
    }
    fingerprint.boolean(audit.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentPressureJumpEnergyAudit buildAudit(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet* volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    const bool staticGeometry = isStaticGeometry(sweep);
    if (volumeRates == nullptr && !staticGeometry) {
        throw std::invalid_argument(
            "planar regional pressure-jump energy audit requires static "
            "geometry");
    }
    if (volumeRates != nullptr) {
        validatePlanarPressureRegionFragmentVolumeRates(
            *volumeRates, grid, sweep, fragments, topology,
            limits.volumeRateLimits);
        if (volumeRates->durationSeconds != settings.timeStepSeconds) {
            throw std::invalid_argument(
                "planar regional pressure-jump energy duration does not "
                "match its volume rates");
        }
    }
    validatePlanarPressureRegionFragmentVelocityState(
        velocityState, grid, sweep, fragments, topology, metric,
        limits.velocityStateLimits);
    if (topology.pressureLayerWallLinkCount > limits.maximumLayers
        || topology.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional pressure-jump energy entity limit exceeded");
    }

    PlanarPressureRegionFragmentPressureJumpEnergyAudit result;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.sourceVelocityStateFingerprint = velocityState.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.volumeRateFingerprint =
        volumeRates == nullptr ? 0 : volumeRates->fingerprint;
    result.staticGeometry = staticGeometry;
    result.usesMovingVolumeRates = volumeRates != nullptr;
    result.settings = settings;
    result.ownedStorageBytes = ownedStorageBytes(
        topology.pressureLayerWallLinkCount, topology.components.size());
    result.workingStorageBytes = checkedProduct(
        checkedProduct(topology.links.size(), sizeof(std::size_t)), 2);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional pressure-jump energy storage limit exceeded");
    }
    result.layers.reserve(topology.pressureLayerWallLinkCount);
    result.components.resize(topology.components.size());
    for (std::size_t index = 0; index < topology.components.size(); ++index) {
        const auto& source = topology.components[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.regionStableId = source.regionStableId;
    }

    constexpr std::size_t absent = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> minusTrace(topology.links.size(), absent);
    std::vector<std::size_t> plusTrace(topology.links.size(), absent);
    for (const auto& dof : metric.dofs) {
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            continue;
        }
        if (dof.sourceFaceLinkIndex >= topology.links.size()
            || topology.links[dof.sourceFaceLinkIndex].kind
                != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            throw std::invalid_argument(
                "planar regional pressure-jump trace source is invalid");
        }
        auto& destination = dof.kind
                == PlanarPressureRegionFragmentVelocityDofKind::
                    PressureLayerMinusTrace
            ? minusTrace[dof.sourceFaceLinkIndex]
            : plusTrace[dof.sourceFaceLinkIndex];
        if (destination != absent) {
            throw std::invalid_argument(
                "planar regional pressure-jump trace is duplicated");
        }
        destination = dof.dofIndex;
    }

    double maximumLayerForceNewtons = 0.0;
    double maximumLayerWorkJoules = 0.0;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        const std::size_t minusDofIndex = minusTrace[link.linkIndex];
        const std::size_t plusDofIndex = plusTrace[link.linkIndex];
        if (minusDofIndex == absent || plusDofIndex == absent) {
            throw std::invalid_argument(
                "planar regional pressure-jump trace coverage is incomplete");
        }
        const auto& minusDof = metric.dofs[minusDofIndex];
        const auto& plusDof = metric.dofs[plusDofIndex];
        const auto& minusSample = velocityState.samples[minusDofIndex];
        const auto& plusSample = velocityState.samples[plusDofIndex];
        const auto& minusFragment =
            fragments.fragments[link.minusFragmentIndex];
        const auto& plusFragment =
            fragments.fragments[link.plusFragmentIndex];
        double expectedMinusVelocity = 0.0;
        double expectedPlusVelocity = 0.0;
        if (volumeRates != nullptr) {
            expectedMinusVelocity = volumeRates
                ->fragments[link.minusFragmentIndex]
                .upperBoundaryVelocityMetersPerSecond;
            expectedPlusVelocity = volumeRates
                ->fragments[link.plusFragmentIndex]
                .lowerBoundaryVelocityMetersPerSecond;
        }
        const double minusVelocityResidual =
            minusSample.normalVelocityMetersPerSecond
            - expectedMinusVelocity;
        const double plusVelocityResidual =
            plusSample.normalVelocityMetersPerSecond
            - expectedPlusVelocity;
        const double materialVelocityResidual =
            expectedPlusVelocity - expectedMinusVelocity;
        const double maximumVelocityResidual = std::max({
            std::abs(minusVelocityResidual),
            std::abs(plusVelocityResidual),
            std::abs(materialVelocityResidual)});
        const double velocityTolerance = scaledTolerance(
            settings.absoluteVelocityResidualToleranceMetersPerSecond,
            settings.relativeVelocityResidualTolerance,
            {expectedMinusVelocity, expectedPlusVelocity,
             minusSample.normalVelocityMetersPerSecond,
             plusSample.normalVelocityMetersPerSecond});

        const double reconstructedJump =
            plusFragment.pressurePascals
            - minusFragment.pressurePascals;
        const double pressureResidual =
            reconstructedJump - link.pressureJumpPascals;
        const double pressureTolerance = scaledTolerance(
            settings.absolutePressureResidualTolerancePascals,
            settings.relativePressureResidualTolerance,
            {minusFragment.pressurePascals,
             plusFragment.pressurePascals,
             link.pressureJumpPascals,
             reconstructedJump});
        if (!std::isfinite(maximumVelocityResidual)
            || maximumVelocityResidual > velocityTolerance
            || !std::isfinite(pressureResidual)
            || std::abs(pressureResidual) > pressureTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-jump wall source is inconsistent");
        }

        const Vector3 minusForce = scaledVector(
            link.unitNormalMinusToPlus,
            -minusFragment.pressurePascals * link.areaSquareMeters);
        const Vector3 plusForce = scaledVector(
            link.unitNormalMinusToPlus,
            plusFragment.pressurePascals * link.areaSquareMeters);
        const Vector3 resolvedForce = vectorSum(minusForce, plusForce);
        const Vector3 authoredForce = scaledVector(
            link.unitNormalMinusToPlus,
            link.pressureJumpPascals * link.areaSquareMeters);
        const Vector3 forceClosureResidual = vectorDifference(
            resolvedForce, authoredForce);
        const Vector3 sheetForce = scaledVector(authoredForce, -1.0);
        const Vector3 actionReactionForceResidual = vectorSum(
            authoredForce, sheetForce);
        const Vector3 fluidImpulse = scaledVector(
            authoredForce, settings.timeStepSeconds);
        const Vector3 sheetImpulse = scaledVector(
            sheetForce, settings.timeStepSeconds);
        const Vector3 actionReactionImpulseResidual = vectorSum(
            fluidImpulse, sheetImpulse);
        const double resolvedPower = link.areaSquareMeters
            * (-minusFragment.pressurePascals
                   * minusSample.normalVelocityMetersPerSecond
               + plusFragment.pressurePascals
                   * plusSample.normalVelocityMetersPerSecond);
        const double authoredPower = link.pressureJumpPascals
            * link.areaSquareMeters * expectedMinusVelocity;
        const double powerClosureResidual =
            resolvedPower - authoredPower;
        const double fluidWork =
            authoredPower * settings.timeStepSeconds;
        const double sheetWork =
            -authoredPower * settings.timeStepSeconds;
        const double actionReactionWorkResidual = fluidWork + sheetWork;
        const double forceScale = std::max({
            maximumAbsoluteComponent(minusForce),
            maximumAbsoluteComponent(plusForce),
            maximumAbsoluteComponent(resolvedForce),
            maximumAbsoluteComponent(authoredForce)});
        const double forceTolerance = scaledTolerance(
            settings.absoluteForceResidualToleranceNewtons,
            settings.relativeForceResidualTolerance, {forceScale});
        const double workScale = settings.timeStepSeconds * std::max(
            std::abs(resolvedPower), std::abs(authoredPower));
        const double workTolerance = scaledTolerance(
            settings.absoluteWorkResidualToleranceJoules,
            settings.relativeWorkResidualTolerance,
            {workScale, fluidWork, sheetWork});
        if (!finiteVector(forceClosureResidual)
            || !finiteVector(actionReactionForceResidual)
            || !finiteVector(actionReactionImpulseResidual)
            || maximumAbsoluteComponent(forceClosureResidual)
                > forceTolerance
            || maximumAbsoluteComponent(actionReactionForceResidual)
                > forceTolerance
            || maximumAbsoluteComponent(actionReactionImpulseResidual)
                > forceTolerance * settings.timeStepSeconds
            || !std::isfinite(powerClosureResidual)
            || std::abs(powerClosureResidual) * settings.timeStepSeconds
                > workTolerance
            || !std::isfinite(actionReactionWorkResidual)
            || std::abs(actionReactionWorkResidual) > workTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-jump layer closure failed");
        }

        result.layers.push_back({
            result.layers.size(),
            link.linkIndex,
            link.stableId,
            link.surfaceStableId,
            link.axis,
            link.minusFragmentIndex,
            link.plusFragmentIndex,
            link.minusComponentIndex,
            link.plusComponentIndex,
            link.minusRegionStableId,
            link.plusRegionStableId,
            minusDofIndex,
            plusDofIndex,
            link.areaSquareMeters,
            link.unitNormalMinusToPlus,
            minusFragment.pressurePascals,
            plusFragment.pressurePascals,
            link.pressureJumpPascals,
            reconstructedJump,
            pressureResidual,
            expectedMinusVelocity,
            minusSample.normalVelocityMetersPerSecond,
            plusSample.normalVelocityMetersPerSecond,
            maximumVelocityResidual,
            minusForce,
            plusForce,
            resolvedForce,
            authoredForce,
            forceClosureResidual,
            sheetForce,
            actionReactionForceResidual,
            fluidImpulse,
            sheetImpulse,
            actionReactionImpulseResidual,
            resolvedPower,
            authoredPower,
            powerClosureResidual,
            fluidWork,
            sheetWork,
            actionReactionWorkResidual,
        });
        result.pressureLayerTraceCount += 2;
        result.maximumAbsolutePressureJumpResidualPascals = std::max(
            result.maximumAbsolutePressureJumpResidualPascals,
            std::abs(pressureResidual));
        result.maximumAbsoluteWallVelocityResidualMetersPerSecond =
            std::max(
                result.maximumAbsoluteWallVelocityResidualMetersPerSecond,
                maximumVelocityResidual);
        result.maximumAbsoluteForceClosureResidualNewtons = std::max({
            result.maximumAbsoluteForceClosureResidualNewtons,
            maximumAbsoluteComponent(forceClosureResidual),
            maximumAbsoluteComponent(actionReactionForceResidual)});
        result.maximumAbsoluteWorkClosureResidualJoules = std::max({
            result.maximumAbsoluteWorkClosureResidualJoules,
            std::abs(powerClosureResidual) * settings.timeStepSeconds,
            std::abs(actionReactionWorkResidual)});
        maximumLayerForceNewtons = std::max(
            maximumLayerForceNewtons, forceScale);
        maximumLayerWorkJoules = std::max(
            maximumLayerWorkJoules,
            std::max(std::abs(fluidWork), std::abs(sheetWork)));

        auto& minusComponent =
            result.components[link.minusComponentIndex];
        ++minusComponent.pressureLayerSideCount;
        addVector(
            minusComponent.resolvedPressureForceOnFluidNewtons,
            minusForce);
        minusComponent.resolvedPressurePowerToFluidWatts +=
            -minusFragment.pressurePascals * link.areaSquareMeters
            * minusSample.normalVelocityMetersPerSecond;
        auto& plusComponent =
            result.components[link.plusComponentIndex];
        ++plusComponent.pressureLayerSideCount;
        addVector(
            plusComponent.resolvedPressureForceOnFluidNewtons,
            plusForce);
        plusComponent.resolvedPressurePowerToFluidWatts +=
            plusFragment.pressurePascals * link.areaSquareMeters
            * plusSample.normalVelocityMetersPerSecond;

        addVector(result.resolvedPressureForceOnFluidNewtons, resolvedForce);
        addVector(
            result.authoredPressureJumpForceOnFluidNewtons, authoredForce);
        addVector(result.pressureForceOnSheetNewtons, sheetForce);
        result.resolvedPressurePowerToFluidWatts += resolvedPower;
        result.authoredPressureJumpPowerToFluidWatts += authoredPower;
    }
    if (result.layers.size() != topology.pressureLayerWallLinkCount
        || result.pressureLayerTraceCount
            != metric.pressureLayerTraceDofCount) {
        throw std::invalid_argument(
            "planar regional pressure-jump wall coverage is incomplete");
    }

    if (volumeRates != nullptr) {
        for (std::size_t index = 0;
             index < fragments.fragments.size(); ++index) {
            const double work = -settings.timeStepSeconds
                * fragments.fragments[index].pressurePascals
                * volumeRates->fragments[index]
                    .geometryVolumeChangeRateCubicMetersPerSecond;
            result.components[topology.fragments[index].componentIndex]
                .geometryPressureWorkToFluidJoules += work;
            result.geometryPressureWorkToFluidJoules += work;
        }
    }

    const double forceTolerance = scaledTolerance(
        settings.absoluteForceResidualToleranceNewtons,
        settings.relativeForceResidualTolerance,
        {maximumLayerForceNewtons,
         maximumAbsoluteComponent(result.resolvedPressureForceOnFluidNewtons),
         maximumAbsoluteComponent(
             result.authoredPressureJumpForceOnFluidNewtons)});
    for (auto& component : result.components) {
        component.closedBoundaryForceResidualNewtons =
            component.resolvedPressureForceOnFluidNewtons;
        component.pressureImpulseOnFluidNewtonSeconds = scaledVector(
            component.resolvedPressureForceOnFluidNewtons,
            settings.timeStepSeconds);
        component.pressureWorkToFluidJoules =
            component.resolvedPressurePowerToFluidWatts
            * settings.timeStepSeconds;
        component.workGeometryResidualJoules =
            component.pressureWorkToFluidJoules
            - component.geometryPressureWorkToFluidJoules;
        result.maximumAbsoluteComponentClosedBoundaryForceResidualNewtons =
            std::max(
                result
                    .maximumAbsoluteComponentClosedBoundaryForceResidualNewtons,
                maximumAbsoluteComponent(
                    component.closedBoundaryForceResidualNewtons));
        const double componentWorkTolerance = scaledTolerance(
            settings.absoluteWorkResidualToleranceJoules,
            settings.relativeWorkResidualTolerance,
            {maximumLayerWorkJoules,
             component.pressureWorkToFluidJoules,
             component.geometryPressureWorkToFluidJoules});
        result.maximumAbsoluteWorkClosureResidualJoules = std::max(
            result.maximumAbsoluteWorkClosureResidualJoules,
            std::abs(component.workGeometryResidualJoules));
        if (!finiteVector(component.closedBoundaryForceResidualNewtons)
            || maximumAbsoluteComponent(
                   component.closedBoundaryForceResidualNewtons)
                > forceTolerance
            || !std::isfinite(component.workGeometryResidualJoules)
            || std::abs(component.workGeometryResidualJoules)
                > componentWorkTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-jump component closure failed");
        }
    }

    result.pressureForceClosureResidualNewtons = vectorDifference(
        result.resolvedPressureForceOnFluidNewtons,
        result.authoredPressureJumpForceOnFluidNewtons);
    result.actionReactionForceResidualNewtons = vectorSum(
        result.authoredPressureJumpForceOnFluidNewtons,
        result.pressureForceOnSheetNewtons);
    result.pressureJumpImpulseOnFluidNewtonSeconds = scaledVector(
        result.authoredPressureJumpForceOnFluidNewtons,
        settings.timeStepSeconds);
    result.pressureImpulseOnSheetNewtonSeconds = scaledVector(
        result.pressureForceOnSheetNewtons, settings.timeStepSeconds);
    result.actionReactionImpulseResidualNewtonSeconds = vectorSum(
        result.pressureJumpImpulseOnFluidNewtonSeconds,
        result.pressureImpulseOnSheetNewtonSeconds);
    result.pressurePowerClosureResidualWatts =
        result.resolvedPressurePowerToFluidWatts
        - result.authoredPressureJumpPowerToFluidWatts;
    result.pressureJumpWorkToFluidJoules =
        result.authoredPressureJumpPowerToFluidWatts
        * settings.timeStepSeconds;
    result.pressureWorkToSheetJoules =
        -result.pressureJumpWorkToFluidJoules;
    result.actionReactionWorkResidualJoules =
        result.pressureJumpWorkToFluidJoules
        + result.pressureWorkToSheetJoules;
    result.workGeometryResidualJoules =
        result.pressureJumpWorkToFluidJoules
        - result.geometryPressureWorkToFluidJoules;
    result.maximumAbsoluteForceClosureResidualNewtons = std::max({
        result.maximumAbsoluteForceClosureResidualNewtons,
        maximumAbsoluteComponent(result.pressureForceClosureResidualNewtons),
        maximumAbsoluteComponent(result.actionReactionForceResidualNewtons),
        maximumAbsoluteComponent(
            result.actionReactionImpulseResidualNewtonSeconds)
            / settings.timeStepSeconds});
    result.maximumAbsoluteWorkClosureResidualJoules = std::max({
        result.maximumAbsoluteWorkClosureResidualJoules,
        std::abs(result.pressurePowerClosureResidualWatts)
            * settings.timeStepSeconds,
        std::abs(result.actionReactionWorkResidualJoules),
        std::abs(result.workGeometryResidualJoules)});
    const double workTolerance = scaledTolerance(
        settings.absoluteWorkResidualToleranceJoules,
        settings.relativeWorkResidualTolerance,
        {maximumLayerWorkJoules,
         result.pressureJumpWorkToFluidJoules,
         result.pressureWorkToSheetJoules,
         result.geometryPressureWorkToFluidJoules});
    if (!finiteVector(result.pressureForceClosureResidualNewtons)
        || !finiteVector(result.actionReactionForceResidualNewtons)
        || !finiteVector(result.actionReactionImpulseResidualNewtonSeconds)
        || maximumAbsoluteComponent(
               result.pressureForceClosureResidualNewtons)
            > forceTolerance
        || maximumAbsoluteComponent(
               result.actionReactionForceResidualNewtons)
            > forceTolerance
        || maximumAbsoluteComponent(
               result.authoredPressureJumpForceOnFluidNewtons)
            > forceTolerance
        || maximumAbsoluteComponent(
               result.resolvedPressureForceOnFluidNewtons)
            > forceTolerance
        || maximumAbsoluteComponent(
               result.actionReactionImpulseResidualNewtonSeconds)
            > forceTolerance * settings.timeStepSeconds
        || !std::isfinite(result.pressurePowerClosureResidualWatts)
        || std::abs(result.pressurePowerClosureResidualWatts)
               * settings.timeStepSeconds
            > workTolerance
        || !std::isfinite(result.actionReactionWorkResidualJoules)
        || std::abs(result.actionReactionWorkResidualJoules)
            > workTolerance
        || !std::isfinite(result.workGeometryResidualJoules)
        || std::abs(result.workGeometryResidualJoules) > workTolerance) {
        throw std::invalid_argument(
            "planar regional pressure-jump global closure failed");
    }

    result.accepted = true;
    result.fingerprint = auditFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentPressureJumpEnergyAudit
auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    return buildAudit(
        grid, sweep, fragments, topology, nullptr, metric, velocityState,
        settings, limits);
}

void validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    validateLimits(limits);
    if (audit.layers.size() > limits.maximumLayers
        || audit.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional pressure-jump energy validation limit "
            "exceeded");
    }
    if (audit != buildAudit(
                     grid, sweep, fragments, topology, nullptr, metric,
                     velocityState, audit.settings, limits)) {
        throw std::invalid_argument(
            "planar regional pressure-jump energy audit is corrupted");
    }
}

PlanarPressureRegionFragmentPressureJumpEnergyAudit
auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergySettings& settings,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    return buildAudit(
        grid, sweep, fragments, topology, &volumeRates, metric,
        velocityState, settings, limits);
}

void validateMovingPlanarPressureRegionFragmentPressureJumpEnergyAudit(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& velocityState,
    const PlanarPressureRegionFragmentPressureJumpEnergyLimits& limits) {
    validateLimits(limits);
    if (audit.layers.size() > limits.maximumLayers
        || audit.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "moving planar regional pressure-jump energy validation limit "
            "exceeded");
    }
    if (audit != buildAudit(
                     grid, sweep, fragments, topology, &volumeRates, metric,
                     velocityState, audit.settings, limits)) {
        throw std::invalid_argument(
            "moving planar regional pressure-jump energy audit is "
            "corrupted");
    }
}

} // namespace simwing::fsi::fluid
