#include "fluid/planar_region_fragment_pressure_state.h"

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

void validateLimits(
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    if (limits.maximumControls == 0 || limits.maximumWalls == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional pressure-state limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional pressure-state storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional pressure-state storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t controlCount,
                              const std::size_t wallCount,
                              const std::size_t componentCount) {
    return checkedSum(
        checkedSum(
            checkedProduct(
                controlCount,
                sizeof(PlanarPressureRegionFragmentPressureStateControl)),
            checkedProduct(
                wallCount,
                sizeof(PlanarPressureRegionFragmentPressureStateWall))),
        checkedProduct(
            componentCount,
            sizeof(PlanarPressureRegionFragmentPressureStateComponent)));
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

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentPressureState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceFragmentFingerprint);
    fingerprint.integer(state.sourceTopologyFingerprint);
    fingerprint.integer(state.sourceMetricFingerprint);
    fingerprint.integer(state.sourceProjectionEnergyFingerprint);
    fingerprint.integer(state.sourcePressureJumpEnergyFingerprint);
    fingerprint.integer(state.volumeRateFingerprint);
    fingerprint.boolean(state.staticGeometry);
    fingerprint.boolean(state.usesMovingVolumeRates);
    fingerprint.real(state.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(state.controls.size()));
    for (const auto& control : state.controls) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.fragmentIndex));
        fingerprint.integer(control.fragmentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.integer(control.regionStableId);
        for (const double value : {
                 control.volumeCubicMeters,
                 control.authoredPressurePascals,
                 control.correctionPressurePascals,
                 control.totalPressurePascals}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.walls.size()));
    for (const auto& wall : state.walls) {
        for (const std::size_t value : {
                 wall.wallIndex,
                 wall.sourceFaceLinkIndex,
                 wall.minusFragmentIndex,
                 wall.plusFragmentIndex,
                 wall.minusComponentIndex,
                 wall.plusComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(wall.sourceFaceLinkStableId);
        fingerprint.integer(wall.surfaceStableId);
        fingerprint.enumeration(wall.axis);
        fingerprint.integer(wall.minusRegionStableId);
        fingerprint.integer(wall.plusRegionStableId);
        fingerprint.real(wall.areaSquareMeters);
        fingerprintVector(fingerprint, wall.unitNormalMinusToPlus);
        for (const double value : {
                 wall.minusTotalPressurePascals,
                 wall.plusTotalPressurePascals,
                 wall.authoredPressureJumpPascals,
                 wall.correctionPressureJumpPascals,
                 wall.totalPressureJumpPascals,
                 wall.pressureSplitResidualPascals}) {
            fingerprint.real(value);
        }
        fingerprintVector(
            fingerprint, wall.authoredPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, wall.correctionPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, wall.totalPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, wall.pressureForceSplitResidualNewtons);
        for (const double value : {
                 wall.materialWallVelocityMetersPerSecond,
                 wall.authoredPressureWorkToFluidJoules,
                 wall.correctionPressureWorkToFluidJoules,
                 wall.totalPressureWorkToFluidJoules,
                 wall.pressureWorkSplitResidualJoules,
                 wall.totalPressureWorkToSheetJoules,
                 wall.actionReactionWorkResidualJoules}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        state.components.size()));
    for (const auto& component : state.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        for (const double value : {
                 component.volumeCubicMeters,
                 component.authoredVolumeMeanPressurePascals,
                 component.correctionVolumeMeanPressurePascals,
                 component.totalVolumeMeanPressurePascals,
                 component.pressureMeanSplitResidualPascals,
                 component.authoredGeometryPressureWorkJoules,
                 component.correctionGeometryPressureWorkJoules,
                 component.totalGeometryPressureWorkJoules,
                 component.totalWallPressureWorkToFluidJoules,
                 component.wallGeometryWorkResidualJoules}) {
            fingerprint.real(value);
        }
    }
    for (const double value : {
             state.maximumAbsoluteCorrectionGaugePascals,
             state.maximumAbsolutePressureSplitResidualPascals,
             state.maximumAbsoluteForceSplitResidualNewtons,
             state.maximumAbsoluteWorkResidualJoules}) {
        fingerprint.real(value);
    }
    fingerprintVector(
        fingerprint, state.authoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, state.correctionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, state.totalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, state.pressureForceSplitResidualNewtons);
    for (const double value : {
             state.authoredPressureWorkToFluidJoules,
             state.correctionPressureWorkToFluidJoules,
             state.totalPressureWorkToFluidJoules,
             state.pressureWorkSplitResidualJoules,
             state.totalPressureWorkToSheetJoules,
             state.actionReactionWorkResidualJoules,
             state.totalGeometryPressureWorkJoules,
             state.wallGeometryWorkResidualJoules}) {
        fingerprint.real(value);
    }
    fingerprint.boolean(state.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        state.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentPressureState buildState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet* volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    validateLimits(limits);
    if (volumeRates == nullptr) {
        validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
            projection, grid, sweep, fragments, topology, metric, before,
            after, limits.projectionEnergyLimits);
        validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
            pressureJump, grid, sweep, fragments, topology, metric, after,
            limits.pressureJumpLimits);
    } else {
        validateMovingPlanarPressureRegionFragmentProjectionEnergyAudit(
            projection, grid, sweep, fragments, topology, *volumeRates,
            metric, before, after, limits.projectionEnergyLimits);
        validateMovingPlanarPressureRegionFragmentPressureJumpEnergyAudit(
            pressureJump, grid, sweep, fragments, topology, *volumeRates,
            metric, after, limits.pressureJumpLimits);
    }
    if (!projection.accepted || !pressureJump.accepted
        || projection.settings.timeStepSeconds
            != pressureJump.settings.timeStepSeconds
        || projection.afterVelocityStateFingerprint != after.fingerprint
        || pressureJump.sourceVelocityStateFingerprint != after.fingerprint
        || projection.sourceMetricFingerprint != metric.fingerprint
        || pressureJump.sourceMetricFingerprint != metric.fingerprint
        || projection.sourceTopologyFingerprint != topology.fingerprint
        || pressureJump.sourceTopologyFingerprint != topology.fingerprint
        || projection.sourceFragmentFingerprint != fragments.fingerprint
        || pressureJump.sourceFragmentFingerprint != fragments.fingerprint
        || projection.pressureCorrectionPascals.size()
            != fragments.fragments.size()
        || projection.staticGeometry != pressureJump.staticGeometry
        || projection.usesMovingVolumeRates
            != pressureJump.usesMovingVolumeRates
        || projection.volumeRateFingerprint
            != pressureJump.volumeRateFingerprint) {
        throw std::invalid_argument(
            "planar regional pressure-state sources are incompatible");
    }
    if (fragments.fragments.size() > limits.maximumControls
        || topology.pressureLayerWallLinkCount > limits.maximumWalls
        || topology.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional pressure-state entity limit exceeded");
    }

    PlanarPressureRegionFragmentPressureState result;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.sourceProjectionEnergyFingerprint = projection.fingerprint;
    result.sourcePressureJumpEnergyFingerprint = pressureJump.fingerprint;
    result.volumeRateFingerprint = projection.volumeRateFingerprint;
    result.staticGeometry = projection.staticGeometry;
    result.usesMovingVolumeRates = projection.usesMovingVolumeRates;
    result.timeStepSeconds = projection.settings.timeStepSeconds;
    result.ownedStorageBytes = ownedStorageBytes(
        fragments.fragments.size(), topology.pressureLayerWallLinkCount,
        topology.components.size());
    result.workingStorageBytes = checkedProduct(
        topology.links.size(), sizeof(std::size_t));
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional pressure-state storage limit exceeded");
    }
    result.controls.reserve(fragments.fragments.size());
    result.walls.reserve(topology.pressureLayerWallLinkCount);
    result.components.resize(topology.components.size());
    for (std::size_t index = 0; index < topology.components.size(); ++index) {
        const auto& source = topology.components[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.regionStableId = source.regionStableId;
        component.volumeCubicMeters = source.volumeCubicMeters;
        component.authoredGeometryPressureWorkJoules =
            pressureJump.components[index]
                .geometryPressureWorkToFluidJoules;
        component.correctionGeometryPressureWorkJoules =
            projection.components[index].geometryPressureWorkJoules;
        component.totalGeometryPressureWorkJoules =
            component.authoredGeometryPressureWorkJoules
            + component.correctionGeometryPressureWorkJoules;
    }

    double maximumAbsolutePressure = 0.0;
    for (std::size_t index = 0;
         index < fragments.fragments.size(); ++index) {
        const auto& fragment = fragments.fragments[index];
        const std::size_t componentIndex =
            topology.fragments[index].componentIndex;
        const double correction =
            projection.pressureCorrectionPascals[index];
        const double total = fragment.pressurePascals + correction;
        if (!std::isfinite(total)) {
            throw std::invalid_argument(
                "planar regional pressure-state pressure is not finite");
        }
        result.controls.push_back({
            index,
            fragment.stableId,
            componentIndex,
            fragment.regionStableId,
            fragment.volumeCubicMeters,
            fragment.pressurePascals,
            correction,
            total,
        });
        auto& component = result.components[componentIndex];
        component.authoredVolumeMeanPressurePascals +=
            fragment.pressurePascals * fragment.volumeCubicMeters;
        component.correctionVolumeMeanPressurePascals +=
            correction * fragment.volumeCubicMeters;
        component.totalVolumeMeanPressurePascals +=
            total * fragment.volumeCubicMeters;
        maximumAbsolutePressure = std::max({
            maximumAbsolutePressure,
            std::abs(fragment.pressurePascals),
            std::abs(correction),
            std::abs(total)});
    }
    const double pressureTolerance = scaledTolerance(
        std::max(
            projection.settings.absolutePressureGaugeTolerancePascals,
            pressureJump.settings
                .absolutePressureResidualTolerancePascals),
        std::max(
            projection.settings.relativePressureGaugeTolerance,
            pressureJump.settings.relativePressureResidualTolerance),
        {maximumAbsolutePressure});
    for (std::size_t index = 0; index < result.components.size(); ++index) {
        auto& component = result.components[index];
        component.authoredVolumeMeanPressurePascals /=
            component.volumeCubicMeters;
        component.correctionVolumeMeanPressurePascals /=
            component.volumeCubicMeters;
        component.totalVolumeMeanPressurePascals /=
            component.volumeCubicMeters;
        component.pressureMeanSplitResidualPascals =
            component.totalVolumeMeanPressurePascals
            - (component.authoredVolumeMeanPressurePascals
               + component.correctionVolumeMeanPressurePascals);
        result.maximumAbsoluteCorrectionGaugePascals = std::max(
            result.maximumAbsoluteCorrectionGaugePascals,
            std::abs(component.correctionVolumeMeanPressurePascals));
        result.maximumAbsolutePressureSplitResidualPascals = std::max(
            result.maximumAbsolutePressureSplitResidualPascals,
            std::abs(component.pressureMeanSplitResidualPascals));
        if (!std::isfinite(component.pressureMeanSplitResidualPascals)
            || std::abs(component.correctionVolumeMeanPressurePascals)
                > pressureTolerance
            || std::abs(component.pressureMeanSplitResidualPascals)
                > pressureTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-state component gauge failed");
        }
    }

    constexpr std::size_t absent = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> jumpLayerByLink(topology.links.size(), absent);
    for (std::size_t index = 0; index < pressureJump.layers.size(); ++index) {
        const std::size_t linkIndex =
            pressureJump.layers[index].sourceFaceLinkIndex;
        if (linkIndex >= topology.links.size()
            || jumpLayerByLink[linkIndex] != absent) {
            throw std::invalid_argument(
                "planar regional pressure-state jump mapping is invalid");
        }
        jumpLayerByLink[linkIndex] = index;
    }

    const double forceAbsolute =
        pressureJump.settings.absoluteForceResidualToleranceNewtons;
    const double forceRelative =
        pressureJump.settings.relativeForceResidualTolerance;
    const double workAbsolute = std::max(
        projection.settings.absoluteEnergyResidualToleranceJoules,
        pressureJump.settings.absoluteWorkResidualToleranceJoules);
    const double workRelative = std::max(
        projection.settings.relativeEnergyResidualTolerance,
        pressureJump.settings.relativeWorkResidualTolerance);
    double maximumForceNewtons = 0.0;
    double maximumWorkJoules = 0.0;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        const std::size_t jumpIndex = jumpLayerByLink[link.linkIndex];
        if (jumpIndex == absent) {
            throw std::invalid_argument(
                "planar regional pressure-state wall coverage is incomplete");
        }
        const auto& sourceJump = pressureJump.layers[jumpIndex];
        const auto& minus = result.controls[link.minusFragmentIndex];
        const auto& plus = result.controls[link.plusFragmentIndex];
        const double authoredJump = link.pressureJumpPascals;
        const double correctionJump = plus.correctionPressurePascals
            - minus.correctionPressurePascals;
        const double totalJump = plus.totalPressurePascals
            - minus.totalPressurePascals;
        const double pressureResidual = totalJump
            - (authoredJump + correctionJump);
        const Vector3 authoredSheetForce = scaledVector(
            link.unitNormalMinusToPlus,
            -authoredJump * link.areaSquareMeters);
        const Vector3 correctionSheetForce = scaledVector(
            link.unitNormalMinusToPlus,
            -correctionJump * link.areaSquareMeters);
        const Vector3 totalSheetForce = scaledVector(
            link.unitNormalMinusToPlus,
            -totalJump * link.areaSquareMeters);
        const Vector3 forceResidual = vectorDifference(
            totalSheetForce,
            vectorSum(authoredSheetForce, correctionSheetForce));
        const double wallVelocity =
            sourceJump.materialWallVelocityMetersPerSecond;
        const double authoredWork = authoredJump * link.areaSquareMeters
            * wallVelocity * result.timeStepSeconds;
        const double correctionWork = correctionJump * link.areaSquareMeters
            * wallVelocity * result.timeStepSeconds;
        const double totalWork = totalJump * link.areaSquareMeters
            * wallVelocity * result.timeStepSeconds;
        const double workSplitResidual =
            totalWork - (authoredWork + correctionWork);
        const double sheetWork = -totalWork;
        const double actionReactionWorkResidual = totalWork + sheetWork;
        maximumForceNewtons = std::max({
            maximumForceNewtons,
            maximumAbsoluteComponent(authoredSheetForce),
            maximumAbsoluteComponent(correctionSheetForce),
            maximumAbsoluteComponent(totalSheetForce)});
        maximumWorkJoules = std::max({
            maximumWorkJoules,
            std::abs(authoredWork),
            std::abs(correctionWork),
            std::abs(totalWork)});
        const double forceTolerance = scaledTolerance(
            forceAbsolute, forceRelative, {maximumForceNewtons});
        const double wallWorkTolerance = scaledTolerance(
            workAbsolute, workRelative, {maximumWorkJoules});
        if (!std::isfinite(pressureResidual)
            || std::abs(pressureResidual) > pressureTolerance
            || !finiteVector(forceResidual)
            || maximumAbsoluteComponent(forceResidual) > forceTolerance
            || !std::isfinite(workSplitResidual)
            || std::abs(workSplitResidual) > wallWorkTolerance
            || !std::isfinite(actionReactionWorkResidual)
            || std::abs(actionReactionWorkResidual) > wallWorkTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-state wall split failed");
        }
        result.walls.push_back({
            result.walls.size(),
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
            link.areaSquareMeters,
            link.unitNormalMinusToPlus,
            minus.totalPressurePascals,
            plus.totalPressurePascals,
            authoredJump,
            correctionJump,
            totalJump,
            pressureResidual,
            authoredSheetForce,
            correctionSheetForce,
            totalSheetForce,
            forceResidual,
            wallVelocity,
            authoredWork,
            correctionWork,
            totalWork,
            workSplitResidual,
            sheetWork,
            actionReactionWorkResidual,
        });
        result.maximumAbsolutePressureSplitResidualPascals = std::max(
            result.maximumAbsolutePressureSplitResidualPascals,
            std::abs(pressureResidual));
        result.maximumAbsoluteForceSplitResidualNewtons = std::max(
            result.maximumAbsoluteForceSplitResidualNewtons,
            maximumAbsoluteComponent(forceResidual));
        result.maximumAbsoluteWorkResidualJoules = std::max({
            result.maximumAbsoluteWorkResidualJoules,
            std::abs(workSplitResidual),
            std::abs(actionReactionWorkResidual)});
        addVector(
            result.authoredPressureForceOnSheetNewtons,
            authoredSheetForce);
        addVector(
            result.correctionPressureForceOnSheetNewtons,
            correctionSheetForce);
        addVector(
            result.totalPressureForceOnSheetNewtons,
            totalSheetForce);
        result.authoredPressureWorkToFluidJoules += authoredWork;
        result.correctionPressureWorkToFluidJoules += correctionWork;
        result.totalPressureWorkToFluidJoules += totalWork;

        result.components[link.minusComponentIndex]
            .totalWallPressureWorkToFluidJoules +=
            -minus.totalPressurePascals * link.areaSquareMeters
            * wallVelocity * result.timeStepSeconds;
        result.components[link.plusComponentIndex]
            .totalWallPressureWorkToFluidJoules +=
            plus.totalPressurePascals * link.areaSquareMeters
            * wallVelocity * result.timeStepSeconds;
    }
    if (result.walls.size() != topology.pressureLayerWallLinkCount) {
        throw std::invalid_argument(
            "planar regional pressure-state wall count is incomplete");
    }

    const double globalForceTolerance = scaledTolerance(
        forceAbsolute, forceRelative, {maximumForceNewtons});
    result.pressureForceSplitResidualNewtons = vectorDifference(
        result.totalPressureForceOnSheetNewtons,
        vectorSum(
            result.authoredPressureForceOnSheetNewtons,
            result.correctionPressureForceOnSheetNewtons));
    result.pressureWorkSplitResidualJoules =
        result.totalPressureWorkToFluidJoules
        - (result.authoredPressureWorkToFluidJoules
           + result.correctionPressureWorkToFluidJoules);
    result.totalPressureWorkToSheetJoules =
        -result.totalPressureWorkToFluidJoules;
    result.actionReactionWorkResidualJoules =
        result.totalPressureWorkToFluidJoules
        + result.totalPressureWorkToSheetJoules;
    result.totalGeometryPressureWorkJoules =
        pressureJump.geometryPressureWorkToFluidJoules
        + projection.geometryPressureWorkJoules;
    result.wallGeometryWorkResidualJoules =
        result.totalPressureWorkToFluidJoules
        - result.totalGeometryPressureWorkJoules;
    result.maximumAbsoluteForceSplitResidualNewtons = std::max(
        result.maximumAbsoluteForceSplitResidualNewtons,
        maximumAbsoluteComponent(result.pressureForceSplitResidualNewtons));
    result.maximumAbsoluteWorkResidualJoules = std::max({
        result.maximumAbsoluteWorkResidualJoules,
        std::abs(result.pressureWorkSplitResidualJoules),
        std::abs(result.actionReactionWorkResidualJoules),
        std::abs(result.wallGeometryWorkResidualJoules),
        std::abs(
            result.authoredPressureWorkToFluidJoules
            - pressureJump.pressureJumpWorkToFluidJoules),
        std::abs(
            result.correctionPressureWorkToFluidJoules
            - projection.geometryPressureWorkJoules)});
    const double globalWorkTolerance = scaledTolerance(
        workAbsolute, workRelative,
        {maximumWorkJoules,
         result.authoredPressureWorkToFluidJoules,
         result.correctionPressureWorkToFluidJoules,
         result.totalPressureWorkToFluidJoules,
         result.totalGeometryPressureWorkJoules});
    if (!finiteVector(result.pressureForceSplitResidualNewtons)
        || maximumAbsoluteComponent(result.pressureForceSplitResidualNewtons)
            > globalForceTolerance
        || !std::isfinite(result.pressureWorkSplitResidualJoules)
        || std::abs(result.pressureWorkSplitResidualJoules)
            > globalWorkTolerance
        || !std::isfinite(result.actionReactionWorkResidualJoules)
        || std::abs(result.actionReactionWorkResidualJoules)
            > globalWorkTolerance
        || !std::isfinite(result.wallGeometryWorkResidualJoules)
        || std::abs(result.wallGeometryWorkResidualJoules)
            > globalWorkTolerance
        || std::abs(
               result.authoredPressureWorkToFluidJoules
               - pressureJump.pressureJumpWorkToFluidJoules)
            > globalWorkTolerance
        || std::abs(
               result.correctionPressureWorkToFluidJoules
               - projection.geometryPressureWorkJoules)
            > globalWorkTolerance) {
        throw std::invalid_argument(
            "planar regional pressure-state global work closure failed");
    }

    for (auto& component : result.components) {
        component.wallGeometryWorkResidualJoules =
            component.totalWallPressureWorkToFluidJoules
            - component.totalGeometryPressureWorkJoules;
        result.maximumAbsoluteWorkResidualJoules = std::max(
            result.maximumAbsoluteWorkResidualJoules,
            std::abs(component.wallGeometryWorkResidualJoules));
        const double componentWorkTolerance = scaledTolerance(
            workAbsolute, workRelative,
            {maximumWorkJoules,
             component.totalWallPressureWorkToFluidJoules,
             component.totalGeometryPressureWorkJoules});
        if (!std::isfinite(component.wallGeometryWorkResidualJoules)
            || std::abs(component.wallGeometryWorkResidualJoules)
                > componentWorkTolerance) {
            throw std::invalid_argument(
                "planar regional pressure-state component work closure "
                "failed");
        }
    }

    result.accepted = true;
    result.fingerprint = stateFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentPressureState
composeStaticPlanarPressureRegionFragmentPressureState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    return buildState(
        grid, sweep, fragments, topology, nullptr, metric, before, after,
        projection, pressureJump, limits);
}

void validateStaticPlanarPressureRegionFragmentPressureState(
    const PlanarPressureRegionFragmentPressureState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    validateLimits(limits);
    if (state.controls.size() > limits.maximumControls
        || state.walls.size() > limits.maximumWalls
        || state.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional pressure-state validation limit exceeded");
    }
    if (state != buildState(
                     grid, sweep, fragments, topology, nullptr, metric,
                     before, after, projection, pressureJump, limits)) {
        throw std::invalid_argument(
            "planar regional pressure state is corrupted");
    }
}

PlanarPressureRegionFragmentPressureState
composeMovingPlanarPressureRegionFragmentPressureState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    return buildState(
        grid, sweep, fragments, topology, &volumeRates, metric, before, after,
        projection, pressureJump, limits);
}

void validateMovingPlanarPressureRegionFragmentPressureState(
    const PlanarPressureRegionFragmentPressureState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyAudit& projection,
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& pressureJump,
    const PlanarPressureRegionFragmentPressureStateLimits& limits) {
    validateLimits(limits);
    if (state.controls.size() > limits.maximumControls
        || state.walls.size() > limits.maximumWalls
        || state.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "moving planar regional pressure-state validation limit "
            "exceeded");
    }
    if (state != buildState(
                     grid, sweep, fragments, topology, &volumeRates, metric,
                     before, after, projection, pressureJump, limits)) {
        throw std::invalid_argument(
            "moving planar regional pressure state is corrupted");
    }
}

} // namespace simwing::fsi::fluid
