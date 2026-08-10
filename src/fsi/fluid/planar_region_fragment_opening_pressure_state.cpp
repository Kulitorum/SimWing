#include "fluid/planar_region_fragment_opening_pressure_state.h"

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
    template<typename Unsigned> void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    template<typename Enum> void enumeration(const Enum value) {
        using Underlying = std::underlying_type_t<Enum>;
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

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening pressure-state storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening pressure-state storage overflows");
    }
    return first + second;
}

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
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings) {
    validateTolerancePair(
        settings.absolutePressureResidualTolerancePascals,
        settings.relativePressureResidualTolerance,
        "opening pressure-state pressure tolerances are invalid");
    validateTolerancePair(
        settings.absoluteForceResidualToleranceNewtons,
        settings.relativeForceResidualTolerance,
        "opening pressure-state force tolerances are invalid");
    validateTolerancePair(
        settings.absoluteWorkResidualToleranceJoules,
        settings.relativeWorkResidualTolerance,
        "opening pressure-state work tolerances are invalid");
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits) {
    if (limits.maximumControls == 0 || limits.maximumWalls == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening pressure-state limits are invalid");
    }
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
                       const std::initializer_list<double> values) {
    double scale = 0.0;
    for (const double value : values)
        scale = std::max(scale, std::abs(value));
    return std::max(absolute, relative * scale);
}

Vector3 scaledVector(const Vector3& value, const double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vector3 vectorSum(const Vector3& first, const Vector3& second) {
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

Vector3 vectorDifference(const Vector3& first, const Vector3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

void addVector(Vector3& target, const Vector3& value) {
    target.x += value.x;
    target.y += value.y;
    target.z += value.z;
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({std::abs(value.x), std::abs(value.y),
                     std::abs(value.z)});
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningPressureState& state) {
    return checkedAdd(
        checkedMultiply(
            state.controls.size(),
            sizeof(
                PlanarPressureRegionFragmentOpeningPressureStateControl)),
        checkedAdd(
            checkedMultiply(
                state.walls.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningPressureStateWall)),
            checkedMultiply(
                state.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningPressureStateComponent))));
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings) {
    fingerprint.real(settings.absolutePressureResidualTolerancePascals);
    fingerprint.real(settings.relativePressureResidualTolerance);
    fingerprint.real(settings.absoluteForceResidualToleranceNewtons);
    fingerprint.real(settings.relativeForceResidualTolerance);
    fingerprint.real(settings.absoluteWorkResidualToleranceJoules);
    fingerprint.real(settings.relativeWorkResidualTolerance);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningPressureState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    for (const std::uint64_t value : {
             state.sourceAcceptedStateFingerprint,
             state.sourcePressureOperatorFingerprint,
             state.sourceBasePressureOperatorFingerprint,
             state.sourceOpeningFingerprint,
             state.sourceFragmentFingerprint,
             state.sourceTopologyFingerprint,
             state.sourceVolumeRateFingerprint}) {
        fingerprint.integer(value);
    }
    fingerprintSettings(fingerprint, state.settings);
    fingerprint.boolean(state.staticGeometry);
    fingerprint.boolean(state.usesMovingVolumeRates);
    fingerprint.real(state.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(state.controls.size()));
    for (const auto& control : state.controls) {
        for (const std::size_t value : {
                 control.fragmentIndex,
                 control.baseComponentIndex,
                 control.connectedComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(control.fragmentStableId);
        fingerprint.integer(control.connectedComponentStableId);
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
                 wall.minusBaseComponentIndex,
                 wall.plusBaseComponentIndex,
                 wall.minusConnectedComponentIndex,
                 wall.plusConnectedComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(wall.sourceFaceLinkStableId);
        fingerprint.integer(wall.surfaceStableId);
        fingerprint.enumeration(wall.axis);
        fingerprint.integer(wall.minusRegionStableId);
        fingerprint.integer(wall.plusRegionStableId);
        fingerprint.real(wall.areaSquareMeters);
        fingerprintVector(fingerprint, wall.wrappedCentroidMeters);
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
        fingerprintVector(
            fingerprint, wall.totalPressureImpulseOnSheetNewtonSeconds);
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
    fingerprint.integer(static_cast<std::uint64_t>(state.components.size()));
    for (const auto& component : state.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.baseComponentCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.fragmentCount));
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
    fingerprintVector(
        fingerprint, state.totalPressureImpulseOnSheetNewtonSeconds);
    for (const double value : {
             state.authoredPressureWorkToFluidJoules,
             state.correctionPressureWorkToFluidJoules,
             state.totalPressureWorkToFluidJoules,
             state.pressureWorkSplitResidualJoules,
             state.totalPressureWorkToSheetJoules,
             state.actionReactionWorkResidualJoules,
             state.authoredGeometryPressureWorkJoules,
             state.correctionGeometryPressureWorkJoules,
             state.totalGeometryPressureWorkJoules,
             state.wallGeometryWorkResidualJoules}) {
        fingerprint.real(value);
    }
    fingerprint.boolean(state.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(state.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentOpeningPressureState buildState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits.acceptedStateLimits);
    if (fragments.fragments.size() > limits.maximumControls
        || topology.pressureLayerWallLinkCount > limits.maximumWalls
        || openings.connectedComponents.size() > limits.maximumComponents) {
        throw std::length_error(
            "opening pressure-state entity limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningPressureState result;
    result.sourceAcceptedStateFingerprint = acceptedState.fingerprint;
    result.sourcePressureOperatorFingerprint = pressureOperator.fingerprint;
    result.sourceBasePressureOperatorFingerprint =
        basePressureOperator.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceVolumeRateFingerprint = volumeRates.fingerprint;
    result.settings = settings;
    result.staticGeometry = isStaticGeometry(sweep);
    result.usesMovingVolumeRates = true;
    result.timeStepSeconds =
        acceptedState.settings.projection.timeStepSeconds;
    result.controls.reserve(fragments.fragments.size());
    result.walls.reserve(topology.pressureLayerWallLinkCount);
    result.components.resize(openings.connectedComponents.size());
    result.ownedStorageBytes = checkedAdd(
        checkedMultiply(
            fragments.fragments.size(),
            sizeof(
                PlanarPressureRegionFragmentOpeningPressureStateControl)),
        checkedAdd(
            checkedMultiply(
                topology.pressureLayerWallLinkCount,
                sizeof(
                    PlanarPressureRegionFragmentOpeningPressureStateWall)),
            checkedMultiply(
                openings.connectedComponents.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningPressureStateComponent))));
    result.workingStorageBytes = checkedMultiply(
        fragments.fragments.size(), sizeof(std::size_t));
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening pressure-state storage limit exceeded");
    }

    for (std::size_t index = 0;
         index < openings.connectedComponents.size(); ++index) {
        const auto& source = openings.connectedComponents[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.baseComponentCount = source.baseComponentCount;
        component.fragmentCount = source.fragmentCount;
        component.volumeCubicMeters = source.volumeCubicMeters;
    }
    std::vector<std::size_t> connectedByFragment(
        fragments.fragments.size(), 0);
    double maximumAbsolutePressure = 0.0;
    for (std::size_t index = 0;
         index < fragments.fragments.size(); ++index) {
        const auto& fragment = fragments.fragments[index];
        const auto& rate = volumeRates.fragments[index];
        const std::size_t baseComponentIndex =
            topology.fragments[index].componentIndex;
        if (baseComponentIndex >= openings.baseComponents.size()) {
            throw std::logic_error(
                "opening pressure-state base-component mapping is missing");
        }
        const auto& mapping = openings.baseComponents[baseComponentIndex];
        const std::size_t connectedComponentIndex =
            mapping.connectedComponentIndex;
        if (connectedComponentIndex >= result.components.size()) {
            throw std::logic_error(
                "opening pressure-state connected-component mapping is missing");
        }
        connectedByFragment[index] = connectedComponentIndex;
        const double authored = fragment.pressurePascals;
        const double correction =
            acceptedState.pressureCorrectionPascals[index];
        const double total = authored + correction;
        if (!std::isfinite(total)) {
            throw std::invalid_argument(
                "opening pressure-state control pressure is non-finite");
        }
        result.controls.push_back({
            index,
            fragment.stableId,
            baseComponentIndex,
            connectedComponentIndex,
            result.components[connectedComponentIndex].stableId,
            fragment.regionStableId,
            fragment.volumeCubicMeters,
            authored,
            correction,
            total,
        });
        auto& component = result.components[connectedComponentIndex];
        component.authoredVolumeMeanPressurePascals +=
            authored * fragment.volumeCubicMeters;
        component.correctionVolumeMeanPressurePascals +=
            correction * fragment.volumeCubicMeters;
        component.totalVolumeMeanPressurePascals +=
            total * fragment.volumeCubicMeters;
        const double authoredGeometryWork =
            -result.timeStepSeconds * authored
            * rate.geometryVolumeChangeRateCubicMetersPerSecond;
        const double correctionGeometryWork =
            -result.timeStepSeconds * correction
            * rate.geometryVolumeChangeRateCubicMetersPerSecond;
        component.authoredGeometryPressureWorkJoules +=
            authoredGeometryWork;
        component.correctionGeometryPressureWorkJoules +=
            correctionGeometryWork;
        component.totalGeometryPressureWorkJoules +=
            authoredGeometryWork + correctionGeometryWork;
        result.authoredGeometryPressureWorkJoules +=
            authoredGeometryWork;
        result.correctionGeometryPressureWorkJoules +=
            correctionGeometryWork;
        result.totalGeometryPressureWorkJoules +=
            authoredGeometryWork + correctionGeometryWork;
        maximumAbsolutePressure = std::max({
            maximumAbsolutePressure,
            std::abs(authored), std::abs(correction), std::abs(total)});
    }
    const double pressureTolerance = scaledTolerance(
        settings.absolutePressureResidualTolerancePascals,
        settings.relativePressureResidualTolerance,
        {maximumAbsolutePressure});
    for (auto& component : result.components) {
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
                "opening pressure-state connected gauge failed");
        }
    }

    double maximumForceNewtons = 0.0;
    double maximumWorkJoules = 0.0;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        const auto& minus = result.controls[link.minusFragmentIndex];
        const auto& plus = result.controls[link.plusFragmentIndex];
        if (connectedByFragment[link.minusFragmentIndex]
                != minus.connectedComponentIndex
            || connectedByFragment[link.plusFragmentIndex]
                != plus.connectedComponentIndex) {
            throw std::logic_error(
                "opening pressure-state fragment mapping is inconsistent");
        }
        const double minusWallVelocity = volumeRates
            .fragments[link.minusFragmentIndex]
            .upperBoundaryVelocityMetersPerSecond;
        const double plusWallVelocity = volumeRates
            .fragments[link.plusFragmentIndex]
            .lowerBoundaryVelocityMetersPerSecond;
        if (minusWallVelocity != plusWallVelocity) {
            throw std::invalid_argument(
                "opening pressure-state material wall velocity disagrees");
        }
        const double authoredJump = link.pressureJumpPascals;
        const double correctionJump = plus.correctionPressurePascals
            - minus.correctionPressurePascals;
        const double totalJump = plus.totalPressurePascals
            - minus.totalPressurePascals;
        const double pressureResidual =
            totalJump - (authoredJump + correctionJump);
        const Vector3 authoredForce = scaledVector(
            link.unitNormalMinusToPlus,
            -authoredJump * link.areaSquareMeters);
        const Vector3 correctionForce = scaledVector(
            link.unitNormalMinusToPlus,
            -correctionJump * link.areaSquareMeters);
        const Vector3 totalForce = scaledVector(
            link.unitNormalMinusToPlus,
            -totalJump * link.areaSquareMeters);
        const Vector3 forceResidual = vectorDifference(
            totalForce, vectorSum(authoredForce, correctionForce));
        const Vector3 totalImpulse = scaledVector(
            totalForce, result.timeStepSeconds);
        const double authoredWork = authoredJump * link.areaSquareMeters
            * minusWallVelocity * result.timeStepSeconds;
        const double correctionWork = correctionJump * link.areaSquareMeters
            * minusWallVelocity * result.timeStepSeconds;
        const double totalWork = totalJump * link.areaSquareMeters
            * minusWallVelocity * result.timeStepSeconds;
        const double workSplitResidual =
            totalWork - (authoredWork + correctionWork);
        const double sheetWork = -totalWork;
        const double actionReactionWorkResidual = totalWork + sheetWork;
        const double forceTolerance = scaledTolerance(
            settings.absoluteForceResidualToleranceNewtons,
            settings.relativeForceResidualTolerance,
            {maximumAbsoluteComponent(authoredForce),
             maximumAbsoluteComponent(correctionForce),
             maximumAbsoluteComponent(totalForce)});
        const double workTolerance = scaledTolerance(
            settings.absoluteWorkResidualToleranceJoules,
            settings.relativeWorkResidualTolerance,
            {authoredWork, correctionWork, totalWork});
        if (!std::isfinite(pressureResidual)
            || std::abs(pressureResidual) > pressureTolerance
            || !finiteVector(forceResidual)
            || maximumAbsoluteComponent(forceResidual) > forceTolerance
            || !std::isfinite(workSplitResidual)
            || std::abs(workSplitResidual) > workTolerance
            || !std::isfinite(actionReactionWorkResidual)
            || std::abs(actionReactionWorkResidual) > workTolerance) {
            throw std::invalid_argument(
                "opening pressure-state wall closure failed");
        }
        result.walls.push_back({
            result.walls.size(),
            link.linkIndex,
            link.stableId,
            link.surfaceStableId,
            link.axis,
            link.minusFragmentIndex,
            link.plusFragmentIndex,
            minus.baseComponentIndex,
            plus.baseComponentIndex,
            minus.connectedComponentIndex,
            plus.connectedComponentIndex,
            link.minusRegionStableId,
            link.plusRegionStableId,
            link.areaSquareMeters,
            link.wrappedCentroidMeters,
            link.unitNormalMinusToPlus,
            minus.totalPressurePascals,
            plus.totalPressurePascals,
            authoredJump,
            correctionJump,
            totalJump,
            pressureResidual,
            authoredForce,
            correctionForce,
            totalForce,
            forceResidual,
            totalImpulse,
            minusWallVelocity,
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
        addVector(result.authoredPressureForceOnSheetNewtons, authoredForce);
        addVector(
            result.correctionPressureForceOnSheetNewtons,
            correctionForce);
        addVector(result.totalPressureForceOnSheetNewtons, totalForce);
        addVector(
            result.totalPressureImpulseOnSheetNewtonSeconds,
            totalImpulse);
        result.authoredPressureWorkToFluidJoules += authoredWork;
        result.correctionPressureWorkToFluidJoules += correctionWork;
        result.totalPressureWorkToFluidJoules += totalWork;
        result.components[minus.connectedComponentIndex]
            .totalWallPressureWorkToFluidJoules +=
            -minus.totalPressurePascals * link.areaSquareMeters
            * minusWallVelocity * result.timeStepSeconds;
        result.components[plus.connectedComponentIndex]
            .totalWallPressureWorkToFluidJoules +=
            plus.totalPressurePascals * link.areaSquareMeters
            * minusWallVelocity * result.timeStepSeconds;
        maximumForceNewtons = std::max({
            maximumForceNewtons,
            maximumAbsoluteComponent(authoredForce),
            maximumAbsoluteComponent(correctionForce),
            maximumAbsoluteComponent(totalForce)});
        maximumWorkJoules = std::max({
            maximumWorkJoules,
            std::abs(authoredWork),
            std::abs(correctionWork),
            std::abs(totalWork)});
    }
    if (result.walls.size() != topology.pressureLayerWallLinkCount) {
        throw std::logic_error(
            "opening pressure-state wall coverage is incomplete");
    }

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
            result.correctionGeometryPressureWorkJoules
            - acceptedState.geometryPressureWorkJoules)});
    const double globalForceTolerance = scaledTolerance(
        settings.absoluteForceResidualToleranceNewtons,
        settings.relativeForceResidualTolerance,
        {maximumForceNewtons});
    const double globalWorkTolerance = scaledTolerance(
        settings.absoluteWorkResidualToleranceJoules,
        settings.relativeWorkResidualTolerance,
        {maximumWorkJoules,
         result.authoredPressureWorkToFluidJoules,
         result.correctionPressureWorkToFluidJoules,
         result.totalPressureWorkToFluidJoules,
         result.totalGeometryPressureWorkJoules,
         acceptedState.geometryPressureWorkJoules});
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
               result.correctionGeometryPressureWorkJoules
               - acceptedState.geometryPressureWorkJoules)
            > globalWorkTolerance) {
        throw std::invalid_argument(
            "opening pressure-state global closure failed");
    }
    for (auto& component : result.components) {
        component.wallGeometryWorkResidualJoules =
            component.totalWallPressureWorkToFluidJoules
            - component.totalGeometryPressureWorkJoules;
        result.maximumAbsoluteWorkResidualJoules = std::max(
            result.maximumAbsoluteWorkResidualJoules,
            std::abs(component.wallGeometryWorkResidualJoules));
        const double componentWorkTolerance = scaledTolerance(
            settings.absoluteWorkResidualToleranceJoules,
            settings.relativeWorkResidualTolerance,
            {component.totalWallPressureWorkToFluidJoules,
             component.totalGeometryPressureWorkJoules,
             maximumWorkJoules});
        if (!std::isfinite(component.wallGeometryWorkResidualJoules)
            || std::abs(component.wallGeometryWorkResidualJoules)
                > componentWorkTolerance) {
            throw std::invalid_argument(
                "opening pressure-state component work closure failed");
        }
    }
    if (!finiteVector(result.authoredPressureForceOnSheetNewtons)
        || !finiteVector(result.correctionPressureForceOnSheetNewtons)
        || !finiteVector(result.totalPressureForceOnSheetNewtons)
        || !finiteVector(result.totalPressureImpulseOnSheetNewtonSeconds)
        || result.ownedStorageBytes != ownedStorageBytes(result)) {
        throw std::logic_error(
            "opening pressure-state aggregate ledger is invalid");
    }
    result.accepted = true;
    result.fingerprint = stateFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningPressureState
composePlanarPressureRegionFragmentOpeningPressureState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings,
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits) {
    auto result = buildState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, settings, limits);
    validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(result);
    return result;
}

void validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureState& state) {
    validateSettings(state.settings);
    if (state.version
            != planarPressureRegionFragmentOpeningPressureStateVersion
        || state.fingerprint == 0 || !state.accepted
        || state.sourceAcceptedStateFingerprint == 0
        || state.sourcePressureOperatorFingerprint == 0
        || state.sourceBasePressureOperatorFingerprint == 0
        || state.sourceOpeningFingerprint == 0
        || state.sourceFragmentFingerprint == 0
        || state.sourceTopologyFingerprint == 0
        || state.sourceVolumeRateFingerprint == 0
        || !std::isfinite(state.timeStepSeconds)
        || !(state.timeStepSeconds > 0.0)
        || !finiteVector(state.authoredPressureForceOnSheetNewtons)
        || !finiteVector(state.correctionPressureForceOnSheetNewtons)
        || !finiteVector(state.totalPressureForceOnSheetNewtons)
        || !finiteVector(state.pressureForceSplitResidualNewtons)
        || !finiteVector(state.totalPressureImpulseOnSheetNewtonSeconds)
        || !std::ranges::all_of(
            std::initializer_list<double>{
                state.maximumAbsoluteCorrectionGaugePascals,
                state.maximumAbsolutePressureSplitResidualPascals,
                state.maximumAbsoluteForceSplitResidualNewtons,
                state.maximumAbsoluteWorkResidualJoules,
                state.authoredPressureWorkToFluidJoules,
                state.correctionPressureWorkToFluidJoules,
                state.totalPressureWorkToFluidJoules,
                state.pressureWorkSplitResidualJoules,
                state.totalPressureWorkToSheetJoules,
                state.actionReactionWorkResidualJoules,
                state.authoredGeometryPressureWorkJoules,
                state.correctionGeometryPressureWorkJoules,
                state.totalGeometryPressureWorkJoules,
                state.wallGeometryWorkResidualJoules},
            [](const double value) { return std::isfinite(value); })
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening pressure-state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningPressureState(
    const PlanarPressureRegionFragmentOpeningPressureState& state,
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(state);
    if (state != buildState(
                     acceptedState, pressureOperator, basePressureOperator,
                     grid, sweep, fragments, topology, volumeRates,
                     openingDefinitions, openings, resistanceDefinitions,
                     state.settings, limits)) {
        throw std::invalid_argument(
            "opening pressure state is corrupted");
    }
}

} // namespace simwing::fsi::fluid
