#include "scene_fluid_regional_opening_momentum_wall_input.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <ranges>
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
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
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

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double norm(const fluid::Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

double tolerance(const double reference) {
    return 128.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, std::abs(reference));
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening momentum wall-input storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "regional opening momentum wall-input storage overflows");
    }
    return first * second;
}

std::size_t ownedStorageBytes(const std::size_t controls,
                              const std::size_t samples) {
    return checkedAdd(
        checkedMultiply(
            controls, sizeof(SceneFluidRegionWallControlVolume)),
        checkedMultiply(samples, sizeof(SceneFluidRegionWallSample)));
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits) {
    if (limits.maximumControlVolumes == 0 || limits.maximumSamples == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening momentum wall-input limits are invalid");
    }
}

void fingerprintVector(Fingerprint& fingerprint,
                       const fluid::Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t inputFingerprint(
    const SceneFluidRegionalOpeningMomentumWallInput& input) {
    Fingerprint fingerprint;
    for (const std::uint64_t value : {
             static_cast<std::uint64_t>(input.version),
             input.sourceTransportFingerprint,
             input.sourceAcceptedStateFingerprint,
             input.sourceCurrentFlowStateFingerprint,
             input.sourceTransportMetricFingerprint,
             input.sourcePressureOperatorFingerprint,
             input.sourceBasePressureOperatorFingerprint,
             input.sourceOpeningFingerprint,
             input.sourceFragmentFingerprint,
             input.sourceTopologyFingerprint,
             input.sourceVolumeRateFingerprint,
             input.sourceLoadStateFingerprint,
             input.sourceSamplingFingerprint,
             input.quadratureFingerprint,
             input.surfaceDefinitionFingerprint,
             input.surfaceStateFingerprint,
             input.structureDefinitionFingerprint,
             input.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(input.simulationTimeSeconds);
    fingerprint.real(input.densityKgPerCubicMeter);
    fingerprint.real(input.timeStepSeconds);
    const auto& pressure = input.settings.pressureState;
    fingerprint.real(pressure.absolutePressureResidualTolerancePascals);
    fingerprint.real(pressure.relativePressureResidualTolerance);
    fingerprint.real(pressure.absoluteForceResidualToleranceNewtons);
    fingerprint.real(pressure.relativeForceResidualTolerance);
    fingerprint.real(pressure.absoluteWorkResidualToleranceJoules);
    fingerprint.real(pressure.relativeWorkResidualTolerance);
    fingerprint.integer(
        static_cast<std::uint64_t>(input.activeControlVolumeCount));
    fingerprint.real(input.wallSampleAreaSquareMeters);
    fingerprint.real(input.controlIncidentWallAreaSquareMeters);
    fingerprint.real(input.maximumIncidentWallAreaSquareMeters);
    fingerprint.integer(
        static_cast<std::uint64_t>(input.controlVolumes.size()));
    for (const auto& control : input.controlVolumes) {
        fingerprint.integer(
            static_cast<std::uint64_t>(control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.real(control.volumeCubicMeters);
        fingerprint.real(control.incidentWallAreaSquareMeters);
        fingerprintVector(fingerprint, control.velocityMetersPerSecond);
        fingerprintVector(fingerprint, control.momentumKilogramMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(input.samples.size()));
    for (const auto& sample : input.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(sample.stableId);
        fingerprint.integer(sample.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.negativeSideControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.positiveSideControlVolumeIndex));
        fingerprint.real(sample.areaSquareMeters);
        fingerprintVector(
            fingerprint, sample.unitNormalNegativeToPositive);
        fingerprintVector(
            fingerprint, sample.wallVelocityMetersPerSecond);
        fingerprintVector(
            fingerprint,
            sample.negativeSideFluidImpulseKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint,
            sample.positiveSideFluidImpulseKilogramMetersPerSecond);
        fingerprint.integer(sample.structureTraction.stableId);
        fingerprint.real(sample.structureTraction.tractionPascals.x);
        fingerprint.real(sample.structureTraction.tractionPascals.y);
        fingerprint.real(sample.structureTraction.tractionPascals.z);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(input.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(input.workingStorageBytes));
    return fingerprint.value();
}

SceneFluidRegionalOpeningMomentumWallInput buildInput(
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        currentAcceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits) {
    validateLimits(limits);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        transport);
    auto currentFlow =
        fluid::capturePlanarPressureRegionFragmentOpeningVelocityState(
            currentAcceptedState, pressureOperator, basePressureOperator,
            grid, sweep, fragments, topology, volumeRates,
            openingDefinitions, openings, resistanceDefinitions, baseMetric,
            metric, limits.flowState);
    if (!transport.diagnostics.accepted
        || transport.targetFlowStateFingerprint != currentFlow.fingerprint
        || transport.targetMetricFingerprint != metric.fingerprint
        || transport.targetVolumeRateFingerprint != volumeRates.fingerprint
        || transport.densityKgPerCubicMeter
            != currentFlow.densityKgPerCubicMeter
        || transport.timeStepSeconds != volumeRates.durationSeconds
        || transport.controls.size() != fragments.fragments.size()) {
        throw std::invalid_argument(
            "regional opening momentum wall-input transport epoch is invalid");
    }

    auto loadState =
        fluid::composePlanarPressureRegionFragmentOpeningLoadState(
            currentAcceptedState, pressureOperator, basePressureOperator,
            grid, sweep, fragments, topology, volumeRates,
            openingDefinitions, openings, resistanceDefinitions,
            settings.pressureState, limits.loadState);
    auto regionalSamples = sampleSceneFluidRegionalOpeningPressure(
        loadState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, surface, surfaceState, quadrature,
        limits.loadState, limits.sampling);
    const auto kinematics = sampleSceneFluidQuadratureKinematics(
        surface, surfaceState, quadrature);
    if (regionalSamples.bindings.size() != quadrature.points.size()
        || kinematics.size() != quadrature.points.size()) {
        throw std::logic_error(
            "regional opening momentum wall-input sampling size changed");
    }

    const std::size_t controlCount = transport.controls.size();
    const std::size_t sampleCount = regionalSamples.bindings.size();
    const std::size_t outputOwned = ownedStorageBytes(
        controlCount, sampleCount);
    std::size_t working = outputOwned;
    working = checkedAdd(working, currentFlow.ownedStorageBytes);
    working = checkedAdd(working, currentFlow.workingStorageBytes);
    working = checkedAdd(working, loadState.ownedStorageBytes);
    working = checkedAdd(working, regionalSamples.ownedStorageBytes);
    working = checkedAdd(working, regionalSamples.workingStorageBytes);
    if (controlCount > limits.maximumControlVolumes
        || sampleCount > limits.maximumSamples
        || outputOwned > limits.maximumOwnedBytes
        || working > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening momentum wall-input limit exceeded");
    }

    SceneFluidRegionalOpeningMomentumWallInput result;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.sourceAcceptedStateFingerprint = currentAcceptedState.fingerprint;
    result.sourceCurrentFlowStateFingerprint = currentFlow.fingerprint;
    result.sourceTransportMetricFingerprint = metric.fingerprint;
    result.sourcePressureOperatorFingerprint = pressureOperator.fingerprint;
    result.sourceBasePressureOperatorFingerprint =
        basePressureOperator.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceVolumeRateFingerprint = volumeRates.fingerprint;
    result.sourceLoadStateFingerprint = loadState.fingerprint;
    result.sourceSamplingFingerprint = regionalSamples.fingerprint;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = surfaceState.fingerprint;
    result.structureDefinitionFingerprint =
        surfaceState.structureDefinitionFingerprint;
    result.acceptedStepCount = surfaceState.acceptedStepCount;
    result.simulationTimeSeconds = surfaceState.simulationTimeSeconds;
    result.densityKgPerCubicMeter = transport.densityKgPerCubicMeter;
    result.timeStepSeconds = transport.timeStepSeconds;
    result.settings = settings;
    result.controlVolumes.reserve(controlCount);
    for (std::size_t index = 0; index < controlCount; ++index) {
        const auto& source = transport.controls[index];
        const auto& fragment = fragments.fragments[index];
        const auto& pressure = loadState.pressure.controls[index];
        if (source.fragmentIndex != index || source.stableId == 0
            || source.stableId != fragment.stableId
            || source.stableId != pressure.fragmentStableId
            || source.regionStableId != fragment.regionStableId
            || source.regionStableId != pressure.regionStableId
            || source.volumeCubicMeters != fragment.volumeCubicMeters
            || source.volumeCubicMeters != pressure.volumeCubicMeters) {
            throw std::invalid_argument(
                "regional opening momentum wall-input control is foreign");
        }
        result.controlVolumes.push_back({
            index,
            source.stableId,
            source.volumeCubicMeters,
            0.0,
            source.velocityMetersPerSecond,
            source.momentumKilogramMetersPerSecond,
        });
    }

    result.samples.reserve(sampleCount);
    for (std::size_t index = 0; index < sampleCount; ++index) {
        const auto& binding = regionalSamples.bindings[index];
        const auto& point = quadrature.points[index];
        const auto& motion = kinematics[index];
        if (binding.sampleIndex != index || binding.stableId != point.stableId
            || motion.stableId != point.stableId
            || binding.sourcePressureWallIndex
                >= loadState.pressure.walls.size()) {
            throw std::invalid_argument(
                "regional opening momentum wall-input sample is foreign");
        }
        const auto& wall =
            loadState.pressure.walls[binding.sourcePressureWallIndex];
        if (wall.sourceFaceLinkIndex != binding.sourceFaceLinkIndex
            || wall.sourceFaceLinkStableId
                != binding.sourceFaceLinkStableId
            || wall.minusRegionStableId != binding.negativeSideRegionId
            || wall.plusRegionStableId != binding.positiveSideRegionId
            || wall.minusFragmentIndex >= result.controlVolumes.size()
            || wall.plusFragmentIndex >= result.controlVolumes.size()
            || transport.controls[wall.minusFragmentIndex].regionStableId
                != wall.minusRegionStableId
            || transport.controls[wall.plusFragmentIndex].regionStableId
                != wall.plusRegionStableId) {
            throw std::invalid_argument(
                "regional opening momentum wall-input wall ownership is invalid");
        }
        SceneFluidRegionWallSample sample;
        sample.sampleIndex = index;
        sample.stableId = point.stableId;
        sample.triangleId = point.triangleId;
        sample.negativeSideControlVolumeIndex = wall.minusFragmentIndex;
        sample.positiveSideControlVolumeIndex = wall.plusFragmentIndex;
        sample.areaSquareMeters = point.areaSquareMeters;
        sample.unitNormalNegativeToPositive = wall.unitNormalMinusToPlus;
        sample.wallVelocityMetersPerSecond = {
            motion.velocityMetersPerSecond.x,
            motion.velocityMetersPerSecond.y,
            motion.velocityMetersPerSecond.z,
        };
        sample.structureTraction.stableId = point.stableId;
        result.samples.push_back(sample);
        result.controlVolumes[wall.minusFragmentIndex]
            .incidentWallAreaSquareMeters += point.areaSquareMeters;
        result.controlVolumes[wall.plusFragmentIndex]
            .incidentWallAreaSquareMeters += point.areaSquareMeters;
    }

    for (const auto& control : result.controlVolumes) {
        if (control.incidentWallAreaSquareMeters > 0.0) {
            ++result.activeControlVolumeCount;
        }
        result.controlIncidentWallAreaSquareMeters +=
            control.incidentWallAreaSquareMeters;
        result.maximumIncidentWallAreaSquareMeters = std::max(
            result.maximumIncidentWallAreaSquareMeters,
            control.incidentWallAreaSquareMeters);
    }
    for (const auto& sample : result.samples) {
        result.wallSampleAreaSquareMeters += sample.areaSquareMeters;
    }
    if (std::abs(
            result.controlIncidentWallAreaSquareMeters
            - 2.0 * result.wallSampleAreaSquareMeters)
        > tolerance(result.controlIncidentWallAreaSquareMeters)) {
        throw std::logic_error(
            "regional opening momentum wall-input incident area does not close");
    }
    result.ownedStorageBytes = outputOwned;
    result.workingStorageBytes = working;
    result.fingerprint = inputFingerprint(result);
    validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(result);
    return result;
}

} // namespace

SceneFluidRegionalOpeningMomentumWallInput
captureSceneFluidRegionalOpeningMomentumWallInput(
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        currentAcceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits) {
    return buildInput(
        transport, currentAcceptedState, pressureOperator,
        basePressureOperator, grid, sweep, fragments, topology, volumeRates,
        openingDefinitions, openings, resistanceDefinitions, baseMetric,
        metric, surface, surfaceState, quadrature, settings, limits);
}

void validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(
    const SceneFluidRegionalOpeningMomentumWallInput& input) {
    bool controlsValid = true;
    std::vector<double> incidentArea(input.controlVolumes.size(), 0.0);
    std::size_t activeControls = 0;
    double incidentAreaSum = 0.0;
    double maximumIncidentArea = 0.0;
    for (std::size_t index = 0; index < input.controlVolumes.size(); ++index) {
        const auto& control = input.controlVolumes[index];
        controlsValid = controlsValid
            && control.controlVolumeIndex == index && control.stableId != 0
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && std::isfinite(control.incidentWallAreaSquareMeters)
            && control.incidentWallAreaSquareMeters >= 0.0
            && finite(control.velocityMetersPerSecond)
            && finite(control.momentumKilogramMetersPerSecond);
        if (control.incidentWallAreaSquareMeters > 0.0) {
            ++activeControls;
        }
        incidentAreaSum += control.incidentWallAreaSquareMeters;
        maximumIncidentArea = std::max(
            maximumIncidentArea, control.incidentWallAreaSquareMeters);
    }
    bool samplesValid = true;
    double sampleArea = 0.0;
    for (std::size_t index = 0; index < input.samples.size(); ++index) {
        const auto& sample = input.samples[index];
        samplesValid = samplesValid
            && sample.sampleIndex == index && sample.stableId != 0
            && sample.triangleId != invalidStableId
            && sample.negativeSideControlVolumeIndex
                < input.controlVolumes.size()
            && sample.positiveSideControlVolumeIndex
                < input.controlVolumes.size()
            && std::isfinite(sample.areaSquareMeters)
            && sample.areaSquareMeters > 0.0
            && finite(sample.unitNormalNegativeToPositive)
            && std::abs(norm(sample.unitNormalNegativeToPositive) - 1.0)
                <= 1.0e-12
            && finite(sample.wallVelocityMetersPerSecond)
            && sample.negativeSideFluidImpulseKilogramMetersPerSecond
                == fluid::Vector3{}
            && sample.positiveSideFluidImpulseKilogramMetersPerSecond
                == fluid::Vector3{}
            && sample.structureTraction.stableId == sample.stableId
            && sample.structureTraction.tractionPascals
                == StructureVector3{};
        if (sample.negativeSideControlVolumeIndex
                < incidentArea.size()
            && sample.positiveSideControlVolumeIndex
                < incidentArea.size()) {
            incidentArea[sample.negativeSideControlVolumeIndex] +=
                sample.areaSquareMeters;
            incidentArea[sample.positiveSideControlVolumeIndex] +=
                sample.areaSquareMeters;
        }
        sampleArea += sample.areaSquareMeters;
    }
    bool incidenceValid = incidentArea.size() == input.controlVolumes.size();
    for (std::size_t index = 0;
         incidenceValid && index < incidentArea.size(); ++index) {
        incidenceValid = incidentArea[index]
            == input.controlVolumes[index].incidentWallAreaSquareMeters;
    }
    const auto& pressure = input.settings.pressureState;
    const bool settingsValid =
        std::isfinite(pressure.absolutePressureResidualTolerancePascals)
        && pressure.absolutePressureResidualTolerancePascals >= 0.0
        && std::isfinite(pressure.relativePressureResidualTolerance)
        && pressure.relativePressureResidualTolerance >= 0.0
        && std::isfinite(pressure.absoluteForceResidualToleranceNewtons)
        && pressure.absoluteForceResidualToleranceNewtons >= 0.0
        && std::isfinite(pressure.relativeForceResidualTolerance)
        && pressure.relativeForceResidualTolerance >= 0.0
        && std::isfinite(pressure.absoluteWorkResidualToleranceJoules)
        && pressure.absoluteWorkResidualToleranceJoules >= 0.0
        && std::isfinite(pressure.relativeWorkResidualTolerance)
        && pressure.relativeWorkResidualTolerance >= 0.0;
    if (input.version
            != sceneFluidRegionalOpeningMomentumWallInputVersion
        || input.fingerprint == 0
        || input.sourceTransportFingerprint == 0
        || input.sourceAcceptedStateFingerprint == 0
        || input.sourceCurrentFlowStateFingerprint == 0
        || input.sourceTransportMetricFingerprint == 0
        || input.sourcePressureOperatorFingerprint == 0
        || input.sourceBasePressureOperatorFingerprint == 0
        || input.sourceOpeningFingerprint == 0
        || input.sourceFragmentFingerprint == 0
        || input.sourceTopologyFingerprint == 0
        || input.sourceVolumeRateFingerprint == 0
        || input.sourceLoadStateFingerprint == 0
        || input.sourceSamplingFingerprint == 0
        || input.quadratureFingerprint == 0
        || input.surfaceDefinitionFingerprint == 0
        || input.surfaceStateFingerprint == 0
        || input.structureDefinitionFingerprint == 0
        || !std::isfinite(input.simulationTimeSeconds)
        || !(input.densityKgPerCubicMeter > 0.0)
        || !(input.timeStepSeconds > 0.0)
        || !settingsValid || input.controlVolumes.empty()
        || !controlsValid || !samplesValid || !incidenceValid
        || input.activeControlVolumeCount != activeControls
        || input.wallSampleAreaSquareMeters != sampleArea
        || input.controlIncidentWallAreaSquareMeters != incidentAreaSum
        || input.maximumIncidentWallAreaSquareMeters != maximumIncidentArea
        || std::abs(incidentAreaSum - 2.0 * sampleArea)
            > tolerance(incidentAreaSum)
        || input.ownedStorageBytes != ownedStorageBytes(
            input.controlVolumes.size(), input.samples.size())
        || input.workingStorageBytes < input.ownedStorageBytes
        || input.fingerprint != inputFingerprint(input)) {
        throw std::invalid_argument(
            "regional opening momentum wall-input integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallInput(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        currentAcceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(input);
    const auto expected = buildInput(
        transport, currentAcceptedState, pressureOperator,
        basePressureOperator, grid, sweep, fragments, topology, volumeRates,
        openingDefinitions, openings, resistanceDefinitions, baseMetric,
        metric, surface, surfaceState, quadrature, settings, limits);
    if (input != expected) {
        throw std::invalid_argument(
            "regional opening momentum wall-input sources are foreign");
    }
}

} // namespace simwing::fsi
