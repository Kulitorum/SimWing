#include "scene_fluid_pressure_epoch_transition.h"

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

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "scene fluid pressure epoch-transition storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const SceneFluidPressureEpochTransition& transition) {
    return checkedAdd(
        transition.currentPressureEpoch.ownedStorageBytes,
        transition.topologyTransition.ownedStorageBytes);
}

bool hasStableControlVolumeTopology(
    const SceneFluidPressureTopologyTransition& transition) {
    return transition.appearedControlVolumeCount == 0
        && transition.disappearedControlVolumeCount == 0
        && transition.previousControlVolumeCount
            == transition.currentControlVolumeCount
        && transition.retainedControlVolumeCount
            == transition.currentControlVolumeCount;
}

std::uint64_t transitionFingerprint(
    const SceneFluidPressureEpochTransition& transition) {
    Fingerprint fingerprint;
    fingerprint.integer(transition.version);
    for (const std::uint64_t value : {
             transition.previousPressureEpochFingerprint,
             transition.acceptedCurrentGridEpochFingerprint,
             transition.currentPressureEpochFingerprint,
             transition.topologyTransitionFingerprint,
             transition.surfaceDefinitionFingerprint,
             transition.structureDefinitionFingerprint,
             transition.previousSurfaceStateFingerprint,
             transition.currentSurfaceStateFingerprint,
             transition.previousAcceptedStepCount,
             transition.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(transition.previousSimulationTimeSeconds);
    fingerprint.real(transition.currentSimulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.controlVolumeTopologyStable));
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.ownedStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const SceneFluidPressureEpochTransitionLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "scene fluid pressure epoch-transition limits are invalid");
    }
}

void validateSources(
    const SceneFluidPressureEpoch& previousPressureEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch) {
    validateSceneFluidPressureEpoch(
        previousPressureEpoch, surface, previousSurfaceState, grid,
        transfer, connectivity);
    validateSceneFluidSurfaceState(surface, currentSurfaceState);
    validateSceneFluidGridEpoch(
        acceptedCurrentGridEpoch, surface, currentSurfaceState, grid,
        transfer);
    if (previousSurfaceState.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentSurfaceState.acceptedStepCount
            != previousSurfaceState.acceptedStepCount + 1
        || !(currentSurfaceState.simulationTimeSeconds
            > previousSurfaceState.simulationTimeSeconds)
        || previousSurfaceState.definitionFingerprint
            != currentSurfaceState.definitionFingerprint
        || previousSurfaceState.structureDefinitionFingerprint
            != currentSurfaceState.structureDefinitionFingerprint
        || previousPressureEpoch.surfaceStateFingerprint
            != previousSurfaceState.fingerprint
        || acceptedCurrentGridEpoch.quadrature.surfaceStateFingerprint
            != currentSurfaceState.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure epoch-transition source identity is invalid");
    }
}

SceneFluidPressureEpochTransition buildTransition(
    const SceneFluidPressureEpoch& previousPressureEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidPressureEpochSettings& settings,
    const SceneFluidPressureEpochTransitionLimits& limits) {
    validateLimits(limits);
    validateSources(
        previousPressureEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch);

    SceneFluidPressureEpochTransition result;
    result.previousPressureEpochFingerprint =
        previousPressureEpoch.fingerprint;
    result.acceptedCurrentGridEpochFingerprint =
        acceptedCurrentGridEpoch.fingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.structureDefinitionFingerprint =
        currentSurfaceState.structureDefinitionFingerprint;
    result.previousSurfaceStateFingerprint =
        previousSurfaceState.fingerprint;
    result.currentSurfaceStateFingerprint = currentSurfaceState.fingerprint;
    result.previousAcceptedStepCount =
        previousSurfaceState.acceptedStepCount;
    result.currentAcceptedStepCount = currentSurfaceState.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousSurfaceState.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentSurfaceState.simulationTimeSeconds;
    result.currentPressureEpoch = buildSceneFluidPressureEpoch(
        surface, currentSurfaceState, grid, transfer, connectivity,
        settings, limits.pressureEpoch);
    if (result.currentPressureEpoch.gridEpoch
            != acceptedCurrentGridEpoch) {
        throw std::invalid_argument(
            "scene fluid pressure epoch-transition rebuilt a foreign current grid epoch");
    }
    result.currentPressureEpochFingerprint =
        result.currentPressureEpoch.fingerprint;
    result.topologyTransition = buildSceneFluidPressureTopologyTransition(
        previousPressureEpoch.pressureControlVolumes,
        previousPressureEpoch.pressureFaceLinks,
        result.currentPressureEpoch.pressureControlVolumes,
        result.currentPressureEpoch.pressureFaceLinks,
        limits.topologyTransition);
    result.topologyTransitionFingerprint =
        result.topologyTransition.fingerprint;
    result.controlVolumeTopologyStable =
        hasStableControlVolumeTopology(result.topologyTransition);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure epoch transition exceeds its aggregate limit");
    }
    result.fingerprint = transitionFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureEpochTransition buildSceneFluidPressureEpochTransition(
    const SceneFluidPressureEpoch& previousPressureEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidPressureEpochSettings& settings,
    const SceneFluidPressureEpochTransitionLimits& limits) {
    auto result = buildTransition(
        previousPressureEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    validateSceneFluidPressureEpochTransition(
        result, previousPressureEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    return result;
}

void validateSceneFluidPressureEpochTransitionIntegrity(
    const SceneFluidPressureEpochTransition& transition) {
    validateSceneFluidPressureControlVolumeIntegrity(
        transition.currentPressureEpoch.pressureControlVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(
        transition.currentPressureEpoch.pressureFaceLinks);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        transition.topologyTransition);
    if (transition.version != sceneFluidPressureEpochTransitionVersion
        || transition.fingerprint == 0
        || transition.previousPressureEpochFingerprint == 0
        || transition.acceptedCurrentGridEpochFingerprint == 0
        || transition.currentPressureEpochFingerprint == 0
        || transition.topologyTransitionFingerprint == 0
        || transition.surfaceDefinitionFingerprint == 0
        || transition.structureDefinitionFingerprint == 0
        || transition.previousSurfaceStateFingerprint == 0
        || transition.currentSurfaceStateFingerprint == 0
        || transition.previousAcceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || transition.currentAcceptedStepCount
            != transition.previousAcceptedStepCount + 1
        || !std::isfinite(transition.previousSimulationTimeSeconds)
        || !std::isfinite(transition.currentSimulationTimeSeconds)
        || transition.previousSimulationTimeSeconds < 0.0
        || !(transition.currentSimulationTimeSeconds
            > transition.previousSimulationTimeSeconds)
        || transition.currentPressureEpoch.fingerprint
            != transition.currentPressureEpochFingerprint
        || transition.topologyTransition.fingerprint
            != transition.topologyTransitionFingerprint
        || transition.currentPressureEpoch.surfaceDefinitionFingerprint
            != transition.surfaceDefinitionFingerprint
        || transition.currentPressureEpoch.structureDefinitionFingerprint
            != transition.structureDefinitionFingerprint
        || transition.currentPressureEpoch.surfaceStateFingerprint
            != transition.currentSurfaceStateFingerprint
        || transition.currentPressureEpoch.acceptedStepCount
            != transition.currentAcceptedStepCount
        || transition.currentPressureEpoch.simulationTimeSeconds
            != transition.currentSimulationTimeSeconds
        || transition.topologyTransition.previousAcceptedStepCount
            != transition.previousAcceptedStepCount
        || transition.topologyTransition.currentAcceptedStepCount
            != transition.currentAcceptedStepCount
        || transition.controlVolumeTopologyStable
            != hasStableControlVolumeTopology(
                transition.topologyTransition)
        || transition.ownedStorageBytes != ownedStorageBytes(transition)
        || transition.fingerprint != transitionFingerprint(transition)) {
        throw std::invalid_argument(
            "scene fluid pressure epoch-transition integrity is invalid");
    }
}

void validateSceneFluidPressureEpochTransition(
    const SceneFluidPressureEpochTransition& transition,
    const SceneFluidPressureEpoch& previousPressureEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidPressureEpochSettings& settings,
    const SceneFluidPressureEpochTransitionLimits& limits) {
    validateLimits(limits);
    validateSceneFluidPressureEpochTransitionIntegrity(transition);
    validateSources(
        previousPressureEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch);
    validateSceneFluidPressureEpoch(
        transition.currentPressureEpoch, surface, currentSurfaceState,
        grid, transfer, connectivity);
    validateSceneFluidPressureTopologyTransition(
        transition.topologyTransition,
        previousPressureEpoch.pressureControlVolumes,
        previousPressureEpoch.pressureFaceLinks,
        transition.currentPressureEpoch.pressureControlVolumes,
        transition.currentPressureEpoch.pressureFaceLinks);
    if (transition.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure epoch transition exceeds its aggregate limit");
    }
    const auto expected = buildTransition(
        previousPressureEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    if (transition != expected) {
        throw std::invalid_argument(
            "scene fluid pressure epoch-transition sources are foreign");
    }
}

} // namespace simwing::fsi
