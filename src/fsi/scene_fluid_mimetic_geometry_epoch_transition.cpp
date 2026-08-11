#include "scene_fluid_mimetic_geometry_epoch_transition.h"

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
            "scene fluid mimetic geometry transition storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const SceneFluidMimeticGeometryEpochTransition& transition) {
    return checkedAdd(
        transition.currentGeometryEpoch.ownedStorageBytes,
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
    const SceneFluidMimeticGeometryEpochTransition& transition) {
    Fingerprint fingerprint;
    fingerprint.integer(transition.version);
    for (const std::uint64_t value : {
             transition.previousGeometryEpochFingerprint,
             transition.acceptedCurrentGridEpochFingerprint,
             transition.currentGeometryEpochFingerprint,
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
    fingerprint.integer(static_cast<std::uint8_t>(
        transition.controlVolumeTopologyStable));
    fingerprint.integer(static_cast<std::uint64_t>(
        transition.ownedStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry transition limit is invalid");
    }
}

void validateSources(
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch) {
    validateSceneFluidMimeticGeometryEpoch(
        previousGeometryEpoch, surface, previousSurfaceState, grid,
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
        || previousGeometryEpoch.surfaceStateFingerprint
            != previousSurfaceState.fingerprint
        || acceptedCurrentGridEpoch.quadrature.surfaceStateFingerprint
            != currentSurfaceState.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry transition sources are not consecutive");
    }
}

SceneFluidMimeticGeometryEpochTransition buildTransition(
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidMimeticGeometryEpochSettings& settings,
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits) {
    validateLimits(limits);
    validateSources(
        previousGeometryEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch);

    SceneFluidMimeticGeometryEpochTransition result;
    result.previousGeometryEpochFingerprint =
        previousGeometryEpoch.fingerprint;
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
    result.currentGeometryEpoch = buildSceneFluidMimeticGeometryEpoch(
        surface, currentSurfaceState, grid, transfer, connectivity,
        settings, limits.geometryEpoch);
    if (result.currentGeometryEpoch.gridEpoch
            != acceptedCurrentGridEpoch) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry transition rebuilt a foreign grid epoch");
    }
    result.currentGeometryEpochFingerprint =
        result.currentGeometryEpoch.fingerprint;
    result.topologyTransition = buildSceneFluidPressureTopologyTransition(
        previousGeometryEpoch.pressureControlVolumes,
        previousGeometryEpoch.pressureFaceLinks,
        result.currentGeometryEpoch.pressureControlVolumes,
        result.currentGeometryEpoch.pressureFaceLinks,
        limits.topologyTransition);
    result.topologyTransitionFingerprint =
        result.topologyTransition.fingerprint;
    result.controlVolumeTopologyStable =
        hasStableControlVolumeTopology(result.topologyTransition);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic geometry transition exceeds its aggregate limit");
    }
    result.fingerprint = transitionFingerprint(result);
    return result;
}

} // namespace

SceneFluidMimeticGeometryEpochTransition
buildSceneFluidMimeticGeometryEpochTransition(
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidMimeticGeometryEpochSettings& settings,
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits) {
    auto result = buildTransition(
        previousGeometryEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    validateSceneFluidMimeticGeometryEpochTransition(
        result, previousGeometryEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    return result;
}

void validateSceneFluidMimeticGeometryEpochTransitionIntegrity(
    const SceneFluidMimeticGeometryEpochTransition& transition) {
    validateSceneFluidMimeticGeometryEpochIntegrity(
        transition.currentGeometryEpoch);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        transition.topologyTransition);
    if (transition.version
            != sceneFluidMimeticGeometryEpochTransitionVersion
        || transition.fingerprint == 0
        || transition.previousGeometryEpochFingerprint == 0
        || transition.acceptedCurrentGridEpochFingerprint == 0
        || transition.currentGeometryEpochFingerprint == 0
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
        || transition.currentGeometryEpoch.fingerprint
            != transition.currentGeometryEpochFingerprint
        || transition.currentGeometryEpoch.surfaceDefinitionFingerprint
            != transition.surfaceDefinitionFingerprint
        || transition.currentGeometryEpoch.structureDefinitionFingerprint
            != transition.structureDefinitionFingerprint
        || transition.currentGeometryEpoch.surfaceStateFingerprint
            != transition.currentSurfaceStateFingerprint
        || transition.topologyTransition.fingerprint
            != transition.topologyTransitionFingerprint
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
            "scene fluid mimetic geometry transition integrity is invalid");
    }
}

void validateSceneFluidMimeticGeometryEpochTransition(
    const SceneFluidMimeticGeometryEpochTransition& transition,
    const SceneFluidMimeticGeometryEpoch& previousGeometryEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidMimeticGeometryEpochSettings& settings,
    const SceneFluidMimeticGeometryEpochTransitionLimits& limits) {
    validateLimits(limits);
    validateSceneFluidMimeticGeometryEpochTransitionIntegrity(transition);
    validateSources(
        previousGeometryEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch);
    validateSceneFluidMimeticGeometryEpoch(
        transition.currentGeometryEpoch, surface, currentSurfaceState,
        grid, transfer, connectivity);
    validateSceneFluidPressureTopologyTransition(
        transition.topologyTransition,
        previousGeometryEpoch.pressureControlVolumes,
        previousGeometryEpoch.pressureFaceLinks,
        transition.currentGeometryEpoch.pressureControlVolumes,
        transition.currentGeometryEpoch.pressureFaceLinks);
    if (transition.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic geometry transition exceeds its aggregate limit");
    }
    const auto expected = buildTransition(
        previousGeometryEpoch, surface, previousSurfaceState,
        currentSurfaceState, grid, transfer, connectivity,
        acceptedCurrentGridEpoch, settings, limits);
    if (transition != expected) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry transition sources are foreign");
    }
}

} // namespace simwing::fsi
