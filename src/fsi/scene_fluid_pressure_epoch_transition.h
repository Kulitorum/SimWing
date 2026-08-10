#pragma once

#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_pressure_topology_transition.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureEpochTransitionVersion = 1;

struct SceneFluidPressureEpochTransitionLimits {
    SceneFluidPressureEpochLimits pressureEpoch;
    SceneFluidPressureTopologyTransitionLimits topologyTransition;
    std::size_t maximumOwnedBytes = 12ULL * 1024ULL * 1024ULL * 1024ULL;
};

// One immutable consecutive pressure-geometry transaction. It rebuilds the
// complete current pressure epoch, requires its grid geometry to equal the
// independently accepted current grid epoch, and owns the shared stable-ID
// control-volume transition from the previous pressure epoch.
//
// controlVolumeTopologyStable means only that no pressure control row appeared
// or disappeared. Face geometry, operator coefficients, and every current
// pressure product are still rebuilt and may differ. This receipt samples no
// velocity, rebases no fluid state, solves no pressure, and applies no load.
struct SceneFluidPressureEpochTransition {
    std::uint32_t version = sceneFluidPressureEpochTransitionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t previousPressureEpochFingerprint = 0;
    std::uint64_t acceptedCurrentGridEpochFingerprint = 0;
    std::uint64_t currentPressureEpochFingerprint = 0;
    std::uint64_t topologyTransitionFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t previousSurfaceStateFingerprint = 0;
    std::uint64_t currentSurfaceStateFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double previousSimulationTimeSeconds = 0.0;
    double currentSimulationTimeSeconds = 0.0;
    bool controlVolumeTopologyStable = false;
    SceneFluidPressureEpoch currentPressureEpoch;
    SceneFluidPressureTopologyTransition topologyTransition;
    std::size_t ownedStorageBytes = 0;

    bool operator==(const SceneFluidPressureEpochTransition&) const = default;
};

[[nodiscard]] SceneFluidPressureEpochTransition
buildSceneFluidPressureEpochTransition(
    const SceneFluidPressureEpoch& previousPressureEpoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& previousSurfaceState,
    const SceneFluidSurfaceState& currentSurfaceState,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidGridEpoch& acceptedCurrentGridEpoch,
    const SceneFluidPressureEpochSettings& settings = {},
    const SceneFluidPressureEpochTransitionLimits& limits = {});

void validateSceneFluidPressureEpochTransitionIntegrity(
    const SceneFluidPressureEpochTransition& transition);

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
    const SceneFluidPressureEpochSettings& settings = {},
    const SceneFluidPressureEpochTransitionLimits& limits = {});

} // namespace simwing::fsi
