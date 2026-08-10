#pragma once

#include "fluid/planar_region_fragment_opening_momentum_transport.h"
#include "scene_fluid_region_wall.h"
#include "scene_fluid_regional_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallInputVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallInputSettings {
    fluid::PlanarPressureRegionFragmentOpeningPressureStateSettings
        pressureState;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallInputSettings&) const =
        default;
};

struct SceneFluidRegionalOpeningMomentumWallInputLimits {
    fluid::PlanarPressureRegionFragmentOpeningVelocityStateLimits flowState;
    fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits loadState;
    SceneFluidRegionalPressureSamplingLimits sampling;
    std::size_t maximumControlVolumes = 20'000'000;
    std::size_t maximumSamples = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16384ULL * 1024ULL * 1024ULL;
};

// Read-only same-epoch ownership adapter from accepted opening-aware
// collocated transport to the source-neutral material-wall kernel. The
// transport target flow is rebuilt from the current accepted pressure epoch;
// the connected-gauge pressure/load state then maps every retained-solid scene
// quadrature sample through its exact regional pressure wall to the two
// transported fragment owners. The output descriptors contain zero wall
// impulse and zero traction: no exchange has run yet.
//
// This boundary does not consume the cycle's consecutive next-pressure state,
// mutate fluid or Structure state, rebase topology, or select a worker.
struct SceneFluidRegionalOpeningMomentumWallInput {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallInputVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourceCurrentFlowStateFingerprint = 0;
    std::uint64_t sourceTransportMetricFingerprint = 0;
    std::uint64_t sourcePressureOperatorFingerprint = 0;
    std::uint64_t sourceBasePressureOperatorFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    std::uint64_t sourceLoadStateFingerprint = 0;
    std::uint64_t sourceSamplingFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    SceneFluidRegionalOpeningMomentumWallInputSettings settings;
    std::size_t activeControlVolumeCount = 0;
    double wallSampleAreaSquareMeters = 0.0;
    double controlIncidentWallAreaSquareMeters = 0.0;
    double maximumIncidentWallAreaSquareMeters = 0.0;
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes;
    std::vector<SceneFluidRegionWallSample> samples;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallInput&) const = default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallInput
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits = {});

void validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(
    const SceneFluidRegionalOpeningMomentumWallInput& input);

// Rebuilds the complete current accepted-flow, retained-solid sampling, and
// fragment/quadrature ownership map and requires bit-identical descriptors.
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
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits = {});

} // namespace simwing::fsi
