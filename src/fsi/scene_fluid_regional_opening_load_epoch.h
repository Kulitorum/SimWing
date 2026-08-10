#pragma once

#include "scene_fluid_regional_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningLoadEpochVersion = 1;

struct SceneFluidRegionalOpeningLoadEpochSettings {
    fluid::PlanarPressureRegionFragmentOpeningPressureStateSettings
        pressureState;
    ConservativeTransferSettings transfer;
};

struct SceneFluidRegionalOpeningLoadEpochLimits {
    fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits loadState;
    SceneFluidRegionalPressureSamplingLimits sampling;
    SceneFluidRegionalPressureLoadApplicationLimits application;
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// One explicit accepted-flow-to-pending-load transaction. The adapter
// privately composes the connected-gauge pressure state, full-wall load,
// retained-solid opening partition, authoritative scene samples, and
// conservative nodal application. Only the final application mutates the
// target, and any later validation or aggregate-limit failure restores the
// exact Structure checkpoint captured immediately before that mutation.
//
// This endpoint neither steps Structure nor consumes pending loads. It owns no
// fluid advection, topology rebase, strong-coupling iteration, persistence, or
// production worker selection.
struct SceneFluidRegionalOpeningLoadEpoch {
    std::uint32_t version =
        sceneFluidRegionalOpeningLoadEpochVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourcePressureOperatorFingerprint = 0;
    std::uint64_t sourceBasePressureOperatorFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t targetDefinitionFingerprint = 0;
    std::uint64_t sourceSettingsFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::PlanarPressureRegionFragmentOpeningLoadState loadState;
    SceneFluidRegionalPressureSampleSet samples;
    SceneFluidRegionalPressureLoadApplication application;
    bool applied = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningLoadEpoch&) const = default;
};

[[nodiscard]] SceneFluidRegionalOpeningLoadEpoch
applySceneFluidRegionalOpeningLoadEpoch(
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        acceptedState,
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningLoadEpochLimits& limits = {});

void validateSceneFluidRegionalOpeningLoadEpochIntegrity(
    const SceneFluidRegionalOpeningLoadEpoch& epoch);

// Rebuilds every immutable pre-application stage and re-evaluates the applied
// node loads against the exact scene transfer. Historical pending loads remain
// self-contained receipt data and no Structure is mutated by validation.
void validateSceneFluidRegionalOpeningLoadEpoch(
    const SceneFluidRegionalOpeningLoadEpoch& epoch,
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        acceptedState,
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
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningLoadEpochLimits& limits = {});

} // namespace simwing::fsi
