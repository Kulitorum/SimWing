#pragma once

#include "scene_fluid_mimetic_pressure_state.h"
#include "scene_fluid_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticPressureSamplingVersion = 1;

// Immutable one-sided surface samples from an accepted mimetic pressure
// state. Each quadrature side resolves through the exact pressure-control row
// also owned by the mimetic control topology. Both sides must share one
// pressure component so their difference is gauge-safe. The resulting
// one-sided samples use the existing conservative pressure-traction path.
struct SceneFluidMimeticPressureSampleSet {
    std::uint32_t version =
        sceneFluidMimeticPressureSamplingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t pressureStateFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t controlVolumeCount = 0;
    std::size_t componentCount = 0;
    double maximumAbsolutePressureDifferencePascals = 0.0;
    std::vector<SceneFluidPressureSampleBinding> bindings;
    std::vector<SceneFluidQuadraturePressure> pressures;

    bool operator==(
        const SceneFluidMimeticPressureSampleSet&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureSampleSet
sampleSceneFluidMimeticPressure(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState,
    const SceneFluidPressureSamplingLimits& limits = {});

void validateSceneFluidMimeticPressureSampleIntegrity(
    const SceneFluidMimeticPressureSampleSet& samples);

void validateSceneFluidMimeticPressureSamples(
    const SceneFluidMimeticPressureSampleSet& samples,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticPressureState& pressureState);

[[nodiscard]] ConservativeTransferResult
evaluateSceneFluidMimeticPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticPressureSampleSet& samples,
    const ConservativeTransferSettings& settings = {});

} // namespace simwing::fsi
