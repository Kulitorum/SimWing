#pragma once

#include "scene_fluid_pressure_projection.h"
#include "scene_fluid_pressure_traction.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureSamplingVersion = 1;

struct SceneFluidPressureSamplingLimits {
    std::size_t maximumSamples = 10'000'000;
    std::size_t maximumSamplingBytes =
        1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureSampleBinding {
    std::size_t sampleIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t negativeSideControlVolumeIndex = 0;
    std::size_t positiveSideControlVolumeIndex = 0;
    std::uint64_t negativeSideControlVolumeStableId = 0;
    std::uint64_t positiveSideControlVolumeStableId = 0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    std::size_t componentIndex = 0;
    double pressureDifferencePascals = 0.0;

    bool operator==(const SceneFluidPressureSampleBinding&) const = default;
};

// Immutable one-sided pressure sampling for the conservative surface
// quadrature. Each authored side resolves to the exact cell/region pressure
// unknown recorded by quadrature-v2. Both sides must share one pressure
// component; otherwise their independently gauged absolute values cannot form
// a physical sheet pressure jump. This adapter only samples an already
// accepted projection and does not apply Structure loads.
struct SceneFluidPressureSampleSet {
    std::uint32_t version = sceneFluidPressureSamplingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t pressureProjectionFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    double maximumAbsolutePressureDifferencePascals = 0.0;
    std::vector<SceneFluidPressureSampleBinding> bindings;
    std::vector<SceneFluidQuadraturePressure> pressures;

    bool operator==(const SceneFluidPressureSampleSet&) const = default;
};

[[nodiscard]] SceneFluidPressureSampleSet
sampleSceneFluidProjectedPressure(
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureProjection& projection,
    const SceneFluidPressureSamplingLimits& limits = {});

void validateSceneFluidPressureSampleIntegrity(
    const SceneFluidPressureSampleSet& samples);

void validateSceneFluidProjectedPressureSamples(
    const SceneFluidPressureSampleSet& samples,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureProjection& projection);

[[nodiscard]] ConservativeTransferResult
evaluateSceneFluidProjectedPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureSampleSet& samples,
    const ConservativeTransferSettings& settings = {});

} // namespace simwing::fsi
