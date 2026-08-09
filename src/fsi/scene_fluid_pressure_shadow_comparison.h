#pragma once

#include "scene_fluid_mimetic_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureShadowComparisonVersion = 1;

struct SceneFluidPressureShadowComparisonLimits {
    std::size_t maximumSamples = 10'000'000;
    std::size_t maximumNodes = 10'000'000;
    std::size_t maximumOwnedBytes =
        1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureShadowSampleComparison {
    std::size_t sampleIndex = 0;
    std::uint64_t stableId = 0;
    double referencePressureDifferencePascals = 0.0;
    double shadowPressureDifferencePascals = 0.0;
    double shadowMinusReferencePascals = 0.0;

    bool operator==(
        const SceneFluidPressureShadowSampleComparison&) const = default;
};

struct SceneFluidPressureShadowNodeComparison {
    std::size_t loadIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;
    StructureVector3 referenceForceNewtons;
    StructureVector3 shadowForceNewtons;
    StructureVector3 shadowMinusReferenceNewtons;

    bool operator==(
        const SceneFluidPressureShadowNodeComparison&) const = default;
};

struct SceneFluidPressureShadowComparisonDiagnostics {
    std::size_t sampleCount = 0;
    std::size_t nodeCount = 0;
    double referencePressureDifferenceL2Pascals = 0.0;
    double shadowPressureDifferenceL2Pascals = 0.0;
    double pressureDifferenceDeltaL2Pascals = 0.0;
    double pressureDifferenceDeltaRmsPascals = 0.0;
    double maximumAbsolutePressureDifferenceDeltaPascals = 0.0;
    double relativePressureDifferenceDeltaL2 = 0.0;
    double nodalForceDeltaL2Newtons = 0.0;
    double maximumNodalForceDeltaNewtons = 0.0;
    ConservativeTransferDiagnostics referenceTransfer;
    ConservativeTransferDiagnostics shadowTransfer;
    StructureVector3 shadowMinusReferenceForceNewtons;
    double forceDeltaNormNewtons = 0.0;
    double relativeForceDelta = 0.0;
    StructureVector3 shadowMinusReferenceMomentNewtonMeters;
    double momentDeltaNormNewtonMeters = 0.0;
    double relativeMomentDelta = 0.0;
    double shadowMinusReferencePowerWatts = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureShadowComparisonDiagnostics&) const =
        default;
};

// Immutable diagnostic comparison at the existing conservative pressure-load
// boundary. "Reference" is normally the accepted graph projection and
// "shadow" is the accepted mimetic state. Both overloads retain every
// material pressure-jump delta and every nodal force delta, but neither adds
// loads to Structure or changes solver selection.
struct SceneFluidPressureShadowComparison {
    std::uint32_t version =
        sceneFluidPressureShadowComparisonVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t referenceSampleFingerprint = 0;
    std::uint64_t shadowSampleFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::vector<SceneFluidPressureShadowSampleComparison> samples;
    std::vector<SceneFluidPressureShadowNodeComparison> nodes;
    SceneFluidPressureShadowComparisonDiagnostics diagnostics;

    bool operator==(
        const SceneFluidPressureShadowComparison&) const = default;
};

[[nodiscard]] SceneFluidPressureShadowComparison
compareSceneFluidPressureShadow(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureSampleSet& referenceSamples,
    const SceneFluidMimeticPressureSampleSet& shadowSamples,
    const ConservativeTransferSettings& transferSettings = {},
    const SceneFluidPressureShadowComparisonLimits& limits = {});

// Independent mimetic-vs-mimetic overload used to compare an accepted audit
// endpoint with a separately assembled oracle on geometry where the graph
// operator is deliberately inadmissible.
[[nodiscard]] SceneFluidPressureShadowComparison
compareSceneFluidPressureShadow(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticPressureSampleSet& referenceSamples,
    const SceneFluidMimeticPressureSampleSet& shadowSamples,
    const ConservativeTransferSettings& transferSettings = {},
    const SceneFluidPressureShadowComparisonLimits& limits = {});

void validateSceneFluidPressureShadowComparisonIntegrity(
    const SceneFluidPressureShadowComparison& comparison);

} // namespace simwing::fsi
