#pragma once

#include "scene_fluid_mimetic_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureShadowComparisonVersion = 2;

struct SceneFluidPressureShadowComparisonLimits {
    std::size_t maximumSamples = 10'000'000;
    std::size_t maximumNodes = 10'000'000;
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumComponents = 1'000'000;
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

struct SceneFluidPressureShadowControlSourceComparison {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t componentIndex = 0;
    double referenceGeometryVolumeRateCubicMetersPerSecond = 0.0;
    double shadowGeometryVolumeRateCubicMetersPerSecond = 0.0;
    double geometryVolumeRateDeltaCubicMetersPerSecond = 0.0;
    double referencePredictedNetOutwardVolumeRateCubicMetersPerSecond = 0.0;
    double shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond = 0.0;
    double predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond = 0.0;
    double referenceContinuityResidualCubicMetersPerSecond = 0.0;
    double shadowContinuityResidualCubicMetersPerSecond = 0.0;
    double continuityResidualDeltaCubicMetersPerSecond = 0.0;
    double referenceIntegratedSourcePascalsMeters = 0.0;
    double shadowIntegratedSourcePascalsMeters = 0.0;
    double integratedSourceDeltaPascalsMeters = 0.0;

    bool operator==(
        const SceneFluidPressureShadowControlSourceComparison&) const =
        default;
};

struct SceneFluidPressureShadowScalarComparisonDiagnostics {
    double referenceL2 = 0.0;
    double shadowL2 = 0.0;
    double deltaL2 = 0.0;
    double maximumAbsoluteDelta = 0.0;
    double relativeDeltaL2 = 0.0;
    bool exact = false;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureShadowScalarComparisonDiagnostics&) const =
        default;
};

struct SceneFluidPressureShadowSourceComparisonDiagnostics {
    std::size_t controlVolumeCount = 0;
    std::size_t componentCount = 0;
    SceneFluidPressureShadowScalarComparisonDiagnostics geometryVolumeRate;
    SceneFluidPressureShadowScalarComparisonDiagnostics
        predictedNetOutwardVolumeRate;
    SceneFluidPressureShadowScalarComparisonDiagnostics continuityResidual;
    SceneFluidPressureShadowScalarComparisonDiagnostics integratedSource;
    double maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureShadowSourceComparisonDiagnostics&) const =
        default;
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
    double pressureDifferenceDotProductPascalsSquared = 0.0;
    double bestFitShadowPressureScale = 0.0;
    double pressureDifferenceCosineSimilarity = 0.0;
    double bestFitPressureShapeResidualL2Pascals = 0.0;
    double relativeBestFitPressureShapeResidualL2 = 0.0;
    double nodalForceDeltaL2Newtons = 0.0;
    double maximumNodalForceDeltaNewtons = 0.0;
    double nodalForceDotProductNewtonsSquared = 0.0;
    double bestFitShadowNodalForceScale = 0.0;
    double nodalForceCosineSimilarity = 0.0;
    double bestFitNodalForceShapeResidualL2Newtons = 0.0;
    double relativeBestFitNodalForceShapeResidualL2 = 0.0;
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
    bool includesSourceComparison = false;
    std::uint64_t referencePressureProjectionFingerprint = 0;
    std::uint64_t shadowControlCellFingerprint = 0;
    std::uint64_t shadowPressureSourceFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::vector<SceneFluidPressureShadowSampleComparison> samples;
    std::vector<SceneFluidPressureShadowNodeComparison> nodes;
    std::vector<SceneFluidPressureShadowControlSourceComparison>
        controlSources;
    std::vector<double> componentIntegratedSourceDeltasPascalsMeters;
    SceneFluidPressureShadowSourceComparisonDiagnostics sourceDiagnostics;
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

// Full live graph-vs-mimetic evidence path. In addition to material samples
// and conservative loads, this overload binds every graph projection source
// row to the exact mimetic control/source row so forcing differences can be
// separated from operator differences.
[[nodiscard]] SceneFluidPressureShadowComparison
compareSceneFluidPressureShadow(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureProjection& referenceProjection,
    const SceneFluidPressureSampleSet& referenceSamples,
    const SceneFluidMimeticControlCellSet& shadowControlCells,
    const SceneFluidMimeticPressureSourceSet& shadowSources,
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
