#include "scene_fluid_pressure_shadow_comparison.h"

#include <algorithm>
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

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double corrected = value - correction_;
        const double updated = sum_ + corrected;
        correction_ = (updated - sum_) - corrected;
        sum_ = updated;
    }

    [[nodiscard]] double value() const noexcept { return sum_; }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

bool finiteVector(const StructureVector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

StructureVector3 subtract(
    const StructureVector3& first,
    const StructureVector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

double norm(const StructureVector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

double relativeDifference(
    const double difference,
    const double first,
    const double second) {
    const double scale = std::max(std::abs(first), std::abs(second));
    return scale > 0.0 ? difference / scale : 0.0;
}

std::size_t storageBytes(
    const SceneFluidPressureShadowComparison& comparison) {
    if (comparison.samples.size()
            > std::numeric_limits<std::size_t>::max()
                / sizeof(SceneFluidPressureShadowSampleComparison)
        || comparison.nodes.size()
            > std::numeric_limits<std::size_t>::max()
                / sizeof(SceneFluidPressureShadowNodeComparison)) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison storage overflows");
    }
    const std::size_t sampleBytes = comparison.samples.size()
        * sizeof(SceneFluidPressureShadowSampleComparison);
    const std::size_t nodeBytes = comparison.nodes.size()
        * sizeof(SceneFluidPressureShadowNodeComparison);
    if (nodeBytes > std::numeric_limits<std::size_t>::max() - sampleBytes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison storage overflows");
    }
    return sampleBytes + nodeBytes;
}

void fingerprintVector(
    Fingerprint& fingerprint,
    const StructureVector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void fingerprintTransfer(
    Fingerprint& fingerprint,
    const ConservativeTransferDiagnostics& transfer) {
    fingerprint.integer(static_cast<std::uint64_t>(transfer.nodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(transfer.triangleCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        transfer.quadraturePointCount));
    fingerprint.real(transfer.surfaceAreaSquareMeters);
    fingerprintVector(fingerprint, transfer.momentReferenceMeters);
    fingerprintVector(fingerprint, transfer.integratedSurfaceForceNewtons);
    fingerprintVector(fingerprint, transfer.transferredNodalForceNewtons);
    fingerprintVector(fingerprint, transfer.forceResidualNewtons);
    fingerprint.real(transfer.forceResidualNormNewtons);
    fingerprintVector(
        fingerprint, transfer.integratedSurfaceMomentNewtonMeters);
    fingerprintVector(
        fingerprint, transfer.transferredNodalMomentNewtonMeters);
    fingerprintVector(fingerprint, transfer.momentResidualNewtonMeters);
    fingerprint.real(transfer.momentResidualNormNewtonMeters);
    fingerprint.real(transfer.integratedSurfacePowerWatts);
    fingerprint.real(transfer.transferredNodalPowerWatts);
    fingerprint.real(transfer.powerResidualWatts);
    fingerprint.integer(static_cast<std::uint8_t>(transfer.finite));
}

std::uint64_t productFingerprint(
    const SceneFluidPressureShadowComparison& comparison) {
    Fingerprint fingerprint;
    fingerprint.integer(comparison.version);
    fingerprint.integer(comparison.referenceSampleFingerprint);
    fingerprint.integer(comparison.shadowSampleFingerprint);
    fingerprint.integer(comparison.quadratureFingerprint);
    fingerprint.integer(comparison.surfaceDefinitionFingerprint);
    fingerprint.integer(comparison.couplingSurfaceFingerprint);
    fingerprint.integer(comparison.structureDefinitionFingerprint);
    fingerprint.integer(comparison.acceptedStepCount);
    fingerprint.real(comparison.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        comparison.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        comparison.samples.size()));
    for (const auto& sample : comparison.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(sample.stableId);
        fingerprint.real(sample.referencePressureDifferencePascals);
        fingerprint.real(sample.shadowPressureDifferencePascals);
        fingerprint.real(sample.shadowMinusReferencePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        comparison.nodes.size()));
    for (const auto& node : comparison.nodes) {
        fingerprint.integer(static_cast<std::uint64_t>(node.loadIndex));
        fingerprint.integer(node.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(node.structureNode));
        fingerprintVector(fingerprint, node.referenceForceNewtons);
        fingerprintVector(fingerprint, node.shadowForceNewtons);
        fingerprintVector(fingerprint, node.shadowMinusReferenceNewtons);
    }
    const auto& diagnostics = comparison.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.sampleCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.nodeCount));
    fingerprint.real(diagnostics.referencePressureDifferenceL2Pascals);
    fingerprint.real(diagnostics.shadowPressureDifferenceL2Pascals);
    fingerprint.real(diagnostics.pressureDifferenceDeltaL2Pascals);
    fingerprint.real(diagnostics.pressureDifferenceDeltaRmsPascals);
    fingerprint.real(
        diagnostics.maximumAbsolutePressureDifferenceDeltaPascals);
    fingerprint.real(diagnostics.relativePressureDifferenceDeltaL2);
    fingerprint.real(diagnostics.nodalForceDeltaL2Newtons);
    fingerprint.real(diagnostics.maximumNodalForceDeltaNewtons);
    fingerprintTransfer(fingerprint, diagnostics.referenceTransfer);
    fingerprintTransfer(fingerprint, diagnostics.shadowTransfer);
    fingerprintVector(
        fingerprint, diagnostics.shadowMinusReferenceForceNewtons);
    fingerprint.real(diagnostics.forceDeltaNormNewtons);
    fingerprint.real(diagnostics.relativeForceDelta);
    fingerprintVector(
        fingerprint, diagnostics.shadowMinusReferenceMomentNewtonMeters);
    fingerprint.real(diagnostics.momentDeltaNormNewtonMeters);
    fingerprint.real(diagnostics.relativeMomentDelta);
    fingerprint.real(diagnostics.shadowMinusReferencePowerWatts);
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.finite));
    return fingerprint.value();
}

SceneFluidPressureShadowComparisonDiagnostics summarize(
    const SceneFluidPressureShadowComparison& comparison,
    const ConservativeTransferDiagnostics& referenceTransfer,
    const ConservativeTransferDiagnostics& shadowTransfer) {
    SceneFluidPressureShadowComparisonDiagnostics result;
    result.sampleCount = comparison.samples.size();
    result.nodeCount = comparison.nodes.size();
    result.referenceTransfer = referenceTransfer;
    result.shadowTransfer = shadowTransfer;
    CompensatedSum referencePressureSquared;
    CompensatedSum shadowPressureSquared;
    CompensatedSum pressureDeltaSquared;
    for (const auto& sample : comparison.samples) {
        referencePressureSquared.add(
            sample.referencePressureDifferencePascals
            * sample.referencePressureDifferencePascals);
        shadowPressureSquared.add(
            sample.shadowPressureDifferencePascals
            * sample.shadowPressureDifferencePascals);
        pressureDeltaSquared.add(
            sample.shadowMinusReferencePascals
            * sample.shadowMinusReferencePascals);
        result.maximumAbsolutePressureDifferenceDeltaPascals = std::max(
            result.maximumAbsolutePressureDifferenceDeltaPascals,
            std::abs(sample.shadowMinusReferencePascals));
    }
    result.referencePressureDifferenceL2Pascals = std::sqrt(
        std::max(0.0, referencePressureSquared.value()));
    result.shadowPressureDifferenceL2Pascals = std::sqrt(
        std::max(0.0, shadowPressureSquared.value()));
    result.pressureDifferenceDeltaL2Pascals = std::sqrt(
        std::max(0.0, pressureDeltaSquared.value()));
    result.pressureDifferenceDeltaRmsPascals = result.sampleCount == 0
        ? 0.0
        : result.pressureDifferenceDeltaL2Pascals
            / std::sqrt(static_cast<double>(result.sampleCount));
    result.relativePressureDifferenceDeltaL2 = relativeDifference(
        result.pressureDifferenceDeltaL2Pascals,
        result.referencePressureDifferenceL2Pascals,
        result.shadowPressureDifferenceL2Pascals);

    CompensatedSum nodalDeltaSquared;
    for (const auto& node : comparison.nodes) {
        const double delta = norm(node.shadowMinusReferenceNewtons);
        nodalDeltaSquared.add(delta * delta);
        result.maximumNodalForceDeltaNewtons = std::max(
            result.maximumNodalForceDeltaNewtons, delta);
    }
    result.nodalForceDeltaL2Newtons = std::sqrt(
        std::max(0.0, nodalDeltaSquared.value()));
    result.shadowMinusReferenceForceNewtons = subtract(
        shadowTransfer.integratedSurfaceForceNewtons,
        referenceTransfer.integratedSurfaceForceNewtons);
    result.forceDeltaNormNewtons = norm(
        result.shadowMinusReferenceForceNewtons);
    result.relativeForceDelta = relativeDifference(
        result.forceDeltaNormNewtons,
        norm(referenceTransfer.integratedSurfaceForceNewtons),
        norm(shadowTransfer.integratedSurfaceForceNewtons));
    result.shadowMinusReferenceMomentNewtonMeters = subtract(
        shadowTransfer.integratedSurfaceMomentNewtonMeters,
        referenceTransfer.integratedSurfaceMomentNewtonMeters);
    result.momentDeltaNormNewtonMeters = norm(
        result.shadowMinusReferenceMomentNewtonMeters);
    result.relativeMomentDelta = relativeDifference(
        result.momentDeltaNormNewtonMeters,
        norm(referenceTransfer.integratedSurfaceMomentNewtonMeters),
        norm(shadowTransfer.integratedSurfaceMomentNewtonMeters));
    result.shadowMinusReferencePowerWatts =
        shadowTransfer.integratedSurfacePowerWatts
        - referenceTransfer.integratedSurfacePowerWatts;
    result.finite = referenceTransfer.finite && shadowTransfer.finite
        && std::isfinite(result.referencePressureDifferenceL2Pascals)
        && std::isfinite(result.shadowPressureDifferenceL2Pascals)
        && std::isfinite(result.pressureDifferenceDeltaL2Pascals)
        && std::isfinite(result.pressureDifferenceDeltaRmsPascals)
        && std::isfinite(
            result.maximumAbsolutePressureDifferenceDeltaPascals)
        && std::isfinite(result.relativePressureDifferenceDeltaL2)
        && std::isfinite(result.nodalForceDeltaL2Newtons)
        && std::isfinite(result.maximumNodalForceDeltaNewtons)
        && finiteVector(result.shadowMinusReferenceForceNewtons)
        && std::isfinite(result.forceDeltaNormNewtons)
        && std::isfinite(result.relativeForceDelta)
        && finiteVector(result.shadowMinusReferenceMomentNewtonMeters)
        && std::isfinite(result.momentDeltaNormNewtonMeters)
        && std::isfinite(result.relativeMomentDelta)
        && std::isfinite(result.shadowMinusReferencePowerWatts);
    return result;
}

bool sameBindingIdentity(
    const SceneFluidPressureSampleBinding& first,
    const SceneFluidPressureSampleBinding& second) {
    return first.sampleIndex == second.sampleIndex
        && first.stableId == second.stableId
        && first.negativeSideControlVolumeIndex
            == second.negativeSideControlVolumeIndex
        && first.positiveSideControlVolumeIndex
            == second.positiveSideControlVolumeIndex
        && first.negativeSideControlVolumeStableId
            == second.negativeSideControlVolumeStableId
        && first.positiveSideControlVolumeStableId
            == second.positiveSideControlVolumeStableId
        && first.negativeSideRegionId == second.negativeSideRegionId
        && first.positiveSideRegionId == second.positiveSideRegionId
        && first.componentIndex == second.componentIndex;
}

template<typename ReferenceSamples, typename ShadowSamples>
SceneFluidPressureShadowComparison buildComparison(
    const ReferenceSamples& referenceSamples,
    const ShadowSamples& shadowSamples,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const ConservativeTransferResult& referenceTransfer,
    const ConservativeTransferResult& shadowTransfer,
    const SceneFluidPressureShadowComparisonLimits& limits) {
    if (referenceSamples.bindings.size() > limits.maximumSamples
        || referenceTransfer.nodeLoads().size() > limits.maximumNodes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison exceeds its count limit");
    }
    if (referenceSamples.fingerprint == 0
        || shadowSamples.fingerprint == 0
        || referenceSamples.quadratureFingerprint != quadrature.fingerprint
        || shadowSamples.quadratureFingerprint != quadrature.fingerprint
        || referenceSamples.surfaceDefinitionFingerprint
            != surface.fingerprint
        || shadowSamples.surfaceDefinitionFingerprint != surface.fingerprint
        || referenceSamples.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || shadowSamples.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || referenceSamples.acceptedStepCount != state.acceptedStepCount
        || shadowSamples.acceptedStepCount != state.acceptedStepCount
        || referenceSamples.simulationTimeSeconds
            != state.simulationTimeSeconds
        || shadowSamples.simulationTimeSeconds
            != state.simulationTimeSeconds
        || referenceSamples.bindings.size()
            != shadowSamples.bindings.size()
        || referenceSamples.pressures.size()
            != shadowSamples.pressures.size()
        || referenceSamples.bindings.size()
            != referenceSamples.pressures.size()
        || referenceTransfer.surfaceFingerprint()
            != transfer.couplingSurfaceFingerprint()
        || shadowTransfer.surfaceFingerprint()
            != transfer.couplingSurfaceFingerprint()
        || referenceTransfer.targetDefinitionFingerprint()
            != state.structureDefinitionFingerprint
        || shadowTransfer.targetDefinitionFingerprint()
            != state.structureDefinitionFingerprint
        || referenceTransfer.nodeLoads().size()
            != shadowTransfer.nodeLoads().size()) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow comparison input identity is invalid");
    }

    SceneFluidPressureShadowComparison result;
    result.referenceSampleFingerprint = referenceSamples.fingerprint;
    result.shadowSampleFingerprint = shadowSamples.fingerprint;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.couplingSurfaceFingerprint =
        transfer.couplingSurfaceFingerprint();
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.samples.reserve(referenceSamples.bindings.size());
    for (std::size_t index = 0;
         index < referenceSamples.bindings.size(); ++index) {
        const auto& referenceBinding = referenceSamples.bindings[index];
        const auto& shadowBinding = shadowSamples.bindings[index];
        const auto& referencePressure = referenceSamples.pressures[index];
        const auto& shadowPressure = shadowSamples.pressures[index];
        if (!sameBindingIdentity(referenceBinding, shadowBinding)
            || referenceBinding.sampleIndex != index
            || referencePressure.stableId != referenceBinding.stableId
            || shadowPressure.stableId != shadowBinding.stableId) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow sample binding is invalid");
        }
        result.samples.push_back({
            index,
            referenceBinding.stableId,
            referenceBinding.pressureDifferencePascals,
            shadowBinding.pressureDifferencePascals,
            shadowBinding.pressureDifferencePascals
                - referenceBinding.pressureDifferencePascals,
        });
    }
    result.nodes.reserve(referenceTransfer.nodeLoads().size());
    for (std::size_t index = 0;
         index < referenceTransfer.nodeLoads().size(); ++index) {
        const auto& reference = referenceTransfer.nodeLoads()[index];
        const auto& shadow = shadowTransfer.nodeLoads()[index];
        if (reference.stableId != shadow.stableId
            || reference.structureNode != shadow.structureNode) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow nodal binding is invalid");
        }
        result.nodes.push_back({
            index,
            reference.stableId,
            reference.structureNode,
            reference.forceNewtons,
            shadow.forceNewtons,
            subtract(shadow.forceNewtons, reference.forceNewtons),
        });
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison exceeds its byte limit");
    }
    result.diagnostics = summarize(
        result, referenceTransfer.diagnostics(), shadowTransfer.diagnostics());
    if (!result.diagnostics.finite) {
        throw std::overflow_error(
            "scene fluid pressure-shadow comparison is non-finite");
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidPressureShadowComparisonIntegrity(result);
    return result;
}

} // namespace

SceneFluidPressureShadowComparison compareSceneFluidPressureShadow(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidPressureSampleSet& referenceSamples,
    const SceneFluidMimeticPressureSampleSet& shadowSamples,
    const ConservativeTransferSettings& transferSettings,
    const SceneFluidPressureShadowComparisonLimits& limits) {
    validateSceneFluidPressureSampleIntegrity(referenceSamples);
    validateSceneFluidMimeticPressureSampleIntegrity(shadowSamples);
    const auto referenceTransfer =
        evaluateSceneFluidProjectedPressureQuadrature(
            surface, state, transfer, quadrature, referenceSamples,
            transferSettings);
    const auto shadowTransfer = evaluateSceneFluidMimeticPressureQuadrature(
        surface, state, transfer, quadrature, shadowSamples,
        transferSettings);
    return buildComparison(
        referenceSamples, shadowSamples, surface, state, transfer,
        quadrature, referenceTransfer, shadowTransfer, limits);
}

SceneFluidPressureShadowComparison compareSceneFluidPressureShadow(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticPressureSampleSet& referenceSamples,
    const SceneFluidMimeticPressureSampleSet& shadowSamples,
    const ConservativeTransferSettings& transferSettings,
    const SceneFluidPressureShadowComparisonLimits& limits) {
    validateSceneFluidMimeticPressureSampleIntegrity(referenceSamples);
    validateSceneFluidMimeticPressureSampleIntegrity(shadowSamples);
    const auto referenceTransfer = evaluateSceneFluidMimeticPressureQuadrature(
        surface, state, transfer, quadrature, referenceSamples,
        transferSettings);
    const auto shadowTransfer = evaluateSceneFluidMimeticPressureQuadrature(
        surface, state, transfer, quadrature, shadowSamples,
        transferSettings);
    return buildComparison(
        referenceSamples, shadowSamples, surface, state, transfer,
        quadrature, referenceTransfer, shadowTransfer, limits);
}

void validateSceneFluidPressureShadowComparisonIntegrity(
    const SceneFluidPressureShadowComparison& comparison) {
    if (comparison.version != sceneFluidPressureShadowComparisonVersion
        || comparison.fingerprint == 0
        || comparison.referenceSampleFingerprint == 0
        || comparison.shadowSampleFingerprint == 0
        || comparison.quadratureFingerprint == 0
        || comparison.surfaceDefinitionFingerprint == 0
        || comparison.couplingSurfaceFingerprint == 0
        || comparison.structureDefinitionFingerprint == 0
        || comparison.samples.empty()
        || comparison.nodes.empty()
        || comparison.ownedStorageBytes != storageBytes(comparison)) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow comparison integrity is invalid");
    }
    for (std::size_t index = 0; index < comparison.samples.size(); ++index) {
        const auto& sample = comparison.samples[index];
        if (sample.sampleIndex != index || sample.stableId == 0
            || !std::isfinite(sample.referencePressureDifferencePascals)
            || !std::isfinite(sample.shadowPressureDifferencePascals)
            || sample.shadowMinusReferencePascals
                != sample.shadowPressureDifferencePascals
                    - sample.referencePressureDifferencePascals) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow comparison sample is invalid");
        }
    }
    for (std::size_t index = 0; index < comparison.nodes.size(); ++index) {
        const auto& node = comparison.nodes[index];
        if (node.loadIndex != index || node.stableId == 0
            || !finiteVector(node.referenceForceNewtons)
            || !finiteVector(node.shadowForceNewtons)
            || node.shadowMinusReferenceNewtons
                != subtract(
                    node.shadowForceNewtons, node.referenceForceNewtons)) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow comparison node is invalid");
        }
    }
    const auto expected = summarize(
        comparison, comparison.diagnostics.referenceTransfer,
        comparison.diagnostics.shadowTransfer);
    const auto& reference = comparison.diagnostics.referenceTransfer;
    const auto& shadow = comparison.diagnostics.shadowTransfer;
    if (!comparison.diagnostics.finite
        || reference.nodeCount != comparison.nodes.size()
        || shadow.nodeCount != comparison.nodes.size()
        || reference.quadraturePointCount != comparison.samples.size()
        || shadow.quadraturePointCount != comparison.samples.size()
        || reference.nodeCount != shadow.nodeCount
        || reference.triangleCount != shadow.triangleCount
        || reference.quadraturePointCount != shadow.quadraturePointCount
        || reference.surfaceAreaSquareMeters
            != shadow.surfaceAreaSquareMeters
        || reference.momentReferenceMeters != shadow.momentReferenceMeters
        || comparison.diagnostics != expected
        || comparison.fingerprint != productFingerprint(comparison)) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow comparison summary is invalid");
    }
}

} // namespace simwing::fsi
