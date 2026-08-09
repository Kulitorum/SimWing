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
                / sizeof(SceneFluidPressureShadowNodeComparison)
        || comparison.controlSources.size()
            > std::numeric_limits<std::size_t>::max()
                / sizeof(SceneFluidPressureShadowControlSourceComparison)
        || comparison.componentIntegratedSourceDeltasPascalsMeters.size()
            > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison storage overflows");
    }
    const std::size_t sampleBytes = comparison.samples.size()
        * sizeof(SceneFluidPressureShadowSampleComparison);
    const std::size_t nodeBytes = comparison.nodes.size()
        * sizeof(SceneFluidPressureShadowNodeComparison);
    const std::size_t controlBytes = comparison.controlSources.size()
        * sizeof(SceneFluidPressureShadowControlSourceComparison);
    const std::size_t componentBytes = comparison
        .componentIntegratedSourceDeltasPascalsMeters.size()
        * sizeof(double);
    if (nodeBytes > std::numeric_limits<std::size_t>::max() - sampleBytes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison storage overflows");
    }
    const std::size_t sampleAndNodeBytes = sampleBytes + nodeBytes;
    if (controlBytes
            > std::numeric_limits<std::size_t>::max() - sampleAndNodeBytes
        || componentBytes > std::numeric_limits<std::size_t>::max()
                - sampleAndNodeBytes - controlBytes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison storage overflows");
    }
    return sampleAndNodeBytes + controlBytes + componentBytes;
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

void fingerprintScalarComparison(
    Fingerprint& fingerprint,
    const SceneFluidPressureShadowScalarComparisonDiagnostics& diagnostics) {
    fingerprint.real(diagnostics.referenceL2);
    fingerprint.real(diagnostics.shadowL2);
    fingerprint.real(diagnostics.deltaL2);
    fingerprint.real(diagnostics.maximumAbsoluteDelta);
    fingerprint.real(diagnostics.relativeDeltaL2);
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.exact));
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.finite));
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
    fingerprint.integer(static_cast<std::uint8_t>(
        comparison.includesSourceComparison));
    fingerprint.integer(comparison.referencePressureProjectionFingerprint);
    fingerprint.integer(comparison.shadowControlCellFingerprint);
    fingerprint.integer(comparison.shadowPressureSourceFingerprint);
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
    fingerprint.integer(static_cast<std::uint64_t>(
        comparison.controlSources.size()));
    for (const auto& source : comparison.controlSources) {
        fingerprint.integer(static_cast<std::uint64_t>(
            source.controlVolumeIndex));
        fingerprint.integer(source.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            source.componentIndex));
        fingerprint.real(
            source.referenceGeometryVolumeRateCubicMetersPerSecond);
        fingerprint.real(
            source.shadowGeometryVolumeRateCubicMetersPerSecond);
        fingerprint.real(
            source.geometryVolumeRateDeltaCubicMetersPerSecond);
        fingerprint.real(
            source.referencePredictedNetOutwardVolumeRateCubicMetersPerSecond);
        fingerprint.real(
            source.shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond);
        fingerprint.real(
            source.predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond);
        fingerprint.real(
            source.referenceContinuityResidualCubicMetersPerSecond);
        fingerprint.real(
            source.shadowContinuityResidualCubicMetersPerSecond);
        fingerprint.real(
            source.continuityResidualDeltaCubicMetersPerSecond);
        fingerprint.real(source.referenceIntegratedSourcePascalsMeters);
        fingerprint.real(source.shadowIntegratedSourcePascalsMeters);
        fingerprint.real(source.integratedSourceDeltaPascalsMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        comparison.componentIntegratedSourceDeltasPascalsMeters.size()));
    for (const double delta
         : comparison.componentIntegratedSourceDeltasPascalsMeters) {
        fingerprint.real(delta);
    }
    const auto& sourceDiagnostics = comparison.sourceDiagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        sourceDiagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        sourceDiagnostics.componentCount));
    fingerprintScalarComparison(
        fingerprint, sourceDiagnostics.geometryVolumeRate);
    fingerprintScalarComparison(
        fingerprint, sourceDiagnostics.predictedNetOutwardVolumeRate);
    fingerprintScalarComparison(
        fingerprint, sourceDiagnostics.continuityResidual);
    fingerprintScalarComparison(
        fingerprint, sourceDiagnostics.integratedSource);
    fingerprint.real(sourceDiagnostics
        .maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters);
    fingerprint.integer(static_cast<std::uint8_t>(
        sourceDiagnostics.finite));
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
    fingerprint.real(
        diagnostics.pressureDifferenceDotProductPascalsSquared);
    fingerprint.real(diagnostics.bestFitShadowPressureScale);
    fingerprint.real(diagnostics.pressureDifferenceCosineSimilarity);
    fingerprint.real(
        diagnostics.bestFitPressureShapeResidualL2Pascals);
    fingerprint.real(
        diagnostics.relativeBestFitPressureShapeResidualL2);
    fingerprint.real(diagnostics.nodalForceDeltaL2Newtons);
    fingerprint.real(diagnostics.maximumNodalForceDeltaNewtons);
    fingerprint.real(diagnostics.nodalForceDotProductNewtonsSquared);
    fingerprint.real(diagnostics.bestFitShadowNodalForceScale);
    fingerprint.real(diagnostics.nodalForceCosineSimilarity);
    fingerprint.real(
        diagnostics.bestFitNodalForceShapeResidualL2Newtons);
    fingerprint.real(
        diagnostics.relativeBestFitNodalForceShapeResidualL2);
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
    CompensatedSum pressureDotProduct;
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
        pressureDotProduct.add(
            sample.referencePressureDifferencePascals
            * sample.shadowPressureDifferencePascals);
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
    result.pressureDifferenceDotProductPascalsSquared =
        pressureDotProduct.value();
    const double referencePressureSquaredValue =
        referencePressureSquared.value();
    result.bestFitShadowPressureScale = referencePressureSquaredValue > 0.0
        ? result.pressureDifferenceDotProductPascalsSquared
            / referencePressureSquaredValue
        : 0.0;
    const double pressureNormProduct =
        result.referencePressureDifferenceL2Pascals
        * result.shadowPressureDifferenceL2Pascals;
    result.pressureDifferenceCosineSimilarity = pressureNormProduct > 0.0
        ? result.pressureDifferenceDotProductPascalsSquared
            / pressureNormProduct
        : (result.referencePressureDifferenceL2Pascals == 0.0
               && result.shadowPressureDifferenceL2Pascals == 0.0
            ? 1.0 : 0.0);
    CompensatedSum pressureShapeResidualSquared;
    for (const auto& sample : comparison.samples) {
        const double residual = sample.shadowPressureDifferencePascals
            - result.bestFitShadowPressureScale
                * sample.referencePressureDifferencePascals;
        pressureShapeResidualSquared.add(residual * residual);
    }
    result.bestFitPressureShapeResidualL2Pascals = std::sqrt(
        std::max(0.0, pressureShapeResidualSquared.value()));
    result.relativeBestFitPressureShapeResidualL2 = relativeDifference(
        result.bestFitPressureShapeResidualL2Pascals,
        result.referencePressureDifferenceL2Pascals,
        result.shadowPressureDifferenceL2Pascals);

    CompensatedSum nodalDeltaSquared;
    CompensatedSum referenceNodalSquared;
    CompensatedSum shadowNodalSquared;
    CompensatedSum nodalDotProduct;
    for (const auto& node : comparison.nodes) {
        const double delta = norm(node.shadowMinusReferenceNewtons);
        nodalDeltaSquared.add(delta * delta);
        referenceNodalSquared.add(
            node.referenceForceNewtons.x * node.referenceForceNewtons.x
            + node.referenceForceNewtons.y * node.referenceForceNewtons.y
            + node.referenceForceNewtons.z * node.referenceForceNewtons.z);
        shadowNodalSquared.add(
            node.shadowForceNewtons.x * node.shadowForceNewtons.x
            + node.shadowForceNewtons.y * node.shadowForceNewtons.y
            + node.shadowForceNewtons.z * node.shadowForceNewtons.z);
        nodalDotProduct.add(
            node.referenceForceNewtons.x * node.shadowForceNewtons.x
            + node.referenceForceNewtons.y * node.shadowForceNewtons.y
            + node.referenceForceNewtons.z * node.shadowForceNewtons.z);
        result.maximumNodalForceDeltaNewtons = std::max(
            result.maximumNodalForceDeltaNewtons, delta);
    }
    result.nodalForceDeltaL2Newtons = std::sqrt(
        std::max(0.0, nodalDeltaSquared.value()));
    result.nodalForceDotProductNewtonsSquared = nodalDotProduct.value();
    result.bestFitShadowNodalForceScale = referenceNodalSquared.value() > 0.0
        ? result.nodalForceDotProductNewtonsSquared
            / referenceNodalSquared.value()
        : 0.0;
    const double referenceNodalL2 = std::sqrt(
        std::max(0.0, referenceNodalSquared.value()));
    const double shadowNodalL2 = std::sqrt(
        std::max(0.0, shadowNodalSquared.value()));
    const double nodalNormProduct = referenceNodalL2 * shadowNodalL2;
    result.nodalForceCosineSimilarity = nodalNormProduct > 0.0
        ? result.nodalForceDotProductNewtonsSquared / nodalNormProduct
        : (referenceNodalL2 == 0.0 && shadowNodalL2 == 0.0 ? 1.0 : 0.0);
    CompensatedSum nodalShapeResidualSquared;
    for (const auto& node : comparison.nodes) {
        const StructureVector3 residual{
            node.shadowForceNewtons.x
                - result.bestFitShadowNodalForceScale
                    * node.referenceForceNewtons.x,
            node.shadowForceNewtons.y
                - result.bestFitShadowNodalForceScale
                    * node.referenceForceNewtons.y,
            node.shadowForceNewtons.z
                - result.bestFitShadowNodalForceScale
                    * node.referenceForceNewtons.z,
        };
        const double residualNorm = norm(residual);
        nodalShapeResidualSquared.add(residualNorm * residualNorm);
    }
    result.bestFitNodalForceShapeResidualL2Newtons = std::sqrt(
        std::max(0.0, nodalShapeResidualSquared.value()));
    result.relativeBestFitNodalForceShapeResidualL2 = relativeDifference(
        result.bestFitNodalForceShapeResidualL2Newtons,
        referenceNodalL2, shadowNodalL2);
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
        && std::isfinite(
            result.pressureDifferenceDotProductPascalsSquared)
        && std::isfinite(result.bestFitShadowPressureScale)
        && std::isfinite(result.pressureDifferenceCosineSimilarity)
        && std::isfinite(result.bestFitPressureShapeResidualL2Pascals)
        && std::isfinite(
            result.relativeBestFitPressureShapeResidualL2)
        && std::isfinite(result.nodalForceDeltaL2Newtons)
        && std::isfinite(result.maximumNodalForceDeltaNewtons)
        && std::isfinite(result.nodalForceDotProductNewtonsSquared)
        && std::isfinite(result.bestFitShadowNodalForceScale)
        && std::isfinite(result.nodalForceCosineSimilarity)
        && std::isfinite(result.bestFitNodalForceShapeResidualL2Newtons)
        && std::isfinite(result.relativeBestFitNodalForceShapeResidualL2)
        && finiteVector(result.shadowMinusReferenceForceNewtons)
        && std::isfinite(result.forceDeltaNormNewtons)
        && std::isfinite(result.relativeForceDelta)
        && finiteVector(result.shadowMinusReferenceMomentNewtonMeters)
        && std::isfinite(result.momentDeltaNormNewtonMeters)
        && std::isfinite(result.relativeMomentDelta)
        && std::isfinite(result.shadowMinusReferencePowerWatts);
    return result;
}

SceneFluidPressureShadowScalarComparisonDiagnostics summarizeSourceScalar(
    const std::vector<SceneFluidPressureShadowControlSourceComparison>&
        controls,
    const double SceneFluidPressureShadowControlSourceComparison::*reference,
    const double SceneFluidPressureShadowControlSourceComparison::*shadow,
    const double SceneFluidPressureShadowControlSourceComparison::*delta) {
    SceneFluidPressureShadowScalarComparisonDiagnostics result;
    CompensatedSum referenceSquared;
    CompensatedSum shadowSquared;
    CompensatedSum deltaSquared;
    result.exact = true;
    for (const auto& control : controls) {
        const double referenceValue = control.*reference;
        const double shadowValue = control.*shadow;
        const double deltaValue = control.*delta;
        referenceSquared.add(referenceValue * referenceValue);
        shadowSquared.add(shadowValue * shadowValue);
        deltaSquared.add(deltaValue * deltaValue);
        result.maximumAbsoluteDelta = std::max(
            result.maximumAbsoluteDelta, std::abs(deltaValue));
        result.exact = result.exact && deltaValue == 0.0;
    }
    result.referenceL2 = std::sqrt(std::max(0.0, referenceSquared.value()));
    result.shadowL2 = std::sqrt(std::max(0.0, shadowSquared.value()));
    result.deltaL2 = std::sqrt(std::max(0.0, deltaSquared.value()));
    result.relativeDeltaL2 = relativeDifference(
        result.deltaL2, result.referenceL2, result.shadowL2);
    result.finite = std::isfinite(result.referenceL2)
        && std::isfinite(result.shadowL2)
        && std::isfinite(result.deltaL2)
        && std::isfinite(result.maximumAbsoluteDelta)
        && std::isfinite(result.relativeDeltaL2);
    return result;
}

SceneFluidPressureShadowSourceComparisonDiagnostics summarizeSources(
    const SceneFluidPressureShadowComparison& comparison) {
    SceneFluidPressureShadowSourceComparisonDiagnostics result;
    result.controlVolumeCount = comparison.controlSources.size();
    result.componentCount = comparison
        .componentIntegratedSourceDeltasPascalsMeters.size();
    result.geometryVolumeRate = summarizeSourceScalar(
        comparison.controlSources,
        &SceneFluidPressureShadowControlSourceComparison::
            referenceGeometryVolumeRateCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            shadowGeometryVolumeRateCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            geometryVolumeRateDeltaCubicMetersPerSecond);
    result.predictedNetOutwardVolumeRate = summarizeSourceScalar(
        comparison.controlSources,
        &SceneFluidPressureShadowControlSourceComparison::
            referencePredictedNetOutwardVolumeRateCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond);
    result.continuityResidual = summarizeSourceScalar(
        comparison.controlSources,
        &SceneFluidPressureShadowControlSourceComparison::
            referenceContinuityResidualCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            shadowContinuityResidualCubicMetersPerSecond,
        &SceneFluidPressureShadowControlSourceComparison::
            continuityResidualDeltaCubicMetersPerSecond);
    result.integratedSource = summarizeSourceScalar(
        comparison.controlSources,
        &SceneFluidPressureShadowControlSourceComparison::
            referenceIntegratedSourcePascalsMeters,
        &SceneFluidPressureShadowControlSourceComparison::
            shadowIntegratedSourcePascalsMeters,
        &SceneFluidPressureShadowControlSourceComparison::
            integratedSourceDeltaPascalsMeters);
    for (const double delta
         : comparison.componentIntegratedSourceDeltasPascalsMeters) {
        result.maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters =
            std::max(
                result
                    .maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters,
                std::abs(delta));
    }
    result.finite = result.geometryVolumeRate.finite
        && result.predictedNetOutwardVolumeRate.finite
        && result.continuityResidual.finite
        && result.integratedSource.finite
        && std::isfinite(result
            .maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters);
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

SceneFluidPressureShadowComparison attachSourceComparison(
    SceneFluidPressureShadowComparison comparison,
    const SceneFluidPressureProjection& referenceProjection,
    const SceneFluidMimeticControlCellSet& shadowControlCells,
    const SceneFluidMimeticPressureSourceSet& shadowSources,
    const SceneFluidPressureShadowComparisonLimits& limits) {
    validateSceneFluidPressureProjectionIntegrity(referenceProjection);
    validateSceneFluidMimeticControlCellIntegrity(shadowControlCells);
    validateSceneFluidMimeticPressureSources(
        shadowSources, shadowControlCells);
    const std::size_t controlCount = referenceProjection.controlVolumes.size();
    const std::size_t componentCount = shadowSources.componentCount;
    if (controlCount > limits.maximumControlVolumes
        || componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid pressure-shadow source comparison exceeds its count limit");
    }
    if (referenceProjection.fingerprint == 0
        || shadowControlCells.fingerprint == 0
        || shadowSources.fingerprint == 0) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow source identity is invalid");
    }
    if (referenceProjection.acceptedStepCount
            != comparison.acceptedStepCount
        || shadowControlCells.acceptedStepCount
            != comparison.acceptedStepCount
        || shadowSources.acceptedStepCount != comparison.acceptedStepCount
        || referenceProjection.simulationTimeSeconds
            != comparison.simulationTimeSeconds
        || shadowControlCells.simulationTimeSeconds
            != comparison.simulationTimeSeconds
        || shadowSources.simulationTimeSeconds
            != comparison.simulationTimeSeconds
        || referenceProjection.pressureControlVolumeFingerprint
            != shadowControlCells.pressureControlVolumeFingerprint
        || shadowSources.mimeticControlCellFingerprint
            != shadowControlCells.fingerprint
        || shadowControlCells.structureDefinitionFingerprint
            != comparison.structureDefinitionFingerprint
        || shadowSources.structureDefinitionFingerprint
            != comparison.structureDefinitionFingerprint
        || referenceProjection.settings.densityKgPerCubicMeter
            != shadowSources.settings.densityKgPerCubicMeter
        || referenceProjection.settings.timeStepSeconds
            != shadowSources.settings.timeStepSeconds
        || controlCount != shadowControlCells.controlCells.size()
        || controlCount != shadowSources.controls.size()
        || componentCount == 0) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow source binding is invalid");
    }

    comparison.includesSourceComparison = true;
    comparison.referencePressureProjectionFingerprint =
        referenceProjection.fingerprint;
    comparison.shadowControlCellFingerprint = shadowControlCells.fingerprint;
    comparison.shadowPressureSourceFingerprint = shadowSources.fingerprint;
    comparison.controlSources.reserve(controlCount);
    std::vector<CompensatedSum> componentDeltas(componentCount);
    for (std::size_t index = 0; index < controlCount; ++index) {
        const auto& reference = referenceProjection.controlVolumes[index];
        const auto& shadowCell = shadowControlCells.controlCells[index];
        const auto& shadow = shadowSources.controls[index];
        if (reference.controlVolumeIndex != index
            || shadowCell.controlCellIndex != index
            || shadowCell.controlVolumeIndex != index
            || shadow.controlCellIndex != index
            || shadow.controlVolumeIndex != index
            || reference.stableId != shadowCell.stableId
            || reference.stableId != shadow.stableId
            || reference.componentIndex != shadowCell.componentIndex
            || reference.componentIndex != shadow.componentIndex
            || reference.componentIndex >= componentCount) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow control source is misbound");
        }
        SceneFluidPressureShadowControlSourceComparison record;
        record.controlVolumeIndex = index;
        record.stableId = reference.stableId;
        record.componentIndex = reference.componentIndex;
        record.referenceGeometryVolumeRateCubicMetersPerSecond =
            reference.geometryVolumeChangeRateCubicMetersPerSecond;
        record.shadowGeometryVolumeRateCubicMetersPerSecond =
            shadow.geometryVolumeChangeRateCubicMetersPerSecond;
        record.geometryVolumeRateDeltaCubicMetersPerSecond =
            record.shadowGeometryVolumeRateCubicMetersPerSecond
            - record.referenceGeometryVolumeRateCubicMetersPerSecond;
        record.referencePredictedNetOutwardVolumeRateCubicMetersPerSecond =
            reference.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        record.shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond =
            shadow.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        record.predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond =
            record.shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond
            - record
                .referencePredictedNetOutwardVolumeRateCubicMetersPerSecond;
        record.referenceContinuityResidualCubicMetersPerSecond =
            reference.predictedContinuityResidualCubicMetersPerSecond;
        record.shadowContinuityResidualCubicMetersPerSecond =
            shadow.predictedContinuityResidualCubicMetersPerSecond;
        record.continuityResidualDeltaCubicMetersPerSecond =
            record.shadowContinuityResidualCubicMetersPerSecond
            - record.referenceContinuityResidualCubicMetersPerSecond;
        record.referenceIntegratedSourcePascalsMeters =
            reference.integratedRightHandSidePascalsMeters;
        record.shadowIntegratedSourcePascalsMeters =
            shadow.integratedSourcePascalsMeters;
        record.integratedSourceDeltaPascalsMeters =
            record.shadowIntegratedSourcePascalsMeters
            - record.referenceIntegratedSourcePascalsMeters;
        componentDeltas[record.componentIndex].add(
            record.integratedSourceDeltaPascalsMeters);
        comparison.controlSources.push_back(record);
    }
    comparison.componentIntegratedSourceDeltasPascalsMeters.resize(
        componentCount);
    for (std::size_t component = 0; component < componentCount; ++component) {
        comparison.componentIntegratedSourceDeltasPascalsMeters[component] =
            componentDeltas[component].value();
    }
    comparison.sourceDiagnostics = summarizeSources(comparison);
    if (!comparison.sourceDiagnostics.finite) {
        throw std::overflow_error(
            "scene fluid pressure-shadow source comparison is non-finite");
    }
    comparison.ownedStorageBytes = storageBytes(comparison);
    if (comparison.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid pressure-shadow comparison exceeds its byte limit");
    }
    comparison.fingerprint = productFingerprint(comparison);
    validateSceneFluidPressureShadowComparisonIntegrity(comparison);
    return comparison;
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
    const SceneFluidPressureProjection& referenceProjection,
    const SceneFluidPressureSampleSet& referenceSamples,
    const SceneFluidMimeticControlCellSet& shadowControlCells,
    const SceneFluidMimeticPressureSourceSet& shadowSources,
    const SceneFluidMimeticPressureSampleSet& shadowSamples,
    const ConservativeTransferSettings& transferSettings,
    const SceneFluidPressureShadowComparisonLimits& limits) {
    if (referenceSamples.pressureProjectionFingerprint
        != referenceProjection.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow reference projection is foreign");
    }
    return attachSourceComparison(
        compareSceneFluidPressureShadow(
            surface, state, transfer, quadrature, referenceSamples,
            shadowSamples, transferSettings, limits),
        referenceProjection, shadowControlCells, shadowSources, limits);
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
    if (comparison.includesSourceComparison) {
        if (comparison.referencePressureProjectionFingerprint == 0
            || comparison.shadowControlCellFingerprint == 0
            || comparison.shadowPressureSourceFingerprint == 0
            || comparison.controlSources.empty()
            || comparison
                .componentIntegratedSourceDeltasPascalsMeters.empty()) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow source comparison is incomplete");
        }
        std::vector<CompensatedSum> componentDeltas(
            comparison.componentIntegratedSourceDeltasPascalsMeters.size());
        for (std::size_t index = 0;
             index < comparison.controlSources.size(); ++index) {
            const auto& source = comparison.controlSources[index];
            const bool finiteFields =
                std::isfinite(source
                    .referenceGeometryVolumeRateCubicMetersPerSecond)
                && std::isfinite(source
                    .shadowGeometryVolumeRateCubicMetersPerSecond)
                && std::isfinite(source
                    .geometryVolumeRateDeltaCubicMetersPerSecond)
                && std::isfinite(source
                    .referencePredictedNetOutwardVolumeRateCubicMetersPerSecond)
                && std::isfinite(source
                    .shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond)
                && std::isfinite(source
                    .predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond)
                && std::isfinite(source
                    .referenceContinuityResidualCubicMetersPerSecond)
                && std::isfinite(source
                    .shadowContinuityResidualCubicMetersPerSecond)
                && std::isfinite(source
                    .continuityResidualDeltaCubicMetersPerSecond)
                && std::isfinite(
                    source.referenceIntegratedSourcePascalsMeters)
                && std::isfinite(
                    source.shadowIntegratedSourcePascalsMeters)
                && std::isfinite(
                    source.integratedSourceDeltaPascalsMeters);
            if (source.controlVolumeIndex != index || source.stableId == 0
                || source.componentIndex >= componentDeltas.size()
                || !finiteFields
                || source.geometryVolumeRateDeltaCubicMetersPerSecond
                    != source
                            .shadowGeometryVolumeRateCubicMetersPerSecond
                        - source
                            .referenceGeometryVolumeRateCubicMetersPerSecond
                || source
                        .predictedNetOutwardVolumeRateDeltaCubicMetersPerSecond
                    != source
                            .shadowPredictedNetOutwardVolumeRateCubicMetersPerSecond
                        - source
                            .referencePredictedNetOutwardVolumeRateCubicMetersPerSecond
                || source.continuityResidualDeltaCubicMetersPerSecond
                    != source
                            .shadowContinuityResidualCubicMetersPerSecond
                        - source
                            .referenceContinuityResidualCubicMetersPerSecond
                || source.integratedSourceDeltaPascalsMeters
                    != source.shadowIntegratedSourcePascalsMeters
                        - source.referenceIntegratedSourcePascalsMeters) {
                throw std::invalid_argument(
                    "scene fluid pressure-shadow control source is invalid");
            }
            componentDeltas[source.componentIndex].add(
                source.integratedSourceDeltaPascalsMeters);
        }
        for (std::size_t component = 0;
             component < componentDeltas.size(); ++component) {
            if (comparison
                    .componentIntegratedSourceDeltasPascalsMeters[component]
                != componentDeltas[component].value()) {
                throw std::invalid_argument(
                    "scene fluid pressure-shadow component source is invalid");
            }
        }
        if (!comparison.sourceDiagnostics.finite
            || comparison.sourceDiagnostics != summarizeSources(comparison)) {
            throw std::invalid_argument(
                "scene fluid pressure-shadow source summary is invalid");
        }
    } else if (comparison.referencePressureProjectionFingerprint != 0
               || comparison.shadowControlCellFingerprint != 0
               || comparison.shadowPressureSourceFingerprint != 0
               || !comparison.controlSources.empty()
               || !comparison
                    .componentIntegratedSourceDeltasPascalsMeters.empty()
               || comparison.sourceDiagnostics
                    != SceneFluidPressureShadowSourceComparisonDiagnostics{}) {
        throw std::invalid_argument(
            "scene fluid pressure-shadow source comparison is unexpected");
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
