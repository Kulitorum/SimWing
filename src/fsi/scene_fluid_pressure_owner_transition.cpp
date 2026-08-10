#include "scene_fluid_pressure_owner_transition.h"

#include <bit>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t knownRejectionMask =
    (1ULL << 18U) - 1ULL;

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

void fingerprintPolicy(
    Fingerprint& fingerprint,
    const SceneFluidPressureOwnerTransitionPolicy& policy) {
    fingerprint.integer(static_cast<std::uint8_t>(
        policy.requireSourceComparison));
    fingerprint.real(policy.maximumRelativeSourceDeltaL2);
    fingerprint.real(
        policy.maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters);
    fingerprint.real(policy.maximumRelativePressureDifferenceDeltaL2);
    fingerprint.real(policy.maximumAbsolutePressureScaleDeviation);
    fingerprint.real(policy.maximumRelativePressureShapeResidualL2);
    fingerprint.real(policy.maximumAbsoluteNodalForceScaleDeviation);
    fingerprint.real(policy.maximumRelativeNodalForceShapeResidualL2);
    fingerprint.real(policy.maximumRelativeNetForceDelta);
    fingerprint.real(policy.maximumRelativeNetMomentDelta);
    fingerprint.real(policy.maximumAbsolutePowerDeltaWatts);
    fingerprint.real(policy.maximumTransferForceResidualNewtons);
    fingerprint.real(policy.maximumTransferMomentResidualNewtonMeters);
    fingerprint.real(policy.maximumTransferPowerResidualWatts);
}

std::uint64_t decisionFingerprint(
    const SceneFluidPressureOwnerTransitionDecision& decision) {
    Fingerprint fingerprint;
    fingerprint.integer(decision.version);
    fingerprint.integer(decision.comparisonFingerprint);
    fingerprint.integer(decision.policyFingerprint);
    fingerprintPolicy(fingerprint, decision.policy);
    fingerprint.integer(static_cast<std::uint8_t>(decision.selectedOwner));
    fingerprint.integer(decision.rejectionMask);
    fingerprint.integer(decision.rejectionCount);
    return fingerprint.value();
}

void validateThreshold(const double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

void reject(
    SceneFluidPressureOwnerTransitionDecision& decision,
    const SceneFluidPressureOwnerTransitionRejection rejection) {
    const auto bit = static_cast<std::uint64_t>(rejection);
    if ((decision.rejectionMask & bit) == 0) {
        decision.rejectionMask |= bit;
        ++decision.rejectionCount;
    }
}

bool exceeds(const double value, const double limit) {
    return !std::isfinite(value) || value > limit;
}

bool scaleExceeds(
    const double scale,
    const double deltaL2,
    const double limit) {
    return deltaL2 != 0.0
        && (!std::isfinite(scale) || std::abs(scale - 1.0) > limit);
}

bool transferExceeds(
    const ConservativeTransferDiagnostics& transfer,
    const SceneFluidPressureOwnerTransitionPolicy& policy) {
    return !transfer.finite
        || exceeds(transfer.forceResidualNormNewtons,
                   policy.maximumTransferForceResidualNewtons)
        || exceeds(transfer.momentResidualNormNewtonMeters,
                   policy.maximumTransferMomentResidualNewtonMeters)
        || exceeds(std::abs(transfer.powerResidualWatts),
                   policy.maximumTransferPowerResidualWatts);
}

} // namespace

std::uint64_t sceneFluidPressureOwnerTransitionPolicyFingerprint(
    const SceneFluidPressureOwnerTransitionPolicy& policy) {
    validateSceneFluidPressureOwnerTransitionPolicy(policy);
    Fingerprint fingerprint;
    fingerprintPolicy(fingerprint, policy);
    return fingerprint.value();
}

void validateSceneFluidPressureOwnerTransitionPolicy(
    const SceneFluidPressureOwnerTransitionPolicy& policy) {
    validateThreshold(
        policy.maximumRelativeSourceDeltaL2,
        "pressure-owner source relative tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters,
        "pressure-owner component-source tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumRelativePressureDifferenceDeltaL2,
        "pressure-owner pressure-difference tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumAbsolutePressureScaleDeviation,
        "pressure-owner pressure-scale tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumRelativePressureShapeResidualL2,
        "pressure-owner pressure-shape tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumAbsoluteNodalForceScaleDeviation,
        "pressure-owner nodal-force scale tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumRelativeNodalForceShapeResidualL2,
        "pressure-owner nodal-force shape tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumRelativeNetForceDelta,
        "pressure-owner net-force tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumRelativeNetMomentDelta,
        "pressure-owner net-moment tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumAbsolutePowerDeltaWatts,
        "pressure-owner power tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumTransferForceResidualNewtons,
        "pressure-owner transfer-force tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumTransferMomentResidualNewtonMeters,
        "pressure-owner transfer-moment tolerance must be finite and nonnegative");
    validateThreshold(
        policy.maximumTransferPowerResidualWatts,
        "pressure-owner transfer-power tolerance must be finite and nonnegative");
}

SceneFluidPressureOwnerTransitionDecision
decideSceneFluidPressureOwnerTransition(
    const SceneFluidPressureShadowComparison& comparison,
    const SceneFluidPressureOwnerTransitionPolicy& policy) {
    validateSceneFluidPressureOwnerTransitionPolicy(policy);
    validateSceneFluidPressureShadowComparisonIntegrity(comparison);

    SceneFluidPressureOwnerTransitionDecision decision;
    decision.comparisonFingerprint = comparison.fingerprint;
    decision.policy = policy;
    decision.policyFingerprint =
        sceneFluidPressureOwnerTransitionPolicyFingerprint(policy);

    if (!comparison.diagnostics.finite) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   NonFiniteLoadComparison);
    }

    if (!comparison.includesSourceComparison) {
        if (policy.requireSourceComparison) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       MissingSourceComparison);
        }
    } else {
        const auto& source = comparison.sourceDiagnostics;
        if (!source.finite) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       NonFiniteSourceComparison);
        }
        if (exceeds(source.geometryVolumeRate.relativeDeltaL2,
                    policy.maximumRelativeSourceDeltaL2)) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       GeometryVolumeRateMismatch);
        }
        if (exceeds(source.predictedNetOutwardVolumeRate.relativeDeltaL2,
                    policy.maximumRelativeSourceDeltaL2)) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       PredictedVolumeRateMismatch);
        }
        if (exceeds(source.continuityResidual.relativeDeltaL2,
                    policy.maximumRelativeSourceDeltaL2)) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       ContinuityResidualMismatch);
        }
        if (exceeds(source.integratedSource.relativeDeltaL2,
                    policy.maximumRelativeSourceDeltaL2)) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       IntegratedSourceMismatch);
        }
        if (exceeds(
                source.maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters,
                policy.maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters)) {
            reject(decision,
                   SceneFluidPressureOwnerTransitionRejection::
                       ComponentIntegratedSourceMismatch);
        }
    }

    const auto& diagnostics = comparison.diagnostics;
    if (exceeds(diagnostics.relativePressureDifferenceDeltaL2,
                policy.maximumRelativePressureDifferenceDeltaL2)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   PressureDifferenceMismatch);
    }
    if (scaleExceeds(
            diagnostics.bestFitShadowPressureScale,
            diagnostics.pressureDifferenceDeltaL2Pascals,
            policy.maximumAbsolutePressureScaleDeviation)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   PressureScaleMismatch);
    }
    if (exceeds(diagnostics.relativeBestFitPressureShapeResidualL2,
                policy.maximumRelativePressureShapeResidualL2)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   PressureShapeMismatch);
    }
    if (scaleExceeds(
            diagnostics.bestFitShadowNodalForceScale,
            diagnostics.nodalForceDeltaL2Newtons,
            policy.maximumAbsoluteNodalForceScaleDeviation)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   NodalForceScaleMismatch);
    }
    if (exceeds(diagnostics.relativeBestFitNodalForceShapeResidualL2,
                policy.maximumRelativeNodalForceShapeResidualL2)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   NodalForceShapeMismatch);
    }
    if (exceeds(diagnostics.relativeForceDelta,
                policy.maximumRelativeNetForceDelta)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   NetForceMismatch);
    }
    if (exceeds(diagnostics.relativeMomentDelta,
                policy.maximumRelativeNetMomentDelta)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   NetMomentMismatch);
    }
    if (exceeds(std::abs(diagnostics.shadowMinusReferencePowerWatts),
                policy.maximumAbsolutePowerDeltaWatts)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::PowerMismatch);
    }
    if (transferExceeds(diagnostics.referenceTransfer, policy)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   ReferenceTransferClosure);
    }
    if (transferExceeds(diagnostics.shadowTransfer, policy)) {
        reject(decision,
               SceneFluidPressureOwnerTransitionRejection::
                   ShadowTransferClosure);
    }

    if (decision.rejectionMask == 0) {
        decision.selectedOwner = SceneFluidPressureOwner::ShadowMimetic;
    }
    decision.fingerprint = decisionFingerprint(decision);
    return decision;
}

bool sceneFluidPressureOwnerTransitionRejectedFor(
    const SceneFluidPressureOwnerTransitionDecision& decision,
    const SceneFluidPressureOwnerTransitionRejection rejection) noexcept {
    const auto bit = static_cast<std::uint64_t>(rejection);
    return bit != 0 && (decision.rejectionMask & bit) != 0;
}

void validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
    const SceneFluidPressureOwnerTransitionDecision& decision,
    const SceneFluidPressureShadowComparison& comparison) {
    if (decision.version
        != sceneFluidPressureOwnerTransitionDecisionVersion) {
        throw std::invalid_argument(
            "unsupported scene fluid pressure-owner transition decision version");
    }
    if (decision.comparisonFingerprint == 0
        || decision.comparisonFingerprint != comparison.fingerprint) {
        throw std::invalid_argument(
            "pressure-owner transition comparison provenance mismatch");
    }
    if ((decision.rejectionMask & ~knownRejectionMask) != 0) {
        throw std::invalid_argument(
            "pressure-owner transition has unknown rejection bits");
    }
    const auto expected =
        decideSceneFluidPressureOwnerTransition(comparison, decision.policy);
    if (decision != expected || decision.fingerprint == 0
        || decision.fingerprint != decisionFingerprint(decision)) {
        throw std::invalid_argument(
            "pressure-owner transition decision integrity failure");
    }
}

} // namespace simwing::fsi
