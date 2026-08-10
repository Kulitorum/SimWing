#pragma once

#include "scene_fluid_pressure_shadow_comparison.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureOwnerTransitionDecisionVersion = 1;

enum class SceneFluidPressureOwner : std::uint8_t {
    ReferenceGraph = 0,
    ShadowMimetic = 1,
};

enum class SceneFluidPressureOwnerTransitionRejection : std::uint64_t {
    None = 0,
    MissingSourceComparison = 1ULL << 0U,
    NonFiniteSourceComparison = 1ULL << 1U,
    GeometryVolumeRateMismatch = 1ULL << 2U,
    PredictedVolumeRateMismatch = 1ULL << 3U,
    ContinuityResidualMismatch = 1ULL << 4U,
    IntegratedSourceMismatch = 1ULL << 5U,
    ComponentIntegratedSourceMismatch = 1ULL << 6U,
    NonFiniteLoadComparison = 1ULL << 7U,
    PressureDifferenceMismatch = 1ULL << 8U,
    PressureScaleMismatch = 1ULL << 9U,
    PressureShapeMismatch = 1ULL << 10U,
    NodalForceScaleMismatch = 1ULL << 11U,
    NodalForceShapeMismatch = 1ULL << 12U,
    NetForceMismatch = 1ULL << 13U,
    NetMomentMismatch = 1ULL << 14U,
    PowerMismatch = 1ULL << 15U,
    ReferenceTransferClosure = 1ULL << 16U,
    ShadowTransferClosure = 1ULL << 17U,
};

struct SceneFluidPressureOwnerTransitionPolicy {
    bool requireSourceComparison = true;
    double maximumRelativeSourceDeltaL2 = 1.0e-10;
    double maximumAbsoluteComponentIntegratedSourceDeltaPascalsMeters =
        1.0e-10;
    double maximumRelativePressureDifferenceDeltaL2 = 1.0e-6;
    double maximumAbsolutePressureScaleDeviation = 1.0e-6;
    double maximumRelativePressureShapeResidualL2 = 1.0e-6;
    double maximumAbsoluteNodalForceScaleDeviation = 1.0e-6;
    double maximumRelativeNodalForceShapeResidualL2 = 1.0e-6;
    double maximumRelativeNetForceDelta = 1.0e-6;
    double maximumRelativeNetMomentDelta = 1.0e-6;
    double maximumAbsolutePowerDeltaWatts = 1.0e-8;
    double maximumTransferForceResidualNewtons = 1.0e-8;
    double maximumTransferMomentResidualNewtonMeters = 1.0e-8;
    double maximumTransferPowerResidualWatts = 1.0e-8;

    bool operator==(
        const SceneFluidPressureOwnerTransitionPolicy&) const = default;
};

// Immutable, read-only decision at the existing graph-vs-mimetic comparison
// boundary. Selecting ShadowMimetic only says that the supplied comparison
// satisfies this explicit policy; applying that field remains a separate
// production integration step.
struct SceneFluidPressureOwnerTransitionDecision {
    std::uint32_t version =
        sceneFluidPressureOwnerTransitionDecisionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t comparisonFingerprint = 0;
    std::uint64_t policyFingerprint = 0;
    SceneFluidPressureOwnerTransitionPolicy policy;
    SceneFluidPressureOwner selectedOwner =
        SceneFluidPressureOwner::ReferenceGraph;
    std::uint64_t rejectionMask = 0;
    std::uint32_t rejectionCount = 0;

    bool operator==(
        const SceneFluidPressureOwnerTransitionDecision&) const = default;
};

[[nodiscard]] std::uint64_t
sceneFluidPressureOwnerTransitionPolicyFingerprint(
    const SceneFluidPressureOwnerTransitionPolicy& policy);

[[nodiscard]] SceneFluidPressureOwnerTransitionDecision
decideSceneFluidPressureOwnerTransition(
    const SceneFluidPressureShadowComparison& comparison,
    const SceneFluidPressureOwnerTransitionPolicy& policy = {});

[[nodiscard]] bool sceneFluidPressureOwnerTransitionRejectedFor(
    const SceneFluidPressureOwnerTransitionDecision& decision,
    SceneFluidPressureOwnerTransitionRejection rejection) noexcept;

void validateSceneFluidPressureOwnerTransitionPolicy(
    const SceneFluidPressureOwnerTransitionPolicy& policy);

void validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
    const SceneFluidPressureOwnerTransitionDecision& decision,
    const SceneFluidPressureShadowComparison& comparison);

} // namespace simwing::fsi
