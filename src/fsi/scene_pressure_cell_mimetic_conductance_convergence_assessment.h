#pragma once

#include "scene_pressure_cell_mimetic_conductance_phase_refinement_audit.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    scenePressureCellMimeticConductanceConvergenceAssessmentVersion = 1;

enum class ScenePressureCellMimeticConductanceConvergenceOutcome
    : std::uint8_t {
    InsufficientEvidence = 0,
    ContinuumTrendCandidate = 1,
};

enum class ScenePressureCellMimeticConductanceConvergenceRejection
    : std::uint64_t {
    None = 0,
    IncompletePhasePopulation = 1ULL << 0U,
    NonGeometricRefinement = 1ULL << 1U,
    MeanIncrementDirectionChange = 1ULL << 2U,
    MeanIncrementNotContracting = 1ULL << 3U,
    ApparentOrderOutsidePolicy = 1ULL << 4U,
    ExtrapolationGapOutsidePolicy = 1ULL << 5U,
    FinePhaseVariationOutsidePolicy = 1ULL << 6U,
    PhaseVariationNotContracting = 1ULL << 7U,
    PhaseIncrementDirectionChange = 1ULL << 8U,
    PhaseIncrementNotContracting = 1ULL << 9U,
};

struct ScenePressureCellMimeticConductanceConvergencePolicy {
    double refinementRatioTolerance = 1.0e-12;
    double maximumMeanIncrementContractionRatio = 0.75;
    double minimumApparentOrder = 0.25;
    double maximumApparentOrder = 4.0;
    double maximumRelativeExtrapolationGap = 0.20;
    double maximumFinePhaseCoefficientOfVariation = 0.10;
    double maximumPhaseVariationContractionRatio = 0.75;
    bool requireConsistentPhaseIncrementDirection = true;
    bool requirePhaseIncrementContraction = true;
    double maximumPhaseIncrementContractionRatio = 0.75;

    bool operator==(
        const ScenePressureCellMimeticConductanceConvergencePolicy&) const =
        default;
};

struct ScenePressureCellMimeticConductancePhaseTrajectory {
    std::size_t phaseIndex = 0;
    fluid::Vector3 gridPhaseFraction;
    bool complete = false;
    std::array<double, 3> normalizedConductance{};
    double previousIncrement = 0.0;
    double latestIncrement = 0.0;
    double incrementContractionRatio = 0.0;
    std::uint64_t rejectionMask = 0;

    bool operator==(
        const ScenePressureCellMimeticConductancePhaseTrajectory&) const =
        default;
};

struct ScenePressureCellMimeticConductanceConvergenceAssessmentLimits {
    std::size_t maximumPhaseTrajectories = 32;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL;
};

// Read-only three-level continuum-trend screen over the finest level of three
// independently immutable phase/refinement audits. Passing this policy names
// only a trend candidate; it is neither a convergence proof nor a pressure-
// owner/load-application command.
struct ScenePressureCellMimeticConductanceConvergenceAssessment {
    std::uint32_t version =
        scenePressureCellMimeticConductanceConvergenceAssessmentVersion;
    std::uint64_t fingerprint = 0;
    std::array<std::uint64_t, 3> sourceAuditFingerprints{};
    std::array<std::size_t, 3> sourceLevelIndices{};
    std::uint64_t policyFingerprint = 0;
    ScenePressureCellMimeticConductanceConvergencePolicy policy;
    ScenePressureCellMimeticConductanceConvergenceOutcome outcome =
        ScenePressureCellMimeticConductanceConvergenceOutcome::
            InsufficientEvidence;
    std::uint64_t rejectionMask = 0;
    std::uint32_t rejectionCount = 0;
    std::size_t ownedStorageBytes = 0;
    double refinementRatio = 0.0;
    std::array<double, 3> meanNormalizedConductance{};
    double previousMeanIncrement = 0.0;
    double latestMeanIncrement = 0.0;
    double meanIncrementContractionRatio = 0.0;
    double apparentOrder = 0.0;
    double extrapolatedNormalizedConductance = 0.0;
    double relativeFineToExtrapolatedGap = 0.0;
    std::array<double, 3> phaseCoefficientOfVariation{};
    double latestPhaseVariationContractionRatio = 0.0;
    std::size_t phaseDirectionChangeCount = 0;
    std::size_t phaseNonContractingCount = 0;
    std::vector<ScenePressureCellMimeticConductancePhaseTrajectory>
        phaseTrajectories;

    bool operator==(
        const ScenePressureCellMimeticConductanceConvergenceAssessment&)
        const = default;
};

[[nodiscard]] std::uint64_t
scenePressureCellMimeticConductanceConvergencePolicyFingerprint(
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy);

[[nodiscard]] ScenePressureCellMimeticConductanceConvergenceAssessment
assessScenePressureCellMimeticConductanceConvergence(
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& coarse,
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& middle,
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& fine,
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy = {},
    const ScenePressureCellMimeticConductanceConvergenceAssessmentLimits&
        limits = {});

void validateScenePressureCellMimeticConductanceConvergenceAssessmentIntegrity(
    const ScenePressureCellMimeticConductanceConvergenceAssessment& assessment,
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& coarse,
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& middle,
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& fine);

[[nodiscard]] bool
scenePressureCellMimeticConductanceConvergenceRejectedFor(
    const ScenePressureCellMimeticConductanceConvergenceAssessment& assessment,
    ScenePressureCellMimeticConductanceConvergenceRejection rejection) noexcept;

} // namespace simwing::fsi
