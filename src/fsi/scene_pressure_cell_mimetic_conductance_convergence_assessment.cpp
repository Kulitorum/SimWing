#include "scene_pressure_cell_mimetic_conductance_convergence_assessment.h"

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

using Assessment =
    ScenePressureCellMimeticConductanceConvergenceAssessment;
using Audit = ScenePressureCellMimeticConductancePhaseRefinementAudit;
using Rejection =
    ScenePressureCellMimeticConductanceConvergenceRejection;

void validatePolicy(
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy) {
    if (!std::isfinite(policy.refinementRatioTolerance)
        || policy.refinementRatioTolerance < 0.0
        || !std::isfinite(
            policy.maximumMeanIncrementContractionRatio)
        || !(policy.maximumMeanIncrementContractionRatio > 0.0)
        || !(policy.maximumMeanIncrementContractionRatio < 1.0)
        || !std::isfinite(policy.minimumApparentOrder)
        || policy.minimumApparentOrder < 0.0
        || !std::isfinite(policy.maximumApparentOrder)
        || policy.maximumApparentOrder < policy.minimumApparentOrder
        || !std::isfinite(policy.maximumRelativeExtrapolationGap)
        || policy.maximumRelativeExtrapolationGap < 0.0
        || !std::isfinite(
            policy.maximumFinePhaseCoefficientOfVariation)
        || policy.maximumFinePhaseCoefficientOfVariation < 0.0
        || !std::isfinite(
            policy.maximumPhaseVariationContractionRatio)
        || !(policy.maximumPhaseVariationContractionRatio > 0.0)
        || !(policy.maximumPhaseVariationContractionRatio < 1.0)
        || !std::isfinite(
            policy.maximumPhaseIncrementContractionRatio)
        || !(policy.maximumPhaseIncrementContractionRatio > 0.0)
        || !(policy.maximumPhaseIncrementContractionRatio < 1.0)) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic convergence policy is invalid");
    }
}

void fingerprintPolicy(
    Fingerprint& fingerprint,
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy) {
    for (const double value : {
             policy.refinementRatioTolerance,
             policy.maximumMeanIncrementContractionRatio,
             policy.minimumApparentOrder,
             policy.maximumApparentOrder,
             policy.maximumRelativeExtrapolationGap,
             policy.maximumFinePhaseCoefficientOfVariation,
             policy.maximumPhaseVariationContractionRatio,
             policy.maximumPhaseIncrementContractionRatio,
         }) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint8_t>(
        policy.requireConsistentPhaseIncrementDirection));
    fingerprint.integer(static_cast<std::uint8_t>(
        policy.requirePhaseIncrementContraction));
}

void validateSources(const Audit& coarse,
                     const Audit& middle,
                     const Audit& fine) {
    validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        coarse);
    validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        middle);
    validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        fine);
    if (coarse.levels.empty() || middle.levels.empty() || fine.levels.empty()
        || coarse.structureDefinitionFingerprint
            != middle.structureDefinitionFingerprint
        || coarse.structureDefinitionFingerprint
            != fine.structureDefinitionFingerprint
        || coarse.settings != middle.settings
        || coarse.settings != fine.settings
        || coarse.gridPhaseFractions != middle.gridPhaseFractions
        || coarse.gridPhaseFractions != fine.gridPhaseFractions) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic convergence sources are incompatible");
    }
}

bool sameDirection(const double previous, const double latest) {
    return (previous > 0.0 && latest > 0.0)
        || (previous < 0.0 && latest < 0.0)
        || (previous == 0.0 && latest == 0.0);
}

double contractionRatio(const double previous, const double latest) {
    if (previous != 0.0) {
        return std::abs(latest / previous);
    }
    return latest == 0.0
        ? 0.0 : std::numeric_limits<double>::max();
}

void addRejection(std::uint64_t& mask, const Rejection rejection) {
    mask |= static_cast<std::uint64_t>(rejection);
}

std::size_t ownedStorageBytes(const Assessment& assessment) {
    if (assessment.phaseTrajectories.size()
        > std::numeric_limits<std::size_t>::max()
            / sizeof(ScenePressureCellMimeticConductancePhaseTrajectory)) {
        throw std::length_error(
            "scene pressure-cell mimetic convergence storage overflows");
    }
    return assessment.phaseTrajectories.size()
        * sizeof(ScenePressureCellMimeticConductancePhaseTrajectory);
}

std::uint64_t productFingerprint(const Assessment& assessment) {
    Fingerprint fingerprint;
    fingerprint.integer(assessment.version);
    for (const auto value : assessment.sourceAuditFingerprints) {
        fingerprint.integer(value);
    }
    for (const auto value : assessment.sourceLevelIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(assessment.policyFingerprint);
    fingerprintPolicy(fingerprint, assessment.policy);
    fingerprint.integer(static_cast<std::uint8_t>(assessment.outcome));
    fingerprint.integer(assessment.rejectionMask);
    fingerprint.integer(assessment.rejectionCount);
    fingerprint.integer(static_cast<std::uint64_t>(
        assessment.ownedStorageBytes));
    fingerprint.real(assessment.refinementRatio);
    for (const auto value : assessment.meanNormalizedConductance) {
        fingerprint.real(value);
    }
    fingerprint.real(assessment.previousMeanIncrement);
    fingerprint.real(assessment.latestMeanIncrement);
    fingerprint.real(assessment.meanIncrementContractionRatio);
    fingerprint.real(assessment.apparentOrder);
    fingerprint.real(assessment.extrapolatedNormalizedConductance);
    fingerprint.real(assessment.relativeFineToExtrapolatedGap);
    for (const auto value : assessment.phaseCoefficientOfVariation) {
        fingerprint.real(value);
    }
    fingerprint.real(
        assessment.latestPhaseVariationContractionRatio);
    fingerprint.integer(static_cast<std::uint64_t>(
        assessment.phaseDirectionChangeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        assessment.phaseNonContractingCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        assessment.phaseTrajectories.size()));
    for (const auto& trajectory : assessment.phaseTrajectories) {
        fingerprint.integer(static_cast<std::uint64_t>(
            trajectory.phaseIndex));
        fingerprint.real(trajectory.gridPhaseFraction.x);
        fingerprint.real(trajectory.gridPhaseFraction.y);
        fingerprint.real(trajectory.gridPhaseFraction.z);
        fingerprint.integer(static_cast<std::uint8_t>(trajectory.complete));
        for (const auto value : trajectory.normalizedConductance) {
            fingerprint.real(value);
        }
        fingerprint.real(trajectory.previousIncrement);
        fingerprint.real(trajectory.latestIncrement);
        fingerprint.real(trajectory.incrementContractionRatio);
        fingerprint.integer(trajectory.rejectionMask);
    }
    return fingerprint.value();
}

Assessment buildAssessment(
    const Audit& coarse,
    const Audit& middle,
    const Audit& fine,
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy) {
    const std::array<const Audit*, 3> audits{{&coarse, &middle, &fine}};
    std::array<const ScenePressureCellMimeticConductanceRefinementLevel*, 3>
        levels{};
    Assessment result;
    result.policy = policy;
    result.policyFingerprint =
        scenePressureCellMimeticConductanceConvergencePolicyFingerprint(
            policy);
    for (std::size_t index = 0; index < audits.size(); ++index) {
        result.sourceAuditFingerprints[index] = audits[index]->fingerprint;
        result.sourceLevelIndices[index] = audits[index]->levels.size() - 1;
        levels[index] = &audits[index]->levels.back();
        result.meanNormalizedConductance[index] =
            levels[index]->meanNormalizedConductance;
        result.phaseCoefficientOfVariation[index] =
            levels[index]
                ->normalizedConductanceCoefficientOfVariation;
    }

    const double previousRefinementRatio =
        levels[0]->cellSpacingMeters.x / levels[1]->cellSpacingMeters.x;
    const double latestRefinementRatio =
        levels[1]->cellSpacingMeters.x / levels[2]->cellSpacingMeters.x;
    result.refinementRatio = latestRefinementRatio;
    const double ratioScale = std::max({
        1.0, std::abs(previousRefinementRatio),
        std::abs(latestRefinementRatio),
    });
    const bool geometricRefinement =
        std::isfinite(previousRefinementRatio)
        && std::isfinite(latestRefinementRatio)
        && previousRefinementRatio > 1.0
        && latestRefinementRatio > 1.0
        && std::abs(previousRefinementRatio - latestRefinementRatio)
            <= policy.refinementRatioTolerance * ratioScale;
    if (!geometricRefinement) {
        addRejection(
            result.rejectionMask, Rejection::NonGeometricRefinement);
    }

    result.previousMeanIncrement = result.meanNormalizedConductance[1]
        - result.meanNormalizedConductance[0];
    result.latestMeanIncrement = result.meanNormalizedConductance[2]
        - result.meanNormalizedConductance[1];
    const bool meanDirectionConsistent = sameDirection(
        result.previousMeanIncrement, result.latestMeanIncrement);
    if (!meanDirectionConsistent) {
        addRejection(
            result.rejectionMask, Rejection::MeanIncrementDirectionChange);
    }
    result.meanIncrementContractionRatio = contractionRatio(
        result.previousMeanIncrement, result.latestMeanIncrement);
    if (result.meanIncrementContractionRatio
        > policy.maximumMeanIncrementContractionRatio) {
        addRejection(
            result.rejectionMask, Rejection::MeanIncrementNotContracting);
    }

    const bool usableRichardsonTrend = geometricRefinement
        && meanDirectionConsistent
        && result.meanIncrementContractionRatio > 0.0
        && result.meanIncrementContractionRatio < 1.0;
    if (usableRichardsonTrend) {
        result.apparentOrder = std::log(
            1.0 / result.meanIncrementContractionRatio)
            / std::log(result.refinementRatio);
        const double refinementPower = std::pow(
            result.refinementRatio, result.apparentOrder);
        result.extrapolatedNormalizedConductance =
            result.meanNormalizedConductance[2]
            + result.latestMeanIncrement / (refinementPower - 1.0);
        result.relativeFineToExtrapolatedGap = std::abs(
            result.extrapolatedNormalizedConductance
            - result.meanNormalizedConductance[2])
            / std::max(
                std::abs(result.extrapolatedNormalizedConductance),
                std::numeric_limits<double>::min());
    }
    if (!usableRichardsonTrend || !std::isfinite(result.apparentOrder)
        || result.apparentOrder < policy.minimumApparentOrder
        || result.apparentOrder > policy.maximumApparentOrder) {
        addRejection(
            result.rejectionMask, Rejection::ApparentOrderOutsidePolicy);
    }
    if (!usableRichardsonTrend
        || !std::isfinite(result.extrapolatedNormalizedConductance)
        || !(result.extrapolatedNormalizedConductance > 0.0)
        || !std::isfinite(result.relativeFineToExtrapolatedGap)
        || result.relativeFineToExtrapolatedGap
            > policy.maximumRelativeExtrapolationGap) {
        addRejection(
            result.rejectionMask,
            Rejection::ExtrapolationGapOutsidePolicy);
    }

    if (result.phaseCoefficientOfVariation[2]
        > policy.maximumFinePhaseCoefficientOfVariation) {
        addRejection(
            result.rejectionMask,
            Rejection::FinePhaseVariationOutsidePolicy);
    }
    result.latestPhaseVariationContractionRatio = contractionRatio(
        result.phaseCoefficientOfVariation[1],
        result.phaseCoefficientOfVariation[2]);
    if (result.latestPhaseVariationContractionRatio
        > policy.maximumPhaseVariationContractionRatio) {
        addRejection(
            result.rejectionMask,
            Rejection::PhaseVariationNotContracting);
    }

    result.phaseTrajectories.reserve(coarse.gridPhaseFractions.size());
    for (std::size_t phaseIndex = 0;
         phaseIndex < coarse.gridPhaseFractions.size(); ++phaseIndex) {
        ScenePressureCellMimeticConductancePhaseTrajectory trajectory;
        trajectory.phaseIndex = phaseIndex;
        trajectory.gridPhaseFraction = coarse.gridPhaseFractions[phaseIndex];
        trajectory.complete = true;
        for (std::size_t levelIndex = 0; levelIndex < levels.size();
             ++levelIndex) {
            const auto& sample = levels[levelIndex]->samples[phaseIndex];
            if (sample.status
                != ScenePressureCellMimeticConductancePhaseSampleStatus::
                    Accepted) {
                trajectory.complete = false;
                addRejection(
                    result.rejectionMask,
                    Rejection::IncompletePhasePopulation);
                continue;
            }
            trajectory.normalizedConductance[levelIndex] =
                sample.normalizedConductance;
        }
        if (trajectory.complete) {
            trajectory.previousIncrement =
                trajectory.normalizedConductance[1]
                - trajectory.normalizedConductance[0];
            trajectory.latestIncrement =
                trajectory.normalizedConductance[2]
                - trajectory.normalizedConductance[1];
            if (!sameDirection(
                    trajectory.previousIncrement,
                    trajectory.latestIncrement)) {
                ++result.phaseDirectionChangeCount;
                if (policy.requireConsistentPhaseIncrementDirection) {
                    addRejection(
                        trajectory.rejectionMask,
                        Rejection::PhaseIncrementDirectionChange);
                    addRejection(
                        result.rejectionMask,
                        Rejection::PhaseIncrementDirectionChange);
                }
            }
            trajectory.incrementContractionRatio = contractionRatio(
                trajectory.previousIncrement, trajectory.latestIncrement);
            if (trajectory.incrementContractionRatio
                > policy.maximumPhaseIncrementContractionRatio) {
                ++result.phaseNonContractingCount;
                if (policy.requirePhaseIncrementContraction) {
                    addRejection(
                        trajectory.rejectionMask,
                        Rejection::PhaseIncrementNotContracting);
                    addRejection(
                        result.rejectionMask,
                        Rejection::PhaseIncrementNotContracting);
                }
            }
        }
        result.phaseTrajectories.push_back(trajectory);
    }

    result.rejectionCount = static_cast<std::uint32_t>(
        std::popcount(result.rejectionMask));
    result.outcome = result.rejectionMask == 0
        ? ScenePressureCellMimeticConductanceConvergenceOutcome::
            ContinuumTrendCandidate
        : ScenePressureCellMimeticConductanceConvergenceOutcome::
            InsufficientEvidence;
    result.ownedStorageBytes = ownedStorageBytes(result);
    result.fingerprint = productFingerprint(result);
    return result;
}

} // namespace

std::uint64_t
scenePressureCellMimeticConductanceConvergencePolicyFingerprint(
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy) {
    validatePolicy(policy);
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x6d696d636f6e7670ULL});
    fingerprintPolicy(fingerprint, policy);
    return fingerprint.value();
}

ScenePressureCellMimeticConductanceConvergenceAssessment
assessScenePressureCellMimeticConductanceConvergence(
    const Audit& coarse,
    const Audit& middle,
    const Audit& fine,
    const ScenePressureCellMimeticConductanceConvergencePolicy& policy,
    const ScenePressureCellMimeticConductanceConvergenceAssessmentLimits&
        limits) {
    validatePolicy(policy);
    validateSources(coarse, middle, fine);
    if (coarse.gridPhaseFractions.size()
            > limits.maximumPhaseTrajectories
        || coarse.gridPhaseFractions.size()
            > std::numeric_limits<std::size_t>::max()
                / sizeof(
                    ScenePressureCellMimeticConductancePhaseTrajectory)
        || coarse.gridPhaseFractions.size()
                * sizeof(
                    ScenePressureCellMimeticConductancePhaseTrajectory)
            > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene pressure-cell mimetic convergence limit exceeded");
    }
    auto result = buildAssessment(coarse, middle, fine, policy);
    validateScenePressureCellMimeticConductanceConvergenceAssessmentIntegrity(
        result, coarse, middle, fine);
    return result;
}

void validateScenePressureCellMimeticConductanceConvergenceAssessmentIntegrity(
    const Assessment& assessment,
    const Audit& coarse,
    const Audit& middle,
    const Audit& fine) {
    validatePolicy(assessment.policy);
    validateSources(coarse, middle, fine);
    const auto expected = buildAssessment(
        coarse, middle, fine, assessment.policy);
    if (assessment.version
            != scenePressureCellMimeticConductanceConvergenceAssessmentVersion
        || assessment.fingerprint == 0
        || assessment.policyFingerprint == 0
        || assessment != expected) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic convergence assessment integrity is invalid");
    }
}

bool scenePressureCellMimeticConductanceConvergenceRejectedFor(
    const Assessment& assessment,
    const Rejection rejection) noexcept {
    return (assessment.rejectionMask
            & static_cast<std::uint64_t>(rejection)) != 0;
}

} // namespace simwing::fsi
