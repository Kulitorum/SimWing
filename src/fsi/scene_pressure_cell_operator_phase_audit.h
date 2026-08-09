#pragma once

#include "scene_pressure_cell_operator_refinement_audit.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    scenePressureCellOperatorPhaseAuditVersion = 1;

enum class ScenePressureCellOperatorPhaseSampleStatus : std::uint8_t {
    Accepted,
    RejectedIncompleteFaceOwnership,
};

struct ScenePressureCellOperatorPhaseAuditSettings {
    fluid::GridCellCounts cellCounts{4, 4, 4};
    SceneFluidPressureOperatorResponseAuditSettings response;

    bool operator==(
        const ScenePressureCellOperatorPhaseAuditSettings&) const = default;
};

struct ScenePressureCellOperatorPhaseAuditLimits {
    std::size_t maximumSamples = 32;
    std::size_t maximumGridCellsPerSample = 2'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
    SceneFluidPressureOperatorResponseAuditLimits response;
};

struct ScenePressureCellOperatorPhaseStatistics {
    std::size_t acceptedSampleCount = 0;
    std::size_t rejectedIncompleteFaceOwnershipSampleCount = 0;
    double minimumNormalizedGraphConductance = 0.0;
    double maximumNormalizedGraphConductance = 0.0;
    double meanNormalizedGraphConductance = 0.0;
    double graphCoefficientOfVariation = 0.0;
    double minimumNormalizedShadowConductance = 0.0;
    double maximumNormalizedShadowConductance = 0.0;
    double meanNormalizedShadowConductance = 0.0;
    double shadowCoefficientOfVariation = 0.0;

    bool operator==(
        const ScenePressureCellOperatorPhaseStatistics&) const = default;
};

struct ScenePressureCellOperatorPhaseSample {
    std::size_t sampleIndex = 0;
    fluid::Vector3 gridPhaseFraction;
    ScenePressureCellOperatorPhaseSampleStatus status =
        ScenePressureCellOperatorPhaseSampleStatus::Accepted;
    std::optional<SceneFluidPressureIncompleteFaceOwnership>
        faceOwnershipRejection;
    std::optional<ScenePressureCellOperatorRefinementAudit> acceptedAudit;

    bool operator==(
        const ScenePressureCellOperatorPhaseSample&) const = default;
};

// Fixed-resolution, fixed-geometry placement discriminator. Each phase moves
// the Cartesian grid by a fraction of one cell, rebuilds the complete cut-cell
// and pressure products, and either owns a validated one-sample response audit
// or typed incomplete-face diagnostics. It never enters worker state.
struct ScenePressureCellOperatorPhaseAudit {
    std::uint32_t version = scenePressureCellOperatorPhaseAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    ScenePressureCellOperatorPhaseAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    ScenePressureCellOperatorPhaseStatistics statistics;
    std::vector<ScenePressureCellOperatorPhaseSample> samples;

    bool operator==(
        const ScenePressureCellOperatorPhaseAudit&) const = default;
};

[[nodiscard]] ScenePressureCellOperatorPhaseAudit
auditScenePressureCellOperatorGridPhases(
    std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellOperatorPhaseAuditSettings& settings = {},
    const ScenePressureCellOperatorPhaseAuditLimits& limits = {});

void validateScenePressureCellOperatorPhaseAuditIntegrity(
    const ScenePressureCellOperatorPhaseAudit& audit);

} // namespace simwing::fsi
