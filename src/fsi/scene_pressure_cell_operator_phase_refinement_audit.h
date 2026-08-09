#pragma once

#include "scene_pressure_cell_operator_phase_audit.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    scenePressureCellOperatorPhaseRefinementAuditVersion = 1;

struct ScenePressureCellOperatorPhaseRefinementAuditSettings {
    SceneFluidPressureOperatorResponseAuditSettings response;

    bool operator==(
        const ScenePressureCellOperatorPhaseRefinementAuditSettings&) const =
        default;
};

struct ScenePressureCellOperatorPhaseRefinementAuditLimits {
    std::size_t maximumResolutionCount = 16;
    std::size_t maximumPhaseCount = 32;
    std::size_t maximumAggregatePhaseSamples = 128;
    std::size_t maximumGridCellsPerSample = 2'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
    SceneFluidPressureOperatorResponseAuditLimits response;
};

struct ScenePressureCellOperatorPhaseRefinementLevel {
    std::size_t levelIndex = 0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 cellSpacingMeters;
    double acceptedTopologyFraction = 0.0;
    ScenePressureCellOperatorPhaseAudit phaseAudit;

    bool operator==(
        const ScenePressureCellOperatorPhaseRefinementLevel&) const = default;
};

// Phase-averaged refinement discriminator. Every strictly increasing
// isotropic resolution runs the same canonical phase set and owns its complete
// fixed-resolution product, including typed incomplete-face rejections. This
// is an offline evidence product and never changes worker pressure ownership.
struct ScenePressureCellOperatorPhaseRefinementAudit {
    std::uint32_t version =
        scenePressureCellOperatorPhaseRefinementAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    ScenePressureCellOperatorPhaseRefinementAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::vector<fluid::Vector3> gridPhaseFractions;
    std::vector<ScenePressureCellOperatorPhaseRefinementLevel> levels;

    bool operator==(
        const ScenePressureCellOperatorPhaseRefinementAudit&) const = default;
};

[[nodiscard]] ScenePressureCellOperatorPhaseRefinementAudit
auditScenePressureCellOperatorPhaseRefinement(
    std::span<const fluid::GridCellCounts> resolutions,
    std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellOperatorPhaseRefinementAuditSettings& settings = {},
    const ScenePressureCellOperatorPhaseRefinementAuditLimits& limits = {});

void validateScenePressureCellOperatorPhaseRefinementAuditIntegrity(
    const ScenePressureCellOperatorPhaseRefinementAudit& audit);

} // namespace simwing::fsi
