#pragma once

#include "scene_fluid_pressure_operator_response_audit.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    scenePressureCellOperatorRefinementAuditVersion = 2;

struct ScenePressureCellOperatorRefinementAuditSettings {
    SceneFluidPressureOperatorResponseAuditSettings response;
    // Signed fraction in [-0.5, 0.5) by which the periodic grid lower corner
    // is moved along each axis. Zero preserves the visible worker's grid.
    fluid::Vector3 gridPhaseFraction;

    bool operator==(
        const ScenePressureCellOperatorRefinementAuditSettings&) const =
        default;
};

struct ScenePressureCellOperatorRefinementAuditLimits {
    std::size_t maximumSamples = 16;
    std::size_t maximumGridCellsPerSample = 2'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
    SceneFluidPressureOperatorResponseAuditLimits response;
};

struct ScenePressureCellOperatorRefinementSample {
    std::size_t sampleIndex = 0;
    fluid::GridCellCounts cellCounts;
    std::size_t gridCellCount = 0;
    fluid::Vector3 gridLowerMeters;
    fluid::Vector3 cellSpacingMeters;
    double intakeAreaSquareMeters = 0.0;
    std::size_t controlVolumeCount = 0;
    std::size_t fullTraceCount = 0;
    std::size_t reducedTraceCount = 0;
    double graphConductanceMeters = 0.0;
    double shadowConductanceMeters = 0.0;
    double graphToShadowConductanceRatio = 0.0;
    double normalizedGraphConductance = 0.0;
    double normalizedShadowConductance = 0.0;
    SceneFluidPressureOperatorResponseAudit response;

    bool operator==(
        const ScenePressureCellOperatorRefinementSample&) const = default;
};

// Offline, rest-geometry refinement discriminator for the analytic intake.
// Every requested isotropic grid rebuilds the complete scene-v2 cut-cell
// graph and mixed-hybrid topology, then runs the immutable response audit.
// Conductance is normalized by h/area for comparison across resolutions.
// No worker state, checkpoint, pressure owner, or Structure load is changed.
struct ScenePressureCellOperatorRefinementAudit {
    std::uint32_t version =
        scenePressureCellOperatorRefinementAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    ScenePressureCellOperatorRefinementAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::vector<ScenePressureCellOperatorRefinementSample> samples;

    bool operator==(
        const ScenePressureCellOperatorRefinementAudit&) const = default;
};

[[nodiscard]] ScenePressureCellOperatorRefinementAudit
auditScenePressureCellOperatorRefinement(
    std::span<const fluid::GridCellCounts> resolutions,
    const ScenePressureCellOperatorRefinementAuditSettings& settings = {},
    const ScenePressureCellOperatorRefinementAuditLimits& limits = {});

void validateScenePressureCellOperatorRefinementAuditIntegrity(
    const ScenePressureCellOperatorRefinementAudit& audit);

} // namespace simwing::fsi
