#pragma once

#include "fluid/grid.h"
#include "scene_fluid_mimetic_region_conductance_audit.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    scenePressureCellMimeticConductancePhaseRefinementAuditVersion = 3;

struct ScenePressureCellMimeticConductancePhaseRefinementAuditSettings {
    // Translates both the canonical diagnostic scene and its periodic grid.
    // This leaves their relative geometry unchanged and exists solely to
    // expose coordinate-origin sensitivity in the offline audit.
    fluid::Vector3 geometryTranslationMeters;
    SceneFluidMimeticRegionConductanceAuditSettings conductance;

    bool operator==(
        const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&)
        const = default;
};

struct ScenePressureCellMimeticConductancePhaseRefinementAuditLimits {
    std::size_t maximumResolutionCount = 16;
    std::size_t maximumPhaseCount = 32;
    std::size_t maximumAggregateSamples = 128;
    std::size_t maximumGridCellsPerSample = 2'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
    SceneFluidMimeticRegionConductanceAuditLimits conductance;
};

enum class ScenePressureCellMimeticConductancePhaseSampleStatus
    : std::uint8_t {
    Accepted = 1,
    RejectedLocalCellLinearConsistency = 2,
};

struct ScenePressureCellMimeticConductancePhaseSample {
    std::size_t sampleIndex = 0;
    std::size_t phaseIndex = 0;
    fluid::Vector3 gridPhaseFraction;
    ScenePressureCellMimeticConductancePhaseSampleStatus status =
        ScenePressureCellMimeticConductancePhaseSampleStatus::Accepted;
    std::size_t gridCellCount = 0;
    fluid::Vector3 gridLowerMeters;
    fluid::Vector3 cellSpacingMeters;
    double intakeAreaSquareMeters = 0.0;
    std::size_t controlVolumeCount = 0;
    std::size_t fullTraceCount = 0;
    std::size_t reducedTraceCount = 0;
    std::size_t openingTraceCount = 0;
    double conductanceMeters = 0.0;
    double normalizedConductance = 0.0;
    std::optional<SceneFluidMimeticTraceLocalCellLinearConsistencyFailure>
        localCellLinearConsistencyRejection;
    std::optional<SceneFluidMimeticRegionConductanceAudit> conductanceAudit;

    bool operator==(
        const ScenePressureCellMimeticConductancePhaseSample&) const =
        default;
};

struct ScenePressureCellMimeticConductanceRefinementLevel {
    std::size_t levelIndex = 0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 cellSpacingMeters;
    std::size_t acceptedSampleCount = 0;
    std::size_t rejectedLocalCellLinearConsistencySampleCount = 0;
    double minimumNormalizedConductance = 0.0;
    double maximumNormalizedConductance = 0.0;
    double meanNormalizedConductance = 0.0;
    double normalizedConductanceCoefficientOfVariation = 0.0;
    std::vector<ScenePressureCellMimeticConductancePhaseSample> samples;

    bool operator==(
        const ScenePressureCellMimeticConductanceRefinementLevel&) const =
        default;
};

// Uncensored shadow-only phase/refinement spectrum. Every requested phase
// builds the pressure face ledger but never constructs or applies the graph
// pressure operator. The mixed-hybrid terminal audit therefore remains
// observable when embedded graph opening ownership is incomplete. This is an
// offline immutable evidence product and never changes worker arithmetic. A
// settings-owned common scene/grid translation supplies a coordinate-origin
// invariance oracle without changing relative geometry.
struct ScenePressureCellMimeticConductancePhaseRefinementAudit {
    std::uint32_t version =
        scenePressureCellMimeticConductancePhaseRefinementAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    ScenePressureCellMimeticConductancePhaseRefinementAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::vector<fluid::Vector3> gridPhaseFractions;
    std::vector<ScenePressureCellMimeticConductanceRefinementLevel> levels;

    bool operator==(
        const ScenePressureCellMimeticConductancePhaseRefinementAudit&) const =
        default;
};

[[nodiscard]] ScenePressureCellMimeticConductancePhaseRefinementAudit
auditScenePressureCellMimeticConductancePhaseRefinement(
    std::span<const fluid::GridCellCounts> resolutions,
    std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&
        settings = {},
    const ScenePressureCellMimeticConductancePhaseRefinementAuditLimits&
        limits = {});

void validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& audit);

} // namespace simwing::fsi
