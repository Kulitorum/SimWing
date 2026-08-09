#pragma once

#include "scene_fluid_mimetic_pressure_solve.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticRegionConductanceAuditVersion = 1;

struct SceneFluidMimeticRegionConductanceAuditSettings {
    double terminalIntegratedTransferPascalsMeters = 1.0;
    double absoluteOpeningPairAreaToleranceSquareMeters = 1.0e-12;
    double relativeOpeningPairAreaTolerance = 1.0e-10;
    SceneFluidMimeticTraceSolveSettings solve;

    bool operator==(
        const SceneFluidMimeticRegionConductanceAuditSettings&) const =
        default;
};

struct SceneFluidMimeticRegionConductanceAuditLimits {
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumOpeningTraces = 50'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticRegionConductanceOpening {
    std::size_t openingIndex = 0;
    std::uint64_t traceStableId = 0;
    SceneFluidMimeticHalfFaceKind traceKind =
        SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace;
    std::size_t lowerHalfFaceIndex = 0;
    std::size_t upperHalfFaceIndex = 0;
    std::size_t lowerControlCellIndex = 0;
    std::size_t upperControlCellIndex = 0;
    double areaSquareMeters = 0.0;
    double areaMismatchSquareMeters = 0.0;
    double integratedTransferPascalsMeters = 0.0;

    bool operator==(
        const SceneFluidMimeticRegionConductanceOpening&) const = default;
};

struct SceneFluidMimeticRegionConductanceResponse {
    std::size_t controlCellIndex = 0;
    std::uint64_t stableId = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    double integratedSourcePascalsMeters = 0.0;
    double gaugeAlignedPressurePascals = 0.0;

    bool operator==(
        const SceneFluidMimeticRegionConductanceResponse&) const = default;
};

// Graph-independent two-region Neumann audit. Every permeable cross-region
// trace representing an authored opening is paired by stable identity,
// receives a uniform area-weighted transfer, and drives the mixed-hybrid solve
// directly. A face-aligned opening is a Cartesian trace; an embedded opening
// is an AuthoredOpeningTrace. Material-wall traces are never terminals. Source
// work defines the effective two-terminal conductance. No graph operator or
// Structure load is involved.
struct SceneFluidMimeticRegionConductanceAudit {
    std::uint32_t version =
        sceneFluidMimeticRegionConductanceAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidMimeticRegionConductanceAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t componentCount = 0;
    StableId lowerTerminalRegionId = invalidStableId;
    StableId upperTerminalRegionId = invalidStableId;
    double openingAreaSquareMeters = 0.0;
    double maximumOpeningPairAreaMismatchSquareMeters = 0.0;
    double achievedIntegratedTransferPascalsMeters = 0.0;
    double lowerTerminalIntegratedSourcePascalsMeters = 0.0;
    double upperTerminalIntegratedSourcePascalsMeters = 0.0;
    double componentIntegratedSourcePascalsMeters = 0.0;
    double sourceL2PascalsMeters = 0.0;
    double maximumAbsoluteSourcePascalsMeters = 0.0;
    double pressureL2Pascals = 0.0;
    double sourcePressureWorkPascalsSquaredMeters = 0.0;
    double effectiveTerminalPressureDifferencePascals = 0.0;
    double conductanceMeters = 0.0;
    SceneFluidMimeticPressureSolveDiagnostics solveDiagnostics;
    std::vector<SceneFluidMimeticRegionConductanceOpening> openings;
    std::vector<SceneFluidMimeticRegionConductanceResponse> responses;

    bool operator==(
        const SceneFluidMimeticRegionConductanceAudit&) const = default;
};

[[nodiscard]] SceneFluidMimeticRegionConductanceAudit
auditSceneFluidMimeticRegionConductance(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystem& fullTraceSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedTraceSystem,
    const SceneFluidMimeticRegionConductanceAuditSettings& settings = {},
    const SceneFluidMimeticRegionConductanceAuditLimits& limits = {});

void validateSceneFluidMimeticRegionConductanceAuditIntegrity(
    const SceneFluidMimeticRegionConductanceAudit& audit);

} // namespace simwing::fsi
