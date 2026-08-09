#pragma once

#include "scene_fluid_mimetic_pressure_solve.h"
#include "scene_fluid_pressure_operator.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidPressureOperatorResponseAuditVersion = 3;

enum class SceneFluidPressureOperatorResponseModeKind : std::uint8_t {
    AcceptedSource = 1,
    CoordinateX = 2,
    CoordinateY = 3,
    CoordinateZ = 4,
    MixedCoordinate = 5,
    StableIdPattern = 6,
    RegionContrast = 7,
};

struct SceneFluidPressureOperatorResponseAuditSettings {
    double manufacturedPressureL2Pascals = 100.0;
    SceneFluidPressureSolveSettings graphSolve;
    SceneFluidMimeticTraceSolveSettings shadowSolve;

    bool operator==(
        const SceneFluidPressureOperatorResponseAuditSettings&) const =
        default;
};

struct SceneFluidPressureOperatorResponseAuditLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumModes = 64;
    std::size_t maximumResponseRecords = 200'000'000;
    std::size_t maximumOwnedBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPressureOperatorResponseRecord {
    std::size_t modeIndex = 0;
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    double integratedSourcePascalsMeters = 0.0;
    double graphGaugeAlignedPressurePascals = 0.0;
    double shadowGaugeAlignedPressurePascals = 0.0;
    double shadowMinusBestFitGraphPressurePascals = 0.0;

    bool operator==(
        const SceneFluidPressureOperatorResponseRecord&) const = default;
};

struct SceneFluidPressureOperatorResponseModeDiagnostics {
    std::size_t modeIndex = 0;
    SceneFluidPressureOperatorResponseModeKind kind =
        SceneFluidPressureOperatorResponseModeKind::AcceptedSource;
    std::size_t firstResponse = 0;
    std::size_t responseCount = 0;
    double sourceL2PascalsMeters = 0.0;
    double maximumAbsoluteSourcePascalsMeters = 0.0;
    double graphPressureL2Pascals = 0.0;
    double shadowPressureL2Pascals = 0.0;
    double pressureDotProductPascalsSquared = 0.0;
    double graphSourcePressureWorkPascalsSquaredMeters = 0.0;
    double shadowSourcePressureWorkPascalsSquaredMeters = 0.0;
    double shadowToGraphSourceComplianceRatio = 0.0;
    double bestFitShadowPressureScale = 0.0;
    double pressureCosineSimilarity = 0.0;
    double bestFitShapeResidualL2Pascals = 0.0;
    double relativeBestFitShapeResidualL2 = 0.0;
    double maximumAbsoluteBestFitShapeResidualPascals = 0.0;
    std::size_t graphIterationCount = 0;
    double graphFinalResidualL2PascalsMeters = 0.0;
    double graphFinalResidualMaximumPascalsMeters = 0.0;
    std::size_t shadowIterationCount = 0;
    double shadowFinalResidualL2PascalsMeters = 0.0;
    double shadowFinalResidualMaximumPascalsMeters = 0.0;
    double shadowMaximumCellConservationResidual = 0.0;
    bool hasTwoTerminalConductance = false;
    StableId lowerTerminalRegionId = invalidStableId;
    StableId upperTerminalRegionId = invalidStableId;
    double twoTerminalIntegratedTransferPascalsMeters = 0.0;
    double graphTwoTerminalConductanceMeters = 0.0;
    double shadowTwoTerminalConductanceMeters = 0.0;
    double graphToShadowTwoTerminalConductanceRatio = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureOperatorResponseModeDiagnostics&) const =
        default;
};

// Offline, immutable inverse-response comparison on one already accepted
// cut-cell topology. The optional accepted source is audited first; six
// deterministic, component-compatible manufactured pressure modes then span
// smooth coordinate, high-frequency, and authored-region directions. Each
// solve is independent, gauge aligned per component, and retains source work.
// A one-component/two-region contrast additionally reports the energy-based
// graph and shadow terminal conductances. The audit never changes either live
// pressure owner or Structure loads.
struct SceneFluidPressureOperatorResponseAudit {
    std::uint32_t version =
        sceneFluidPressureOperatorResponseAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t graphPressureOperatorFingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidPressureOperatorResponseAuditSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t controlVolumeCount = 0;
    std::size_t componentCount = 0;
    bool includesAcceptedSource = false;
    std::vector<SceneFluidPressureOperatorResponseModeDiagnostics> modes;
    std::vector<SceneFluidPressureOperatorResponseRecord> responses;

    bool operator==(
        const SceneFluidPressureOperatorResponseAudit&) const = default;
};

[[nodiscard]] SceneFluidPressureOperatorResponseAudit
auditSceneFluidPressureOperatorResponses(
    const SceneFluidPressureOperator& graphOperator,
    const SceneFluidMimeticControlCellSet& mimeticControlCells,
    const SceneFluidMimeticTraceSystem& fullTraceSystem,
    const SceneFluidMimeticCondensedTraceSystem& condensedTraceSystem,
    std::span<const double> acceptedIntegratedSourcePascalsMeters = {},
    const SceneFluidPressureOperatorResponseAuditSettings& settings = {},
    const SceneFluidPressureOperatorResponseAuditLimits& limits = {});

void validateSceneFluidPressureOperatorResponseAuditIntegrity(
    const SceneFluidPressureOperatorResponseAudit& audit);

} // namespace simwing::fsi
