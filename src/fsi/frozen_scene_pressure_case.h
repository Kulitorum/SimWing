#pragma once

#include "fluid/grid.h"
#include "scene.h"
#include "structure.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

inline constexpr const char* frozenScenePressureSolverId =
    "frozen-scene-mimetic-pressure-v5";
inline constexpr const char* movingScenePressureSolverId =
    "moving-scene-mimetic-pressure-v1";
inline constexpr const char* movingSceneWallTracePressureSolverId =
    "moving-scene-wall-trace-pressure-v1";
inline constexpr const char* movingScenePreflowBootstrapSolverId =
    "moving-scene-preflow-bootstrap-v1";
inline constexpr const char* movingScenePreflowWallTraceSolverId =
    "moving-scene-preflow-wall-trace-v1";
inline constexpr const char* frozenSceneAggregateOpeningPressureSolverId =
    "frozen-scene-aggregate-opening-pressure-v1";
inline constexpr const char* movingSceneAggregateOpeningPressureSolverId =
    "moving-scene-aggregate-opening-pressure-v1";
inline constexpr const char*
    movingSceneAggregateOpeningWallTracePressureSolverId =
        "moving-scene-aggregate-opening-wall-trace-v1";
inline constexpr const char*
    movingSceneAggregateOpeningPreflowBootstrapSolverId =
        "moving-scene-aggregate-opening-preflow-v1";
inline constexpr const char*
    movingSceneAggregateOpeningPreflowWallTraceSolverId =
        "moving-scene-aggregate-opening-preflow-wall-trace-v1";

struct FrozenScenePressureCaseSettings {
    fluid::GridCellCounts cellCounts{2, 2, 2};
    double domainPaddingMeters = 0.5;
    bool useExplicitDomain = false;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    fluid::Vector3 backgroundWindMetersPerSecond{0.0, 0.85, 0.0};
    double windRampSeconds = 0.5;
    bool useCorrectedTraceFlowContinuation = false;
    bool useRegionalTransportFlowPrediction = false;
    bool useMovingGeometryFsi = false;
    bool useMaterialWallTracePressureLoads = false;
    bool aggregateCoplanarOpeningPatches = false;
    std::size_t preflowBootstrapSteps = 0;
    std::size_t movingStructureSubsteps = 8;
    double movingLoadRampSeconds = 0.0;
    double movingPressureRelativeResidualTolerance = 1.0e-5;
    double movingPressureReconstructedResidualTolerancePascalsMeters = 0.0;
    double diagnosticPerturbationSpeedMetersPerSecond = 0.0;
    double timeStepSeconds = 1.0 / 60.0;
    double densityKgPerCubicMeter = 1.225;
};

struct FrozenScenePressureCaseDiagnostics {
    fluid::GridCellCounts gridCellCounts;
    fluid::Vector3 gridLowerMeters;
    fluid::Vector3 gridUpperMeters;
    fluid::Vector3 backgroundWindMetersPerSecond;
    double windRampSeconds = 0.0;
    double windRampFraction = 0.0;
    bool usesCorrectedTraceFlowContinuation = false;
    bool usesRegionalTransportFlowPrediction = false;
    bool usesMovingGeometryFsi = false;
    bool usesMaterialWallTracePressureLoads = false;
    bool usesCoplanarOpeningAggregation = false;
    bool usesPreflowBootstrap = false;
    std::size_t preflowBootstrapStepCount = 0;
    std::size_t completedPreflowBootstrapStepCount = 0;
    bool usesConsecutivePressureWarmStart = false;
    bool usesRegionWallPrediction = false;
    bool pressureControlTopologyStable = false;
    std::size_t geometryAdvanceCount = 0;
    double movingLoadActivationFraction = 0.0;
    double maximumGeometryDisplacementMeters = 0.0;
    double maximumSuspensionResidualMeters = 0.0;
    double maximumCarriedTraceCorrectionCubicMetersPerSecond = 0.0;
    double maximumTraceBulkIncrementCubicMetersPerSecond = 0.0;
    double maximumRegionalTransportFlowDifferenceFromBulkBaselineCubicMetersPerSecond =
        0.0;
    double diagnosticPerturbationSpeedMetersPerSecond = 0.0;
    std::size_t pressureControlCount = 0;
    std::size_t sharedTraceCount = 0;
    std::size_t pressureIterationCount = 0;
    std::size_t extrapolatedZeroVolumePressureSideCount = 0;
    std::size_t reconstructedMaterialWallPressureSideCount = 0;
    std::size_t fallbackMaterialWallPressureSideCount = 0;
    double maximumPressureExtrapolationDistanceMeters = 0.0;
    double maximumAbsolutePressureDifferencePascals = 0.0;
    double maximumAbsoluteMaterialWallTracePressureDifferencePascals = 0.0;
    double maximumAbsoluteMaterialWallTraceDifferenceFromControlSamplePascals =
        0.0;
    double controlSampleIntegratedAbsolutePressureForceNewtons = 0.0;
    double controlSampleMaximumQuadraturePressureForceNewtons = 0.0;
    double materialWallTraceIntegratedAbsolutePressureForceNewtons = 0.0;
    double materialWallTraceMaximumQuadraturePressureForceNewtons = 0.0;
    double materialWallTraceMaximumNodalPressureForceNewtons = 0.0;
    double maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond =
        0.0;
    double correctedContinuityToleranceCubicMetersPerSecond = 0.0;
    double maximumCollapsedMacVelocityMetersPerSecond = 0.0;
    double maximumCollapsedSubfaceVelocityDeviationMetersPerSecond = 0.0;
    double maximumRegionalVelocityMetersPerSecond = 0.0;
    double maximumRegionalLinkVelocityResidualMetersPerSecond = 0.0;
    double regionalKineticEnergyJoules = 0.0;
    std::size_t regionalTransportSubstepCount = 0;
    double regionalTransportMaximumVelocityChangeMetersPerSecond = 0.0;
    double regionalTransportMomentumResidualKilogramMetersPerSecond = 0.0;
    double regionalTransportAdvectiveEnergyLossJoules = 0.0;
    double regionalTransportViscousEnergyLossJoules = 0.0;
    std::size_t embeddedOpeningTraceCount = 0;
    std::size_t flowAdvanceCount = 0;
    std::size_t bulkFlowSubstepCount = 0;
    double bulkFlowMaximumVelocityChangeMetersPerSecond = 0.0;
    fluid::Vector3 meanVelocityBeforePumpMetersPerSecond;
    fluid::Vector3 windPumpIncrementMetersPerSecond;
    double meanStreamwiseVelocityBeforePumpMetersPerSecond = 0.0;
    double streamwisePumpIncrementMetersPerSecond = 0.0;
    double bulkProjectionDivergenceBeforePerSecond = 0.0;
    double bulkProjectionDivergenceAfterPerSecond = 0.0;
    StructureVector3 pressureForceNewtons;
    StructureVector3 pressureMomentNewtonMeters;
    StructureVector3 materialWallTracePressureForceNewtons;
    double materialWallTracePressureTransferForceResidualNewtons = 0.0;
    double materialWallTracePressureTransferMomentResidualNewtonMeters = 0.0;
    StructureVector3 wallForceNewtons;
    StructureVector3 totalFluidForceNewtons;
    double wallMomentumResidualKilogramMetersPerSecond = 0.0;
    double totalFluidTransferForceResidualNewtons = 0.0;
    double totalFluidTransferMomentResidualNewtonMeters = 0.0;
    double transferForceResidualNewtons = 0.0;
    double transferMomentResidualNewtonMeters = 0.0;
    bool finite = false;

    bool operator==(
        const FrozenScenePressureCaseDiagnostics&) const = default;
};

// Loads one already authoritative scene-v2 payload and freezes its Structure
// geometry. The first frame evaluates a deterministic mixed-hybrid pressure
// projection from a prescribed periodic MAC predictor. Later frames advance a
// bounded target-wind ramp, project the collapsed continuation field, advance
// periodic bulk transport, and reapply the fixed-geometry mimetic pressure
// boundary. One opt-in continuation retains the exact last corrected trace
// flow and applies only the bulk predictor delta. A mutually exclusive opt-in
// instead projects the accepted conservative regional transport state onto the
// next pressure predictor. A separate opt-in moving mode applies the accepted
// pressure-plus-wall load, advances XPBD once, rebuilds the graph-free geometry,
// rebases regional momentum, and accepts the next consecutive pressure endpoint
// atomically. A further moving-only opt-in substitutes reconstructed material
// wall trace pressures for cell-constant pressure samples through the same
// conservative transfer while retaining both ledgers. Every path commits only
// after pressure, continuity, and
// conservative load transfer accept. This is an integration probe, not a
// validated external-flow wake, lift polar, or aerodynamic model.
class FrozenScenePressureCase final {
public:
    explicit FrozenScenePressureCase(
        Scene scene,
        const FrozenScenePressureCaseSettings& settings = {});
    ~FrozenScenePressureCase();

    FrozenScenePressureCase(const FrozenScenePressureCase&) = delete;
    FrozenScenePressureCase& operator=(
        const FrozenScenePressureCase&) = delete;
    FrozenScenePressureCase(FrozenScenePressureCase&&) noexcept;
    FrozenScenePressureCase& operator=(
        FrozenScenePressureCase&&) noexcept;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] const FrozenScenePressureCaseDiagnostics& diagnostics()
        const noexcept;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace simwing::fsi
