#pragma once

#include "coupling.h"
#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_pressure_sampling.h"
#include "scene_fluid_region_link_flow.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureCouplingVersion = 1;
inline constexpr std::uint32_t
    sceneFluidPressureCouplingCheckpointVersion = 2;
inline constexpr std::uint32_t
    sceneFluidPressureMacVelocityCollapseVersion = 1;

struct SceneFluidPressureCouplingSettings {
    SceneFluidPressureCouplingSettings();

    StructureStepSettings structure;
    AitkenRelaxationSettings relaxation;
    CouplingConvergenceSettings convergence;
    SceneFluidPressureEpochSettings pressureEpoch;
    SceneFluidPressureProjectionSettings pressureProjection;
    SceneFluidRegionWallSettings regionWall;
    ConservativeTransferSettings transfer;

};

struct SceneFluidPressureCouplingLimits {
    SceneFluidRegionConnectivityLimits connectivity;
    SceneFluidPressureEpochLimits pressureEpoch;
    SceneFluidPressureVolumeRateLimits volumeRates;
    SceneFluidOpeningFluxLimits openingFlux;
    SceneFluidPressureProjectionLimits pressureProjection;
    SceneFluidPressureSamplingLimits pressureSampling;
    SceneFluidRegionLinkFlowLimits regionLinkFlow;
    SceneFluidRegionRebaseLimits regionRebase;
    SceneFluidRegionWallLimits regionWall;
    std::size_t maximumCouplingNodes = 10'000'000;
    std::size_t maximumInterfaceBytes =
        1024ULL * 1024ULL * 1024ULL;

};

struct SceneFluidPressureCouplingStepDiagnostics {
    std::uint32_t version = sceneFluidPressureCouplingVersion;
    bool accepted = false;
    std::uint64_t solverRunCount = 0;
    std::uint64_t previousPressureEpochFingerprint = 0;
    std::uint64_t currentPressureEpochFingerprint = 0;
    StrongCouplingIterationResult iteration;
    StructureDiagnostics structure;
    SceneFluidPressureProjectionDiagnostics pressureProjection;
    SceneFluidRegionWallDiagnostics regionWall;
    bool usesRegionWall = false;
    SceneFluidRegionRebaseDiagnostics regionRebase;
    bool usesRegionRebase = false;
    ConservativeTransferDiagnostics pressureTransfer;
    ConservativeTransferDiagnostics totalFluidTransfer;
    double interfaceForceClosureNewtons = 0.0;
    double interfaceForceReferenceNewtons = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureCouplingStepDiagnostics&) const = default;
};

struct SceneFluidPressureCouplingCheckpoint {
    std::uint32_t version = sceneFluidPressureCouplingCheckpointVersion;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidPressureCouplingSettings settings;
    StructureCheckpoint structure;
    std::optional<SceneFluidPressureProjection> pressureProjection;
    std::optional<SceneFluidAcceptedWallTractionSet> wallTractions;
};

struct SceneFluidPressureMacVelocityCollapseDiagnostics {
    std::uint32_t version =
        sceneFluidPressureMacVelocityCollapseVersion;
    std::uint64_t pressureProjectionFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t faceCount = 0;
    std::size_t linkCount = 0;
    std::size_t openingLinkCount = 0;
    std::size_t multiLinkFaceCount = 0;
    double maximumAbsoluteVelocityMetersPerSecond = 0.0;
    double maximumSubfaceVelocityDeviationMetersPerSecond = 0.0;
    double maximumVolumeFlowClosureCubicMetersPerSecond = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidPressureMacVelocityCollapseDiagnostics&) const =
        default;
};

struct SceneFluidPressureMacVelocityCollapse {
    explicit SceneFluidPressureMacVelocityCollapse(
        const fluid::PeriodicCartesianGrid& grid)
        : velocityMetersPerSecond(grid) {}

    fluid::MacVelocityField velocityMetersPerSecond;
    SceneFluidPressureMacVelocityCollapseDiagnostics diagnostics;

    bool operator==(
        const SceneFluidPressureMacVelocityCollapse&) const = default;
};

// First strong feedback owner for the scene-v2 pressure path.
// The interface iterate is the canonical end-of-step total fluid nodal load.
// Every solve restores the same Structure baseline, applies the trapezoidal
// average of the accepted start load and relaxed end-load guess, advances
// XPBD, rebuilds one atomic pressure epoch, projects moving-volume flow, and
// returns the resulting conservative pressure-plus-wall load to Aitken
// relaxation.
// Only a converged physical iterate becomes the next accepted baseline.
// Exhaustion, topology change, projection failure, and exceptions restore the
// caller's exact Structure checkpoint and leave this owner unchanged.
//
// The predicted MAC field is held fixed across the nonlinear iteration. The
// region-transport overload likewise holds one accepted transport fixed,
// remaps it to every current geometry iterate, exchanges tangential momentum
// with the material wall, and projects the adjusted link flow. Its bounded
// first crossing subset admits appeared controls using current same-region
// one-ring donors, and retires a disappeared control only to one unique
// previous same-region neighbour. This is not a general swept-volume topology
// remap.
class SceneFluidPressureCoupling final {
public:
    SceneFluidPressureCoupling(
        SceneFluidSurfaceDefinition surface,
        SceneStructureMappings structureMappings,
        const Structure& target,
        fluid::PeriodicCartesianGrid grid,
        const SceneFluidPressureCouplingSettings& settings = {},
        const SceneFluidPressureCouplingLimits& limits = {});

    SceneFluidPressureCoupling(
        const SceneFluidPressureCoupling&) = delete;
    SceneFluidPressureCoupling& operator=(
        const SceneFluidPressureCoupling&) = delete;
    SceneFluidPressureCoupling(
        SceneFluidPressureCoupling&&) noexcept = default;
    SceneFluidPressureCoupling& operator=(
        SceneFluidPressureCoupling&&) noexcept = default;

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const SceneFluidSurfaceState& acceptedSurfaceState()
        const noexcept;
    [[nodiscard]] const SceneFluidPressureEpoch& acceptedPressureEpoch()
        const noexcept;
    [[nodiscard]] const ConservativeTransferResult& acceptedPressureTransfer()
        const noexcept;
    [[nodiscard]] const ConservativeTransferResult&
    acceptedPressureOnlyTransfer() const noexcept;
    [[nodiscard]] const SceneFluidPressureProjection*
    acceptedPressureProjection() const noexcept;
    [[nodiscard]] const SceneFluidPressureSampleSet*
    acceptedPressureSamples() const noexcept;
    [[nodiscard]] const SceneFluidAcceptedWallTractionSet*
    acceptedWallTractions() const noexcept;
    // Collapses accepted link-resolved corrected flow back onto one bulk MAC
    // normal velocity per Cartesian face. Opening links first restore cap
    // sweep so the MAC value is absolute fluid velocity. Mixed-region
    // subface differences are area-averaged and reported, not hidden. This is
    // a deterministic continuation state, not cut-cell advection/diffusion.
    [[nodiscard]] SceneFluidPressureMacVelocityCollapse
    acceptedPressureCorrectedMacVelocity() const;
    [[nodiscard]] const SceneFluidPressureCouplingSettings& settings()
        const noexcept;
    [[nodiscard]] SceneFluidPressureCouplingCheckpoint checkpoint(
        const Structure& target) const;
    void restore(
        Structure& target,
        const SceneFluidPressureCouplingCheckpoint& checkpoint);

    [[nodiscard]] SceneFluidPressureCouplingStepDiagnostics advance(
        Structure& target,
        const fluid::MacVelocityField& predictedVelocityMetersPerSecond);

    // Uses one already-accepted fixed-epoch region transport as the immutable
    // fluid predictor across every structural iterate of this macro step.
    [[nodiscard]] SceneFluidPressureCouplingStepDiagnostics advance(
        Structure& target,
        const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
        const SceneFluidRegionTransport& transportedRegionMomentum);

private:
    [[nodiscard]] SceneFluidPressureCouplingStepDiagnostics advanceImpl(
        Structure& target,
        const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
        const SceneFluidRegionTransport* transportedRegionMomentum);

    SceneFluidSurfaceDefinition surface_;
    SceneStructureMappings structureMappings_;
    fluid::PeriodicCartesianGrid grid_;
    SceneFluidRegionConnectivity connectivity_;
    SceneFluidSurfaceTransfer transfer_;
    ConservativeMacroStepCoupling macroCoupling_;
    SceneFluidPressureCouplingSettings settings_;
    SceneFluidPressureCouplingLimits limits_;
    SceneFluidSurfaceState acceptedSurfaceState_;
    SceneFluidPressureEpoch acceptedPressureEpoch_;
    std::vector<double> acceptedPressurePascals_;
    std::optional<ConservativeTransferResult> acceptedPressureTransfer_;
    std::optional<ConservativeTransferResult> acceptedPressureOnlyTransfer_;
    std::optional<SceneFluidPressureProjection> acceptedPressureProjection_;
    std::optional<SceneFluidPressureSampleSet> acceptedPressureSamples_;
    std::optional<SceneFluidAcceptedWallTractionSet> acceptedWallTractions_;
};

} // namespace simwing::fsi
