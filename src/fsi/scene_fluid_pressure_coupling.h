#pragma once

#include "coupling.h"
#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_pressure_sampling.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidPressureCouplingVersion = 1;
inline constexpr std::uint32_t
    sceneFluidPressureCouplingCheckpointVersion = 1;

struct SceneFluidPressureCouplingSettings {
    SceneFluidPressureCouplingSettings();

    StructureStepSettings structure;
    AitkenRelaxationSettings relaxation;
    CouplingConvergenceSettings convergence;
    SceneFluidPressureEpochSettings pressureEpoch;
    SceneFluidPressureProjectionSettings pressureProjection;
    ConservativeTransferSettings transfer;

};

struct SceneFluidPressureCouplingLimits {
    SceneFluidRegionConnectivityLimits connectivity;
    SceneFluidPressureEpochLimits pressureEpoch;
    SceneFluidPressureVolumeRateLimits volumeRates;
    SceneFluidOpeningFluxLimits openingFlux;
    SceneFluidPressureProjectionLimits pressureProjection;
    SceneFluidPressureSamplingLimits pressureSampling;
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
    ConservativeTransferDiagnostics pressureTransfer;
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
};

// First topology-stable strong feedback owner for the scene-v2 pressure path.
// The interface iterate is the canonical end-of-step pressure nodal load.
// Every solve restores the same Structure baseline, applies the trapezoidal
// average of the accepted start load and relaxed end-load guess, advances
// XPBD, rebuilds one atomic pressure epoch, projects moving-volume flow, and
// returns the resulting conservative pressure load to Aitken relaxation.
// Only a converged physical iterate becomes the next accepted baseline.
// Exhaustion, topology change, projection failure, and exceptions restore the
// caller's exact Structure checkpoint and leave this owner unchanged.
//
// The predicted MAC field is held fixed across the nonlinear iteration. This
// closes structural/pressure feedback but is not yet a momentum-evolving or
// topology-rebasing CFD macro-step.
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
    [[nodiscard]] const SceneFluidPressureProjection*
    acceptedPressureProjection() const noexcept;
    [[nodiscard]] const SceneFluidPressureSampleSet*
    acceptedPressureSamples() const noexcept;
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

private:
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
    std::optional<SceneFluidPressureProjection> acceptedPressureProjection_;
    std::optional<SceneFluidPressureSampleSet> acceptedPressureSamples_;
};

} // namespace simwing::fsi
