#pragma once

#include "scene_fluid_region_momentum.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionTransportVersion = 1;

struct SceneFluidRegionTransportSettings {
    double timeStepSeconds = 1.0 / 60.0;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    double maximumOutgoingCourantNumber = 0.8;
    double maximumViscousNumber = 0.5;
    std::size_t maximumSubsteps = 1024;
    double absoluteMomentumToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-11;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-11;

    bool operator==(const SceneFluidRegionTransportSettings&) const = default;
};

struct SceneFluidRegionTransportLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumLinks = 100'000'000;
    std::size_t maximumTransportBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

enum class SceneFluidRegionTransportFailureStage : std::uint8_t {
    None = 0,
    FlowContinuity = 1,
    SubstepLimit = 2,
    AdvectionEnergy = 3,
    ViscosityEnergy = 4,
    Conservation = 5,
    NonFinite = 6,
};

struct SceneFluidRegionTransportControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    double volumeCubicMeters = 0.0;
    fluid::Vector3 velocityMetersPerSecond;
    fluid::Vector3 momentumKilogramMetersPerSecond;

    bool operator==(
        const SceneFluidRegionTransportControlVolume&) const = default;
};

struct SceneFluidRegionTransportDiagnostics {
    std::size_t controlVolumeCount = 0;
    std::size_t linkCount = 0;
    std::size_t openingLinkCount = 0;
    std::size_t substepCount = 0;
    double maximumCorrectedContinuityResidualCubicMetersPerSecond = 0.0;
    double maximumFullStepOutgoingCourantNumber = 0.0;
    double maximumAcceptedSubstepOutgoingCourantNumber = 0.0;
    double maximumFullStepViscousNumber = 0.0;
    double maximumAcceptedSubstepViscousNumber = 0.0;
    fluid::Vector3 momentumBeforeKilogramMetersPerSecond;
    fluid::Vector3 momentumAfterKilogramMetersPerSecond;
    fluid::Vector3 momentumResidualKilogramMetersPerSecond;
    double momentumResidualNormKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterAdvectionJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double advectiveEnergyLossJoules = 0.0;
    double viscousEnergyLossJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    SceneFluidRegionTransportFailureStage failureStage =
        SceneFluidRegionTransportFailureStage::None;
    bool finite = false;
    bool accepted = false;

    bool operator==(
        const SceneFluidRegionTransportDiagnostics&) const = default;
};

// Fixed-epoch conservative cell/region momentum advance. Corrected relative
// pressure-link flow carries donor-cell vector momentum between its two
// control-volume owners; graph viscosity then exchanges equal-and-opposite
// impulse over the same open fluid connection. Deterministic subcycling bounds
// both outgoing volume Courant number and the explicit pair diffusion number.
// Moving volumes, topology rebase, material-wall shear, and a pressure update
// remain outside this first transport boundary.
struct SceneFluidRegionTransport {
    std::uint32_t version = sceneFluidRegionTransportVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceMomentumFingerprint = 0;
    std::uint64_t pressureProjectionFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double sourceSimulationTimeSeconds = 0.0;
    double targetSimulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    SceneFluidRegionTransportSettings settings;
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionTransportDiagnostics diagnostics;
    std::vector<SceneFluidRegionTransportControlVolume> controlVolumes;

    bool operator==(const SceneFluidRegionTransport&) const = default;
};

[[nodiscard]] SceneFluidRegionTransport advanceSceneFluidRegionMomentum(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureProjection& correctedProjection,
    const SceneFluidRegionTransportSettings& settings = {},
    const SceneFluidRegionTransportLimits& limits = {});

void validateSceneFluidRegionTransportIntegrity(
    const SceneFluidRegionTransport& transport);

void validateSceneFluidRegionTransport(
    const SceneFluidRegionTransport& transport,
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureProjection& correctedProjection);

} // namespace simwing::fsi
