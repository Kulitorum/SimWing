#pragma once

#include "scene_fluid_region_transport.h"
#include "scene_fluid_quadrature.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionWallExchangeVersion = 1;
inline constexpr std::uint32_t sceneFluidAcceptedWallTractionVersion = 1;

struct SceneFluidRegionWallSettings {
    double timeStepSeconds = 1.0 / 60.0;
    double kinematicViscositySquareMetersPerSecond = 1.5e-5;
    double minimumWallDistanceMeters = 1.0e-5;
    double maximumViscousNumber = 0.5;
    std::size_t maximumSubsteps = 1024;
    double absoluteMomentumToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-11;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-11;

    bool operator==(const SceneFluidRegionWallSettings&) const = default;
};

struct SceneFluidRegionWallLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumQuadraturePoints = 100'000'000;
    std::size_t maximumWallBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

enum class SceneFluidRegionWallFailureStage : std::uint8_t {
    None = 0,
    SubstepLimit = 1,
    Energy = 2,
    Conservation = 3,
    NonFinite = 4,
};

struct SceneFluidRegionWallControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    double volumeCubicMeters = 0.0;
    double incidentWallAreaSquareMeters = 0.0;
    fluid::Vector3 velocityMetersPerSecond;
    fluid::Vector3 momentumKilogramMetersPerSecond;

    bool operator==(
        const SceneFluidRegionWallControlVolume&) const = default;
};

struct SceneFluidRegionWallSample {
    std::size_t sampleIndex = 0;
    std::uint64_t stableId = 0;
    StableId triangleId = invalidStableId;
    std::size_t negativeSideControlVolumeIndex = 0;
    std::size_t positiveSideControlVolumeIndex = 0;
    double areaSquareMeters = 0.0;
    fluid::Vector3 unitNormalNegativeToPositive;
    fluid::Vector3 wallVelocityMetersPerSecond;
    fluid::Vector3 negativeSideFluidImpulseKilogramMetersPerSecond;
    fluid::Vector3 positiveSideFluidImpulseKilogramMetersPerSecond;
    SceneFluidQuadratureTraction structureTraction;

    bool operator==(const SceneFluidRegionWallSample&) const = default;
};

struct SceneFluidRegionWallDiagnostics {
    std::size_t controlVolumeCount = 0;
    std::size_t quadraturePointCount = 0;
    std::size_t substepCount = 0;
    double maximumFullStepViscousNumber = 0.0;
    double maximumAcceptedSubstepViscousNumber = 0.0;
    double maximumWallDistanceMeters = 0.0;
    double maximumRelativeTangentialSpeedMetersPerSecond = 0.0;
    fluid::Vector3 fluidMomentumBeforeKilogramMetersPerSecond;
    fluid::Vector3 fluidMomentumAfterKilogramMetersPerSecond;
    fluid::Vector3 fluidImpulseKilogramMetersPerSecond;
    fluid::Vector3 structureImpulseKilogramMetersPerSecond;
    fluid::Vector3 momentumResidualKilogramMetersPerSecond;
    double momentumResidualNormKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double wallWorkOnFluidJoules = 0.0;
    double viscousDissipationJoules = 0.0;
    SceneFluidRegionWallFailureStage failureStage =
        SceneFluidRegionWallFailureStage::None;
    bool finite = false;
    bool accepted = false;

    bool operator==(const SceneFluidRegionWallDiagnostics&) const = default;
};

// Two-sided tangential exchange between a transported cell/region state and
// the current material surface. Each quadrature side uses its exact authored
// cell/region owner. Explicit subcycling bounds the aggregate viscous row per
// control volume. Fluid impulse and the equal-and-opposite structure traction
// are published together, with wall work separated from nonnegative viscous
// dissipation. Normal exchange remains owned by pressure projection.
struct SceneFluidRegionWallExchange {
    std::uint32_t version = sceneFluidRegionWallExchangeVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t currentPressureControlVolumeFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    SceneFluidRegionWallSettings settings;
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionWallDiagnostics diagnostics;
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes;
    std::vector<SceneFluidRegionWallSample> samples;

    bool operator==(const SceneFluidRegionWallExchange&) const = default;
};

// Minimal accepted endpoint needed by strong-coupling restart. The full wall
// exchange owns transient adjusted fluid controls; this state retains only the
// conservative Structure traction with exact wall/quadrature provenance.
struct SceneFluidAcceptedWallTractionSet {
    std::uint32_t version = sceneFluidAcceptedWallTractionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t wallExchangeFingerprint = 0;
    std::uint64_t quadratureFingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::vector<SceneFluidQuadratureTraction> tractions;

    bool operator==(const SceneFluidAcceptedWallTractionSet&) const = default;
};

[[nodiscard]] SceneFluidRegionWallExchange
exchangeSceneFluidRegionWallMomentum(
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionWallSettings& settings = {},
    const SceneFluidRegionWallLimits& limits = {});

void validateSceneFluidRegionWallExchangeIntegrity(
    const SceneFluidRegionWallExchange& exchange);

void validateSceneFluidRegionWallExchange(
    const SceneFluidRegionWallExchange& exchange,
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature);

[[nodiscard]] SceneFluidAcceptedWallTractionSet
captureSceneFluidAcceptedWallTractions(
    const SceneFluidRegionWallExchange& exchange);

void validateSceneFluidAcceptedWallTractionSetIntegrity(
    const SceneFluidAcceptedWallTractionSet& tractions);

void validateSceneFluidAcceptedWallTractions(
    const SceneFluidAcceptedWallTractionSet& tractions,
    const SceneFluidQuadratureDefinition& quadrature,
    std::uint64_t expectedWallExchangeFingerprint);

} // namespace simwing::fsi
