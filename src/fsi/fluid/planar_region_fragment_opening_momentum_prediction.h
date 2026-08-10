#pragma once

#include "fluid/planar_region_fragment_opening_momentum_transport.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumPredictionVersion = 1;

struct PlanarPressureRegionFragmentOpeningMomentumPredictionLimits {
    PlanarPressureRegionFragmentOpeningVelocityStateLimits stateLimits;
    PlanarPressureRegionFragmentVolumeRateLimits volumeRateLimits;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumDofs = 140'000'000;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 8192ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningMomentumPredictionDiagnostics {
    std::size_t fragmentCount = 0;
    std::size_t dofCount = 0;
    std::size_t sharedRegionGridDofCount = 0;
    std::size_t solidWallTraceDofCount = 0;
    std::size_t openingDofCount = 0;
    Vector3 sourceMomentumKilogramMetersPerSecond;
    Vector3 remappedMomentumKilogramMetersPerSecond;
    Vector3 geometricMomentumChangeKilogramMetersPerSecond;
    Vector3 predictedStateMomentumKilogramMetersPerSecond;
    Vector3 reconstructionMomentumChangeKilogramMetersPerSecond;
    double sourceKineticEnergyJoules = 0.0;
    double remappedKineticEnergyJoules = 0.0;
    double predictedStateKineticEnergyJoules = 0.0;
    double maximumAbsoluteVolumeChangeCubicMeters = 0.0;
    double maximumSourceSpeedMetersPerSecond = 0.0;
    double maximumEndpointNormalVelocityJumpMetersPerSecond = 0.0;
    double maximumAbsolutePredictedRelativeVelocityMetersPerSecond = 0.0;
    double maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond =
        0.0;
    bool finite = false;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumPredictionDiagnostics&)
        const = default;
};

// First-order face/aperture predictor for one consecutive topology-stable
// pressure epoch. Accepted transported fragment velocity is retained while
// density times the new physical volume diagnoses the geometric momentum
// remap. Each current same-region or aperture degree receives the arithmetic
// mean of its two endpoint vector components. Fixed grid faces therefore keep
// that absolute value as relative flow; retained solid traces take exact
// material motion; apertures subtract exact material motion from the averaged
// absolute fluid velocity.
//
// The nested immutable velocity state is ready to seed a later pressure
// transaction. This product does not perform pressure projection, viscosity,
// wall shear, topology rebase, or production-worker mutation.
struct PlanarPressureRegionFragmentOpeningMomentumPrediction {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumPredictionVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceTransportFingerprint = 0;
    std::uint64_t sourceTransportMetricFingerprint = 0;
    std::uint64_t currentMetricFingerprint = 0;
    std::uint64_t currentVolumeRateFingerprint = 0;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    PlanarPressureRegionFragmentOpeningMomentumPredictionDiagnostics
        diagnostics;
    PlanarPressureRegionFragmentOpeningVelocityState predictedVelocityState;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumPrediction&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningMomentumPrediction
predictPlanarPressureRegionFragmentOpeningMomentum(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction);

void validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
    const PlanarPressureRegionFragmentOpeningMomentumPrediction& prediction,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportTargetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentVelocityMetric& currentBaseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& currentMetric,
    const PlanarPressureRegionFragmentOpeningMomentumPredictionLimits& limits =
        {});

} // namespace simwing::fsi::fluid
