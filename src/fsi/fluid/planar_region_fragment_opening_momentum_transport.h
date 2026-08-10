#pragma once

#include "fluid/planar_region_fragment_opening_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningMomentumTransportVersion = 1;

struct PlanarPressureRegionFragmentOpeningMomentumTransportSettings {
    double maximumOutgoingCourantNumber = 0.8;
    std::size_t maximumSubsteps = 1024;
    double absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeContinuityTolerance = 1.0e-10;
    double absoluteMomentumToleranceKilogramMetersPerSecond = 1.0e-12;
    double relativeMomentumTolerance = 1.0e-11;
    double absoluteEnergyToleranceJoules = 1.0e-12;
    double relativeEnergyTolerance = 1.0e-11;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningMomentumTransportLimits {
    PlanarPressureRegionFragmentOpeningVelocityStateLimits stateLimits;
    PlanarPressureRegionFragmentVolumeRateLimits volumeRateLimits;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumDofs = 140'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 8192ULL * 1024ULL * 1024ULL;
};

enum class PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage
    : std::uint8_t {
    None = 0,
    FlowContinuity = 1,
    SubstepLimit = 2,
    AdvectionEnergy = 3,
    Conservation = 4,
    GeometryVolume = 5,
    NonFinite = 6,
};

struct PlanarPressureRegionFragmentOpeningMomentumTransportControl {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t connectedComponentIndex = 0;
    double volumeCubicMeters = 0.0;
    Vector3 velocityMetersPerSecond;
    Vector3 momentumKilogramMetersPerSecond;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumTransportControl&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics {
    std::size_t fragmentCount = 0;
    std::size_t transportDofCount = 0;
    std::size_t openingDofCount = 0;
    std::size_t substepCount = 0;
    double maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond =
        0.0;
    double maximumContinuityResidualCubicMetersPerSecond = 0.0;
    double continuityToleranceCubicMetersPerSecond = 0.0;
    double maximumFullStepOutgoingCourantNumber = 0.0;
    double maximumAcceptedSubstepOutgoingCourantNumber = 0.0;
    Vector3 momentumBeforeKilogramMetersPerSecond;
    Vector3 momentumAfterKilogramMetersPerSecond;
    Vector3 momentumResidualKilogramMetersPerSecond;
    double momentumResidualNormKilogramMetersPerSecond = 0.0;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double advectiveKineticEnergyLossJoules = 0.0;
    double maximumVelocityChangeMetersPerSecond = 0.0;
    PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage
        failureStage =
            PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                None;
    bool finite = false;
    bool accepted = false;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics&)
        const = default;
};

// First-order conservative ALE advance of complete collocated fragment
// momentum. The source state supplies one vector momentum per previous
// physical fragment. The target velocity state supplies only the corrected
// relative normal flow on fixed-grid and aperture degrees; retained solid
// traces carry no inter-fragment mass. Donor selection transports the complete
// vector through each connection while target fragment volumes advance
// linearly through the accepted geometry rate.
//
// Corrected continuity makes a uniform vector field an exact discrete-GCL
// solution. Deterministic subcycling bounds outgoing volume Courant number;
// pair impulses conserve global momentum and donor-cell mixing may only reduce
// collocated kinetic energy. This opt-in artifact does not apply pressure,
// viscosity, wall shear, topology rebase, or production-worker state.
struct PlanarPressureRegionFragmentOpeningMomentumTransport {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningMomentumTransportVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceStateFingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t targetFlowStateFingerprint = 0;
    std::uint64_t targetMetricFingerprint = 0;
    std::uint64_t targetVolumeRateFingerprint = 0;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    PlanarPressureRegionFragmentOpeningMomentumTransportSettings settings;
    PlanarPressureRegionFragmentOpeningMomentumTransportDiagnostics
        diagnostics;
    std::vector<
        PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        controls;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningMomentumTransport&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningMomentumTransport
advancePlanarPressureRegionFragmentOpeningMomentum(
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityState& targetFlowState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& targetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& targetSweep,
    const PlanarPressureRegionFragmentSet& targetFragments,
    const PlanarPressureRegionFragmentTopology& targetTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& targetVolumeRates,
    const PlanarPressureRegionFragmentOpeningMomentumTransportSettings&
        settings = {},
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport);

void validatePlanarPressureRegionFragmentOpeningMomentumTransport(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const PlanarPressureRegionFragmentOpeningVelocityState& sourceState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& sourceMetric,
    const PlanarPressureRegionFragmentOpeningVelocityState& targetFlowState,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& targetMetric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& targetSweep,
    const PlanarPressureRegionFragmentSet& targetFragments,
    const PlanarPressureRegionFragmentTopology& targetTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& targetVolumeRates,
    const PlanarPressureRegionFragmentOpeningMomentumTransportLimits& limits =
        {});

} // namespace simwing::fsi::fluid
