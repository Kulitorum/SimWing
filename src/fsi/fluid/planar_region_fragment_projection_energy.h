#pragma once

#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_velocity_state.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentProjectionEnergyVersion = 1;

struct PlanarPressureRegionFragmentProjectionEnergySettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    double absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeContinuityTolerance = 1.0e-10;
    double absoluteVelocityResidualToleranceMetersPerSecond = 1.0e-12;
    double relativeVelocityResidualTolerance = 1.0e-10;
    double absoluteMomentumResidualToleranceKilogramMetersPerSecond =
        1.0e-12;
    double relativeMomentumResidualTolerance = 1.0e-10;
    double absoluteEnergyResidualToleranceJoules = 1.0e-12;
    double relativeEnergyResidualTolerance = 1.0e-10;
    double absolutePressureGaugeTolerancePascals = 1.0e-12;
    double relativePressureGaugeTolerance = 1.0e-10;

    bool operator==(
        const PlanarPressureRegionFragmentProjectionEnergySettings&) const =
        default;
};

struct PlanarPressureRegionFragmentProjectionEnergyCorrection {
    std::size_t correctionIndex = 0;
    std::size_t dofIndex = 0;
    std::uint64_t dofStableId = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t componentIndex = 0;
    std::uint64_t regionStableId = 0;
    double areaSquareMeters = 0.0;
    double centerDistanceMeters = 0.0;
    double dualVolumeCubicMeters = 0.0;
    double diagonalMassKilograms = 0.0;
    double pressureDifferencePascals = 0.0;
    double velocityBeforeMetersPerSecond = 0.0;
    double velocityAfterMetersPerSecond = 0.0;
    double expectedVelocityChangeMetersPerSecond = 0.0;
    double velocityChangeResidualMetersPerSecond = 0.0;
    double momentumChangeKilogramMetersPerSecond = 0.0;
    double pressureImpulseKilogramMetersPerSecond = 0.0;
    double momentumImpulseResidualKilogramMetersPerSecond = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double midpointPressureWorkJoules = 0.0;
    double workEnergyResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentProjectionEnergyCorrection&) const =
        default;
};

struct PlanarPressureRegionFragmentProjectionEnergyComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t correctionCount = 0;
    double pressureCorrectionVolumeMeanPascals = 0.0;
    double predictedNetOutwardFlowCubicMetersPerSecond = 0.0;
    double correctedNetOutwardFlowCubicMetersPerSecond = 0.0;
    Vector3 momentumChangeKilogramMetersPerSecond;
    Vector3 pressureImpulseKilogramMetersPerSecond;
    Vector3 momentumImpulseResidualKilogramMetersPerSecond;
    double kineticEnergyChangeJoules = 0.0;
    double midpointPressureWorkJoules = 0.0;
    double workEnergyResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentProjectionEnergyComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentProjectionEnergyLimits {
    PlanarPressureRegionFragmentVelocityStateLimits velocityStateLimits;
    std::size_t maximumCorrections = 120'000'000;
    std::size_t maximumPressureSamples = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Independent certificate for one accepted static regional pressure
// correction. The before/after states use the same diagonal face metric.
// Same-region velocity change must equal dt/rho times the correction-pressure
// difference over center distance. The matching pressure impulse must equal
// momentum change, and midpoint pressure work must equal kinetic-energy
// change. Both one-sided pressure-layer traces remain exactly zero and the
// corrected static continuity residual closes per fragment.
//
// This audit covers the correction potential only, not the separate authored
// static pressure jump. It does not prescribe moving-wall traces, certify an
// affine moving-volume projection, transport momentum, or enter production.
struct PlanarPressureRegionFragmentProjectionEnergyAudit {
    std::uint32_t version =
        planarPressureRegionFragmentProjectionEnergyVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceMetricFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t beforeVelocityStateFingerprint = 0;
    std::uint64_t afterVelocityStateFingerprint = 0;
    PlanarPressureRegionFragmentProjectionEnergySettings settings;
    std::vector<double> pressureCorrectionPascals;
    std::vector<PlanarPressureRegionFragmentProjectionEnergyCorrection>
        corrections;
    std::vector<PlanarPressureRegionFragmentProjectionEnergyComponent>
        components;
    std::size_t sealedPressureLayerTraceCount = 0;
    double predictedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double predictedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double correctedContinuityResidualL2CubicMetersPerSecond = 0.0;
    double correctedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    double continuityToleranceCubicMetersPerSecond = 0.0;
    double maximumAbsoluteVelocityChangeResidualMetersPerSecond = 0.0;
    double maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        0.0;
    double maximumAbsolutePressureGaugePascals = 0.0;
    Vector3 momentumBeforeKilogramMetersPerSecond;
    Vector3 momentumAfterKilogramMetersPerSecond;
    Vector3 momentumChangeKilogramMetersPerSecond;
    Vector3 pressureImpulseKilogramMetersPerSecond;
    Vector3 momentumImpulseResidualKilogramMetersPerSecond;
    double kineticEnergyBeforeJoules = 0.0;
    double kineticEnergyAfterJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double midpointPressureWorkJoules = 0.0;
    double workEnergyResidualJoules = 0.0;
    double kineticEnergyRemovedJoules = 0.0;
    bool nonIncreasingKineticEnergy = false;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentProjectionEnergyAudit&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentProjectionEnergyAudit
auditStaticPlanarPressureRegionFragmentProjectionEnergy(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentProjectionEnergySettings& settings = {},
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits = {});

void validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
    const PlanarPressureRegionFragmentProjectionEnergyAudit& audit,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityState& before,
    const PlanarPressureRegionFragmentVelocityState& after,
    const PlanarPressureRegionFragmentProjectionEnergyLimits& limits = {});

} // namespace simwing::fsi::fluid
