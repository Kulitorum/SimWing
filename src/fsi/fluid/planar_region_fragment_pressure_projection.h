#pragma once

#include "fluid/planar_region_fragment_pressure_solve.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

struct PlanarPressureRegionFragmentPressureProjectionSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;
    double absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeContinuityTolerance = 1.0e-10;
    PlanarPressureRegionFragmentPressureSolveSettings pressureSolve;

    bool operator==(
        const PlanarPressureRegionFragmentPressureProjectionSettings&) const =
        default;
};

struct PlanarPressureRegionFragmentPressureProjectionLimits {
    PlanarPressureRegionFragmentPressureOperatorLimits pressureOperatorLimits;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentPressureProjectionDiagnostics {
    bool accepted = false;
    bool finite = false;
    bool staticGeometry = false;
    std::uint64_t pressureOperatorFingerprint = 0;
    std::uint64_t topologyFingerprint = 0;
    std::uint64_t fragmentFingerprint = 0;
    std::size_t fragmentCount = 0;
    std::size_t linkCount = 0;
    std::size_t projectedSameRegionGridLinkCount = 0;
    std::size_t sealedPressureLayerWallLinkCount = 0;
    std::size_t workingStorageBytes = 0;
    double predictedNetOutwardFlowL2CubicMetersPerSecond = 0.0;
    double predictedNetOutwardFlowMaximumCubicMetersPerSecond = 0.0;
    double maximumAbsolutePredictedComponentBalanceCubicMetersPerSecond =
        0.0;
    double correctedNetOutwardFlowL2CubicMetersPerSecond = 0.0;
    double correctedNetOutwardFlowMaximumCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectedComponentBalanceCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteVelocityCorrectionMetersPerSecond = 0.0;
    double continuityToleranceCubicMetersPerSecond = 0.0;
    PlanarPressureRegionFragmentPressureSolveDiagnostics pressureSolve;

    bool operator==(
        const PlanarPressureRegionFragmentPressureProjectionDiagnostics&)
        const = default;
};

// Projects one oriented normal velocity per fragment-topology link on a
// static accepted geometry. A positive velocity points from the link's minus
// fragment to its plus fragment. Same-region grid links contribute their
// area-weighted flow to the integrated pressure RHS and receive the matching
// pressure-gradient correction. Pressure-layer walls must enter with exactly
// zero velocity and remain sealed; their separate authored static pressure
// jump is not part of the correction solve.
//
// This first velocity slice deliberately rejects moving layer geometry until
// the separate local fragment volume-rate product is explicitly composed into
// its continuity RHS. It owns no face momentum mass, kinetic-energy claim,
// opening conductance, or production worker state. Both link velocities and
// pressure correction are committed only after the solve and an explicit
// corrected-continuity check succeed.
[[nodiscard]] PlanarPressureRegionFragmentPressureProjectionDiagnostics
projectStaticPlanarPressureRegionFragmentFaceVelocities(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::vector<double>& orientedNormalVelocityMetersPerSecond,
    std::vector<double>& pressureCorrectionPascals,
    const PlanarPressureRegionFragmentPressureProjectionSettings& settings =
        {},
    const PlanarPressureRegionFragmentPressureProjectionLimits& limits = {});

} // namespace simwing::fsi::fluid
