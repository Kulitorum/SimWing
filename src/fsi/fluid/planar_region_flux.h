#pragma once

#include "fluid/planar_region_sweep.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarPressureRegionFluxVersion = 1;

struct PlanarPressureRegionFluxSettings {
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-12;
    double relativeVelocityTolerance = 1.0e-12;
    double absoluteVolumeToleranceCubicMeters = 1.0e-12;
    double relativeVolumeTolerance = 1.0e-12;

    bool operator==(const PlanarPressureRegionFluxSettings&) const = default;
};

struct PlanarPressureRegionFluxLimits {
    std::size_t maximumIntervals = 1'000'000;
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumOwnedBytes = 256ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionIntervalFlux {
    std::uint64_t lowerSurfaceStableId = 0;
    std::uint64_t upperSurfaceStableId = 0;
    std::uint64_t regionStableId = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double lowerSurfaceVelocityMetersPerSecond = 0.0;
    double upperSurfaceVelocityMetersPerSecond = 0.0;
    double leastSquaresFluidVelocityMetersPerSecond = 0.0;
    double lowerOutwardRelativeVelocityMetersPerSecond = 0.0;
    double upperOutwardRelativeVelocityMetersPerSecond = 0.0;
    double lowerOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double upperOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double outwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double integratedOutwardRelativeVolumeCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double maximumAbsoluteInterfaceSlipMetersPerSecond = 0.0;
    double rmsInterfaceSlipMetersPerSecond = 0.0;
    double totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond = 0.0;
    double velocityToleranceMetersPerSecond = 0.0;
    double continuityToleranceCubicMeters = 0.0;
    bool impermeableWithinTolerance = false;
    bool continuityWithinTolerance = false;

    bool operator==(const PlanarPressureRegionIntervalFlux&) const = default;
};

struct PlanarPressureRegionFluxSummary {
    std::uint64_t regionStableId = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double outwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double integratedOutwardRelativeVolumeCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double maximumAbsoluteInterfaceSlipMetersPerSecond = 0.0;
    double totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond = 0.0;
    double continuityToleranceCubicMeters = 0.0;
    bool impermeableWithinTolerance = false;
    bool continuityWithinTolerance = false;

    bool operator==(const PlanarPressureRegionFluxSummary&) const = default;
};

// Minimum one-dimensional material-relative flux implied by one uniform axial
// fluid velocity per planar interval. The least-squares velocity is the mean
// of the two boundary velocities. Rigid motion is exactly impermeable; unequal
// boundary velocities expose unavoidable one-sided slip. Continuity still
// closes because that relative flux is integrated with the regional volume
// change. This is an offline compatibility screen, not a fluid solve.
struct PlanarPressureRegionFluxCompatibility {
    std::uint32_t version = planarPressureRegionFluxVersion;
    std::uint32_t sourceSweepVersion = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    double durationSeconds = 0.0;
    double crossSectionAreaSquareMeters = 0.0;
    PlanarPressureRegionFluxSettings settings;
    std::vector<PlanarPressureRegionIntervalFlux> intervals;
    std::vector<PlanarPressureRegionFluxSummary> regions;
    std::size_t failedImpermeableIntervalCount = 0;
    std::size_t failedContinuityIntervalCount = 0;
    std::size_t failedImpermeableRegionCount = 0;
    std::size_t failedContinuityRegionCount = 0;
    double maximumAbsoluteInterfaceSlipMetersPerSecond = 0.0;
    double maximumAbsoluteContinuityResidualCubicMeters = 0.0;
    double globalGeometryVolumeChangeCubicMeters = 0.0;
    double globalIntegratedOutwardRelativeVolumeCubicMeters = 0.0;
    double globalContinuityResidualCubicMeters = 0.0;
    bool allIntervalsImpermeableWithinTolerance = false;
    bool allIntervalsContinuousWithinTolerance = false;
    bool allRegionsImpermeableWithinTolerance = false;
    bool allRegionsContinuousWithinTolerance = false;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFluxCompatibility&) const = default;
};

[[nodiscard]] PlanarPressureRegionFluxCompatibility
assessPlanarPressureRegionFluxCompatibility(
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFluxSettings& settings = {},
    const PlanarPressureRegionFluxLimits& limits = {});

} // namespace simwing::fsi::fluid
