#pragma once

#include "fluid/planar_region_sweep.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarPressureRegionOpeningFlowVersion = 1;

struct PlanarPressureRegionOpeningDefinition {
    std::uint64_t openingStableId = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    double areaSquareMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionOpeningDefinition&) const = default;
};

struct PlanarPressureRegionOpeningFlowSettings {
    double absoluteFlowRateToleranceCubicMetersPerSecond = 1.0e-12;
    double relativeFlowRateTolerance = 1.0e-12;
    double absoluteVolumeToleranceCubicMeters = 1.0e-12;
    double relativeVolumeTolerance = 1.0e-12;
    double relativeCholeskyPivotTolerance = 1.0e-13;

    bool operator==(
        const PlanarPressureRegionOpeningFlowSettings&) const = default;
};

struct PlanarPressureRegionOpeningFlowLimits {
    std::size_t maximumIntervals = 16'384;
    std::size_t maximumRegions = 4096;
    std::size_t maximumOpenings = 16'384;
    std::size_t maximumOwnedBytes = 256ULL * 1024ULL * 1024ULL;
    std::size_t maximumFactorizationBytes = 256ULL * 1024ULL * 1024ULL;
    std::uint64_t maximumFactorizationWork = 500'000'000ULL;
};

struct PlanarPressureRegionOpeningFlow {
    std::uint64_t openingStableId = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    std::uint64_t componentStableId = 0;
    double areaSquareMeters = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;

    bool operator==(const PlanarPressureRegionOpeningFlow&) const = default;
};

struct PlanarPressureRegionOpeningBalance {
    std::uint64_t regionStableId = 0;
    std::uint64_t componentStableId = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double requestedOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double solvedOutwardRelativeFlowRateCubicMetersPerSecond = 0.0;
    double integratedOutwardRelativeVolumeCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double continuityToleranceCubicMeters = 0.0;
    bool withinTolerance = false;

    bool operator==(
        const PlanarPressureRegionOpeningBalance&) const = default;
};

struct PlanarPressureRegionOpeningComponent {
    std::uint64_t componentStableId = 0;
    std::uint64_t correctionRegionStableId = 0;
    std::size_t regionCount = 0;
    std::size_t openingCount = 0;
    double requestedOutwardFlowRateSumCubicMetersPerSecond = 0.0;
    double appliedCorrectionFlowRateCubicMetersPerSecond = 0.0;
    double compatibilityToleranceCubicMetersPerSecond = 0.0;
    double maximumAbsoluteFlowBalanceResidualCubicMetersPerSecond = 0.0;
    bool sourceCompatible = false;
    bool linearSolveWithinTolerance = false;
    bool feasible = false;

    bool operator==(
        const PlanarPressureRegionOpeningComponent&) const = default;
};

// Positive opening flow travels from the authored negative region to the
// positive region. The oracle minimizes sum(flow^2 / area), so parallel
// openings carry equal normal velocity. Per-region continuity is
//
//   volume change + duration * outward relative opening flow = 0.
//
// Every connected opening component must balance independently. A sealed
// moving component is reported infeasible rather than repaired through fabric.
// This is a bounded offline topology/continuity oracle, not a fluid state or
// pressure/velocity solve.
struct PlanarPressureRegionOpeningFlowAllocation {
    std::uint32_t version = planarPressureRegionOpeningFlowVersion;
    std::uint64_t fingerprint = 0;
    std::uint32_t sourceSweepVersion = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    double durationSeconds = 0.0;
    PlanarPressureRegionOpeningFlowSettings settings;
    std::vector<PlanarPressureRegionOpeningFlow> openings;
    std::vector<PlanarPressureRegionOpeningBalance> regions;
    std::vector<PlanarPressureRegionOpeningComponent> components;
    std::size_t failedComponentCount = 0;
    std::size_t failedRegionCount = 0;
    double maximumAbsoluteOpeningNormalVelocityMetersPerSecond = 0.0;
    double maximumAbsoluteContinuityResidualCubicMeters = 0.0;
    double globalGeometryVolumeChangeCubicMeters = 0.0;
    double globalIntegratedOutwardRelativeVolumeCubicMeters = 0.0;
    double globalContinuityResidualCubicMeters = 0.0;
    bool allComponentsFeasible = false;
    bool allRegionsWithinTolerance = false;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionOpeningFlowAllocation&) const = default;
};

[[nodiscard]] PlanarPressureRegionOpeningFlowAllocation
solvePlanarPressureRegionOpeningFlow(
    const PlanarPressureRegionSweepLedger& sweep,
    std::span<const PlanarPressureRegionOpeningDefinition> openings,
    const PlanarPressureRegionOpeningFlowSettings& settings = {},
    const PlanarPressureRegionOpeningFlowLimits& limits = {});

void validatePlanarPressureRegionOpeningFlow(
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionSweepLedger& sweep,
    std::span<const PlanarPressureRegionOpeningDefinition> openings,
    const PlanarPressureRegionOpeningFlowLimits& limits = {});

} // namespace simwing::fsi::fluid
