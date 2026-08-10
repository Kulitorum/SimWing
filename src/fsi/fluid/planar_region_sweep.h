#pragma once

#include "fluid/planar_pressure_jump.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarPressureRegionSweepVersion = 1;

struct PlanarPressureRegionSweepLimits {
    std::size_t maximumLayers = 1'000'000;
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumOwnedBytes = 256ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionIntervalSweep {
    std::uint64_t lowerSurfaceStableId = 0;
    std::uint64_t upperSurfaceStableId = 0;
    std::uint64_t regionStableId = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double lowerSurfaceDisplacementMeters = 0.0;
    double upperSurfaceDisplacementMeters = 0.0;
    double lowerSurfaceVelocityMetersPerSecond = 0.0;
    double upperSurfaceVelocityMetersPerSecond = 0.0;
    double boundarySweptVolumeCubicMeters = 0.0;
    double surfaceGeometryResidualCubicMeters = 0.0;

    bool operator==(const PlanarPressureRegionIntervalSweep&) const = default;
};

struct PlanarPressureRegionSweepSummary {
    std::uint64_t regionStableId = 0;
    double previousVolumeCubicMeters = 0.0;
    double currentVolumeCubicMeters = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double boundarySweptVolumeCubicMeters = 0.0;
    double surfaceGeometryResidualCubicMeters = 0.0;

    bool operator==(const PlanarPressureRegionSweepSummary&) const = default;
};

// Owning two-epoch geometric conservation ledger for one complete
// axis-aligned layer chain. Each endpoint has a closed static regional pressure
// potential, stable layer identity and strict topology placement. Every layer
// may move through at most one topology segment. This owns boundary kinematics
// only: it supplies no Eulerian regional flux or velocity unknown.
struct PlanarPressureRegionSweepLedger {
    std::uint32_t version = planarPressureRegionSweepVersion;
    GridFaceAxis axis = GridFaceAxis::X;
    double durationSeconds = 0.0;
    double crossSectionAreaSquareMeters = 0.0;
    StaticPlanarPressureRegionProfile previousProfile;
    StaticPlanarPressureRegionProfile currentProfile;
    std::vector<PlanarPressureRegionIntervalSweep> intervals;
    std::vector<PlanarPressureRegionSweepSummary> regions;
    double maximumAbsoluteSurfaceGeometryResidualCubicMeters = 0.0;
    double globalGeometryVolumeChangeCubicMeters = 0.0;
    double globalBoundarySweptVolumeCubicMeters = 0.0;
    double globalSurfaceGeometryResidualCubicMeters = 0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(const PlanarPressureRegionSweepLedger&) const = default;
};

// Builds a pure candidate ledger without mutating either authored endpoint.
// Layer order, stable identities, region sides and pressure jumps must remain
// unchanged. Exact-boundary, skipped-segment, crossing, or foreign-topology
// motion rejects before a result is published.
[[nodiscard]] PlanarPressureRegionSweepLedger
makePlanarPressureRegionSweepLedger(
    const PeriodicCartesianGrid& grid,
    std::span<const PlanarPressureJumpLayerDefinition> previousLayers,
    std::span<const PlanarPressureJumpLayerDefinition> currentLayers,
    double durationSeconds,
    const PlanarPressureRegionSweepLimits& limits = {});

} // namespace simwing::fsi::fluid
