#pragma once

#include "fluid/interface_jump.h"
#include "fluid/planar_face_topology.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

// One complete axis-aligned pressure-jump plane in an unwrapped periodic
// topology epoch. Multiple definitions form a globally ordered region chain;
// a closed chain can represent a thin pocket even when several layers occupy
// one cell-centre segment.
struct PlanarPressureJumpLayerDefinition {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    MovingPlanarFaceTopology topology;
    double physicalPlaneCoordinateMeters = 0.0;
    double pressureJumpPascals = 0.0;

    bool operator==(const PlanarPressureJumpLayerDefinition&) const = default;
};

struct PlanarPressureJumpLayerTranslation {
    std::vector<PlanarPressureJumpLayerDefinition> layers;
    std::vector<PlanarTopologyRebaseDirection> rebaseDirections;

    bool operator==(const PlanarPressureJumpLayerTranslation&) const = default;
};

// One exact interval between consecutive complete pressure layers in an
// unwrapped one-period window. The last interval ends at the first layer in
// the next periodic image. Its static pressure is a region potential, not a
// cell-centre interpolation.
struct PlanarPressureRegionInterval {
    std::uint64_t lowerSurfaceStableId = 0;
    std::uint64_t upperSurfaceStableId = 0;
    std::uint64_t regionStableId = 0;
    double lowerCoordinateMeters = 0.0;
    double upperCoordinateMeters = 0.0;
    double volumeCubicMeters = 0.0;
    double pressurePascals = 0.0;

    bool operator==(const PlanarPressureRegionInterval&) const = default;
};

struct PlanarPressureRegionSummary {
    std::uint64_t regionStableId = 0;
    double volumeCubicMeters = 0.0;
    double pressurePascals = 0.0;

    bool operator==(const PlanarPressureRegionSummary&) const = default;
};

struct StaticPlanarPressureRegionProfile {
    GridFaceAxis axis = GridFaceAxis::X;
    double windowLowerCoordinateMeters = 0.0;
    double windowUpperCoordinateMeters = 0.0;
    double geometricDomainVolumeCubicMeters = 0.0;
    double intervalVolumeCubicMeters = 0.0;
    double volumeClosureResidualCubicMeters = 0.0;
    double requestedVolumeMeanPressurePascals = 0.0;
    double achievedVolumeMeanPressurePascals = 0.0;
    std::vector<PlanarPressureRegionInterval> intervals;
    std::vector<PlanarPressureRegionSummary> regions;

    bool operator==(const StaticPlanarPressureRegionProfile&) const = default;
};

// Canonicalizes one closed periodic region chain by unwrapped physical
// coordinate and expands every layer into a complete transverse plane. Every
// definition must use the same axis, unique stable surface identity, strict
// topology segment ownership, and a span smaller than one domain period.
[[nodiscard]] SharpPressureJumpField makePlanarPressureJumpField(
    const PeriodicCartesianGrid& grid,
    std::span<const PlanarPressureJumpLayerDefinition> layers);

// Rigidly translates a complete layer chain. Each layer may retain or rebase
// by exactly one dual-cell segment, including a periodic wrap. The current and
// candidate complete fields are both validated before the owning result is
// returned, so invalid, exact-boundary, or skipped-segment motion publishes no
// partial topology.
[[nodiscard]] PlanarPressureJumpLayerTranslation
translatePlanarPressureJumpLayers(
    const PeriodicCartesianGrid& grid,
    std::span<const PlanarPressureJumpLayerDefinition> currentLayers,
    double displacementMeters);

// Reconstructs the exact axis-aligned interval volumes and one consistent
// piecewise-constant static pressure potential from a complete closed layer
// chain. The jump sum around the period and every repeated region potential
// must close exactly. A single gauge shift gives the requested volume-weighted
// mean pressure. This is a verification oracle only: it creates no regional
// velocity degrees of freedom and does not change the dense projection.
[[nodiscard]] StaticPlanarPressureRegionProfile
makeStaticPlanarPressureRegionProfile(
    const PeriodicCartesianGrid& grid,
    std::span<const PlanarPressureJumpLayerDefinition> layers,
    double volumeMeanPressurePascals = 0.0);

} // namespace simwing::fsi::fluid
