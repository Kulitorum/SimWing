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

} // namespace simwing::fsi::fluid
