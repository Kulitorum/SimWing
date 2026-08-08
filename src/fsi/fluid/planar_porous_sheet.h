#pragma once

#include "fluid/porous_interface.h"
#include "fluid/porous_topology.h"

#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

// Authored data for one complete, axis-aligned porous plane in a periodic
// topology epoch. The physical coordinate is unwrapped; topology owns the
// wrapped MAC face and signed periodic image used to locate it.
struct PlanarPorousSheetDefinition {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    MovingPorousFaceTopology topology;
    double physicalPlaneCoordinateMeters = 0.0;
    double surfaceNormalVelocityMetersPerSecond = 0.0;
    DarcyForchheimerResistance resistance;

    bool operator==(const PlanarPorousSheetDefinition&) const = default;
};

// Expands one topology-bound physical sheet into exactly one porous crossing
// for every transverse MAC tile. Output order is deterministic and canonical
// for X, Y, and Z. Invalid identity, material, kinematics, topology, or exact
// segment-boundary placement is rejected before an output vector is returned.
[[nodiscard]] std::vector<PorousGridFaceCrossing>
makePlanarPorousSheetCrossings(
    const PeriodicCartesianGrid& grid,
    const PlanarPorousSheetDefinition& definition);

} // namespace simwing::fsi::fluid
