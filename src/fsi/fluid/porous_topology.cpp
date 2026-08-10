#include "fluid/porous_topology.h"

namespace simwing::fsi::fluid {

double movingPorousCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& topology,
    const double physicalPlaneCoordinateMeters) {
    return movingPlanarCrossingFraction(
        grid, topology, physicalPlaneCoordinateMeters);
}

MovingPorousTopologySelection selectMovingPorousTopology(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& current,
    const double physicalPlaneCoordinateMeters) {
    return selectMovingPlanarTopology(
        grid, current, physicalPlaneCoordinateMeters);
}

} // namespace simwing::fsi::fluid
