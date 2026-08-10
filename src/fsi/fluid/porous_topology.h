#pragma once

#include "fluid/planar_face_topology.h"

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t movingPorousFaceTopologyVersion =
    movingPlanarFaceTopologyVersion;

// One unwrapped planar porous-sheet epoch on a periodic MAC grid. The face
// coordinate owns the wrapped fluid topology; periodicImage disambiguates the
// physical plane after a domain wrap without changing any field index.
using MovingPorousFaceTopology = MovingPlanarFaceTopology;

using PorousTopologyRebaseDirection = PlanarTopologyRebaseDirection;

using MovingPorousTopologySelection = MovingPlanarTopologySelection;

// Returns the strict (0, 1) crossing coordinate for a physical plane already
// owned by the supplied epoch. Exact MAC-plane placement is intentionally not
// representable by the sharp-jump stencil and is rejected.
[[nodiscard]] double movingPorousCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& topology,
    double physicalPlaneCoordinateMeters);

// Selects the current, immediately previous, or immediately next dual-cell
// segment for one physical plane. It is pure and supports X/Y/Z, both motion
// directions, and periodic wraps. Crossing more than one segment in a single
// call is rejected so callers cannot silently skip topology events.
[[nodiscard]] MovingPorousTopologySelection selectMovingPorousTopology(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& current,
    double physicalPlaneCoordinateMeters);

} // namespace simwing::fsi::fluid
