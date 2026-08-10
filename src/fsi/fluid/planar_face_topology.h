#pragma once

#include "fluid/interface_jump.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t movingPlanarFaceTopologyVersion = 1;

// One unwrapped axis-aligned plane epoch on a periodic MAC grid. The face
// coordinate owns the wrapped fluid topology; periodicImage disambiguates the
// physical plane after a domain wrap without changing any field index.
struct MovingPlanarFaceTopology {
    std::uint32_t version = movingPlanarFaceTopologyVersion;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t faceCoordinate = 0;
    std::int64_t periodicImage = 0;

    bool operator==(const MovingPlanarFaceTopology&) const = default;
};

enum class PlanarTopologyRebaseDirection : std::int8_t {
    Negative = -1,
    None = 0,
    Positive = 1,
};

struct MovingPlanarTopologySelection {
    MovingPlanarFaceTopology topology;
    double crossingFraction = 0.5;
    PlanarTopologyRebaseDirection rebaseDirection =
        PlanarTopologyRebaseDirection::None;

    bool operator==(const MovingPlanarTopologySelection&) const = default;
};

// Returns the strict (0, 1) crossing coordinate for a physical plane already
// owned by the supplied epoch. Exact MAC-plane placement is intentionally not
// representable by the sharp-jump stencil and is rejected.
[[nodiscard]] double movingPlanarCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPlanarFaceTopology& topology,
    double physicalPlaneCoordinateMeters);

// Selects the current, immediately previous, or immediately next dual-cell
// segment for one physical plane. It is pure and supports X/Y/Z, both motion
// directions, and periodic wraps. Crossing more than one segment in a single
// call is rejected so callers cannot silently skip topology events.
[[nodiscard]] MovingPlanarTopologySelection selectMovingPlanarTopology(
    const PeriodicCartesianGrid& grid,
    const MovingPlanarFaceTopology& current,
    double physicalPlaneCoordinateMeters);

} // namespace simwing::fsi::fluid
