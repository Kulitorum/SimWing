#pragma once

#include "fluid/interface_jump.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t movingPorousFaceTopologyVersion = 1;

// One unwrapped planar porous-sheet epoch on a periodic MAC grid. The face
// coordinate owns the wrapped fluid topology; periodicImage disambiguates the
// physical plane after a domain wrap without changing any field index.
struct MovingPorousFaceTopology {
    std::uint32_t version = movingPorousFaceTopologyVersion;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t faceCoordinate = 0;
    std::int64_t periodicImage = 0;

    bool operator==(const MovingPorousFaceTopology&) const = default;
};

enum class PorousTopologyRebaseDirection : std::int8_t {
    Negative = -1,
    None = 0,
    Positive = 1,
};

struct MovingPorousTopologySelection {
    MovingPorousFaceTopology topology;
    double crossingFraction = 0.5;
    PorousTopologyRebaseDirection rebaseDirection =
        PorousTopologyRebaseDirection::None;

    bool operator==(const MovingPorousTopologySelection&) const = default;
};

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
