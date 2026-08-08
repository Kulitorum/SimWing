#include "fluid/porous_topology.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

struct AxisGeometry {
    std::size_t faceCount = 0;
    double spacingMeters = 0.0;
    double domainLengthMeters = 0.0;
    double wrappedFaceCenterMeters = 0.0;
};

AxisGeometry axisGeometry(const PeriodicCartesianGrid& grid,
                          const MovingPorousFaceTopology& topology) {
    if (topology.version != movingPorousFaceTopologyVersion) {
        throw std::invalid_argument(
            "moving porous topology version is unsupported");
    }
    const GridCellCounts counts = grid.cellCounts();
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    AxisGeometry result;
    switch (topology.axis) {
    case GridFaceAxis::X:
        result = {
            counts.x, spacing.x, upper.x - lower.x,
            topology.faceCoordinate < counts.x
                ? grid.xFaceCenterMeters(
                      topology.faceCoordinate, 0, 0).x
                : 0.0,
        };
        break;
    case GridFaceAxis::Y:
        result = {
            counts.y, spacing.y, upper.y - lower.y,
            topology.faceCoordinate < counts.y
                ? grid.yFaceCenterMeters(
                      0, topology.faceCoordinate, 0).y
                : 0.0,
        };
        break;
    case GridFaceAxis::Z:
        result = {
            counts.z, spacing.z, upper.z - lower.z,
            topology.faceCoordinate < counts.z
                ? grid.zFaceCenterMeters(
                      0, 0, topology.faceCoordinate).z
                : 0.0,
        };
        break;
    default:
        throw std::invalid_argument(
            "moving porous topology axis is invalid");
    }
    if (topology.faceCoordinate >= result.faceCount) {
        throw std::out_of_range(
            "moving porous topology face coordinate is outside the grid");
    }
    return result;
}

double rawCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& topology,
    const double physicalPlaneCoordinateMeters) {
    if (!std::isfinite(physicalPlaneCoordinateMeters)) {
        throw std::invalid_argument(
            "moving porous physical plane coordinate is non-finite");
    }
    const AxisGeometry geometry = axisGeometry(grid, topology);
    const double unwrappedFaceCenterMeters =
        geometry.wrappedFaceCenterMeters
        + static_cast<double>(topology.periodicImage)
            * geometry.domainLengthMeters;
    const double fraction = 0.5
        + (physicalPlaneCoordinateMeters - unwrappedFaceCenterMeters)
            / geometry.spacingMeters;
    if (!std::isfinite(fraction)) {
        throw std::invalid_argument(
            "moving porous crossing fraction is non-finite");
    }
    return fraction;
}

MovingPorousFaceTopology adjacentTopology(
    const PeriodicCartesianGrid& grid,
    MovingPorousFaceTopology topology,
    const PorousTopologyRebaseDirection direction) {
    const AxisGeometry geometry = axisGeometry(grid, topology);
    if (direction == PorousTopologyRebaseDirection::Positive) {
        if (topology.faceCoordinate + 1 == geometry.faceCount) {
            if (topology.periodicImage
                == std::numeric_limits<std::int64_t>::max()) {
                throw std::overflow_error(
                    "moving porous periodic image overflow");
            }
            topology.faceCoordinate = 0;
            ++topology.periodicImage;
        } else {
            ++topology.faceCoordinate;
        }
    } else if (direction == PorousTopologyRebaseDirection::Negative) {
        if (topology.faceCoordinate == 0) {
            if (topology.periodicImage
                == std::numeric_limits<std::int64_t>::min()) {
                throw std::overflow_error(
                    "moving porous periodic image underflow");
            }
            topology.faceCoordinate = geometry.faceCount - 1;
            --topology.periodicImage;
        } else {
            --topology.faceCoordinate;
        }
    }
    return topology;
}

} // namespace

double movingPorousCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& topology,
    const double physicalPlaneCoordinateMeters) {
    const double fraction = rawCrossingFraction(
        grid, topology, physicalPlaneCoordinateMeters);
    if (!(fraction > 0.0) || !(fraction < 1.0)) {
        throw std::runtime_error(
            "moving porous plane lies outside its topology segment");
    }
    return fraction;
}

MovingPorousTopologySelection selectMovingPorousTopology(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFaceTopology& current,
    const double physicalPlaneCoordinateMeters) {
    const double fraction = rawCrossingFraction(
        grid, current, physicalPlaneCoordinateMeters);
    if (fraction > 0.0 && fraction < 1.0) {
        return {current, fraction,
                PorousTopologyRebaseDirection::None};
    }

    PorousTopologyRebaseDirection direction;
    if (fraction >= 1.0 && fraction < 2.0) {
        direction = PorousTopologyRebaseDirection::Positive;
    } else if (fraction <= 0.0 && fraction > -1.0) {
        direction = PorousTopologyRebaseDirection::Negative;
    } else {
        throw std::runtime_error(
            "moving porous plane crossed more than one topology segment");
    }
    const MovingPorousFaceTopology rebased = adjacentTopology(
        grid, current, direction);
    return {
        rebased,
        movingPorousCrossingFraction(
            grid, rebased, physicalPlaneCoordinateMeters),
        direction,
    };
}

} // namespace simwing::fsi::fluid
