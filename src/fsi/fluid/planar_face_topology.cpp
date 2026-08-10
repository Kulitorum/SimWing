#include "fluid/planar_face_topology.h"

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
                          const MovingPlanarFaceTopology& topology) {
    if (topology.version != movingPlanarFaceTopologyVersion) {
        throw std::invalid_argument(
            "moving planar topology version is unsupported");
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
            "moving planar topology axis is invalid");
    }
    if (topology.faceCoordinate >= result.faceCount) {
        throw std::out_of_range(
            "moving planar topology face coordinate is outside the grid");
    }
    return result;
}

double rawCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPlanarFaceTopology& topology,
    const double physicalPlaneCoordinateMeters) {
    if (!std::isfinite(physicalPlaneCoordinateMeters)) {
        throw std::invalid_argument(
            "moving planar physical plane coordinate is non-finite");
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
            "moving planar crossing fraction is non-finite");
    }
    return fraction;
}

MovingPlanarFaceTopology adjacentTopology(
    const PeriodicCartesianGrid& grid,
    MovingPlanarFaceTopology topology,
    const PlanarTopologyRebaseDirection direction) {
    const AxisGeometry geometry = axisGeometry(grid, topology);
    if (direction == PlanarTopologyRebaseDirection::Positive) {
        if (topology.faceCoordinate + 1 == geometry.faceCount) {
            if (topology.periodicImage
                == std::numeric_limits<std::int64_t>::max()) {
                throw std::overflow_error(
                    "moving planar periodic image overflow");
            }
            topology.faceCoordinate = 0;
            ++topology.periodicImage;
        } else {
            ++topology.faceCoordinate;
        }
    } else if (direction == PlanarTopologyRebaseDirection::Negative) {
        if (topology.faceCoordinate == 0) {
            if (topology.periodicImage
                == std::numeric_limits<std::int64_t>::min()) {
                throw std::overflow_error(
                    "moving planar periodic image underflow");
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

double movingPlanarCrossingFraction(
    const PeriodicCartesianGrid& grid,
    const MovingPlanarFaceTopology& topology,
    const double physicalPlaneCoordinateMeters) {
    const double fraction = rawCrossingFraction(
        grid, topology, physicalPlaneCoordinateMeters);
    if (!(fraction > 0.0) || !(fraction < 1.0)) {
        throw std::runtime_error(
            "moving planar plane lies outside its topology segment");
    }
    return fraction;
}

MovingPlanarTopologySelection selectMovingPlanarTopology(
    const PeriodicCartesianGrid& grid,
    const MovingPlanarFaceTopology& current,
    const double physicalPlaneCoordinateMeters) {
    const double fraction = rawCrossingFraction(
        grid, current, physicalPlaneCoordinateMeters);
    if (fraction > 0.0 && fraction < 1.0) {
        return {current, fraction, PlanarTopologyRebaseDirection::None};
    }

    PlanarTopologyRebaseDirection direction;
    if (fraction >= 1.0 && fraction < 2.0) {
        direction = PlanarTopologyRebaseDirection::Positive;
    } else if (fraction <= 0.0 && fraction > -1.0) {
        direction = PlanarTopologyRebaseDirection::Negative;
    } else {
        throw std::runtime_error(
            "moving planar plane crossed more than one topology segment");
    }
    const MovingPlanarFaceTopology rebased = adjacentTopology(
        grid, current, direction);
    return {
        rebased,
        movingPlanarCrossingFraction(
            grid, rebased, physicalPlaneCoordinateMeters),
        direction,
    };
}

} // namespace simwing::fsi::fluid
