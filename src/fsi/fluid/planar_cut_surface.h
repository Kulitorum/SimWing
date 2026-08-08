#pragma once

#include "fluid/moving_control_volume.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarCutSurfacePressureVersion = 2;

struct PlanarCutSurfacePressureSettings {
    Vector3 momentReferenceMeters;
    double absolutePositionToleranceMeters = 1.0e-12;
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-11;
    double relativeVelocityTolerance = 1.0e-11;
    double absoluteAreaToleranceSquareMeters = 1.0e-11;
    double relativeAreaTolerance = 1.0e-11;
    double absoluteForceToleranceNewtons = 1.0e-10;
    double relativeForceTolerance = 1.0e-11;
    double absolutePowerToleranceWatts = 1.0e-10;
    double relativePowerTolerance = 1.0e-11;
};

struct PlanarCutSurfacePressureFaceDiagnostics {
    std::uint64_t surfaceStableId = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    Vector3 gridLowerCornerMeters;
    Vector3 gridUpperCornerMeters;
    Vector3 physicalLowerCornerMeters;
    Vector3 physicalUpperCornerMeters;
    double areaSquareMeters = 0.0;
    double normalVelocityMetersPerSecond = 0.0;
    Vector3 pressureTractionPascals;
    Vector3 pressureForceNewtons;
    double pressurePowerWatts = 0.0;

    bool operator==(
        const PlanarCutSurfacePressureFaceDiagnostics&) const = default;
};

// Accepted physical geometry and reaction ledger for one planar moving
// control-volume surface. The pressure reaction is the face-aligned constraint
// reaction produced by the projection; this operator translates its application
// geometry from the Eulerian MAC plane to the bounded physical cut plane. It
// does not interpolate a new cell pressure or claim general cut-cell metrics.
struct PlanarCutSurfacePressureDiagnostics {
    std::uint32_t version = planarCutSurfacePressureVersion;
    std::uint32_t sourceInterfaceVersion = 0;
    std::uint64_t surfaceStableId = 0;
    std::uint64_t fluidRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t movingPlaneCoordinate = 0;
    std::size_t faceCount = 0;
    Vector3 momentReferenceMeters;
    double surfaceOffsetMeters = 0.0;
    double gridPlaneCoordinateMeters = 0.0;
    double physicalPlaneCoordinateMeters = 0.0;
    double periodicPositionResidualMeters = 0.0;
    double normalVelocityMetersPerSecond = 0.0;
    double maximumNormalVelocitySpreadMetersPerSecond = 0.0;
    bool kinematicsResampled = false;
    double reactionSourcePhysicalPlaneCoordinateMeters = 0.0;
    double reactionSourceNormalVelocityMetersPerSecond = 0.0;
    double areaSquareMeters = 0.0;
    double sourceAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;
    Vector3 pressureForceNewtons;
    Vector3 sourcePressureForceNewtons;
    Vector3 forceResidualNewtons;
    double forceResidualNormNewtons = 0.0;
    Vector3 pressureMomentNewtonMeters;
    double pressurePowerWatts = 0.0;
    double sourcePressurePowerWatts = 0.0;
    double powerResidualWatts = 0.0;
    std::vector<PlanarCutSurfacePressureFaceDiagnostics> faces;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PlanarCutSurfacePressureDiagnostics&) const = default;
};

// Bounded planar physical-surface evaluation for the current control-volume
// epoch. surfaceOffsetMeters must lie in the positive-axis partial cell.
// physicalPlaneCoordinateMeters may be an unwrapped periodic image, but it must
// be congruent with grid-plane + offset. The function owns no mutable state.
[[nodiscard]] PlanarCutSurfacePressureDiagnostics
evaluatePlanarCutSurfacePressure(
    const PeriodicCartesianGrid& grid,
    const PlanarMovingControlVolume& controlVolume,
    const MovingInterfaceProjectionDiagnostics& interfaceDiagnostics,
    double surfaceOffsetMeters,
    double physicalPlaneCoordinateMeters,
    const PlanarCutSurfacePressureSettings& settings = {});

// Retains an already accepted face-resolved pressure reaction while sampling
// its power on another congruent physical plane and rigid normal velocity in
// the same topology epoch. This is the explicit temporal adapter for a
// projection reaction that represents a macro-step average rather than an
// endpoint force. Traction and force remain unchanged; geometry and power are
// fully revalidated through evaluatePlanarCutSurfacePressure.
[[nodiscard]] PlanarCutSurfacePressureDiagnostics
resamplePlanarCutSurfaceReaction(
    const PeriodicCartesianGrid& grid,
    const PlanarMovingControlVolume& controlVolume,
    const PlanarCutSurfacePressureDiagnostics& acceptedReaction,
    double surfaceOffsetMeters,
    double physicalPlaneCoordinateMeters,
    double normalVelocityMetersPerSecond,
    const PlanarCutSurfacePressureSettings& settings = {});

} // namespace simwing::fsi::fluid
