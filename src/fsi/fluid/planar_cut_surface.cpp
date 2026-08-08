#include "fluid/planar_cut_surface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

[[nodiscard]] bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] Vector3 add(const Vector3& first, const Vector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

[[nodiscard]] Vector3 subtract(const Vector3& first,
                               const Vector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] Vector3 scale(const Vector3& value, const double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

[[nodiscard]] Vector3 cross(const Vector3& first, const Vector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

[[nodiscard]] double length(const Vector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

[[nodiscard]] double dot(const Vector3& first, const Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

[[nodiscard]] Vector3 axisVector(const GridFaceAxis axis,
                                 const double value) {
    switch (axis) {
    case GridFaceAxis::X:
        return {value, 0.0, 0.0};
    case GridFaceAxis::Y:
        return {0.0, value, 0.0};
    case GridFaceAxis::Z:
        return {0.0, 0.0, value};
    }
    throw std::invalid_argument(
        "planar cut surface has an unknown face axis");
}

[[nodiscard]] double normalCoordinate(const Vector3& value,
                                      const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return value.x;
    case GridFaceAxis::Y:
        return value.y;
    case GridFaceAxis::Z:
        return value.z;
    }
    throw std::invalid_argument(
        "planar cut surface has an unknown face axis");
}

void setNormalCoordinate(Vector3& value,
                         const GridFaceAxis axis,
                         const double coordinate) {
    switch (axis) {
    case GridFaceAxis::X:
        value.x = coordinate;
        return;
    case GridFaceAxis::Y:
        value.y = coordinate;
        return;
    case GridFaceAxis::Z:
        value.z = coordinate;
        return;
    }
    throw std::invalid_argument(
        "planar cut surface has an unknown face axis");
}

[[nodiscard]] std::size_t normalIndex(
    const MovingInterfaceFaceDiagnostics& face) {
    switch (face.axis) {
    case GridFaceAxis::X:
        return face.i;
    case GridFaceAxis::Y:
        return face.j;
    case GridFaceAxis::Z:
        return face.k;
    }
    throw std::invalid_argument(
        "planar cut surface has an unknown face axis");
}

[[nodiscard]] std::size_t faceLinearIndex(
    const GridCellCounts counts,
    const MovingInterfaceFaceDiagnostics& face) {
    if (face.i >= counts.x || face.j >= counts.y || face.k >= counts.z) {
        throw std::invalid_argument(
            "planar cut surface face index is outside the grid");
    }
    return face.i + counts.x * (face.j + counts.y * face.k);
}

[[nodiscard]] std::size_t transverseFaceCount(
    const GridCellCounts counts,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return counts.y * counts.z;
    case GridFaceAxis::Y:
        return counts.x * counts.z;
    case GridFaceAxis::Z:
        return counts.x * counts.y;
    }
    throw std::invalid_argument(
        "planar cut surface has an unknown face axis");
}

[[nodiscard]] std::pair<Vector3, Vector3> expectedFaceCorners(
    const PeriodicCartesianGrid& grid,
    const MovingInterfaceFaceDiagnostics& face) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    Vector3 first{
        lower.x + static_cast<double>(face.i) * spacing.x,
        lower.y + static_cast<double>(face.j) * spacing.y,
        lower.z + static_cast<double>(face.k) * spacing.z,
    };
    Vector3 second = first;
    switch (face.axis) {
    case GridFaceAxis::X:
        second.y += spacing.y;
        second.z += spacing.z;
        break;
    case GridFaceAxis::Y:
        second.x += spacing.x;
        second.z += spacing.z;
        break;
    case GridFaceAxis::Z:
        second.x += spacing.x;
        second.y += spacing.y;
        break;
    default:
        throw std::invalid_argument(
            "planar cut surface has an unknown face axis");
    }
    return {first, second};
}

[[nodiscard]] double combinedTolerance(const double absoluteTolerance,
                                       const double relativeTolerance,
                                       const double firstScale,
                                       const double secondScale) {
    return absoluteTolerance
        + relativeTolerance * std::max(firstScale, secondScale);
}

void validateSettings(const PlanarCutSurfacePressureSettings& settings) {
    const std::array nonnegative{
        settings.absolutePositionToleranceMeters,
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeAreaTolerance,
        settings.absoluteForceToleranceNewtons,
        settings.relativeForceTolerance,
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
    };
    if (!finite(settings.momentReferenceMeters)
        || !std::ranges::all_of(
            nonnegative,
            [](const double value) {
                return std::isfinite(value) && value >= 0.0;
            })) {
        throw std::invalid_argument(
            "planar cut-surface pressure settings are invalid");
    }
}

} // namespace

PlanarCutSurfacePressureDiagnostics evaluatePlanarCutSurfacePressure(
    const PeriodicCartesianGrid& grid,
    const PlanarMovingControlVolume& controlVolume,
    const MovingInterfaceProjectionDiagnostics& interfaceDiagnostics,
    const double surfaceOffsetMeters,
    const double physicalPlaneCoordinateMeters,
    const PlanarCutSurfacePressureSettings& settings) {
    validateSettings(settings);
    if (!controlVolume.matches(grid)) {
        throw std::invalid_argument(
            "planar cut surface does not match its control-volume grid");
    }
    if (interfaceDiagnostics.interfaceVersion
            != faceAlignedMovingInterfaceVersion
        || !interfaceDiagnostics.projection.converged
        || !interfaceDiagnostics.finite
        || interfaceDiagnostics.interfaceFaceCount
            != interfaceDiagnostics.faces.size()) {
        throw std::invalid_argument(
            "planar cut surface requires accepted interface diagnostics");
    }
    const double spacing = controlVolume.normalCellSpacingMeters();
    if (!std::isfinite(surfaceOffsetMeters)
        || !std::isfinite(physicalPlaneCoordinateMeters)
        || surfaceOffsetMeters < 0.0
        || surfaceOffsetMeters > spacing) {
        throw std::invalid_argument(
            "planar cut surface lies outside the current partial cell");
    }

    std::vector<const MovingInterfaceFaceDiagnostics*> selectedFaces;
    for (const auto& face : interfaceDiagnostics.faces) {
        if (face.surfaceStableId == controlVolume.movingSurfaceStableId()) {
            selectedFaces.push_back(&face);
        }
    }
    if (selectedFaces.empty()) {
        throw std::invalid_argument(
            "planar cut surface cannot find its fluid pressure reaction");
    }

    const GridFaceAxis axis = controlVolume.axis();
    const GridCellCounts counts = grid.cellCounts();
    if (selectedFaces.size() != interfaceDiagnostics.faces.size()
        || selectedFaces.size() != transverseFaceCount(counts, axis)) {
        throw std::invalid_argument(
            "planar cut surface requires one complete fluid face plane");
    }
    const double gridPlaneCoordinateMeters = normalCoordinate(
        selectedFaces.front()->lowerCornerMeters, axis);
    const double periodMeters = normalCoordinate(
        grid.upperMeters(), axis) - normalCoordinate(
        grid.lowerMeters(), axis);
    const double expectedPhysicalCoordinateMeters =
        gridPlaneCoordinateMeters + surfaceOffsetMeters;
    const double periodicDisplacement = physicalPlaneCoordinateMeters
        - expectedPhysicalCoordinateMeters;
    const double periodicImage = std::round(
        periodicDisplacement / periodMeters);

    PlanarCutSurfacePressureDiagnostics diagnostics;
    diagnostics.sourceInterfaceVersion =
        interfaceDiagnostics.interfaceVersion;
    diagnostics.surfaceStableId = controlVolume.movingSurfaceStableId();
    diagnostics.fluidRegionStableId = controlVolume.fluidRegionStableId();
    diagnostics.axis = axis;
    diagnostics.movingPlaneCoordinate =
        controlVolume.movingPlaneCoordinate();
    diagnostics.faceCount = selectedFaces.size();
    diagnostics.momentReferenceMeters = settings.momentReferenceMeters;
    diagnostics.surfaceOffsetMeters = surfaceOffsetMeters;
    diagnostics.gridPlaneCoordinateMeters = gridPlaneCoordinateMeters;
    diagnostics.physicalPlaneCoordinateMeters =
        physicalPlaneCoordinateMeters;
    diagnostics.periodicPositionResidualMeters = std::abs(
        periodicDisplacement - periodicImage * periodMeters);
    diagnostics.faces.reserve(selectedFaces.size());

    std::size_t previousFaceIndex = 0;
    bool havePreviousFaceIndex = false;
    for (const auto* source : selectedFaces) {
        const std::size_t currentFaceIndex = faceLinearIndex(counts, *source);
        const auto [expectedGridLower, expectedGridUpper] =
            expectedFaceCorners(grid, *source);
        const double tileArea = source->axis == GridFaceAxis::X
            ? grid.cellSpacingMeters().y * grid.cellSpacingMeters().z
            : (source->axis == GridFaceAxis::Y
                ? grid.cellSpacingMeters().x
                    * grid.cellSpacingMeters().z
                : grid.cellSpacingMeters().x
                    * grid.cellSpacingMeters().y);
        const double tileAreaTolerance = combinedTolerance(
            settings.absoluteAreaToleranceSquareMeters,
            settings.relativeAreaTolerance,
            source->areaSquareMeters, tileArea);
        const Vector3 axialTraction = axisVector(
            axis, normalCoordinate(source->pressureTractionPascals, axis));
        const Vector3 expectedForce = scale(
            source->pressureTractionPascals,
            source->areaSquareMeters);
        const double faceForceTolerance = combinedTolerance(
            settings.absoluteForceToleranceNewtons,
            settings.relativeForceTolerance,
            length(expectedForce), length(source->pressureForceNewtons));
        const double expectedPower = dot(
            source->pressureForceNewtons,
            axisVector(axis, source->normalVelocityMetersPerSecond));
        const double facePowerTolerance = combinedTolerance(
            settings.absolutePowerToleranceWatts,
            settings.relativePowerTolerance,
            std::abs(expectedPower), std::abs(source->pressurePowerWatts));
        if (source->minusRegionStableId
                != controlVolume.fluidRegionStableId()
            || source->plusRegionStableId
                != controlVolume.fluidRegionStableId()
            || source->axis != axis
            || normalIndex(*source)
                != controlVolume.movingPlaneCoordinate()
            || !finite(source->lowerCornerMeters)
            || !finite(source->upperCornerMeters)
            || !finite(source->pressureTractionPascals)
            || !finite(source->pressureForceNewtons)
            || !std::isfinite(source->areaSquareMeters)
            || source->areaSquareMeters <= 0.0
            || !std::isfinite(source->normalVelocityMetersPerSecond)
            || !std::isfinite(source->pressurePowerWatts)
            || (havePreviousFaceIndex
                && currentFaceIndex <= previousFaceIndex)
            || length(subtract(
                   source->lowerCornerMeters, expectedGridLower))
                > settings.absolutePositionToleranceMeters
            || length(subtract(
                   source->upperCornerMeters, expectedGridUpper))
                > settings.absolutePositionToleranceMeters
            || std::abs(source->areaSquareMeters - tileArea)
                > tileAreaTolerance
            || length(subtract(
                   source->pressureTractionPascals, axialTraction))
                > 0.0
            || length(subtract(
                   source->pressureForceNewtons, expectedForce))
                > faceForceTolerance
            || std::abs(source->pressurePowerWatts - expectedPower)
                > facePowerTolerance
            || std::abs(normalCoordinate(
                    source->lowerCornerMeters, axis)
                - gridPlaneCoordinateMeters)
                > settings.absolutePositionToleranceMeters
            || std::abs(normalCoordinate(
                    source->upperCornerMeters, axis)
                - gridPlaneCoordinateMeters)
                > settings.absolutePositionToleranceMeters) {
            throw std::invalid_argument(
                "fluid pressure faces do not match the planar cut surface");
        }
        previousFaceIndex = currentFaceIndex;
        havePreviousFaceIndex = true;
        Vector3 physicalLower = source->lowerCornerMeters;
        Vector3 physicalUpper = source->upperCornerMeters;
        setNormalCoordinate(
            physicalLower, axis, physicalPlaneCoordinateMeters);
        setNormalCoordinate(
            physicalUpper, axis, physicalPlaneCoordinateMeters);
        diagnostics.faces.push_back({
            source->surfaceStableId,
            source->minusRegionStableId,
            source->plusRegionStableId,
            source->axis,
            source->i,
            source->j,
            source->k,
            source->lowerCornerMeters,
            source->upperCornerMeters,
            physicalLower,
            physicalUpper,
            source->areaSquareMeters,
            source->normalVelocityMetersPerSecond,
            source->pressureTractionPascals,
            source->pressureForceNewtons,
            source->pressurePowerWatts,
        });
        diagnostics.areaSquareMeters += source->areaSquareMeters;
        diagnostics.pressureForceNewtons = add(
            diagnostics.pressureForceNewtons,
            source->pressureForceNewtons);
        const Vector3 physicalCenter = scale(
            add(physicalLower, physicalUpper), 0.5);
        diagnostics.pressureMomentNewtonMeters = add(
            diagnostics.pressureMomentNewtonMeters,
            cross(subtract(
                      physicalCenter,
                      settings.momentReferenceMeters),
                  source->pressureForceNewtons));
        diagnostics.pressurePowerWatts += source->pressurePowerWatts;
    }

    const auto surface = std::lower_bound(
        interfaceDiagnostics.surfaces.begin(),
        interfaceDiagnostics.surfaces.end(),
        controlVolume.movingSurfaceStableId(),
        [](const MovingInterfaceSurfaceDiagnostics& candidate,
           const std::uint64_t stableId) {
            return candidate.stableId < stableId;
        });
    if (surface == interfaceDiagnostics.surfaces.end()
        || surface->stableId != controlVolume.movingSurfaceStableId()
        || surface->faceCount != selectedFaces.size()
        || !std::isfinite(surface->areaSquareMeters)
        || !finite(surface->pressureForceNewtons)
        || !std::isfinite(surface->pressurePowerWatts)) {
        throw std::invalid_argument(
            "planar cut surface is missing its pressure aggregate");
    }
    diagnostics.sourceAreaSquareMeters = surface->areaSquareMeters;
    diagnostics.areaResidualSquareMeters =
        diagnostics.areaSquareMeters - surface->areaSquareMeters;
    diagnostics.sourcePressureForceNewtons =
        surface->pressureForceNewtons;
    diagnostics.forceResidualNewtons = subtract(
        diagnostics.pressureForceNewtons,
        surface->pressureForceNewtons);
    diagnostics.forceResidualNormNewtons = length(
        diagnostics.forceResidualNewtons);
    diagnostics.sourcePressurePowerWatts = surface->pressurePowerWatts;
    diagnostics.powerResidualWatts = diagnostics.pressurePowerWatts
        - surface->pressurePowerWatts;

    diagnostics.finite =
        std::isfinite(diagnostics.surfaceOffsetMeters)
        && std::isfinite(diagnostics.gridPlaneCoordinateMeters)
        && std::isfinite(diagnostics.physicalPlaneCoordinateMeters)
        && std::isfinite(diagnostics.periodicPositionResidualMeters)
        && std::isfinite(diagnostics.areaSquareMeters)
        && std::isfinite(diagnostics.sourceAreaSquareMeters)
        && std::isfinite(diagnostics.areaResidualSquareMeters)
        && finite(diagnostics.pressureForceNewtons)
        && finite(diagnostics.sourcePressureForceNewtons)
        && finite(diagnostics.forceResidualNewtons)
        && std::isfinite(diagnostics.forceResidualNormNewtons)
        && finite(diagnostics.pressureMomentNewtonMeters)
        && std::isfinite(diagnostics.pressurePowerWatts)
        && std::isfinite(diagnostics.sourcePressurePowerWatts)
        && std::isfinite(diagnostics.powerResidualWatts);
    const double areaTolerance = combinedTolerance(
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeAreaTolerance,
        diagnostics.areaSquareMeters,
        std::max(surface->areaSquareMeters,
                 controlVolume.crossSectionAreaSquareMeters()));
    const double forceTolerance = combinedTolerance(
        settings.absoluteForceToleranceNewtons,
        settings.relativeForceTolerance,
        length(diagnostics.pressureForceNewtons),
        length(surface->pressureForceNewtons));
    const double powerTolerance = combinedTolerance(
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
        std::abs(diagnostics.pressurePowerWatts),
        std::abs(surface->pressurePowerWatts));
    diagnostics.accepted = diagnostics.finite
        && diagnostics.periodicPositionResidualMeters
            <= settings.absolutePositionToleranceMeters
        && std::abs(diagnostics.areaSquareMeters
                    - controlVolume.crossSectionAreaSquareMeters())
            <= areaTolerance
        && std::abs(diagnostics.areaResidualSquareMeters) <= areaTolerance
        && diagnostics.forceResidualNormNewtons <= forceTolerance
        && std::abs(diagnostics.powerResidualWatts) <= powerTolerance;
    return diagnostics;
}

} // namespace simwing::fsi::fluid
