#include "fluid/moving_control_volume.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

[[nodiscard]] std::size_t faceIndex(
    const GridCellCounts counts,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) {
    return i + counts.x * (j + counts.y * k);
}

[[nodiscard]] std::size_t axisCount(
    const GridCellCounts counts,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return counts.x;
    case GridFaceAxis::Y:
        return counts.y;
    case GridFaceAxis::Z:
        return counts.z;
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] std::size_t faceCoordinate(
    const GridFaceMovingInterface& face) {
    switch (face.axis) {
    case GridFaceAxis::X:
        return face.i;
    case GridFaceAxis::Y:
        return face.j;
    case GridFaceAxis::Z:
        return face.k;
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] double normalSpacing(
    const Vector3 spacing,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return spacing.x;
    case GridFaceAxis::Y:
        return spacing.y;
    case GridFaceAxis::Z:
        return spacing.z;
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] double tileArea(
    const Vector3 spacing,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return spacing.y * spacing.z;
    case GridFaceAxis::Y:
        return spacing.z * spacing.x;
    case GridFaceAxis::Z:
        return spacing.x * spacing.y;
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] std::size_t transverseFaceCount(
    const GridCellCounts counts,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return counts.y * counts.z;
    case GridFaceAxis::Y:
        return counts.z * counts.x;
    case GridFaceAxis::Z:
        return counts.x * counts.y;
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] std::span<const std::uint8_t> constraints(
    const FaceAlignedMovingInterface& interfaces,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return interfaces.xFaceConstraints();
    case GridFaceAxis::Y:
        return interfaces.yFaceConstraints();
    case GridFaceAxis::Z:
        return interfaces.zFaceConstraints();
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] std::span<const double> velocities(
    const MacVelocityField& velocity,
    const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return velocity.xFaces();
    case GridFaceAxis::Y:
        return velocity.yFaces();
    case GridFaceAxis::Z:
        return velocity.zFaces();
    }
    throw std::invalid_argument(
        "planar control volume has an unknown face axis");
}

[[nodiscard]] double combinedTolerance(
    const double absoluteTolerance,
    const double relativeTolerance,
    const double firstMagnitude,
    const double secondMagnitude) {
    return absoluteTolerance
        + relativeTolerance * std::max(firstMagnitude, secondMagnitude);
}

void validateSettings(const PlanarControlVolumeSettings& settings) {
    const double nonnegative[] = {
        settings.absoluteVelocityToleranceMetersPerSecond,
        settings.relativeVelocityTolerance,
        settings.absolutePositionToleranceMeters,
        settings.absoluteAreaToleranceSquareMeters,
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance,
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
        settings.minimumRemainingCellLengthMeters,
    };
    if (!std::ranges::all_of(nonnegative, [](const double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument(
            "planar control-volume settings are invalid");
    }
}

[[nodiscard]] bool sameGrid(
    const PeriodicCartesianGrid& grid,
    const GridCellCounts counts,
    const Vector3 lower,
    const Vector3 upper) {
    return grid.cellCounts() == counts
        && grid.lowerMeters() == lower
        && grid.upperMeters() == upper;
}

} // namespace

PlanarMovingControlVolume::PlanarMovingControlVolume(
    const PeriodicCartesianGrid& grid,
    const FaceAlignedMovingInterface& interfaces,
    const std::uint64_t movingSurfaceStableId,
    const std::size_t openingPlaneCoordinate)
    : cellCounts_(grid.cellCounts()),
      lowerMeters_(grid.lowerMeters()),
      upperMeters_(grid.upperMeters()),
      movingSurfaceStableId_(movingSurfaceStableId),
      openingPlaneCoordinate_(openingPlaneCoordinate) {
    if (!interfaces.matches(grid) || movingSurfaceStableId_ == 0) {
        throw std::invalid_argument(
            "planar control volume requires matching grid topology and a surface ID");
    }
    std::vector<GridFaceMovingInterface> selectedFaces;
    for (const auto& face : interfaces.faces()) {
        if (face.surfaceStableId == movingSurfaceStableId_) {
            selectedFaces.push_back(face);
        }
    }
    if (selectedFaces.empty()) {
        throw std::invalid_argument(
            "planar control volume cannot find its moving surface");
    }
    if (selectedFaces.size() != interfaces.faceCount()
        || interfaces.regionCount() != 1) {
        throw std::invalid_argument(
            "planar control volume requires one nonseparating surface and one fluid region");
    }

    axis_ = selectedFaces.front().axis;
    movingPlaneCoordinate_ = faceCoordinate(selectedFaces.front());
    fluidRegionStableId_ = selectedFaces.front().minusRegionStableId;
    const std::size_t normalCount = axisCount(cellCounts_, axis_);
    if (openingPlaneCoordinate_ >= normalCount
        || openingPlaneCoordinate_ == movingPlaneCoordinate_) {
        throw std::invalid_argument(
            "planar control-volume opening coordinate is invalid");
    }
    const std::size_t expectedFaceCount = transverseFaceCount(
        cellCounts_, axis_);
    if (selectedFaces.size() != expectedFaceCount) {
        throw std::invalid_argument(
            "planar moving surface does not cover one complete grid plane");
    }
    surfaceFaces_.reserve(selectedFaces.size());
    for (const auto& face : selectedFaces) {
        if (face.axis != axis_
            || faceCoordinate(face) != movingPlaneCoordinate_
            || face.minusRegionStableId != fluidRegionStableId_
            || face.plusRegionStableId != fluidRegionStableId_) {
            throw std::invalid_argument(
                "planar moving surface is not one nonseparating fluid plane");
        }
        surfaceFaces_.push_back({
            face.minusRegionStableId,
            face.plusRegionStableId,
            face.i,
            face.j,
            face.k,
        });
    }

    const auto blocked = constraints(interfaces, axis_);
    openingFaceIndices_.reserve(expectedFaceCount);
    for (std::size_t k = 0; k < cellCounts_.z; ++k) {
        for (std::size_t j = 0; j < cellCounts_.y; ++j) {
            for (std::size_t i = 0; i < cellCounts_.x; ++i) {
                const std::size_t coordinate = axis_ == GridFaceAxis::X
                    ? i : (axis_ == GridFaceAxis::Y ? j : k);
                if (coordinate != openingPlaneCoordinate_) {
                    continue;
                }
                const std::size_t index = faceIndex(cellCounts_, i, j, k);
                if (blocked[index] != 0) {
                    throw std::invalid_argument(
                        "planar control-volume opening is constrained");
                }
                openingFaceIndices_.push_back(index);
            }
        }
    }
    if (openingFaceIndices_.size() != expectedFaceCount) {
        throw std::logic_error(
            "planar control-volume opening enumeration failed");
    }

    const Vector3 spacing = grid.cellSpacingMeters();
    normalCellSpacingMeters_ = normalSpacing(spacing, axis_);
    tileAreaSquareMeters_ = tileArea(spacing, axis_);
    crossSectionAreaSquareMeters_ = tileAreaSquareMeters_
        * static_cast<double>(expectedFaceCount);
    const std::size_t layerCount =
        (movingPlaneCoordinate_ + normalCount - openingPlaneCoordinate_)
        % normalCount;
    if (layerCount == 0) {
        throw std::invalid_argument(
            "planar control volume has no full reference layer");
    }
    referenceVolumeCubicMeters_ = crossSectionAreaSquareMeters_
        * normalCellSpacingMeters_ * static_cast<double>(layerCount);
}

std::uint64_t
PlanarMovingControlVolume::movingSurfaceStableId() const noexcept {
    return movingSurfaceStableId_;
}

std::uint64_t
PlanarMovingControlVolume::fluidRegionStableId() const noexcept {
    return fluidRegionStableId_;
}

bool PlanarMovingControlVolume::matches(
    const PeriodicCartesianGrid& grid) const noexcept {
    return sameGrid(grid, cellCounts_, lowerMeters_, upperMeters_);
}

GridFaceAxis PlanarMovingControlVolume::axis() const noexcept {
    return axis_;
}

std::size_t
PlanarMovingControlVolume::movingPlaneCoordinate() const noexcept {
    return movingPlaneCoordinate_;
}

std::size_t
PlanarMovingControlVolume::openingPlaneCoordinate() const noexcept {
    return openingPlaneCoordinate_;
}

double PlanarMovingControlVolume::crossSectionAreaSquareMeters()
    const noexcept {
    return crossSectionAreaSquareMeters_;
}

double PlanarMovingControlVolume::normalCellSpacingMeters()
    const noexcept {
    return normalCellSpacingMeters_;
}

double PlanarMovingControlVolume::referenceVolumeCubicMeters()
    const noexcept {
    return referenceVolumeCubicMeters_;
}

PlanarControlVolumeDiagnostics PlanarMovingControlVolume::evaluate(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& projectedVelocityMetersPerSecond,
    const MovingInterfaceProjectionDiagnostics& interfaceDiagnostics,
    const PlanarControlVolumeStep& step,
    const PlanarControlVolumeSettings& settings) const {
    validateSettings(settings);
    if (!sameGrid(grid, cellCounts_, lowerMeters_, upperMeters_)
        || !projectedVelocityMetersPerSecond.matches(grid)
        || !isFinite(projectedVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "planar control-volume fields do not match their bound grid");
    }
    if (interfaceDiagnostics.interfaceVersion
            != faceAlignedMovingInterfaceVersion
        || !interfaceDiagnostics.projection.converged
        || !interfaceDiagnostics.finite
        || interfaceDiagnostics.interfaceFaceCount
            != interfaceDiagnostics.faces.size()
        || interfaceDiagnostics.faces.size() != surfaceFaces_.size()) {
        throw std::invalid_argument(
            "planar control volume requires accepted interface diagnostics");
    }
    if (!std::isfinite(step.startSurfaceOffsetMeters)
        || !std::isfinite(step.endSurfaceOffsetMeters)
        || !std::isfinite(step.durationSeconds)
        || step.startSurfaceOffsetMeters < 0.0
        || step.endSurfaceOffsetMeters < 0.0
        || !(step.durationSeconds > 0.0)
        || step.startSurfaceOffsetMeters >= normalCellSpacingMeters_
        || step.startSurfaceOffsetMeters
            > normalCellSpacingMeters_
                - settings.minimumRemainingCellLengthMeters
        || (step.endsAtCellBoundary
                ? step.endSurfaceOffsetMeters
                    != normalCellSpacingMeters_
                : step.endSurfaceOffsetMeters >= normalCellSpacingMeters_
                    || step.endSurfaceOffsetMeters
                        > normalCellSpacingMeters_
                            - settings.minimumRemainingCellLengthMeters)) {
        throw std::invalid_argument(
            "planar control-volume step leaves its supported cut-cell interval");
    }

    PlanarControlVolumeDiagnostics diagnostics;
    diagnostics.movingSurfaceStableId = movingSurfaceStableId_;
    diagnostics.fluidRegionStableId = fluidRegionStableId_;
    diagnostics.axis = axis_;
    diagnostics.movingPlaneCoordinate = movingPlaneCoordinate_;
    diagnostics.openingPlaneCoordinate = openingPlaneCoordinate_;
    diagnostics.movingSurfaceFaceCount = surfaceFaces_.size();
    diagnostics.openingFaceCount = openingFaceIndices_.size();
    diagnostics.crossSectionAreaSquareMeters =
        crossSectionAreaSquareMeters_;
    diagnostics.referenceVolumeCubicMeters = referenceVolumeCubicMeters_;
    diagnostics.startCutCellVolumeCubicMeters =
        crossSectionAreaSquareMeters_ * step.startSurfaceOffsetMeters;
    diagnostics.endCutCellVolumeCubicMeters =
        crossSectionAreaSquareMeters_ * step.endSurfaceOffsetMeters;
    diagnostics.startCutCellVolumeFraction =
        step.startSurfaceOffsetMeters / normalCellSpacingMeters_;
    diagnostics.endCutCellVolumeFraction =
        step.endSurfaceOffsetMeters / normalCellSpacingMeters_;
    diagnostics.startVolumeCubicMeters = referenceVolumeCubicMeters_
        + diagnostics.startCutCellVolumeCubicMeters;
    diagnostics.endVolumeCubicMeters = referenceVolumeCubicMeters_
        + diagnostics.endCutCellVolumeCubicMeters;
    diagnostics.geometryVolumeChangeCubicMeters =
        diagnostics.endVolumeCubicMeters
        - diagnostics.startVolumeCubicMeters;
    const double geometryVelocity =
        (step.endSurfaceOffsetMeters - step.startSurfaceOffsetMeters)
        / step.durationSeconds;

    std::size_t selectedFaceIndex = 0;
    double faceReactionPowerWatts = 0.0;
    const auto projectedNormalVelocity = velocities(
        projectedVelocityMetersPerSecond, axis_);
    for (const auto& face : interfaceDiagnostics.faces) {
        if (face.surfaceStableId != movingSurfaceStableId_) {
            continue;
        }
        if (selectedFaceIndex >= surfaceFaces_.size()) {
            throw std::invalid_argument(
                "interface diagnostics contain extra piston faces");
        }
        const auto& expected = surfaceFaces_[selectedFaceIndex++];
        if (face.minusRegionStableId != expected.minusRegionStableId
            || face.plusRegionStableId != expected.plusRegionStableId
            || face.axis != axis_
            || face.i != expected.i
            || face.j != expected.j
            || face.k != expected.k
            || !std::isfinite(face.areaSquareMeters)
            || !std::isfinite(face.normalVelocityMetersPerSecond)
            || !std::isfinite(face.pressurePowerWatts)
            || !std::isfinite(face.constraintReactionPowerWatts)) {
            throw std::invalid_argument(
                "interface diagnostics do not match the bound piston plane");
        }
        const double areaTolerance = combinedTolerance(
            settings.absoluteAreaToleranceSquareMeters,
            settings.relativeVolumeTolerance,
            face.areaSquareMeters, tileAreaSquareMeters_);
        if (std::abs(face.areaSquareMeters - tileAreaSquareMeters_)
            > areaTolerance) {
            throw std::invalid_argument(
                "piston face area changed incompatibly");
        }
        const std::size_t index = faceIndex(
            cellCounts_, face.i, face.j, face.k);
        diagnostics.maximumSurfaceVelocityErrorMetersPerSecond = std::max({
            diagnostics.maximumSurfaceVelocityErrorMetersPerSecond,
            std::abs(face.normalVelocityMetersPerSecond - geometryVelocity),
            std::abs(projectedNormalVelocity[index]
                     - face.normalVelocityMetersPerSecond),
        });
        diagnostics.surfaceSweptVolumeCubicMeters +=
            face.normalVelocityMetersPerSecond
            * face.areaSquareMeters * step.durationSeconds;
        faceReactionPowerWatts += face.constraintReactionPowerWatts;
    }
    if (selectedFaceIndex != surfaceFaces_.size()) {
        throw std::invalid_argument(
            "interface diagnostics omit piston faces");
    }

    for (const std::size_t index : openingFaceIndices_) {
        diagnostics.openingTransportVolumeCubicMeters +=
            projectedNormalVelocity[index] * tileAreaSquareMeters_
            * step.durationSeconds;
    }
    diagnostics.surfaceGeometryResidualCubicMeters =
        diagnostics.surfaceSweptVolumeCubicMeters
        - diagnostics.geometryVolumeChangeCubicMeters;
    diagnostics.continuityResidualCubicMeters =
        diagnostics.openingTransportVolumeCubicMeters
        - diagnostics.geometryVolumeChangeCubicMeters;

    const auto surface = std::lower_bound(
        interfaceDiagnostics.surfaces.begin(),
        interfaceDiagnostics.surfaces.end(), movingSurfaceStableId_,
        [](const MovingInterfaceSurfaceDiagnostics& candidate,
           const std::uint64_t stableId) {
            return candidate.stableId < stableId;
        });
    if (surface == interfaceDiagnostics.surfaces.end()
        || surface->stableId != movingSurfaceStableId_
        || surface->faceCount != surfaceFaces_.size()
        || !std::isfinite(surface->constraintReactionPowerWatts)) {
        throw std::invalid_argument(
            "piston surface aggregate is absent or invalid");
    }
    diagnostics.surfacePressurePowerWatts =
        surface->constraintReactionPowerWatts;
    diagnostics.rectangularSurfacePressureWorkJoules =
        surface->constraintReactionPowerWatts
        * step.durationSeconds;
    const double powerTolerance = combinedTolerance(
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
        std::abs(faceReactionPowerWatts),
        std::abs(surface->constraintReactionPowerWatts));
    if (std::abs(
            faceReactionPowerWatts
            - surface->constraintReactionPowerWatts)
        > powerTolerance) {
        throw std::invalid_argument(
            "piston face power does not match its surface aggregate");
    }

    diagnostics.finite =
        std::isfinite(diagnostics.crossSectionAreaSquareMeters)
        && std::isfinite(diagnostics.referenceVolumeCubicMeters)
        && std::isfinite(diagnostics.startVolumeCubicMeters)
        && std::isfinite(diagnostics.endVolumeCubicMeters)
        && std::isfinite(diagnostics.startCutCellVolumeCubicMeters)
        && std::isfinite(diagnostics.endCutCellVolumeCubicMeters)
        && std::isfinite(diagnostics.startCutCellVolumeFraction)
        && std::isfinite(diagnostics.endCutCellVolumeFraction)
        && std::isfinite(diagnostics.geometryVolumeChangeCubicMeters)
        && std::isfinite(diagnostics.surfaceSweptVolumeCubicMeters)
        && std::isfinite(diagnostics.openingTransportVolumeCubicMeters)
        && std::isfinite(diagnostics.surfaceGeometryResidualCubicMeters)
        && std::isfinite(diagnostics.continuityResidualCubicMeters)
        && std::isfinite(
            diagnostics.maximumSurfaceVelocityErrorMetersPerSecond)
        && std::isfinite(diagnostics.surfacePressurePowerWatts)
        && std::isfinite(
            diagnostics.rectangularSurfacePressureWorkJoules);
    const double volumeTolerance = combinedTolerance(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance,
        std::abs(diagnostics.geometryVolumeChangeCubicMeters),
        std::max(std::abs(diagnostics.surfaceSweptVolumeCubicMeters),
                 std::abs(diagnostics.openingTransportVolumeCubicMeters)));
    const double velocityTolerance = combinedTolerance(
        settings.absoluteVelocityToleranceMetersPerSecond,
        settings.relativeVelocityTolerance,
        std::abs(geometryVelocity),
        std::abs(geometryVelocity)
            + diagnostics.maximumSurfaceVelocityErrorMetersPerSecond);
    diagnostics.accepted = diagnostics.finite
        && diagnostics.maximumSurfaceVelocityErrorMetersPerSecond
            <= velocityTolerance
        && std::abs(diagnostics.surfaceGeometryResidualCubicMeters)
            <= volumeTolerance
        && std::abs(diagnostics.continuityResidualCubicMeters)
            <= volumeTolerance;
    return diagnostics;
}

PlanarControlVolumeRebaseResult rebasePlanarMovingControlVolume(
    const PeriodicCartesianGrid& grid,
    const PlanarMovingControlVolume& current,
    const FaceAlignedMovingInterface& rebasedInterfaces,
    const PlanarControlVolumeDiagnostics& terminalDiagnostics,
    const PlanarControlVolumeSettings& settings) {
    validateSettings(settings);
    if (!current.matches(grid)) {
        throw std::invalid_argument(
            "planar control-volume rebase grid does not match "
            "its current epoch");
    }
    if (terminalDiagnostics.version != planarMovingControlVolumeVersion
        || !terminalDiagnostics.finite
        || !terminalDiagnostics.accepted
        || terminalDiagnostics.movingSurfaceStableId
            != current.movingSurfaceStableId()
        || terminalDiagnostics.fluidRegionStableId
            != current.fluidRegionStableId()
        || terminalDiagnostics.axis != current.axis()
        || terminalDiagnostics.movingPlaneCoordinate
            != current.movingPlaneCoordinate()
        || terminalDiagnostics.openingPlaneCoordinate
            != current.openingPlaneCoordinate()) {
        throw std::invalid_argument(
            "planar control-volume rebase requires its accepted terminal ledger");
    }

    PlanarMovingControlVolume rebased(
        grid, rebasedInterfaces, current.movingSurfaceStableId(),
        current.openingPlaneCoordinate());
    const std::size_t expectedPlane =
        (current.movingPlaneCoordinate()
         + 1) % axisCount(grid.cellCounts(), current.axis());
    if (rebased.movingPlaneCoordinate() != expectedPlane
        || rebased.axis() != current.axis()
        || rebased.fluidRegionStableId() != current.fluidRegionStableId()
        || rebased.openingPlaneCoordinate()
            != current.openingPlaneCoordinate()
        || rebased.crossSectionAreaSquareMeters()
            != current.crossSectionAreaSquareMeters()
        || rebased.normalCellSpacingMeters()
            != current.normalCellSpacingMeters()) {
        throw std::invalid_argument(
            "planar control-volume rebase must advance exactly "
            "one compatible MAC plane");
    }

    PlanarControlVolumeRebaseDiagnostics diagnostics;
    diagnostics.movingSurfaceStableId = current.movingSurfaceStableId();
    diagnostics.fluidRegionStableId = current.fluidRegionStableId();
    diagnostics.axis = current.axis();
    diagnostics.previousMovingPlaneCoordinate =
        current.movingPlaneCoordinate();
    diagnostics.rebasedMovingPlaneCoordinate =
        rebased.movingPlaneCoordinate();
    diagnostics.openingPlaneCoordinate = current.openingPlaneCoordinate();
    diagnostics.completedCellOffsetMeters =
        terminalDiagnostics.endCutCellVolumeFraction
        * current.normalCellSpacingMeters();
    diagnostics.previousTerminalVolumeCubicMeters =
        terminalDiagnostics.endVolumeCubicMeters;
    diagnostics.rebasedReferenceVolumeCubicMeters =
        rebased.referenceVolumeCubicMeters();
    diagnostics.volumeContinuityResidualCubicMeters =
        diagnostics.rebasedReferenceVolumeCubicMeters
        - diagnostics.previousTerminalVolumeCubicMeters;

    const double expectedCompletedVolume =
        current.crossSectionAreaSquareMeters()
        * current.normalCellSpacingMeters();
    const double expectedTerminalVolume =
        current.referenceVolumeCubicMeters() + expectedCompletedVolume;
    const double volumeTolerance = combinedTolerance(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance,
        std::abs(expectedTerminalVolume),
        std::max(
            std::abs(terminalDiagnostics.endVolumeCubicMeters),
            std::abs(rebased.referenceVolumeCubicMeters())));
    const double areaTolerance = combinedTolerance(
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeVolumeTolerance,
        current.crossSectionAreaSquareMeters(),
        std::abs(terminalDiagnostics.crossSectionAreaSquareMeters));
    const double positionTolerance = combinedTolerance(
        settings.absolutePositionToleranceMeters,
        settings.relativeVolumeTolerance,
        current.normalCellSpacingMeters(),
        std::abs(diagnostics.completedCellOffsetMeters));
    diagnostics.finite =
        std::isfinite(diagnostics.completedCellOffsetMeters)
        && std::isfinite(diagnostics.previousTerminalVolumeCubicMeters)
        && std::isfinite(diagnostics.rebasedReferenceVolumeCubicMeters)
        && std::isfinite(diagnostics.volumeContinuityResidualCubicMeters);
    diagnostics.accepted = diagnostics.finite
        && std::abs(diagnostics.completedCellOffsetMeters
                    - current.normalCellSpacingMeters())
            <= positionTolerance
        && std::abs(terminalDiagnostics.referenceVolumeCubicMeters
                    - current.referenceVolumeCubicMeters())
            <= volumeTolerance
        && std::abs(terminalDiagnostics.crossSectionAreaSquareMeters
                    - current.crossSectionAreaSquareMeters())
            <= areaTolerance
        && std::abs(terminalDiagnostics.endCutCellVolumeCubicMeters
                    - expectedCompletedVolume)
            <= volumeTolerance
        && std::abs(terminalDiagnostics.endVolumeCubicMeters
                    - expectedTerminalVolume)
            <= volumeTolerance
        && std::abs(diagnostics.volumeContinuityResidualCubicMeters)
            <= volumeTolerance;
    return {std::move(rebased), diagnostics};
}

} // namespace simwing::fsi::fluid
