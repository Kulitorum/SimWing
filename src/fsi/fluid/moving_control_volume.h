#pragma once

#include "fluid/moving_interface.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarMovingControlVolumeVersion = 1;

struct PlanarControlVolumeStep {
    // Positive offsets move the complete planar surface from its reference MAC
    // face into the plus-side cell. This first slice stops before the next grid
    // face, where a future topology-rebase operation is required.
    double startSurfaceOffsetMeters = 0.0;
    double endSurfaceOffsetMeters = 0.0;
    double durationSeconds = 0.0;

    bool operator==(const PlanarControlVolumeStep&) const = default;
};

struct PlanarControlVolumeSettings {
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-11;
    double relativeVelocityTolerance = 1.0e-11;
    double absoluteAreaToleranceSquareMeters = 1.0e-12;
    double absoluteVolumeToleranceCubicMeters = 1.0e-11;
    double relativeVolumeTolerance = 1.0e-11;
    double absolutePowerToleranceWatts = 1.0e-10;
    double relativePowerTolerance = 1.0e-11;
    double minimumRemainingCellLengthMeters = 1.0e-12;
};

// Independent geometry and transport ledgers for one open, axis-aligned
// piston chamber. The opening is a complete unconstrained MAC plane and the
// moving surface is a complete nonseparating plane whose two sides retain one
// connected fluid-region ID. Positive opening transport enters the chamber.
struct PlanarControlVolumeDiagnostics {
    std::uint32_t version = planarMovingControlVolumeVersion;
    std::uint64_t movingSurfaceStableId = 0;
    std::uint64_t fluidRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t movingSurfaceFaceCount = 0;
    std::size_t openingFaceCount = 0;
    double crossSectionAreaSquareMeters = 0.0;
    double referenceVolumeCubicMeters = 0.0;
    double startVolumeCubicMeters = 0.0;
    double endVolumeCubicMeters = 0.0;
    double startCutCellVolumeCubicMeters = 0.0;
    double endCutCellVolumeCubicMeters = 0.0;
    double startCutCellVolumeFraction = 0.0;
    double endCutCellVolumeFraction = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double surfaceSweptVolumeCubicMeters = 0.0;
    double openingTransportVolumeCubicMeters = 0.0;
    double surfaceGeometryResidualCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double maximumSurfaceVelocityErrorMetersPerSecond = 0.0;
    double surfacePressurePowerWatts = 0.0;
    double rectangularSurfacePressureWorkJoules = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PlanarControlVolumeDiagnostics&) const = default;
};

// Fixed-topology, one-cell-motion verification operator. It owns no fluid
// state and performs no pressure solve. Instead it binds the moving-interface
// topology to one explicit open control volume and checks the geometric volume
// update independently against both interface sweep and projected opening
// transport. General cut-cell metrics, topology rebasing, and folded surfaces
// remain future operators.
class PlanarMovingControlVolume final {
public:
    PlanarMovingControlVolume(
        const PeriodicCartesianGrid& grid,
        const FaceAlignedMovingInterface& interfaces,
        std::uint64_t movingSurfaceStableId,
        std::size_t openingPlaneCoordinate);

    [[nodiscard]] std::uint64_t movingSurfaceStableId() const noexcept;
    [[nodiscard]] std::uint64_t fluidRegionStableId() const noexcept;
    [[nodiscard]] GridFaceAxis axis() const noexcept;
    [[nodiscard]] std::size_t movingPlaneCoordinate() const noexcept;
    [[nodiscard]] std::size_t openingPlaneCoordinate() const noexcept;
    [[nodiscard]] double crossSectionAreaSquareMeters() const noexcept;
    [[nodiscard]] double referenceVolumeCubicMeters() const noexcept;

    [[nodiscard]] PlanarControlVolumeDiagnostics evaluate(
        const PeriodicCartesianGrid& grid,
        const MacVelocityField& projectedVelocityMetersPerSecond,
        const MovingInterfaceProjectionDiagnostics& interfaceDiagnostics,
        const PlanarControlVolumeStep& step,
        const PlanarControlVolumeSettings& settings = {}) const;

private:
    struct SurfaceFace {
        std::uint64_t minusRegionStableId = 0;
        std::uint64_t plusRegionStableId = 0;
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t k = 0;
    };

    GridCellCounts cellCounts_;
    Vector3 lowerMeters_;
    Vector3 upperMeters_;
    std::uint64_t movingSurfaceStableId_ = 0;
    std::uint64_t fluidRegionStableId_ = 0;
    GridFaceAxis axis_ = GridFaceAxis::X;
    std::size_t movingPlaneCoordinate_ = 0;
    std::size_t openingPlaneCoordinate_ = 0;
    double normalCellSpacingMeters_ = 0.0;
    double tileAreaSquareMeters_ = 0.0;
    double crossSectionAreaSquareMeters_ = 0.0;
    double referenceVolumeCubicMeters_ = 0.0;
    std::vector<SurfaceFace> surfaceFaces_;
    std::vector<std::size_t> openingFaceIndices_;
};

} // namespace simwing::fsi::fluid
