#pragma once

#include "fluid/moving_interface.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarMovingControlVolumeVersion = 2;
inline constexpr std::uint32_t planarControlVolumeRebaseVersion = 1;

struct PlanarControlVolumeStep {
    // Positive offsets move the complete planar surface from its reference MAC
    // face into the plus-side cell. An ordinary step stops before the next grid
    // face; a terminal step reaches it exactly and must then use the explicit
    // rebase operation below.
    double startSurfaceOffsetMeters = 0.0;
    double endSurfaceOffsetMeters = 0.0;
    double durationSeconds = 0.0;
    // The ordinary supported interval retains a positive remainder in the
    // plus-side cell. This explicit flag permits an accepted step to end
    // exactly on the next MAC face, after which the topology must rebase
    // before another step is evaluated.
    bool endsAtCellBoundary = false;

    bool operator==(const PlanarControlVolumeStep&) const = default;
};

struct PlanarControlVolumeSettings {
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-11;
    double relativeVelocityTolerance = 1.0e-11;
    double absolutePositionToleranceMeters = 1.0e-12;
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
    std::size_t movingPlaneCoordinate = 0;
    std::size_t openingPlaneCoordinate = 0;
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
    // Complete fluid-on-surface constraint-reaction power, including direct
    // MAC velocity enforcement as well as adjacent pressure traction.
    double surfacePressurePowerWatts = 0.0;
    double rectangularSurfacePressureWorkJoules = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PlanarControlVolumeDiagnostics&) const = default;
};

// Transaction ledger for advancing one accepted planar surface epoch to the
// next positive-axis MAC plane. A rebase moves no geometry and advances no
// time: the completed old partial cell becomes one full reference layer in
// the new topology, so the two independently computed volumes must agree.
struct PlanarControlVolumeRebaseDiagnostics {
    std::uint32_t version = planarControlVolumeRebaseVersion;
    std::uint64_t movingSurfaceStableId = 0;
    std::uint64_t fluidRegionStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t previousMovingPlaneCoordinate = 0;
    std::size_t rebasedMovingPlaneCoordinate = 0;
    std::size_t openingPlaneCoordinate = 0;
    double completedCellOffsetMeters = 0.0;
    double previousTerminalVolumeCubicMeters = 0.0;
    double rebasedReferenceVolumeCubicMeters = 0.0;
    double volumeContinuityResidualCubicMeters = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const PlanarControlVolumeRebaseDiagnostics&) const = default;
};

// One-epoch, one-cell-motion verification operator. It owns no fluid
// state and performs no pressure solve. Instead it binds the moving-interface
// topology to one explicit open control volume and checks the geometric volume
// update independently against both interface sweep and projected opening
// transport. The separate operation below builds an exact next-plane epoch;
// general cut-cell metrics, nonplanar events, and folded surfaces remain future
// operators.
class PlanarMovingControlVolume final {
public:
    PlanarMovingControlVolume(
        const PeriodicCartesianGrid& grid,
        const FaceAlignedMovingInterface& interfaces,
        std::uint64_t movingSurfaceStableId,
        std::size_t openingPlaneCoordinate);

    [[nodiscard]] std::uint64_t movingSurfaceStableId() const noexcept;
    [[nodiscard]] std::uint64_t fluidRegionStableId() const noexcept;
    [[nodiscard]] bool matches(
        const PeriodicCartesianGrid& grid) const noexcept;
    [[nodiscard]] GridFaceAxis axis() const noexcept;
    [[nodiscard]] std::size_t movingPlaneCoordinate() const noexcept;
    [[nodiscard]] std::size_t openingPlaneCoordinate() const noexcept;
    [[nodiscard]] double crossSectionAreaSquareMeters() const noexcept;
    [[nodiscard]] double normalCellSpacingMeters() const noexcept;
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

struct PlanarControlVolumeRebaseResult {
    PlanarMovingControlVolume controlVolume;
    PlanarControlVolumeRebaseDiagnostics diagnostics;
};

// Builds and validates a candidate next-plane control volume without mutating
// the current epoch. The caller commits the returned candidate only after its
// own fluid/structure transaction also succeeds.
[[nodiscard]] PlanarControlVolumeRebaseResult
rebasePlanarMovingControlVolume(
    const PeriodicCartesianGrid& grid,
    const PlanarMovingControlVolume& current,
    const FaceAlignedMovingInterface& rebasedInterfaces,
    const PlanarControlVolumeDiagnostics& terminalDiagnostics,
    const PlanarControlVolumeSettings& settings = {});

} // namespace simwing::fsi::fluid
