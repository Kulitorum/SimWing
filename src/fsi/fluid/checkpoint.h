#pragma once

#include "fluid/moving_interface.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t movingInterfaceFluidCheckpointVersion = 1;

struct MovingInterfaceFluidState {
    MacVelocityField velocityMetersPerSecond;
    CellScalarField pressurePascals;
    FaceAlignedMovingInterface interfaces;
    MovingInterfaceProjectionDiagnostics diagnostics;
};

// In-memory committed state for one fixed moving-interface topology epoch.
// The public metadata is suitable for orchestration and persistence adapters;
// the immutable private payload prevents a caller from changing fields without
// invalidating the checkpoint boundary. The topology fingerprint excludes
// prescribed velocity but includes the grid, stable identities, orientation,
// face coordinates, regions, and cell-to-region partition.
struct MovingInterfaceFluidCheckpoint {
    std::uint32_t version = movingInterfaceFluidCheckpointVersion;
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::uint64_t topologyFingerprint = 0;

private:
    friend MovingInterfaceFluidCheckpoint
    checkpointMovingInterfaceFluidState(
        const PeriodicCartesianGrid&,
        const MacVelocityField&,
        const CellScalarField&,
        const FaceAlignedMovingInterface&,
        const MovingInterfaceProjectionDiagnostics&);
    friend MovingInterfaceFluidState
    restoreMovingInterfaceFluidState(
        const PeriodicCartesianGrid&,
        const MovingInterfaceFluidCheckpoint&);

    struct Detail;
    std::shared_ptr<const Detail> detail;
};

// Captures accepted, finite velocity, pressure, interface, and diagnostic
// state. Unaccepted projection output is never checkpointable.
[[nodiscard]] MovingInterfaceFluidCheckpoint
checkpointMovingInterfaceFluidState(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocityMetersPerSecond,
    const CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionDiagnostics& diagnostics);

// Fully validates and copies the immutable payload before returning it. The
// caller may therefore assemble a larger transaction and commit the returned
// candidate only after its structural and topology-specific checks succeed.
[[nodiscard]] MovingInterfaceFluidState
restoreMovingInterfaceFluidState(
    const PeriodicCartesianGrid& grid,
    const MovingInterfaceFluidCheckpoint& checkpoint);

} // namespace simwing::fsi::fluid
