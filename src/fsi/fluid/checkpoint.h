#pragma once

#include "fluid/porous_interface.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t movingInterfaceFluidCheckpointVersion = 1;
inline constexpr std::uint32_t movingPorousFluidCheckpointVersion = 1;
inline constexpr std::uint16_t
    movingInterfaceFluidCheckpointProtocolVersion = 1;

struct MovingInterfaceFluidCheckpointLimits {
    std::uint64_t maximumBytes = 64ULL * 1024ULL * 1024ULL;
    std::uint64_t maximumScalarSamples = 5'000'000;
    std::uint64_t maximumInterfaceFaces = 1'000'000;
    std::uint64_t maximumFluidRegions = 1'000'000;
    std::uint64_t maximumDiagnosticSurfaces = 1'000'000;
};

enum class MovingInterfaceFluidCheckpointErrorCode {
    None,
    InvalidData,
    InvalidMagic,
    UnsupportedVersion,
    LimitExceeded,
    Truncated,
    TrailingData,
    ChecksumMismatch,
};

struct MovingInterfaceFluidCheckpointError {
    MovingInterfaceFluidCheckpointErrorCode code =
        MovingInterfaceFluidCheckpointErrorCode::None;
    std::string message;

    [[nodiscard]] explicit operator bool() const noexcept {
        return code != MovingInterfaceFluidCheckpointErrorCode::None;
    }
};

struct MovingInterfaceFluidState {
    MacVelocityField velocityMetersPerSecond;
    CellScalarField pressurePascals;
    FaceAlignedMovingInterface interfaces;
    MovingInterfaceProjectionDiagnostics diagnostics;
};

struct MovingPorousFluidState {
    // Predicted field supplied to the nonlinear projection. It is retained as
    // diagnostic provenance so midpoint constitutive samples and prescribed
    // jump power can be reconstructed exactly after persistence.
    MacVelocityField predictedVelocityMetersPerSecond;
    MacVelocityField velocityMetersPerSecond;
    CellScalarField pressurePascals;
    FaceAlignedMovingInterface interfaces;
    std::vector<PorousGridFaceCrossing> porousCrossings;
    SharpPressureJumpField prescribedPressureJumps;
    MovingPorousProjectionDiagnostics diagnostics;
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

// In-memory committed state for one coupled moving/porous topology epoch.
// Porous definitions are stored in canonical crossing order. The topology
// fingerprint includes moving, porous, and prescribed-jump geometry and stable
// identities, while prescribed velocities, resistance coefficients, and jump
// values remain mutable state within that fixed topology.
struct MovingPorousFluidCheckpoint {
    std::uint32_t version = movingPorousFluidCheckpointVersion;
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::uint64_t topologyFingerprint = 0;

private:
    friend MovingPorousFluidCheckpoint checkpointMovingPorousFluidState(
        const PeriodicCartesianGrid&,
        const MacVelocityField&,
        const MacVelocityField&,
        const CellScalarField&,
        const FaceAlignedMovingInterface&,
        const std::vector<PorousGridFaceCrossing>&,
        const SharpPressureJumpField&,
        const MovingPorousProjectionDiagnostics&);
    friend MovingPorousFluidState restoreMovingPorousFluidState(
        const PeriodicCartesianGrid&,
        const MovingPorousFluidCheckpoint&);

    struct Detail;
    std::shared_ptr<const Detail> detail;
};

// Deterministic bounded little-endian persistence for the complete accepted
// state. Decode reconstructs the grid/interface objects and commits the output
// only after checkpointMovingInterfaceFluidState accepts all fields and nested
// diagnostics, including the recomputed topology fingerprint.
[[nodiscard]] bool serializeMovingInterfaceFluidCheckpoint(
    const MovingInterfaceFluidCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    MovingInterfaceFluidCheckpointError* error = nullptr,
    const MovingInterfaceFluidCheckpointLimits& limits = {});

[[nodiscard]] bool deserializeMovingInterfaceFluidCheckpoint(
    std::span<const std::uint8_t> bytes,
    MovingInterfaceFluidCheckpoint& checkpoint,
    MovingInterfaceFluidCheckpointError* error = nullptr,
    const MovingInterfaceFluidCheckpointLimits& limits = {});

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

// Captures only a completely accepted moving/porous solve. Validation binds
// the constitutive samples to the calibrated resistance law at their declared
// endpoint or midpoint time, and binds the nested projection diagnostics to
// the committed moving-interface state.
[[nodiscard]] MovingPorousFluidCheckpoint checkpointMovingPorousFluidState(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& predictedVelocityMetersPerSecond,
    const MacVelocityField& velocityMetersPerSecond,
    const CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const std::vector<PorousGridFaceCrossing>& porousCrossings,
    const SharpPressureJumpField& prescribedPressureJumps,
    const MovingPorousProjectionDiagnostics& diagnostics);

[[nodiscard]] MovingPorousFluidState restoreMovingPorousFluidState(
    const PeriodicCartesianGrid& grid,
    const MovingPorousFluidCheckpoint& checkpoint);

} // namespace simwing::fsi::fluid
