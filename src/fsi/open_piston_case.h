#pragma once

#include "coupling.h"
#include "face_resolved_bridge.h"
#include "fluid/checkpoint.h"
#include "fluid/moving_control_volume.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

inline constexpr std::uint32_t openPistonConservationVersion = 1;
inline constexpr std::uint32_t openPistonCaseCheckpointVersion = 1;
inline constexpr std::uint64_t openPistonCaseDefinitionFingerprint =
    0x6f70656e70737436ull;

inline constexpr char openPistonCaseChecksum[] =
    "sha256:simwing-open-control-volume-piston-case-v6";
inline constexpr char openPistonCaseSolverId[] =
    "simwing-fsi-open-control-volume-piston-worker-v6";

struct OpenPistonConservationDiagnostics {
    std::uint32_t version = openPistonConservationVersion;
    StructureVector3 structureMomentumChangeNewtonSeconds;
    StructureVector3 fluidMomentumChangeNewtonSeconds;
    StructureVector3 pressureImpulseToStructureNewtonSeconds;
    StructureVector3 actuatorImpulseNewtonSeconds;
    StructureVector3 structureMomentumResidualNewtonSeconds;
    double structureMomentumResidualNormNewtonSeconds = 0.0;
    StructureVector3 fluidMomentumResidualNewtonSeconds;
    double fluidMomentumResidualNormNewtonSeconds = 0.0;
    StructureVector3 systemMomentumResidualNewtonSeconds;
    double systemMomentumResidualNormNewtonSeconds = 0.0;
    double structureKineticEnergyChangeJoules = 0.0;
    double fluidKineticEnergyChangeJoules = 0.0;
    double pressureWorkToStructureJoules = 0.0;
    double actuatorWorkJoules = 0.0;
    double structureEnergyResidualJoules = 0.0;
    double fluidEnergyResidualJoules = 0.0;
    double systemEnergyResidualJoules = 0.0;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const OpenPistonConservationDiagnostics&) const = default;
};

// Composite accepted macro-step checkpoint. The public epoch metadata lets a
// worker persist and route checkpoints without exposing mutable numerical
// payloads. The private payload composes the complete Structure checkpoint
// with a topology-bound moving-interface fluid checkpoint and every committed
// open-piston diagnostic/accounting value.
struct OpenPistonCaseCheckpoint {
    // Publicly nameable only for the separate persistence translation unit.
    // The definition remains private to src/fsi.
    struct Detail;

    std::uint32_t version = openPistonCaseCheckpointVersion;
    std::uint64_t caseDefinitionFingerprint =
        openPistonCaseDefinitionFingerprint;
    std::uint64_t acceptedStepCount = 0;
    std::uint64_t topologyRebaseCount = 0;
    std::size_t movingPlaneCoordinate = 0;
    double surfaceOffsetMeters = 0.0;

private:
    friend class OpenPistonCase;
    friend struct OpenPistonCheckpointPersistenceAccess;
    std::shared_ptr<const Detail> detail;
};

// Visible verification harness for an accelerating planar piston in one
// connected periodic fluid region. A complete moving sheet is nonseparating:
// fluid is projected around the remaining grid path and crosses an explicit
// open control-volume section. An actuator supplies the prescribed initial
// structural impulse while CFD supplies the resisting pressure load. The
// piston then coasts, exercising partial-cell volume and opening-flux GCL
// ledgers on every accepted step. At an exact MAC-face crossing it validates
// volume continuity, remaps the constraint by one face without a material
// velocity jump, and commits the new topology only with the complete frame.
// The pressure transfer retains face-resolved material patches while requiring
// the structural plate to match the unwrapped physical cut-surface plane. The
// projection's complete macro-step-average constraint reaction is sampled at
// both endpoint kinematics; fluid, structure, actuator, and system momentum and
// kinetic-energy ledgers must all close before commit.
class OpenPistonCase final {
public:
    OpenPistonCase();

    OpenPistonCase(const OpenPistonCase&) = delete;
    OpenPistonCase& operator=(const OpenPistonCase&) = delete;
    OpenPistonCase(OpenPistonCase&&) = delete;
    OpenPistonCase& operator=(OpenPistonCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] OpenPistonCaseCheckpoint checkpoint() const;
    void restore(const OpenPistonCaseCheckpoint& checkpoint);

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] const fluid::PlanarControlVolumeDiagnostics&
    controlVolumeDiagnostics() const noexcept;
    [[nodiscard]] const fluid::PlanarControlVolumeRebaseDiagnostics&
    lastRebaseDiagnostics() const noexcept;
    [[nodiscard]] const PlanarFaceResolvedBridgeDiagnostics&
    bridgeDiagnostics() const noexcept;
    [[nodiscard]] const fluid::PlanarCutSurfacePressureDiagnostics&
    cutSurfaceDiagnostics() const noexcept;
    [[nodiscard]] const OpenPistonConservationDiagnostics&
    conservationDiagnostics() const noexcept;
    [[nodiscard]] double surfaceOffsetMeters() const noexcept;
    [[nodiscard]] std::size_t movingPlaneCoordinate() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] double lastRebaseVelocityResidualMetersPerSecond()
        const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::FaceAlignedMovingInterface interfaces_;
    fluid::MacVelocityField fluidVelocity_;
    fluid::CellScalarField fluidPressure_;
    fluid::MovingInterfaceProjectionDiagnostics fluidDiagnostics_;
    Structure structure_;
    PlanarFaceResolvedFluidStructureBridge bridge_;
    ConservativeMacroStepCoupling coupling_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    fluid::PlanarMovingControlVolume controlVolume_;
    fluid::PlanarControlVolumeDiagnostics controlVolumeDiagnostics_;
    fluid::PlanarControlVolumeRebaseDiagnostics lastRebaseDiagnostics_;
    fluid::PlanarCutSurfacePressureDiagnostics cutSurfaceDiagnostics_;
    PlanarFaceResolvedBridgeDiagnostics bridgeDiagnostics_;
    OpenPistonConservationDiagnostics conservationDiagnostics_;
    double surfaceOffsetMeters_ = 0.0;
    double lastRebaseVelocityResidualMetersPerSecond_ = 0.0;
    std::uint64_t topologyRebaseCount_ = 0;
};

} // namespace simwing::fsi
