#pragma once

#include "coupling.h"
#include "face_resolved_bridge.h"
#include "fluid/moving_control_volume.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>

namespace simwing::fsi {

inline constexpr char openPistonCaseChecksum[] =
    "sha256:simwing-open-control-volume-piston-case-v3";
inline constexpr char openPistonCaseSolverId[] =
    "simwing-fsi-open-control-volume-piston-worker-v3";

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
// the structural plate to match the unwrapped physical cut-surface plane.
class OpenPistonCase final {
public:
    OpenPistonCase();

    OpenPistonCase(const OpenPistonCase&) = delete;
    OpenPistonCase& operator=(const OpenPistonCase&) = delete;
    OpenPistonCase(OpenPistonCase&&) = delete;
    OpenPistonCase& operator=(OpenPistonCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] const fluid::PlanarControlVolumeDiagnostics&
    controlVolumeDiagnostics() const noexcept;
    [[nodiscard]] const fluid::PlanarControlVolumeRebaseDiagnostics&
    lastRebaseDiagnostics() const noexcept;
    [[nodiscard]] const PlanarFaceResolvedBridgeDiagnostics&
    bridgeDiagnostics() const noexcept;
    [[nodiscard]] double surfaceOffsetMeters() const noexcept;
    [[nodiscard]] std::size_t movingPlaneCoordinate() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] double lastRebaseVelocityResidualMetersPerSecond()
        const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
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
    PlanarFaceResolvedBridgeDiagnostics bridgeDiagnostics_;
    double surfaceOffsetMeters_ = 0.0;
    double lastRebaseVelocityResidualMetersPerSecond_ = 0.0;
    std::uint64_t topologyRebaseCount_ = 0;
};

} // namespace simwing::fsi
