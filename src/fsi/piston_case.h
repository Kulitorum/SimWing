#pragma once

#include "coupling.h"
#include "face_resolved_bridge.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char coupledPistonCaseChecksum[] =
    "sha256:simwing-uniform-coupled-piston-case-v1";
inline constexpr char coupledPistonCaseSolverId[] =
    "simwing-fsi-uniform-piston-worker-v1";

// Visible end-to-end verification harness for the current fixed-topology
// numerical foundations. A uniform face-aligned pressure solve is bridged by
// stable ID, integrated in time, accepted through XPBD, and then copied into
// an immutable viewer frame. It is intentionally not a general moving cut-cell
// or strongly iterated FSI algorithm.
class CoupledPistonCase final {
public:
    CoupledPistonCase();

    CoupledPistonCase(const CoupledPistonCase&) = delete;
    CoupledPistonCase& operator=(const CoupledPistonCase&) = delete;
    CoupledPistonCase(CoupledPistonCase&&) = delete;
    CoupledPistonCase& operator=(CoupledPistonCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;

private:
    Structure structure_;
    PlanarFaceResolvedFluidStructureBridge bridge_;
    ConservativeMacroStepCoupling coupling_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
};

} // namespace simwing::fsi
