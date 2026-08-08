#pragma once

#include "coupled_state.h"
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
inline constexpr char strongCoupledPistonCaseSolverId[] =
    "simwing-fsi-strong-light-piston-v1";

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

struct StrongCoupledPistonStepDiagnostics {
    StrongCouplingMacroStepRunResult coupling;
    double startSpeedMetersPerSecond = 0.0;
    double acceptedSpeedMetersPerSecond = 0.0;
    double acceptedInterfaceSpeedMetersPerSecond = 0.0;
    double velocityClosureMetersPerSecond = 0.0;
    bool finite = true;

    bool operator==(
        const StrongCoupledPistonStepDiagnostics&) const = default;
};

// Added-mass canonical for the generic strong-coupling loop. The light rigid
// XPBD plate is driven by the actual moving-interface pressure projection.
// Each fixed-point solve starts from the same accepted Structure/fluid epoch;
// only a converged macro-step becomes the next persistent state.
class StrongCoupledPistonCase final {
public:
    StrongCoupledPistonCase();

    StrongCoupledPistonCase(const StrongCoupledPistonCase&) = delete;
    StrongCoupledPistonCase& operator=(
        const StrongCoupledPistonCase&) = delete;
    StrongCoupledPistonCase(StrongCoupledPistonCase&&) = delete;
    StrongCoupledPistonCase& operator=(
        StrongCoupledPistonCase&&) = delete;

    [[nodiscard]] StrongCoupledPistonStepDiagnostics advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const fluid::MovingInterfaceFluidState&
    fluidState() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;

private:
    Structure structure_;
    fluid::PeriodicCartesianGrid grid_;
    fluid::MovingInterfaceFluidState fluidState_;
    PlanarFaceResolvedFluidStructureBridge bridge_;
    ConservativeMacroStepCoupling coupling_;
    StructureStepSettings stepSettings_;
    AitkenRelaxationSettings relaxationSettings_;
    CouplingConvergenceSettings convergenceSettings_;
    CouplingMacroStepRetrySettings retrySettings_;
};

} // namespace simwing::fsi
