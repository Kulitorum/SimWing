#pragma once

#include "fluid/interface_jump.h"
#include "fluid/projection.h"
#include "pressure_jump_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char pressureJumpCaseChecksum[] =
    "sha256:simwing-split-static-pressure-jump-case-v1";
inline constexpr char pressureJumpCaseSolverId[] =
    "simwing-fsi-pressure-jump-worker-v1";

// Visible static sharp-interface verification harness. Each slab boundary is
// represented by two ordered subcell region transitions on every crossed grid
// face. Projection must preserve the exact piecewise-constant pressure and zero
// velocity while immutable frames expose every separate interface layer. This
// is not moving folded topology or a cut-cell evolution case.
class PressureJumpCase final {
public:
    PressureJumpCase();

    PressureJumpCase(const PressureJumpCase&) = delete;
    PressureJumpCase& operator=(const PressureJumpCase&) = delete;
    PressureJumpCase(PressureJumpCase&&) = delete;
    PressureJumpCase& operator=(PressureJumpCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::SharpPressureJumpField&
    pressureJumps() const noexcept;
    [[nodiscard]] const fluid::ProjectionSettings&
    stepSettings() const noexcept;
    [[nodiscard]] const fluid::ProjectionDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::SharpPressureJumpField pressureJumps_;
    fluid::ProjectionSettings stepSettings_;
    fluid::ProjectionDiagnostics diagnostics_;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
