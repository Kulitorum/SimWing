#pragma once

#include "fluid/evolution.h"
#include "fluid/porous_flow.h"
#include "fluid/projection.h"
#include "pressure_jump_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char porousFlowCaseChecksum[] =
    "sha256:simwing-pressure-driven-porous-plug-case-v1";
inline constexpr char porousFlowCaseSolverId[] =
    "simwing-fsi-porous-flow-worker-v1";

// Visible pressure-driven porous-flow verification harness. A midpoint plug
// momentum solve determines the uniform periodic flux through a calibrated
// porous plane. A gauge plane cancels the endpoint jump so the periodic
// pressure field is single-valued; the physical driving rise remains separate
// in the plug impulse/work ledger. The accepted endpoint jump is then verified
// by grid projection and published as immutable cell/layer diagnostics. This
// is not general implicit porous coupling on a nonuniform or moving interface.
class PorousFlowCase final {
public:
    PorousFlowCase();

    PorousFlowCase(const PorousFlowCase&) = delete;
    PorousFlowCase& operator=(const PorousFlowCase&) = delete;
    PorousFlowCase(PorousFlowCase&&) = delete;
    PorousFlowCase& operator=(PorousFlowCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::SharpPressureJumpField&
    pressureJumps() const noexcept;
    [[nodiscard]] const fluid::ProjectionSettings&
    stepSettings() const noexcept;
    [[nodiscard]] const fluid::PorousPlugFlowSettings&
    plugSettings() const noexcept;
    [[nodiscard]] const fluid::ProjectionDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSspRk2Diagnostics&
    flowDiagnostics() const noexcept;
    [[nodiscard]] const fluid::PorousPlugFlowDiagnostics&
    plugDiagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::SharpPressureJumpField pressureJumps_;
    fluid::ProjectionSettings stepSettings_;
    fluid::PeriodicFlowStrangSspRk2Settings flowSettings_;
    fluid::PorousPlugFlowSettings plugSettings_;
    fluid::ProjectionDiagnostics diagnostics_;
    fluid::PeriodicFlowStrangSspRk2Diagnostics flowDiagnostics_;
    fluid::PorousPlugFlowDiagnostics plugDiagnostics_;
    double flowVelocityMetersPerSecond_ = 0.0;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
