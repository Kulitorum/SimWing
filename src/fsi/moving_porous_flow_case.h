#pragma once

#include "fluid/evolution.h"
#include "pressure_jump_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char movingPorousFlowCaseChecksum[] =
    "sha256:simwing-moving-planar-porous-flow-case-v1";
inline constexpr char movingPorousFlowCaseSolverId[] =
    "simwing-fsi-moving-planar-porous-flow-worker-v1";

// Visible Qt-free oracle for one prescribed translating porous plane inside a
// complete periodic flow step. Each advance samples the sheet at both porous
// half-step midpoints, admits retained or adjacent topology, and publishes the
// final wrapped crossing plane plus its unwrapped topology diagnostics. This
// is prescribed planar source motion, not fluid/structure coupling.
class MovingPorousFlowCase final {
public:
    MovingPorousFlowCase();

    MovingPorousFlowCase(const MovingPorousFlowCase&) = delete;
    MovingPorousFlowCase& operator=(const MovingPorousFlowCase&) = delete;
    MovingPorousFlowCase(MovingPorousFlowCase&&) = delete;
    MovingPorousFlowCase& operator=(MovingPorousFlowCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::SharpPressureJumpField&
    pressureJumps() const noexcept;
    [[nodiscard]] const fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] const fluid::MovingPorousFaceTopology&
    porousTopology() const noexcept;
    [[nodiscard]] double sheetPositionMeters() const noexcept;
    [[nodiscard]] double sheetVelocityMetersPerSecond() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::SharpPressureJumpField pressureJumps_;
    fluid::PorousIterationSettings porousIteration_;
    fluid::PeriodicFlowStrangSspRk2Settings flowSettings_;
    fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics diagnostics_;
    fluid::MovingPorousFaceTopology porousTopology_;
    double sheetPositionMeters_ = 3.48;
    double sheetVelocityMetersPerSecond_ = 0.4;
    std::uint64_t topologyRebaseCount_ = 0;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
