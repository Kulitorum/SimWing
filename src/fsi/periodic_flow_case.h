#pragma once

#include "fluid/evolution.h"
#include "fluid_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char periodicFlowCaseChecksum[] =
    "sha256:simwing-periodic-taylor-green-flow-case-v1";
inline constexpr char periodicFlowCaseSolverId[] =
    "simwing-fsi-periodic-flow-worker-v1";

// A visible, deterministic CFD verification case. It advances a smooth
// Galilean-shifted Taylor-Green field through the bounded periodic
// Strang/SSPRK2 subcycler and publishes only accepted cell-centred snapshots.
// It is deliberately not a canopy, cut-cell, or aerodynamic-truth case.
class PeriodicFlowCase final {
public:
    PeriodicFlowCase();

    PeriodicFlowCase(const PeriodicFlowCase&) = delete;
    PeriodicFlowCase& operator=(const PeriodicFlowCase&) = delete;
    PeriodicFlowCase(PeriodicFlowCase&&) = delete;
    PeriodicFlowCase& operator=(PeriodicFlowCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSubcyclingSettings&
    stepSettings() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSubcyclingDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::PeriodicFlowStrangSubcyclingSettings stepSettings_;
    fluid::PeriodicFlowStrangSubcyclingDiagnostics diagnostics_;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
