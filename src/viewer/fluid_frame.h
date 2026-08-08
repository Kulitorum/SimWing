#pragma once

#include "fluid/evolution.h"
#include "viewer_protocol.h"

#include <cstdint>
#include <string>

namespace simwing::viewer {

struct PeriodicFluidFrameContext {
    std::string sceneChecksum;
    std::string solverCommit;
    std::uint64_t step = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
};

// Copies one committed periodic-fluid state into an owning diagnostic frame.
// Each cell becomes one stable, unconnected DiagnosticVertex at its exact cell
// centre. Pressure is copied directly; the published velocity is the standard
// arithmetic average of the two bounding MAC faces for that cell. No frame
// storage aliases the solver fields. Only an accepted complete subcycled outer
// interval may be published.
[[nodiscard]] DiagnosticFrame buildPeriodicFluidFrame(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond,
    const fsi::fluid::CellScalarField& pressurePascals,
    const fsi::fluid::PeriodicFlowStrangSubcyclingDiagnostics& diagnostics,
    const PeriodicFluidFrameContext& context);

} // namespace simwing::viewer
