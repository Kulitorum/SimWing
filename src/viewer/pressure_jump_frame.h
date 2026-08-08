#pragma once

#include "fluid/interface_jump.h"
#include "fluid/projection.h"
#include "viewer_protocol.h"

#include <cstdint>
#include <string>

namespace simwing::viewer {

struct PressureJumpFrameContext {
    std::string sceneChecksum;
    std::string solverCommit;
    std::uint64_t step = 0;
    double simulationTimeSeconds = 0.0;
    double timeStepSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
};

// Builds an owning diagnostic frame for one accepted sharp-jump projection.
// Cell-centre vertices retain pressure and interpolated MAC velocity. Each
// authored crossing becomes a separate oriented quad at its normal segment
// fraction, including multiple layers on one grid face. Quad pressure samples
// reconstruct the smooth endpoint variation plus half of the local jump; they
// are diagnostic samples, not a cut-cell pressure interpolation claim.
[[nodiscard]] DiagnosticFrame buildPressureJumpFrame(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond,
    const fsi::fluid::CellScalarField& pressurePascals,
    const fsi::fluid::SharpPressureJumpField& pressureJumps,
    const fsi::fluid::ProjectionDiagnostics& diagnostics,
    const PressureJumpFrameContext& context);

} // namespace simwing::viewer
