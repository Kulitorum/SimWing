#pragma once

#include "fluid/evolution.h"
#include "viewer_protocol.h"

#include <cstdint>
#include <string>
#include <vector>

namespace simwing::viewer {

struct PeriodicFluidFrameContext {
    std::string sceneChecksum;
    std::string solverCommit;
    std::uint64_t step = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
};

struct PeriodicFluidCellFields {
    std::vector<Vec3d> velocityMetersPerSecond;
    std::vector<double> speedMetersPerSecond;
    std::vector<double> divergencePerSecond;
    std::vector<Vec3d> vorticityPerSecond;
    std::vector<double> vorticityMagnitudePerSecond;
    double maximumAbsoluteDivergencePerSecond = 0.0;
    double maximumVorticityPerSecond = 0.0;
};

// Samples owning cell-centred diagnostics from one finite periodic MAC field.
// Divergence is the exact finite-volume operator used by projection. Velocity
// averages the two bounding faces per component; vorticity is the centred
// periodic curl of that published velocity.
[[nodiscard]] PeriodicFluidCellFields buildPeriodicFluidCellFields(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond);

// Copies one committed periodic-fluid state into an owning diagnostic frame.
// Each cell becomes one stable, unconnected DiagnosticVertex at its exact cell
// centre. Pressure is copied directly; the published velocity is the standard
// arithmetic average of the two bounding MAC faces for that cell. Divergence
// uses the solver's finite-volume MAC operator. Vorticity is a diagnostic curl
// of that published cell-centred velocity. No frame storage aliases the solver
// fields. Only an accepted complete subcycled outer interval may be published.
[[nodiscard]] DiagnosticFrame buildPeriodicFluidFrame(
    const fsi::fluid::PeriodicCartesianGrid& grid,
    const fsi::fluid::MacVelocityField& velocityMetersPerSecond,
    const fsi::fluid::CellScalarField& pressurePascals,
    const fsi::fluid::PeriodicFlowStrangSubcyclingDiagnostics& diagnostics,
    const PeriodicFluidFrameContext& context);

} // namespace simwing::viewer
