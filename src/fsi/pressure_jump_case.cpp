#include "pressure_jump_case.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

fluid::PeriodicCartesianGrid makeGrid() {
    return {{16, 4, 3}, {}, {4.0, 2.0, 1.5}};
}

std::vector<fluid::GridFacePressureJump> makePressureJumps(
    const fluid::PeriodicCartesianGrid& grid) {
    std::vector<fluid::GridFacePressureJump> result;
    const auto counts = grid.cellCounts();
    result.reserve(4 * counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                110, 1, 3, fluid::GridFaceAxis::X,
                4, j, k, 100.0, 0.2});
            result.push_back({
                120, 3, 2, fluid::GridFaceAxis::X,
                4, j, k, 150.0, 0.8});
            result.push_back({
                130, 2, 4, fluid::GridFaceAxis::X,
                12, j, k, -150.0, 0.2});
            result.push_back({
                140, 4, 1, fluid::GridFaceAxis::X,
                12, j, k, -100.0, 0.8});
        }
    }
    return result;
}

fluid::ProjectionSettings makeStepSettings() {
    fluid::ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 1.0 / 60.0;
    settings.absoluteResidualTolerance = 1.0e-10;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 2000;
    return settings;
}

} // namespace

PressureJumpCase::PressureJumpCase()
    : grid_(makeGrid()),
      velocity_(grid_),
      pressure_(grid_),
      pressureJumps_(grid_, makePressureJumps(grid_)),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader PressureJumpCase::traceHeader() const {
    return {pressureJumpCaseChecksum, pressureJumpCaseSolverId};
}

viewer::DiagnosticFrame PressureJumpCase::advance() {
    // This is a static oracle, so every published interval repeats the same
    // fresh projection instead of accumulating solver-tolerance drift.
    fluid::MacVelocityField candidateVelocity(grid_);
    fluid::CellScalarField candidatePressure(grid_);
    const fluid::ProjectionDiagnostics candidateDiagnostics =
        fluid::projectVelocityWithPressureJumps(
            grid_, candidateVelocity, candidatePressure,
            pressureJumps_, stepSettings_);
    if (!candidateDiagnostics.converged) {
        throw std::runtime_error(
            "pressure-jump case projection did not converge");
    }
    const std::uint64_t candidateStep = acceptedStepCount_ + 1;
    const double candidateTime = simulationTimeSeconds_
        + stepSettings_.timeStepSeconds;
    viewer::PressureJumpFrameContext context;
    context.sceneChecksum = pressureJumpCaseChecksum;
    context.solverCommit = pressureJumpCaseSolverId;
    context.step = candidateStep;
    context.simulationTimeSeconds = candidateTime;
    context.timeStepSeconds = stepSettings_.timeStepSeconds;
    context.densityKgPerCubicMeter =
        stepSettings_.densityKgPerCubicMeter;
    viewer::DiagnosticFrame frame = viewer::buildPressureJumpFrame(
        grid_, candidateVelocity, candidatePressure,
        pressureJumps_, candidateDiagnostics, context);
    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    diagnostics_ = candidateDiagnostics;
    acceptedStepCount_ = candidateStep;
    simulationTimeSeconds_ = candidateTime;
    return frame;
}

const fluid::PeriodicCartesianGrid&
PressureJumpCase::grid() const noexcept {
    return grid_;
}

const fluid::MacVelocityField&
PressureJumpCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField&
PressureJumpCase::pressure() const noexcept {
    return pressure_;
}

const fluid::SharpPressureJumpField&
PressureJumpCase::pressureJumps() const noexcept {
    return pressureJumps_;
}

const fluid::ProjectionSettings&
PressureJumpCase::stepSettings() const noexcept {
    return stepSettings_;
}

const fluid::ProjectionDiagnostics&
PressureJumpCase::diagnostics() const noexcept {
    return diagnostics_;
}

std::uint64_t PressureJumpCase::acceptedStepCount() const noexcept {
    return acceptedStepCount_;
}

double PressureJumpCase::simulationTimeSeconds() const noexcept {
    return simulationTimeSeconds_;
}

} // namespace simwing::fsi
