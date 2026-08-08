#include "periodic_flow_case.h"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace simwing::fsi {
namespace {

fluid::PeriodicCartesianGrid makeGrid() {
    const double twoPi = 2.0 * std::numbers::pi;
    return {{18, 18, 2}, {}, {twoPi, twoPi, 1.0}};
}

fluid::MacVelocityField makeInitialVelocity(
    const fluid::PeriodicCartesianGrid& grid) {
    constexpr double backgroundX = 0.35;
    constexpr double backgroundY = -0.2;
    fluid::MacVelocityField result(grid);
    const fluid::GridCellCounts counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const fluid::Vector3 xFace =
                    grid.xFaceCenterMeters(i, j, k);
                const fluid::Vector3 yFace =
                    grid.yFaceCenterMeters(i, j, k);
                result.xFaces()[index] = backgroundX
                    + std::sin(xFace.x) * std::cos(xFace.y);
                result.yFaces()[index] = backgroundY
                    - std::cos(yFace.x) * std::sin(yFace.y);
            }
        }
    }
    return result;
}

fluid::PeriodicFlowStrangSubcyclingSettings makeStepSettings() {
    fluid::PeriodicFlowStrangSubcyclingSettings result;
    result.flow.densityKgPerCubicMeter = 1.225;
    result.flow.kinematicViscositySquareMetersPerSecond = 0.02;
    result.flow.timeStepSeconds = 1.0 / 60.0;
    result.flow.advectionReconstruction =
        fluid::VariableMacReconstruction::MonotonizedCentral;
    result.flow.advectionAbsoluteDivergenceTolerancePerSecond = 2.0e-10;
    result.flow.advectionRelativeDivergenceTolerance = 1.0e-12;
    result.flow.projectionAbsoluteResidualTolerance = 1.0e-10;
    result.flow.projectionRelativeResidualTolerance = 1.0e-12;
    result.flow.projectionMaximumIterations = 2000;
    result.flow.absoluteMomentumToleranceNewtonSeconds = 1.0e-10;
    result.flow.relativeMomentumTolerance = 1.0e-12;
    result.flow.absoluteEnergyToleranceJoules = 1.0e-10;
    result.flow.relativeEnergyTolerance = 1.0e-12;
    result.maximumSubsteps = 64;
    return result;
}

} // namespace

PeriodicFlowCase::PeriodicFlowCase()
    : grid_(makeGrid()),
      velocity_(makeInitialVelocity(grid_)),
      pressure_(grid_),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader PeriodicFlowCase::traceHeader() const {
    return {periodicFlowCaseChecksum, periodicFlowCaseSolverId};
}

viewer::DiagnosticFrame PeriodicFlowCase::advance() {
    diagnostics_ = fluid::advancePeriodicFlowStrangSspRk2Subcycled(
        grid_, velocity_, pressure_, stepSettings_);
    if (!diagnostics_.accepted) {
        throw std::runtime_error(
            "periodic flow case rejected its requested outer interval");
    }
    ++acceptedStepCount_;
    simulationTimeSeconds_ += stepSettings_.flow.timeStepSeconds;
    viewer::PeriodicFluidFrameContext context;
    context.sceneChecksum = periodicFlowCaseChecksum;
    context.solverCommit = periodicFlowCaseSolverId;
    context.step = acceptedStepCount_;
    context.simulationTimeSeconds = simulationTimeSeconds_;
    context.densityKgPerCubicMeter =
        stepSettings_.flow.densityKgPerCubicMeter;
    return viewer::buildPeriodicFluidFrame(
        grid_, velocity_, pressure_, diagnostics_, context);
}

const fluid::PeriodicCartesianGrid& PeriodicFlowCase::grid() const noexcept {
    return grid_;
}

const fluid::MacVelocityField& PeriodicFlowCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField& PeriodicFlowCase::pressure() const noexcept {
    return pressure_;
}

const fluid::PeriodicFlowStrangSubcyclingSettings&
PeriodicFlowCase::stepSettings() const noexcept {
    return stepSettings_;
}

const fluid::PeriodicFlowStrangSubcyclingDiagnostics&
PeriodicFlowCase::diagnostics() const noexcept {
    return diagnostics_;
}

std::uint64_t PeriodicFlowCase::acceptedStepCount() const noexcept {
    return acceptedStepCount_;
}

double PeriodicFlowCase::simulationTimeSeconds() const noexcept {
    return simulationTimeSeconds_;
}

} // namespace simwing::fsi
