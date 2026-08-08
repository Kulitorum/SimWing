#include "periodic_flow_case.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

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

bool gridMetadataMatches(const fluid::PeriodicCartesianGrid& grid,
                         const fluid::GridCellCounts cellCounts,
                         const fluid::Vector3 lowerMeters,
                         const fluid::Vector3 upperMeters) noexcept {
    return grid.cellCounts() == cellCounts
        && grid.lowerMeters() == lowerMeters
        && grid.upperMeters() == upperMeters;
}

bool validCommittedState(
    const std::uint64_t acceptedStepCount,
    const double simulationTimeSeconds,
    const fluid::PeriodicFlowStrangSubcyclingSettings& settings,
    const fluid::PeriodicFlowStrangSubcyclingDiagnostics& diagnostics) {
    const double expectedTime = static_cast<double>(acceptedStepCount)
        * settings.flow.timeStepSeconds;
    if (!std::isfinite(simulationTimeSeconds)
        || simulationTimeSeconds < 0.0
        || std::abs(simulationTimeSeconds - expectedTime)
            > 1.0e-12 * std::max(1.0, std::abs(expectedTime))) {
        return false;
    }
    if (acceptedStepCount == 0) {
        return diagnostics
            == fluid::PeriodicFlowStrangSubcyclingDiagnostics{};
    }
    return diagnostics.version == fluid::periodicFlowStrangSubcyclingVersion
        && diagnostics.requestedIntervalSeconds
            == settings.flow.timeStepSeconds
        && std::isfinite(diagnostics.substepSeconds)
        && diagnostics.substepSeconds > 0.0
        && diagnostics.plannedSubstepCount > 0
        && diagnostics.completedSubstepCount
            == diagnostics.plannedSubstepCount
        && diagnostics.substeps.size()
            == diagnostics.plannedSubstepCount
        && std::ranges::all_of(
            diagnostics.substeps,
            [](const auto& substep) {
                return substep.finite && substep.accepted;
            })
        && diagnostics.failureStage
            == fluid::PeriodicFlowStrangSubcyclingFailureStage::None
        && diagnostics.finite
        && diagnostics.accepted;
}

} // namespace

struct PeriodicFlowCaseCheckpoint::Detail {
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    fluid::MacVelocityField velocityMetersPerSecond;
    fluid::CellScalarField pressurePascals;
    fluid::PeriodicFlowStrangSubcyclingDiagnostics diagnostics;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
};

PeriodicFlowCase::PeriodicFlowCase()
    : grid_(makeGrid()),
      velocity_(makeInitialVelocity(grid_)),
      pressure_(grid_),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader PeriodicFlowCase::traceHeader() const {
    return {periodicFlowCaseChecksum, periodicFlowCaseSolverId};
}

viewer::DiagnosticFrame PeriodicFlowCase::advance() {
    fluid::MacVelocityField candidateVelocity = velocity_;
    fluid::CellScalarField candidatePressure = pressure_;
    fluid::PeriodicFlowStrangSubcyclingDiagnostics candidateDiagnostics =
        fluid::advancePeriodicFlowStrangSspRk2Subcycled(
            grid_, candidateVelocity, candidatePressure, stepSettings_);
    if (!candidateDiagnostics.accepted) {
        throw std::runtime_error(
            "periodic flow case rejected its requested outer interval");
    }
    const std::uint64_t candidateStep = acceptedStepCount_ + 1;
    const double candidateTime = simulationTimeSeconds_
        + stepSettings_.flow.timeStepSeconds;
    viewer::PeriodicFluidFrameContext context;
    context.sceneChecksum = periodicFlowCaseChecksum;
    context.solverCommit = periodicFlowCaseSolverId;
    context.step = candidateStep;
    context.simulationTimeSeconds = candidateTime;
    context.densityKgPerCubicMeter =
        stepSettings_.flow.densityKgPerCubicMeter;
    viewer::DiagnosticFrame frame = viewer::buildPeriodicFluidFrame(
        grid_, candidateVelocity, candidatePressure,
        candidateDiagnostics, context);
    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    diagnostics_ = std::move(candidateDiagnostics);
    acceptedStepCount_ = candidateStep;
    simulationTimeSeconds_ = candidateTime;
    return frame;
}

PeriodicFlowCaseCheckpoint PeriodicFlowCase::checkpoint() const {
    if (!fluid::isFinite(velocity_) || !fluid::isFinite(pressure_)
        || !validCommittedState(
            acceptedStepCount_, simulationTimeSeconds_,
            stepSettings_, diagnostics_)) {
        throw std::logic_error(
            "periodic flow case cannot checkpoint invalid committed state");
    }

    PeriodicFlowCaseCheckpoint result;
    result.cellCounts = grid_.cellCounts();
    result.lowerMeters = grid_.lowerMeters();
    result.upperMeters = grid_.upperMeters();
    result.scalarSampleCount = grid_.cellCount();
    result.acceptedStepCount = acceptedStepCount_;
    result.simulationTimeSeconds = simulationTimeSeconds_;
    result.detail = std::make_shared<PeriodicFlowCaseCheckpoint::Detail>(
        PeriodicFlowCaseCheckpoint::Detail{
            result.cellCounts,
            result.lowerMeters,
            result.upperMeters,
            velocity_,
            pressure_,
            diagnostics_,
            acceptedStepCount_,
            simulationTimeSeconds_,
        });
    return result;
}

void PeriodicFlowCase::restore(
    const PeriodicFlowCaseCheckpoint& checkpointValue) {
    if (checkpointValue.version != periodicFlowCaseCheckpointVersion
        || checkpointValue.caseDefinitionFingerprint
            != periodicFlowCaseDefinitionFingerprint
        || !checkpointValue.detail
        || checkpointValue.cellCounts != checkpointValue.detail->cellCounts
        || checkpointValue.lowerMeters != checkpointValue.detail->lowerMeters
        || checkpointValue.upperMeters != checkpointValue.detail->upperMeters
        || checkpointValue.scalarSampleCount != grid_.cellCount()
        || checkpointValue.acceptedStepCount
            != checkpointValue.detail->acceptedStepCount
        || checkpointValue.simulationTimeSeconds
            != checkpointValue.detail->simulationTimeSeconds
        || !gridMetadataMatches(
            grid_, checkpointValue.cellCounts,
            checkpointValue.lowerMeters, checkpointValue.upperMeters)
        || !validCommittedState(
            checkpointValue.detail->acceptedStepCount,
            checkpointValue.detail->simulationTimeSeconds,
            stepSettings_, checkpointValue.detail->diagnostics)) {
        throw std::invalid_argument(
            "periodic flow case checkpoint metadata is invalid");
    }

    fluid::MacVelocityField candidateVelocity =
        checkpointValue.detail->velocityMetersPerSecond;
    fluid::CellScalarField candidatePressure =
        checkpointValue.detail->pressurePascals;
    auto candidateDiagnostics = checkpointValue.detail->diagnostics;
    if (!candidateVelocity.matches(grid_)
        || !candidatePressure.matches(grid_)
        || !fluid::isFinite(candidateVelocity)
        || !fluid::isFinite(candidatePressure)) {
        throw std::invalid_argument(
            "periodic flow case checkpoint payload is invalid");
    }

    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    diagnostics_ = std::move(candidateDiagnostics);
    acceptedStepCount_ = checkpointValue.acceptedStepCount;
    simulationTimeSeconds_ = checkpointValue.simulationTimeSeconds;
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
