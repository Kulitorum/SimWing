#include "porous_flow_case.h"
#include "fluid/planar_porous_sheet.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

fluid::PeriodicCartesianGrid makeGrid() {
    return {{16, 4, 3}, {}, {4.0, 2.0, 1.5}};
}

fluid::ProjectionSettings makeProjectionSettings() {
    fluid::ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 1.0 / 60.0;
    settings.absoluteResidualTolerance = 1.0e-10;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 2000;
    return settings;
}

fluid::PeriodicFlowStrangSspRk2Settings makeFlowSettings() {
    fluid::PeriodicFlowStrangSspRk2Settings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.kinematicViscositySquareMetersPerSecond = 1.5e-5;
    settings.timeStepSeconds = 1.0 / 60.0;
    settings.advectionReconstruction =
        fluid::VariableMacReconstruction::DonorCell;
    settings.maximumLocalOutgoingCourantNumber = 1.0;
    settings.advectionAbsoluteDivergenceTolerancePerSecond = 1.0e-11;
    settings.advectionRelativeDivergenceTolerance = 1.0e-12;
    settings.maximumDiffusionNumber = 0.5;
    settings.projectionAbsoluteResidualTolerance = 1.0e-10;
    settings.projectionRelativeResidualTolerance = 1.0e-13;
    settings.projectionMaximumIterations = 2000;
    settings.absoluteMomentumToleranceNewtonSeconds = 2.0e-11;
    settings.relativeMomentumTolerance = 1.0e-12;
    settings.absoluteEnergyToleranceJoules = 2.0e-11;
    settings.relativeEnergyTolerance = 1.0e-12;
    return settings;
}

fluid::PorousPlugFlowSettings makePlugSettings() {
    fluid::PorousPlugFlowSettings settings;
    settings.resistance = {100.0, 25.0};
    settings.densityKgPerCubicMeter = 1.2;
    settings.flowLengthMeters = 4.0;
    settings.crossSectionAreaSquareMeters = 3.0;
    settings.drivingPressureRisePascals = 250.0;
    settings.timeStepSeconds = 1.0 / 60.0;
    return settings;
}

std::vector<fluid::PorousGridFaceCrossing> makePorousPlane(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::DarcyForchheimerResistance& resistance) {
    return fluid::makePlanarPorousSheetCrossings(
        grid,
        {
            110,
            1,
            2,
            {
                fluid::movingPorousFaceTopologyVersion,
                fluid::GridFaceAxis::X,
                4,
                0,
            },
            1.0,
            0.0,
            resistance,
        });
}

std::vector<fluid::GridFacePressureJump> completePressureCircuit(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PorousPressureJumpField& porous) {
    if (porous.samples().empty()) {
        throw std::invalid_argument(
            "porous pressure circuit requires resolved fabric tiles");
    }
    const double pumpJump = -porous.samples().front()
        .pressureJump.pressureJumpPascals;
    for (const auto& sample : porous.samples()) {
        if (sample.pressureJump.pressureJumpPascals != -pumpJump) {
            throw std::invalid_argument(
                "porous pressure circuit requires one uniform endpoint loss");
        }
    }
    std::vector<fluid::GridFacePressureJump> result(
        porous.pressureJumps().faces().begin(),
        porous.pressureJumps().faces().end());
    const auto counts = grid.cellCounts();
    result.reserve(result.size() + counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                120, 2, 1, fluid::GridFaceAxis::X,
                12, j, k, pumpJump,
                0.5});
        }
    }
    return result;
}

void addGlobalScalar(viewer::DiagnosticFrame& frame,
                     const char* name,
                     const char* unit,
                     const double value) {
    frame.scalarFields.push_back({
        name, unit, viewer::FieldAssociation::Global, {value}});
}

} // namespace

PorousFlowCase::PorousFlowCase()
    : grid_(makeGrid()),
      velocity_(grid_),
      pressure_(grid_),
      pressureJumps_(grid_),
      stepSettings_(makeProjectionSettings()),
      flowSettings_(makeFlowSettings()),
      plugSettings_(makePlugSettings()) {
    if (stepSettings_.densityKgPerCubicMeter
            != plugSettings_.densityKgPerCubicMeter
        || stepSettings_.timeStepSeconds
            != plugSettings_.timeStepSeconds
        || flowSettings_.densityKgPerCubicMeter
            != plugSettings_.densityKgPerCubicMeter
        || flowSettings_.timeStepSeconds
            != plugSettings_.timeStepSeconds) {
        throw std::logic_error(
            "porous worker projection and plug settings disagree");
    }
}

viewer::TraceHeader PorousFlowCase::traceHeader() const {
    return {porousFlowCaseChecksum, porousFlowCaseSolverId};
}

viewer::DiagnosticFrame PorousFlowCase::advance() {
    double candidateFlowVelocity = flowVelocityMetersPerSecond_;
    const fluid::PorousPlugFlowDiagnostics candidatePlugDiagnostics =
        fluid::advancePorousPlugFlow(
            candidateFlowVelocity, plugSettings_);

    fluid::MacVelocityField candidateVelocity(grid_);
    std::ranges::fill(
        candidateVelocity.xFaces(), candidateFlowVelocity);
    const fluid::PorousPressureJumpField porous(
        grid_, candidateVelocity,
        makePorousPlane(grid_, plugSettings_.resistance));
    fluid::SharpPressureJumpField candidatePressureJumps(
        grid_, completePressureCircuit(grid_, porous));
    fluid::CellScalarField candidatePressure(grid_);
    const fluid::PeriodicFlowStrangSspRk2Diagnostics
        candidateFlowDiagnostics =
            fluid::advancePeriodicFlowStrangSspRk2(
                grid_, candidateVelocity, candidatePressure,
                candidatePressureJumps, flowSettings_);
    if (!candidateFlowDiagnostics.accepted) {
        throw std::runtime_error(
            "porous-flow case endpoint fluid step was not accepted");
    }
    const fluid::ProjectionDiagnostics candidateDiagnostics =
        candidateFlowDiagnostics.projectedAdvection.secondProjection;
    double maximumVelocityError = 0.0;
    for (const double sample : candidateVelocity.xFaces()) {
        maximumVelocityError = std::max(
            maximumVelocityError,
            std::abs(sample - candidateFlowVelocity));
    }
    fluid::CellScalarField candidateDivergence(grid_);
    fluid::computeDivergence(
        grid_, candidateVelocity, candidateDivergence);
    if (maximumVelocityError > 2.0e-13
        || fluid::maximumAbsoluteValue(candidateDivergence) > 1.0e-12) {
        throw std::runtime_error(
            "porous-flow endpoint projection changed the uniform plug");
    }

    const std::uint64_t candidateStep = acceptedStepCount_ + 1;
    const double candidateTime = simulationTimeSeconds_
        + stepSettings_.timeStepSeconds;
    viewer::PressureJumpFrameContext context;
    context.sceneChecksum = porousFlowCaseChecksum;
    context.solverCommit = porousFlowCaseSolverId;
    context.step = candidateStep;
    context.simulationTimeSeconds = candidateTime;
    context.timeStepSeconds = stepSettings_.timeStepSeconds;
    context.densityKgPerCubicMeter =
        stepSettings_.densityKgPerCubicMeter;
    viewer::DiagnosticFrame frame = viewer::buildPressureJumpFrame(
        grid_, candidateVelocity, candidatePressure,
        candidatePressureJumps, candidateDiagnostics, context);
    addGlobalScalar(
        frame, "porous normal velocity", "m/s",
        candidateFlowVelocity);
    addGlobalScalar(
        frame, "porous endpoint pressure drop", "Pa",
        candidatePlugDiagnostics.endpointPressureDropPascals);
    addGlobalScalar(
        frame, "porous midpoint pressure drop", "Pa",
        candidatePlugDiagnostics.midpointPressureDropPascals);
    addGlobalScalar(
        frame, "driving pressure rise", "Pa",
        plugSettings_.drivingPressureRisePascals);
    addGlobalScalar(
        frame, "porous dissipation", "J",
        candidatePlugDiagnostics.porousDissipationJoules);
    addGlobalScalar(
        frame, "driving pressure work", "J",
        candidatePlugDiagnostics.drivingPressureWorkJoules);
    addGlobalScalar(
        frame, "porous momentum residual", "N*s",
        candidatePlugDiagnostics.momentumResidualNewtonSeconds);
    addGlobalScalar(
        frame, "porous energy residual", "J",
        candidatePlugDiagnostics.energyResidualJoules);
    viewer::ProtocolError frameError;
    if (!viewer::validateFrame(frame, &frameError)) {
        throw std::runtime_error(
            "porous-flow enriched frame is invalid: "
            + frameError.message);
    }

    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    pressureJumps_ = std::move(candidatePressureJumps);
    diagnostics_ = candidateDiagnostics;
    flowDiagnostics_ = candidateFlowDiagnostics;
    plugDiagnostics_ = candidatePlugDiagnostics;
    flowVelocityMetersPerSecond_ = candidateFlowVelocity;
    acceptedStepCount_ = candidateStep;
    simulationTimeSeconds_ = candidateTime;
    return frame;
}

const fluid::PeriodicCartesianGrid&
PorousFlowCase::grid() const noexcept {
    return grid_;
}

const fluid::MacVelocityField&
PorousFlowCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField&
PorousFlowCase::pressure() const noexcept {
    return pressure_;
}

const fluid::SharpPressureJumpField&
PorousFlowCase::pressureJumps() const noexcept {
    return pressureJumps_;
}

const fluid::ProjectionSettings&
PorousFlowCase::stepSettings() const noexcept {
    return stepSettings_;
}

const fluid::PorousPlugFlowSettings&
PorousFlowCase::plugSettings() const noexcept {
    return plugSettings_;
}

const fluid::ProjectionDiagnostics&
PorousFlowCase::diagnostics() const noexcept {
    return diagnostics_;
}

const fluid::PeriodicFlowStrangSspRk2Diagnostics&
PorousFlowCase::flowDiagnostics() const noexcept {
    return flowDiagnostics_;
}

const fluid::PorousPlugFlowDiagnostics&
PorousFlowCase::plugDiagnostics() const noexcept {
    return plugDiagnostics_;
}

std::uint64_t PorousFlowCase::acceptedStepCount() const noexcept {
    return acceptedStepCount_;
}

double PorousFlowCase::simulationTimeSeconds() const noexcept {
    return simulationTimeSeconds_;
}

} // namespace simwing::fsi
