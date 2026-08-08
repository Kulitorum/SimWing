#include "moving_porous_flow_case.h"

#include "fluid/planar_porous_sheet.h"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t sheetSurfaceStableId = 100;
constexpr std::uint64_t pumpSurfaceStableId = 200;
constexpr std::uint64_t minusRegionStableId = 1;
constexpr std::uint64_t plusRegionStableId = 2;
constexpr double pumpPressureJumpPascals = 20.0;
constexpr fluid::DarcyForchheimerResistance sheetResistance{
    10.0, 0.0};

fluid::PeriodicCartesianGrid makeGrid() {
    return {{4, 3, 2}, {}, {4.0, 3.0, 2.0}};
}

fluid::PorousIterationSettings makePorousIteration() {
    fluid::PorousIterationSettings result;
    result.constitutiveEvaluation =
        fluid::PorousConstitutiveEvaluation::Midpoint;
    result.absoluteNormalVelocityToleranceMetersPerSecond = 1.0e-13;
    result.relativeNormalVelocityTolerance = 1.0e-13;
    result.absolutePressureJumpTolerancePascals = 1.0e-12;
    result.relativePressureJumpTolerance = 1.0e-13;
    result.relaxation = 0.5;
    result.maximumNonlinearIterations = 200;
    return result;
}

fluid::PeriodicFlowStrangSspRk2Settings makeFlowSettings() {
    fluid::PeriodicFlowStrangSspRk2Settings result;
    result.densityKgPerCubicMeter = 1.0;
    result.kinematicViscositySquareMetersPerSecond = 0.0;
    result.timeStepSeconds = 0.1;
    result.advectionReconstruction =
        fluid::VariableMacReconstruction::DonorCell;
    result.maximumLocalOutgoingCourantNumber = 1.0;
    result.advectionAbsoluteDivergenceTolerancePerSecond = 1.0e-11;
    result.advectionRelativeDivergenceTolerance = 1.0e-12;
    result.maximumDiffusionNumber = 0.5;
    result.projectionAbsoluteResidualTolerance = 1.0e-12;
    result.projectionRelativeResidualTolerance = 1.0e-13;
    result.projectionMaximumIterations = 1000;
    result.absoluteMomentumToleranceNewtonSeconds = 2.0e-10;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 2.0e-10;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

fluid::SharpPressureJumpField makePump(
    const fluid::PeriodicCartesianGrid& grid) {
    const auto counts = grid.cellCounts();
    std::vector<fluid::GridFacePressureJump> faces;
    faces.reserve(counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                pumpSurfaceStableId,
                plusRegionStableId,
                minusRegionStableId,
                fluid::GridFaceAxis::X,
                2,
                j,
                k,
                pumpPressureJumpPascals,
                0.5,
            });
        }
    }
    return {grid, std::move(faces)};
}

fluid::PlanarPorousSheetDefinition makeSheet(
    const fluid::MovingPorousFaceTopology& topology,
    const double positionMeters,
    const double velocityMetersPerSecond) {
    return {
        sheetSurfaceStableId,
        minusRegionStableId,
        plusRegionStableId,
        topology,
        positionMeters,
        velocityMetersPerSecond,
        sheetResistance,
    };
}

fluid::SharpPressureJumpField makeAcceptedJumps(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PorousProjectionDiagnostics& porous,
    const fluid::SharpPressureJumpField& pump) {
    std::vector<fluid::GridFacePressureJump> faces;
    faces.reserve(porous.samples.size() + pump.faceCount());
    for (const auto& sample : porous.samples) {
        faces.push_back(sample.pressureJump);
    }
    faces.insert(
        faces.end(), pump.faces().begin(), pump.faces().end());
    return {grid, std::move(faces)};
}

void addGlobalScalar(viewer::DiagnosticFrame& frame,
                     const char* name,
                     const char* unit,
                     const double value) {
    frame.scalarFields.push_back({
        name, unit, viewer::FieldAssociation::Global, {value}});
}

} // namespace

MovingPorousFlowCase::MovingPorousFlowCase()
    : grid_(makeGrid()),
      velocity_(grid_),
      pressure_(grid_),
      pressureJumps_(grid_),
      porousIteration_(makePorousIteration()),
      flowSettings_(makeFlowSettings()),
      porousTopology_{
          fluid::movingPorousFaceTopologyVersion,
          fluid::GridFaceAxis::X,
          3,
          0,
      } {}

viewer::TraceHeader MovingPorousFlowCase::traceHeader() const {
    return {movingPorousFlowCaseChecksum, movingPorousFlowCaseSolverId};
}

viewer::DiagnosticFrame MovingPorousFlowCase::advance() {
    const double timeStep = flowSettings_.timeStepSeconds;
    const double firstPosition = sheetPositionMeters_
        + 0.25 * timeStep * sheetVelocityMetersPerSecond_;
    const double secondPosition = sheetPositionMeters_
        + 0.75 * timeStep * sheetVelocityMetersPerSecond_;
    const auto firstSelection = fluid::selectMovingPorousTopology(
        grid_, porousTopology_, firstPosition);
    const auto secondSelection = fluid::selectMovingPorousTopology(
        grid_, firstSelection.topology, secondPosition);
    const fluid::MovingPlanarPorousSheetStrangStages stages{
        makeSheet(
            firstSelection.topology, firstPosition,
            sheetVelocityMetersPerSecond_),
        makeSheet(
            secondSelection.topology, secondPosition,
            sheetVelocityMetersPerSecond_),
    };
    const fluid::SharpPressureJumpField pump = makePump(grid_);
    auto candidateVelocity = velocity_;
    auto candidatePressure = pressure_;
    const auto candidateDiagnostics =
        fluid::advancePeriodicFlowStrangSspRk2WithMovingPlanarPorousSheet(
            grid_, candidateVelocity, candidatePressure,
            stages, pump, porousIteration_, flowSettings_);
    if (!candidateDiagnostics.accepted) {
        throw std::runtime_error(
            "moving porous-flow case fluid step was not accepted");
    }
    fluid::SharpPressureJumpField candidateJumps = makeAcceptedJumps(
        grid_, candidateDiagnostics.flow.secondHalfPorous, pump);
    const double candidateSheetPosition = sheetPositionMeters_
        + timeStep * sheetVelocityMetersPerSecond_;
    if (!std::isfinite(candidateSheetPosition)) {
        throw std::runtime_error(
            "moving porous-flow sheet position is not finite");
    }
    const std::uint64_t candidateStep = acceptedStepCount_ + 1;
    const double candidateTime = simulationTimeSeconds_ + timeStep;
    viewer::PressureJumpFrameContext context;
    context.sceneChecksum = movingPorousFlowCaseChecksum;
    context.solverCommit = movingPorousFlowCaseSolverId;
    context.step = candidateStep;
    context.simulationTimeSeconds = candidateTime;
    context.timeStepSeconds = timeStep;
    context.densityKgPerCubicMeter =
        flowSettings_.densityKgPerCubicMeter;
    viewer::DiagnosticFrame frame = viewer::buildPressureJumpFrame(
        grid_, candidateVelocity, candidatePressure, candidateJumps,
        candidateDiagnostics.flow.secondHalfPorous.projection, context);
    frame.conservation.totalMomentumNewtonSeconds = {
        candidateDiagnostics.flow.momentumAfterNewtonSeconds.x,
        candidateDiagnostics.flow.momentumAfterNewtonSeconds.y,
        candidateDiagnostics.flow.momentumAfterNewtonSeconds.z,
    };
    frame.conservation.totalEnergyJoules =
        candidateDiagnostics.flow.kineticEnergyAfterJoules;
    addGlobalScalar(
        frame, "sheet position", "m", candidateSheetPosition);
    addGlobalScalar(
        frame, "sheet velocity", "m/s",
        sheetVelocityMetersPerSecond_);
    addGlobalScalar(
        frame, "first porous face", "1",
        static_cast<double>(firstSelection.topology.faceCoordinate));
    addGlobalScalar(
        frame, "first porous image", "1",
        static_cast<double>(firstSelection.topology.periodicImage));
    addGlobalScalar(
        frame, "second porous face", "1",
        static_cast<double>(secondSelection.topology.faceCoordinate));
    addGlobalScalar(
        frame, "second porous image", "1",
        static_cast<double>(secondSelection.topology.periodicImage));
    addGlobalScalar(
        frame, "porous kinematic residual", "m",
        candidateDiagnostics.kinematicResidualMeters);
    addGlobalScalar(
        frame, "porous dissipation", "J",
        candidateDiagnostics.flow.porousDissipationJoules);
    addGlobalScalar(
        frame, "pressure-jump work", "J",
        candidateDiagnostics.flow.pressureJumpWorkToFluidJoules);
    addGlobalScalar(
        frame, "flow momentum residual", "N*s",
        candidateDiagnostics.flow.momentumResidualNormNewtonSeconds);
    viewer::ProtocolError frameError;
    if (!viewer::validateFrame(frame, &frameError)) {
        throw std::runtime_error(
            "moving porous-flow frame is invalid: "
            + frameError.message);
    }

    const auto didRebase = [](const auto& selection) {
        return selection.rebaseDirection
            != fluid::PorousTopologyRebaseDirection::None;
    };
    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    pressureJumps_ = std::move(candidateJumps);
    diagnostics_ = candidateDiagnostics;
    porousTopology_ = secondSelection.topology;
    sheetPositionMeters_ = candidateSheetPosition;
    topologyRebaseCount_ += static_cast<std::uint64_t>(
        didRebase(firstSelection))
        + static_cast<std::uint64_t>(didRebase(secondSelection));
    acceptedStepCount_ = candidateStep;
    simulationTimeSeconds_ = candidateTime;
    return frame;
}

const fluid::PeriodicCartesianGrid&
MovingPorousFlowCase::grid() const noexcept {
    return grid_;
}

const fluid::MacVelocityField&
MovingPorousFlowCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField&
MovingPorousFlowCase::pressure() const noexcept {
    return pressure_;
}

const fluid::SharpPressureJumpField&
MovingPorousFlowCase::pressureJumps() const noexcept {
    return pressureJumps_;
}

const fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics&
MovingPorousFlowCase::diagnostics() const noexcept {
    return diagnostics_;
}

const fluid::MovingPorousFaceTopology&
MovingPorousFlowCase::porousTopology() const noexcept {
    return porousTopology_;
}

double MovingPorousFlowCase::sheetPositionMeters() const noexcept {
    return sheetPositionMeters_;
}

double MovingPorousFlowCase::sheetVelocityMetersPerSecond() const noexcept {
    return sheetVelocityMetersPerSecond_;
}

std::uint64_t MovingPorousFlowCase::topologyRebaseCount() const noexcept {
    return topologyRebaseCount_;
}

std::uint64_t MovingPorousFlowCase::acceptedStepCount() const noexcept {
    return acceptedStepCount_;
}

double MovingPorousFlowCase::simulationTimeSeconds() const noexcept {
    return simulationTimeSeconds_;
}

} // namespace simwing::fsi
