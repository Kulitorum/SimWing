#include "moving_porous_flow_case.h"

#include "fluid/planar_porous_sheet.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {

struct MovingPorousFlowCaseCheckpoint::Detail {
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::size_t pressureJumpCount = 0;
    fluid::MacVelocityField velocityMetersPerSecond;
    fluid::CellScalarField pressurePascals;
    fluid::SharpPressureJumpField pressureJumps;
    fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics diagnostics;
    fluid::MovingPorousFaceTopology porousTopology;
    double sheetPositionMeters = 0.0;
    double sheetVelocityMetersPerSecond = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
};

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

bool validCommittedState(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    const fluid::CellScalarField& pressure,
    const fluid::SharpPressureJumpField& pressureJumps,
    const fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics& diagnostics,
    const fluid::MovingPorousFaceTopology& topology,
    const double sheetPositionMeters,
    const double sheetVelocityMetersPerSecond,
    const std::uint64_t topologyRebaseCount,
    const std::uint64_t acceptedStepCount,
    const double simulationTimeSeconds,
    const fluid::PeriodicFlowStrangSspRk2Settings& settings) {
    if (!velocity.matches(grid) || !pressure.matches(grid)
        || !pressureJumps.matches(grid)
        || !fluid::isFinite(velocity) || !fluid::isFinite(pressure)
        || !std::isfinite(sheetPositionMeters)
        || sheetVelocityMetersPerSecond
            != movingPorousFlowSheetVelocityMetersPerSecond
        || !std::isfinite(simulationTimeSeconds)
        || simulationTimeSeconds < 0.0
        || topology.version != fluid::movingPorousFaceTopologyVersion
        || topology.axis != fluid::GridFaceAxis::X
        || topology.faceCoordinate >= grid.cellCounts().x
        || topology.periodicImage < 0) {
        return false;
    }
    const double expectedTime = static_cast<double>(acceptedStepCount)
        * settings.timeStepSeconds;
    const double expectedPosition = movingPorousFlowInitialSheetPositionMeters
        + expectedTime * sheetVelocityMetersPerSecond;
    if (!std::isfinite(expectedTime) || !std::isfinite(expectedPosition)) {
        return false;
    }
    const double tolerance = 4.0e-14
        + 8.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, static_cast<double>(acceptedStepCount))
            * std::max({1.0, std::abs(expectedTime),
                        std::abs(expectedPosition)});
    if (std::abs(simulationTimeSeconds - expectedTime) > tolerance
        || std::abs(sheetPositionMeters - expectedPosition) > tolerance) {
        return false;
    }
    const auto image = static_cast<std::uint64_t>(topology.periodicImage);
    const auto counts = grid.cellCounts();
    if (image > (std::numeric_limits<std::uint64_t>::max()
                 - topology.faceCoordinate) / counts.x) {
        return false;
    }
    const std::uint64_t unwrappedFace = image * counts.x
        + topology.faceCoordinate;
    if (unwrappedFace < movingPorousFlowInitialFaceCoordinate
        || unwrappedFace - movingPorousFlowInitialFaceCoordinate
            != topologyRebaseCount) {
        return false;
    }
    if (acceptedStepCount == 0) {
        return topologyRebaseCount == 0
            && topology.faceCoordinate
                == movingPorousFlowInitialFaceCoordinate
            && topology.periodicImage == 0
            && pressureJumps.empty()
            && diagnostics
                == fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics{}
            && velocity == fluid::MacVelocityField(grid)
            && pressure == fluid::CellScalarField(grid);
    }
    if (pressureJumps.faceCount() != 2 * counts.y * counts.z
        || diagnostics.version
            != fluid::movingPlanarPorousFlowStrangSspRk2Version
        || !diagnostics.finite || !diagnostics.accepted
        || !diagnostics.flow.finite || !diagnostics.flow.accepted
        || diagnostics.secondHalfSheet.topology != topology
        || diagnostics.firstHalfSheet.surfaceNormalVelocityMetersPerSecond
            != sheetVelocityMetersPerSecond
        || diagnostics.secondHalfSheet.surfaceNormalVelocityMetersPerSecond
            != sheetVelocityMetersPerSecond) {
        return false;
    }
    try {
        static_cast<void>(fluid::movingPorousCrossingFraction(
            grid, diagnostics.firstHalfSheet.topology,
            diagnostics.firstHalfSheet.physicalPlaneCoordinateMeters));
        static_cast<void>(fluid::movingPorousCrossingFraction(
            grid, diagnostics.secondHalfSheet.topology,
            diagnostics.secondHalfSheet.physicalPlaneCoordinateMeters));
    } catch (const std::exception&) {
        return false;
    }
    return pressureJumps == makeAcceptedJumps(
        grid, diagnostics.flow.secondHalfPorous, makePump(grid));
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
          movingPorousFlowInitialFaceCoordinate,
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

MovingPorousFlowCaseCheckpoint MovingPorousFlowCase::checkpoint() const {
    if (!validCommittedState(
            grid_, velocity_, pressure_, pressureJumps_, diagnostics_,
            porousTopology_, sheetPositionMeters_,
            sheetVelocityMetersPerSecond_, topologyRebaseCount_,
            acceptedStepCount_, simulationTimeSeconds_, flowSettings_)) {
        throw std::logic_error(
            "moving porous-flow case cannot checkpoint invalid state");
    }
    MovingPorousFlowCaseCheckpoint result;
    result.cellCounts = grid_.cellCounts();
    result.lowerMeters = grid_.lowerMeters();
    result.upperMeters = grid_.upperMeters();
    result.scalarSampleCount = grid_.cellCount();
    result.pressureJumpCount = pressureJumps_.faceCount();
    result.acceptedStepCount = acceptedStepCount_;
    result.simulationTimeSeconds = simulationTimeSeconds_;
    result.sheetPositionMeters = sheetPositionMeters_;
    result.sheetVelocityMetersPerSecond = sheetVelocityMetersPerSecond_;
    result.topologyRebaseCount = topologyRebaseCount_;
    result.porousTopology = porousTopology_;
    result.detail = std::make_shared<
        MovingPorousFlowCaseCheckpoint::Detail>(
            MovingPorousFlowCaseCheckpoint::Detail{
                result.cellCounts,
                result.lowerMeters,
                result.upperMeters,
                result.scalarSampleCount,
                result.pressureJumpCount,
                velocity_,
                pressure_,
                pressureJumps_,
                diagnostics_,
                porousTopology_,
                sheetPositionMeters_,
                sheetVelocityMetersPerSecond_,
                topologyRebaseCount_,
                acceptedStepCount_,
                simulationTimeSeconds_,
            });
    return result;
}

void MovingPorousFlowCase::restore(
    const MovingPorousFlowCaseCheckpoint& checkpointValue) {
    if (checkpointValue.version
            != movingPorousFlowCaseCheckpointVersion
        || checkpointValue.caseDefinitionFingerprint
            != movingPorousFlowCaseDefinitionFingerprint
        || !checkpointValue.detail) {
        throw std::invalid_argument(
            "moving porous-flow checkpoint metadata is invalid");
    }
    const auto& detail = *checkpointValue.detail;
    if (checkpointValue.cellCounts != detail.cellCounts
        || checkpointValue.lowerMeters != detail.lowerMeters
        || checkpointValue.upperMeters != detail.upperMeters
        || checkpointValue.scalarSampleCount != detail.scalarSampleCount
        || checkpointValue.pressureJumpCount != detail.pressureJumpCount
        || checkpointValue.acceptedStepCount != detail.acceptedStepCount
        || checkpointValue.simulationTimeSeconds
            != detail.simulationTimeSeconds
        || checkpointValue.sheetPositionMeters
            != detail.sheetPositionMeters
        || checkpointValue.sheetVelocityMetersPerSecond
            != detail.sheetVelocityMetersPerSecond
        || checkpointValue.topologyRebaseCount
            != detail.topologyRebaseCount
        || checkpointValue.porousTopology != detail.porousTopology
        || checkpointValue.cellCounts != grid_.cellCounts()
        || checkpointValue.lowerMeters != grid_.lowerMeters()
        || checkpointValue.upperMeters != grid_.upperMeters()
        || checkpointValue.scalarSampleCount != grid_.cellCount()
        || checkpointValue.pressureJumpCount
            != detail.pressureJumps.faceCount()
        || !validCommittedState(
            grid_, detail.velocityMetersPerSecond,
            detail.pressurePascals, detail.pressureJumps,
            detail.diagnostics, detail.porousTopology,
            detail.sheetPositionMeters,
            detail.sheetVelocityMetersPerSecond,
            detail.topologyRebaseCount, detail.acceptedStepCount,
            detail.simulationTimeSeconds, flowSettings_)) {
        throw std::invalid_argument(
            "moving porous-flow checkpoint state is invalid");
    }

    auto candidateVelocity = detail.velocityMetersPerSecond;
    auto candidatePressure = detail.pressurePascals;
    auto candidateJumps = detail.pressureJumps;
    auto candidateDiagnostics = detail.diagnostics;
    velocity_ = std::move(candidateVelocity);
    pressure_ = std::move(candidatePressure);
    pressureJumps_ = std::move(candidateJumps);
    diagnostics_ = std::move(candidateDiagnostics);
    porousTopology_ = checkpointValue.porousTopology;
    sheetPositionMeters_ = checkpointValue.sheetPositionMeters;
    sheetVelocityMetersPerSecond_ =
        checkpointValue.sheetVelocityMetersPerSecond;
    topologyRebaseCount_ = checkpointValue.topologyRebaseCount;
    acceptedStepCount_ = checkpointValue.acceptedStepCount;
    simulationTimeSeconds_ = checkpointValue.simulationTimeSeconds;
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

const fluid::PeriodicFlowStrangSspRk2Settings&
MovingPorousFlowCase::stepSettings() const noexcept {
    return flowSettings_;
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
