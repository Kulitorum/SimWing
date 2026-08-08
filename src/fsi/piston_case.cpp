#include "piston_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t rightPistonSurfaceStableId = 200;
constexpr double pistonPressurePascals = 110.0;
constexpr double fluidDensityKgPerCubicMeter = 1.2;
constexpr double fluidMassKilograms = 28.8;
constexpr double structuralMassKilograms = 6000.0;
constexpr double strongStructuralMassKilograms = 6.0;

StructureDefinition makeDefinition(
    const double totalMassKilograms = structuralMassKilograms) {
    StructureDefinition definition;
    // The two-triangle barycentric load shares are 2:1:2:1. Mass follows the
    // same ratio so uniform pressure preserves rigid piston translation.
    definition.nodes = {
        {{3.0, 0.0, 0.0}, totalMassKilograms / 3.0, false},
        {{3.0, 2.0, 0.0}, totalMassKilograms / 6.0, false},
        {{3.0, 2.0, 3.0}, totalMassKilograms / 3.0, false},
        {{3.0, 0.0, 3.0}, totalMassKilograms / 6.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

std::vector<CouplingSurfaceNodeDefinition> makeCouplingNodes() {
    return {{40, 3}, {10, 0}, {30, 2}, {20, 1}};
}

std::vector<CouplingSurfaceTriangleDefinition> makeCouplingTriangles() {
    return {
        {2000, {10, 30, 40}},
        {1000, {10, 20, 30}},
    };
}

viewer::StructureFrameMappingDefinition makeFrameMapping() {
    viewer::StructureFrameMappingDefinition mapping;
    mapping.vertexStableIds = {10, 20, 30, 40};
    mapping.triangles = {{1000, 2, 1}, {2000, 2, 1}};
    return mapping;
}

StructureStepSettings makeStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 120.0;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.cableConstraintSweepPairs = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    settings.workerThreads = 0;
    return settings;
}

fluid::PeriodicCartesianGrid pistonGrid() {
    return fluid::PeriodicCartesianGrid(
        {8, 2, 3}, {}, {4.0, 2.0, 3.0});
}

std::vector<fluid::GridFaceMovingInterface> slabFaces(
    const fluid::PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond) {
    std::vector<fluid::GridFaceMovingInterface> faces;
    const auto counts = grid.cellCounts();
    faces.reserve(2 * counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                100, 1, 2, fluid::GridFaceAxis::X, 2, j, k,
                speedMetersPerSecond});
            faces.push_back({
                rightPistonSurfaceStableId, 2, 1,
                fluid::GridFaceAxis::X, 6, j, k,
                speedMetersPerSecond});
        }
    }
    return faces;
}

fluid::MovingInterfaceProjectionDiagnostics solveFluidSample(
    const double speedMetersPerSecond,
    const double timeStepSeconds) {
    const auto grid = pistonGrid();
    fluid::FaceAlignedMovingInterface interfaces(
        grid, slabFaces(grid, speedMetersPerSecond));
    fluid::MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), speedMetersPerSecond);
    fluid::CellScalarField pressure(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                pressure.values()[grid.cellIndex(i, j, k)] =
                    i >= 2 && i < 6 ? pistonPressurePascals : 0.0;
            }
        }
    }
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter =
        fluidDensityKgPerCubicMeter;
    settings.projection.timeStepSeconds = timeStepSeconds;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    const auto diagnostics = fluid::projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, settings);
    if (!diagnostics.projection.converged || !diagnostics.finite) {
        throw std::runtime_error(
            "coupled piston fluid projection was not accepted");
    }
    return diagnostics;
}

fluid::MovingInterfaceFluidState initialFluidState(
    const fluid::PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond,
    const double timeStepSeconds) {
    fluid::MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), speedMetersPerSecond);
    fluid::CellScalarField pressure(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                pressure.values()[grid.cellIndex(i, j, k)] =
                    i >= 2 && i < 6 ? pistonPressurePascals : 0.0;
            }
        }
    }
    fluid::FaceAlignedMovingInterface interfaces(
        grid, slabFaces(grid, speedMetersPerSecond));
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter =
        fluidDensityKgPerCubicMeter;
    settings.projection.timeStepSeconds = timeStepSeconds;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    const auto diagnostics = fluid::projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, settings);
    if (!diagnostics.projection.converged || !diagnostics.finite) {
        throw std::runtime_error(
            "strong piston initial fluid projection was not accepted");
    }
    return {
        std::move(velocity),
        std::move(pressure),
        std::move(interfaces),
        diagnostics,
    };
}

fluid::MovingInterfaceProjectionSettings fluidProjectionSettings(
    const double timeStepSeconds) {
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter =
        fluidDensityKgPerCubicMeter;
    settings.projection.timeStepSeconds = timeStepSeconds;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return settings;
}

const fluid::MovingInterfaceSurfaceDiagnostics& rightSurface(
    const fluid::MovingInterfaceProjectionDiagnostics& diagnostics) {
    const auto found = std::ranges::lower_bound(
        diagnostics.surfaces, rightPistonSurfaceStableId, {},
        &fluid::MovingInterfaceSurfaceDiagnostics::stableId);
    if (found == diagnostics.surfaces.end()
        || found->stableId != rightPistonSurfaceStableId) {
        throw std::logic_error(
            "coupled piston fluid result lost its right surface");
    }
    return *found;
}

std::vector<fluid::MovingInterfaceFaceDiagnostics> referenceFluidFaces(
    const double timeStepSeconds) {
    return solveFluidSample(0.0, timeStepSeconds).faces;
}

double maximumVelocityDifference(
    const std::vector<StructureNodeState>& states,
    const double expectedSpeedMetersPerSecond) {
    double maximum = 0.0;
    for (const auto& state : states) {
        maximum = std::max(maximum, std::abs(
            state.velocityMetersPerSecond.x
                - expectedSpeedMetersPerSecond));
        maximum = std::max(maximum,
                           std::abs(state.velocityMetersPerSecond.y));
        maximum = std::max(maximum,
                           std::abs(state.velocityMetersPerSecond.z));
    }
    return maximum;
}

double acceptanceDisplacementResidual(
    const std::vector<StructureNodeState>& before,
    const std::vector<StructureNodeState>& accepted,
    const double expectedSpeedMetersPerSecond,
    const double timeStepSeconds) {
    if (before.size() != accepted.size()) {
        throw std::logic_error(
            "coupled piston acceptance changed the node count");
    }
    double maximum = maximumVelocityDifference(
        accepted, expectedSpeedMetersPerSecond) * timeStepSeconds;
    for (std::size_t index = 0; index < before.size(); ++index) {
        maximum = std::max(maximum, std::abs(
            accepted[index].positionMeters.x
                - before[index].positionMeters.x
                - expectedSpeedMetersPerSecond * timeStepSeconds));
        maximum = std::max(maximum, std::abs(
            accepted[index].positionMeters.y
                - before[index].positionMeters.y));
        maximum = std::max(maximum, std::abs(
            accepted[index].positionMeters.z
                - before[index].positionMeters.z));
    }
    return maximum;
}

viewer::Vec3d toViewer(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

StructureVector3 scale(const StructureVector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

void appendPistonFields(
    viewer::DiagnosticFrame& frame,
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const PlanarFaceResolvedBridgeDiagnostics& bridgeDiagnostics,
    const TimeIntegratedTransferDiagnostics& integratedDiagnostics) {
    frame.scalarFields.push_back({
        "interface.pressure_traction", "Pa",
        viewer::FieldAssociation::Triangle,
        {pistonPressurePascals, pistonPressurePascals}});
    frame.scalarFields.push_back({
        "fluid.outside_pressure_mean", "Pa",
        viewer::FieldAssociation::Global,
        {fluidDiagnostics.regions[0].pressureMeanAfterPascals}});
    frame.scalarFields.push_back({
        "fluid.inside_pressure_mean", "Pa",
        viewer::FieldAssociation::Global,
        {fluidDiagnostics.regions[1].pressureMeanAfterPascals}});
    frame.scalarFields.push_back({
        "interface.step_work", "J", viewer::FieldAssociation::Global,
        {integratedDiagnostics.integratedSurfaceWorkJoules}});
    frame.scalarFields.push_back({
        "fluid.divergence_l2", "1/s", viewer::FieldAssociation::Global,
        {fluidDiagnostics.projection.divergenceL2AfterPerSecond}});
    frame.vectorFields.push_back({
        "interface.uniform_pressure_traction", "Pa",
        viewer::FieldAssociation::Global,
        {toViewer(scale(
            bridgeDiagnostics.fluidPressureForceNewtons,
            1.0 / bridgeDiagnostics.fluidAreaSquareMeters))}});
    frame.vectorFields.push_back({
        "interface.step_impulse", "N*s", viewer::FieldAssociation::Global,
        {toViewer(
            integratedDiagnostics.integratedSurfaceImpulseNewtonSeconds)}});
}

} // namespace

CoupledPistonCase::CoupledPistonCase()
    : structure_(makeDefinition()),
      bridge_(structure_, rightPistonSurfaceStableId,
              makeCouplingNodes(), makeCouplingTriangles(),
              referenceFluidFaces(makeStepSettings().timeStepSeconds)),
      coupling_(bridge_.transfer()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader CoupledPistonCase::traceHeader() const {
    return {coupledPistonCaseChecksum, coupledPistonCaseSolverId};
}

viewer::DiagnosticFrame CoupledPistonCase::advance() {
    const double timeStep = stepSettings_.timeStepSeconds;
    const auto beforeStates = structure_.nodeStates();
    if (beforeStates.size() != 4) {
        throw std::logic_error("coupled piston topology changed unexpectedly");
    }
    const double startSpeed = beforeStates.front()
        .velocityMetersPerSecond.x;
    if (maximumVelocityDifference(beforeStates, startSpeed) > 1.0e-12) {
        throw std::runtime_error(
            "coupled piston no longer has rigid translational velocity");
    }

    const auto startFluid = solveFluidSample(startSpeed, timeStep);
    const double forceNewtons = rightSurface(startFluid)
        .pressureForceNewtons.x;
    const double predictedEndSpeed = startSpeed
        + forceNewtons * timeStep / structuralMassKilograms;
    const auto endFluid = solveFluidSample(predictedEndSpeed, timeStep);

    auto startKinematics = bridge_.transfer().captureKinematics(structure_);
    auto endKinematics = startKinematics;
    for (auto& node : endKinematics) {
        node.positionMeters.x += predictedEndSpeed * timeStep;
        node.velocityMetersPerSecond.x = predictedEndSpeed;
    }
    const auto startTransfer = bridge_.evaluate(
        startFluid, startKinematics);
    const auto endTransfer = bridge_.evaluate(
        endFluid, endKinematics);
    const std::array<double, 2> offsets{0.0, timeStep};
    const std::array transferSamples{
        startTransfer.transferResult(), endTransfer.transferResult()};
    const TimeIntegratedTransferResult integrated = coupling_.integrate(
        offsets, transferSamples);
    const StructureDiagnostics structureDiagnostics =
        coupling_.advanceStructure(structure_, integrated, stepSettings_);
    if (!structureDiagnostics.finite) {
        throw std::runtime_error(
            "coupled piston structural step produced non-finite diagnostics");
    }

    const auto acceptedStates = structure_.nodeStates();
    const double displacementResidual = acceptanceDisplacementResidual(
        beforeStates, acceptedStates, predictedEndSpeed, timeStep);
    if (displacementResidual > 1.0e-12) {
        throw std::runtime_error(
            "coupled piston structural acceptance missed its analytic predictor");
    }

    const auto& bridgeDiagnostics = endTransfer.diagnostics();
    const auto& integratedDiagnostics = integrated.diagnostics();
    viewer::StructureFrameContext context;
    context.sceneChecksum = coupledPistonCaseChecksum;
    context.solverCommit = coupledPistonCaseSolverId;
    context.timeStepSeconds = timeStep;
    context.couplingIteration = 1;
    context.couplingResiduals.displacementMetres = displacementResidual;
    context.couplingResiduals.tractionNewtons =
        bridgeDiagnostics.forceResidualNormNewtons;
    context.couplingResiduals.fluid =
        endFluid.projection.divergenceL2AfterPerSecond;
    context.couplingResiduals.structure =
        structureDiagnostics.maximumMembraneResidual;
    context.couplingResiduals.interfacePowerWatts =
        std::abs(bridgeDiagnostics.powerResidualWatts);
    context.conservation.fluidMassKilograms = fluidMassKilograms;
    context.conservation.totalMomentumNewtonSeconds = {
        structureDiagnostics.linearMomentumKgMetersPerSecond.x
            + fluidMassKilograms * predictedEndSpeed,
        structureDiagnostics.linearMomentumKgMetersPerSecond.y,
        structureDiagnostics.linearMomentumKgMetersPerSecond.z,
    };
    context.conservation.totalEnergyJoules =
        structureDiagnostics.kineticEnergyJoules
        + endFluid.projection.kineticEnergyAfterJoules;
    context.conservation.interfaceForceResidualNewtons =
        toViewer(bridgeDiagnostics.forceResidualNewtons);
    context.conservation.interfaceMomentResidualNewtonMetres =
        toViewer(endTransfer.transferResult().diagnostics()
                     .momentResidualNewtonMeters);
    context.conservation.interfacePowerResidualWatts =
        bridgeDiagnostics.powerResidualWatts;

    viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
        structure_, frameMapping_, context);
    appendPistonFields(
        frame, endFluid, bridgeDiagnostics, integratedDiagnostics);
    viewer::ProtocolError error;
    if (!viewer::validateFrame(frame, &error)) {
        throw std::runtime_error(
            "coupled piston produced an invalid diagnostic frame: "
            + error.message);
    }
    return frame;
}

const Structure& CoupledPistonCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& CoupledPistonCase::stepSettings() const noexcept {
    return stepSettings_;
}

StrongCoupledPistonCase::StrongCoupledPistonCase()
    : structure_(makeDefinition(strongStructuralMassKilograms)),
      grid_(pistonGrid()),
      fluidState_(initialFluidState(
          grid_, 0.0, makeStepSettings().timeStepSeconds)),
      bridge_(structure_, rightPistonSurfaceStableId,
              makeCouplingNodes(), makeCouplingTriangles(),
              fluidState_.diagnostics.faces),
      coupling_(bridge_.transfer()),
      stepSettings_(makeStepSettings()) {
    relaxationSettings_.initialRelaxation = 0.25;
    relaxationSettings_.minimumRelaxation = 0.02;
    relaxationSettings_.maximumRelaxation = 1.0;
    convergenceSettings_.minimumIterations = 3;
    convergenceSettings_.maximumIterations = 10;
    convergenceSettings_.absoluteDisplacementToleranceMetres = 1.0e-10;
    convergenceSettings_.relativeDisplacementTolerance = 1.0e-9;
    convergenceSettings_.absoluteVelocityToleranceMetersPerSecond = 1.0e-9;
    convergenceSettings_.relativeVelocityTolerance = 1.0e-9;
    convergenceSettings_.absoluteTractionToleranceNewtons = 1.0e-7;
    convergenceSettings_.relativeTractionTolerance = 1.0e-9;
    retrySettings_.maximumRetries = 2;
    retrySettings_.reductionFactor = 0.5;
    retrySettings_.minimumTimeStepSeconds =
        stepSettings_.timeStepSeconds * 0.25;
}

StrongCoupledPistonStepDiagnostics StrongCoupledPistonCase::advance() {
    const auto baselineKinematics =
        bridge_.transfer().captureKinematics(structure_);
    if (baselineKinematics.empty()) {
        throw std::logic_error(
            "strong piston has no coupling kinematics");
    }
    const double startSpeed =
        baselineKinematics.front().velocityMetersPerSecond.x;
    if (maximumVelocityDifference(
            structure_.nodeStates(), startSpeed) > 1.0e-12) {
        throw std::runtime_error(
            "strong piston no longer has rigid translational velocity");
    }
    const auto baselineTransfer = bridge_.evaluate(
        fluidState_.diagnostics, baselineKinematics);
    const std::array<double, 1> initialInterface{startSpeed};
    const std::uint64_t fingerprint = coupling_.surfaceFingerprint();
    StrongCouplingIteration iteration(
        fingerprint,
        initialInterface,
        relaxationSettings_,
        convergenceSettings_);
    StrongCouplingRollbackState rollback(
        fingerprint,
        std::move(structure_),
        grid_,
        std::move(fluidState_),
        std::move(iteration));
    StrongCouplingMacroStepState macroStep(
        std::move(rollback),
        stepSettings_.timeStepSeconds,
        retrySettings_);

    auto previousKinematics = baselineKinematics;
    ConservativeTransferResult previousTraction =
        baselineTransfer.transferResult();
    double activeTimeStepSeconds = 0.0;
    const StrongCouplingSolverCallback solve =
        [&](Structure& structure,
            const fluid::PeriodicCartesianGrid& grid,
            fluid::MovingInterfaceFluidState& fluidState,
            const std::span<const double> relaxedInterface,
            const double timeStepSeconds) {
            if (relaxedInterface.size() != 1) {
                throw std::logic_error(
                    "strong piston interface vector changed size");
            }
            if (timeStepSeconds != activeTimeStepSeconds) {
                activeTimeStepSeconds = timeStepSeconds;
                previousKinematics = baselineKinematics;
                previousTraction = baselineTransfer.transferResult();
            }
            const double endSpeed = relaxedInterface.front();
            fluidState.interfaces = fluid::FaceAlignedMovingInterface(
                grid, slabFaces(grid, endSpeed));
            fluidState.diagnostics =
                fluid::projectVelocityWithMovingInterfaces(
                    grid,
                    fluidState.velocityMetersPerSecond,
                    fluidState.pressurePascals,
                    fluidState.interfaces,
                    fluidProjectionSettings(timeStepSeconds));
            if (!fluidState.diagnostics.projection.converged
                || !fluidState.diagnostics.finite) {
                throw std::runtime_error(
                    "strong piston fluid projection was not accepted");
            }

            auto guessedKinematics = baselineKinematics;
            for (auto& node : guessedKinematics) {
                node.positionMeters.x += endSpeed * timeStepSeconds;
                node.velocityMetersPerSecond = {endSpeed, 0.0, 0.0};
            }
            const auto endTransfer = bridge_.evaluate(
                fluidState.diagnostics, guessedKinematics);
            const std::array<double, 2> offsets{
                0.0, timeStepSeconds};
            const std::array transferSamples{
                baselineTransfer.transferResult(),
                endTransfer.transferResult()};
            const auto integrated = coupling_.integrate(
                offsets, transferSamples);
            StructureStepSettings structuralSettings = stepSettings_;
            structuralSettings.timeStepSeconds = timeStepSeconds;
            const auto structureDiagnostics = coupling_.advanceStructure(
                structure, integrated, structuralSettings);
            if (!structureDiagnostics.finite) {
                throw std::runtime_error(
                    "strong piston structural solve was not accepted");
            }
            const auto currentKinematics =
                bridge_.transfer().captureKinematics(structure);
            const CouplingResidualNorms residuals =
                coupling_.measureResiduals(
                    baselineKinematics,
                    previousKinematics,
                    currentKinematics,
                    previousTraction,
                    endTransfer.transferResult());
            previousKinematics = currentKinematics;
            previousTraction = endTransfer.transferResult();

            const double solvedSpeed =
                currentKinematics.front().velocityMetersPerSecond.x;
            for (const auto& node : currentKinematics) {
                if (std::abs(node.velocityMetersPerSecond.x - solvedSpeed)
                        > 1.0e-12
                    || std::abs(node.velocityMetersPerSecond.y) > 1.0e-12
                    || std::abs(node.velocityMetersPerSecond.z) > 1.0e-12) {
                    throw std::runtime_error(
                        "strong piston structural solve lost rigid translation");
                }
            }
            return StrongCouplingSolverResult{{solvedSpeed}, residuals};
        };

    StrongCouplingMacroStepRunResult run;
    try {
        run = runStrongCouplingMacroStep(macroStep, solve);
    } catch (...) {
        structure_ = std::move(macroStep.rollbackState().structure());
        fluidState_ = std::move(macroStep.rollbackState().fluidState());
        throw;
    }
    structure_ = std::move(macroStep.rollbackState().structure());
    fluidState_ = std::move(macroStep.rollbackState().fluidState());

    const double acceptedSpeed =
        structure_.nodeStates().front().velocityMetersPerSecond.x;
    const double interfaceSpeed = fluidState_.interfaces.faces().front()
        .normalVelocityMetersPerSecond;
    StrongCoupledPistonStepDiagnostics result;
    result.coupling = std::move(run);
    result.startSpeedMetersPerSecond = startSpeed;
    result.acceptedSpeedMetersPerSecond = acceptedSpeed;
    result.acceptedInterfaceSpeedMetersPerSecond = interfaceSpeed;
    result.velocityClosureMetersPerSecond =
        std::abs(acceptedSpeed - interfaceSpeed);
    result.finite = std::isfinite(acceptedSpeed)
        && std::isfinite(interfaceSpeed)
        && std::isfinite(result.velocityClosureMetersPerSecond);
    return result;
}

const Structure& StrongCoupledPistonCase::structure() const noexcept {
    return structure_;
}

const fluid::MovingInterfaceFluidState&
StrongCoupledPistonCase::fluidState() const noexcept {
    return fluidState_;
}

const StructureStepSettings&
StrongCoupledPistonCase::stepSettings() const noexcept {
    return stepSettings_;
}

} // namespace simwing::fsi
