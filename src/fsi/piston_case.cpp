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

StructureDefinition makeDefinition() {
    StructureDefinition definition;
    // The two-triangle barycentric load shares are 2:1:2:1. Mass follows the
    // same ratio so uniform pressure preserves rigid piston translation.
    definition.nodes = {
        {{3.0, 0.0, 0.0}, 2000.0, false},
        {{3.0, 2.0, 0.0}, 1000.0, false},
        {{3.0, 2.0, 3.0}, 2000.0, false},
        {{3.0, 0.0, 3.0}, 1000.0, false},
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
    const fluid::FaceAlignedMovingInterface interfaces(
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

} // namespace simwing::fsi
