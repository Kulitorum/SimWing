#include "open_piston_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t pistonSurfaceStableId = 300;
constexpr std::uint64_t connectedFluidRegionStableId = 9;
constexpr double fluidDensityKgPerCubicMeter = 1.2;
constexpr double structuralMassKilograms = 6000.0;
constexpr double coastSpeedMetersPerSecond = 0.05;
constexpr double rebasePositionToleranceMeters = 1.0e-12;
constexpr double rebaseVelocityToleranceMetersPerSecond = 2.0e-12;

StructureDefinition makeDefinition() {
    StructureDefinition definition;
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
    mapping.triangles = {
        {1000, connectedFluidRegionStableId, connectedFluidRegionStableId},
        {2000, connectedFluidRegionStableId, connectedFluidRegionStableId},
    };
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

fluid::PeriodicCartesianGrid makeGrid() {
    return fluid::PeriodicCartesianGrid(
        {8, 2, 3}, {}, {4.0, 2.0, 3.0});
}

std::vector<fluid::GridFaceMovingInterface> makeInterfaceFaces(
    const fluid::PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond,
    const std::size_t movingPlaneCoordinate) {
    std::vector<fluid::GridFaceMovingInterface> faces;
    const auto counts = grid.cellCounts();
    faces.reserve(counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                pistonSurfaceStableId,
                connectedFluidRegionStableId,
                connectedFluidRegionStableId,
                fluid::GridFaceAxis::X,
                movingPlaneCoordinate,
                j,
                k,
                speedMetersPerSecond,
            });
        }
    }
    return faces;
}

fluid::FaceAlignedMovingInterface makeInterfaces(
    const fluid::PeriodicCartesianGrid& grid,
    const double speedMetersPerSecond,
    const std::size_t movingPlaneCoordinate = 6) {
    return fluid::FaceAlignedMovingInterface(
        grid, makeInterfaceFaces(
            grid, speedMetersPerSecond, movingPlaneCoordinate));
}

fluid::MovingInterfaceProjectionSettings projectionSettings(
    const double timeStepSeconds) {
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter =
        fluidDensityKgPerCubicMeter;
    settings.projection.timeStepSeconds = timeStepSeconds;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond =
        1.0e-12;
    return settings;
}

std::vector<CouplingNodeKinematics> predictedKinematics(
    const ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const double endSpeedMetersPerSecond,
    const double timeStepSeconds) {
    auto result = transfer.captureKinematics(structure);
    for (auto& node : result) {
        node.positionMeters.x +=
            endSpeedMetersPerSecond * timeStepSeconds;
        node.velocityMetersPerSecond = {
            endSpeedMetersPerSecond, 0.0, 0.0};
    }
    return result;
}

double maximumRigidVelocityError(
    const std::vector<StructureNodeState>& states,
    const double expectedSpeedMetersPerSecond) {
    double result = 0.0;
    for (const auto& state : states) {
        result = std::max({
            result,
            std::abs(state.velocityMetersPerSecond.x
                     - expectedSpeedMetersPerSecond),
            std::abs(state.velocityMetersPerSecond.y),
            std::abs(state.velocityMetersPerSecond.z),
        });
    }
    return result;
}

double maximumFluidVelocityDifference(
    const fluid::MacVelocityField& first,
    const fluid::MacVelocityField& second) {
    double result = 0.0;
    const auto accumulate = [&](const auto firstValues,
                                const auto secondValues) {
        for (std::size_t index = 0; index < firstValues.size(); ++index) {
            result = std::max(
                result, std::abs(firstValues[index] - secondValues[index]));
        }
    };
    accumulate(first.xFaces(), second.xFaces());
    accumulate(first.yFaces(), second.yFaces());
    accumulate(first.zFaces(), second.zFaces());
    return result;
}

double acceptanceResidual(
    const std::vector<StructureNodeState>& before,
    const std::vector<StructureNodeState>& after,
    const double expectedSpeedMetersPerSecond,
    const double timeStepSeconds) {
    if (before.size() != after.size()) {
        throw std::logic_error(
            "open piston acceptance changed the node count");
    }
    double result = maximumRigidVelocityError(
        after, expectedSpeedMetersPerSecond) * timeStepSeconds;
    for (std::size_t index = 0; index < before.size(); ++index) {
        result = std::max({
            result,
            std::abs(after[index].positionMeters.x
                     - before[index].positionMeters.x
                     - expectedSpeedMetersPerSecond * timeStepSeconds),
            std::abs(after[index].positionMeters.y
                     - before[index].positionMeters.y),
            std::abs(after[index].positionMeters.z
                     - before[index].positionMeters.z),
        });
    }
    return result;
}

viewer::Vec3d toViewer(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

void appendFields(
    viewer::DiagnosticFrame& frame,
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const fluid::PlanarControlVolumeDiagnostics& controlDiagnostics,
    const UniformFluidStructureBridgeDiagnostics& bridgeDiagnostics,
    const TimeIntegratedTransferDiagnostics& integratedDiagnostics,
    const double actuatorImpulseNewtonSeconds,
    const std::uint64_t topologyRebaseCount,
    const double rebaseVolumeResidualCubicMeters,
    const double rebaseVelocityResidualMetersPerSecond) {
    frame.scalarFields.push_back({
        "interface.pressure_traction", "Pa",
        viewer::FieldAssociation::Triangle,
        {bridgeDiagnostics.uniformPressureTractionPascals.x,
         bridgeDiagnostics.uniformPressureTractionPascals.x}});
    frame.scalarFields.push_back({
        "fluid.chamber_volume", "m^3", viewer::FieldAssociation::Global,
        {controlDiagnostics.endVolumeCubicMeters}});
    frame.scalarFields.push_back({
        "fluid.cut_cell_fraction", "1", viewer::FieldAssociation::Global,
        {controlDiagnostics.endCutCellVolumeFraction}});
    frame.scalarFields.push_back({
        "fluid.opening_transport", "m^3", viewer::FieldAssociation::Global,
        {controlDiagnostics.openingTransportVolumeCubicMeters}});
    frame.scalarFields.push_back({
        "fluid.gcl_residual", "m^3", viewer::FieldAssociation::Global,
        {controlDiagnostics.continuityResidualCubicMeters}});
    frame.scalarFields.push_back({
        "interface.step_work", "J", viewer::FieldAssociation::Global,
        {integratedDiagnostics.integratedSurfaceWorkJoules}});
    frame.scalarFields.push_back({
        "actuator.step_impulse", "N*s", viewer::FieldAssociation::Global,
        {actuatorImpulseNewtonSeconds}});
    frame.vectorFields.push_back({
        "interface.step_impulse", "N*s", viewer::FieldAssociation::Global,
        {toViewer(
            integratedDiagnostics.integratedSurfaceImpulseNewtonSeconds)}});
    frame.vectorFields.push_back({
        "interface.pressure_force", "N", viewer::FieldAssociation::Global,
        {toViewer(bridgeDiagnostics.fluidPressureForceNewtons)}});
    frame.scalarFields.push_back({
        "fluid.divergence_l2", "1/s", viewer::FieldAssociation::Global,
        {fluidDiagnostics.projection.divergenceL2AfterPerSecond}});
    frame.scalarFields.push_back({
        "fluid.topology_rebase_count", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(topologyRebaseCount)}});
    frame.scalarFields.push_back({
        "fluid.rebase_volume_residual", "m^3",
        viewer::FieldAssociation::Global,
        {rebaseVolumeResidualCubicMeters}});
    frame.scalarFields.push_back({
        "fluid.rebase_velocity_residual", "m/s",
        viewer::FieldAssociation::Global,
        {rebaseVelocityResidualMetersPerSecond}});
}

} // namespace

OpenPistonCase::OpenPistonCase()
    : grid_(makeGrid()),
      fluidVelocity_(grid_),
      fluidPressure_(grid_),
      structure_(makeDefinition()),
      bridge_(structure_, pistonSurfaceStableId,
              makeCouplingNodes(), makeCouplingTriangles()),
      coupling_(bridge_.transfer()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()),
      controlVolume_(grid_, makeInterfaces(grid_, 0.0),
                     pistonSurfaceStableId, 2) {
    const auto interfaces = makeInterfaces(grid_, 0.0);
    fluidDiagnostics_ = fluid::projectVelocityWithMovingInterfaces(
        grid_, fluidVelocity_, fluidPressure_, interfaces,
        projectionSettings(stepSettings_.timeStepSeconds));
    if (!fluidDiagnostics_.projection.converged
        || !fluidDiagnostics_.finite) {
        throw std::runtime_error(
            "open piston failed to initialize its accepted fluid state");
    }
}

viewer::TraceHeader OpenPistonCase::traceHeader() const {
    return {openPistonCaseChecksum, openPistonCaseSolverId};
}

viewer::DiagnosticFrame OpenPistonCase::advance() {
    const double timeStep = stepSettings_.timeStepSeconds;
    const auto beforeStates = structure_.nodeStates();
    if (beforeStates.size() != 4) {
        throw std::logic_error(
            "open piston structure topology changed unexpectedly");
    }
    const double startSpeed = beforeStates.front()
        .velocityMetersPerSecond.x;
    if (maximumRigidVelocityError(beforeStates, startSpeed) > 1.0e-12) {
        throw std::runtime_error(
            "open piston no longer has rigid translational velocity");
    }
    const double endSpeed = coastSpeedMetersPerSecond;
    const double rawEndOffset =
        surfaceOffsetMeters_ + endSpeed * timeStep;
    const double cellSpacing =
        controlVolume_.normalCellSpacingMeters();
    const bool endsAtCellBoundary =
        rawEndOffset >= cellSpacing - rebasePositionToleranceMeters;
    if (rawEndOffset > cellSpacing + rebasePositionToleranceMeters) {
        throw std::runtime_error(
            "open piston crossed a MAC face away from a macro-step boundary");
    }
    const double endOffset = endsAtCellBoundary
        ? cellSpacing : rawEndOffset;

    auto candidateVelocity = fluidVelocity_;
    auto candidatePressure = fluidPressure_;
    const auto endInterfaces = makeInterfaces(
        grid_, endSpeed, controlVolume_.movingPlaneCoordinate());
    const auto endFluid = fluid::projectVelocityWithMovingInterfaces(
        grid_, candidateVelocity, candidatePressure, endInterfaces,
        projectionSettings(timeStep));
    if (!endFluid.projection.converged || !endFluid.finite) {
        throw std::runtime_error(
            "open piston fluid projection was not accepted");
    }
    const auto controlDiagnostics = controlVolume_.evaluate(
        grid_, candidateVelocity, endFluid,
        {surfaceOffsetMeters_, endOffset, timeStep,
         endsAtCellBoundary});
    if (!controlDiagnostics.accepted) {
        throw std::runtime_error(
            "open piston geometric-conservation ledger did not close");
    }

    auto acceptedFluid = endFluid;
    std::optional<fluid::PlanarControlVolumeRebaseResult> rebase;
    double rebaseVolumeResidual = 0.0;
    double rebaseVelocityResidual = 0.0;
    if (endsAtCellBoundary) {
        const std::size_t rebasedPlane =
            (controlVolume_.movingPlaneCoordinate() + 1)
            % grid_.cellCounts().x;
        const auto rebasedInterfaces = makeInterfaces(
            grid_, endSpeed, rebasedPlane);
        rebase.emplace(fluid::rebasePlanarMovingControlVolume(
            grid_, controlVolume_, rebasedInterfaces,
            controlDiagnostics));
        if (!rebase->diagnostics.accepted) {
            throw std::runtime_error(
                "open piston topology-rebase volume ledger did not close");
        }

        auto rebasedVelocity = candidateVelocity;
        auto rebasedPressure = candidatePressure;
        const auto rebasedFluid =
            fluid::projectVelocityWithMovingInterfaces(
                grid_, rebasedVelocity, rebasedPressure,
                rebasedInterfaces, projectionSettings(timeStep));
        rebaseVelocityResidual = maximumFluidVelocityDifference(
            candidateVelocity, rebasedVelocity);
        if (!rebasedFluid.projection.converged
            || !rebasedFluid.finite
            || rebaseVelocityResidual
                > rebaseVelocityToleranceMetersPerSecond) {
            throw std::runtime_error(
                "open piston topology rebase changed its accepted "
                "fluid velocity");
        }
        rebaseVolumeResidual = rebase->diagnostics
            .volumeContinuityResidualCubicMeters;
        candidateVelocity = std::move(rebasedVelocity);
        candidatePressure = std::move(rebasedPressure);
        acceptedFluid = rebasedFluid;
    }

    const auto startKinematics = bridge_.transfer().captureKinematics(
        structure_);
    const auto endKinematics = predictedKinematics(
        bridge_.transfer(), structure_, endSpeed, timeStep);
    const auto startTransfer = bridge_.evaluate(
        fluidDiagnostics_, startKinematics);
    const auto endTransfer = bridge_.evaluate(
        endFluid, endKinematics);
    const std::array<double, 2> offsets{0.0, timeStep};
    const std::array samples{
        startTransfer.transferResult(), endTransfer.transferResult()};
    const auto integrated = coupling_.integrate(offsets, samples);

    const double desiredStructuralImpulse = structuralMassKilograms
        * (endSpeed - startSpeed);
    const double fluidImpulse = integrated.diagnostics()
        .integratedSurfaceImpulseNewtonSeconds.x;
    const double actuatorImpulse = desiredStructuralImpulse - fluidImpulse;
    const StructureCheckpoint beforeActuator = structure_.checkpoint();
    StructureDiagnostics structureDiagnostics;
    double displacementResidual = 0.0;
    try {
        for (std::size_t node = 0;
             node < structure_.definition().nodes.size(); ++node) {
            const double fraction =
                structure_.definition().nodes[node].massKg
                / structuralMassKilograms;
            structure_.addExternalForce(
                node,
                {fraction * actuatorImpulse / timeStep, 0.0, 0.0});
        }
        structureDiagnostics = coupling_.advanceStructure(
            structure_, integrated, stepSettings_);
        if (!structureDiagnostics.finite) {
            throw std::runtime_error(
                "open piston structural step produced non-finite diagnostics");
        }
        const auto acceptedStates = structure_.nodeStates();
        displacementResidual = acceptanceResidual(
            beforeStates, acceptedStates, endSpeed, timeStep);
        if (displacementResidual > 1.0e-12) {
            throw std::runtime_error(
                "open piston structural acceptance missed its prescribed drive");
        }
    } catch (...) {
        structure_.restore(beforeActuator);
        throw;
    }

    const auto& bridgeDiagnostics = endTransfer.diagnostics();
    const auto& integratedDiagnostics = integrated.diagnostics();
    viewer::StructureFrameContext context;
    context.sceneChecksum = openPistonCaseChecksum;
    context.solverCommit = openPistonCaseSolverId;
    context.timeStepSeconds = timeStep;
    context.couplingIteration = 1;
    context.couplingResiduals.displacementMetres = displacementResidual;
    context.couplingResiduals.tractionNewtons =
        bridgeDiagnostics.forceResidualNormNewtons;
    context.couplingResiduals.fluid =
        acceptedFluid.projection.divergenceL2AfterPerSecond;
    context.couplingResiduals.structure =
        structureDiagnostics.maximumMembraneResidual;
    context.couplingResiduals.interfacePowerWatts =
        std::abs(bridgeDiagnostics.powerResidualWatts);
    const double domainVolume = grid_.cellVolumeCubicMeters()
        * static_cast<double>(grid_.cellCount());
    const double fluidMass = fluidDensityKgPerCubicMeter * domainVolume;
    context.conservation.fluidMassKilograms = fluidMass;
    context.conservation.totalMomentumNewtonSeconds = {
        structureDiagnostics.linearMomentumKgMetersPerSecond.x
            + fluidMass * endSpeed,
        structureDiagnostics.linearMomentumKgMetersPerSecond.y,
        structureDiagnostics.linearMomentumKgMetersPerSecond.z,
    };
    context.conservation.totalEnergyJoules =
        structureDiagnostics.kineticEnergyJoules
        + acceptedFluid.projection.kineticEnergyAfterJoules;
    context.conservation.interfaceForceResidualNewtons =
        toViewer(bridgeDiagnostics.forceResidualNewtons);
    context.conservation.interfaceMomentResidualNewtonMetres =
        toViewer(endTransfer.transferResult().diagnostics()
                     .momentResidualNewtonMeters);
    context.conservation.interfacePowerResidualWatts =
        bridgeDiagnostics.powerResidualWatts;

    viewer::DiagnosticFrame frame;
    try {
        frame = viewer::buildStructureFrame(
            structure_, frameMapping_, context);
        appendFields(
            frame, acceptedFluid, controlDiagnostics, bridgeDiagnostics,
            integratedDiagnostics, actuatorImpulse,
            topologyRebaseCount_ + (endsAtCellBoundary ? 1 : 0),
            rebaseVolumeResidual, rebaseVelocityResidual);
        viewer::ProtocolError error;
        if (!viewer::validateFrame(frame, &error)) {
            throw std::runtime_error(
                "open piston produced an invalid diagnostic frame: "
                + error.message);
        }
    } catch (...) {
        structure_.restore(beforeActuator);
        throw;
    }

    fluidVelocity_ = std::move(candidateVelocity);
    fluidPressure_ = std::move(candidatePressure);
    fluidDiagnostics_ = acceptedFluid;
    controlVolumeDiagnostics_ = controlDiagnostics;
    if (rebase.has_value()) {
        controlVolume_ = std::move(rebase->controlVolume);
        lastRebaseDiagnostics_ = rebase->diagnostics;
        lastRebaseVelocityResidualMetersPerSecond_ =
            rebaseVelocityResidual;
        ++topologyRebaseCount_;
        surfaceOffsetMeters_ = 0.0;
    } else {
        surfaceOffsetMeters_ = endOffset;
    }
    return frame;
}

const Structure& OpenPistonCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& OpenPistonCase::stepSettings() const noexcept {
    return stepSettings_;
}

const fluid::PlanarControlVolumeDiagnostics&
OpenPistonCase::controlVolumeDiagnostics() const noexcept {
    return controlVolumeDiagnostics_;
}

const fluid::PlanarControlVolumeRebaseDiagnostics&
OpenPistonCase::lastRebaseDiagnostics() const noexcept {
    return lastRebaseDiagnostics_;
}

double OpenPistonCase::surfaceOffsetMeters() const noexcept {
    return surfaceOffsetMeters_;
}

std::size_t OpenPistonCase::movingPlaneCoordinate() const noexcept {
    return controlVolume_.movingPlaneCoordinate();
}

std::uint64_t OpenPistonCase::topologyRebaseCount() const noexcept {
    return topologyRebaseCount_;
}

double OpenPistonCase::lastRebaseVelocityResidualMetersPerSecond()
    const noexcept {
    return lastRebaseVelocityResidualMetersPerSecond_;
}

} // namespace simwing::fsi
