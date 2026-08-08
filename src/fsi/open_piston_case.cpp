#include "open_piston_case.h"

#include "open_piston_checkpoint_detail.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
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

fluid::MovingInterfaceProjectionDiagnostics initializeFluid(
    const fluid::PeriodicCartesianGrid& grid,
    fluid::MacVelocityField& velocity,
    fluid::CellScalarField& pressure,
    const fluid::FaceAlignedMovingInterface& interfaces,
    const double timeStepSeconds) {
    const auto diagnostics = fluid::projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces,
        projectionSettings(timeStepSeconds));
    if (!diagnostics.projection.converged || !diagnostics.finite) {
        throw std::runtime_error(
            "open piston failed to initialize its accepted fluid state");
    }
    return diagnostics;
}

PlanarFaceResolvedBridgeSettings movingBridgeSettings() {
    PlanarFaceResolvedBridgeSettings settings;
    settings.correspondenceMode =
        PlanarFaceCorrespondenceMode::RigidNormalTranslation;
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

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

StructureVector3 scale(const StructureVector3& value,
                       const double factor) {
    return {value.x * factor,
            value.y * factor,
            value.z * factor};
}

double length(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

StructureVector3 fluidMomentumNewtonSeconds(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    StructureVector3 result;
    for (const double value : velocity.xFaces()) {
        result.x += value;
    }
    for (const double value : velocity.yFaces()) {
        result.y += value;
    }
    for (const double value : velocity.zFaces()) {
        result.z += value;
    }
    return scale(
        result,
        fluidDensityKgPerCubicMeter
            * grid.cellVolumeCubicMeters());
}

double combinedTolerance(const double absoluteTolerance,
                         const double relativeTolerance,
                         const double firstScale,
                         const double secondScale) {
    return absoluteTolerance
        + relativeTolerance * std::max(firstScale, secondScale);
}

OpenPistonConservationDiagnostics makeConservationDiagnostics(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& startFluidVelocity,
    const fluid::MacVelocityField& endFluidVelocity,
    const StructureDiagnostics& startStructure,
    const StructureDiagnostics& endStructure,
    const TimeIntegratedTransferDiagnostics& pressureTransfer,
    const double actuatorImpulseNewtonSeconds,
    const double startSpeedMetersPerSecond,
    const double endSpeedMetersPerSecond) {
    OpenPistonConservationDiagnostics result;
    result.structureMomentumChangeNewtonSeconds = subtract(
        endStructure.linearMomentumKgMetersPerSecond,
        startStructure.linearMomentumKgMetersPerSecond);
    result.fluidMomentumChangeNewtonSeconds = subtract(
        fluidMomentumNewtonSeconds(grid, endFluidVelocity),
        fluidMomentumNewtonSeconds(grid, startFluidVelocity));
    result.pressureImpulseToStructureNewtonSeconds =
        pressureTransfer.integratedSurfaceImpulseNewtonSeconds;
    result.actuatorImpulseNewtonSeconds = {
        actuatorImpulseNewtonSeconds, 0.0, 0.0};
    result.structureMomentumResidualNewtonSeconds = subtract(
        result.structureMomentumChangeNewtonSeconds,
        add(result.pressureImpulseToStructureNewtonSeconds,
            result.actuatorImpulseNewtonSeconds));
    result.structureMomentumResidualNormNewtonSeconds = length(
        result.structureMomentumResidualNewtonSeconds);
    result.fluidMomentumResidualNewtonSeconds = add(
        result.fluidMomentumChangeNewtonSeconds,
        result.pressureImpulseToStructureNewtonSeconds);
    result.fluidMomentumResidualNormNewtonSeconds = length(
        result.fluidMomentumResidualNewtonSeconds);
    result.systemMomentumResidualNewtonSeconds = subtract(
        add(result.structureMomentumChangeNewtonSeconds,
            result.fluidMomentumChangeNewtonSeconds),
        result.actuatorImpulseNewtonSeconds);
    result.systemMomentumResidualNormNewtonSeconds = length(
        result.systemMomentumResidualNewtonSeconds);

    result.structureKineticEnergyChangeJoules =
        endStructure.kineticEnergyJoules
        - startStructure.kineticEnergyJoules;
    result.fluidKineticEnergyChangeJoules =
        fluid::kineticEnergyJoules(
            grid, endFluidVelocity, fluidDensityKgPerCubicMeter)
        - fluid::kineticEnergyJoules(
            grid, startFluidVelocity, fluidDensityKgPerCubicMeter);
    result.pressureWorkToStructureJoules =
        pressureTransfer.integratedSurfaceWorkJoules;
    result.actuatorWorkJoules = actuatorImpulseNewtonSeconds
        * 0.5 * (startSpeedMetersPerSecond
                 + endSpeedMetersPerSecond);
    result.structureEnergyResidualJoules =
        result.structureKineticEnergyChangeJoules
        - result.pressureWorkToStructureJoules
        - result.actuatorWorkJoules;
    result.fluidEnergyResidualJoules =
        result.fluidKineticEnergyChangeJoules
        + result.pressureWorkToStructureJoules;
    result.systemEnergyResidualJoules =
        result.structureKineticEnergyChangeJoules
        + result.fluidKineticEnergyChangeJoules
        - result.actuatorWorkJoules;

    result.finite = startStructure.finite && endStructure.finite
        && std::isfinite(result.structureMomentumResidualNormNewtonSeconds)
        && std::isfinite(result.fluidMomentumResidualNormNewtonSeconds)
        && std::isfinite(result.systemMomentumResidualNormNewtonSeconds)
        && std::isfinite(result.structureKineticEnergyChangeJoules)
        && std::isfinite(result.fluidKineticEnergyChangeJoules)
        && std::isfinite(result.pressureWorkToStructureJoules)
        && std::isfinite(result.actuatorWorkJoules)
        && std::isfinite(result.structureEnergyResidualJoules)
        && std::isfinite(result.fluidEnergyResidualJoules)
        && std::isfinite(result.systemEnergyResidualJoules);
    const double momentumTolerance = combinedTolerance(
        1.0e-8, 1.0e-11,
        length(result.structureMomentumChangeNewtonSeconds),
        std::max(
            length(result.fluidMomentumChangeNewtonSeconds),
            length(result.actuatorImpulseNewtonSeconds)));
    const double energyTolerance = combinedTolerance(
        2.0e-9, 1.0e-11,
        std::abs(result.structureKineticEnergyChangeJoules),
        std::max(
            std::abs(result.fluidKineticEnergyChangeJoules),
            std::abs(result.actuatorWorkJoules)));
    result.accepted = result.finite
        && result.structureMomentumResidualNormNewtonSeconds
            <= momentumTolerance
        && result.fluidMomentumResidualNormNewtonSeconds
            <= momentumTolerance
        && result.systemMomentumResidualNormNewtonSeconds
            <= momentumTolerance
        && std::abs(result.structureEnergyResidualJoules)
            <= energyTolerance
        && std::abs(result.fluidEnergyResidualJoules)
            <= energyTolerance
        && std::abs(result.systemEnergyResidualJoules)
            <= energyTolerance;
    return result;
}

viewer::Vec3d toViewer(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

void appendFields(
    viewer::DiagnosticFrame& frame,
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const fluid::PlanarControlVolumeDiagnostics& controlDiagnostics,
    const fluid::PlanarCutSurfacePressureDiagnostics& cutDiagnostics,
    const PlanarFaceResolvedBridgeDiagnostics& bridgeDiagnostics,
    const TimeIntegratedTransferDiagnostics& integratedDiagnostics,
    const OpenPistonConservationDiagnostics& conservation,
    const double actuatorImpulseNewtonSeconds,
    const std::uint64_t topologyRebaseCount,
    const double rebaseVolumeResidualCubicMeters,
    const double rebaseVelocityResidualMetersPerSecond) {
    frame.scalarFields.push_back({
        "interface.pressure_traction", "Pa",
        viewer::FieldAssociation::Triangle,
        {cutDiagnostics.pressureForceNewtons.x
             / cutDiagnostics.areaSquareMeters,
         cutDiagnostics.pressureForceNewtons.x
             / cutDiagnostics.areaSquareMeters}});
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
    frame.scalarFields.push_back({
        "actuator.step_work", "J", viewer::FieldAssociation::Global,
        {conservation.actuatorWorkJoules}});
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
    frame.scalarFields.push_back({
        "interface.grid_plane", "m",
        viewer::FieldAssociation::Global,
        {cutDiagnostics.gridPlaneCoordinateMeters}});
    frame.scalarFields.push_back({
        "interface.physical_plane", "m",
        viewer::FieldAssociation::Global,
        {cutDiagnostics.physicalPlaneCoordinateMeters}});
    frame.scalarFields.push_back({
        "interface.normal_translation", "m",
        viewer::FieldAssociation::Global,
        {bridgeDiagnostics.normalTranslationFromReferenceMeters}});
    frame.scalarFields.push_back({
        "interface.correspondence_residual", "m",
        viewer::FieldAssociation::Global,
        {bridgeDiagnostics.maximumRigidPositionResidualMeters}});
    frame.scalarFields.push_back({
        "interface.correspondence_velocity_residual", "m/s",
        viewer::FieldAssociation::Global,
        {bridgeDiagnostics
             .maximumRigidVelocityResidualMetersPerSecond}});
    frame.scalarFields.push_back({
        "fluid.cut_surface_periodic_residual", "m",
        viewer::FieldAssociation::Global,
        {cutDiagnostics.periodicPositionResidualMeters}});
    frame.scalarFields.push_back({
        "fluid.cut_surface_force_residual", "N",
        viewer::FieldAssociation::Global,
        {cutDiagnostics.forceResidualNormNewtons}});
    frame.scalarFields.push_back({
        "fluid.cut_surface_power_residual", "W",
        viewer::FieldAssociation::Global,
        {cutDiagnostics.powerResidualWatts}});
    frame.scalarFields.push_back({
        "conservation.structure_momentum_residual", "N*s",
        viewer::FieldAssociation::Global,
        {conservation.structureMomentumResidualNormNewtonSeconds}});
    frame.scalarFields.push_back({
        "conservation.fluid_momentum_residual", "N*s",
        viewer::FieldAssociation::Global,
        {conservation.fluidMomentumResidualNormNewtonSeconds}});
    frame.scalarFields.push_back({
        "conservation.system_momentum_residual", "N*s",
        viewer::FieldAssociation::Global,
        {conservation.systemMomentumResidualNormNewtonSeconds}});
    frame.scalarFields.push_back({
        "conservation.structure_energy_residual", "J",
        viewer::FieldAssociation::Global,
        {conservation.structureEnergyResidualJoules}});
    frame.scalarFields.push_back({
        "conservation.fluid_energy_residual", "J",
        viewer::FieldAssociation::Global,
        {conservation.fluidEnergyResidualJoules}});
    frame.scalarFields.push_back({
        "conservation.system_energy_residual", "J",
        viewer::FieldAssociation::Global,
        {conservation.systemEnergyResidualJoules}});
}

bool closeValue(const double first,
                const double second,
                const double absoluteTolerance = 1.0e-10,
                const double relativeTolerance = 1.0e-11) {
    return std::isfinite(first) && std::isfinite(second)
        && std::abs(first - second)
            <= absoluteTolerance
                + relativeTolerance
                    * std::max(std::abs(first), std::abs(second));
}

bool closeVector(const StructureVector3& first,
                 const StructureVector3& second,
                 const double tolerance = 1.0e-9) {
    return closeValue(first.x, second.x, tolerance)
        && closeValue(first.y, second.y, tolerance)
        && closeValue(first.z, second.z, tolerance);
}

bool closeVector(const fluid::Vector3& first,
                 const fluid::Vector3& second,
                 const double tolerance = 1.0e-9) {
    return closeValue(first.x, second.x, tolerance)
        && closeValue(first.y, second.y, tolerance)
        && closeValue(first.z, second.z, tolerance);
}

const char* invalidStoredOpenPistonDiagnostics(
    const OpenPistonCaseCheckpoint& checkpoint,
    const OpenPistonCaseCheckpoint::Detail& detail) {
    // These ledgers are restart state, not trusted commentary. Reconstruct
    // their canonical relationships before allowing restore to publish them.
    if (checkpoint.acceptedStepCount == 0) {
        const bool valid = checkpoint.topologyRebaseCount == 0
            && checkpoint.movingPlaneCoordinate == 6
            && checkpoint.surfaceOffsetMeters == 0.0
            && detail.controlVolumeDiagnostics
                == fluid::PlanarControlVolumeDiagnostics{}
            && detail.lastRebaseDiagnostics
                == fluid::PlanarControlVolumeRebaseDiagnostics{}
            && detail.cutSurfaceDiagnostics
                == fluid::PlanarCutSurfacePressureDiagnostics{}
            && detail.bridgeDiagnostics
                == PlanarFaceResolvedBridgeDiagnostics{}
            && detail.conservationDiagnostics
                == OpenPistonConservationDiagnostics{}
            && detail.lastRebaseVelocityResidualMetersPerSecond == 0.0;
        return valid ? nullptr : "initial diagnostic state";
    }

    constexpr std::uint64_t stepsPerCell = 1200;
    constexpr double cellSpacingMeters = 0.5;
    constexpr double initialGridPlaneMeters = 3.0;
    constexpr double initialReferenceVolumeCubicMeters = 12.0;
    constexpr double crossSectionAreaSquareMeters = 6.0;
    constexpr std::size_t faceCount = 6;
    const double timeStepSeconds = makeStepSettings().timeStepSeconds;
    const double stepTranslationMeters =
        coastSpeedMetersPerSecond * timeStepSeconds;
    const std::uint64_t expectedRebaseCount =
        checkpoint.acceptedStepCount / stepsPerCell;
    const std::uint64_t stepsInEpoch =
        checkpoint.acceptedStepCount % stepsPerCell;
    const bool crossingCheckpoint = stepsInEpoch == 0;
    const double expectedStoredOffset =
        static_cast<double>(stepsInEpoch) * stepTranslationMeters;
    if (checkpoint.topologyRebaseCount != expectedRebaseCount
        || !closeValue(checkpoint.surfaceOffsetMeters,
                       expectedStoredOffset, 2.0e-12)) {
        return "worker epoch state";
    }

    const std::uint64_t diagnosticEpoch = checkpoint.topologyRebaseCount
        - (crossingCheckpoint ? 1 : 0);
    const std::size_t diagnosticPlane = static_cast<std::size_t>(
        (6 + diagnosticEpoch) % 8);
    const double gridPlaneMeters =
        static_cast<double>(diagnosticPlane) * cellSpacingMeters;
    const double physicalPlaneMeters = initialGridPlaneMeters
        + static_cast<double>(checkpoint.topologyRebaseCount)
            * cellSpacingMeters
        + checkpoint.surfaceOffsetMeters;
    const double endOffsetMeters = crossingCheckpoint
        ? cellSpacingMeters : checkpoint.surfaceOffsetMeters;
    const double startOffsetMeters =
        endOffsetMeters - stepTranslationMeters;
    const double referenceVolumeCubicMeters =
        initialReferenceVolumeCubicMeters
        + static_cast<double>(diagnosticEpoch)
            * crossSectionAreaSquareMeters * cellSpacingMeters;

    const auto& control = detail.controlVolumeDiagnostics;
    if (control.version != fluid::planarMovingControlVolumeVersion
        || control.movingSurfaceStableId != pistonSurfaceStableId
        || control.fluidRegionStableId != connectedFluidRegionStableId
        || control.axis != fluid::GridFaceAxis::X
        || control.movingPlaneCoordinate != diagnosticPlane
        || control.openingPlaneCoordinate != 2
        || control.movingSurfaceFaceCount != faceCount
        || control.openingFaceCount != faceCount
        || !control.finite || !control.accepted
        || !closeValue(control.crossSectionAreaSquareMeters,
                       crossSectionAreaSquareMeters)
        || !closeValue(control.referenceVolumeCubicMeters,
                       referenceVolumeCubicMeters)
        || !closeValue(control.startCutCellVolumeCubicMeters,
                       crossSectionAreaSquareMeters * startOffsetMeters)
        || !closeValue(control.endCutCellVolumeCubicMeters,
                       crossSectionAreaSquareMeters * endOffsetMeters)
        || !closeValue(control.startCutCellVolumeFraction,
                       startOffsetMeters / cellSpacingMeters)
        || !closeValue(control.endCutCellVolumeFraction,
                       endOffsetMeters / cellSpacingMeters)
        || !closeValue(control.startVolumeCubicMeters,
                       referenceVolumeCubicMeters
                           + crossSectionAreaSquareMeters
                               * startOffsetMeters)
        || !closeValue(control.endVolumeCubicMeters,
                       referenceVolumeCubicMeters
                           + crossSectionAreaSquareMeters
                               * endOffsetMeters)
        || !closeValue(control.geometryVolumeChangeCubicMeters,
                       crossSectionAreaSquareMeters
                           * stepTranslationMeters)
        || !closeValue(control.surfaceSweptVolumeCubicMeters,
                       control.geometryVolumeChangeCubicMeters)
        || !closeValue(control.openingTransportVolumeCubicMeters,
                       control.geometryVolumeChangeCubicMeters)
        || !std::isfinite(control.surfaceGeometryResidualCubicMeters)
        || std::abs(control.surfaceGeometryResidualCubicMeters) > 1.0e-11
        || !std::isfinite(control.continuityResidualCubicMeters)
        || std::abs(control.continuityResidualCubicMeters) > 1.0e-11
        || !std::isfinite(
            control.maximumSurfaceVelocityErrorMetersPerSecond)
        || control.maximumSurfaceVelocityErrorMetersPerSecond < 0.0
        || control.maximumSurfaceVelocityErrorMetersPerSecond > 1.0e-11
        || !closeValue(control.rectangularSurfacePressureWorkJoules,
                       control.surfacePressurePowerWatts
                           * timeStepSeconds, 1.0e-9)) {
        return "moving control-volume ledger";
    }

    const auto& rebase = detail.lastRebaseDiagnostics;
    if (checkpoint.topologyRebaseCount == 0) {
        if (rebase != fluid::PlanarControlVolumeRebaseDiagnostics{}
            || detail.lastRebaseVelocityResidualMetersPerSecond != 0.0) {
            return "initial topology-rebase ledger";
        }
    } else {
        const std::size_t previousPlane = static_cast<std::size_t>(
            (6 + checkpoint.topologyRebaseCount - 1) % 8);
        const double rebasedVolume = initialReferenceVolumeCubicMeters
            + static_cast<double>(checkpoint.topologyRebaseCount)
                * crossSectionAreaSquareMeters * cellSpacingMeters;
        if (rebase.version != fluid::planarControlVolumeRebaseVersion
            || rebase.movingSurfaceStableId != pistonSurfaceStableId
            || rebase.fluidRegionStableId != connectedFluidRegionStableId
            || rebase.axis != fluid::GridFaceAxis::X
            || rebase.previousMovingPlaneCoordinate != previousPlane
            || rebase.rebasedMovingPlaneCoordinate
                != checkpoint.movingPlaneCoordinate
            || rebase.openingPlaneCoordinate != 2
            || !rebase.finite || !rebase.accepted
            || !closeValue(rebase.completedCellOffsetMeters,
                           cellSpacingMeters)
            || !closeValue(rebase.previousTerminalVolumeCubicMeters,
                           rebasedVolume)
            || !closeValue(rebase.rebasedReferenceVolumeCubicMeters,
                           rebasedVolume)
            || !std::isfinite(
                rebase.volumeContinuityResidualCubicMeters)
            || std::abs(rebase.volumeContinuityResidualCubicMeters)
                > 1.0e-11
            || detail.lastRebaseVelocityResidualMetersPerSecond < 0.0
            || detail.lastRebaseVelocityResidualMetersPerSecond
                > rebaseVelocityToleranceMetersPerSecond) {
            return "topology-rebase ledger";
        }
    }

    const auto& cut = detail.cutSurfaceDiagnostics;
    if (cut.version != fluid::planarCutSurfacePressureVersion
        || cut.sourceInterfaceVersion
            != fluid::faceAlignedMovingInterfaceVersion
        || cut.surfaceStableId != pistonSurfaceStableId
        || cut.fluidRegionStableId != connectedFluidRegionStableId
        || cut.axis != fluid::GridFaceAxis::X
        || cut.movingPlaneCoordinate != diagnosticPlane
        || cut.faceCount != faceCount || cut.faces.size() != faceCount
        || cut.kinematicsResampled || !cut.finite || !cut.accepted
        || !closeVector(cut.momentReferenceMeters, {})
        || !closeValue(cut.surfaceOffsetMeters, endOffsetMeters)
        || !closeValue(cut.gridPlaneCoordinateMeters, gridPlaneMeters)
        || !closeValue(cut.physicalPlaneCoordinateMeters,
                       physicalPlaneMeters)
        || !std::isfinite(cut.periodicPositionResidualMeters)
        || cut.periodicPositionResidualMeters < 0.0
        || cut.periodicPositionResidualMeters > 1.0e-11
        || !closeValue(cut.normalVelocityMetersPerSecond,
                       coastSpeedMetersPerSecond)
        || !std::isfinite(
            cut.maximumNormalVelocitySpreadMetersPerSecond)
        || cut.maximumNormalVelocitySpreadMetersPerSecond < 0.0
        || cut.maximumNormalVelocitySpreadMetersPerSecond > 1.0e-11
        || !closeValue(
            cut.reactionSourcePhysicalPlaneCoordinateMeters,
            physicalPlaneMeters)
        || !closeValue(
            cut.reactionSourceNormalVelocityMetersPerSecond,
            coastSpeedMetersPerSecond)) {
        return "cut-surface header";
    }

    fluid::Vector3 faceForce;
    fluid::Vector3 faceMoment;
    double faceArea = 0.0;
    double facePower = 0.0;
    std::array<bool, faceCount> seenFaces{};
    for (const auto& face : cut.faces) {
        if (face.surfaceStableId != pistonSurfaceStableId
            || face.minusRegionStableId != connectedFluidRegionStableId
            || face.plusRegionStableId != connectedFluidRegionStableId
            || face.axis != fluid::GridFaceAxis::X
            || face.i != diagnosticPlane || face.j >= 2 || face.k >= 3) {
            return "cut-surface face topology";
        }
        const std::size_t faceIndex = face.k * 2 + face.j;
        if (seenFaces[faceIndex]) return "duplicate cut-surface face";
        seenFaces[faceIndex] = true;
        const fluid::Vector3 expectedGridLower{
            gridPlaneMeters, static_cast<double>(face.j),
            static_cast<double>(face.k)};
        const fluid::Vector3 expectedGridUpper{
            gridPlaneMeters, static_cast<double>(face.j + 1),
            static_cast<double>(face.k + 1)};
        const fluid::Vector3 expectedPhysicalLower{
            physicalPlaneMeters, expectedGridLower.y,
            expectedGridLower.z};
        const fluid::Vector3 expectedPhysicalUpper{
            physicalPlaneMeters, expectedGridUpper.y,
            expectedGridUpper.z};
        const fluid::Vector3 expectedForce{
            face.pressureTractionPascals.x * face.areaSquareMeters,
            face.pressureTractionPascals.y * face.areaSquareMeters,
            face.pressureTractionPascals.z * face.areaSquareMeters};
        if (!closeVector(face.gridLowerCornerMeters, expectedGridLower)
            || !closeVector(face.gridUpperCornerMeters, expectedGridUpper)
            || !closeVector(face.physicalLowerCornerMeters,
                            expectedPhysicalLower)
            || !closeVector(face.physicalUpperCornerMeters,
                            expectedPhysicalUpper)
            || !closeValue(face.areaSquareMeters, 1.0)
            || !closeValue(face.normalVelocityMetersPerSecond,
                           coastSpeedMetersPerSecond)
            || !closeValue(face.pressureTractionPascals.y, 0.0)
            || !closeValue(face.pressureTractionPascals.z, 0.0)
            || !closeVector(face.pressureForceNewtons, expectedForce)
            || !closeValue(face.pressurePowerWatts,
                           face.pressureForceNewtons.x
                               * coastSpeedMetersPerSecond, 1.0e-9)) {
            return "cut-surface face geometry or traction";
        }
        faceArea += face.areaSquareMeters;
        faceForce.x += face.pressureForceNewtons.x;
        faceForce.y += face.pressureForceNewtons.y;
        faceForce.z += face.pressureForceNewtons.z;
        facePower += face.pressurePowerWatts;
        const double centerY = static_cast<double>(face.j) + 0.5;
        const double centerZ = static_cast<double>(face.k) + 0.5;
        faceMoment.x += centerY * face.pressureForceNewtons.z
            - centerZ * face.pressureForceNewtons.y;
        faceMoment.y += centerZ * face.pressureForceNewtons.x
            - physicalPlaneMeters * face.pressureForceNewtons.z;
        faceMoment.z += physicalPlaneMeters * face.pressureForceNewtons.y
            - centerY * face.pressureForceNewtons.x;
    }
    if (!std::ranges::all_of(seenFaces, [](const bool seen) { return seen; })
        || !closeValue(cut.areaSquareMeters, faceArea)
        || !closeValue(cut.sourceAreaSquareMeters, faceArea)
        || !closeValue(cut.areaResidualSquareMeters,
                       cut.areaSquareMeters - cut.sourceAreaSquareMeters)
        || !closeVector(cut.pressureForceNewtons, faceForce)
        || !closeVector(cut.sourcePressureForceNewtons, faceForce)
        || !closeVector(cut.forceResidualNewtons,
                        {cut.pressureForceNewtons.x
                             - cut.sourcePressureForceNewtons.x,
                         cut.pressureForceNewtons.y
                             - cut.sourcePressureForceNewtons.y,
                         cut.pressureForceNewtons.z
                             - cut.sourcePressureForceNewtons.z})
        || !closeValue(cut.forceResidualNormNewtons,
                       std::sqrt(
                           cut.forceResidualNewtons.x
                               * cut.forceResidualNewtons.x
                           + cut.forceResidualNewtons.y
                               * cut.forceResidualNewtons.y
                           + cut.forceResidualNewtons.z
                               * cut.forceResidualNewtons.z))
        || !closeVector(cut.pressureMomentNewtonMeters, faceMoment)
        || !closeValue(cut.pressurePowerWatts, facePower)
        || !closeValue(cut.sourcePressurePowerWatts, facePower)
        || !closeValue(cut.powerResidualWatts,
                       cut.pressurePowerWatts
                           - cut.sourcePressurePowerWatts)
        || !closeValue(control.surfacePressurePowerWatts,
                       cut.pressurePowerWatts, 1.0e-9)) {
        return "cut-surface aggregate ledger";
    }

    const auto& bridge = detail.bridgeDiagnostics;
    const auto bridgeSettings = movingBridgeSettings();
    const double bridgeAreaTolerance =
        bridgeSettings.absoluteAreaToleranceSquareMeters
        + bridgeSettings.relativeAreaTolerance * std::max(
            std::abs(bridge.fluidAreaSquareMeters),
            std::abs(bridge.referenceStructureAreaSquareMeters));
    const double bridgeForceTolerance =
        bridgeSettings.absoluteForceToleranceNewtons
        + bridgeSettings.relativeForceTolerance * std::max(
            length(bridge.fluidPressureForceNewtons),
            length(bridge.structureSurfaceForceNewtons));
    const double bridgeMomentTolerance =
        bridgeSettings.absoluteMomentToleranceNewtonMeters
        + bridgeSettings.relativeMomentTolerance * std::max(
            length(bridge.fluidPressureMomentNewtonMeters),
            length(bridge.structureSurfaceMomentNewtonMeters));
    const double bridgePowerTolerance =
        bridgeSettings.absolutePowerToleranceWatts
        + bridgeSettings.relativePowerTolerance * std::max(
            std::abs(bridge.fluidPressurePowerWatts),
            std::abs(bridge.structureSurfacePowerWatts));
    if (bridge.version != planarFaceResolvedBridgeVersion
        || bridge.fluidSurfaceStableId != pistonSurfaceStableId
        || bridge.fluidFaceCount != faceCount
        || bridge.structureTriangleCount != 2
        || bridge.overlapPatchCount != 10
        || bridge.correspondenceMode
            != PlanarFaceCorrespondenceMode::RigidNormalTranslation
        || !bridge.finite
        || !closeValue(bridge.fluidAreaSquareMeters,
                       crossSectionAreaSquareMeters)
        || !closeValue(bridge.referenceStructureAreaSquareMeters,
                       crossSectionAreaSquareMeters)
        || !closeValue(bridge.areaResidualSquareMeters,
                       bridge.referenceStructureAreaSquareMeters
                           - bridge.fluidAreaSquareMeters)
        || std::abs(bridge.areaResidualSquareMeters)
            > bridgeAreaTolerance
        || !closeVector(bridge.fluidPressureForceNewtons,
                       {cut.pressureForceNewtons.x,
                        cut.pressureForceNewtons.y,
                        cut.pressureForceNewtons.z})
        || !closeVector(bridge.forceResidualNewtons,
                       subtract(bridge.structureSurfaceForceNewtons,
                                bridge.fluidPressureForceNewtons))
        || !closeValue(bridge.forceResidualNormNewtons,
                       length(bridge.forceResidualNewtons))
        || bridge.forceResidualNormNewtons > bridgeForceTolerance
        || !closeVector(bridge.fluidPressureMomentNewtonMeters,
                       {cut.pressureMomentNewtonMeters.x,
                        cut.pressureMomentNewtonMeters.y,
                        cut.pressureMomentNewtonMeters.z})
        || !closeVector(bridge.momentResidualNewtonMeters,
                       subtract(bridge.structureSurfaceMomentNewtonMeters,
                                bridge.fluidPressureMomentNewtonMeters))
        || !closeValue(bridge.momentResidualNormNewtonMeters,
                       length(bridge.momentResidualNewtonMeters))
        || bridge.momentResidualNormNewtonMeters > bridgeMomentTolerance
        || !closeValue(bridge.fluidPressurePowerWatts,
                       cut.pressurePowerWatts, 1.0e-9)
        || !closeValue(bridge.powerResidualWatts,
                       bridge.structureSurfacePowerWatts
                           - bridge.fluidPressurePowerWatts, 1.0e-9)
        || std::abs(bridge.powerResidualWatts) > bridgePowerTolerance
        || !std::isfinite(bridge.maximumFacePowerResidualWatts)
        || bridge.maximumFacePowerResidualWatts < 0.0
        || bridge.maximumFacePowerResidualWatts > bridgePowerTolerance
        || !closeValue(bridge.gridPlaneCoordinateMeters, gridPlaneMeters)
        || !closeValue(bridge.physicalPlaneCoordinateMeters,
                       physicalPlaneMeters)
        || !closeValue(bridge.normalTranslationFromReferenceMeters,
                       physicalPlaneMeters - initialGridPlaneMeters)
        || !std::isfinite(bridge.maximumRigidPositionResidualMeters)
        || bridge.maximumRigidPositionResidualMeters < 0.0
        || bridge.maximumRigidPositionResidualMeters > 1.0e-10
        || !std::isfinite(
            bridge.maximumRigidVelocityResidualMetersPerSecond)
        || bridge.maximumRigidVelocityResidualMetersPerSecond < 0.0
        || bridge.maximumRigidVelocityResidualMetersPerSecond > 1.0e-10) {
        return "structure-fluid bridge ledger";
    }

    const auto& conservation = detail.conservationDiagnostics;
    const double startSpeed = checkpoint.acceptedStepCount == 1
        ? 0.0 : coastSpeedMetersPerSecond;
    const StructureVector3 expectedPressureImpulse{
        cut.pressureForceNewtons.x * timeStepSeconds,
        cut.pressureForceNewtons.y * timeStepSeconds,
        cut.pressureForceNewtons.z * timeStepSeconds};
    const StructureVector3 expectedStructureMomentum{
        structuralMassKilograms
            * (coastSpeedMetersPerSecond - startSpeed), 0.0, 0.0};
    const StructureVector3 expectedFluidMomentum =
        scale(expectedPressureImpulse, -1.0);
    const StructureVector3 expectedActuatorImpulse = subtract(
        expectedStructureMomentum, expectedPressureImpulse);
    const double averageSpeed =
        0.5 * (startSpeed + coastSpeedMetersPerSecond);
    const double expectedPressureWork =
        cut.pressureForceNewtons.x * averageSpeed * timeStepSeconds;
    const double expectedActuatorWork =
        expectedActuatorImpulse.x * averageSpeed;
    const double expectedStructureEnergy = 0.5
        * structuralMassKilograms
        * (coastSpeedMetersPerSecond * coastSpeedMetersPerSecond
           - startSpeed * startSpeed);
    const double expectedFluidEnergy = -expectedPressureWork;
    if (conservation.version != openPistonConservationVersion
        || !conservation.finite || !conservation.accepted
        || !closeVector(conservation.structureMomentumChangeNewtonSeconds,
                       expectedStructureMomentum, 1.0e-8)
        || !closeVector(conservation.fluidMomentumChangeNewtonSeconds,
                       expectedFluidMomentum, 1.0e-8)
        || !closeVector(
            conservation.pressureImpulseToStructureNewtonSeconds,
            expectedPressureImpulse, 1.0e-8)
        || !closeVector(conservation.actuatorImpulseNewtonSeconds,
                       expectedActuatorImpulse, 1.0e-8)
        || !closeVector(conservation.structureMomentumResidualNewtonSeconds,
                       subtract(expectedStructureMomentum,
                                add(expectedPressureImpulse,
                                    expectedActuatorImpulse)), 1.0e-8)
        || !closeValue(
            conservation.structureMomentumResidualNormNewtonSeconds,
            length(conservation.structureMomentumResidualNewtonSeconds),
            1.0e-8)
        || !closeVector(conservation.fluidMomentumResidualNewtonSeconds,
                       add(expectedFluidMomentum,
                           expectedPressureImpulse), 1.0e-8)
        || !closeValue(
            conservation.fluidMomentumResidualNormNewtonSeconds,
            length(conservation.fluidMomentumResidualNewtonSeconds),
            1.0e-8)
        || !closeVector(conservation.systemMomentumResidualNewtonSeconds,
                       subtract(add(expectedStructureMomentum,
                                    expectedFluidMomentum),
                                expectedActuatorImpulse), 1.0e-8)
        || !closeValue(
            conservation.systemMomentumResidualNormNewtonSeconds,
            length(conservation.systemMomentumResidualNewtonSeconds),
            1.0e-8)
        || !closeValue(conservation.structureKineticEnergyChangeJoules,
                       expectedStructureEnergy, 2.0e-9)
        || !closeValue(conservation.fluidKineticEnergyChangeJoules,
                       expectedFluidEnergy, 2.0e-9)
        || !closeValue(conservation.pressureWorkToStructureJoules,
                       expectedPressureWork, 2.0e-9)
        || !closeValue(conservation.actuatorWorkJoules,
                       expectedActuatorWork, 2.0e-9)
        || !closeValue(conservation.structureEnergyResidualJoules,
                       expectedStructureEnergy - expectedPressureWork
                           - expectedActuatorWork, 2.0e-9)
        || !closeValue(conservation.fluidEnergyResidualJoules,
                       expectedFluidEnergy + expectedPressureWork, 2.0e-9)
        || !closeValue(conservation.systemEnergyResidualJoules,
                       expectedStructureEnergy + expectedFluidEnergy
                           - expectedActuatorWork, 2.0e-9)) {
        return "coupled conservation ledger";
    }
    return nullptr;
}

} // namespace

OpenPistonCase::OpenPistonCase()
    : grid_(makeGrid()),
      interfaces_(makeInterfaces(grid_, 0.0)),
      fluidVelocity_(grid_),
      fluidPressure_(grid_),
      fluidDiagnostics_(initializeFluid(
          grid_, fluidVelocity_, fluidPressure_, interfaces_,
          makeStepSettings().timeStepSeconds)),
      structure_(makeDefinition()),
      bridge_(structure_, pistonSurfaceStableId,
              makeCouplingNodes(), makeCouplingTriangles(),
              fluidDiagnostics_.faces, movingBridgeSettings()),
      coupling_(bridge_.transfer()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()),
      controlVolume_(grid_, interfaces_,
                     pistonSurfaceStableId, 2) {}

viewer::TraceHeader OpenPistonCase::traceHeader() const {
    return {openPistonCaseChecksum, openPistonCaseSolverId};
}

viewer::DiagnosticFrame OpenPistonCase::advance() {
    const double timeStep = stepSettings_.timeStepSeconds;
    const auto beforeStates = structure_.nodeStates();
    const StructureDiagnostics beforeStructureDiagnostics =
        structure_.diagnostics();
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
    const double startPhysicalPlaneMeters =
        beforeStates.front().positionMeters.x;
    const double endPhysicalPlaneMeters =
        startPhysicalPlaneMeters + endSpeed * timeStep;
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
    auto acceptedInterfaces = endInterfaces;
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
        acceptedInterfaces = rebasedInterfaces;
    }

    const auto startKinematics = bridge_.transfer().captureKinematics(
        structure_);
    const auto endKinematics = predictedKinematics(
        bridge_.transfer(), structure_, endSpeed, timeStep);
    const auto endCutSurface =
        fluid::evaluatePlanarCutSurfacePressure(
            grid_, controlVolume_, endFluid,
            endOffset, endPhysicalPlaneMeters);
    const auto startCutSurface =
        fluid::resamplePlanarCutSurfaceReaction(
            grid_, controlVolume_, endCutSurface,
            surfaceOffsetMeters_, startPhysicalPlaneMeters,
            startSpeed);
    if (!startCutSurface.accepted || !endCutSurface.accepted) {
        throw std::runtime_error(
            "open piston physical cut-surface pressure was not accepted");
    }
    const auto startTransfer = bridge_.evaluateCutSurface(
        startCutSurface, startKinematics);
    const auto endTransfer = bridge_.evaluateCutSurface(
        endCutSurface, endKinematics);
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
    OpenPistonConservationDiagnostics conservation;
    try {
        conservation = makeConservationDiagnostics(
            grid_, fluidVelocity_, candidateVelocity,
            beforeStructureDiagnostics, structureDiagnostics,
            integratedDiagnostics, actuatorImpulse,
            startSpeed, endSpeed);
        if (!conservation.accepted) {
            throw std::runtime_error(
                "open piston coupled momentum or energy ledger did not close: "
                "structure momentum="
                + std::to_string(
                    conservation
                        .structureMomentumResidualNormNewtonSeconds)
                + ", fluid momentum="
                + std::to_string(
                    conservation.fluidMomentumResidualNormNewtonSeconds)
                + ", system momentum="
                + std::to_string(
                    conservation.systemMomentumResidualNormNewtonSeconds)
                + ", structure energy="
                + std::to_string(
                    conservation.structureEnergyResidualJoules)
                + ", fluid energy="
                + std::to_string(
                    conservation.fluidEnergyResidualJoules)
                + ", system energy="
                + std::to_string(
                    conservation.systemEnergyResidualJoules));
        }
    } catch (...) {
        structure_.restore(beforeActuator);
        throw;
    }
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
    const StructureVector3 endFluidMomentum =
        fluidMomentumNewtonSeconds(grid_, candidateVelocity);
    context.conservation.fluidMassKilograms = fluidMass;
    context.conservation.totalMomentumNewtonSeconds = toViewer(add(
        structureDiagnostics.linearMomentumKgMetersPerSecond,
        endFluidMomentum));
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
            frame, acceptedFluid, controlDiagnostics, endCutSurface,
            bridgeDiagnostics, integratedDiagnostics, conservation,
            actuatorImpulse,
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
    interfaces_ = std::move(acceptedInterfaces);
    fluidDiagnostics_ = acceptedFluid;
    controlVolumeDiagnostics_ = controlDiagnostics;
    cutSurfaceDiagnostics_ = endCutSurface;
    bridgeDiagnostics_ = bridgeDiagnostics;
    conservationDiagnostics_ = conservation;
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

OpenPistonCaseCheckpoint OpenPistonCase::checkpoint() const {
    OpenPistonCaseCheckpoint result;
    const StructureCheckpoint structureCheckpoint = structure_.checkpoint();
    result.acceptedStepCount = structureCheckpoint.acceptedStepCount;
    result.topologyRebaseCount = topologyRebaseCount_;
    result.movingPlaneCoordinate = controlVolume_.movingPlaneCoordinate();
    result.surfaceOffsetMeters = surfaceOffsetMeters_;
    auto detail = std::make_shared<OpenPistonCaseCheckpoint::Detail>(
        OpenPistonCaseCheckpoint::Detail{
            structureCheckpoint,
            fluid::checkpointMovingInterfaceFluidState(
                grid_, fluidVelocity_, fluidPressure_,
                interfaces_, fluidDiagnostics_),
            controlVolumeDiagnostics_,
            lastRebaseDiagnostics_,
            cutSurfaceDiagnostics_,
            bridgeDiagnostics_,
            conservationDiagnostics_,
            lastRebaseVelocityResidualMetersPerSecond_,
            topologyRebaseCount_,
            controlVolume_.movingPlaneCoordinate(),
            surfaceOffsetMeters_,
        });
    result.detail = std::move(detail);
    return result;
}

void OpenPistonCase::restore(
    const OpenPistonCaseCheckpoint& checkpointValue) {
    if (checkpointValue.version != openPistonCaseCheckpointVersion
        || checkpointValue.caseDefinitionFingerprint
            != openPistonCaseDefinitionFingerprint
        || !checkpointValue.detail
        || checkpointValue.acceptedStepCount
            != checkpointValue.detail->structure.acceptedStepCount
        || checkpointValue.topologyRebaseCount
            != checkpointValue.detail->topologyRebaseCount
        || checkpointValue.movingPlaneCoordinate
            != checkpointValue.detail->movingPlaneCoordinate
        || checkpointValue.surfaceOffsetMeters
            != checkpointValue.detail->surfaceOffsetMeters
        || checkpointValue.surfaceOffsetMeters < 0.0
        || checkpointValue.surfaceOffsetMeters
            >= controlVolume_.normalCellSpacingMeters()) {
        throw std::invalid_argument(
            "open piston checkpoint metadata is invalid");
    }

    auto candidateFluid = fluid::restoreMovingInterfaceFluidState(
        grid_, checkpointValue.detail->fluid);
    fluid::PlanarMovingControlVolume candidateControlVolume(
        grid_, candidateFluid.interfaces,
        pistonSurfaceStableId, 2);
    const std::size_t expectedPlane = static_cast<std::size_t>(
        (6 + checkpointValue.topologyRebaseCount)
        % grid_.cellCounts().x);
    const double cellSpacing =
        candidateControlVolume.normalCellSpacingMeters();
    const double expectedPhysicalPlane = 3.0
        + static_cast<double>(checkpointValue.topologyRebaseCount)
            * cellSpacing
        + checkpointValue.surfaceOffsetMeters;
    const auto& structureCheckpoint = checkpointValue.detail->structure;
    bool structureMatchesEpoch =
        structureCheckpoint.version == structureCheckpointVersion
        && structureCheckpoint.definitionFingerprint
            == structure_.definitionFingerprint()
        && structureCheckpoint.nodes.size()
            == structure_.definition().nodes.size()
        && std::abs(
            structureCheckpoint.simulationTimeSeconds
            - static_cast<double>(structureCheckpoint.acceptedStepCount)
                * stepSettings_.timeStepSeconds) <= 1.0e-12;
    double structuralSpeed = 0.0;
    if (!structureCheckpoint.nodes.empty()) {
        structuralSpeed = structureCheckpoint.nodes.front()
            .velocityMetersPerSecond.x;
    }
    for (const auto& node : structureCheckpoint.nodes) {
        structureMatchesEpoch = structureMatchesEpoch
            && std::abs(node.positionMeters.x - expectedPhysicalPlane)
                <= 2.0e-12
            && std::abs(
                node.velocityMetersPerSecond.x - structuralSpeed)
                <= 1.0e-12
            && std::abs(node.velocityMetersPerSecond.y) <= 1.0e-12
            && std::abs(node.velocityMetersPerSecond.z) <= 1.0e-12;
    }
    bool interfacesMatchEpoch = candidateFluid.interfaces.faceCount() == 6;
    for (const auto& face : candidateFluid.interfaces.faces()) {
        interfacesMatchEpoch = interfacesMatchEpoch
            && face.surfaceStableId == pistonSurfaceStableId
            && face.minusRegionStableId == connectedFluidRegionStableId
            && face.plusRegionStableId == connectedFluidRegionStableId
            && face.axis == fluid::GridFaceAxis::X
            && face.i == expectedPlane
            && std::abs(
                face.normalVelocityMetersPerSecond - structuralSpeed)
                <= 1.0e-12;
    }
    if (checkpointValue.movingPlaneCoordinate != expectedPlane
        || candidateControlVolume.movingPlaneCoordinate() != expectedPlane
        || candidateControlVolume.fluidRegionStableId()
            != connectedFluidRegionStableId
        || candidateControlVolume.openingPlaneCoordinate() != 2
        || checkpointValue.detail->lastRebaseVelocityResidualMetersPerSecond
            < 0.0
        || !std::isfinite(
            checkpointValue.detail
                ->lastRebaseVelocityResidualMetersPerSecond)
        || !structureMatchesEpoch || !interfacesMatchEpoch) {
        throw std::invalid_argument(
            "open piston checkpoint does not match its topology epoch");
    }
    if (const char* diagnosticError = invalidStoredOpenPistonDiagnostics(
            checkpointValue, *checkpointValue.detail)) {
        throw std::invalid_argument(
            std::string("open piston checkpoint diagnostics are inconsistent: ")
                + diagnosticError);
    }

    auto controlDiagnostics =
        checkpointValue.detail->controlVolumeDiagnostics;
    auto rebaseDiagnostics =
        checkpointValue.detail->lastRebaseDiagnostics;
    auto cutDiagnostics = checkpointValue.detail->cutSurfaceDiagnostics;
    auto bridgeDiagnostics = checkpointValue.detail->bridgeDiagnostics;
    auto conservationDiagnostics =
        checkpointValue.detail->conservationDiagnostics;
    const double rebaseVelocityResidual = checkpointValue.detail
        ->lastRebaseVelocityResidualMetersPerSecond;

    structure_.restore(structureCheckpoint);
    fluidVelocity_ = std::move(candidateFluid.velocityMetersPerSecond);
    fluidPressure_ = std::move(candidateFluid.pressurePascals);
    interfaces_ = std::move(candidateFluid.interfaces);
    fluidDiagnostics_ = std::move(candidateFluid.diagnostics);
    controlVolume_ = std::move(candidateControlVolume);
    controlVolumeDiagnostics_ = std::move(controlDiagnostics);
    lastRebaseDiagnostics_ = std::move(rebaseDiagnostics);
    cutSurfaceDiagnostics_ = std::move(cutDiagnostics);
    bridgeDiagnostics_ = std::move(bridgeDiagnostics);
    conservationDiagnostics_ = std::move(conservationDiagnostics);
    surfaceOffsetMeters_ = checkpointValue.surfaceOffsetMeters;
    lastRebaseVelocityResidualMetersPerSecond_ = rebaseVelocityResidual;
    topologyRebaseCount_ = checkpointValue.topologyRebaseCount;
}

std::uint64_t OpenPistonCase::acceptedStepCount() const noexcept {
    return structure_.acceptedStepCount();
}

double OpenPistonCase::simulationTimeSeconds() const noexcept {
    return structure_.simulationTimeSeconds();
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

const PlanarFaceResolvedBridgeDiagnostics&
OpenPistonCase::bridgeDiagnostics() const noexcept {
    return bridgeDiagnostics_;
}

const fluid::PlanarCutSurfacePressureDiagnostics&
OpenPistonCase::cutSurfaceDiagnostics() const noexcept {
    return cutSurfaceDiagnostics_;
}

const OpenPistonConservationDiagnostics&
OpenPistonCase::conservationDiagnostics() const noexcept {
    return conservationDiagnostics_;
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
