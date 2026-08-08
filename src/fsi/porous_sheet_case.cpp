#include "porous_sheet_case.h"
#include "porous_sheet_checkpoint_detail.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t porousSheetSurfaceStableId = 300;
constexpr std::uint64_t pumpSurfaceStableId = 400;
constexpr std::uint64_t minusRegionStableId = 10;
constexpr std::uint64_t plusRegionStableId = 11;
constexpr std::size_t pumpFaceCoordinate = 5;
constexpr double initialSheetPositionMeters = 1.45;
constexpr double fluidDensityKgPerCubicMeter = 1.2;
constexpr double sheetMassKilograms = 60.0;
constexpr double pumpPressureJumpPascals = 1.5;
constexpr double linearResistancePascalSecondsPerMeter = 10.0;

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

StructureVector3 scale(const StructureVector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double length(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

viewer::Vec3d toViewer(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

StructureVector3 toStructure(const fluid::Vector3& value) {
    return {value.x, value.y, value.z};
}

fluid::PeriodicCartesianGrid makeGrid() {
    return {{8, 2, 3}, {}, {4.0, 2.0, 3.0}};
}

double fluidMassKilograms(const fluid::PeriodicCartesianGrid& grid) {
    return fluidDensityKgPerCubicMeter
        * static_cast<double>(grid.cellCount())
        * grid.cellVolumeCubicMeters();
}

double sheetAreaSquareMeters(const fluid::PeriodicCartesianGrid& grid) {
    const auto spacing = grid.cellSpacingMeters();
    const auto counts = grid.cellCounts();
    return spacing.y * spacing.z
        * static_cast<double>(counts.y * counts.z);
}

StructureDefinition makeDefinition() {
    StructureDefinition definition;
    definition.nodes = {
        {{initialSheetPositionMeters, 0.0, 0.0}, 20.0, false},
        {{initialSheetPositionMeters, 2.0, 0.0}, 10.0, false},
        {{initialSheetPositionMeters, 2.0, 3.0}, 20.0, false},
        {{initialSheetPositionMeters, 0.0, 3.0}, 10.0, false},
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
        {1000, minusRegionStableId, plusRegionStableId},
        {2000, minusRegionStableId, plusRegionStableId},
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

fluid::PorousProjectionSettings makeProjectionSettings() {
    fluid::PorousProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter =
        fluidDensityKgPerCubicMeter;
    settings.projection.timeStepSeconds =
        makeStepSettings().timeStepSeconds;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    settings.iteration.constitutiveEvaluation =
        fluid::PorousConstitutiveEvaluation::Midpoint;
    settings.iteration.absoluteNormalVelocityToleranceMetersPerSecond =
        1.0e-12;
    settings.iteration.relativeNormalVelocityTolerance = 1.0e-12;
    settings.iteration.absolutePressureJumpTolerancePascals = 1.0e-11;
    settings.iteration.relativePressureJumpTolerance = 1.0e-12;
    settings.iteration.relaxation = 0.5;
    settings.iteration.maximumNonlinearIterations = 100;
    return settings;
}

double crossingFraction(const fluid::PeriodicCartesianGrid& grid,
                        const std::size_t faceCoordinate,
                        const double sheetPositionMeters) {
    const double faceCenter = grid.xFaceCenterMeters(
        faceCoordinate, 0, 0).x;
    const double fraction = 0.5
        + (sheetPositionMeters - faceCenter)
            / grid.cellSpacingMeters().x;
    if (!std::isfinite(fraction)
        || !(fraction > 0.0) || !(fraction < 1.0)) {
        throw std::runtime_error(
            "coupled porous sheet reached its fixed-topology segment boundary");
    }
    return fraction;
}

struct PorousTopologySelection {
    std::size_t faceCoordinate =
        coupledPorousSheetInitialFaceCoordinate;
    bool rebased = false;
};

PorousTopologySelection selectPorousTopology(
    const fluid::PeriodicCartesianGrid& grid,
    const std::size_t currentFaceCoordinate,
    const double sheetPositionMeters) {
    const double faceCenter = grid.xFaceCenterMeters(
        currentFaceCoordinate, 0, 0).x;
    const double fraction = 0.5
        + (sheetPositionMeters - faceCenter)
            / grid.cellSpacingMeters().x;
    if (!std::isfinite(fraction)) {
        throw std::runtime_error(
            "coupled porous sheet topology position is non-finite");
    }
    if (fraction > 0.0 && fraction < 1.0) {
        return {currentFaceCoordinate, false};
    }
    if (!(fraction >= 1.0) || !(fraction < 2.0)) {
        throw std::runtime_error(
            "coupled porous sheet crossed more than one topology segment");
    }
    const std::size_t rebasedCoordinate =
        (currentFaceCoordinate + 1) % grid.cellCounts().x;
    if (rebasedCoordinate == pumpFaceCoordinate) {
        throw std::runtime_error(
            "coupled porous sheet reached the pump-surface topology");
    }
    static_cast<void>(crossingFraction(
        grid, rebasedCoordinate, sheetPositionMeters));
    return {rebasedCoordinate, true};
}

std::vector<fluid::PorousGridFaceCrossing> makePorousSheet(
    const fluid::PeriodicCartesianGrid& grid,
    const std::size_t faceCoordinate,
    const double sheetPositionMeters,
    const double sheetVelocityMetersPerSecond) {
    const double fraction = crossingFraction(
        grid, faceCoordinate, sheetPositionMeters);
    const auto counts = grid.cellCounts();
    std::vector<fluid::PorousGridFaceCrossing> result;
    result.reserve(counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                porousSheetSurfaceStableId,
                minusRegionStableId,
                plusRegionStableId,
                fluid::GridFaceAxis::X,
                faceCoordinate,
                j,
                k,
                fraction,
                sheetVelocityMetersPerSecond,
                {linearResistancePascalSecondsPerMeter, 0.0},
            });
        }
    }
    return result;
}

fluid::SharpPressureJumpField makePump(
    const fluid::PeriodicCartesianGrid& grid) {
    const auto counts = grid.cellCounts();
    std::vector<fluid::GridFacePressureJump> result;
    result.reserve(counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                pumpSurfaceStableId,
                plusRegionStableId,
                minusRegionStableId,
                fluid::GridFaceAxis::X,
                pumpFaceCoordinate,
                j,
                k,
                pumpPressureJumpPascals,
                0.5,
            });
        }
    }
    return {grid, std::move(result)};
}

std::vector<fluid::PorousFaceTractionDiagnostics> referenceFaces(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PorousProjectionSettings& settings) {
    fluid::MacVelocityField velocity(grid);
    fluid::CellScalarField pressure(grid);
    const auto porous = makePorousSheet(
        grid, coupledPorousSheetInitialFaceCoordinate,
        initialSheetPositionMeters, 0.0);
    const auto pump = makePump(grid);
    const auto projected = fluid::projectVelocityWithPorousInterfaces(
        grid, velocity, pressure, porous, pump, settings);
    const auto traction = fluid::evaluatePorousSurfaceTraction(
        grid, projected, settings.projection.timeStepSeconds);
    if (!traction.accepted) {
        throw std::runtime_error(
            "coupled porous sheet reference projection was not accepted");
    }
    return traction.faces;
}

PlanarFaceResolvedBridgeSettings makeBridgeSettings() {
    PlanarFaceResolvedBridgeSettings settings;
    settings.correspondenceMode =
        PlanarFaceCorrespondenceMode::RigidNormalTranslation;
    return settings;
}

double rigidSheetState(const std::vector<StructureNodeState>& states,
                       const bool velocity,
                       double& maximumResidual) {
    if (states.size() != 4) {
        throw std::logic_error(
            "coupled porous sheet topology changed unexpectedly");
    }
    const double value = velocity
        ? states.front().velocityMetersPerSecond.x
        : states.front().positionMeters.x;
    maximumResidual = 0.0;
    for (const auto& state : states) {
        const auto& vector = velocity
            ? state.velocityMetersPerSecond : state.positionMeters;
        maximumResidual = std::max(
            maximumResidual, std::abs(vector.x - value));
        if (velocity) {
            maximumResidual = std::max(
                maximumResidual, std::abs(vector.y));
            maximumResidual = std::max(
                maximumResidual, std::abs(vector.z));
        }
    }
    return value;
}

double uniformFluidVelocity(
    const fluid::MacVelocityField& velocity,
    double& maximumResidual) {
    if (velocity.xFaces().empty()) {
        throw std::logic_error("coupled porous sheet fluid field is empty");
    }
    const double value = velocity.xFaces().front();
    maximumResidual = 0.0;
    for (const double sample : velocity.xFaces()) {
        maximumResidual = std::max(
            maximumResidual, std::abs(sample - value));
    }
    for (const double sample : velocity.yFaces()) {
        maximumResidual = std::max(maximumResidual, std::abs(sample));
    }
    for (const double sample : velocity.zFaces()) {
        maximumResidual = std::max(maximumResidual, std::abs(sample));
    }
    return value;
}

std::vector<CouplingNodeKinematics> midpointKinematics(
    const ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const double positionMeters,
    const double velocityMetersPerSecond) {
    auto result = transfer.captureKinematics(structure);
    for (auto& node : result) {
        node.positionMeters.x = positionMeters;
        node.velocityMetersPerSecond = {
            velocityMetersPerSecond, 0.0, 0.0};
    }
    return result;
}

void appendFrameFields(
    viewer::DiagnosticFrame& frame,
    const CoupledPorousSheetStepDiagnostics& diagnostics) {
    const double jump = diagnostics.porousTraction.faces.empty()
        ? 0.0
        : diagnostics.porousTraction.faces.front().pressureJumpPascals;
    const double relativeVelocity = diagnostics.porousTraction.faces.empty()
        ? 0.0
        : diagnostics.porousTraction.faces.front()
            .relativeNormalVelocityMetersPerSecond;
    frame.scalarFields.push_back({
        "porous.pressure_jump", "Pa",
        viewer::FieldAssociation::Triangle, {jump, jump}});
    frame.scalarFields.push_back({
        "porous.relative_velocity", "m/s",
        viewer::FieldAssociation::Global, {relativeVelocity}});
    frame.scalarFields.push_back({
        "porous.dissipated_energy", "J",
        viewer::FieldAssociation::Global,
        {diagnostics.porousDissipatedEnergyJoules}});
    frame.scalarFields.push_back({
        "pump.step_work", "J",
        viewer::FieldAssociation::Global,
        {diagnostics.pumpWorkJoules}});
    frame.scalarFields.push_back({
        "coupled.step_energy_residual", "J",
        viewer::FieldAssociation::Global,
        {diagnostics.energyResidualJoules}});
    frame.scalarFields.push_back({
        "fluid.topology_rebase_count", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(diagnostics.topologyRebaseCount)}});
    frame.scalarFields.push_back({
        "porous.face_coordinate", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(diagnostics.porousFaceCoordinate)}});
    frame.vectorFields.push_back({
        "porous.sheet_impulse", "N*s",
        viewer::FieldAssociation::Global,
        {toViewer(diagnostics.porousImpulseOnSheetNewtonSeconds)}});
    frame.vectorFields.push_back({
        "pump.fluid_impulse", "N*s",
        viewer::FieldAssociation::Global,
        {toViewer(diagnostics.pumpImpulseNewtonSeconds)}});
}

} // namespace

CoupledPorousSheetCase::CoupledPorousSheetCase()
    : grid_(makeGrid()),
      velocity_(grid_),
      pressure_(grid_),
      projectionSettings_(makeProjectionSettings()),
      structure_(makeDefinition()),
      bridge_(
          structure_, porousSheetSurfaceStableId,
          makeCouplingNodes(), makeCouplingTriangles(),
          referenceFaces(grid_, projectionSettings_),
          makeBridgeSettings()),
      coupling_(bridge_.transfer()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()) {
    if (projectionSettings_.projection.timeStepSeconds
        != stepSettings_.timeStepSeconds) {
        throw std::logic_error(
            "coupled porous sheet fluid and structure time steps disagree");
    }
}

viewer::TraceHeader CoupledPorousSheetCase::traceHeader() const {
    return {
        coupledPorousSheetCaseChecksum,
        coupledPorousSheetCaseSolverId,
    };
}

viewer::DiagnosticFrame CoupledPorousSheetCase::advance() {
    const StructureCheckpoint structureBefore = structure_.checkpoint();
    try {
        const auto statesBefore = structure_.nodeStates();
        const StructureDiagnostics structureDiagnosticsBefore =
            structure_.diagnostics();
        double sheetVelocityResidual = 0.0;
        double sheetPositionResidual = 0.0;
        const double sheetVelocityBefore = rigidSheetState(
            statesBefore, true, sheetVelocityResidual);
        const double sheetPositionBefore = rigidSheetState(
            statesBefore, false, sheetPositionResidual);
        double fluidUniformityBefore = 0.0;
        const double fluidVelocityBefore = uniformFluidVelocity(
            velocity_, fluidUniformityBefore);
        if (sheetVelocityResidual > 1.0e-12
            || sheetPositionResidual > 1.0e-12
            || fluidUniformityBefore > 1.0e-12) {
            throw std::runtime_error(
                "coupled porous sheet state is not rigid and uniform");
        }

        const double timeStep = stepSettings_.timeStepSeconds;
        const double fluidMass = fluidMassKilograms(grid_);
        const double area = sheetAreaSquareMeters(grid_);
        const double relativeVelocityBefore =
            fluidVelocityBefore - sheetVelocityBefore;
        const double denominator = 1.0
            + 0.5 * timeStep * area
                * linearResistancePascalSecondsPerMeter
                * (1.0 / fluidMass + 1.0 / sheetMassKilograms);
        const double relativeVelocityAtConstitutiveTime =
            (relativeVelocityBefore
             + 0.5 * timeStep * area
                 * pumpPressureJumpPascals / fluidMass)
            / denominator;
        const double predictedSheetForce = area
            * linearResistancePascalSecondsPerMeter
            * relativeVelocityAtConstitutiveTime;
        const double sheetVelocityAtConstitutiveTime =
            sheetVelocityBefore
            + 0.5 * timeStep * predictedSheetForce
                / sheetMassKilograms;
        const double predictedSheetVelocityAfter =
            sheetVelocityBefore
            + timeStep * predictedSheetForce / sheetMassKilograms;
        const double sheetPositionAtConstitutiveTime =
            sheetPositionBefore
            + 0.5 * timeStep * predictedSheetVelocityAfter;
        const PorousTopologySelection topology = selectPorousTopology(
            grid_, porousFaceCoordinate_,
            sheetPositionAtConstitutiveTime);

        fluid::MacVelocityField candidateVelocity = velocity_;
        fluid::CellScalarField candidatePressure = pressure_;
        const auto porous = makePorousSheet(
            grid_, topology.faceCoordinate,
            sheetPositionAtConstitutiveTime,
            sheetVelocityAtConstitutiveTime);
        const auto pump = makePump(grid_);
        const fluid::PorousProjectionDiagnostics projection =
            fluid::projectVelocityWithPorousInterfaces(
                grid_, candidateVelocity, candidatePressure,
                porous, pump, projectionSettings_);
        if (!projection.accepted || !projection.finite) {
            throw std::runtime_error(
                "coupled porous sheet fluid projection was not accepted");
        }
        const fluid::PorousSurfaceTractionDiagnostics traction =
            fluid::evaluatePorousSurfaceTraction(
                grid_, projection, timeStep);
        const auto kinematics = midpointKinematics(
            bridge_.transfer(), structure_,
            sheetPositionAtConstitutiveTime,
            sheetVelocityAtConstitutiveTime);
        const PorousFaceResolvedTransferResult mapped =
            bridge_.evaluateMovingPorousSurface(
                traction, kinematics,
                sheetPositionAtConstitutiveTime);
        const std::array<double, 2> sampleOffsets{0.0, timeStep};
        const std::array transferSamples{
            mapped.transferResult(), mapped.transferResult()};
        const TimeIntegratedTransferResult integrated =
            coupling_.integrate(sampleOffsets, transferSamples);
        const StructureDiagnostics structureDiagnosticsAfter =
            coupling_.advanceStructure(
                structure_, integrated, stepSettings_);
        if (!structureDiagnosticsAfter.finite) {
            throw std::runtime_error(
                "coupled porous sheet structural step was not accepted");
        }

        const auto statesAfter = structure_.nodeStates();
        double sheetVelocityAfterResidual = 0.0;
        double sheetPositionAfterResidual = 0.0;
        const double sheetVelocityAfter = rigidSheetState(
            statesAfter, true, sheetVelocityAfterResidual);
        const double sheetPositionAfter = rigidSheetState(
            statesAfter, false, sheetPositionAfterResidual);
        double fluidUniformityAfter = 0.0;
        const double fluidVelocityAfter = uniformFluidVelocity(
            candidateVelocity, fluidUniformityAfter);
        const double fluidVelocityAtConstitutiveTime =
            0.5 * (fluidVelocityBefore + fluidVelocityAfter);

        CoupledPorousSheetStepDiagnostics candidate;
        candidate.acceptedStepCount = structure_.acceptedStepCount();
        candidate.simulationTimeSeconds =
            structure_.simulationTimeSeconds();
        candidate.sheetPositionBeforeMeters = sheetPositionBefore;
        candidate.sheetPositionAtConstitutiveTimeMeters =
            sheetPositionAtConstitutiveTime;
        candidate.sheetPositionAfterMeters = sheetPositionAfter;
        candidate.fluidVelocityBeforeMetersPerSecond = fluidVelocityBefore;
        candidate.fluidVelocityAtConstitutiveTimeMetersPerSecond =
            fluidVelocityAtConstitutiveTime;
        candidate.fluidVelocityAfterMetersPerSecond = fluidVelocityAfter;
        candidate.sheetVelocityBeforeMetersPerSecond = sheetVelocityBefore;
        candidate.sheetVelocityAtConstitutiveTimeMetersPerSecond =
            sheetVelocityAtConstitutiveTime;
        candidate.sheetVelocityAfterMetersPerSecond = sheetVelocityAfter;
        candidate.actualFluidImpulseNewtonSeconds = {
            fluidMass * (fluidVelocityAfter - fluidVelocityBefore),
            0.0,
            0.0,
        };
        candidate.pressureJumpImpulseOnFluidNewtonSeconds = toStructure(
            projection.totalPressureJumpImpulseOnFluidNewtonSeconds);
        candidate.actualSheetImpulseNewtonSeconds = subtract(
            structureDiagnosticsAfter.linearMomentumKgMetersPerSecond,
            structureDiagnosticsBefore.linearMomentumKgMetersPerSecond);
        candidate.porousImpulseOnSheetNewtonSeconds =
            mapped.diagnostics().pressureImpulseOnSurfaceNewtonSeconds;
        candidate.pumpImpulseNewtonSeconds = {
            area * pumpPressureJumpPascals * timeStep, 0.0, 0.0};
        candidate.momentumResidualNewtonSeconds = subtract(
            add(candidate.actualFluidImpulseNewtonSeconds,
                candidate.actualSheetImpulseNewtonSeconds),
            candidate.pumpImpulseNewtonSeconds);
        candidate.momentumResidualNormNewtonSeconds = length(
            candidate.momentumResidualNewtonSeconds);
        candidate.fluidKineticEnergyChangeJoules =
            fluid::kineticEnergyJoules(
                grid_, candidateVelocity, fluidDensityKgPerCubicMeter)
            - fluid::kineticEnergyJoules(
                grid_, velocity_, fluidDensityKgPerCubicMeter);
        candidate.sheetKineticEnergyChangeJoules =
            structureDiagnosticsAfter.kineticEnergyJoules
            - structureDiagnosticsBefore.kineticEnergyJoules;
        candidate.pumpWorkJoules =
            projection.totalPressureJumpWorkToFluidJoules
            - traction.totalPressureWorkToFluidJoules;
        candidate.porousDissipatedEnergyJoules =
            traction.totalDissipatedEnergyJoules;
        candidate.energyResidualJoules =
            candidate.fluidKineticEnergyChangeJoules
            + candidate.sheetKineticEnergyChangeJoules
            + candidate.porousDissipatedEnergyJoules
            - candidate.pumpWorkJoules;
        candidate.maximumFluidUniformityResidualMetersPerSecond =
            std::max(fluidUniformityBefore, fluidUniformityAfter);
        candidate.maximumSheetRigidMotionResidualMeters = std::max({
            sheetVelocityResidual * timeStep,
            sheetPositionResidual,
            sheetVelocityAfterResidual * timeStep,
            sheetPositionAfterResidual,
            std::abs(sheetVelocityAfter
                     - predictedSheetVelocityAfter) * timeStep,
            std::abs(sheetPositionAfter
                     - (sheetPositionBefore
                        + timeStep * predictedSheetVelocityAfter)),
        });
        candidate.topologyRebaseCount = topologyRebaseCount_
            + (topology.rebased ? 1 : 0);
        candidate.porousFaceCoordinate = topology.faceCoordinate;
        candidate.topologyRebasedThisStep = topology.rebased;
        candidate.fluidProjection = projection;
        candidate.porousTraction = traction;
        candidate.bridge = mapped.diagnostics();
        candidate.transfer = integrated.diagnostics();
        const StructureVector3 fluidImpulseResidual = subtract(
            candidate.actualFluidImpulseNewtonSeconds,
            candidate.pressureJumpImpulseOnFluidNewtonSeconds);
        const StructureVector3 sheetImpulseResidual = subtract(
            candidate.actualSheetImpulseNewtonSeconds,
            candidate.porousImpulseOnSheetNewtonSeconds);
        const double energyScale =
            std::abs(candidate.fluidKineticEnergyChangeJoules)
            + std::abs(candidate.sheetKineticEnergyChangeJoules)
            + std::abs(candidate.porousDissipatedEnergyJoules)
            + std::abs(candidate.pumpWorkJoules);
        candidate.finite =
            std::isfinite(candidate.sheetPositionBeforeMeters)
            && std::isfinite(
                candidate.sheetPositionAtConstitutiveTimeMeters)
            && std::isfinite(candidate.sheetPositionAfterMeters)
            && std::isfinite(candidate.fluidVelocityBeforeMetersPerSecond)
            && std::isfinite(
                candidate.fluidVelocityAtConstitutiveTimeMetersPerSecond)
            && std::isfinite(candidate.fluidVelocityAfterMetersPerSecond)
            && std::isfinite(candidate.sheetVelocityBeforeMetersPerSecond)
            && std::isfinite(
                candidate.sheetVelocityAtConstitutiveTimeMetersPerSecond)
            && std::isfinite(candidate.sheetVelocityAfterMetersPerSecond)
            && finite(candidate.actualFluidImpulseNewtonSeconds)
            && finite(candidate.pressureJumpImpulseOnFluidNewtonSeconds)
            && finite(candidate.actualSheetImpulseNewtonSeconds)
            && finite(candidate.porousImpulseOnSheetNewtonSeconds)
            && finite(candidate.pumpImpulseNewtonSeconds)
            && finite(candidate.momentumResidualNewtonSeconds)
            && std::isfinite(candidate.momentumResidualNormNewtonSeconds)
            && std::isfinite(candidate.fluidKineticEnergyChangeJoules)
            && std::isfinite(candidate.sheetKineticEnergyChangeJoules)
            && std::isfinite(candidate.pumpWorkJoules)
            && std::isfinite(candidate.porousDissipatedEnergyJoules)
            && std::isfinite(candidate.energyResidualJoules)
            && std::isfinite(
                candidate.maximumFluidUniformityResidualMetersPerSecond)
            && std::isfinite(
                candidate.maximumSheetRigidMotionResidualMeters)
            && projection.finite && traction.finite
            && mapped.diagnostics().finite
            && integrated.diagnostics().finite;
        candidate.accepted = candidate.finite
            && projection.accepted && traction.accepted
            && mapped.diagnostics().accepted
            && length(fluidImpulseResidual) <= 2.0e-10
            && length(sheetImpulseResidual) <= 2.0e-10
            && candidate.momentumResidualNormNewtonSeconds <= 3.0e-10
            && std::abs(candidate.energyResidualJoules)
                <= 2.0e-10 + 1.0e-11 * energyScale
            && candidate.maximumFluidUniformityResidualMetersPerSecond
                <= 2.0e-12
            && candidate.maximumSheetRigidMotionResidualMeters <= 2.0e-12;
        if (!candidate.accepted) {
            throw std::runtime_error(
                "coupled porous sheet momentum, energy, or rigid-motion ledger did not close");
        }

        viewer::StructureFrameContext context;
        context.sceneChecksum = coupledPorousSheetCaseChecksum;
        context.solverCommit = coupledPorousSheetCaseSolverId;
        context.timeStepSeconds = timeStep;
        context.couplingIteration = 1;
        context.couplingResiduals.displacementMetres =
            candidate.maximumSheetRigidMotionResidualMeters;
        context.couplingResiduals.tractionNewtons =
            candidate.bridge.mapping.forceResidualNormNewtons;
        context.couplingResiduals.fluid =
            projection.projection.divergenceL2AfterPerSecond;
        context.couplingResiduals.structure =
            structureDiagnosticsAfter.maximumMembraneResidual;
        context.couplingResiduals.interfacePowerWatts =
            std::abs(candidate.bridge.mapping.powerResidualWatts);
        context.conservation.fluidMassKilograms = fluidMass;
        context.conservation.totalMomentumNewtonSeconds = {
            structureDiagnosticsAfter.linearMomentumKgMetersPerSecond.x
                + fluidMass * fluidVelocityAfter,
            structureDiagnosticsAfter.linearMomentumKgMetersPerSecond.y,
            structureDiagnosticsAfter.linearMomentumKgMetersPerSecond.z,
        };
        context.conservation.totalEnergyJoules =
            structureDiagnosticsAfter.kineticEnergyJoules
            + fluid::kineticEnergyJoules(
                grid_, candidateVelocity, fluidDensityKgPerCubicMeter);
        context.conservation.interfaceForceResidualNewtons =
            toViewer(candidate.bridge.mapping.forceResidualNewtons);
        context.conservation.interfaceMomentResidualNewtonMetres =
            toViewer(candidate.bridge.mapping.momentResidualNewtonMeters);
        context.conservation.interfacePowerResidualWatts =
            candidate.bridge.mapping.powerResidualWatts;
        viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
            structure_, frameMapping_, context);
        appendFrameFields(frame, candidate);
        viewer::ProtocolError error;
        if (!viewer::validateFrame(frame, &error)) {
            throw std::runtime_error(
                "coupled porous sheet produced an invalid frame: "
                + error.message);
        }

        velocity_ = std::move(candidateVelocity);
        pressure_ = std::move(candidatePressure);
        porousFaceCoordinate_ = topology.faceCoordinate;
        topologyRebaseCount_ = candidate.topologyRebaseCount;
        diagnostics_ = std::move(candidate);
        return frame;
    } catch (...) {
        structure_.restore(structureBefore);
        throw;
    }
}

CoupledPorousSheetCheckpoint
CoupledPorousSheetCase::checkpoint() const {
    CoupledPorousSheetCheckpoint result;
    result.acceptedStepCount = structure_.acceptedStepCount();
    result.simulationTimeSeconds = structure_.simulationTimeSeconds();
    result.topologyRebaseCount = topologyRebaseCount_;
    result.porousFaceCoordinate = porousFaceCoordinate_;
    result.detail = std::make_shared<
        CoupledPorousSheetCheckpoint::Detail>(
            structure_.checkpoint(), velocity_, pressure_, diagnostics_);
    return result;
}

void CoupledPorousSheetCase::restore(
    const CoupledPorousSheetCheckpoint& checkpointValue) {
    if (checkpointValue.version != coupledPorousSheetCheckpointVersion
        || checkpointValue.caseFingerprint
            != coupledPorousSheetCaseFingerprint
        || !std::isfinite(checkpointValue.simulationTimeSeconds)
        || checkpointValue.simulationTimeSeconds < 0.0
        || checkpointValue.porousFaceCoordinate
            >= grid_.cellCounts().x
        || checkpointValue.porousFaceCoordinate == pumpFaceCoordinate
        || checkpointValue.topologyRebaseCount > 1
        || checkpointValue.porousFaceCoordinate
            != coupledPorousSheetInitialFaceCoordinate
                + checkpointValue.topologyRebaseCount
        || !checkpointValue.detail) {
        throw std::invalid_argument(
            "coupled porous sheet checkpoint metadata is invalid");
    }
    const auto& detail = *checkpointValue.detail;
    if (!detail.velocity.matches(grid_)
        || !detail.pressure.matches(grid_)
        || !fluid::isFinite(detail.velocity)
        || !fluid::isFinite(detail.pressure)
        || detail.structure.acceptedStepCount
            != checkpointValue.acceptedStepCount
        || detail.structure.simulationTimeSeconds
            != checkpointValue.simulationTimeSeconds) {
        throw std::invalid_argument(
            "coupled porous sheet checkpoint fields are invalid");
    }

    Structure candidateStructure(makeDefinition());
    candidateStructure.restore(detail.structure);
    const auto candidateStates = candidateStructure.nodeStates();
    double positionResidual = 0.0;
    double velocityResidual = 0.0;
    const double sheetPosition = rigidSheetState(
        candidateStates, false, positionResidual);
    const double sheetVelocity = rigidSheetState(
        candidateStates, true, velocityResidual);
    double fluidUniformity = 0.0;
    const double fluidVelocity = uniformFluidVelocity(
        detail.velocity, fluidUniformity);
    // The accepted fluid/traction topology is evaluated at the constitutive
    // midpoint. A fast-moving sheet endpoint may already lie beyond that
    // segment while the accepted step is still valid; restoring such a
    // terminal safe point must reproduce the next topology rejection.
    const double topologyPosition =
        checkpointValue.acceptedStepCount == 0
        ? sheetPosition
        : detail.diagnostics.sheetPositionAtConstitutiveTimeMeters;
    bool canonicalTopologyPosition = true;
    try {
        static_cast<void>(crossingFraction(
            grid_, checkpointValue.porousFaceCoordinate,
            topologyPosition));
    } catch (const std::exception&) {
        canonicalTopologyPosition = false;
    }
    const double expectedTime = static_cast<double>(
        checkpointValue.acceptedStepCount)
        * stepSettings_.timeStepSeconds;
    if (positionResidual > 1.0e-12
        || velocityResidual > 1.0e-12
        || fluidUniformity > 2.0e-12
        || !canonicalTopologyPosition
        || std::abs(checkpointValue.simulationTimeSeconds - expectedTime)
            > 2.0e-14) {
        throw std::invalid_argument(
            "coupled porous sheet checkpoint state is not a canonical epoch");
    }

    if (checkpointValue.acceptedStepCount == 0) {
        if (detail.diagnostics
                != CoupledPorousSheetStepDiagnostics{}
            || sheetPosition != initialSheetPositionMeters
            || sheetVelocity != 0.0 || fluidVelocity != 0.0
            || checkpointValue.topologyRebaseCount != 0
            || checkpointValue.porousFaceCoordinate
                != coupledPorousSheetInitialFaceCoordinate) {
            throw std::invalid_argument(
                "coupled porous sheet initial checkpoint is inconsistent");
        }
    } else {
        const auto& diagnostics = detail.diagnostics;
        const auto structureDiagnostics = candidateStructure.diagnostics();
        const double fluidMass = fluidMassKilograms(grid_);
        const double expectedSystemMomentum =
            sheetAreaSquareMeters(grid_)
            * pumpPressureJumpPascals
            * checkpointValue.simulationTimeSeconds;
        const double actualSystemMomentum =
            structureDiagnostics.linearMomentumKgMetersPerSecond.x
            + fluidMass * fluidVelocity;
        const double fieldKineticEnergy = fluid::kineticEnergyJoules(
            grid_, detail.velocity, fluidDensityKgPerCubicMeter);
        if (!diagnostics.accepted || !diagnostics.finite
            || diagnostics.version
                != coupledPorousSheetDiagnosticsVersion
            || diagnostics.acceptedStepCount
                != checkpointValue.acceptedStepCount
            || diagnostics.simulationTimeSeconds
                != checkpointValue.simulationTimeSeconds
            || diagnostics.sheetPositionAfterMeters != sheetPosition
            || diagnostics.sheetVelocityAfterMetersPerSecond
                != sheetVelocity
            || diagnostics.fluidVelocityAfterMetersPerSecond
                != fluidVelocity
            || diagnostics.topologyRebaseCount
                != checkpointValue.topologyRebaseCount
            || diagnostics.porousFaceCoordinate
                != checkpointValue.porousFaceCoordinate
            || !diagnostics.fluidProjection.accepted
            || !diagnostics.porousTraction.accepted
            || !diagnostics.bridge.accepted
            || !diagnostics.transfer.finite
            || diagnostics.bridge.mapping.gridPlaneCoordinateMeters
                != diagnostics.sheetPositionAtConstitutiveTimeMeters
            || std::ranges::any_of(
                diagnostics.porousTraction.faces,
                [&](const auto& face) {
                    return face.surfaceStableId
                            == porousSheetSurfaceStableId
                        && face.i
                            != checkpointValue.porousFaceCoordinate;
                })
            || std::abs(
                diagnostics.fluidProjection.projection
                    .kineticEnergyAfterJoules
                - fieldKineticEnergy) > 2.0e-10
            || std::abs(actualSystemMomentum - expectedSystemMomentum)
                > 2.0e-8) {
            throw std::invalid_argument(
                "coupled porous sheet checkpoint diagnostics are inconsistent");
        }
    }

    structure_.restore(detail.structure);
    velocity_ = detail.velocity;
    pressure_ = detail.pressure;
    diagnostics_ = detail.diagnostics;
    topologyRebaseCount_ = checkpointValue.topologyRebaseCount;
    porousFaceCoordinate_ = checkpointValue.porousFaceCoordinate;
}

const Structure& CoupledPorousSheetCase::structure() const noexcept {
    return structure_;
}

const fluid::PeriodicCartesianGrid&
CoupledPorousSheetCase::grid() const noexcept {
    return grid_;
}

const fluid::MacVelocityField&
CoupledPorousSheetCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField&
CoupledPorousSheetCase::pressure() const noexcept {
    return pressure_;
}

const StructureStepSettings&
CoupledPorousSheetCase::stepSettings() const noexcept {
    return stepSettings_;
}

const CoupledPorousSheetStepDiagnostics&
CoupledPorousSheetCase::diagnostics() const noexcept {
    return diagnostics_;
}

std::uint64_t
CoupledPorousSheetCase::acceptedStepCount() const noexcept {
    return structure_.acceptedStepCount();
}

double CoupledPorousSheetCase::simulationTimeSeconds() const noexcept {
    return structure_.simulationTimeSeconds();
}

std::uint64_t
CoupledPorousSheetCase::topologyRebaseCount() const noexcept {
    return topologyRebaseCount_;
}

std::size_t CoupledPorousSheetCase::porousFaceCoordinate() const noexcept {
    return porousFaceCoordinate_;
}

} // namespace simwing::fsi
