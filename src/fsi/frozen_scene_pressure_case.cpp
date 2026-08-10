#include "frozen_scene_pressure_case.h"

#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_cell_volume.h"
#include "scene_fluid_mimetic_pressure_audit.h"
#include "scene_fluid_mimetic_pressure_sampling.h"
#include "scene_fluid_opening_cap.h"
#include "scene_fluid_opening_face_crossing.h"
#include "scene_fluid_opening_flux.h"
#include "scene_fluid_opening_patch.h"
#include "scene_fluid_opening_quadrature.h"
#include "scene_fluid_pressure_control_volume.h"
#include "scene_fluid_pressure_face_link.h"
#include "scene_fluid_region_connectivity.h"
#include "scene_fluid_surface.h"
#include "scene_fluid_surface_transfer.h"
#include "scene_structure.h"
#include "structure_frame.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

void validateSettings(const FrozenScenePressureCaseSettings& settings) {
    const auto counts = settings.cellCounts;
    if (counts.x < 2 || counts.y < 2 || counts.z < 2
        || !std::isfinite(settings.domainPaddingMeters)
        || !(settings.domainPaddingMeters > 0.0)
        || (settings.useExplicitDomain
            && (!std::isfinite(settings.lowerMeters.x)
                || !std::isfinite(settings.lowerMeters.y)
                || !std::isfinite(settings.lowerMeters.z)
                || !std::isfinite(settings.upperMeters.x)
                || !std::isfinite(settings.upperMeters.y)
                || !std::isfinite(settings.upperMeters.z)
                || !(settings.upperMeters.x > settings.lowerMeters.x)
                || !(settings.upperMeters.y > settings.lowerMeters.y)
                || !(settings.upperMeters.z > settings.lowerMeters.z)))
        || !std::isfinite(settings.backgroundWindMetersPerSecond.x)
        || !std::isfinite(settings.backgroundWindMetersPerSecond.y)
        || !std::isfinite(settings.backgroundWindMetersPerSecond.z)
        || !std::isfinite(
            settings.diagnosticPerturbationSpeedMetersPerSecond)
        || settings.diagnosticPerturbationSpeedMetersPerSecond < 0.0
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument(
            "frozen scene pressure settings are invalid");
    }
}

std::string assemblyError(const SceneStructureAssembly& assembly) {
    if (assembly.diagnostics.empty()) {
        return "scene has no supported structural assembly";
    }
    return assembly.diagnostics.front().message;
}

std::string surfaceError(const SceneFluidSurfaceAssembly& surface) {
    if (surface.diagnostics.empty()) {
        return "scene has no supported fluid surface";
    }
    return surface.diagnostics.front().message;
}

fluid::PeriodicCartesianGrid makeGrid(
    const SceneFluidSurfaceState& state,
    const FrozenScenePressureCaseSettings& settings) {
    if (state.vertices.empty()) {
        throw std::invalid_argument(
            "frozen scene pressure surface is empty");
    }
    if (settings.useExplicitDomain) {
        return {
            settings.cellCounts, settings.lowerMeters,
            settings.upperMeters,
        };
    }
    fluid::Vector3 lower{
        state.vertices.front().positionMeters.x,
        state.vertices.front().positionMeters.y,
        state.vertices.front().positionMeters.z,
    };
    fluid::Vector3 upper = lower;
    for (const auto& vertex : state.vertices) {
        lower.x = std::min(lower.x, vertex.positionMeters.x);
        lower.y = std::min(lower.y, vertex.positionMeters.y);
        lower.z = std::min(lower.z, vertex.positionMeters.z);
        upper.x = std::max(upper.x, vertex.positionMeters.x);
        upper.y = std::max(upper.y, vertex.positionMeters.y);
        upper.z = std::max(upper.z, vertex.positionMeters.z);
    }
    lower.x -= settings.domainPaddingMeters;
    lower.y -= settings.domainPaddingMeters;
    lower.z -= settings.domainPaddingMeters;
    upper.x += settings.domainPaddingMeters;
    upper.y += settings.domainPaddingMeters;
    upper.z += settings.domainPaddingMeters;
    return {settings.cellCounts, lower, upper};
}

fluid::MacVelocityField makePredictor(
    const fluid::PeriodicCartesianGrid& grid,
    const FrozenScenePressureCaseSettings& settings) {
    fluid::MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const double sample = static_cast<double>(index + 1);
                result.xFaces()[index] =
                    settings.backgroundWindMetersPerSecond.x
                    + settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.01 * sample;
                result.yFaces()[index] =
                    settings.backgroundWindMetersPerSecond.y
                    + settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.02 * sample;
                result.zFaces()[index] =
                    settings.backgroundWindMetersPerSecond.z
                    - settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.015 * sample;
            }
        }
    }
    return result;
}

struct BuiltCase {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly assembly;
    Structure structure;
    SceneFluidSurfaceState surfaceState;
    SceneFluidSurfaceTransfer transfer;
    fluid::PeriodicCartesianGrid grid;
    SceneFluidGridEpoch gridEpoch;
    SceneFluidMimeticPressureAuditEndpoint pressure;
    ConservativeTransferResult pressureTransfer;
    viewer::StructureFrameMapping frameMapping;
};

BuiltCase buildCase(
    Scene scene,
    const FrozenScenePressureCaseSettings& settings) {
    validateSettings(settings);
    const auto validation = validateScene(scene);
    if (!validation.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure input is not a valid scene-v2 payload");
    }
    auto surface = assembleSceneFluidSurface(scene);
    if (!surface.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure surface failed: "
            + surfaceError(surface));
    }
    auto assembly = assembleSceneStructure(scene);
    if (!assembly.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure structure failed: "
            + assemblyError(assembly));
    }
    Structure structure(assembly.definition);
    auto state = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    SceneFluidSurfaceTransfer transfer(
        surface.definition, assembly.mappings, structure);
    auto grid = makeGrid(state, settings);
    auto gridEpoch = buildSceneFluidGridEpoch(
        surface.definition, state, grid, transfer);
    const auto openingCaps = buildSceneFluidOpeningCaps(
        surface.definition, state);
    const auto openingQuadrature = buildSceneFluidOpeningQuadrature(
        surface.definition, state, openingCaps);
    const auto openingPatches = buildSceneFluidOpeningGridPatches(
        surface.definition, state, openingCaps, openingQuadrature, grid);
    const auto openingFaceCrossings =
        buildSceneFluidOpeningFaceCrossings(
            surface.definition, state, openingCaps, openingQuadrature,
            openingPatches, grid);
    const auto cappedFacePartitions =
        buildSceneFluidCappedFacePartitions(
            surface.definition, state, grid, transfer, gridEpoch,
            openingCaps, openingQuadrature, openingPatches,
            openingFaceCrossings);
    const auto cellVolumes = buildSceneFluidCellVolumes(
        surface.definition, state, grid, transfer, gridEpoch);
    const auto connectivity = buildSceneFluidRegionConnectivity(
        surface.definition);
    const auto pressureVolumes = buildSceneFluidPressureControlVolumes(
        surface.definition, cellVolumes, connectivity);
    const auto pressureFaceLinks = buildSceneFluidPressureFaceLinks(
        surface.definition, state, grid, transfer, gridEpoch,
        openingCaps, openingQuadrature, openingPatches,
        openingFaceCrossings, cappedFacePartitions, cellVolumes,
        connectivity, pressureVolumes);
    const auto predictor = makePredictor(grid, settings);
    const auto openingFlux = evaluateSceneFluidOpeningFlux(
        surface.definition, state, openingCaps, openingQuadrature,
        openingPatches, grid, predictor);
    SceneFluidMimeticPressureAuditSettings pressureSettings;
    pressureSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    pressureSettings.timeStepSeconds = settings.timeStepSeconds;
    pressureSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    pressureSettings.pressureSolve.relativeResidualTolerance = 1.0e-5;
    pressureSettings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-8;
    pressureSettings.pressureSolve.maximumIterations = 1000;
    auto pressure = buildSceneFluidMimeticPressureAuditEndpoint(
        surface.definition, state, grid, gridEpoch, openingCaps,
        openingQuadrature, openingPatches, pressureVolumes,
        pressureFaceLinks, openingFlux, predictor, pressureSettings);
    if (!pressure.pressureEpoch.diagnostics.accepted) {
        throw std::runtime_error(
            "frozen scene mimetic pressure solve was not accepted");
    }
    auto pressureTransfer = evaluateSceneFluidMimeticPressureQuadrature(
        surface.definition, state, transfer, gridEpoch.quadrature,
        pressure.pressureEpoch.acceptedPressureSamples);
    auto frameMapping = viewer::makeStructureFrameMapping(
        scene, assembly, structure);
    return {
        std::move(scene), std::move(surface), std::move(assembly),
        std::move(structure), std::move(state), std::move(transfer),
        std::move(grid), std::move(gridEpoch), std::move(pressure),
        std::move(pressureTransfer), std::move(frameMapping),
    };
}

FrozenScenePressureCaseDiagnostics makeDiagnostics(
    const BuiltCase& built) {
    FrozenScenePressureCaseDiagnostics result;
    result.pressureControlCount =
        built.pressure.controlCells.controlCells.size();
    result.sharedTraceCount =
        built.pressure.fullTraceSystem.sharedTraceCount;
    result.pressureIterationCount = built.pressure.pressureEpoch
        .diagnostics.pressureSolve.reducedTraceSolve.iterationCount;
    result.maximumAbsolutePressureDifferencePascals =
        built.pressure.pressureEpoch.acceptedPressureSamples
            .maximumAbsolutePressureDifferencePascals;
    result.maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
        built.pressure.pressureSources
            .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond;
    const auto& transfer = built.pressureTransfer.diagnostics();
    result.pressureForceNewtons = transfer.transferredNodalForceNewtons;
    result.pressureMomentNewtonMeters =
        transfer.transferredNodalMomentNewtonMeters;
    result.transferForceResidualNewtons = transfer.forceResidualNormNewtons;
    result.transferMomentResidualNewtonMeters =
        transfer.momentResidualNormNewtonMeters;
    result.finite = built.pressure.pressureEpoch.diagnostics.accepted
        && transfer.finite
        && std::isfinite(result.maximumAbsolutePressureDifferencePascals)
        && std::isfinite(
            result.maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond)
        && std::isfinite(result.transferForceResidualNewtons)
        && std::isfinite(result.transferMomentResidualNewtonMeters);
    return result;
}

} // namespace

struct FrozenScenePressureCase::Implementation {
    Implementation(
        BuiltCase builtValue,
        const FrozenScenePressureCaseSettings& settingsValue)
        : built(std::move(builtValue)), settings(settingsValue),
          diagnostics(makeDiagnostics(built)) {
        stepSettings.timeStepSeconds = settings.timeStepSeconds;
        stepSettings.gravityMetersPerSecondSquared = {};
        if (!diagnostics.finite) {
            throw std::runtime_error(
                "frozen scene pressure diagnostics are non-finite");
        }
    }

    BuiltCase built;
    FrozenScenePressureCaseSettings settings;
    FrozenScenePressureCaseDiagnostics diagnostics;
    StructureStepSettings stepSettings;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
};

FrozenScenePressureCase::FrozenScenePressureCase(
    Scene scene,
    const FrozenScenePressureCaseSettings& settings)
    : implementation_(std::make_unique<Implementation>(
          buildCase(std::move(scene), settings), settings)) {}

FrozenScenePressureCase::~FrozenScenePressureCase() = default;
FrozenScenePressureCase::FrozenScenePressureCase(
    FrozenScenePressureCase&&) noexcept = default;
FrozenScenePressureCase& FrozenScenePressureCase::operator=(
    FrozenScenePressureCase&&) noexcept = default;

viewer::TraceHeader FrozenScenePressureCase::traceHeader() const {
    return {
        implementation_->built.scene.metadata.designChecksum,
        frozenScenePressureSolverId,
    };
}

viewer::DiagnosticFrame FrozenScenePressureCase::advance() {
    auto& state = *implementation_;
    const std::uint64_t nextStep = state.acceptedStepCount + 1;
    const double nextTime = state.simulationTimeSeconds
        + state.settings.timeStepSeconds;
    viewer::StructureFrameContext context;
    context.sceneChecksum =
        state.built.scene.metadata.designChecksum;
    context.solverCommit = frozenScenePressureSolverId;
    context.timeStepSeconds = state.settings.timeStepSeconds;
    context.couplingResiduals.fluid = state.diagnostics
        .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond;
    context.couplingResiduals.tractionNewtons =
        state.diagnostics.transferForceResidualNewtons;
    context.conservation.interfaceForceResidualNewtons = {
        state.built.pressureTransfer.diagnostics().forceResidualNewtons.x,
        state.built.pressureTransfer.diagnostics().forceResidualNewtons.y,
        state.built.pressureTransfer.diagnostics().forceResidualNewtons.z,
    };
    context.conservation.interfaceMomentResidualNewtonMetres = {
        state.built.pressureTransfer.diagnostics()
            .momentResidualNewtonMeters.x,
        state.built.pressureTransfer.diagnostics()
            .momentResidualNewtonMeters.y,
        state.built.pressureTransfer.diagnostics()
            .momentResidualNewtonMeters.z,
    };
    auto frame = viewer::buildStructureFrame(
        state.built.structure, state.built.frameMapping, context);
    frame.step = nextStep;
    frame.simulationTimeSeconds = nextTime;

    std::map<StableId, std::size_t> triangleIndices;
    for (std::size_t index = 0;
         index < state.built.assembly.mappings.triangleIds.size();
         ++index) {
        triangleIndices.emplace(
            state.built.assembly.mappings.triangleIds[index], index);
    }
    std::vector<double> trianglePressure(
        state.built.structure.definition().triangles.size(), 0.0);
    std::vector<double> triangleArea(trianglePressure.size(), 0.0);
    const auto& quadrature = state.built.gridEpoch.quadrature;
    const auto& samples = state.built.pressure.pressureEpoch
        .acceptedPressureSamples;
    if (samples.bindings.size() != quadrature.points.size()) {
        throw std::logic_error(
            "frozen scene pressure sample count changed");
    }
    for (std::size_t index = 0; index < quadrature.points.size(); ++index) {
        const auto found = triangleIndices.find(
            quadrature.points[index].triangleId);
        if (found == triangleIndices.end()) {
            throw std::logic_error(
                "frozen scene pressure triangle is missing");
        }
        const double area = quadrature.points[index].areaSquareMeters;
        trianglePressure[found->second] +=
            samples.bindings[index].pressureDifferencePascals * area;
        triangleArea[found->second] += area;
    }
    for (std::size_t index = 0; index < trianglePressure.size(); ++index) {
        if (!(triangleArea[index] > 0.0)) {
            throw std::logic_error(
                "frozen scene triangle has no pressure quadrature");
        }
        trianglePressure[index] /= triangleArea[index];
    }
    frame.scalarFields.push_back({
        "frozen_scene.pressure_jump", "Pa",
        viewer::FieldAssociation::Triangle, std::move(trianglePressure),
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_pressure_jump", "Pa",
        viewer::FieldAssociation::Global,
        {state.diagnostics.maximumAbsolutePressureDifferencePascals},
    });
    frame.scalarFields.push_back({
        "frozen_scene.pressure_iterations", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.pressureIterationCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.pressure_controls", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.pressureControlCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.shared_traces", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.sharedTraceCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.component_continuity_residual", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond},
    });

    std::vector<viewer::Vec3d> nodalPressureForces(
        frame.vertices.size());
    for (const auto& load : state.built.pressureTransfer.nodeLoads()) {
        nodalPressureForces[load.structureNode] = {
            load.forceNewtons.x, load.forceNewtons.y,
            load.forceNewtons.z,
        };
    }
    frame.vectorFields.push_back({
        "frozen_scene.nodal_pressure_force", "N",
        viewer::FieldAssociation::Vertex,
        std::move(nodalPressureForces),
    });
    frame.vectorFields.push_back({
        "frozen_scene.total_pressure_force", "N",
        viewer::FieldAssociation::Global,
        {{state.diagnostics.pressureForceNewtons.x,
          state.diagnostics.pressureForceNewtons.y,
          state.diagnostics.pressureForceNewtons.z}},
    });
    frame.vectorFields.push_back({
        "frozen_scene.total_pressure_moment", "N m",
        viewer::FieldAssociation::Global,
        {{state.diagnostics.pressureMomentNewtonMeters.x,
          state.diagnostics.pressureMomentNewtonMeters.y,
          state.diagnostics.pressureMomentNewtonMeters.z}},
    });
    frame.vectorFields.push_back({
        "frozen_scene.background_wind", "m/s",
        viewer::FieldAssociation::Global,
        {{state.settings.backgroundWindMetersPerSecond.x,
          state.settings.backgroundWindMetersPerSecond.y,
          state.settings.backgroundWindMetersPerSecond.z}},
    });

    viewer::ProtocolError error;
    if (!viewer::validateFrame(frame, &error)) {
        throw std::runtime_error(
            "frozen scene pressure frame is invalid: " + error.message);
    }
    state.acceptedStepCount = nextStep;
    state.simulationTimeSeconds = nextTime;
    return frame;
}

const StructureStepSettings&
FrozenScenePressureCase::stepSettings() const noexcept {
    return implementation_->stepSettings;
}

std::uint64_t FrozenScenePressureCase::acceptedStepCount() const noexcept {
    return implementation_->acceptedStepCount;
}

double FrozenScenePressureCase::simulationTimeSeconds() const noexcept {
    return implementation_->simulationTimeSeconds;
}

const FrozenScenePressureCaseDiagnostics&
FrozenScenePressureCase::diagnostics() const noexcept {
    return implementation_->diagnostics;
}

} // namespace simwing::fsi
