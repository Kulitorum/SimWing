#include "scene_pressure_cell_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

constexpr StableId apexVertexId = 10;
constexpr double targetMeanWindMetersPerSecond = -0.85;

template<std::size_t VertexCount>
std::array<Vec2, 3> intrinsicChart(
    const std::array<Vec3, VertexCount>& positions,
    const std::array<std::size_t, 3>& vertices) {
    const Vec3& first = positions[vertices[0]];
    const Vec3& second = positions[vertices[1]];
    const Vec3& third = positions[vertices[2]];
    const Vec3 edge{second.x - first.x,
                    second.y - first.y,
                    second.z - first.z};
    const Vec3 diagonal{third.x - first.x,
                        third.y - first.y,
                        third.z - first.z};
    const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
    const double projected = (diagonal.x * edge.x
                              + diagonal.y * edge.y
                              + diagonal.z * edge.z) / edgeLength;
    const double diagonalSquared = diagonal.x * diagonal.x
        + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
    return {{{0.0, 0.0},
             {edgeLength, 0.0},
             {projected, std::sqrt(std::max(
                 0.0, diagonalSquared - projected * projected))}}};
}

Scene makeScene() {
    Scene scene;
    scene.metadata.designChecksum = scenePressureCellCaseChecksum;
    scene.metadata.exporterVersion = scenePressureCellCaseSolverId;
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "soft diagnostic fabric", 8.0, 6.0, 2.0, 0.0002,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 4> positions{{
        {1.2, 1.5, 1.45},
        {2.0, 1.2, 1.15},
        {2.0, 1.8, 1.15},
        {2.0, 1.5, 1.75},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 3> faces{{
        {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        scene.triangles.push_back({
            500 + face,
            {10 + faces[face][0], 10 + faces[face][1],
             10 + faces[face][2]},
            intrinsicChart(positions, faces[face]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

SceneStructureAssembly makeAssembly(const Scene& scene) {
    auto assembly = assembleSceneStructure(scene);
    for (std::size_t node = 0;
         node < assembly.mappings.nodeVertexIds.size(); ++node) {
        if (assembly.mappings.nodeVertexIds[node] != apexVertexId) {
            assembly.definition.nodes[node].fixed = true;
        }
    }
    return assembly;
}

fluid::PeriodicCartesianGrid makeGrid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

SceneFluidPressureCouplingSettings makeSettings() {
    SceneFluidPressureCouplingSettings settings;
    settings.structure.timeStepSeconds = 1.0 / 60.0;
    settings.structure.substeps = 4;
    settings.structure.constraintIterations = 10;
    settings.structure.gravityMetersPerSecondSquared = {};
    settings.structure.velocityDampingPerSecond = 0.5;
    settings.pressureProjection.timeStepSeconds =
        settings.structure.timeStepSeconds;
    settings.pressureProjection
        .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond = 2.0e-11;
    settings.pressureProjection.relativeCorrectedVolumeRateTolerance =
        1.0e-10;
    settings.pressureProjection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.pressureProjection.pressureSolve.relativeResidualTolerance =
        1.0e-12;
    settings.pressureProjection.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    settings.relaxation.initialRelaxation = 0.5;
    settings.relaxation.minimumRelaxation = 0.02;
    settings.relaxation.maximumRelaxation = 1.0;
    settings.convergence.minimumIterations = 2;
    settings.convergence.maximumIterations = 20;
    settings.convergence.absoluteDisplacementToleranceMetres = 1.0e-7;
    settings.convergence.relativeDisplacementTolerance = 1.0e-4;
    settings.convergence.absoluteVelocityToleranceMetersPerSecond = 1.0e-6;
    settings.convergence.relativeVelocityTolerance = 1.0e-4;
    settings.convergence.absoluteTractionToleranceNewtons = 1.0e-6;
    settings.convergence.relativeTractionTolerance = 1.0e-4;
    return settings;
}

viewer::StructureFrameMappingDefinition makeFrameMapping(
    const Scene& scene,
    const SceneStructureAssembly& assembly) {
    viewer::StructureFrameMappingDefinition result;
    result.vertexStableIds = assembly.mappings.nodeVertexIds;
    std::map<StableId, const Triangle*> triangles;
    for (const auto& triangle : scene.triangles) {
        triangles.emplace(triangle.id, &triangle);
    }
    for (const StableId id : assembly.mappings.triangleIds) {
        const auto found = triangles.find(id);
        if (found == triangles.end()) {
            throw std::logic_error(
                "scene pressure cell frame triangle is missing");
        }
        result.triangles.push_back({
            id, found->second->negativeSideRegionId,
            found->second->positiveSideRegionId,
        });
    }
    return result;
}

fluid::MacVelocityField makeInitialVelocity(
    const fluid::PeriodicCartesianGrid& grid) {
    fluid::MacVelocityField result(grid);
    std::ranges::fill(
        result.xFaces(), targetMeanWindMetersPerSecond);
    return result;
}

double meanXVelocity(const fluid::MacVelocityField& velocity) {
    double result = 0.0;
    for (const double value : velocity.xFaces()) {
        result += value;
    }
    return result / static_cast<double>(velocity.xFaces().size());
}

double gridVolume(const fluid::PeriodicCartesianGrid& grid) {
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    return (upper.x - lower.x) * (upper.y - lower.y)
        * (upper.z - lower.z);
}

double vectorNorm(const StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

double maximumDisplacement(const Structure& structure) {
    const auto states = structure.nodeStates();
    const auto& definition = structure.definition();
    double maximum = 0.0;
    for (std::size_t index = 0; index < states.size(); ++index) {
        const StructureVector3 difference{
            states[index].positionMeters.x
                - definition.nodes[index].positionMeters.x,
            states[index].positionMeters.y
                - definition.nodes[index].positionMeters.y,
            states[index].positionMeters.z
                - definition.nodes[index].positionMeters.z,
        };
        maximum = std::max(maximum, vectorNorm(difference));
    }
    return maximum;
}

bool finite(const ScenePressureCellDiagnostics& diagnostics) {
    return diagnostics.coupling.finite
        && diagnostics.macVelocity.finite
        && diagnostics.bulkAdvection.finite
        && diagnostics.bulkAdvection.accepted
        && std::isfinite(diagnostics.targetMeanWindMetersPerSecond)
        && std::isfinite(diagnostics.meanWindBeforePumpMetersPerSecond)
        && std::isfinite(diagnostics.flowPumpForceNewtons)
        && std::isfinite(diagnostics.pressureForceNewtons.x)
        && std::isfinite(diagnostics.pressureForceNewtons.y)
        && std::isfinite(diagnostics.pressureForceNewtons.z)
        && std::isfinite(diagnostics.maximumAbsolutePressurePascals)
        && std::isfinite(diagnostics.maximumDisplacementMeters);
}

} // namespace

ScenePressureCellCase::ScenePressureCellCase()
    : scene_(makeScene()),
      surface_(assembleSceneFluidSurface(scene_)),
      assembly_(makeAssembly(scene_)),
      structure_(assembly_.definition),
      coupling_(surface_.definition, assembly_.mappings, structure_,
                makeGrid(), makeSettings()),
      predictedVelocity_(makeInitialVelocity(coupling_.grid())),
      frameMapping_(structure_, makeFrameMapping(scene_, assembly_)) {}

viewer::TraceHeader ScenePressureCellCase::traceHeader() const {
    return {scenePressureCellCaseChecksum, scenePressureCellCaseSolverId};
}

viewer::DiagnosticFrame ScenePressureCellCase::advance() {
    const double meanWind = meanXVelocity(predictedVelocity_);
    const double pumpVelocityIncrement =
        targetMeanWindMetersPerSecond - meanWind;
    for (double& value : predictedVelocity_.xFaces()) {
        value += pumpVelocityIncrement;
    }
    fluid::CellScalarField bulkPressure(coupling_.grid());
    fluid::ProjectedMacAdvectionSspRk2Settings bulkSettings;
    bulkSettings.densityKgPerCubicMeter = coupling_.settings()
        .pressureProjection.densityKgPerCubicMeter;
    bulkSettings.timeStepSeconds =
        coupling_.settings().structure.timeStepSeconds;
    const auto bulk = fluid::advectVelocityProjectedSspRk2(
        coupling_.grid(), predictedVelocity_, bulkPressure, bulkSettings);
    if (!bulk.accepted) {
        throw std::runtime_error(
            "scene pressure cell bulk advection rejected its collapsed MAC predictor");
    }
    const auto coupled = coupling_.advance(structure_, predictedVelocity_);
    if (!coupled.accepted) {
        throw std::runtime_error(
            "scene pressure cell exhausted its coupling iteration budget");
    }
    const auto correctedMac =
        coupling_.acceptedPressureCorrectedMacVelocity();
    predictedVelocity_ = correctedMac.velocityMetersPerSecond;

    ScenePressureCellDiagnostics nextDiagnostics;
    nextDiagnostics.coupling = coupled;
    nextDiagnostics.macVelocity = correctedMac.diagnostics;
    nextDiagnostics.bulkAdvection = bulk;
    nextDiagnostics.targetMeanWindMetersPerSecond =
        targetMeanWindMetersPerSecond;
    nextDiagnostics.meanWindBeforePumpMetersPerSecond = meanWind;
    nextDiagnostics.flowPumpForceNewtons =
        bulkSettings.densityKgPerCubicMeter
        * gridVolume(coupling_.grid()) * pumpVelocityIncrement
        / bulkSettings.timeStepSeconds;
    nextDiagnostics.pressureForceNewtons =
        coupling_.acceptedPressureTransfer().diagnostics()
            .transferredNodalForceNewtons;
    const auto* projection = coupling_.acceptedPressureProjection();
    for (const double pressure : projection->pressurePascals) {
        nextDiagnostics.maximumAbsolutePressurePascals = std::max(
            nextDiagnostics.maximumAbsolutePressurePascals,
            std::abs(pressure));
    }
    nextDiagnostics.maximumDisplacementMeters =
        maximumDisplacement(structure_);
    nextDiagnostics.finite = finite(nextDiagnostics);
    if (!nextDiagnostics.finite) {
        throw std::runtime_error(
            "scene pressure cell diagnostics are non-finite");
    }
    diagnostics_ = nextDiagnostics;

    viewer::StructureFrameContext context;
    context.sceneChecksum = scenePressureCellCaseChecksum;
    context.solverCommit = scenePressureCellCaseSolverId;
    context.timeStepSeconds =
        coupling_.settings().structure.timeStepSeconds;
    context.couplingIteration = static_cast<std::uint32_t>(
        coupled.iteration.convergence.iteration);
    context.couplingResiduals.displacementMetres =
        coupled.iteration.convergence.residuals.displacementMetres;
    context.couplingResiduals.tractionNewtons =
        coupled.interfaceForceClosureNewtons;
    context.couplingResiduals.fluid = coupled.pressureProjection
        .correctedContinuityResidualMaximumCubicMetersPerSecond;
    context.couplingResiduals.structure =
        coupled.structure.maximumMembraneResidual;
    context.conservation.interfaceForceResidualNewtons = {
        coupled.pressureTransfer.forceResidualNewtons.x,
        coupled.pressureTransfer.forceResidualNewtons.y,
        coupled.pressureTransfer.forceResidualNewtons.z,
    };
    context.conservation.interfaceMomentResidualNewtonMetres = {
        coupled.pressureTransfer.momentResidualNewtonMeters.x,
        coupled.pressureTransfer.momentResidualNewtonMeters.y,
        coupled.pressureTransfer.momentResidualNewtonMeters.z,
    };
    viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
        structure_, frameMapping_, context);

    const auto states = structure_.nodeStates();
    std::vector<double> displacement;
    displacement.reserve(states.size());
    for (std::size_t index = 0; index < states.size(); ++index) {
        const auto& reference = structure_.definition().nodes[index]
            .positionMeters;
        displacement.push_back(vectorNorm({
            states[index].positionMeters.x - reference.x,
            states[index].positionMeters.y - reference.y,
            states[index].positionMeters.z - reference.z,
        }));
    }
    frame.scalarFields.push_back({
        "pressure_cell.displacement", "m",
        viewer::FieldAssociation::Vertex, std::move(displacement),
    });

    std::map<StableId, std::size_t> triangleIndices;
    for (std::size_t index = 0;
         index < assembly_.mappings.triangleIds.size(); ++index) {
        triangleIndices.emplace(
            assembly_.mappings.triangleIds[index], index);
    }
    std::vector<double> trianglePressure(
        structure_.definition().triangles.size(), 0.0);
    std::vector<double> triangleArea(trianglePressure.size(), 0.0);
    const auto* samples = coupling_.acceptedPressureSamples();
    const auto& quadrature = coupling_.acceptedPressureEpoch()
        .gridEpoch.quadrature;
    for (std::size_t index = 0;
         index < quadrature.points.size(); ++index) {
        const auto found = triangleIndices.find(
            quadrature.points[index].triangleId);
        if (found == triangleIndices.end()) {
            throw std::logic_error(
                "scene pressure cell sample triangle is missing");
        }
        const double area = quadrature.points[index].areaSquareMeters;
        trianglePressure[found->second] +=
            samples->bindings[index].pressureDifferencePascals * area;
        triangleArea[found->second] += area;
    }
    for (std::size_t index = 0; index < trianglePressure.size(); ++index) {
        if (!(triangleArea[index] > 0.0)) {
            throw std::logic_error(
                "scene pressure cell triangle has no pressure area");
        }
        trianglePressure[index] /= triangleArea[index];
    }
    frame.scalarFields.push_back({
        "pressure_cell.pressure_jump", "Pa",
        viewer::FieldAssociation::Triangle, std::move(trianglePressure),
    });
    frame.scalarFields.push_back({
        "pressure_cell.maximum_pressure", "Pa",
        viewer::FieldAssociation::Global,
        {diagnostics_.maximumAbsolutePressurePascals},
    });
    frame.scalarFields.push_back({
        "pressure_cell.coupling_iterations", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(diagnostics_.coupling.solverRunCount)},
    });
    frame.scalarFields.push_back({
        "pressure_cell.maximum_mac_speed", "m/s",
        viewer::FieldAssociation::Global,
        {diagnostics_.macVelocity
             .maximumAbsoluteVelocityMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "pressure_cell.mac_subface_deviation", "m/s",
        viewer::FieldAssociation::Global,
        {diagnostics_.macVelocity
             .maximumSubfaceVelocityDeviationMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "pressure_cell.bulk_advection_change", "m/s",
        viewer::FieldAssociation::Global,
        {diagnostics_.bulkAdvection
             .maximumVelocityChangeMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "pressure_cell.bulk_divergence", "1/s",
        viewer::FieldAssociation::Global,
        {diagnostics_.bulkAdvection.finalDivergenceL2PerSecond},
    });

    std::vector<viewer::Vec3d> nodalPressureForces(
        structure_.definition().nodes.size());
    for (const auto& load : coupling_.acceptedPressureTransfer().nodeLoads()) {
        nodalPressureForces[load.structureNode] = {
            load.forceNewtons.x, load.forceNewtons.y, load.forceNewtons.z,
        };
    }
    frame.vectorFields.push_back({
        "pressure_cell.nodal_pressure_force", "N",
        viewer::FieldAssociation::Vertex, std::move(nodalPressureForces),
    });
    frame.vectorFields.push_back({
        "pressure_cell.total_pressure_force", "N",
        viewer::FieldAssociation::Global,
        {{diagnostics_.pressureForceNewtons.x,
          diagnostics_.pressureForceNewtons.y,
          diagnostics_.pressureForceNewtons.z}},
    });
    frame.vectorFields.push_back({
        "pressure_cell.flow_pump_force", "N",
        viewer::FieldAssociation::Global,
        {{diagnostics_.flowPumpForceNewtons, 0.0, 0.0}},
    });

    viewer::ProtocolError error;
    if (!viewer::validateFrame(frame, &error)) {
        throw std::runtime_error(
            "scene pressure cell frame is invalid: " + error.message);
    }
    return frame;
}

const Structure& ScenePressureCellCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings&
ScenePressureCellCase::stepSettings() const noexcept {
    return coupling_.settings().structure;
}

std::uint64_t ScenePressureCellCase::acceptedStepCount() const noexcept {
    return structure_.acceptedStepCount();
}

double ScenePressureCellCase::simulationTimeSeconds() const noexcept {
    return structure_.simulationTimeSeconds();
}

const fluid::MacVelocityField&
ScenePressureCellCase::predictedVelocity() const noexcept {
    return predictedVelocity_;
}

const ScenePressureCellDiagnostics&
ScenePressureCellCase::diagnostics() const noexcept {
    return diagnostics_;
}

ScenePressureCellCheckpoint ScenePressureCellCase::checkpoint() const {
    ScenePressureCellCheckpoint result;
    result.coupling = coupling_.checkpoint(structure_);
    return result;
}

void ScenePressureCellCase::restore(
    const ScenePressureCellCheckpoint& checkpointValue) {
    if (checkpointValue.version != scenePressureCellCheckpointVersion) {
        throw std::invalid_argument(
            "scene pressure cell checkpoint version is invalid");
    }
    coupling_.restore(structure_, checkpointValue.coupling);
    predictedVelocity_ = makeInitialVelocity(coupling_.grid());
    if (coupling_.acceptedPressureProjection() != nullptr) {
        predictedVelocity_ = coupling_
            .acceptedPressureCorrectedMacVelocity()
            .velocityMetersPerSecond;
    }
    diagnostics_ = {};
}

} // namespace simwing::fsi
