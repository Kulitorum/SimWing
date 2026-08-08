#include "projected_flag_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::size_t gridCellsX = 12;
constexpr std::size_t gridCellsY = 8;
constexpr std::size_t gridCellsZ = 8;
constexpr std::size_t panelFaceX = 6;
constexpr std::size_t panelFaceYBegin = 2;
constexpr std::size_t panelFaceZBegin = 2;
constexpr double panelPlaneXMeters = 1.5;
constexpr double panelMinimumYMeters = -0.5;
constexpr double panelMinimumZMeters = 0.5;
constexpr double panelSpacingMeters = 0.25;
constexpr double fabricArealDensityKgPerSquareMeter = 0.08;
constexpr double gustAmplitudeMetersPerSecond = 1.8;
constexpr double gustFrequencyHertz = 0.65;

[[nodiscard]] StructureVector3 subtract(const StructureVector3& first,
                                        const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] StructureVector3 cross(const StructureVector3& first,
                                     const StructureVector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

[[nodiscard]] double dot(const StructureVector3& first,
                         const StructureVector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

[[nodiscard]] double length(const StructureVector3& value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] std::size_t panelNode(const std::size_t y,
                                    const std::size_t z) {
    return z * projectedFlagNodesPerSide + y;
}

[[nodiscard]] fluid::PeriodicCartesianGrid makeGrid() {
    return {{gridCellsX, gridCellsY, gridCellsZ},
            {0.0, -1.0, 0.0},
            {3.0, 1.0, 2.0}};
}

[[nodiscard]] std::vector<fluid::GridFaceMovingInterface>
makeInterfaceFaces() {
    std::vector<fluid::GridFaceMovingInterface> faces;
    faces.reserve(projectedFlagTilesPerSide * projectedFlagTilesPerSide);
    for (std::size_t z = 0; z < projectedFlagTilesPerSide; ++z) {
        for (std::size_t y = 0; y < projectedFlagTilesPerSide; ++y) {
            faces.push_back({
                projectedFlagSurfaceStableId,
                1,
                1,
                fluid::GridFaceAxis::X,
                panelFaceX,
                panelFaceYBegin + y,
                panelFaceZBegin + z,
                0.0,
            });
        }
    }
    return faces;
}

[[nodiscard]] fluid::MovingInterfaceProjectionSettings
makeFluidSettings() {
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.225;
    settings.projection.timeStepSeconds = 1.0 / 120.0;
    settings.projection.absoluteResidualTolerance = 1.0e-10;
    settings.projection.relativeResidualTolerance = 1.0e-12;
    settings.projection.maximumIterations = 2000;
    settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond =
        1.0e-11;
    return settings;
}

[[nodiscard]] fluid::MovingInterfaceProjectionDiagnostics
makeReferenceFluidDiagnostics(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::FaceAlignedMovingInterface& interfaceValue,
    const fluid::MovingInterfaceProjectionSettings& settings) {
    fluid::MacVelocityField velocity(grid);
    fluid::CellScalarField pressure(grid);
    const auto diagnostics = fluid::projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaceValue, settings);
    if (!diagnostics.projection.converged || !diagnostics.finite
        || diagnostics.faces.size()
            != projectedFlagTilesPerSide * projectedFlagTilesPerSide) {
        throw std::logic_error(
            "projected flag could not construct its reference CFD surface");
    }
    return diagnostics;
}

[[nodiscard]] StructureMembraneMaterial membraneMaterial() {
    StructureMembraneMaterial material;
    material.warpStiffnessNewtonsPerMeter = 650.0;
    material.weftStiffnessNewtonsPerMeter = 500.0;
    material.couplingStiffnessNewtonsPerMeter = 80.0;
    material.shearStiffnessNewtonsPerMeter = 140.0;
    material.dampingSeconds = 0.01;
    material.compressionStiffnessRatio = 0.08;
    return material;
}

struct EdgeIncidence {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t opposite = 0;
};

void addRestShapeDihedrals(StructureDefinition& definition) {
    using Edge = std::pair<std::size_t, std::size_t>;
    std::map<Edge, std::vector<EdgeIncidence>> incidences;
    for (const auto& triangle : definition.triangles) {
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const std::size_t from = triangle.nodes[edge];
            const std::size_t to = triangle.nodes[(edge + 1) % 3];
            incidences[{std::min(from, to), std::max(from, to)}].push_back(
                {from, to, triangle.nodes[(edge + 2) % 3]});
        }
    }
    for (const auto& [edge, adjacent] : incidences) {
        static_cast<void>(edge);
        if (adjacent.size() == 1) {
            continue;
        }
        if (adjacent.size() != 2
            || adjacent[0].from != adjacent[1].to
            || adjacent[0].to != adjacent[1].from) {
            throw std::logic_error(
                "projected flag mesh is not an oriented manifold");
        }
        definition.dihedrals.push_back({
            {adjacent[0].from, adjacent[0].to,
             adjacent[0].opposite, adjacent[1].opposite},
            0.0,
            0.15,
        });
    }
}

[[nodiscard]] StructureDefinition makeDefinition() {
    StructureDefinition definition;
    definition.nodes.reserve(
        projectedFlagNodesPerSide * projectedFlagNodesPerSide);
    for (std::size_t z = 0; z < projectedFlagNodesPerSide; ++z) {
        for (std::size_t y = 0; y < projectedFlagNodesPerSide; ++y) {
            definition.nodes.push_back({
                {panelPlaneXMeters,
                 panelMinimumYMeters + panelSpacingMeters * y,
                panelMinimumZMeters + panelSpacingMeters * z},
                0.0,
                y <= 1,
            });
        }
    }

    definition.triangles.reserve(
        2 * projectedFlagTilesPerSide * projectedFlagTilesPerSide);
    for (std::size_t z = 0; z < projectedFlagTilesPerSide; ++z) {
        for (std::size_t y = 0; y < projectedFlagTilesPerSide; ++y) {
            const std::size_t lower = panelNode(y, z);
            const std::size_t right = panelNode(y + 1, z);
            const std::size_t upperRight = panelNode(y + 1, z + 1);
            const std::size_t upper = panelNode(y, z + 1);
            definition.triangles.push_back({{lower, right, upperRight}});
            definition.triangles.push_back({{lower, upperRight, upper}});
        }
    }

    const StructureMembraneMaterial material = membraneMaterial();
    definition.membranes.reserve(definition.triangles.size());
    for (std::size_t index = 0;
         index < definition.triangles.size(); ++index) {
        const auto& triangle = definition.triangles[index];
        const auto& first =
            definition.nodes[triangle.nodes[0]].positionMeters;
        const auto& second =
            definition.nodes[triangle.nodes[1]].positionMeters;
        const auto& third =
            definition.nodes[triangle.nodes[2]].positionMeters;
        const double area = 0.5 * length(cross(
            subtract(second, first), subtract(third, first)));
        for (const std::size_t node : triangle.nodes) {
            definition.nodes[node].massKg +=
                fabricArealDensityKgPerSquareMeter * area / 3.0;
        }
        definition.membranes.push_back({
            index,
            {StructureVector2{first.y, first.z},
             StructureVector2{second.y, second.z},
             StructureVector2{third.y, third.z}},
            material,
            StructureMaterialRole::Bulk,
        });
    }
    addRestShapeDihedrals(definition);
    return definition;
}

[[nodiscard]] std::vector<CouplingSurfaceNodeDefinition>
makeCouplingNodes() {
    std::vector<CouplingSurfaceNodeDefinition> nodes;
    nodes.reserve(projectedFlagNodesPerSide * projectedFlagNodesPerSide);
    for (std::size_t index = 0;
         index < projectedFlagNodesPerSide * projectedFlagNodesPerSide;
         ++index) {
        nodes.push_back({100'000 + index, index});
    }
    return nodes;
}

[[nodiscard]] std::vector<CouplingSurfaceTriangleDefinition>
makeCouplingTriangles() {
    std::vector<CouplingSurfaceTriangleDefinition> triangles;
    triangles.reserve(
        2 * projectedFlagTilesPerSide * projectedFlagTilesPerSide);
    for (std::size_t z = 0; z < projectedFlagTilesPerSide; ++z) {
        for (std::size_t y = 0; y < projectedFlagTilesPerSide; ++y) {
            const std::size_t lower = panelNode(y, z);
            const std::size_t right = panelNode(y + 1, z);
            const std::size_t upperRight = panelNode(y + 1, z + 1);
            const std::size_t upper = panelNode(y, z + 1);
            triangles.push_back(
                {200'000 + triangles.size(),
                 {100'000 + lower, 100'000 + right,
                  100'000 + upperRight}});
            triangles.push_back(
                {200'000 + triangles.size(),
                 {100'000 + lower, 100'000 + upperRight,
                  100'000 + upper}});
        }
    }
    return triangles;
}

[[nodiscard]] viewer::StructureFrameMappingDefinition makeFrameMapping() {
    viewer::StructureFrameMappingDefinition mapping;
    mapping.vertexStableIds.reserve(
        projectedFlagNodesPerSide * projectedFlagNodesPerSide);
    for (std::size_t index = 0;
         index < projectedFlagNodesPerSide * projectedFlagNodesPerSide;
         ++index) {
        mapping.vertexStableIds.push_back(100'000 + index);
    }
    mapping.triangles.reserve(
        2 * projectedFlagTilesPerSide * projectedFlagTilesPerSide);
    for (std::size_t index = 0;
         index < 2 * projectedFlagTilesPerSide
                       * projectedFlagTilesPerSide;
         ++index) {
        mapping.triangles.push_back({200'000 + index, 1, 1});
    }
    return mapping;
}

[[nodiscard]] StructureStepSettings makeStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 120.0;
    settings.substeps = 4;
    settings.constraintIterations = 20;
    settings.cableConstraintSweepPairs = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 1.5;
    settings.workerThreads = 0;
    return settings;
}

[[nodiscard]] const fluid::MovingInterfaceSurfaceDiagnostics*
findSurface(const fluid::MovingInterfaceProjectionDiagnostics& diagnostics) {
    const auto found = std::ranges::find_if(
        diagnostics.surfaces,
        [](const auto& surface) {
            return surface.stableId == projectedFlagSurfaceStableId;
        });
    return found == diagnostics.surfaces.end() ? nullptr : &*found;
}

} // namespace

ProjectedGustFlagCase::ProjectedGustFlagCase()
    : grid_(makeGrid()),
      interface_(grid_, makeInterfaceFaces()),
      velocity_(grid_),
      pressure_(grid_),
      structure_(makeDefinition()),
      fluidSettings_(makeFluidSettings()),
      bridge_(
          structure_, projectedFlagSurfaceStableId,
          makeCouplingNodes(), makeCouplingTriangles(),
          makeReferenceFluidDiagnostics(
              grid_, interface_, fluidSettings_).faces),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()),
      referenceKinematics_(
          bridge_.transfer().captureKinematics(structure_)),
      lastNodalForcesNewtons_(
          projectedFlagNodesPerSide * projectedFlagNodesPerSide) {}

viewer::TraceHeader ProjectedGustFlagCase::traceHeader() const {
    return {projectedGustFlagCaseChecksum, projectedGustFlagCaseSolverId};
}

viewer::DiagnosticFrame ProjectedGustFlagCase::advance() {
    const StructureCheckpoint structureBefore = structure_.checkpoint();
    const fluid::MacVelocityField velocityBefore = velocity_;
    const fluid::CellScalarField pressureBefore = pressure_;
    const auto fluidDiagnosticsBefore = fluidDiagnostics_;
    const auto transferDiagnosticsBefore = transferDiagnostics_;
    const auto nodalForcesBefore = lastNodalForcesNewtons_;
    const double gustBefore = gustSpeedMetersPerSecond_;
    try {
        const double nextTime = structure_.simulationTimeSeconds()
            + stepSettings_.timeStepSeconds;
        const double nextGust = gustAmplitudeMetersPerSecond
            * std::sin(2.0 * std::numbers::pi
                       * gustFrequencyHertz * nextTime);
        const double gustIncrement =
            nextGust - gustSpeedMetersPerSecond_;
        for (double& xFaceVelocity : velocity_.xFaces()) {
            xFaceVelocity += gustIncrement;
        }

        const auto projected = fluid::projectVelocityWithMovingInterfaces(
            grid_, velocity_, pressure_, interface_, fluidSettings_);
        if (!projected.projection.converged || !projected.finite) {
            throw std::runtime_error(
                "projected flag CFD projection did not converge");
        }
        const auto transferred = bridge_.evaluateConstraintReaction(
            projected, referenceKinematics_);
        bridge_.transfer().addLoadsTo(
            structure_, transferred.transferResult());

        std::ranges::fill(
            lastNodalForcesNewtons_, StructureVector3{});
        for (const auto& load : transferred.transferResult().nodeLoads()) {
            lastNodalForcesNewtons_[load.structureNode] =
                load.forceNewtons;
        }
        const StructureDiagnostics structureDiagnostics =
            structure_.step(stepSettings_);
        if (!structureDiagnostics.finite) {
            throw std::runtime_error(
                "projected flag XPBD step produced non-finite diagnostics");
        }

        fluidDiagnostics_ = projected;
        transferDiagnostics_ = transferred.diagnostics();
        gustSpeedMetersPerSecond_ = nextGust;

        viewer::StructureFrameContext context;
        context.sceneChecksum = projectedGustFlagCaseChecksum;
        context.solverCommit = projectedGustFlagCaseSolverId;
        context.timeStepSeconds = stepSettings_.timeStepSeconds;
        context.couplingIteration = 0;
        context.couplingResiduals.tractionNewtons =
            transferDiagnostics_.forceResidualNormNewtons;
        context.couplingResiduals.structure =
            structureDiagnostics.maximumMembraneResidual;
        context.couplingResiduals.interfacePowerWatts =
            transferDiagnostics_.powerResidualWatts;
        context.conservation.interfaceForceResidualNewtons = {
            transferDiagnostics_.forceResidualNewtons.x,
            transferDiagnostics_.forceResidualNewtons.y,
            transferDiagnostics_.forceResidualNewtons.z,
        };
        context.conservation.interfaceMomentResidualNewtonMetres = {
            transferDiagnostics_.momentResidualNewtonMeters.x,
            transferDiagnostics_.momentResidualNewtonMeters.y,
            transferDiagnostics_.momentResidualNewtonMeters.z,
        };
        context.conservation.interfacePowerResidualWatts =
            transferDiagnostics_.powerResidualWatts;
        viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
            structure_, frameMapping_, context);

        std::vector<double> normalDisplacements;
        normalDisplacements.reserve(frame.vertices.size());
        for (const auto& vertex : frame.vertices) {
            normalDisplacements.push_back(
                vertex.positionMetres.x - panelPlaneXMeters);
        }
        frame.scalarFields.push_back({
            "flag.normal_displacement", "m",
            viewer::FieldAssociation::Vertex,
            std::move(normalDisplacements),
        });

        const std::size_t triangleCount =
            2 * projectedFlagTilesPerSide * projectedFlagTilesPerSide;
        std::vector<double> pressureTraction(
            triangleCount, 0.0);
        std::vector<double> directConstraintTraction(
            triangleCount, 0.0);
        std::vector<double> completeReactionTraction(
            triangleCount, 0.0);
        for (const auto& face : fluidDiagnostics_.faces) {
            if (face.surfaceStableId != projectedFlagSurfaceStableId) {
                continue;
            }
            if (face.j < panelFaceYBegin || face.k < panelFaceZBegin
                || face.j >= panelFaceYBegin + projectedFlagTilesPerSide
                || face.k >= panelFaceZBegin + projectedFlagTilesPerSide) {
                throw std::logic_error(
                    "projected flag CFD face escaped its reference panel");
            }
            const std::size_t tile =
                (face.k - panelFaceZBegin) * projectedFlagTilesPerSide
                + (face.j - panelFaceYBegin);
            for (const std::size_t triangle :
                 std::array<std::size_t, 2>{2 * tile, 2 * tile + 1}) {
                pressureTraction[triangle] =
                    face.pressureTractionPascals.x;
                directConstraintTraction[triangle] =
                    face.directConstraintForceNewtons.x
                    / face.areaSquareMeters;
                completeReactionTraction[triangle] =
                    face.constraintReactionTractionPascals.x;
            }
        }
        frame.scalarFields.push_back({
            "flag.cfd_pressure_traction", "Pa",
            viewer::FieldAssociation::Triangle,
            std::move(pressureTraction),
        });
        frame.scalarFields.push_back({
            "flag.cfd_direct_constraint_traction", "Pa",
            viewer::FieldAssociation::Triangle,
            std::move(directConstraintTraction),
        });
        frame.scalarFields.push_back({
            "flag.cfd_complete_reaction_traction", "Pa",
            viewer::FieldAssociation::Triangle,
            std::move(completeReactionTraction),
        });
        frame.scalarFields.push_back({
            "flag.gust_speed", "m/s",
            viewer::FieldAssociation::Global,
            {gustSpeedMetersPerSecond_},
        });
        frame.scalarFields.push_back({
            "flag.fluid_divergence_l2", "1/s",
            viewer::FieldAssociation::Global,
            {fluidDiagnostics_.projection.divergenceL2AfterPerSecond},
        });
        const auto* surface = findSurface(fluidDiagnostics_);
        if (surface == nullptr) {
            throw std::logic_error(
                "projected flag CFD surface aggregate disappeared");
        }
        frame.scalarFields.push_back({
            "flag.pressure_force_x", "N",
            viewer::FieldAssociation::Global,
            {surface->pressureForceNewtons.x},
        });
        frame.scalarFields.push_back({
            "flag.direct_constraint_force_x", "N",
            viewer::FieldAssociation::Global,
            {surface->directConstraintForceNewtons.x},
        });
        frame.scalarFields.push_back({
            "flag.complete_reaction_force_x", "N",
            viewer::FieldAssociation::Global,
            {surface->constraintReactionForceNewtons.x},
        });

        std::vector<viewer::Vec3d> nodalForces;
        nodalForces.reserve(lastNodalForcesNewtons_.size());
        for (const auto& force : lastNodalForcesNewtons_) {
            nodalForces.push_back({force.x, force.y, force.z});
        }
        frame.vectorFields.push_back({
            "flag.cfd_nodal_force", "N",
            viewer::FieldAssociation::Vertex,
            std::move(nodalForces),
        });

        viewer::ProtocolError error;
        if (!viewer::validateFrame(frame, &error)) {
            throw std::logic_error(
                "projected flag frame is invalid: " + error.message);
        }
        return frame;
    } catch (...) {
        structure_.restore(structureBefore);
        velocity_ = velocityBefore;
        pressure_ = pressureBefore;
        fluidDiagnostics_ = fluidDiagnosticsBefore;
        transferDiagnostics_ = transferDiagnosticsBefore;
        lastNodalForcesNewtons_ = nodalForcesBefore;
        gustSpeedMetersPerSecond_ = gustBefore;
        throw;
    }
}

const Structure& ProjectedGustFlagCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& ProjectedGustFlagCase::stepSettings() const
    noexcept {
    return stepSettings_;
}

const fluid::MacVelocityField& ProjectedGustFlagCase::velocity() const
    noexcept {
    return velocity_;
}

const fluid::CellScalarField& ProjectedGustFlagCase::pressure() const
    noexcept {
    return pressure_;
}

const fluid::MovingInterfaceProjectionDiagnostics&
ProjectedGustFlagCase::fluidDiagnostics() const noexcept {
    return fluidDiagnostics_;
}

const PlanarFaceResolvedBridgeDiagnostics&
ProjectedGustFlagCase::transferDiagnostics() const noexcept {
    return transferDiagnostics_;
}

double ProjectedGustFlagCase::gustSpeedMetersPerSecond() const noexcept {
    return gustSpeedMetersPerSecond_;
}

double ProjectedGustFlagCase::maximumNormalDisplacementMeters() const {
    double maximum = 0.0;
    for (const auto& node : structure_.nodeStates()) {
        maximum = std::max(
            maximum, std::abs(node.positionMeters.x - panelPlaneXMeters));
    }
    return maximum;
}

double ProjectedGustFlagCase::maximumFreeEdgeDisplacementMeters() const {
    const auto states = structure_.nodeStates();
    const auto& definition = structure_.definition();
    double maximum = 0.0;
    for (std::size_t z = 0; z < projectedFlagNodesPerSide; ++z) {
        const std::size_t node =
            panelNode(projectedFlagNodesPerSide - 1, z);
        maximum = std::max(
            maximum,
            length(subtract(
                states[node].positionMeters,
                definition.nodes[node].positionMeters)));
    }
    return maximum;
}

} // namespace simwing::fsi
