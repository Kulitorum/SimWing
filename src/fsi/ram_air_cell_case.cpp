#include "ram_air_cell_case.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::size_t gridCellsX = 16;
constexpr std::size_t gridCellsY = 12;
constexpr std::size_t gridCellsZ = 12;
constexpr std::size_t cellFrontX = 4;
constexpr std::size_t cellBackX = 8;
constexpr std::size_t cellLeftY = 4;
constexpr std::size_t cellRightY = 8;
constexpr std::size_t cellBottomZ = 4;
constexpr std::size_t cellTopZ = 8;
constexpr double cellMinimumXMeters = 1.0;
constexpr double cellMinimumYMeters = -0.5;
constexpr double cellMinimumZMeters = 1.0;
constexpr double cellSpacingMeters = 0.25;
constexpr double fabricArealDensityKgPerSquareMeter = 0.08;
constexpr double gustAmplitudeMetersPerSecond = 2.5;
constexpr double gustFrequencyHertz = 0.45;
constexpr std::size_t trianglesPerPanel =
    2 * ramAirCellTilesPerEdge * ramAirCellTilesPerEdge;
constexpr std::size_t cellTriangleCount =
    ramAirCellPanelCount * trianglesPerPanel;

enum class PanelKind : std::size_t {
    Back = 0,
    Left = 1,
    Right = 2,
    Bottom = 3,
    Top = 4,
};

struct PanelDescriptor {
    PanelKind kind;
    std::uint64_t surfaceStableId;
    fluid::GridFaceAxis axis;
    std::size_t planeIndex;
};

constexpr std::array<PanelDescriptor, ramAirCellPanelCount> panels{{
    {PanelKind::Back, 51'000, fluid::GridFaceAxis::X, cellBackX},
    {PanelKind::Left, 52'000, fluid::GridFaceAxis::Y, cellLeftY},
    {PanelKind::Right, 53'000, fluid::GridFaceAxis::Y, cellRightY},
    {PanelKind::Bottom, 54'000, fluid::GridFaceAxis::Z, cellBottomZ},
    {PanelKind::Top, 55'000, fluid::GridFaceAxis::Z, cellTopZ},
}};

using NodeKey = std::array<std::size_t, 3>;
using NodeMap = std::map<NodeKey, std::size_t>;

[[nodiscard]] StructureVector3 add(const StructureVector3& first,
                                   const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

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

[[nodiscard]] StructureVector3 scaled(const StructureVector3& value,
                                      const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

[[nodiscard]] double dot(const StructureVector3& first,
                         const StructureVector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

[[nodiscard]] double length(const StructureVector3& value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] bool isSurfaceNode(const std::size_t x,
                                 const std::size_t y,
                                 const std::size_t z) {
    return x == ramAirCellTilesPerEdge
        || y == 0 || y == ramAirCellTilesPerEdge
        || z == 0 || z == ramAirCellTilesPerEdge;
}

[[nodiscard]] NodeMap makeNodeMap() {
    NodeMap result;
    for (std::size_t x = 0; x < ramAirCellNodesPerEdge; ++x) {
        for (std::size_t y = 0; y < ramAirCellNodesPerEdge; ++y) {
            for (std::size_t z = 0; z < ramAirCellNodesPerEdge; ++z) {
                if (isSurfaceNode(x, y, z)) {
                    result.emplace(NodeKey{x, y, z}, result.size());
                }
            }
        }
    }
    return result;
}

[[nodiscard]] std::size_t nodeIndex(
    const NodeMap& nodes,
    const std::size_t x,
    const std::size_t y,
    const std::size_t z) {
    const auto found = nodes.find({x, y, z});
    if (found == nodes.end()) {
        throw std::logic_error(
            "ram-air cell panel references a missing shell node");
    }
    return found->second;
}

[[nodiscard]] std::size_t panelIndex(const PanelKind kind) {
    return static_cast<std::size_t>(kind);
}

[[nodiscard]] std::size_t panelTriangleBegin(const PanelKind kind) {
    return panelIndex(kind) * trianglesPerPanel;
}

void addPanelQuad(StructureDefinition& definition,
                  const std::array<std::size_t, 4>& corners) {
    definition.triangles.push_back(
        {{corners[0], corners[1], corners[2]}});
    definition.triangles.push_back(
        {{corners[0], corners[2], corners[3]}});
}

void addPanelTriangles(StructureDefinition& definition,
                       const NodeMap& nodes,
                       const PanelKind kind) {
    switch (kind) {
    case PanelKind::Back:
        for (std::size_t z = 0; z < ramAirCellTilesPerEdge; ++z) {
            for (std::size_t y = 0; y < ramAirCellTilesPerEdge; ++y) {
                addPanelQuad(definition, {
                    nodeIndex(nodes, 4, y, z),
                    nodeIndex(nodes, 4, y + 1, z),
                    nodeIndex(nodes, 4, y + 1, z + 1),
                    nodeIndex(nodes, 4, y, z + 1),
                });
            }
        }
        return;
    case PanelKind::Left:
    case PanelKind::Right: {
        const std::size_t y = kind == PanelKind::Left ? 0 : 4;
        for (std::size_t z = 0; z < ramAirCellTilesPerEdge; ++z) {
            for (std::size_t x = 0; x < ramAirCellTilesPerEdge; ++x) {
                // The bridge's positive-Y material chart is (z,x).
                addPanelQuad(definition, {
                    nodeIndex(nodes, x, y, z),
                    nodeIndex(nodes, x, y, z + 1),
                    nodeIndex(nodes, x + 1, y, z + 1),
                    nodeIndex(nodes, x + 1, y, z),
                });
            }
        }
        return;
    }
    case PanelKind::Bottom:
    case PanelKind::Top: {
        const std::size_t z = kind == PanelKind::Bottom ? 0 : 4;
        for (std::size_t y = 0; y < ramAirCellTilesPerEdge; ++y) {
            for (std::size_t x = 0; x < ramAirCellTilesPerEdge; ++x) {
                addPanelQuad(definition, {
                    nodeIndex(nodes, x, y, z),
                    nodeIndex(nodes, x + 1, y, z),
                    nodeIndex(nodes, x + 1, y + 1, z),
                    nodeIndex(nodes, x, y + 1, z),
                });
            }
        }
        return;
    }
    }
    throw std::logic_error("ram-air cell panel kind is invalid");
}

[[nodiscard]] StructureMembraneMaterial membraneMaterial() {
    StructureMembraneMaterial material;
    material.warpStiffnessNewtonsPerMeter = 600.0;
    material.weftStiffnessNewtonsPerMeter = 480.0;
    material.couplingStiffnessNewtonsPerMeter = 80.0;
    material.shearStiffnessNewtonsPerMeter = 130.0;
    material.dampingSeconds = 0.012;
    material.compressionStiffnessRatio = 0.08;
    return material;
}

[[nodiscard]] std::array<StructureVector2, 3> intrinsicChart(
    const StructureDefinition& definition,
    const StructureTriangleDefinition& triangle) {
    const auto& first = definition.nodes[triangle.nodes[0]].positionMeters;
    const auto& second = definition.nodes[triangle.nodes[1]].positionMeters;
    const auto& third = definition.nodes[triangle.nodes[2]].positionMeters;
    const StructureVector3 edge = subtract(second, first);
    const StructureVector3 diagonal = subtract(third, first);
    const double edgeLength = length(edge);
    const double thirdX = dot(diagonal, edge) / edgeLength;
    const double thirdYSquared = std::max(
        0.0, dot(diagonal, diagonal) - thirdX * thirdX);
    return {
        StructureVector2{0.0, 0.0},
        StructureVector2{edgeLength, 0.0},
        StructureVector2{thirdX, std::sqrt(thirdYSquared)},
    };
}

struct EdgeIncidence {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t opposite = 0;
};

void addPanelDihedrals(StructureDefinition& definition,
                       const PanelKind kind) {
    using Edge = std::pair<std::size_t, std::size_t>;
    std::map<Edge, std::vector<EdgeIncidence>> incidences;
    const std::size_t begin = panelTriangleBegin(kind);
    const std::size_t end = begin + trianglesPerPanel;
    for (std::size_t index = begin; index < end; ++index) {
        const auto& triangle = definition.triangles[index];
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
                "ram-air cell panel is not an oriented manifold");
        }
        const auto& first = adjacent[0];
        const auto& second = adjacent[1];
        const auto& a = definition.nodes[first.from].positionMeters;
        const auto& b = definition.nodes[first.to].positionMeters;
        const auto& c = definition.nodes[first.opposite].positionMeters;
        const auto& d = definition.nodes[second.opposite].positionMeters;
        const StructureVector3 edgeVector = subtract(b, a);
        const StructureVector3 firstArea = cross(
            edgeVector, subtract(c, a));
        const StructureVector3 secondArea = cross(
            subtract(d, a), edgeVector);
        const StructureVector3 firstNormal = scaled(
            firstArea, 1.0 / length(firstArea));
        const StructureVector3 secondNormal = scaled(
            secondArea, 1.0 / length(secondArea));
        const StructureVector3 edgeUnit = scaled(
            edgeVector, 1.0 / length(edgeVector));
        const double cosine = std::clamp(
            dot(firstNormal, secondNormal), -1.0, 1.0);
        const double sine = dot(
            edgeUnit, cross(firstNormal, secondNormal));
        definition.dihedrals.push_back({
            {first.from, first.to, first.opposite, second.opposite},
            std::atan2(sine, cosine),
            0.35,
        });
    }
}

[[nodiscard]] StructureDefinition makeDefinition() {
    StructureDefinition definition;
    const NodeMap nodes = makeNodeMap();
    definition.nodes.resize(nodes.size());
    for (const auto& [key, index] : nodes) {
        const bool mouthPerimeter = key[1] == 0 || key[1] == 4
            || key[2] == 0 || key[2] == 4;
        definition.nodes[index] = {
            {cellMinimumXMeters + cellSpacingMeters * key[0],
             cellMinimumYMeters + cellSpacingMeters * key[1],
             cellMinimumZMeters + cellSpacingMeters * key[2]},
            0.0,
            key[0] <= 1 && mouthPerimeter,
        };
    }

    definition.triangles.reserve(cellTriangleCount);
    for (const auto& panel : panels) {
        addPanelTriangles(definition, nodes, panel.kind);
    }
    if (definition.triangles.size() != cellTriangleCount) {
        throw std::logic_error(
            "ram-air cell did not construct every panel triangle");
    }

    const auto material = membraneMaterial();
    definition.membranes.reserve(definition.triangles.size());
    for (std::size_t triangleIndex = 0;
         triangleIndex < definition.triangles.size(); ++triangleIndex) {
        const auto& triangle = definition.triangles[triangleIndex];
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
            triangleIndex,
            intrinsicChart(definition, triangle),
            material,
            StructureMaterialRole::Bulk,
        });
    }
    for (const auto& panel : panels) {
        addPanelDihedrals(definition, panel.kind);
    }
    return definition;
}

[[nodiscard]] fluid::PeriodicCartesianGrid makeGrid() {
    return {{gridCellsX, gridCellsY, gridCellsZ},
            {0.0, -1.5, 0.0},
            {4.0, 1.5, 3.0}};
}

[[nodiscard]] std::vector<fluid::GridFaceMovingInterface>
makeInterfaceFaces() {
    std::vector<fluid::GridFaceMovingInterface> result;
    result.reserve(
        ramAirCellPanelCount
        * ramAirCellTilesPerEdge * ramAirCellTilesPerEdge);
    for (const auto& panel : panels) {
        switch (panel.kind) {
        case PanelKind::Back:
            for (std::size_t z = 0; z < 4; ++z) {
                for (std::size_t y = 0; y < 4; ++y) {
                    result.push_back({
                        panel.surfaceStableId, 1, 1,
                        fluid::GridFaceAxis::X,
                        panel.planeIndex, cellLeftY + y, cellBottomZ + z,
                        0.0,
                    });
                }
            }
            break;
        case PanelKind::Left:
        case PanelKind::Right:
            for (std::size_t z = 0; z < 4; ++z) {
                for (std::size_t x = 0; x < 4; ++x) {
                    result.push_back({
                        panel.surfaceStableId, 1, 1,
                        fluid::GridFaceAxis::Y,
                        cellFrontX + x, panel.planeIndex, cellBottomZ + z,
                        0.0,
                    });
                }
            }
            break;
        case PanelKind::Bottom:
        case PanelKind::Top:
            for (std::size_t y = 0; y < 4; ++y) {
                for (std::size_t x = 0; x < 4; ++x) {
                    result.push_back({
                        panel.surfaceStableId, 1, 1,
                        fluid::GridFaceAxis::Z,
                        cellFrontX + x, cellLeftY + y, panel.planeIndex,
                        0.0,
                    });
                }
            }
            break;
        }
    }
    return result;
}

[[nodiscard]] fluid::MovingInterfaceProjectionSettings
makeFluidSettings() {
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.225;
    settings.projection.timeStepSeconds = 1.0 / 120.0;
    settings.projection.absoluteResidualTolerance = 1.0e-9;
    settings.projection.relativeResidualTolerance = 1.0e-11;
    settings.projection.maximumIterations = 4000;
    settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond =
        1.0e-10;
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
            != ramAirCellPanelCount
                * ramAirCellTilesPerEdge * ramAirCellTilesPerEdge) {
        throw std::logic_error(
            "ram-air cell could not construct its reference CFD cavity");
    }
    return diagnostics;
}

[[nodiscard]] std::vector<CouplingSurfaceNodeDefinition>
makeCouplingNodes(const StructureDefinition& definition,
                  const PanelKind kind) {
    std::set<std::size_t> selected;
    const std::size_t begin = panelTriangleBegin(kind);
    for (std::size_t index = begin;
         index < begin + trianglesPerPanel; ++index) {
        selected.insert(definition.triangles[index].nodes.begin(),
                        definition.triangles[index].nodes.end());
    }
    std::vector<CouplingSurfaceNodeDefinition> result;
    result.reserve(selected.size());
    for (const std::size_t node : selected) {
        result.push_back({100'000 + node, node});
    }
    return result;
}

[[nodiscard]] std::vector<CouplingSurfaceTriangleDefinition>
makeCouplingTriangles(const StructureDefinition& definition,
                      const PanelKind kind) {
    std::vector<CouplingSurfaceTriangleDefinition> result;
    result.reserve(trianglesPerPanel);
    const std::size_t begin = panelTriangleBegin(kind);
    for (std::size_t index = begin;
         index < begin + trianglesPerPanel; ++index) {
        const auto& triangle = definition.triangles[index];
        result.push_back({
            200'000 + index,
            {100'000 + triangle.nodes[0],
             100'000 + triangle.nodes[1],
             100'000 + triangle.nodes[2]},
        });
    }
    return result;
}

[[nodiscard]] viewer::StructureFrameMappingDefinition makeFrameMapping() {
    viewer::StructureFrameMappingDefinition mapping;
    const NodeMap nodes = makeNodeMap();
    mapping.vertexStableIds.reserve(nodes.size());
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        mapping.vertexStableIds.push_back(100'000 + index);
    }
    mapping.triangles.reserve(cellTriangleCount);
    for (std::size_t index = 0; index < cellTriangleCount; ++index) {
        mapping.triangles.push_back({200'000 + index, 1, 1});
    }
    return mapping;
}

[[nodiscard]] StructureStepSettings makeStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 120.0;
    settings.substeps = 4;
    settings.constraintIterations = 24;
    settings.cableConstraintSweepPairs = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 1.25;
    settings.workerThreads = 0;
    return settings;
}

[[nodiscard]] std::size_t panelIndexForSurface(
    const std::uint64_t surfaceStableId) {
    for (std::size_t index = 0; index < panels.size(); ++index) {
        if (panels[index].surfaceStableId == surfaceStableId) {
            return index;
        }
    }
    throw std::logic_error(
        "ram-air cell fluid face has an unknown panel ID");
}

[[nodiscard]] std::size_t tileIndex(
    const PanelDescriptor& panel,
    const fluid::MovingInterfaceFaceDiagnostics& face) {
    switch (panel.kind) {
    case PanelKind::Back:
        return (face.k - cellBottomZ) * ramAirCellTilesPerEdge
            + (face.j - cellLeftY);
    case PanelKind::Left:
    case PanelKind::Right:
        return (face.k - cellBottomZ) * ramAirCellTilesPerEdge
            + (face.i - cellFrontX);
    case PanelKind::Bottom:
    case PanelKind::Top:
        return (face.j - cellLeftY) * ramAirCellTilesPerEdge
            + (face.i - cellFrontX);
    }
    throw std::logic_error("ram-air cell panel kind is invalid");
}

[[nodiscard]] double axisComponent(
    const fluid::Vector3& value,
    const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return value.x;
    case fluid::GridFaceAxis::Y: return value.y;
    case fluid::GridFaceAxis::Z: return value.z;
    }
    throw std::logic_error("ram-air cell face axis is invalid");
}

[[nodiscard]] double openingMeanVelocity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    double sum = 0.0;
    for (std::size_t z = 0; z < ramAirCellTilesPerEdge; ++z) {
        for (std::size_t y = 0; y < ramAirCellTilesPerEdge; ++y) {
            sum += velocity.xFaces()[grid.cellIndex(
                cellFrontX, cellLeftY + y, cellBottomZ + z)];
        }
    }
    return sum / static_cast<double>(
        ramAirCellTilesPerEdge * ramAirCellTilesPerEdge);
}

[[nodiscard]] double openingRootMeanSquareVelocity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    double sumSquares = 0.0;
    for (std::size_t z = 0; z < ramAirCellTilesPerEdge; ++z) {
        for (std::size_t y = 0; y < ramAirCellTilesPerEdge; ++y) {
            const double sample = velocity.xFaces()[grid.cellIndex(
                cellFrontX, cellLeftY + y, cellBottomZ + z)];
            sumSquares += sample * sample;
        }
    }
    return std::sqrt(sumSquares / static_cast<double>(
        ramAirCellTilesPerEdge * ramAirCellTilesPerEdge));
}

} // namespace

RamAirCellCase::RamAirCellCase()
    : grid_(makeGrid()),
      interface_(grid_, makeInterfaceFaces()),
      velocity_(grid_),
      pressure_(grid_),
      structure_(makeDefinition()),
      fluidSettings_(makeFluidSettings()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()),
      lastNodalForcesNewtons_(structure_.definition().nodes.size()) {
    const auto reference = makeReferenceFluidDiagnostics(
        grid_, interface_, fluidSettings_);
    bridges_.reserve(ramAirCellPanelCount);
    referenceKinematics_.reserve(ramAirCellPanelCount);
    for (const auto& panel : panels) {
        auto bridge = std::make_unique<
            PlanarFaceResolvedFluidStructureBridge>(
                structure_, panel.surfaceStableId,
                makeCouplingNodes(structure_.definition(), panel.kind),
                makeCouplingTriangles(structure_.definition(), panel.kind),
                reference.faces);
        referenceKinematics_.push_back(
            bridge->transfer().captureKinematics(structure_));
        bridges_.push_back(std::move(bridge));
    }
}

RamAirCellCase::~RamAirCellCase() = default;

viewer::TraceHeader RamAirCellCase::traceHeader() const {
    return {ramAirCellCaseChecksum, ramAirCellCaseSolverId};
}

viewer::DiagnosticFrame RamAirCellCase::advance() {
    const StructureCheckpoint structureBefore = structure_.checkpoint();
    const fluid::MacVelocityField velocityBefore = velocity_;
    const fluid::CellScalarField pressureBefore = pressure_;
    const auto fluidDiagnosticsBefore = fluidDiagnostics_;
    const auto diagnosticsBefore = diagnostics_;
    const auto nodalForcesBefore = lastNodalForcesNewtons_;
    const double gustBefore = gustSpeedMetersPerSecond_;
    try {
        const double nextTime = structure_.simulationTimeSeconds()
            + stepSettings_.timeStepSeconds;
        const double nextGust = gustAmplitudeMetersPerSecond
            * std::sin(2.0 * std::numbers::pi
                       * gustFrequencyHertz * nextTime);
        const double increment = nextGust - gustSpeedMetersPerSecond_;
        for (double& xFace : velocity_.xFaces()) {
            xFace += increment;
        }

        const auto projected = fluid::projectVelocityWithMovingInterfaces(
            grid_, velocity_, pressure_, interface_, fluidSettings_);
        if (!projected.projection.converged || !projected.finite) {
            throw std::runtime_error(
                "ram-air cell CFD projection did not converge");
        }

        RamAirCellDiagnostics nextDiagnostics;
        nextDiagnostics.panels.reserve(bridges_.size());
        std::ranges::fill(
            lastNodalForcesNewtons_, StructureVector3{});
        for (std::size_t panel = 0; panel < bridges_.size(); ++panel) {
            const auto transferred =
                bridges_[panel]->evaluateConstraintReaction(
                    projected, referenceKinematics_[panel]);
            bridges_[panel]->transfer().addLoadsTo(
                structure_, transferred.transferResult());
            for (const auto& load
                 : transferred.transferResult().nodeLoads()) {
                lastNodalForcesNewtons_[load.structureNode] = add(
                    lastNodalForcesNewtons_[load.structureNode],
                    load.forceNewtons);
            }
            const auto& mapping = transferred.diagnostics();
            nextDiagnostics.panels.push_back(mapping);
            nextDiagnostics.fluidPressureForceNewtons = add(
                nextDiagnostics.fluidPressureForceNewtons,
                mapping.fluidPressureForceNewtons);
            nextDiagnostics.fluidReactionForceNewtons = add(
                nextDiagnostics.fluidReactionForceNewtons,
                mapping.fluidLoadForceNewtons);
            nextDiagnostics.transferredForceNewtons = add(
                nextDiagnostics.transferredForceNewtons,
                mapping.structureSurfaceForceNewtons);
            nextDiagnostics.forceResidualNewtons = add(
                nextDiagnostics.forceResidualNewtons,
                mapping.forceResidualNewtons);
            nextDiagnostics.fluidReactionMomentNewtonMeters = add(
                nextDiagnostics.fluidReactionMomentNewtonMeters,
                mapping.fluidLoadMomentNewtonMeters);
            nextDiagnostics.transferredMomentNewtonMeters = add(
                nextDiagnostics.transferredMomentNewtonMeters,
                mapping.structureSurfaceMomentNewtonMeters);
            nextDiagnostics.momentResidualNewtonMeters = add(
                nextDiagnostics.momentResidualNewtonMeters,
                mapping.momentResidualNewtonMeters);
            nextDiagnostics.maximumPanelForceResidualNewtons = std::max(
                nextDiagnostics.maximumPanelForceResidualNewtons,
                mapping.forceResidualNormNewtons);
            nextDiagnostics.maximumPanelMomentResidualNewtonMeters = std::max(
                nextDiagnostics.maximumPanelMomentResidualNewtonMeters,
                mapping.momentResidualNormNewtonMeters);
        }

        const StructureDiagnostics structureDiagnostics =
            structure_.step(stepSettings_);
        if (!structureDiagnostics.finite) {
            throw std::runtime_error(
                "ram-air cell XPBD step produced non-finite diagnostics");
        }
        nextDiagnostics.fluidDivergenceL2PerSecond =
            projected.projection.divergenceL2AfterPerSecond;
        nextDiagnostics.finite =
            finite(nextDiagnostics.fluidPressureForceNewtons)
            && finite(nextDiagnostics.fluidReactionForceNewtons)
            && finite(nextDiagnostics.transferredForceNewtons)
            && finite(nextDiagnostics.forceResidualNewtons)
            && finite(nextDiagnostics.fluidReactionMomentNewtonMeters)
            && finite(nextDiagnostics.transferredMomentNewtonMeters)
            && finite(nextDiagnostics.momentResidualNewtonMeters)
            && std::isfinite(
                nextDiagnostics.maximumPanelForceResidualNewtons)
            && std::isfinite(
                nextDiagnostics.maximumPanelMomentResidualNewtonMeters)
            && std::isfinite(
                nextDiagnostics.fluidDivergenceL2PerSecond);
        if (!nextDiagnostics.finite) {
            throw std::runtime_error(
                "ram-air cell coupling diagnostics are non-finite");
        }

        fluidDiagnostics_ = projected;
        diagnostics_ = std::move(nextDiagnostics);
        gustSpeedMetersPerSecond_ = nextGust;

        viewer::StructureFrameContext context;
        context.sceneChecksum = ramAirCellCaseChecksum;
        context.solverCommit = ramAirCellCaseSolverId;
        context.timeStepSeconds = stepSettings_.timeStepSeconds;
        context.couplingResiduals.tractionNewtons =
            diagnostics_.maximumPanelForceResidualNewtons;
        context.couplingResiduals.structure =
            structureDiagnostics.maximumMembraneResidual;
        context.conservation.interfaceForceResidualNewtons = {
            diagnostics_.forceResidualNewtons.x,
            diagnostics_.forceResidualNewtons.y,
            diagnostics_.forceResidualNewtons.z,
        };
        context.conservation.interfaceMomentResidualNewtonMetres = {
            diagnostics_.momentResidualNewtonMeters.x,
            diagnostics_.momentResidualNewtonMeters.y,
            diagnostics_.momentResidualNewtonMeters.z,
        };
        viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
            structure_, frameMapping_, context);

        const auto states = structure_.nodeStates();
        const auto& definition = structure_.definition();
        std::vector<double> displacements;
        displacements.reserve(states.size());
        for (std::size_t index = 0; index < states.size(); ++index) {
            displacements.push_back(length(subtract(
                states[index].positionMeters,
                definition.nodes[index].positionMeters)));
        }
        frame.scalarFields.push_back({
            "cell.displacement", "m", viewer::FieldAssociation::Vertex,
            std::move(displacements),
        });

        std::vector<double> panelIds(cellTriangleCount, 0.0);
        std::vector<double> pressureTractions(cellTriangleCount, 0.0);
        std::vector<double> reactionTractions(cellTriangleCount, 0.0);
        for (std::size_t panel = 0; panel < panels.size(); ++panel) {
            const std::size_t begin = panel * trianglesPerPanel;
            std::fill_n(
                panelIds.begin() + static_cast<std::ptrdiff_t>(begin),
                trianglesPerPanel, static_cast<double>(panel + 1));
        }
        for (const auto& face : fluidDiagnostics_.faces) {
            const std::size_t panel = panelIndexForSurface(
                face.surfaceStableId);
            const std::size_t tile = tileIndex(panels[panel], face);
            if (tile >= ramAirCellTilesPerEdge
                            * ramAirCellTilesPerEdge) {
                throw std::logic_error(
                    "ram-air cell CFD tile escaped its structural panel");
            }
            const std::size_t triangle =
                panelTriangleBegin(panels[panel].kind) + 2 * tile;
            for (const std::size_t offset : {std::size_t{0}, std::size_t{1}}) {
                pressureTractions[triangle + offset] = axisComponent(
                    face.pressureTractionPascals, face.axis);
                reactionTractions[triangle + offset] = axisComponent(
                    face.constraintReactionTractionPascals, face.axis);
            }
        }
        frame.scalarFields.push_back({
            "cell.panel", "1", viewer::FieldAssociation::Triangle,
            std::move(panelIds),
        });
        frame.scalarFields.push_back({
            "cell.cfd_pressure_traction", "Pa",
            viewer::FieldAssociation::Triangle,
            std::move(pressureTractions),
        });
        frame.scalarFields.push_back({
            "cell.cfd_complete_reaction_traction", "Pa",
            viewer::FieldAssociation::Triangle,
            std::move(reactionTractions),
        });
        frame.scalarFields.push_back({
            "cell.gust_speed", "m/s", viewer::FieldAssociation::Global,
            {gustSpeedMetersPerSecond_},
        });
        frame.scalarFields.push_back({
            "cell.opening_mean_velocity", "m/s",
            viewer::FieldAssociation::Global,
            {openingMeanVelocity(grid_, velocity_)},
        });
        frame.scalarFields.push_back({
            "cell.opening_rms_velocity", "m/s",
            viewer::FieldAssociation::Global,
            {openingRootMeanSquareVelocity(grid_, velocity_)},
        });
        frame.scalarFields.push_back({
            "cell.fluid_divergence_l2", "1/s",
            viewer::FieldAssociation::Global,
            {diagnostics_.fluidDivergenceL2PerSecond},
        });

        std::vector<viewer::Vec3d> nodalForces;
        nodalForces.reserve(lastNodalForcesNewtons_.size());
        for (const auto& force : lastNodalForcesNewtons_) {
            nodalForces.push_back({force.x, force.y, force.z});
        }
        frame.vectorFields.push_back({
            "cell.cfd_nodal_force", "N",
            viewer::FieldAssociation::Vertex,
            std::move(nodalForces),
        });
        frame.vectorFields.push_back({
            "cell.total_cfd_reaction", "N",
            viewer::FieldAssociation::Global,
            {{diagnostics_.fluidReactionForceNewtons.x,
              diagnostics_.fluidReactionForceNewtons.y,
              diagnostics_.fluidReactionForceNewtons.z}},
        });

        viewer::ProtocolError error;
        if (!viewer::validateFrame(frame, &error)) {
            throw std::logic_error(
                "ram-air cell frame is invalid: " + error.message);
        }
        return frame;
    } catch (...) {
        structure_.restore(structureBefore);
        velocity_ = velocityBefore;
        pressure_ = pressureBefore;
        fluidDiagnostics_ = fluidDiagnosticsBefore;
        diagnostics_ = diagnosticsBefore;
        lastNodalForcesNewtons_ = nodalForcesBefore;
        gustSpeedMetersPerSecond_ = gustBefore;
        throw;
    }
}

const Structure& RamAirCellCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& RamAirCellCase::stepSettings() const noexcept {
    return stepSettings_;
}

const fluid::MacVelocityField& RamAirCellCase::velocity() const noexcept {
    return velocity_;
}

const fluid::CellScalarField& RamAirCellCase::pressure() const noexcept {
    return pressure_;
}

const fluid::MovingInterfaceProjectionDiagnostics&
RamAirCellCase::fluidDiagnostics() const noexcept {
    return fluidDiagnostics_;
}

const RamAirCellDiagnostics& RamAirCellCase::diagnostics() const noexcept {
    return diagnostics_;
}

double RamAirCellCase::gustSpeedMetersPerSecond() const noexcept {
    return gustSpeedMetersPerSecond_;
}

double RamAirCellCase::openingMeanVelocityMetersPerSecond() const {
    return openingMeanVelocity(grid_, velocity_);
}

double RamAirCellCase::openingRmsVelocityMetersPerSecond() const {
    return openingRootMeanSquareVelocity(grid_, velocity_);
}

double RamAirCellCase::maximumDisplacementMeters() const {
    const auto states = structure_.nodeStates();
    const auto& definition = structure_.definition();
    double maximum = 0.0;
    for (std::size_t index = 0; index < states.size(); ++index) {
        maximum = std::max(maximum, length(subtract(
            states[index].positionMeters,
            definition.nodes[index].positionMeters)));
    }
    return maximum;
}

double RamAirCellCase::maximumOutwardInflationMeters() const {
    const auto states = structure_.nodeStates();
    const auto& definition = structure_.definition();
    double maximum = 0.0;
    for (std::size_t index = 0; index < states.size(); ++index) {
        const auto& reference = definition.nodes[index].positionMeters;
        const auto& current = states[index].positionMeters;
        if (reference.x == cellMinimumXMeters + 1.0) {
            maximum = std::max(maximum, current.x - reference.x);
        }
        if (reference.y == cellMinimumYMeters) {
            maximum = std::max(maximum, reference.y - current.y);
        }
        if (reference.y == cellMinimumYMeters + 1.0) {
            maximum = std::max(maximum, current.y - reference.y);
        }
        if (reference.z == cellMinimumZMeters) {
            maximum = std::max(maximum, reference.z - current.z);
        }
        if (reference.z == cellMinimumZMeters + 1.0) {
            maximum = std::max(maximum, current.z - reference.z);
        }
    }
    return maximum;
}

} // namespace simwing::fsi
