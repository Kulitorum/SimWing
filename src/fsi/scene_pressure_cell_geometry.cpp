#include "scene_pressure_cell_geometry.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace simwing::fsi {
namespace {

constexpr StableId apexVertexId = 10;

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

Scene makeGeometry(const std::array<Vec3, 4>& positions,
                   const char* designChecksum) {
    Scene scene;
    scene.metadata.designChecksum = designChecksum;
    scene.metadata.exporterVersion = scenePressureCellCaseSolverId;
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "soft diagnostic fabric", 8.0, 6.0, 2.0, 0.0002,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
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

} // namespace

Scene makeScenePressureCellGeometry() {
    const std::array<Vec3, 4> positions{{
        {1.2, 1.5, 1.45},
        {2.0, 1.2, 1.15},
        {2.0, 1.8, 1.15},
        {2.0, 1.5, 1.75},
    }};
    return makeGeometry(positions, scenePressureCellCaseChecksum);
}

Scene makeScenePressureCellRefinementGeometry() {
    const std::array<Vec3, 4> positions{{
        {1.173, 1.483, 1.437},
        {2.613, 1.187, 1.123},
        {2.817, 1.817, 1.169},
        {2.697, 1.503, 1.793},
    }};
    return makeGeometry(positions, scenePressureCellRefinementChecksum);
}

SceneStructureAssembly makeScenePressureCellAssembly(const Scene& scene) {
    auto assembly = assembleSceneStructure(scene);
    for (std::size_t node = 0;
         node < assembly.mappings.nodeVertexIds.size(); ++node) {
        if (assembly.mappings.nodeVertexIds[node] != apexVertexId) {
            assembly.definition.nodes[node].fixed = true;
        }
    }
    return assembly;
}

fluid::PeriodicCartesianGrid makeScenePressureCellGrid(
    const fluid::GridCellCounts cellCounts,
    const fluid::Vector3 lowerMeters) {
    return {
        cellCounts,
        lowerMeters,
        {lowerMeters.x + 4.0,
         lowerMeters.y + 4.0,
         lowerMeters.z + 4.0},
    };
}

} // namespace simwing::fsi
