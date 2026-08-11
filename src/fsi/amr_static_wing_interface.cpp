#include "amr_static_wing_interface.h"

#include "fluid/scene_surface_clipping.h"
#include "fluid/scene_surface_intersection.h"
#include "scene_fluid_surface.h"
#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace simwing::fsi::amr {
namespace {

struct LevelGeometry {
    StaticWingInterfaceLevel level;
    fluid::SceneFluidGridCandidateSet candidates;
    fluid::SceneFluidGridPatchSet patches;
};

std::size_t checkedCellCount(const fluid::GridCellCounts counts) {
    if (counts.x == 0 || counts.y == 0 || counts.z == 0
        || counts.x > std::numeric_limits<std::size_t>::max() / counts.y
        || counts.x * counts.y
               > std::numeric_limits<std::size_t>::max() / counts.z) {
        throw std::overflow_error(
            "static wing interface grid cell count overflowed");
    }
    return counts.x * counts.y * counts.z;
}

std::size_t cellIndex(const fluid::GridCellCounts counts,
                      const std::size_t i,
                      const std::size_t j,
                      const std::size_t k) {
    return i + counts.x * (j + counts.y * k);
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

double length(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

std::array<double, 6> triangleNormalTensor(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const std::size_t triangleIndex) {
    const auto& triangle = surface.triangles.at(triangleIndex);
    const Vec3 first = state.vertices.at(triangle.vertexIndices[0]).positionMeters;
    const Vec3 second = state.vertices.at(triangle.vertexIndices[1]).positionMeters;
    const Vec3 third = state.vertices.at(triangle.vertexIndices[2]).positionMeters;
    Vec3 normal = cross(subtract(second, first), subtract(third, first));
    const double magnitude = length(normal);
    if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
        throw std::invalid_argument(
            "static wing interface contains a degenerate triangle");
    }
    normal.x /= magnitude;
    normal.y /= magnitude;
    normal.z /= magnitude;
    return {
        normal.x * normal.x,
        normal.x * normal.y,
        normal.x * normal.z,
        normal.y * normal.y,
        normal.y * normal.z,
        normal.z * normal.z,
    };
}

void addCellPatch(StaticWingInterfaceCell& cell,
                  const std::array<double, 6>& tensor,
                  const double areaSquareMeters) {
    if (!std::isfinite(areaSquareMeters) || !(areaSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "static wing interface patch area is invalid");
    }
    for (std::size_t component = 0; component < tensor.size(); ++component) {
        cell.normalProjectionTensor[component] +=
            areaSquareMeters * tensor[component];
    }
    cell.surfaceAreaSquareMeters += areaSquareMeters;
}

LevelGeometry buildLevel(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::GridCellCounts counts,
    const fluid::Vector3 lower,
    const fluid::Vector3 upper) {
    const fluid::PeriodicCartesianGrid grid(counts, lower, upper);
    LevelGeometry result;
    result.level.cellCounts = counts;
    result.level.cells.resize(checkedCellCount(counts));
    result.candidates = fluid::buildSceneFluidGridCandidates(
        surface, state, grid);
    const auto intersections = fluid::intersectSceneFluidSurfaceWithGrid(
        surface, state, grid, result.candidates);
    result.patches = fluid::clipSceneFluidSurfaceToCells(
        surface, state, grid, result.candidates, intersections);
    std::vector<std::array<double, 6>> triangleTensors;
    triangleTensors.reserve(surface.triangles.size());
    for (std::size_t triangle = 0;
         triangle < surface.triangles.size(); ++triangle) {
        triangleTensors.push_back(
            triangleNormalTensor(surface, state, triangle));
    }
    for (const auto& patch : result.patches.patches) {
        if (patch.dimension != fluid::SceneFluidPatchDimension::Area
            || !(patch.areaSquareMeters > 0.0)) {
            continue;
        }
        addCellPatch(
            result.level.cells.at(patch.cellIndex),
            triangleTensors.at(patch.triangleIndex),
            patch.areaSquareMeters);
    }
    for (auto& cell : result.level.cells) {
        if (!cell.active()) {
            continue;
        }
        for (double& component : cell.normalProjectionTensor) {
            component /= cell.surfaceAreaSquareMeters;
        }
    }
    return result;
}

bool inRefinedPatch(const WindTunnelGridSettings& settings,
                    const fluid::GridCellCounts levelCounts,
                    const std::size_t i,
                    const std::size_t j,
                    const std::size_t k) {
    const bool fine = levelCounts != settings.coarseCellCounts;
    const std::size_t scale = fine ? settings.refinementRatio : 1;
    return i >= scale * settings.refinedCoarseCellLower[0]
        && i < scale * settings.refinedCoarseCellUpperExclusive[0]
        && j >= scale * settings.refinedCoarseCellLower[1]
        && j < scale * settings.refinedCoarseCellUpperExclusive[1]
        && k >= scale * settings.refinedCoarseCellLower[2]
        && k < scale * settings.refinedCoarseCellUpperExclusive[2];
}

std::size_t activeCellCount(const StaticWingInterfaceLevel& level) {
    return static_cast<std::size_t>(std::ranges::count_if(
        level.cells, [](const StaticWingInterfaceCell& cell) {
            return cell.active();
        }));
}

} // namespace

StaticWingInterface::StaticWingInterface(
    const Scene& scene,
    const WindTunnelGridSettings& gridSettings) {
    const SceneFluidSurfaceAssembly surface =
        assembleSceneFluidSurface(scene);
    if (!surface.ok()) {
        const std::string detail = surface.diagnostics.empty()
            ? std::string("unknown surface assembly error")
            : surface.diagnostics.front().message;
        throw std::invalid_argument(
            "static wing scene surface assembly failed: " + detail);
    }
    const SceneStructureAssembly structureAssembly =
        assembleSceneStructure(scene);
    if (!structureAssembly.ok()) {
        const std::string detail = structureAssembly.diagnostics.empty()
            ? std::string("unknown structure assembly error")
            : structureAssembly.diagnostics.front().message;
        throw std::invalid_argument(
            "static wing structure assembly failed: " + detail);
    }
    Structure structure(structureAssembly.definition);
    const SceneFluidSurfaceState state = captureSceneFluidSurfaceState(
        surface.definition, structureAssembly.mappings, structure);

    if (state.vertices.empty() || surface.definition.triangles.empty()) {
        throw std::invalid_argument(
            "static wing interface requires surface geometry");
    }
    diagnostics_.lowerMeters = state.vertices.front().positionMeters;
    diagnostics_.upperMeters = diagnostics_.lowerMeters;
    vertices_.reserve(state.vertices.size());
    for (const auto& vertex : state.vertices) {
        diagnostics_.lowerMeters.x = std::min(
            diagnostics_.lowerMeters.x, vertex.positionMeters.x);
        diagnostics_.lowerMeters.y = std::min(
            diagnostics_.lowerMeters.y, vertex.positionMeters.y);
        diagnostics_.lowerMeters.z = std::min(
            diagnostics_.lowerMeters.z, vertex.positionMeters.z);
        diagnostics_.upperMeters.x = std::max(
            diagnostics_.upperMeters.x, vertex.positionMeters.x);
        diagnostics_.upperMeters.y = std::max(
            diagnostics_.upperMeters.y, vertex.positionMeters.y);
        diagnostics_.upperMeters.z = std::max(
            diagnostics_.upperMeters.z, vertex.positionMeters.z);
        vertices_.push_back({vertex.id, vertex.positionMeters});
    }
    triangles_.reserve(surface.definition.triangles.size());
    for (const auto& triangle : surface.definition.triangles) {
        triangles_.push_back({
            triangle.id,
            triangle.vertexIndices,
            surface.definition.regions.at(
                triangle.negativeSideRegionIndex).id,
            surface.definition.regions.at(
                triangle.positiveSideRegionIndex).id,
            triangle.role,
        });
    }

    const bool inside =
        diagnostics_.lowerMeters.x >= gridSettings.lowerMeters.x
        && diagnostics_.lowerMeters.y >= gridSettings.lowerMeters.y
        && diagnostics_.lowerMeters.z >= gridSettings.lowerMeters.z
        && diagnostics_.upperMeters.x <= gridSettings.upperMeters.x
        && diagnostics_.upperMeters.y <= gridSettings.upperMeters.y
        && diagnostics_.upperMeters.z <= gridSettings.upperMeters.z;
    if (!inside) {
        std::ostringstream message;
        message.precision(17);
        message << "static wing surface bounds ["
                << diagnostics_.lowerMeters.x << ", "
                << diagnostics_.lowerMeters.y << ", "
                << diagnostics_.lowerMeters.z << "] to ["
                << diagnostics_.upperMeters.x << ", "
                << diagnostics_.upperMeters.y << ", "
                << diagnostics_.upperMeters.z
                << "] lie outside wind-tunnel bounds ["
                << gridSettings.lowerMeters.x << ", "
                << gridSettings.lowerMeters.y << ", "
                << gridSettings.lowerMeters.z << "] to ["
                << gridSettings.upperMeters.x << ", "
                << gridSettings.upperMeters.y << ", "
                << gridSettings.upperMeters.z << "]";
        throw std::invalid_argument(message.str());
    }

    const fluid::GridCellCounts fineCounts{
        gridSettings.coarseCellCounts.x * gridSettings.refinementRatio,
        gridSettings.coarseCellCounts.y * gridSettings.refinementRatio,
        gridSettings.coarseCellCounts.z * gridSettings.refinementRatio,
    };
    LevelGeometry coarse = buildLevel(
        surface.definition, state, gridSettings.coarseCellCounts,
        gridSettings.lowerMeters, gridSettings.upperMeters);
    LevelGeometry fine = buildLevel(
        surface.definition, state, fineCounts,
        gridSettings.lowerMeters, gridSettings.upperMeters);
    levels_[0] = std::move(coarse.level);
    levels_[1] = std::move(fine.level);

    diagnostics_.vertexCount = vertices_.size();
    diagnostics_.triangleCount = triangles_.size();
    diagnostics_.coarseCandidateCount = coarse.candidates.candidates.size();
    diagnostics_.refinedCandidateCount = fine.candidates.candidates.size();
    diagnostics_.coarseCutCellCount = activeCellCount(levels_[0]);
    diagnostics_.refinedCutCellCount = activeCellCount(levels_[1]);
    for (std::size_t triangleIndex = 0;
         triangleIndex < surface.definition.triangles.size();
         ++triangleIndex) {
        const auto& triangle = surface.definition.triangles[triangleIndex];
        const Vec3 first = state.vertices[
            triangle.vertexIndices[0]].positionMeters;
        const Vec3 second = state.vertices[
            triangle.vertexIndices[1]].positionMeters;
        const Vec3 third = state.vertices[
            triangle.vertexIndices[2]].positionMeters;
        diagnostics_.surfaceAreaSquareMeters += 0.5 * length(
            cross(subtract(second, first), subtract(third, first)));
    }

    const fluid::Vector3 coarseSpacing{
        (gridSettings.upperMeters.x - gridSettings.lowerMeters.x)
            / static_cast<double>(gridSettings.coarseCellCounts.x),
        (gridSettings.upperMeters.y - gridSettings.lowerMeters.y)
            / static_cast<double>(gridSettings.coarseCellCounts.y),
        (gridSettings.upperMeters.z - gridSettings.lowerMeters.z)
            / static_cast<double>(gridSettings.coarseCellCounts.z),
    };
    const Vec3 refinedLower{
        gridSettings.lowerMeters.x
            + static_cast<double>(
                  gridSettings.refinedCoarseCellLower[0])
                * coarseSpacing.x,
        gridSettings.lowerMeters.y
            + static_cast<double>(
                  gridSettings.refinedCoarseCellLower[1])
                * coarseSpacing.y,
        gridSettings.lowerMeters.z
            + static_cast<double>(
                  gridSettings.refinedCoarseCellLower[2])
                * coarseSpacing.z,
    };
    const Vec3 refinedUpper{
        gridSettings.lowerMeters.x
            + static_cast<double>(
                  gridSettings.refinedCoarseCellUpperExclusive[0])
                * coarseSpacing.x,
        gridSettings.lowerMeters.y
            + static_cast<double>(
                  gridSettings.refinedCoarseCellUpperExclusive[1])
                * coarseSpacing.y,
        gridSettings.lowerMeters.z
            + static_cast<double>(
                  gridSettings.refinedCoarseCellUpperExclusive[2])
                * coarseSpacing.z,
    };
    for (const auto& triangle : surface.definition.triangles) {
        const std::array<Vec3, 3> positions{
            state.vertices[triangle.vertexIndices[0]].positionMeters,
            state.vertices[triangle.vertexIndices[1]].positionMeters,
            state.vertices[triangle.vertexIndices[2]].positionMeters,
        };
        const auto clipped = fluid::clipSceneFluidTriangleToAxisAlignedBox(
            positions, refinedLower, refinedUpper);
        if (clipped
            && clipped->dimension == fluid::SceneFluidPatchDimension::Area) {
            diagnostics_.refinedPatchSurfaceAreaSquareMeters +=
                clipped->areaSquareMeters;
        }
    }

    for (std::size_t k = 0; k < fineCounts.z; ++k) {
        for (std::size_t j = 0; j < fineCounts.y; ++j) {
            for (std::size_t i = 0; i < fineCounts.x; ++i) {
                if (!inRefinedPatch(
                        gridSettings, fineCounts,
                        i, j, k)) {
                    continue;
                }
                const auto& cell = levels_[1].cells[
                    cellIndex(fineCounts, i, j, k)];
                if (cell.active()) {
                    ++diagnostics_.activeCompositeCutCellCount;
                }
            }
        }
    }
    const auto coarseCounts = gridSettings.coarseCellCounts;
    for (std::size_t k = 0; k < coarseCounts.z; ++k) {
        for (std::size_t j = 0; j < coarseCounts.y; ++j) {
            for (std::size_t i = 0; i < coarseCounts.x; ++i) {
                if (inRefinedPatch(
                        gridSettings, coarseCounts, i, j, k)) {
                    continue;
                }
                if (levels_[0].cells[
                        cellIndex(coarseCounts, i, j, k)].active()) {
                    ++diagnostics_.activeCompositeCutCellCount;
                }
            }
        }
    }

    diagnostics_.entirelyInsideWindTunnel =
        diagnostics_.lowerMeters.x >= gridSettings.lowerMeters.x
        && diagnostics_.lowerMeters.y >= gridSettings.lowerMeters.y
        && diagnostics_.lowerMeters.z >= gridSettings.lowerMeters.z
        && diagnostics_.upperMeters.x <= gridSettings.upperMeters.x
        && diagnostics_.upperMeters.y <= gridSettings.upperMeters.y
        && diagnostics_.upperMeters.z <= gridSettings.upperMeters.z;
    diagnostics_.finite =
        std::isfinite(diagnostics_.surfaceAreaSquareMeters)
        && std::isfinite(
            diagnostics_.refinedPatchSurfaceAreaSquareMeters)
        && std::isfinite(diagnostics_.lowerMeters.x)
        && std::isfinite(diagnostics_.lowerMeters.y)
        && std::isfinite(diagnostics_.lowerMeters.z)
        && std::isfinite(diagnostics_.upperMeters.x)
        && std::isfinite(diagnostics_.upperMeters.y)
        && std::isfinite(diagnostics_.upperMeters.z);
    diagnostics_.accepted = diagnostics_.finite
        && diagnostics_.entirelyInsideWindTunnel
        && diagnostics_.vertexCount > 0
        && diagnostics_.triangleCount > 0
        && diagnostics_.coarseCandidateCount > 0
        && diagnostics_.refinedCandidateCount > 0
        && diagnostics_.coarseCutCellCount > 0
        && diagnostics_.refinedCutCellCount > 0
        && diagnostics_.activeCompositeCutCellCount > 0
        && diagnostics_.surfaceAreaSquareMeters > 0.0;
    if (!diagnostics_.accepted) {
        throw std::runtime_error(
            "static wing interface grid binding was rejected");
    }
}

const StaticWingInterfaceDiagnostics&
StaticWingInterface::diagnostics() const noexcept {
    return diagnostics_;
}

const StaticWingInterfaceLevel& StaticWingInterface::level(
    const std::size_t levelIndex) const {
    return levels_.at(levelIndex);
}

const std::vector<StaticWingInterfaceVertex>&
StaticWingInterface::vertices() const noexcept {
    return vertices_;
}

const std::vector<StaticWingInterfaceTriangle>&
StaticWingInterface::triangles() const noexcept {
    return triangles_;
}

} // namespace simwing::fsi::amr
