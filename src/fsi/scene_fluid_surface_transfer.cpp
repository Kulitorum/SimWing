#include "scene_fluid_surface_transfer.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

bool finite(const Vec3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

} // namespace

SceneFluidSurfaceTransfer::Topology SceneFluidSurfaceTransfer::makeTopology(
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const Structure& target) {
    // This performs the complete scene-surface/reference-geometry/winding
    // binding before any coupling topology is published.
    static_cast<void>(captureSceneFluidSurfaceState(
        surface, structureMappings, target));

    SceneFluidSurfaceTransfer::Topology result;
    result.nodes.reserve(surface.vertices.size());
    for (const auto& vertex : surface.vertices) {
        const auto node = structureMappings.nodeIndex(vertex.id);
        if (!node) {
            throw std::invalid_argument(
                "scene fluid coupling vertex has no Structure node");
        }
        result.nodes.push_back({vertex.id, *node});
    }
    result.triangles.reserve(surface.triangles.size());
    for (const auto& triangle : surface.triangles) {
        CouplingSurfaceTriangleDefinition converted;
        converted.stableId = triangle.id;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            if (triangle.vertexIndices[corner] >= surface.vertices.size()) {
                throw std::invalid_argument(
                    "scene fluid coupling triangle vertex is out of range");
            }
            converted.nodeStableIds[corner] =
                surface.vertices[triangle.vertexIndices[corner]].id;
        }
        result.triangles.push_back(std::move(converted));
    }
    return result;
}

SceneFluidSurfaceTransfer::SceneFluidSurfaceTransfer(
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const Structure& target)
    : SceneFluidSurfaceTransfer(
          surface.fingerprint,
          target,
          makeTopology(surface, structureMappings, target)) {}

SceneFluidSurfaceTransfer::SceneFluidSurfaceTransfer(
    const std::uint64_t surfaceDefinitionFingerprint,
    const Structure& target,
    Topology topology)
    : surfaceDefinitionFingerprint_(surfaceDefinitionFingerprint),
      transfer_(target,
                std::move(topology.nodes),
                std::move(topology.triangles)) {}

std::uint64_t
SceneFluidSurfaceTransfer::surfaceDefinitionFingerprint() const noexcept {
    return surfaceDefinitionFingerprint_;
}

std::uint64_t
SceneFluidSurfaceTransfer::couplingSurfaceFingerprint() const noexcept {
    return transfer_.fingerprint();
}

std::uint64_t
SceneFluidSurfaceTransfer::targetDefinitionFingerprint() const noexcept {
    return transfer_.targetDefinitionFingerprint();
}

std::span<const CouplingSurfaceNodeDefinition>
SceneFluidSurfaceTransfer::nodes() const noexcept {
    return transfer_.nodes();
}

std::span<const CouplingSurfaceTriangleDefinition>
SceneFluidSurfaceTransfer::triangles() const noexcept {
    return transfer_.triangles();
}

std::vector<CouplingNodeKinematics> SceneFluidSurfaceTransfer::kinematics(
    const SceneFluidSurfaceState& state) const {
    if (state.version != sceneFluidSurfaceStateVersion
        || state.definitionFingerprint != surfaceDefinitionFingerprint_
        || state.structureDefinitionFingerprint
            != transfer_.targetDefinitionFingerprint()
        || !std::isfinite(state.simulationTimeSeconds)
        || state.simulationTimeSeconds < 0.0
        || state.vertices.size() != transfer_.nodes().size()) {
        throw std::invalid_argument(
            "scene fluid coupling state identity is invalid");
    }

    std::vector<CouplingNodeKinematics> result;
    result.reserve(state.vertices.size());
    for (std::size_t index = 0; index < state.vertices.size(); ++index) {
        const auto& source = state.vertices[index];
        if (source.id != transfer_.nodes()[index].stableId
            || !finite(source.positionMeters)
            || !finite(source.velocityMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid coupling kinematics are non-finite or out of order");
        }
        result.push_back({
            source.id,
            {source.positionMeters.x,
             source.positionMeters.y,
             source.positionMeters.z},
            {source.velocityMetersPerSecond.x,
             source.velocityMetersPerSecond.y,
             source.velocityMetersPerSecond.z},
        });
    }
    return result;
}

ConservativeTransferResult SceneFluidSurfaceTransfer::evaluate(
    const SceneFluidSurfaceState& state,
    const std::span<const CouplingTriangleTraction> triangleTractions,
    const ConservativeTransferSettings& settings) const {
    return transfer_.evaluate(
        kinematics(state), triangleTractions, settings);
}

ConservativeTransferResult SceneFluidSurfaceTransfer::evaluateQuadrature(
    const SceneFluidSurfaceState& state,
    const std::span<const CouplingTriangleTractionQuadrature> quadrature,
    const ConservativeTransferSettings& settings) const {
    return transfer_.evaluateQuadrature(
        kinematics(state), quadrature, settings);
}

void SceneFluidSurfaceTransfer::addLoadsTo(
    Structure& target,
    const ConservativeTransferResult& result) const {
    transfer_.addLoadsTo(target, result);
}

} // namespace simwing::fsi
