#include "scene_fluid_pressure_traction.h"

#include <cmath>
#include <stdexcept>

namespace simwing::fsi {
namespace {

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

Vec3 cross(const Vec3& first, const Vec3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

StructureVector3 unitNormal(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const std::size_t triangleIndex) {
    const auto& triangle = surface.triangles[triangleIndex];
    const Vec3& first =
        state.vertices[triangle.vertexIndices[0]].positionMeters;
    const Vec3& second =
        state.vertices[triangle.vertexIndices[1]].positionMeters;
    const Vec3& third =
        state.vertices[triangle.vertexIndices[2]].positionMeters;
    const Vec3 normal = cross(subtract(second, first), subtract(third, first));
    const double magnitude = std::hypot(normal.x, normal.y, normal.z);
    if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
        throw std::invalid_argument(
            "scene fluid pressure traction has a degenerate current triangle");
    }
    return {normal.x / magnitude,
            normal.y / magnitude,
            normal.z / magnitude};
}

} // namespace

std::vector<SceneFluidQuadratureTraction> buildSceneFluidPressureTractions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& quadrature,
    const std::span<const SceneFluidQuadraturePressure> pressures) {
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidQuadratureDefinition(quadrature);
    if (quadrature.surfaceDefinitionFingerprint != surface.fingerprint
        || quadrature.surfaceStateFingerprint != state.fingerprint
        || quadrature.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || quadrature.acceptedStepCount != state.acceptedStepCount
        || quadrature.simulationTimeSeconds != state.simulationTimeSeconds
        || pressures.size() != quadrature.points.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-traction binding is invalid");
    }

    std::vector<SceneFluidQuadratureTraction> result;
    result.reserve(pressures.size());
    StableId previousTriangleId = invalidStableId;
    StructureVector3 normal;
    for (std::size_t index = 0; index < pressures.size(); ++index) {
        const auto& point = quadrature.points[index];
        const auto& pressure = pressures[index];
        if (pressure.stableId != point.stableId
            || !std::isfinite(pressure.negativeSidePressurePascals)
            || !std::isfinite(pressure.positiveSidePressurePascals)) {
            throw std::invalid_argument(
                "scene fluid pressure samples are non-finite or out of order");
        }
        if (point.triangleId != previousTriangleId) {
            const auto triangleIndex =
                surface.mappings.triangleIndex(point.triangleId);
            if (!triangleIndex) {
                throw std::invalid_argument(
                    "scene fluid pressure sample references a foreign triangle");
            }
            normal = unitNormal(surface, state, *triangleIndex);
            previousTriangleId = point.triangleId;
        }
        const double pressureDifference =
            pressure.negativeSidePressurePascals
            - pressure.positiveSidePressurePascals;
        const StructureVector3 traction{
            pressureDifference * normal.x,
            pressureDifference * normal.y,
            pressureDifference * normal.z,
        };
        if (!std::isfinite(traction.x)
            || !std::isfinite(traction.y)
            || !std::isfinite(traction.z)) {
            throw std::overflow_error(
                "scene fluid pressure traction is not finite");
        }
        result.push_back({point.stableId, traction});
    }
    return result;
}

ConservativeTransferResult evaluateSceneFluidPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const std::span<const SceneFluidQuadraturePressure> pressures,
    const ConservativeTransferSettings& settings) {
    const auto tractions = buildSceneFluidPressureTractions(
        surface, state, quadrature, pressures);
    return evaluateSceneFluidQuadrature(
        transfer, state, quadrature, tractions, settings);
}

} // namespace simwing::fsi
