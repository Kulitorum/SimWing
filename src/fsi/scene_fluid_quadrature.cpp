#include "scene_fluid_quadrature.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    template<typename Enum>
    void enumeration(const Enum value) {
        using Underlying = std::underlying_type_t<Enum>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

std::uint64_t pointStableId(
    const SceneFluidQuadratureOwnerKind kind,
    const StableId triangleId,
    const std::array<std::uint64_t, 4>& ownerCoordinates) {
    Fingerprint fingerprint;
    fingerprint.enumeration(kind);
    fingerprint.integer(triangleId);
    for (const std::uint64_t coordinate : ownerCoordinates) {
        fingerprint.integer(coordinate);
    }
    return fingerprint.value();
}

std::uint64_t definitionFingerprint(
    const SceneFluidQuadratureDefinition& definition) {
    Fingerprint fingerprint;
    fingerprint.integer(definition.version);
    fingerprint.integer(definition.surfaceDefinitionFingerprint);
    fingerprint.integer(definition.surfaceStateFingerprint);
    fingerprint.integer(definition.patchOwnershipFingerprint);
    fingerprint.integer(definition.couplingSurfaceFingerprint);
    fingerprint.integer(definition.structureDefinitionFingerprint);
    fingerprint.integer(definition.acceptedStepCount);
    fingerprint.real(definition.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(definition.points.size()));
    for (const auto& point : definition.points) {
        fingerprint.integer(point.stableId);
        fingerprint.integer(point.triangleId);
        for (const double coordinate : point.barycentricCoordinates) {
            fingerprint.real(coordinate);
        }
        fingerprint.real(point.areaSquareMeters);
        fingerprint.integer(point.negativeSideRegionId);
        fingerprint.integer(point.positiveSideRegionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            point.negativeSideCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            point.positiveSideCellIndex));
        fingerprint.integer(point.materialId);
        fingerprint.integer(point.sheetId);
        fingerprint.enumeration(point.role);
        fingerprint.enumeration(point.ownerKind);
    }
    return fingerprint.value();
}

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

SceneFluidQuadraturePoint pointFromPatch(
    const SceneFluidSurfaceDefinition& surface,
    const fluid::SceneFluidGridPatchSet& patches,
    const std::size_t sourcePatchIndex,
    const SceneFluidQuadratureOwnerKind ownerKind,
    const std::uint64_t stableId,
    const StableId negativeRegion,
    const StableId positiveRegion,
    const std::size_t negativeCellIndex,
    const std::size_t positiveCellIndex) {
    const auto& patch = patches.patches[sourcePatchIndex];
    const auto& triangle = surface.triangles[patch.triangleIndex];
    return {
        stableId,
        patch.triangleId,
        patch.centroidBarycentricCoordinates,
        patch.areaSquareMeters,
        negativeRegion,
        positiveRegion,
        negativeCellIndex,
        positiveCellIndex,
        surface.materials[triangle.materialIndex].id,
        triangle.sheetId,
        triangle.role,
        ownerKind,
    };
}

std::size_t lowerCellIndex(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::GridFaceAxis axis,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) {
    const auto counts = grid.cellCounts();
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return grid.cellIndex(
            (i + counts.x - 1) % counts.x, j, k);
    case fluid::GridFaceAxis::Y:
        return grid.cellIndex(
            i, (j + counts.y - 1) % counts.y, k);
    case fluid::GridFaceAxis::Z:
        return grid.cellIndex(
            i, j, (k + counts.z - 1) % counts.z);
    }
    throw std::invalid_argument(
        "scene fluid quadrature has an invalid face axis");
}

} // namespace

void validateSceneFluidQuadratureDefinition(
    const SceneFluidQuadratureDefinition& definition) {
    if (definition.version != sceneFluidQuadratureVersion
        || definition.fingerprint == 0
        || definition.surfaceDefinitionFingerprint == 0
        || definition.surfaceStateFingerprint == 0
        || definition.patchOwnershipFingerprint == 0
        || definition.couplingSurfaceFingerprint == 0
        || definition.structureDefinitionFingerprint == 0
        || !std::isfinite(definition.simulationTimeSeconds)
        || definition.simulationTimeSeconds < 0.0
        || definition.points.empty()) {
        throw std::invalid_argument(
            "scene fluid quadrature identity is invalid");
    }
    std::set<std::uint64_t> stableIds;
    std::pair<StableId, std::uint64_t> previousKey{};
    bool havePrevious = false;
    for (const auto& point : definition.points) {
        const std::pair key{point.triangleId, point.stableId};
        double barycentricSum = 0.0;
        for (const double coordinate : point.barycentricCoordinates) {
            if (!std::isfinite(coordinate)
                || coordinate < -1.0e-12
                || coordinate > 1.0 + 1.0e-12) {
                throw std::invalid_argument(
                    "scene fluid quadrature barycentric coordinate is invalid");
            }
            barycentricSum += coordinate;
        }
        if (point.stableId == 0 || point.triangleId == invalidStableId
            || !stableIds.insert(point.stableId).second
            || (havePrevious && !(previousKey < key))
            || !std::isfinite(barycentricSum)
            || std::abs(barycentricSum - 1.0) > 1.0e-12
            || !std::isfinite(point.areaSquareMeters)
            || !(point.areaSquareMeters > 0.0)
            || point.negativeSideRegionId == invalidStableId
            || point.positiveSideRegionId == invalidStableId
            || point.materialId == invalidStableId
            || point.sheetId == invalidStableId
            || (point.ownerKind != SceneFluidQuadratureOwnerKind::Cell
                && point.ownerKind != SceneFluidQuadratureOwnerKind::Face)) {
            throw std::invalid_argument(
                "scene fluid quadrature point is invalid or out of order");
        }
        havePrevious = true;
        previousKey = key;
    }
    if (definition.fingerprint != definitionFingerprint(definition)) {
        throw std::invalid_argument(
            "scene fluid quadrature fingerprint does not match its payload");
    }
}

SceneFluidQuadratureDefinition buildSceneFluidQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::SceneFluidGridCandidateSet& candidates,
    const fluid::SceneFluidGridIntersectionSet& intersections,
    const fluid::SceneFluidGridPatchSet& patches,
    const fluid::SceneFluidPatchOwnership& ownership,
    const SceneFluidSurfaceTransfer& transfer) {
    fluid::validateSceneFluidPatchOwnership(
        ownership, surface, state, grid, candidates, intersections, patches);
    if (transfer.surfaceDefinitionFingerprint() != surface.fingerprint
        || transfer.targetDefinitionFingerprint()
            != state.structureDefinitionFingerprint) {
        throw std::invalid_argument(
            "scene fluid quadrature transfer binding is foreign");
    }

    SceneFluidQuadratureDefinition result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.patchOwnershipFingerprint = ownership.fingerprint;
    result.couplingSurfaceFingerprint =
        transfer.couplingSurfaceFingerprint();
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.points.reserve(
        ownership.cellPatches.size() + ownership.facePatches.size());
    for (const auto& owned : ownership.cellPatches) {
        const std::uint64_t stableId = pointStableId(
            SceneFluidQuadratureOwnerKind::Cell,
            owned.triangleId,
            {static_cast<std::uint64_t>(owned.cellIndex), 0, 0, 0});
        result.points.push_back(pointFromPatch(
            surface, patches, owned.sourcePatchIndex,
            SceneFluidQuadratureOwnerKind::Cell, stableId,
            owned.negativeSideRegionId, owned.positiveSideRegionId,
            owned.cellIndex, owned.cellIndex));
    }
    for (const auto& owned : ownership.facePatches) {
        const std::uint64_t stableId = pointStableId(
            SceneFluidQuadratureOwnerKind::Face,
            owned.triangleId,
            {static_cast<std::uint64_t>(owned.axis),
             static_cast<std::uint64_t>(owned.i),
             static_cast<std::uint64_t>(owned.j),
             static_cast<std::uint64_t>(owned.k)});
        const std::size_t lower = lowerCellIndex(
            grid, owned.axis, owned.i, owned.j, owned.k);
        const std::size_t upper = grid.cellIndex(
            owned.i, owned.j, owned.k);
        const std::size_t negativeCell =
            owned.triangleNormalAxisSign > 0 ? lower : upper;
        const std::size_t positiveCell =
            owned.triangleNormalAxisSign > 0 ? upper : lower;
        result.points.push_back(pointFromPatch(
            surface, patches, owned.lowerCellSourcePatchIndex,
            SceneFluidQuadratureOwnerKind::Face, stableId,
            owned.negativeSideRegionId, owned.positiveSideRegionId,
            negativeCell, positiveCell));
    }
    std::sort(result.points.begin(), result.points.end(),
              [](const auto& first, const auto& second) {
                  return std::pair{first.triangleId, first.stableId}
                      < std::pair{second.triangleId, second.stableId};
              });
    result.fingerprint = definitionFingerprint(result);
    validateSceneFluidQuadratureDefinition(result);
    return result;
}

ConservativeTransferResult evaluateSceneFluidQuadrature(
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& definition,
    const std::span<const SceneFluidQuadratureTraction> tractions,
    const ConservativeTransferSettings& settings) {
    validateSceneFluidSurfaceState(state);
    validateSceneFluidQuadratureDefinition(definition);
    if (definition.surfaceDefinitionFingerprint
            != transfer.surfaceDefinitionFingerprint()
        || definition.surfaceStateFingerprint != state.fingerprint
        || definition.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || definition.structureDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || definition.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || definition.acceptedStepCount != state.acceptedStepCount
        || definition.simulationTimeSeconds != state.simulationTimeSeconds
        || tractions.size() != definition.points.size()) {
        throw std::invalid_argument(
            "scene fluid quadrature evaluation binding is invalid");
    }

    std::vector<CouplingTriangleTractionQuadrature> coupling;
    coupling.reserve(definition.points.size());
    for (std::size_t index = 0; index < definition.points.size(); ++index) {
        const auto& point = definition.points[index];
        const auto& traction = tractions[index];
        if (traction.stableId != point.stableId
            || !finite(traction.tractionPascals)) {
            throw std::invalid_argument(
                "scene fluid quadrature traction is non-finite or out of order");
        }
        coupling.push_back({
            point.stableId,
            point.triangleId,
            point.barycentricCoordinates,
            point.areaSquareMeters,
            traction.tractionPascals,
        });
    }
    return transfer.evaluateQuadrature(state, coupling, settings);
}

std::vector<SceneFluidQuadratureKinematics>
sampleSceneFluidQuadratureKinematics(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& definition) {
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidQuadratureDefinition(definition);
    if (definition.surfaceDefinitionFingerprint != surface.fingerprint
        || definition.surfaceStateFingerprint != state.fingerprint
        || definition.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || definition.acceptedStepCount != state.acceptedStepCount
        || definition.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid quadrature kinematics binding is invalid");
    }

    std::vector<SceneFluidQuadratureKinematics> result;
    result.reserve(definition.points.size());
    for (const auto& point : definition.points) {
        const auto triangleIndex = surface.mappings.triangleIndex(
            point.triangleId);
        if (!triangleIndex) {
            throw std::invalid_argument(
                "scene fluid quadrature kinematics references a foreign triangle");
        }
        const auto& triangle = surface.triangles[*triangleIndex];
        SceneFluidQuadratureKinematics sample;
        sample.stableId = point.stableId;
        for (std::size_t corner = 0; corner < 3; ++corner) {
            const double weight = point.barycentricCoordinates[corner];
            const auto& vertex = state.vertices[
                triangle.vertexIndices[corner]];
            sample.positionMeters.x += weight * vertex.positionMeters.x;
            sample.positionMeters.y += weight * vertex.positionMeters.y;
            sample.positionMeters.z += weight * vertex.positionMeters.z;
            sample.velocityMetersPerSecond.x +=
                weight * vertex.velocityMetersPerSecond.x;
            sample.velocityMetersPerSecond.y +=
                weight * vertex.velocityMetersPerSecond.y;
            sample.velocityMetersPerSecond.z +=
                weight * vertex.velocityMetersPerSecond.z;
        }
        if (!finite(sample.positionMeters)
            || !finite(sample.velocityMetersPerSecond)) {
            throw std::overflow_error(
                "scene fluid quadrature kinematics is not finite");
        }
        result.push_back(sample);
    }
    return result;
}

} // namespace simwing::fsi
