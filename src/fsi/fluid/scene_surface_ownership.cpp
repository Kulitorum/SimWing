#include "fluid/scene_surface_ownership.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::size_t missingPatch = std::numeric_limits<std::size_t>::max();

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

bool checkedMultiply(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

struct FaceLocation {
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    bool sourceIsLowerCell = false;
};

FaceLocation faceLocation(const SceneFluidCellPatch& patch,
                          const GridCellCounts counts) {
    const std::uint8_t mask = patch.coincidentBoundaryPlanes;
    if (mask == 0 || (mask & (mask - 1U)) != 0) {
        throw std::invalid_argument(
            "positive-area scene patch must be coincident with at most one cell plane");
    }
    if (mask == CellBoundaryXMinus) {
        if (patch.cell.i == 0) {
            throw std::invalid_argument(
                "scene patch on periodic X boundary has no explicit image owner");
        }
        return {GridFaceAxis::X,
                patch.cell.i, patch.cell.j, patch.cell.k, false};
    }
    if (mask == CellBoundaryXPlus) {
        if (patch.cell.i + 1 >= counts.x) {
            throw std::invalid_argument(
                "scene patch on periodic X boundary has no explicit image owner");
        }
        return {GridFaceAxis::X,
                patch.cell.i + 1, patch.cell.j, patch.cell.k, true};
    }
    if (mask == CellBoundaryYMinus) {
        if (patch.cell.j == 0) {
            throw std::invalid_argument(
                "scene patch on periodic Y boundary has no explicit image owner");
        }
        return {GridFaceAxis::Y,
                patch.cell.i, patch.cell.j, patch.cell.k, false};
    }
    if (mask == CellBoundaryYPlus) {
        if (patch.cell.j + 1 >= counts.y) {
            throw std::invalid_argument(
                "scene patch on periodic Y boundary has no explicit image owner");
        }
        return {GridFaceAxis::Y,
                patch.cell.i, patch.cell.j + 1, patch.cell.k, true};
    }
    if (mask == CellBoundaryZMinus) {
        if (patch.cell.k == 0) {
            throw std::invalid_argument(
                "scene patch on periodic Z boundary has no explicit image owner");
        }
        return {GridFaceAxis::Z,
                patch.cell.i, patch.cell.j, patch.cell.k, false};
    }
    if (mask == CellBoundaryZPlus) {
        if (patch.cell.k + 1 >= counts.z) {
            throw std::invalid_argument(
                "scene patch on periodic Z boundary has no explicit image owner");
        }
        return {GridFaceAxis::Z,
                patch.cell.i, patch.cell.j, patch.cell.k + 1, true};
    }
    throw std::invalid_argument(
        "scene patch has an unknown coincident cell plane");
}

double axisNormalComponent(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const std::size_t triangleIndex,
    const GridFaceAxis axis) {
    const auto& triangle = surface.triangles[triangleIndex];
    const Vec3& first =
        state.vertices[triangle.vertexIndices[0]].positionMeters;
    const Vec3& second =
        state.vertices[triangle.vertexIndices[1]].positionMeters;
    const Vec3& third =
        state.vertices[triangle.vertexIndices[2]].positionMeters;
    const Vec3 firstEdge{
        second.x - first.x,
        second.y - first.y,
        second.z - first.z,
    };
    const Vec3 secondEdge{
        third.x - first.x,
        third.y - first.y,
        third.z - first.z,
    };
    const Vec3 normal{
        firstEdge.y * secondEdge.z - firstEdge.z * secondEdge.y,
        firstEdge.z * secondEdge.x - firstEdge.x * secondEdge.z,
        firstEdge.x * secondEdge.y - firstEdge.y * secondEdge.x,
    };
    if (axis == GridFaceAxis::X) {
        return normal.x;
    }
    if (axis == GridFaceAxis::Y) {
        return normal.y;
    }
    return normal.z;
}

using FaceKey = std::tuple<std::uint8_t,
                           std::size_t,
                           std::size_t,
                           std::size_t,
                           std::size_t>;

struct FacePair {
    FaceLocation location;
    std::size_t lowerPatch = missingPatch;
    std::size_t upperPatch = missingPatch;
};

bool samePatchGeometry(const SceneFluidGridPatchSet& patches,
                       const SceneFluidCellPatch& first,
                       const SceneFluidCellPatch& second) {
    if (first.triangleIndex != second.triangleIndex
        || first.triangleId != second.triangleId
        || first.dimension != second.dimension
        || first.areaSquareMeters != second.areaSquareMeters
        || first.centroidMeters.x != second.centroidMeters.x
        || first.centroidMeters.y != second.centroidMeters.y
        || first.centroidMeters.z != second.centroidMeters.z
        || first.centroidBarycentricCoordinates
            != second.centroidBarycentricCoordinates) {
        return false;
    }
    const auto firstVertices = patches.verticesForPatch(first);
    const auto secondVertices = patches.verticesForPatch(second);
    return std::equal(firstVertices.begin(), firstVertices.end(),
                      secondVertices.begin(), secondVertices.end());
}

std::uint64_t ownershipFingerprint(
    const SceneFluidPatchOwnership& ownership) {
    Fingerprint fingerprint;
    fingerprint.integer(ownership.version);
    fingerprint.integer(ownership.surfaceDefinitionFingerprint);
    fingerprint.integer(ownership.surfaceStateFingerprint);
    fingerprint.integer(ownership.patchSetFingerprint);
    fingerprint.integer(ownership.structureDefinitionFingerprint);
    fingerprint.integer(ownership.acceptedStepCount);
    fingerprint.real(ownership.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        ownership.pointContactPatchCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        ownership.segmentContactPatchCount));
    fingerprint.real(ownership.ownedAreaSquareMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        ownership.cellPatches.size()));
    for (const auto& patch : ownership.cellPatches) {
        fingerprint.integer(static_cast<std::uint64_t>(patch.sourcePatchIndex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(patch.triangleIndex));
        fingerprint.integer(patch.triangleId);
        fingerprint.integer(patch.negativeSideRegionId);
        fingerprint.integer(patch.positiveSideRegionId);
        fingerprint.real(patch.areaSquareMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        ownership.facePatches.size()));
    for (const auto& patch : ownership.facePatches) {
        fingerprint.integer(static_cast<std::uint64_t>(
            patch.lowerCellSourcePatchIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            patch.upperCellSourcePatchIndex));
        fingerprint.enumeration(patch.axis);
        fingerprint.integer(static_cast<std::uint64_t>(patch.i));
        fingerprint.integer(static_cast<std::uint64_t>(patch.j));
        fingerprint.integer(static_cast<std::uint64_t>(patch.k));
        fingerprint.integer(static_cast<std::uint64_t>(patch.triangleIndex));
        fingerprint.integer(patch.triangleId);
        fingerprint.integer(patch.negativeSideRegionId);
        fingerprint.integer(patch.positiveSideRegionId);
        fingerprint.integer(static_cast<std::uint8_t>(
            patch.triangleNormalAxisSign));
        fingerprint.real(patch.areaSquareMeters);
    }
    return fingerprint.value();
}

SceneFluidPatchOwnership buildOwnership(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnershipLimits& limits) {
    SceneFluidPatchOwnership result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.patchSetFingerprint = patches.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;

    std::map<FaceKey, FacePair> facePairs;
    for (std::size_t patchIndex = 0;
         patchIndex < patches.patches.size(); ++patchIndex) {
        const auto& patch = patches.patches[patchIndex];
        if (patch.dimension == SceneFluidPatchDimension::Point) {
            ++result.pointContactPatchCount;
            continue;
        }
        if (patch.dimension == SceneFluidPatchDimension::Segment) {
            ++result.segmentContactPatchCount;
            continue;
        }
        const auto& triangle = surface.triangles[patch.triangleIndex];
        const StableId negativeRegion =
            surface.regions[triangle.negativeSideRegionIndex].id;
        const StableId positiveRegion =
            surface.regions[triangle.positiveSideRegionIndex].id;
        if (patch.coincidentBoundaryPlanes == CellBoundaryNone) {
            if (result.cellPatches.size() + result.facePatches.size()
                == limits.maximumOwnedAreaPatches) {
                throw std::length_error(
                    "scene fluid owned patches exceed their count limit");
            }
            result.cellPatches.push_back({
                patchIndex,
                patch.cellIndex,
                patch.cell,
                patch.triangleIndex,
                patch.triangleId,
                negativeRegion,
                positiveRegion,
                patch.areaSquareMeters,
            });
            result.ownedAreaSquareMeters += patch.areaSquareMeters;
            continue;
        }
        const FaceLocation location = faceLocation(
            patch, grid.cellCounts());
        const FaceKey key{
            static_cast<std::uint8_t>(location.axis),
            location.i, location.j, location.k,
            patch.triangleIndex,
        };
        FacePair& pair = facePairs[key];
        pair.location = location;
        std::size_t& destination = location.sourceIsLowerCell
            ? pair.lowerPatch : pair.upperPatch;
        if (destination != missingPatch) {
            throw std::invalid_argument(
                "scene fluid face ownership has a duplicate cell-side patch");
        }
        destination = patchIndex;
    }

    for (const auto& [key, pair] : facePairs) {
        static_cast<void>(key);
        if (pair.lowerPatch == missingPatch
            || pair.upperPatch == missingPatch) {
            throw std::invalid_argument(
                "scene fluid face ownership requires both adjacent cell patches");
        }
        if (result.cellPatches.size() + result.facePatches.size()
            == limits.maximumOwnedAreaPatches) {
            throw std::length_error(
                "scene fluid owned patches exceed their count limit");
        }
        const auto& lower = patches.patches[pair.lowerPatch];
        const auto& upper = patches.patches[pair.upperPatch];
        if (!samePatchGeometry(patches, lower, upper)) {
            throw std::invalid_argument(
                "scene fluid adjacent face patches do not share exact geometry");
        }
        const auto& triangle = surface.triangles[lower.triangleIndex];
        const double normalComponent = axisNormalComponent(
            surface, state, lower.triangleIndex, pair.location.axis);
        if (!std::isfinite(normalComponent) || normalComponent == 0.0) {
            throw std::invalid_argument(
                "scene fluid face patch winding is not axis-resolvable");
        }
        result.facePatches.push_back({
            pair.lowerPatch,
            pair.upperPatch,
            pair.location.axis,
            pair.location.i,
            pair.location.j,
            pair.location.k,
            lower.triangleIndex,
            lower.triangleId,
            surface.regions[triangle.negativeSideRegionIndex].id,
            surface.regions[triangle.positiveSideRegionIndex].id,
            static_cast<std::int8_t>(normalComponent > 0.0 ? 1 : -1),
            lower.areaSquareMeters,
        });
        result.ownedAreaSquareMeters += lower.areaSquareMeters;
    }

    std::size_t cellBytes = 0;
    std::size_t faceBytes = 0;
    if (!checkedMultiply(result.cellPatches.size(),
                         sizeof(SceneFluidOwnedCellPatch), cellBytes)
        || !checkedMultiply(result.facePatches.size(),
                            sizeof(SceneFluidOwnedFacePatch), faceBytes)
        || faceBytes > std::numeric_limits<std::size_t>::max() - cellBytes
        || cellBytes + faceBytes > limits.maximumOwnershipBytes
        || !std::isfinite(result.ownedAreaSquareMeters)) {
        throw std::length_error(
            "scene fluid patch ownership exceeds its storage or area limit");
    }
    result.fingerprint = ownershipFingerprint(result);
    return result;
}

} // namespace

SceneFluidPatchOwnership ownSceneFluidSurfacePatches(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnershipLimits& limits) {
    validateSceneFluidGridPatches(
        patches, surface, state, grid, candidates, intersections);
    SceneFluidPatchOwnership result = buildOwnership(
        surface, state, grid, patches, limits);
    validateSceneFluidPatchOwnership(
        result, surface, state, grid, candidates, intersections, patches);
    return result;
}

void validateSceneFluidPatchOwnership(
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches) {
    validateSceneFluidGridPatches(
        patches, surface, state, grid, candidates, intersections);
    if (ownership.version != sceneFluidPatchOwnershipVersion
        || ownership.fingerprint == 0
        || ownership.surfaceDefinitionFingerprint != surface.fingerprint
        || ownership.surfaceStateFingerprint != state.fingerprint
        || ownership.patchSetFingerprint != patches.fingerprint
        || ownership.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || ownership.acceptedStepCount != state.acceptedStepCount
        || ownership.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid patch-ownership identity is invalid");
    }
    const SceneFluidPatchOwnership expected = buildOwnership(
        surface,
        state,
        grid,
        patches,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (ownership != expected) {
        throw std::invalid_argument(
            "scene fluid patch ownership does not match its source patches");
    }
}

} // namespace simwing::fsi::fluid
