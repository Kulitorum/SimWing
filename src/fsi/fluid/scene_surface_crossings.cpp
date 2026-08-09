#include "fluid/scene_surface_crossings.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace simwing::fsi::fluid {
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

double coordinate(const Vec3& value, const std::size_t axis) {
    if (axis == 0) {
        return value.x;
    }
    if (axis == 1) {
        return value.y;
    }
    return value.z;
}

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

double squaredNorm(const Vec3& value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}

double norm(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

bool vertexLess(const SceneFluidClippedVertex& first,
                const SceneFluidClippedVertex& second) {
    return std::tie(first.positionMeters.x,
                    first.positionMeters.y,
                    first.positionMeters.z,
                    first.barycentricCoordinates)
        < std::tie(second.positionMeters.x,
                   second.positionMeters.y,
                   second.positionMeters.z,
                   second.barycentricCoordinates);
}

struct FaceLocation {
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    bool sourceIsLowerCell = false;
};

struct BoundaryPlane {
    FaceLocation location;
    std::size_t coordinateAxis = 0;
    double coordinateMeters = 0.0;
};

std::array<BoundaryPlane, 6> boundaryPlanes(
    const SceneFluidCellPatch& patch,
    const PeriodicCartesianGrid& grid) {
    const Vector3 gridLower = grid.lowerMeters();
    const Vector3 spacing = grid.cellSpacingMeters();
    const Vec3 lower{
        gridLower.x + static_cast<double>(patch.cell.i) * spacing.x,
        gridLower.y + static_cast<double>(patch.cell.j) * spacing.y,
        gridLower.z + static_cast<double>(patch.cell.k) * spacing.z,
    };
    const Vec3 upper{
        gridLower.x
            + static_cast<double>(patch.cell.i + 1) * spacing.x,
        gridLower.y
            + static_cast<double>(patch.cell.j + 1) * spacing.y,
        gridLower.z
            + static_cast<double>(patch.cell.k + 1) * spacing.z,
    };
    return {{
        {{GridFaceAxis::X, patch.cell.i, patch.cell.j, patch.cell.k, false},
         0, lower.x},
        {{GridFaceAxis::X, patch.cell.i + 1, patch.cell.j, patch.cell.k, true},
         0, upper.x},
        {{GridFaceAxis::Y, patch.cell.i, patch.cell.j, patch.cell.k, false},
         1, lower.y},
        {{GridFaceAxis::Y, patch.cell.i, patch.cell.j + 1, patch.cell.k, true},
         1, upper.y},
        {{GridFaceAxis::Z, patch.cell.i, patch.cell.j, patch.cell.k, false},
         2, lower.z},
        {{GridFaceAxis::Z, patch.cell.i, patch.cell.j, patch.cell.k + 1, true},
         2, upper.z},
    }};
}

struct CandidateSegment {
    FaceLocation location;
    std::size_t sourcePatchIndex = 0;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    SceneFluidClippedVertex first;
    SceneFluidClippedVertex second;
    bool gridEdgeAligned = false;
};

std::optional<CandidateSegment> candidateSegment(
    const SceneFluidCellPatch& patch,
    const std::span<const SceneFluidClippedVertex> vertices,
    const BoundaryPlane& plane) {
    std::vector<SceneFluidClippedVertex> onPlane;
    onPlane.reserve(vertices.size());
    for (const auto& vertex : vertices) {
        if (coordinate(vertex.positionMeters, plane.coordinateAxis)
            == plane.coordinateMeters) {
            onPlane.push_back(vertex);
        }
    }
    if (onPlane.size() < 2) {
        return std::nullopt;
    }

    double greatestDistanceSquared = 0.0;
    std::size_t firstIndex = 0;
    std::size_t secondIndex = 0;
    for (std::size_t first = 0; first < onPlane.size(); ++first) {
        for (std::size_t second = first + 1;
             second < onPlane.size(); ++second) {
            const double distanceSquared = squaredNorm(subtract(
                onPlane[second].positionMeters,
                onPlane[first].positionMeters));
            if (distanceSquared > greatestDistanceSquared) {
                greatestDistanceSquared = distanceSquared;
                firstIndex = first;
                secondIndex = second;
            }
        }
    }
    if (!(greatestDistanceSquared > 0.0)
        || !std::isfinite(greatestDistanceSquared)) {
        return std::nullopt;
    }

    CandidateSegment result;
    result.location = plane.location;
    result.triangleIndex = patch.triangleIndex;
    result.triangleId = patch.triangleId;
    result.first = onPlane[firstIndex];
    result.second = onPlane[secondIndex];
    if (vertexLess(result.second, result.first)) {
        std::swap(result.first, result.second);
    }

    return result;
}

std::size_t commonBoundaryPlaneCount(
    const SceneFluidCellPatch& patch,
    const PeriodicCartesianGrid& grid,
    const SceneFluidClippedVertex& first,
    const SceneFluidClippedVertex& second) {
    std::size_t count = 0;
    for (const auto& plane : boundaryPlanes(patch, grid)) {
        if (coordinate(first.positionMeters, plane.coordinateAxis)
                == plane.coordinateMeters
            && coordinate(second.positionMeters, plane.coordinateAxis)
                == plane.coordinateMeters) {
            ++count;
        }
    }
    return count;
}

bool near(const double first,
          const double second,
          const double tolerance) {
    return std::isfinite(first) && std::isfinite(second)
        && std::abs(first - second) <= tolerance;
}

bool sameVertexWithinRoundoff(
    const SceneFluidClippedVertex& first,
    const SceneFluidClippedVertex& second,
    const double positionTolerance) {
    constexpr double barycentricTolerance =
        256.0 * std::numeric_limits<double>::epsilon();
    return near(first.positionMeters.x, second.positionMeters.x,
                positionTolerance)
        && near(first.positionMeters.y, second.positionMeters.y,
                positionTolerance)
        && near(first.positionMeters.z, second.positionMeters.z,
                positionTolerance)
        && near(first.barycentricCoordinates[0],
                second.barycentricCoordinates[0], barycentricTolerance)
        && near(first.barycentricCoordinates[1],
                second.barycentricCoordinates[1], barycentricTolerance)
        && near(first.barycentricCoordinates[2],
                second.barycentricCoordinates[2], barycentricTolerance);
}

bool sameSegmentWithinRoundoff(const CandidateSegment& first,
                               const CandidateSegment& second) {
    if (first.triangleIndex != second.triangleIndex
        || first.triangleId != second.triangleId) {
        return false;
    }
    const double scale = std::max({
        1.0,
        std::abs(first.first.positionMeters.x),
        std::abs(first.first.positionMeters.y),
        std::abs(first.first.positionMeters.z),
        std::abs(first.second.positionMeters.x),
        std::abs(first.second.positionMeters.y),
        std::abs(first.second.positionMeters.z),
        std::abs(second.first.positionMeters.x),
        std::abs(second.first.positionMeters.y),
        std::abs(second.first.positionMeters.z),
        std::abs(second.second.positionMeters.x),
        std::abs(second.second.positionMeters.y),
        std::abs(second.second.positionMeters.z),
    });
    const double positionTolerance =
        256.0 * std::numeric_limits<double>::epsilon() * scale;
    const bool direct = sameVertexWithinRoundoff(
                            first.first, second.first, positionTolerance)
        && sameVertexWithinRoundoff(
            first.second, second.second, positionTolerance);
    const bool reversed = sameVertexWithinRoundoff(
                              first.first, second.second, positionTolerance)
        && sameVertexWithinRoundoff(
            first.second, second.first, positionTolerance);
    return direct || reversed;
}

using FaceKey = std::tuple<std::uint8_t,
                           std::size_t,
                           std::size_t,
                           std::size_t,
                           std::size_t>;

struct FacePair {
    FaceLocation location;
    std::optional<CandidateSegment> lower;
    std::optional<CandidateSegment> upper;
};

std::uint64_t crossingStableId(const FaceLocation& location,
                               const StableId triangleId) {
    Fingerprint fingerprint;
    fingerprint.integer(sceneFluidFaceCrossingVersion);
    fingerprint.enumeration(location.axis);
    fingerprint.integer(static_cast<std::uint64_t>(location.i));
    fingerprint.integer(static_cast<std::uint64_t>(location.j));
    fingerprint.integer(static_cast<std::uint64_t>(location.k));
    fingerprint.integer(triangleId);
    return fingerprint.value();
}

Vec3 inFacePositiveDirection(
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
    Vec3 normal = cross(subtract(second, first), subtract(third, first));
    if (axis == GridFaceAxis::X) {
        normal.x = 0.0;
    } else if (axis == GridFaceAxis::Y) {
        normal.y = 0.0;
    } else {
        normal.z = 0.0;
    }
    const double magnitude = norm(normal);
    if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
        throw std::invalid_argument(
            "scene fluid transverse crossing has no in-face orientation");
    }
    return {normal.x / magnitude, normal.y / magnitude, normal.z / magnitude};
}

std::uint64_t crossingsFingerprint(
    const SceneFluidFaceCrossingSet& crossings) {
    Fingerprint fingerprint;
    fingerprint.integer(crossings.version);
    fingerprint.integer(crossings.surfaceDefinitionFingerprint);
    fingerprint.integer(crossings.surfaceStateFingerprint);
    fingerprint.integer(crossings.patchOwnershipFingerprint);
    fingerprint.integer(crossings.structureDefinitionFingerprint);
    fingerprint.integer(crossings.acceptedStepCount);
    fingerprint.real(crossings.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.candidateSegmentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.unpairedContactSegmentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.coplanarAreaPatchCount));
    fingerprint.real(crossings.crossingLengthMeters);
    fingerprint.integer(static_cast<std::uint64_t>(crossings.crossings.size()));
    for (const auto& crossing : crossings.crossings) {
        fingerprint.integer(crossing.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.lowerCellSourcePatchIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.upperCellSourcePatchIndex));
        fingerprint.enumeration(crossing.axis);
        fingerprint.integer(static_cast<std::uint64_t>(crossing.i));
        fingerprint.integer(static_cast<std::uint64_t>(crossing.j));
        fingerprint.integer(static_cast<std::uint64_t>(crossing.k));
        fingerprint.integer(static_cast<std::uint64_t>(crossing.triangleIndex));
        fingerprint.integer(crossing.triangleId);
        for (const auto* vertex : {&crossing.first, &crossing.second}) {
            fingerprint.real(vertex->positionMeters.x);
            fingerprint.real(vertex->positionMeters.y);
            fingerprint.real(vertex->positionMeters.z);
            for (const double value : vertex->barycentricCoordinates) {
                fingerprint.real(value);
            }
        }
        fingerprint.real(crossing.midpointMeters.x);
        fingerprint.real(crossing.midpointMeters.y);
        fingerprint.real(crossing.midpointMeters.z);
        fingerprint.real(crossing.negativeToPositiveDirectionInFace.x);
        fingerprint.real(crossing.negativeToPositiveDirectionInFace.y);
        fingerprint.real(crossing.negativeToPositiveDirectionInFace.z);
        fingerprint.real(crossing.lengthMeters);
        fingerprint.integer(crossing.negativeSideRegionId);
        fingerprint.integer(crossing.positiveSideRegionId);
        fingerprint.integer(crossing.materialId);
        fingerprint.integer(crossing.sheetId);
        fingerprint.enumeration(crossing.role);
    }
    return fingerprint.value();
}

SceneFluidFaceCrossingSet buildCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingLimits& limits) {
    SceneFluidFaceCrossingSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.patchOwnershipFingerprint = ownership.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.coplanarAreaPatchCount = ownership.facePatches.size();

    std::map<FaceKey, FacePair> pairs;
    for (const auto& owned : ownership.cellPatches) {
        const auto& patch = patches.patches[owned.sourcePatchIndex];
        const auto vertices = patches.verticesForPatch(patch);
        for (const auto& plane : boundaryPlanes(patch, grid)) {
            auto candidate = candidateSegment(patch, vertices, plane);
            if (!candidate) {
                continue;
            }
            if (result.candidateSegmentCount
                == limits.maximumCandidateSegments) {
                throw std::length_error(
                    "scene fluid face-crossing candidates exceed their limit");
            }
            ++result.candidateSegmentCount;
            candidate->sourcePatchIndex = owned.sourcePatchIndex;
            candidate->gridEdgeAligned = commonBoundaryPlaneCount(
                patch, grid, candidate->first, candidate->second) > 1;
            const FaceKey key{
                static_cast<std::uint8_t>(plane.location.axis),
                plane.location.i,
                plane.location.j,
                plane.location.k,
                patch.triangleIndex,
            };
            FacePair& pair = pairs[key];
            pair.location = plane.location;
            auto& side = plane.location.sourceIsLowerCell
                ? pair.lower : pair.upper;
            if (side) {
                throw std::invalid_argument(
                    "scene fluid face crossing has a duplicate cell-side segment");
            }
            side = std::move(candidate);
        }
    }

    std::set<std::uint64_t> stableIds;
    for (const auto& [key, pair] : pairs) {
        static_cast<void>(key);
        if (!pair.lower || !pair.upper) {
            ++result.unpairedContactSegmentCount;
            continue;
        }
        if (pair.lower->gridEdgeAligned || pair.upper->gridEdgeAligned) {
            throw std::invalid_argument(
                "scene fluid crossing on a grid edge needs explicit edge ownership");
        }
        if (!sameSegmentWithinRoundoff(*pair.lower, *pair.upper)) {
            throw std::invalid_argument(
                "scene fluid adjacent face crossings disagree beyond roundoff");
        }
        if (result.crossings.size() == limits.maximumCrossings) {
            throw std::length_error(
                "scene fluid face crossings exceed their count limit");
        }
        const auto& triangle =
            surface.triangles[pair.lower->triangleIndex];
        const Vec3 delta = subtract(pair.lower->second.positionMeters,
                                    pair.lower->first.positionMeters);
        const double length = norm(delta);
        if (!std::isfinite(length) || !(length > 0.0)) {
            throw std::invalid_argument(
                "scene fluid face crossing length is invalid");
        }
        const std::uint64_t stableId = crossingStableId(
            pair.location, pair.lower->triangleId);
        if (!stableIds.insert(stableId).second) {
            throw std::invalid_argument(
                "scene fluid face crossing stable-ID collision");
        }
        SceneFluidFaceCrossing crossing;
        crossing.stableId = stableId;
        crossing.lowerCellSourcePatchIndex = pair.lower->sourcePatchIndex;
        crossing.upperCellSourcePatchIndex = pair.upper->sourcePatchIndex;
        crossing.axis = pair.location.axis;
        crossing.i = pair.location.i;
        crossing.j = pair.location.j;
        crossing.k = pair.location.k;
        crossing.triangleIndex = pair.lower->triangleIndex;
        crossing.triangleId = pair.lower->triangleId;
        crossing.first = pair.lower->first;
        crossing.second = pair.lower->second;
        crossing.midpointMeters = {
            0.5 * (crossing.first.positionMeters.x
                   + crossing.second.positionMeters.x),
            0.5 * (crossing.first.positionMeters.y
                   + crossing.second.positionMeters.y),
            0.5 * (crossing.first.positionMeters.z
                   + crossing.second.positionMeters.z),
        };
        crossing.negativeToPositiveDirectionInFace = inFacePositiveDirection(
            surface, state, crossing.triangleIndex, crossing.axis);
        crossing.lengthMeters = length;
        crossing.negativeSideRegionId =
            surface.regions[triangle.negativeSideRegionIndex].id;
        crossing.positiveSideRegionId =
            surface.regions[triangle.positiveSideRegionIndex].id;
        crossing.materialId = surface.materials[triangle.materialIndex].id;
        crossing.sheetId = triangle.sheetId;
        crossing.role = triangle.role;
        result.crossingLengthMeters += length;
        result.crossings.push_back(std::move(crossing));
    }

    std::size_t crossingBytes = 0;
    if (!checkedMultiply(result.crossings.size(),
                         sizeof(SceneFluidFaceCrossing), crossingBytes)
        || crossingBytes > limits.maximumCrossingBytes
        || !std::isfinite(result.crossingLengthMeters)) {
        throw std::length_error(
            "scene fluid face crossings exceed their storage or length limit");
    }
    result.fingerprint = crossingsFingerprint(result);
    return result;
}

} // namespace

SceneFluidFaceCrossingSet buildSceneFluidFaceCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingLimits& limits) {
    validateSceneFluidPatchOwnership(
        ownership, surface, state, grid, candidates, intersections, patches);
    SceneFluidFaceCrossingSet result = buildCrossings(
        surface, state, grid, patches, ownership, limits);
    validateSceneFluidFaceCrossings(
        result, surface, state, grid, candidates, intersections, patches,
        ownership);
    return result;
}

void validateSceneFluidFaceCrossings(
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership) {
    validateSceneFluidPatchOwnership(
        ownership, surface, state, grid, candidates, intersections, patches);
    if (crossings.version != sceneFluidFaceCrossingVersion
        || crossings.fingerprint == 0
        || crossings.surfaceDefinitionFingerprint != surface.fingerprint
        || crossings.surfaceStateFingerprint != state.fingerprint
        || crossings.patchOwnershipFingerprint != ownership.fingerprint
        || crossings.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || crossings.acceptedStepCount != state.acceptedStepCount
        || crossings.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid face-crossing identity is invalid");
    }
    const SceneFluidFaceCrossingSet expected = buildCrossings(
        surface,
        state,
        grid,
        patches,
        ownership,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (crossings != expected) {
        throw std::invalid_argument(
            "scene fluid face crossings do not match their owned patches");
    }
}

} // namespace simwing::fsi::fluid
