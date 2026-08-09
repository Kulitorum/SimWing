#include "scene_fluid_opening_face_crossing.h"

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
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t crossingIdentityDomain = 0x6f70656e78666163ULL;

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
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
    if (axis == 0) return value.x;
    if (axis == 1) return value.y;
    return value.z;
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

double squaredNorm(const Vec3& value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}

double norm(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

bool vertexLess(const fluid::SceneFluidClippedVertex& first,
                const fluid::SceneFluidClippedVertex& second) {
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
    fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    bool sourceIsLowerCell = false;
};

struct BoundaryPlane {
    FaceLocation location;
    std::size_t coordinateAxis = 0;
    double coordinateMeters = 0.0;
    bool internal = false;
};

std::array<BoundaryPlane, 6> boundaryPlanes(
    const SceneFluidOpeningGridPatch& patch,
    const fluid::PeriodicCartesianGrid& grid) {
    const auto lowerGrid = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    const auto counts = grid.cellCounts();
    const Vec3 lower{
        lowerGrid.x + static_cast<double>(patch.cell.i) * spacing.x,
        lowerGrid.y + static_cast<double>(patch.cell.j) * spacing.y,
        lowerGrid.z + static_cast<double>(patch.cell.k) * spacing.z,
    };
    const Vec3 upper{
        lower.x + spacing.x,
        lower.y + spacing.y,
        lower.z + spacing.z,
    };
    return {{
        {{fluid::GridFaceAxis::X, patch.cell.i,
          patch.cell.j, patch.cell.k, false},
         0, lower.x, patch.cell.i != 0},
        {{fluid::GridFaceAxis::X, patch.cell.i + 1,
          patch.cell.j, patch.cell.k, true},
         0, upper.x, patch.cell.i + 1 < counts.x},
        {{fluid::GridFaceAxis::Y, patch.cell.i,
          patch.cell.j, patch.cell.k, false},
         1, lower.y, patch.cell.j != 0},
        {{fluid::GridFaceAxis::Y, patch.cell.i,
          patch.cell.j + 1, patch.cell.k, true},
         1, upper.y, patch.cell.j + 1 < counts.y},
        {{fluid::GridFaceAxis::Z, patch.cell.i,
          patch.cell.j, patch.cell.k, false},
         2, lower.z, patch.cell.k != 0},
        {{fluid::GridFaceAxis::Z, patch.cell.i,
          patch.cell.j, patch.cell.k + 1, true},
         2, upper.z, patch.cell.k + 1 < counts.z},
    }};
}

struct CandidateSegment {
    FaceLocation location;
    std::size_t sourcePatchIndex = 0;
    std::uint64_t sourcePointStableId = 0;
    fluid::SceneFluidClippedVertex first;
    fluid::SceneFluidClippedVertex second;
    bool gridEdgeAligned = false;
};

std::optional<CandidateSegment> candidateSegment(
    const SceneFluidOpeningGridPatch& patch,
    const std::span<const fluid::SceneFluidClippedVertex> vertices,
    const BoundaryPlane& plane) {
    std::vector<fluid::SceneFluidClippedVertex> onPlane;
    onPlane.reserve(vertices.size());
    for (const auto& vertex : vertices) {
        if (coordinate(vertex.positionMeters, plane.coordinateAxis)
            == plane.coordinateMeters) {
            onPlane.push_back(vertex);
        }
    }
    if (onPlane.size() < 2) return std::nullopt;

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
    result.sourcePointStableId = patch.sourcePointStableId;
    result.first = onPlane[firstIndex];
    result.second = onPlane[secondIndex];
    if (vertexLess(result.second, result.first)) {
        std::swap(result.first, result.second);
    }
    return result;
}

std::size_t commonBoundaryPlaneCount(
    const SceneFluidOpeningGridPatch& patch,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::SceneFluidClippedVertex& first,
    const fluid::SceneFluidClippedVertex& second) {
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

bool sameVector(const Vec3& first, const Vec3& second) {
    return first.x == second.x
        && first.y == second.y
        && first.z == second.z;
}

bool sameVertexWithinRoundoff(
    const fluid::SceneFluidClippedVertex& first,
    const fluid::SceneFluidClippedVertex& second,
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
    if (first.sourcePointStableId != second.sourcePointStableId) return false;
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
    const double tolerance =
        256.0 * std::numeric_limits<double>::epsilon() * scale;
    const bool direct = sameVertexWithinRoundoff(
                            first.first, second.first, tolerance)
        && sameVertexWithinRoundoff(
            first.second, second.second, tolerance);
    const bool reversed = sameVertexWithinRoundoff(
                              first.first, second.second, tolerance)
        && sameVertexWithinRoundoff(
            first.second, second.first, tolerance);
    return direct || reversed;
}

using FaceKey = std::tuple<std::uint8_t, std::size_t, std::size_t,
                           std::size_t, std::uint64_t>;

struct FacePair {
    FaceLocation location;
    std::optional<CandidateSegment> lower;
    std::optional<CandidateSegment> upper;
};

std::uint64_t crossingStableId(const FaceLocation& location,
                               const std::uint64_t sourcePointStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(crossingIdentityDomain);
    fingerprint.enumeration(location.axis);
    fingerprint.integer(static_cast<std::uint64_t>(location.i));
    fingerprint.integer(static_cast<std::uint64_t>(location.j));
    fingerprint.integer(static_cast<std::uint64_t>(location.k));
    fingerprint.integer(sourcePointStableId);
    return fingerprint.value();
}

Vec3 projectedDirection(const Vec3& normal,
                        const fluid::GridFaceAxis axis) {
    Vec3 result = normal;
    if (axis == fluid::GridFaceAxis::X) result.x = 0.0;
    else if (axis == fluid::GridFaceAxis::Y) result.y = 0.0;
    else result.z = 0.0;
    const double magnitude = norm(result);
    if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
        throw std::invalid_argument(
            "scene fluid opening cap is not transverse to its grid face");
    }
    result.x /= magnitude;
    result.y /= magnitude;
    result.z /= magnitude;
    return result;
}

void appendVector(Fingerprint& fingerprint, const Vec3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void appendVector(Fingerprint& fingerprint,
                  const fluid::Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t crossingsFingerprint(
    const SceneFluidOpeningFaceCrossingSet& crossings) {
    Fingerprint fingerprint;
    fingerprint.integer(crossings.version);
    fingerprint.integer(crossings.surfaceDefinitionFingerprint);
    fingerprint.integer(crossings.surfaceStateFingerprint);
    fingerprint.integer(crossings.openingCapFingerprint);
    fingerprint.integer(crossings.openingQuadratureFingerprint);
    fingerprint.integer(crossings.openingPatchFingerprint);
    fingerprint.integer(crossings.structureDefinitionFingerprint);
    fingerprint.integer(crossings.acceptedStepCount);
    fingerprint.real(crossings.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(crossings.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(crossings.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(crossings.cellCounts.z));
    appendVector(fingerprint, crossings.lowerMeters);
    appendVector(fingerprint, crossings.upperMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.candidateSegmentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.unpairedContactSegmentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.faceOwnedPatchCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.ownedStorageBytes));
    fingerprint.real(crossings.crossingLengthMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        crossings.crossings.size()));
    for (const auto& crossing : crossings.crossings) {
        fingerprint.integer(crossing.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.lowerCellPatchIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.upperCellPatchIndex));
        fingerprint.enumeration(crossing.axis);
        fingerprint.integer(static_cast<std::uint64_t>(crossing.i));
        fingerprint.integer(static_cast<std::uint64_t>(crossing.j));
        fingerprint.integer(static_cast<std::uint64_t>(crossing.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.openingIndex));
        fingerprint.integer(crossing.openingId);
        fingerprint.integer(static_cast<std::uint64_t>(
            crossing.triangleOrdinal));
        fingerprint.integer(crossing.sourcePointStableId);
        for (const auto* vertex : {&crossing.first, &crossing.second}) {
            appendVector(fingerprint, vertex->positionMeters);
            for (const double value : vertex->barycentricCoordinates) {
                fingerprint.real(value);
            }
        }
        appendVector(fingerprint, crossing.midpointMeters);
        appendVector(fingerprint,
                     crossing.negativeToPositiveDirectionInFace);
        fingerprint.real(crossing.lengthMeters);
        fingerprint.integer(crossing.negativeSideRegionId);
        fingerprint.integer(crossing.positiveSideRegionId);
        fingerprint.enumeration(crossing.role);
    }
    return fingerprint.value();
}

std::size_t storageBytes(
    const SceneFluidOpeningFaceCrossingSet& crossings) {
    std::size_t result = 0;
    if (!checkedMultiply(crossings.crossings.size(),
                         sizeof(SceneFluidOpeningFaceCrossing), result)) {
        throw std::length_error(
            "scene fluid opening-face-crossing storage overflows");
    }
    return result;
}

bool sameGrid(const SceneFluidOpeningFaceCrossingSet& crossings,
              const fluid::PeriodicCartesianGrid& grid) {
    return crossings.cellCounts == grid.cellCounts()
        && crossings.lowerMeters == grid.lowerMeters()
        && crossings.upperMeters == grid.upperMeters();
}

SceneFluidOpeningFaceCrossingSet buildCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningFaceCrossingLimits& limits) {
    SceneFluidOpeningFaceCrossingSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.openingCapFingerprint = caps.fingerprint;
    result.openingQuadratureFingerprint = quadrature.fingerprint;
    result.openingPatchFingerprint = patches.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();

    std::map<FaceKey, FacePair> pairs;
    for (std::size_t patchIndex = 0;
         patchIndex < patches.patches.size(); ++patchIndex) {
        const auto& patch = patches.patches[patchIndex];
        if (patch.ownerKind == SceneFluidOpeningPatchOwnerKind::Face) {
            ++result.faceOwnedPatchCount;
            continue;
        }
        const auto vertices = patches.verticesForPatch(patch);
        for (const auto& plane : boundaryPlanes(patch, grid)) {
            if (!plane.internal) continue;
            auto candidate = candidateSegment(patch, vertices, plane);
            if (!candidate) continue;
            if (result.candidateSegmentCount
                == limits.maximumCandidateSegments) {
                throw std::length_error(
                    "scene fluid opening-face candidates exceed their limit");
            }
            ++result.candidateSegmentCount;
            candidate->sourcePatchIndex = patchIndex;
            candidate->gridEdgeAligned = commonBoundaryPlaneCount(
                patch, grid, candidate->first, candidate->second) > 1;
            const auto key = FaceKey{
                static_cast<std::uint8_t>(plane.location.axis),
                plane.location.i, plane.location.j, plane.location.k,
                patch.sourcePointStableId};
            auto& pair = pairs[key];
            pair.location = plane.location;
            auto& destination = plane.location.sourceIsLowerCell
                ? pair.lower : pair.upper;
            if (destination) {
                throw std::invalid_argument(
                    "scene fluid opening face has duplicate cap contact");
            }
            destination = std::move(*candidate);
        }
    }

    std::set<std::uint64_t> stableIds;
    for (const auto& [key, pair] : pairs) {
        static_cast<void>(key);
        if (!pair.lower || !pair.upper) {
            ++result.unpairedContactSegmentCount;
            continue;
        }
        if (!sameSegmentWithinRoundoff(*pair.lower, *pair.upper)) {
            throw std::invalid_argument(
                "scene fluid opening face has mismatched adjacent cap clips");
        }
        if (pair.lower->gridEdgeAligned || pair.upper->gridEdgeAligned) {
            throw std::invalid_argument(
                "scene fluid opening face crossing lies on a grid edge");
        }
        if (result.crossings.size() == limits.maximumCrossings) {
            throw std::length_error(
                "scene fluid opening-face crossings exceed their limit");
        }
        std::size_t prospectiveBytes = 0;
        if (!checkedMultiply(
                result.crossings.size() + 1,
                sizeof(SceneFluidOpeningFaceCrossing), prospectiveBytes)
            || prospectiveBytes > limits.maximumCrossingBytes) {
            throw std::length_error(
                "scene fluid opening-face crossings exceed their byte limit");
        }
        const auto& lowerPatch = patches.patches[
            pair.lower->sourcePatchIndex];
        const auto& upperPatch = patches.patches[
            pair.upper->sourcePatchIndex];
        if (lowerPatch.sourcePointStableId
                != upperPatch.sourcePointStableId
            || lowerPatch.openingIndex != upperPatch.openingIndex
            || lowerPatch.openingId != upperPatch.openingId
            || lowerPatch.triangleOrdinal != upperPatch.triangleOrdinal
            || lowerPatch.negativeSideRegionId
                != upperPatch.negativeSideRegionId
            || lowerPatch.positiveSideRegionId
                != upperPatch.positiveSideRegionId
            || lowerPatch.role != upperPatch.role
            || !sameVector(lowerPatch.unitNormalNegativeToPositive,
                           upperPatch.unitNormalNegativeToPositive)) {
            throw std::invalid_argument(
                "scene fluid opening face cap metadata disagrees");
        }
        SceneFluidOpeningFaceCrossing crossing;
        crossing.stableId = crossingStableId(
            pair.location, lowerPatch.sourcePointStableId);
        crossing.lowerCellPatchIndex = pair.lower->sourcePatchIndex;
        crossing.upperCellPatchIndex = pair.upper->sourcePatchIndex;
        crossing.axis = pair.location.axis;
        crossing.i = pair.location.i;
        crossing.j = pair.location.j;
        crossing.k = pair.location.k;
        crossing.openingIndex = lowerPatch.openingIndex;
        crossing.openingId = lowerPatch.openingId;
        crossing.triangleOrdinal = lowerPatch.triangleOrdinal;
        crossing.sourcePointStableId = lowerPatch.sourcePointStableId;
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
        crossing.negativeToPositiveDirectionInFace = projectedDirection(
            lowerPatch.unitNormalNegativeToPositive, crossing.axis);
        crossing.lengthMeters = norm(subtract(
            crossing.second.positionMeters,
            crossing.first.positionMeters));
        crossing.negativeSideRegionId = lowerPatch.negativeSideRegionId;
        crossing.positiveSideRegionId = lowerPatch.positiveSideRegionId;
        crossing.role = lowerPatch.role;
        if (!std::isfinite(crossing.lengthMeters)
            || !(crossing.lengthMeters > 0.0)
            || crossing.negativeSideRegionId == invalidStableId
            || crossing.positiveSideRegionId == invalidStableId
            || crossing.negativeSideRegionId
                == crossing.positiveSideRegionId
            || !stableIds.insert(crossing.stableId).second) {
            throw std::invalid_argument(
                "scene fluid opening face crossing is invalid");
        }
        result.crossingLengthMeters += crossing.lengthMeters;
        result.crossings.push_back(crossing);
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumCrossingBytes) {
        throw std::length_error(
            "scene fluid opening-face crossings exceed their byte limit");
    }
    result.fingerprint = crossingsFingerprint(result);
    return result;
}

} // namespace

bool SceneFluidOpeningFaceCrossing::operator==(
    const SceneFluidOpeningFaceCrossing& other) const {
    return stableId == other.stableId
        && lowerCellPatchIndex == other.lowerCellPatchIndex
        && upperCellPatchIndex == other.upperCellPatchIndex
        && axis == other.axis
        && i == other.i && j == other.j && k == other.k
        && openingIndex == other.openingIndex
        && openingId == other.openingId
        && triangleOrdinal == other.triangleOrdinal
        && sourcePointStableId == other.sourcePointStableId
        && first == other.first && second == other.second
        && sameVector(midpointMeters, other.midpointMeters)
        && sameVector(negativeToPositiveDirectionInFace,
                      other.negativeToPositiveDirectionInFace)
        && lengthMeters == other.lengthMeters
        && negativeSideRegionId == other.negativeSideRegionId
        && positiveSideRegionId == other.positiveSideRegionId
        && role == other.role;
}

SceneFluidOpeningFaceCrossingSet buildSceneFluidOpeningFaceCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningFaceCrossingLimits& limits) {
    validateSceneFluidOpeningGridPatches(
        patches, surface, state, caps, quadrature, grid);
    auto result = buildCrossings(
        surface, state, caps, quadrature, patches, grid, limits);
    validateSceneFluidOpeningFaceCrossings(
        result, surface, state, caps, quadrature, patches, grid);
    return result;
}

void validateSceneFluidOpeningFaceCrossings(
    const SceneFluidOpeningFaceCrossingSet& crossings,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid) {
    validateSceneFluidOpeningGridPatches(
        patches, surface, state, caps, quadrature, grid);
    if (crossings.version != sceneFluidOpeningFaceCrossingVersion
        || crossings.fingerprint == 0
        || crossings.surfaceDefinitionFingerprint != surface.fingerprint
        || crossings.surfaceStateFingerprint != state.fingerprint
        || crossings.openingCapFingerprint != caps.fingerprint
        || crossings.openingQuadratureFingerprint != quadrature.fingerprint
        || crossings.openingPatchFingerprint != patches.fingerprint
        || crossings.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || crossings.acceptedStepCount != state.acceptedStepCount
        || crossings.simulationTimeSeconds != state.simulationTimeSeconds
        || !sameGrid(crossings, grid)) {
        throw std::invalid_argument(
            "scene fluid opening-face-crossing identity is invalid");
    }
    const SceneFluidOpeningFaceCrossingLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildCrossings(
        surface, state, caps, quadrature, patches, grid, unlimited);
    if (crossings != expected
        || crossings.ownedStorageBytes != storageBytes(crossings)
        || crossings.fingerprint != crossingsFingerprint(crossings)) {
        throw std::invalid_argument(
            "scene fluid opening-face-crossing payload is invalid");
    }
}

} // namespace simwing::fsi
