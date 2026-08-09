#include "scene_fluid_opening_patch.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t patchIdentityDomain = 0x6f70656e70617463ULL;

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

struct AxisCellRange {
    std::size_t first = 0;
    std::size_t last = 0;
};

struct FaceLocation {
    SceneFluidOpeningPatchFaceAxis axis =
        SceneFluidOpeningPatchFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    bool sourceIsLowerCell = false;
};

struct CandidateClip {
    std::size_t cellIndex = 0;
    fluid::GridCellCoordinate cell;
    fluid::SceneFluidTriangleBoxClip clip;
};

struct FacePair {
    std::optional<CandidateClip> lowerCell;
    std::optional<CandidateClip> upperCell;
};

struct OwnedClip {
    SceneFluidOpeningPatchOwnerKind ownerKind =
        SceneFluidOpeningPatchOwnerKind::Cell;
    std::size_t cellIndex = 0;
    fluid::GridCellCoordinate cell;
    SceneFluidOpeningPatchFaceAxis faceAxis =
        SceneFluidOpeningPatchFaceAxis::X;
    std::size_t faceI = 0;
    std::size_t faceJ = 0;
    std::size_t faceK = 0;
    fluid::SceneFluidTriangleBoxClip clip;
};

using FaceKey = std::tuple<SceneFluidOpeningPatchFaceAxis,
                           std::size_t, std::size_t, std::size_t>;

bool same(const Vec3& first, const Vec3& second) {
    return first.x == second.x
        && first.y == second.y
        && first.z == second.z;
}

bool finite(const Vec3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
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

double boundedCoordinate(const double value,
                         const double lower,
                         const double upper,
                         const double tolerance) {
    if (!std::isfinite(value)
        || value < lower - tolerance
        || value > upper + tolerance) {
        throw std::invalid_argument(
            "scene fluid opening patch lies outside the grid domain");
    }
    return std::clamp(value, lower, upper);
}

AxisCellRange axisCellRange(const double minimum,
                            const double maximum,
                            const double domainLower,
                            const double domainUpper,
                            const double spacing,
                            const std::size_t cellCount,
                            const double padding) {
    const double paddedMinimum = std::max(domainLower, minimum - padding);
    const double paddedMaximum = std::min(domainUpper, maximum + padding);
    if (!std::isfinite(paddedMinimum) || !std::isfinite(paddedMaximum)
        || paddedMinimum > paddedMaximum) {
        throw std::invalid_argument(
            "scene fluid opening-patch padding produced invalid bounds");
    }
    const double firstCoordinate = std::nextafter(
        (paddedMinimum - domainLower) / spacing,
        -std::numeric_limits<double>::infinity());
    const double lastCoordinate = std::nextafter(
        (paddedMaximum - domainLower) / spacing,
        std::numeric_limits<double>::infinity());
    const auto clampedIndex = [cellCount](const double value) {
        if (value <= 0.0) {
            return std::size_t{0};
        }
        if (value >= static_cast<double>(cellCount)) {
            return cellCount - 1;
        }
        return static_cast<std::size_t>(std::floor(value));
    };
    return {clampedIndex(firstCoordinate), clampedIndex(lastCoordinate)};
}

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

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

void validateSettings(const SceneFluidOpeningGridPatchSettings& settings) {
    if (!std::isfinite(settings.boundingPaddingMeters)
        || settings.boundingPaddingMeters < 0.0
        || !std::isfinite(settings.domainBoundaryToleranceMeters)
        || settings.domainBoundaryToleranceMeters < 0.0
        || !std::isfinite(settings.absoluteAreaToleranceSquareMeters)
        || settings.absoluteAreaToleranceSquareMeters < 0.0
        || !std::isfinite(settings.relativeAreaTolerance)
        || settings.relativeAreaTolerance < 0.0
        || !std::isfinite(
            settings.absoluteSweepRateToleranceCubicMetersPerSecond)
        || settings.absoluteSweepRateToleranceCubicMetersPerSecond < 0.0
        || !std::isfinite(settings.relativeSweepRateTolerance)
        || settings.relativeSweepRateTolerance < 0.0) {
        throw std::invalid_argument(
            "scene fluid opening-patch settings are invalid");
    }
}

bool withinTolerance(const double first,
                     const double second,
                     const double absoluteTolerance,
                     const double relativeTolerance) {
    const double tolerance = absoluteTolerance
        + relativeTolerance * std::max(std::abs(first), std::abs(second));
    return std::isfinite(first) && std::isfinite(second)
        && std::abs(first - second) <= tolerance;
}

FaceLocation faceLocation(const CandidateClip& candidate,
                          const fluid::GridCellCounts counts) {
    const std::uint8_t mask = candidate.clip.coincidentBoundaryPlanes;
    if (mask == 0 || (mask & (mask - 1U)) != 0) {
        throw std::invalid_argument(
            "positive-area opening patch must coincide with at most one cell plane");
    }
    if (mask == fluid::CellBoundaryXMinus) {
        if (candidate.cell.i == 0) {
            throw std::invalid_argument(
                "opening patch on periodic X boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::X,
                candidate.cell.i, candidate.cell.j, candidate.cell.k, false};
    }
    if (mask == fluid::CellBoundaryXPlus) {
        if (candidate.cell.i + 1 >= counts.x) {
            throw std::invalid_argument(
                "opening patch on periodic X boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::X,
                candidate.cell.i + 1, candidate.cell.j, candidate.cell.k,
                true};
    }
    if (mask == fluid::CellBoundaryYMinus) {
        if (candidate.cell.j == 0) {
            throw std::invalid_argument(
                "opening patch on periodic Y boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::Y,
                candidate.cell.i, candidate.cell.j, candidate.cell.k, false};
    }
    if (mask == fluid::CellBoundaryYPlus) {
        if (candidate.cell.j + 1 >= counts.y) {
            throw std::invalid_argument(
                "opening patch on periodic Y boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::Y,
                candidate.cell.i, candidate.cell.j + 1, candidate.cell.k,
                true};
    }
    if (mask == fluid::CellBoundaryZMinus) {
        if (candidate.cell.k == 0) {
            throw std::invalid_argument(
                "opening patch on periodic Z boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::Z,
                candidate.cell.i, candidate.cell.j, candidate.cell.k, false};
    }
    if (mask == fluid::CellBoundaryZPlus) {
        if (candidate.cell.k + 1 >= counts.z) {
            throw std::invalid_argument(
                "opening patch on periodic Z boundary has no explicit image owner");
        }
        return {SceneFluidOpeningPatchFaceAxis::Z,
                candidate.cell.i, candidate.cell.j, candidate.cell.k + 1,
                true};
    }
    throw std::invalid_argument(
        "opening patch has an unknown coincident cell plane");
}

bool equivalentFaceClips(
    const fluid::SceneFluidTriangleBoxClip& first,
    const fluid::SceneFluidTriangleBoxClip& second) {
    return first.dimension == second.dimension
        && first.areaSquareMeters == second.areaSquareMeters
        && same(first.centroidMeters, second.centroidMeters)
        && first.centroidBarycentricCoordinates
            == second.centroidBarycentricCoordinates
        && first.vertices == second.vertices;
}

std::uint64_t patchStableId(
    const std::uint64_t sourcePointStableId,
    const OwnedClip& owner) {
    Fingerprint fingerprint;
    fingerprint.integer(patchIdentityDomain);
    fingerprint.integer(sourcePointStableId);
    fingerprint.enumeration(owner.ownerKind);
    if (owner.ownerKind == SceneFluidOpeningPatchOwnerKind::Cell) {
        fingerprint.integer(static_cast<std::uint64_t>(owner.cellIndex));
    } else {
        fingerprint.enumeration(owner.faceAxis);
        fingerprint.integer(static_cast<std::uint64_t>(owner.faceI));
        fingerprint.integer(static_cast<std::uint64_t>(owner.faceJ));
        fingerprint.integer(static_cast<std::uint64_t>(owner.faceK));
    }
    return fingerprint.value();
}

auto ownerKey(const OwnedClip& owner) {
    return std::tuple{
        owner.ownerKind,
        owner.ownerKind == SceneFluidOpeningPatchOwnerKind::Cell
            ? owner.cellIndex : std::size_t{0},
        owner.faceAxis, owner.faceI, owner.faceJ, owner.faceK};
}

Vec3 interpolatedVelocity(
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningQuadraturePoint& point,
    const std::array<double, 3>& barycentric) {
    Vec3 result;
    for (std::size_t corner = 0; corner < 3; ++corner) {
        const Vec3& velocity = state.vertices[
            point.vertexIndices[corner]].velocityMetersPerSecond;
        result.x += barycentric[corner] * velocity.x;
        result.y += barycentric[corner] * velocity.y;
        result.z += barycentric[corner] * velocity.z;
    }
    return result;
}

std::size_t storageBytesForCounts(const std::size_t rangeCount,
                                  const std::size_t patchCount,
                                  const std::size_t vertexCount) {
    std::size_t rangeBytes = 0;
    std::size_t patchBytes = 0;
    std::size_t vertexBytes = 0;
    std::size_t firstTotal = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            rangeCount,
            sizeof(SceneFluidOpeningGridPointRange), rangeBytes)
        || !checkedMultiply(
            patchCount,
            sizeof(SceneFluidOpeningGridPatch), patchBytes)
        || !checkedMultiply(
            vertexCount,
            sizeof(fluid::SceneFluidClippedVertex), vertexBytes)
        || !checkedAdd(rangeBytes, patchBytes, firstTotal)
        || !checkedAdd(firstTotal, vertexBytes, total)) {
        throw std::length_error(
            "scene fluid opening-patch storage size overflows");
    }
    return total;
}

std::size_t storageBytes(
    const SceneFluidOpeningGridPatchSet& patches) {
    return storageBytesForCounts(
        patches.pointRanges.size(), patches.patches.size(),
        patches.vertices.size());
}

void appendVector(Fingerprint& fingerprint, const Vec3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t patchFingerprint(
    const SceneFluidOpeningGridPatchSet& patches) {
    Fingerprint fingerprint;
    fingerprint.integer(patches.version);
    fingerprint.integer(patches.surfaceDefinitionFingerprint);
    fingerprint.integer(patches.surfaceStateFingerprint);
    fingerprint.integer(patches.openingCapFingerprint);
    fingerprint.integer(patches.openingQuadratureFingerprint);
    fingerprint.integer(patches.structureDefinitionFingerprint);
    fingerprint.integer(patches.acceptedStepCount);
    fingerprint.real(patches.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(patches.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(patches.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(patches.cellCounts.z));
    for (const double value : {
             patches.lowerMeters.x, patches.lowerMeters.y,
             patches.lowerMeters.z, patches.upperMeters.x,
             patches.upperMeters.y, patches.upperMeters.z,
             patches.settings.boundingPaddingMeters,
             patches.settings.domainBoundaryToleranceMeters,
             patches.settings.absoluteAreaToleranceSquareMeters,
             patches.settings.relativeAreaTolerance,
             patches.settings.absoluteSweepRateToleranceCubicMetersPerSecond,
             patches.settings.relativeSweepRateTolerance,
             patches.totalAreaSquareMeters,
             patches.totalSurfaceSweepRateCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        patches.candidateCellCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        patches.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        patches.pointRanges.size()));
    for (const auto& range : patches.pointRanges) {
        fingerprint.integer(range.sourcePointStableId);
        fingerprint.integer(static_cast<std::uint64_t>(range.firstPatch));
        fingerprint.integer(static_cast<std::uint64_t>(range.patchCount));
        fingerprint.real(range.areaSquareMeters);
        fingerprint.real(range.surfaceSweepRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(patches.patches.size()));
    for (const auto& patch : patches.patches) {
        fingerprint.integer(patch.stableId);
        fingerprint.integer(patch.sourcePointStableId);
        fingerprint.integer(static_cast<std::uint64_t>(patch.openingIndex));
        fingerprint.integer(patch.openingId);
        fingerprint.integer(static_cast<std::uint64_t>(
            patch.triangleOrdinal));
        fingerprint.integer(patch.negativeSideRegionId);
        fingerprint.integer(patch.positiveSideRegionId);
        fingerprint.enumeration(patch.role);
        fingerprint.enumeration(patch.ownerKind);
        fingerprint.integer(static_cast<std::uint64_t>(patch.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(patch.cell.k));
        fingerprint.enumeration(patch.faceAxis);
        fingerprint.integer(static_cast<std::uint64_t>(patch.faceI));
        fingerprint.integer(static_cast<std::uint64_t>(patch.faceJ));
        fingerprint.integer(static_cast<std::uint64_t>(patch.faceK));
        fingerprint.integer(static_cast<std::uint64_t>(patch.firstVertex));
        fingerprint.integer(static_cast<std::uint64_t>(patch.vertexCount));
        appendVector(fingerprint, patch.centroidMeters);
        for (const double value : patch.centroidBarycentricCoordinates) {
            fingerprint.real(value);
        }
        appendVector(fingerprint, patch.velocityMetersPerSecond);
        appendVector(fingerprint, patch.unitNormalNegativeToPositive);
        fingerprint.real(patch.areaSquareMeters);
        fingerprint.real(patch.surfaceSweepRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(patches.vertices.size()));
    for (const auto& vertex : patches.vertices) {
        appendVector(fingerprint, vertex.positionMeters);
        for (const double value : vertex.barycentricCoordinates) {
            fingerprint.real(value);
        }
    }
    return fingerprint.value();
}

SceneFluidOpeningGridPatchSet buildPatches(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningGridPatchSettings& settings,
    const SceneFluidOpeningGridPatchLimits& limits) {
    validateSettings(settings);
    SceneFluidOpeningGridPatchSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.openingCapFingerprint = caps.fingerprint;
    result.openingQuadratureFingerprint = quadrature.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    if (storageBytesForCounts(quadrature.points.size(), 0, 0)
        > limits.maximumPatchBytes) {
        throw std::length_error(
            "scene fluid opening-patch point ranges exceed their byte limit");
    }
    result.pointRanges.reserve(quadrature.points.size());

    const fluid::Vector3 spacing = grid.cellSpacingMeters();
    std::set<std::uint64_t> stableIds;
    for (const auto& point : quadrature.points) {
        std::array<Vec3, 3> positions;
        const double infinity = std::numeric_limits<double>::infinity();
        Vec3 minimum{infinity, infinity, infinity};
        Vec3 maximum{-infinity, -infinity, -infinity};
        for (std::size_t corner = 0; corner < 3; ++corner) {
            positions[corner] = state.vertices[
                point.vertexIndices[corner]].positionMeters;
            positions[corner].x = boundedCoordinate(
                positions[corner].x, result.lowerMeters.x,
                result.upperMeters.x,
                settings.domainBoundaryToleranceMeters);
            positions[corner].y = boundedCoordinate(
                positions[corner].y, result.lowerMeters.y,
                result.upperMeters.y,
                settings.domainBoundaryToleranceMeters);
            positions[corner].z = boundedCoordinate(
                positions[corner].z, result.lowerMeters.z,
                result.upperMeters.z,
                settings.domainBoundaryToleranceMeters);
            minimum.x = std::min(minimum.x, positions[corner].x);
            minimum.y = std::min(minimum.y, positions[corner].y);
            minimum.z = std::min(minimum.z, positions[corner].z);
            maximum.x = std::max(maximum.x, positions[corner].x);
            maximum.y = std::max(maximum.y, positions[corner].y);
            maximum.z = std::max(maximum.z, positions[corner].z);
        }
        const AxisCellRange xRange = axisCellRange(
            minimum.x, maximum.x, result.lowerMeters.x,
            result.upperMeters.x, spacing.x, result.cellCounts.x,
            settings.boundingPaddingMeters);
        const AxisCellRange yRange = axisCellRange(
            minimum.y, maximum.y, result.lowerMeters.y,
            result.upperMeters.y, spacing.y, result.cellCounts.y,
            settings.boundingPaddingMeters);
        const AxisCellRange zRange = axisCellRange(
            minimum.z, maximum.z, result.lowerMeters.z,
            result.upperMeters.z, spacing.z, result.cellCounts.z,
            settings.boundingPaddingMeters);
        const std::size_t xCount = xRange.last - xRange.first + 1;
        const std::size_t yCount = yRange.last - yRange.first + 1;
        const std::size_t zCount = zRange.last - zRange.first + 1;
        std::size_t xyCount = 0;
        std::size_t candidateCount = 0;
        std::size_t newCandidateTotal = 0;
        if (!checkedMultiply(xCount, yCount, xyCount)
            || !checkedMultiply(xyCount, zCount, candidateCount)
            || !checkedAdd(result.candidateCellCount, candidateCount,
                           newCandidateTotal)
            || newCandidateTotal > limits.maximumCandidateCells) {
            throw std::length_error(
                "scene fluid opening-patch candidates exceed their limit");
        }
        result.candidateCellCount = newCandidateTotal;

        std::vector<OwnedClip> owned;
        std::map<FaceKey, FacePair> facePairs;
        for (std::size_t k = zRange.first; k <= zRange.last; ++k) {
            for (std::size_t j = yRange.first; j <= yRange.last; ++j) {
                for (std::size_t i = xRange.first; i <= xRange.last; ++i) {
                    const Vec3 lower{
                        result.lowerMeters.x + static_cast<double>(i)
                            * spacing.x,
                        result.lowerMeters.y + static_cast<double>(j)
                            * spacing.y,
                        result.lowerMeters.z + static_cast<double>(k)
                            * spacing.z,
                    };
                    const Vec3 upper{
                        lower.x + spacing.x,
                        lower.y + spacing.y,
                        lower.z + spacing.z,
                    };
                    auto clip = fluid::clipSceneFluidTriangleToAxisAlignedBox(
                        positions, lower, upper);
                    if (!clip
                        || clip->dimension
                            != fluid::SceneFluidPatchDimension::Area) {
                        continue;
                    }
                    CandidateClip candidate{
                        grid.cellIndex(i, j, k), {i, j, k},
                        std::move(*clip)};
                    if (candidate.clip.coincidentBoundaryPlanes == 0) {
                        OwnedClip owner;
                        owner.cellIndex = candidate.cellIndex;
                        owner.cell = candidate.cell;
                        owner.clip = std::move(candidate.clip);
                        owned.push_back(std::move(owner));
                        continue;
                    }
                    const FaceLocation location = faceLocation(
                        candidate, result.cellCounts);
                    auto& pair = facePairs[{
                        location.axis, location.i, location.j, location.k}];
                    auto& destination = location.sourceIsLowerCell
                        ? pair.lowerCell : pair.upperCell;
                    if (destination) {
                        throw std::invalid_argument(
                            "scene fluid opening face patch has duplicate ownership");
                    }
                    destination = std::move(candidate);
                }
            }
        }
        for (auto& [key, pair] : facePairs) {
            if (!pair.lowerCell || !pair.upperCell
                || !equivalentFaceClips(
                    pair.lowerCell->clip, pair.upperCell->clip)) {
                throw std::invalid_argument(
                    "scene fluid opening face patch lacks a matching adjacent owner");
            }
            OwnedClip owner;
            owner.ownerKind = SceneFluidOpeningPatchOwnerKind::Face;
            owner.cellIndex = pair.lowerCell->cellIndex;
            owner.cell = pair.lowerCell->cell;
            owner.faceAxis = std::get<0>(key);
            owner.faceI = std::get<1>(key);
            owner.faceJ = std::get<2>(key);
            owner.faceK = std::get<3>(key);
            owner.clip = std::move(pair.lowerCell->clip);
            owned.push_back(std::move(owner));
        }
        std::ranges::sort(owned, [](const auto& first, const auto& second) {
            return ownerKey(first) < ownerKey(second);
        });

        SceneFluidOpeningGridPointRange range;
        range.sourcePointStableId = point.stableId;
        range.firstPatch = result.patches.size();
        for (auto& owner : owned) {
            if (result.patches.size() == limits.maximumPatches) {
                throw std::length_error(
                    "scene fluid opening patches exceed their count limit");
            }
            std::size_t newPatchCount = 0;
            std::size_t newVertexCount = 0;
            if (!checkedAdd(result.patches.size(), 1, newPatchCount)
                || !checkedAdd(result.vertices.size(), owner.clip.vertices.size(),
                            newVertexCount)
                || newVertexCount > limits.maximumVertices) {
                throw std::length_error(
                    "scene fluid opening-patch vertices exceed their limit");
            }
            if (storageBytesForCounts(
                    quadrature.points.size(), newPatchCount,
                    newVertexCount) > limits.maximumPatchBytes) {
                throw std::length_error(
                    "scene fluid opening-patch result exceeds its byte limit");
            }
            const Vec3 velocity = interpolatedVelocity(
                state, point, owner.clip.centroidBarycentricCoordinates);
            const double sweepRate = owner.clip.areaSquareMeters
                * dot(velocity, point.unitNormalNegativeToPositive);
            const std::uint64_t stableId = patchStableId(
                point.stableId, owner);
            if (!finite(velocity) || !std::isfinite(sweepRate)
                || !stableIds.insert(stableId).second) {
                throw std::invalid_argument(
                    "scene fluid opening patch kinematics or identity is invalid");
            }
            SceneFluidOpeningGridPatch patch;
            patch.stableId = stableId;
            patch.sourcePointStableId = point.stableId;
            patch.openingIndex = point.openingIndex;
            patch.openingId = point.openingId;
            patch.triangleOrdinal = point.triangleOrdinal;
            patch.negativeSideRegionId = point.negativeSideRegionId;
            patch.positiveSideRegionId = point.positiveSideRegionId;
            patch.role = point.role;
            patch.ownerKind = owner.ownerKind;
            patch.cellIndex = owner.cellIndex;
            patch.cell = owner.cell;
            patch.faceAxis = owner.faceAxis;
            patch.faceI = owner.faceI;
            patch.faceJ = owner.faceJ;
            patch.faceK = owner.faceK;
            patch.firstVertex = result.vertices.size();
            patch.vertexCount = owner.clip.vertices.size();
            patch.centroidMeters = owner.clip.centroidMeters;
            patch.centroidBarycentricCoordinates =
                owner.clip.centroidBarycentricCoordinates;
            patch.velocityMetersPerSecond = velocity;
            patch.unitNormalNegativeToPositive =
                point.unitNormalNegativeToPositive;
            patch.areaSquareMeters = owner.clip.areaSquareMeters;
            patch.surfaceSweepRateCubicMetersPerSecond = sweepRate;
            result.vertices.insert(
                result.vertices.end(),
                owner.clip.vertices.begin(), owner.clip.vertices.end());
            result.patches.push_back(patch);
            range.areaSquareMeters += patch.areaSquareMeters;
            range.surfaceSweepRateCubicMetersPerSecond += sweepRate;
        }
        range.patchCount = result.patches.size() - range.firstPatch;
        if (range.patchCount == 0
            || !withinTolerance(
                range.areaSquareMeters, point.areaSquareMeters,
                settings.absoluteAreaToleranceSquareMeters,
                settings.relativeAreaTolerance)
            || !withinTolerance(
                range.surfaceSweepRateCubicMetersPerSecond,
                point.surfaceSweepRateCubicMetersPerSecond,
                settings.absoluteSweepRateToleranceCubicMetersPerSecond,
                settings.relativeSweepRateTolerance)) {
            throw std::invalid_argument(
                "scene fluid opening-patch area or sweep does not close its source triangle");
        }
        result.totalAreaSquareMeters += range.areaSquareMeters;
        result.totalSurfaceSweepRateCubicMetersPerSecond +=
            range.surfaceSweepRateCubicMetersPerSecond;
        result.pointRanges.push_back(range);
    }
    if (!withinTolerance(
            result.totalAreaSquareMeters,
            quadrature.totalAreaSquareMeters,
            settings.absoluteAreaToleranceSquareMeters,
            settings.relativeAreaTolerance)
        || !withinTolerance(
            result.totalSurfaceSweepRateCubicMetersPerSecond,
            quadrature.totalSurfaceSweepRateCubicMetersPerSecond,
            settings.absoluteSweepRateToleranceCubicMetersPerSecond,
            settings.relativeSweepRateTolerance)) {
        throw std::invalid_argument(
            "scene fluid opening-patch global area or sweep does not close");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumPatchBytes) {
        throw std::length_error(
            "scene fluid opening-patch result exceeds its byte limit");
    }
    result.fingerprint = patchFingerprint(result);
    return result;
}

bool sameGrid(const SceneFluidOpeningGridPatchSet& patches,
              const fluid::PeriodicCartesianGrid& grid) {
    return patches.cellCounts == grid.cellCounts()
        && patches.lowerMeters == grid.lowerMeters()
        && patches.upperMeters == grid.upperMeters();
}

} // namespace

bool SceneFluidOpeningGridPatch::operator==(
    const SceneFluidOpeningGridPatch& other) const {
    return stableId == other.stableId
        && sourcePointStableId == other.sourcePointStableId
        && openingIndex == other.openingIndex
        && openingId == other.openingId
        && triangleOrdinal == other.triangleOrdinal
        && negativeSideRegionId == other.negativeSideRegionId
        && positiveSideRegionId == other.positiveSideRegionId
        && role == other.role
        && ownerKind == other.ownerKind
        && cellIndex == other.cellIndex
        && cell == other.cell
        && faceAxis == other.faceAxis
        && faceI == other.faceI
        && faceJ == other.faceJ
        && faceK == other.faceK
        && firstVertex == other.firstVertex
        && vertexCount == other.vertexCount
        && same(centroidMeters, other.centroidMeters)
        && centroidBarycentricCoordinates
            == other.centroidBarycentricCoordinates
        && same(velocityMetersPerSecond, other.velocityMetersPerSecond)
        && same(unitNormalNegativeToPositive,
                other.unitNormalNegativeToPositive)
        && areaSquareMeters == other.areaSquareMeters
        && surfaceSweepRateCubicMetersPerSecond
            == other.surfaceSweepRateCubicMetersPerSecond;
}

std::span<const SceneFluidOpeningGridPatch>
SceneFluidOpeningGridPatchSet::patchesForPoint(
    const SceneFluidOpeningGridPointRange& range) const {
    if (range.firstPatch > patches.size()
        || range.patchCount > patches.size() - range.firstPatch) {
        throw std::out_of_range(
            "scene fluid opening-patch range is out of bounds");
    }
    return std::span<const SceneFluidOpeningGridPatch>(patches).subspan(
        range.firstPatch, range.patchCount);
}

std::span<const fluid::SceneFluidClippedVertex>
SceneFluidOpeningGridPatchSet::verticesForPatch(
    const SceneFluidOpeningGridPatch& patch) const {
    if (patch.firstVertex > vertices.size()
        || patch.vertexCount > vertices.size() - patch.firstVertex) {
        throw std::out_of_range(
            "scene fluid opening-patch vertex range is out of bounds");
    }
    return std::span<const fluid::SceneFluidClippedVertex>(vertices).subspan(
        patch.firstVertex, patch.vertexCount);
}

SceneFluidOpeningGridPatchSet buildSceneFluidOpeningGridPatches(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningGridPatchSettings& settings,
    const SceneFluidOpeningGridPatchLimits& limits) {
    validateSceneFluidOpeningQuadrature(
        quadrature, surface, state, caps);
    auto result = buildPatches(
        surface, state, caps, quadrature, grid, settings, limits);
    validateSceneFluidOpeningGridPatches(
        result, surface, state, caps, quadrature, grid);
    return result;
}

void validateSceneFluidOpeningGridPatches(
    const SceneFluidOpeningGridPatchSet& patches,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const fluid::PeriodicCartesianGrid& grid) {
    validateSceneFluidOpeningQuadrature(
        quadrature, surface, state, caps);
    validateSettings(patches.settings);
    if (patches.version != sceneFluidOpeningGridPatchVersion
        || patches.fingerprint == 0
        || patches.surfaceDefinitionFingerprint != surface.fingerprint
        || patches.surfaceStateFingerprint != state.fingerprint
        || patches.openingCapFingerprint != caps.fingerprint
        || patches.openingQuadratureFingerprint != quadrature.fingerprint
        || patches.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || patches.acceptedStepCount != state.acceptedStepCount
        || patches.simulationTimeSeconds != state.simulationTimeSeconds
        || !sameGrid(patches, grid)) {
        throw std::invalid_argument(
            "scene fluid opening-patch identity is invalid");
    }
    const SceneFluidOpeningGridPatchLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildPatches(
        surface, state, caps, quadrature, grid,
        patches.settings, unlimited);
    if (patches != expected
        || patches.ownedStorageBytes != storageBytes(patches)
        || patches.fingerprint != patchFingerprint(patches)) {
        throw std::invalid_argument(
            "scene fluid opening-patch payload is invalid");
    }
}

} // namespace simwing::fsi
