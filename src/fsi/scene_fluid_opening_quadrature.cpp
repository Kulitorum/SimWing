#include "scene_fluid_opening_quadrature.h"

#include <bit>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t pointIdentityDomain = 0x6f70656e71756164ULL;

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

std::uint64_t pointStableId(
    const StableId openingId,
    const std::size_t triangleOrdinal) {
    Fingerprint fingerprint;
    fingerprint.integer(pointIdentityDomain);
    fingerprint.integer(openingId);
    fingerprint.integer(static_cast<std::uint64_t>(triangleOrdinal));
    return fingerprint.value();
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

std::size_t storageBytes(
    const SceneFluidOpeningQuadratureSet& quadrature) {
    std::size_t openingBytes = 0;
    std::size_t pointBytes = 0;
    std::size_t totalBytes = 0;
    if (!checkedMultiply(
            quadrature.openings.size(),
            sizeof(SceneFluidOpeningQuadrature), openingBytes)
        || !checkedMultiply(
            quadrature.points.size(),
            sizeof(SceneFluidOpeningQuadraturePoint), pointBytes)
        || !checkedAdd(openingBytes, pointBytes, totalBytes)) {
        throw std::length_error(
            "scene fluid opening-quadrature storage size overflows");
    }
    return totalBytes;
}

void appendVector(Fingerprint& fingerprint, const Vec3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t quadratureFingerprint(
    const SceneFluidOpeningQuadratureSet& quadrature) {
    Fingerprint fingerprint;
    fingerprint.integer(quadrature.version);
    fingerprint.integer(quadrature.surfaceDefinitionFingerprint);
    fingerprint.integer(quadrature.surfaceStateFingerprint);
    fingerprint.integer(quadrature.openingCapFingerprint);
    fingerprint.integer(quadrature.structureDefinitionFingerprint);
    fingerprint.integer(quadrature.acceptedStepCount);
    fingerprint.real(quadrature.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        quadrature.ownedStorageBytes));
    fingerprint.real(quadrature.totalAreaSquareMeters);
    fingerprint.real(
        quadrature.totalSurfaceSweepRateCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint64_t>(
        quadrature.openings.size()));
    for (const auto& opening : quadrature.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.openingIndex));
        fingerprint.integer(opening.openingId);
        fingerprint.integer(opening.negativeSideRegionId);
        fingerprint.integer(opening.positiveSideRegionId);
        fingerprint.enumeration(opening.role);
        fingerprint.integer(static_cast<std::uint64_t>(opening.firstPoint));
        fingerprint.integer(static_cast<std::uint64_t>(opening.pointCount));
        appendVector(fingerprint, opening.unitNormalNegativeToPositive);
        appendVector(fingerprint, opening.centroidMeters);
        fingerprint.real(opening.areaSquareMeters);
        fingerprint.real(
            opening.surfaceSweepRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        quadrature.points.size()));
    for (const auto& point : quadrature.points) {
        fingerprint.integer(point.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(point.openingIndex));
        fingerprint.integer(point.openingId);
        fingerprint.integer(static_cast<std::uint64_t>(
            point.triangleOrdinal));
        for (const std::size_t vertex : point.vertexIndices) {
            fingerprint.integer(static_cast<std::uint64_t>(vertex));
        }
        fingerprint.integer(point.negativeSideRegionId);
        fingerprint.integer(point.positiveSideRegionId);
        fingerprint.enumeration(point.role);
        appendVector(fingerprint, point.positionMeters);
        appendVector(fingerprint, point.velocityMetersPerSecond);
        appendVector(fingerprint, point.unitNormalNegativeToPositive);
        fingerprint.real(point.areaSquareMeters);
        fingerprint.real(point.surfaceSweepRateCubicMetersPerSecond);
    }
    return fingerprint.value();
}

SceneFluidOpeningQuadratureSet buildQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureLimits& limits) {
    if (caps.caps.size() > limits.maximumOpenings
        || caps.triangles.size() > limits.maximumPoints) {
        throw std::length_error(
            "scene fluid opening quadrature exceeds its count limit");
    }

    SceneFluidOpeningQuadratureSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.openingCapFingerprint = caps.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.openings.reserve(caps.caps.size());
    result.points.reserve(caps.triangles.size());

    std::set<std::uint64_t> stablePointIds;
    for (const auto& cap : caps.caps) {
        SceneFluidOpeningQuadrature opening;
        opening.openingIndex = cap.openingIndex;
        opening.openingId = cap.openingId;
        opening.negativeSideRegionId =
            surface.regions[cap.negativeSideRegionIndex].id;
        opening.positiveSideRegionId =
            surface.regions[cap.positiveSideRegionIndex].id;
        opening.role = cap.role;
        opening.firstPoint = result.points.size();
        opening.pointCount = cap.triangleCount;
        opening.unitNormalNegativeToPositive =
            cap.unitNormalNegativeToPositive;
        opening.centroidMeters = cap.centroidMeters;
        opening.areaSquareMeters = cap.areaSquareMeters;

        for (std::size_t localIndex = 0;
             localIndex < cap.triangleCount; ++localIndex) {
            const auto& triangle = caps.triangles[
                cap.firstTriangle + localIndex];
            Vec3 velocity;
            for (const std::size_t vertex : triangle.vertexIndices) {
                const auto& source = state.vertices[vertex]
                    .velocityMetersPerSecond;
                velocity.x += source.x / 3.0;
                velocity.y += source.y / 3.0;
                velocity.z += source.z / 3.0;
            }
            const double sweepRate = triangle.areaSquareMeters
                * dot(velocity, triangle.unitNormalNegativeToPositive);
            const std::uint64_t stableId = pointStableId(
                cap.openingId, triangle.triangleOrdinal);
            if (!finite(velocity) || !std::isfinite(sweepRate)
                || !stablePointIds.insert(stableId).second) {
                throw std::invalid_argument(
                    "scene fluid opening quadrature kinematics or identity is invalid");
            }
            result.points.push_back({
                stableId,
                cap.openingIndex,
                cap.openingId,
                triangle.triangleOrdinal,
                triangle.vertexIndices,
                opening.negativeSideRegionId,
                opening.positiveSideRegionId,
                cap.role,
                triangle.centroidMeters,
                velocity,
                triangle.unitNormalNegativeToPositive,
                triangle.areaSquareMeters,
                sweepRate,
            });
            opening.surfaceSweepRateCubicMetersPerSecond += sweepRate;
        }
        if (!std::isfinite(
                opening.surfaceSweepRateCubicMetersPerSecond)) {
            throw std::overflow_error(
                "scene fluid opening sweep rate is not finite");
        }
        result.totalAreaSquareMeters += opening.areaSquareMeters;
        result.totalSurfaceSweepRateCubicMetersPerSecond +=
            opening.surfaceSweepRateCubicMetersPerSecond;
        result.openings.push_back(opening);
    }
    if (!std::isfinite(result.totalAreaSquareMeters)
        || !std::isfinite(
            result.totalSurfaceSweepRateCubicMetersPerSecond)
        || result.totalAreaSquareMeters != caps.totalAreaSquareMeters) {
        throw std::overflow_error(
            "scene fluid opening quadrature aggregate is invalid");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumQuadratureBytes) {
        throw std::length_error(
            "scene fluid opening quadrature exceeds its byte limit");
    }
    result.fingerprint = quadratureFingerprint(result);
    return result;
}

} // namespace

bool SceneFluidOpeningQuadraturePoint::operator==(
    const SceneFluidOpeningQuadraturePoint& other) const {
    return stableId == other.stableId
        && openingIndex == other.openingIndex
        && openingId == other.openingId
        && triangleOrdinal == other.triangleOrdinal
        && vertexIndices == other.vertexIndices
        && negativeSideRegionId == other.negativeSideRegionId
        && positiveSideRegionId == other.positiveSideRegionId
        && role == other.role
        && same(positionMeters, other.positionMeters)
        && same(velocityMetersPerSecond, other.velocityMetersPerSecond)
        && same(unitNormalNegativeToPositive,
                other.unitNormalNegativeToPositive)
        && areaSquareMeters == other.areaSquareMeters
        && surfaceSweepRateCubicMetersPerSecond
            == other.surfaceSweepRateCubicMetersPerSecond;
}

bool SceneFluidOpeningQuadrature::operator==(
    const SceneFluidOpeningQuadrature& other) const {
    return openingIndex == other.openingIndex
        && openingId == other.openingId
        && negativeSideRegionId == other.negativeSideRegionId
        && positiveSideRegionId == other.positiveSideRegionId
        && role == other.role
        && firstPoint == other.firstPoint
        && pointCount == other.pointCount
        && same(unitNormalNegativeToPositive,
                other.unitNormalNegativeToPositive)
        && same(centroidMeters, other.centroidMeters)
        && areaSquareMeters == other.areaSquareMeters
        && surfaceSweepRateCubicMetersPerSecond
            == other.surfaceSweepRateCubicMetersPerSecond;
}

SceneFluidOpeningQuadratureSet buildSceneFluidOpeningQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureLimits& limits) {
    validateSceneFluidOpeningCaps(caps, surface, state);
    auto result = buildQuadrature(surface, state, caps, limits);
    validateSceneFluidOpeningQuadrature(result, surface, state, caps);
    return result;
}

void validateSceneFluidOpeningQuadrature(
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps) {
    validateSceneFluidOpeningCaps(caps, surface, state);
    if (quadrature.version != sceneFluidOpeningQuadratureVersion
        || quadrature.fingerprint == 0
        || quadrature.surfaceDefinitionFingerprint != surface.fingerprint
        || quadrature.surfaceStateFingerprint != state.fingerprint
        || quadrature.openingCapFingerprint != caps.fingerprint
        || quadrature.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || quadrature.acceptedStepCount != state.acceptedStepCount
        || quadrature.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid opening-quadrature identity is invalid");
    }
    const SceneFluidOpeningQuadratureLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildQuadrature(surface, state, caps, unlimited);
    if (quadrature != expected
        || quadrature.ownedStorageBytes != storageBytes(quadrature)
        || quadrature.fingerprint != quadratureFingerprint(quadrature)) {
        throw std::invalid_argument(
            "scene fluid opening-quadrature payload is invalid");
    }
}

} // namespace simwing::fsi
