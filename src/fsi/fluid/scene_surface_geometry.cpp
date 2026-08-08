#include "fluid/scene_surface_geometry.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
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

bool checkedAdd(std::size_t& total, const std::size_t value) {
    if (value > std::numeric_limits<std::size_t>::max() - total) {
        return false;
    }
    total += value;
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

bool checkedBytes(std::size_t& total,
                  const std::size_t count,
                  const std::size_t itemBytes) {
    std::size_t bytes = 0;
    return checkedMultiply(count, itemBytes, bytes)
        && checkedAdd(total, bytes);
}

bool sameGrid(const SceneFluidGridCandidateSet& candidateSet,
              const PeriodicCartesianGrid& grid) {
    return candidateSet.cellCounts == grid.cellCounts()
        && candidateSet.lowerMeters == grid.lowerMeters()
        && candidateSet.upperMeters == grid.upperMeters();
}

void validateSettings(const SceneFluidGridCandidateSettings& settings) {
    if (!std::isfinite(settings.boundingPaddingMeters)
        || settings.boundingPaddingMeters < 0.0
        || !std::isfinite(settings.domainBoundaryToleranceMeters)
        || settings.domainBoundaryToleranceMeters < 0.0) {
        throw std::invalid_argument(
            "scene fluid grid-candidate settings are invalid");
    }
}

double boundedCoordinate(
    const double value,
    const double lower,
    const double upper,
    const double tolerance) {
    if (value < lower) {
        if (lower - value > tolerance) {
            throw std::invalid_argument(
                "scene fluid surface lies below the grid domain");
        }
        return lower;
    }
    if (value > upper) {
        if (value - upper > tolerance) {
            throw std::invalid_argument(
                "scene fluid surface lies above the grid domain");
        }
        return upper;
    }
    return value;
}

struct AxisCellRange {
    std::size_t first = 0;
    std::size_t last = 0;
};

AxisCellRange axisCellRange(
    const double minimum,
    const double maximum,
    const double domainLower,
    const double domainUpper,
    const double spacing,
    const std::size_t cellCount,
    const double padding) {
    const double paddedMinimum = std::max(
        domainLower, minimum - padding);
    const double paddedMaximum = std::min(
        domainUpper, maximum + padding);
    if (!std::isfinite(paddedMinimum) || !std::isfinite(paddedMaximum)
        || paddedMinimum > paddedMaximum) {
        throw std::invalid_argument(
            "scene fluid surface padding produced invalid bounds");
    }

    const double firstCoordinate = std::nextafter(
        (paddedMinimum - domainLower) / spacing,
        -std::numeric_limits<double>::infinity());
    const double lastCoordinate = std::nextafter(
        (paddedMaximum - domainLower) / spacing,
        std::numeric_limits<double>::infinity());
    const auto clampedIndex = [cellCount](const double coordinate) {
        if (coordinate <= 0.0) {
            return std::size_t{0};
        }
        if (coordinate >= static_cast<double>(cellCount)) {
            return cellCount - 1;
        }
        return static_cast<std::size_t>(std::floor(coordinate));
    };
    return {clampedIndex(firstCoordinate), clampedIndex(lastCoordinate)};
}

std::uint64_t candidateFingerprint(
    const SceneFluidGridCandidateSet& candidateSet) {
    Fingerprint fingerprint;
    fingerprint.integer(candidateSet.version);
    fingerprint.integer(candidateSet.surfaceDefinitionFingerprint);
    fingerprint.integer(candidateSet.surfaceStateFingerprint);
    fingerprint.integer(candidateSet.structureDefinitionFingerprint);
    fingerprint.integer(candidateSet.acceptedStepCount);
    fingerprint.real(candidateSet.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        candidateSet.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(
        candidateSet.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(
        candidateSet.cellCounts.z));
    fingerprint.real(candidateSet.lowerMeters.x);
    fingerprint.real(candidateSet.lowerMeters.y);
    fingerprint.real(candidateSet.lowerMeters.z);
    fingerprint.real(candidateSet.upperMeters.x);
    fingerprint.real(candidateSet.upperMeters.y);
    fingerprint.real(candidateSet.upperMeters.z);
    fingerprint.real(candidateSet.settings.boundingPaddingMeters);
    fingerprint.real(
        candidateSet.settings.domainBoundaryToleranceMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        candidateSet.triangleBounds.size()));
    for (const auto& bounds : candidateSet.triangleBounds) {
        fingerprint.integer(static_cast<std::uint64_t>(
            bounds.triangleIndex));
        fingerprint.integer(bounds.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(bounds.firstCell.i));
        fingerprint.integer(static_cast<std::uint64_t>(bounds.firstCell.j));
        fingerprint.integer(static_cast<std::uint64_t>(bounds.firstCell.k));
        fingerprint.integer(static_cast<std::uint64_t>(bounds.lastCell.i));
        fingerprint.integer(static_cast<std::uint64_t>(bounds.lastCell.j));
        fingerprint.integer(static_cast<std::uint64_t>(bounds.lastCell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            bounds.candidateCount));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        candidateSet.candidates.size()));
    for (const auto& candidate : candidateSet.candidates) {
        fingerprint.integer(static_cast<std::uint64_t>(candidate.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(candidate.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(candidate.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(candidate.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            candidate.triangleIndex));
        fingerprint.integer(candidate.triangleId);
    }
    return fingerprint.value();
}

auto candidateKey(const SceneFluidCellCandidate& candidate) {
    return std::tuple{candidate.cellIndex, candidate.triangleIndex};
}

} // namespace

std::span<const SceneFluidCellCandidate>
SceneFluidGridCandidateSet::candidatesForCell(
    const std::size_t cellIndex) const noexcept {
    const auto first = std::lower_bound(
        candidates.begin(), candidates.end(), cellIndex,
        [](const SceneFluidCellCandidate& candidate,
           const std::size_t expected) {
            return candidate.cellIndex < expected;
        });
    const auto last = std::upper_bound(
        first, candidates.end(), cellIndex,
        [](const std::size_t expected,
           const SceneFluidCellCandidate& candidate) {
            return expected < candidate.cellIndex;
        });
    return {first, last};
}

SceneFluidGridCandidateSet buildSceneFluidGridCandidates(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSettings& settings,
    const SceneFluidGridCandidateLimits& limits) {
    validateSceneFluidSurfaceState(surface, state);
    validateSettings(settings);

    std::size_t candidateBytes = 0;
    if (!checkedBytes(candidateBytes,
                      surface.triangles.size(),
                      sizeof(SceneFluidTriangleCandidateBounds))
        || candidateBytes > limits.maximumCandidateBytes) {
        throw std::length_error(
            "scene fluid triangle bounds exceed their byte limit");
    }

    SceneFluidGridCandidateSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    result.triangleBounds.reserve(surface.triangles.size());

    const Vector3 spacing = grid.cellSpacingMeters();
    std::size_t totalCandidates = 0;
    for (std::size_t triangleIndex = 0;
         triangleIndex < surface.triangles.size(); ++triangleIndex) {
        const auto& triangle = surface.triangles[triangleIndex];
        const double infinity = std::numeric_limits<double>::infinity();
        Vec3 minimum{infinity, infinity, infinity};
        Vec3 maximum{-infinity, -infinity, -infinity};
        for (const std::size_t vertexIndex : triangle.vertexIndices) {
            const Vec3& source = state.vertices[vertexIndex].positionMeters;
            const double x = boundedCoordinate(
                source.x, result.lowerMeters.x, result.upperMeters.x,
                settings.domainBoundaryToleranceMeters);
            const double y = boundedCoordinate(
                source.y, result.lowerMeters.y, result.upperMeters.y,
                settings.domainBoundaryToleranceMeters);
            const double z = boundedCoordinate(
                source.z, result.lowerMeters.z, result.upperMeters.z,
                settings.domainBoundaryToleranceMeters);
            minimum.x = std::min(minimum.x, x);
            minimum.y = std::min(minimum.y, y);
            minimum.z = std::min(minimum.z, z);
            maximum.x = std::max(maximum.x, x);
            maximum.y = std::max(maximum.y, y);
            maximum.z = std::max(maximum.z, z);
        }
        const AxisCellRange xRange = axisCellRange(
            minimum.x, maximum.x,
            result.lowerMeters.x, result.upperMeters.x,
            spacing.x, result.cellCounts.x,
            settings.boundingPaddingMeters);
        const AxisCellRange yRange = axisCellRange(
            minimum.y, maximum.y,
            result.lowerMeters.y, result.upperMeters.y,
            spacing.y, result.cellCounts.y,
            settings.boundingPaddingMeters);
        const AxisCellRange zRange = axisCellRange(
            minimum.z, maximum.z,
            result.lowerMeters.z, result.upperMeters.z,
            spacing.z, result.cellCounts.z,
            settings.boundingPaddingMeters);

        const std::size_t xCount = xRange.last - xRange.first + 1;
        const std::size_t yCount = yRange.last - yRange.first + 1;
        const std::size_t zCount = zRange.last - zRange.first + 1;
        std::size_t xyCount = 0;
        std::size_t triangleCandidateCount = 0;
        if (!checkedMultiply(xCount, yCount, xyCount)
            || !checkedMultiply(
                xyCount, zCount, triangleCandidateCount)
            || !checkedAdd(totalCandidates, triangleCandidateCount)
            || totalCandidates > limits.maximumCandidates) {
            throw std::length_error(
                "scene fluid grid candidates exceed their count limit");
        }
        result.triangleBounds.push_back({
            triangleIndex,
            triangle.id,
            {xRange.first, yRange.first, zRange.first},
            {xRange.last, yRange.last, zRange.last},
            triangleCandidateCount,
        });
    }

    if (!checkedBytes(candidateBytes,
                      totalCandidates,
                      sizeof(SceneFluidCellCandidate))
        || candidateBytes > limits.maximumCandidateBytes) {
        throw std::length_error(
            "scene fluid grid candidates exceed their byte limit");
    }

    result.candidates.reserve(totalCandidates);
    for (const auto& bounds : result.triangleBounds) {
        for (std::size_t k = bounds.firstCell.k;
             k <= bounds.lastCell.k; ++k) {
            for (std::size_t j = bounds.firstCell.j;
                 j <= bounds.lastCell.j; ++j) {
                for (std::size_t i = bounds.firstCell.i;
                     i <= bounds.lastCell.i; ++i) {
                    result.candidates.push_back({
                        grid.cellIndex(i, j, k),
                        {i, j, k},
                        bounds.triangleIndex,
                        bounds.triangleId,
                    });
                }
            }
        }
    }
    std::sort(result.candidates.begin(), result.candidates.end(),
              [](const auto& first, const auto& second) {
                  return candidateKey(first) < candidateKey(second);
              });
    result.fingerprint = candidateFingerprint(result);
    validateSceneFluidGridCandidates(result, surface, state, grid);
    return result;
}

void validateSceneFluidGridCandidates(
    const SceneFluidGridCandidateSet& candidateSet,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid) {
    validateSceneFluidSurfaceState(surface, state);
    validateSettings(candidateSet.settings);
    if (candidateSet.version != sceneFluidGridCandidateVersion
        || candidateSet.fingerprint == 0
        || candidateSet.surfaceDefinitionFingerprint != surface.fingerprint
        || candidateSet.surfaceStateFingerprint != state.fingerprint
        || candidateSet.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || candidateSet.acceptedStepCount != state.acceptedStepCount
        || candidateSet.simulationTimeSeconds
            != state.simulationTimeSeconds
        || !sameGrid(candidateSet, grid)
        || candidateSet.triangleBounds.size()
            != surface.triangles.size()) {
        throw std::invalid_argument(
            "scene fluid grid-candidate identity is invalid");
    }

    std::vector<std::size_t> observedCounts(
        surface.triangles.size(), 0);
    std::tuple<std::size_t, std::size_t> previousKey{};
    bool havePrevious = false;
    for (const auto& candidate : candidateSet.candidates) {
        if (candidate.cell.i >= candidateSet.cellCounts.x
            || candidate.cell.j >= candidateSet.cellCounts.y
            || candidate.cell.k >= candidateSet.cellCounts.z
            || candidate.cellIndex != grid.cellIndex(
                candidate.cell.i, candidate.cell.j, candidate.cell.k)
            || candidate.triangleIndex >= surface.triangles.size()
            || candidate.triangleId
                != surface.triangles[candidate.triangleIndex].id) {
            throw std::invalid_argument(
                "scene fluid grid candidate is invalid");
        }
        const auto& bounds =
            candidateSet.triangleBounds[candidate.triangleIndex];
        if (candidate.cell.i < bounds.firstCell.i
            || candidate.cell.i > bounds.lastCell.i
            || candidate.cell.j < bounds.firstCell.j
            || candidate.cell.j > bounds.lastCell.j
            || candidate.cell.k < bounds.firstCell.k
            || candidate.cell.k > bounds.lastCell.k) {
            throw std::invalid_argument(
                "scene fluid grid candidate lies outside its triangle bounds");
        }
        const auto key = candidateKey(candidate);
        if (havePrevious && !(previousKey < key)) {
            throw std::invalid_argument(
                "scene fluid grid candidates are not canonical and unique");
        }
        havePrevious = true;
        previousKey = key;
        ++observedCounts[candidate.triangleIndex];
    }
    for (std::size_t triangleIndex = 0;
         triangleIndex < candidateSet.triangleBounds.size();
         ++triangleIndex) {
        const auto& bounds = candidateSet.triangleBounds[triangleIndex];
        std::size_t xyCount = 0;
        std::size_t expectedCount = 0;
        const bool boundsOrdered =
            bounds.firstCell.i <= bounds.lastCell.i
            && bounds.firstCell.j <= bounds.lastCell.j
            && bounds.firstCell.k <= bounds.lastCell.k;
        const bool countValid = boundsOrdered
            && checkedMultiply(
                bounds.lastCell.i - bounds.firstCell.i + 1,
                bounds.lastCell.j - bounds.firstCell.j + 1,
                xyCount)
            && checkedMultiply(
                xyCount,
                bounds.lastCell.k - bounds.firstCell.k + 1,
                expectedCount);
        if (bounds.triangleIndex != triangleIndex
            || bounds.triangleId != surface.triangles[triangleIndex].id
            || !boundsOrdered
            || bounds.lastCell.i >= candidateSet.cellCounts.x
            || bounds.lastCell.j >= candidateSet.cellCounts.y
            || bounds.lastCell.k >= candidateSet.cellCounts.z
            || !countValid
            || bounds.candidateCount != expectedCount
            || bounds.candidateCount != observedCounts[triangleIndex]) {
            throw std::invalid_argument(
                "scene fluid triangle candidate bounds are invalid");
        }
    }
    if (candidateSet.fingerprint != candidateFingerprint(candidateSet)) {
        throw std::invalid_argument(
            "scene fluid grid-candidate fingerprint does not match its payload");
    }
}

} // namespace simwing::fsi::fluid
