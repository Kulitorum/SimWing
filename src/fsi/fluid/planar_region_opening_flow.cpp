#include "fluid/planar_region_opening_flow.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
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
    std::uint64_t value_ = fnvOffsetBasis;
};

class DisjointSet final {
public:
    explicit DisjointSet(const std::size_t count)
        : parents_(count), ranks_(count, 0) {
        for (std::size_t index = 0; index < count; ++index) {
            parents_[index] = index;
        }
    }

    std::size_t root(const std::size_t index) {
        if (parents_[index] != index) {
            parents_[index] = root(parents_[index]);
        }
        return parents_[index];
    }

    void join(const std::size_t first, const std::size_t second) {
        std::size_t a = root(first);
        std::size_t b = root(second);
        if (a == b) return;
        if (ranks_[a] < ranks_[b]) std::swap(a, b);
        parents_[b] = a;
        if (ranks_[a] == ranks_[b]) ++ranks_[a];
    }

private:
    std::vector<std::size_t> parents_;
    std::vector<std::uint8_t> ranks_;
};

struct ComponentWork {
    std::uint64_t stableId = 0;
    std::vector<std::size_t> regions;
    std::vector<std::size_t> openings;
};

void validateSettings(
    const PlanarPressureRegionOpeningFlowSettings& settings) {
    if (!std::isfinite(
            settings.absoluteFlowRateToleranceCubicMetersPerSecond)
        || settings.absoluteFlowRateToleranceCubicMetersPerSecond < 0.0
        || !std::isfinite(settings.relativeFlowRateTolerance)
        || settings.relativeFlowRateTolerance < 0.0
        || (settings.absoluteFlowRateToleranceCubicMetersPerSecond == 0.0
            && settings.relativeFlowRateTolerance == 0.0)
        || !std::isfinite(settings.absoluteVolumeToleranceCubicMeters)
        || settings.absoluteVolumeToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeVolumeTolerance)
        || settings.relativeVolumeTolerance < 0.0
        || (settings.absoluteVolumeToleranceCubicMeters == 0.0
            && settings.relativeVolumeTolerance == 0.0)
        || !std::isfinite(settings.relativeCholeskyPivotTolerance)
        || !(settings.relativeCholeskyPivotTolerance > 0.0)
        || settings.relativeCholeskyPivotTolerance >= 1.0) {
        throw std::invalid_argument(
            "planar pressure region opening-flow settings are invalid");
    }
}

void validateLimits(const PlanarPressureRegionOpeningFlowLimits& limits) {
    if (limits.maximumIntervals == 0 || limits.maximumRegions == 0
        || limits.maximumOpenings == 0 || limits.maximumOwnedBytes == 0
        || limits.maximumFactorizationBytes == 0
        || limits.maximumFactorizationWork == 0) {
        throw std::invalid_argument(
            "planar pressure region opening-flow limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t count,
                           const std::size_t itemBytes,
                           const char* message) {
    if (count != 0
        && itemBytes > std::numeric_limits<std::size_t>::max() / count) {
        throw std::length_error(message);
    }
    return count * itemBytes;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second,
                       const char* message) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(message);
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t openingCount,
                              const std::size_t regionCount,
                              const std::size_t componentCount) {
    constexpr const char* message =
        "planar pressure region opening-flow storage size overflows";
    const std::size_t openings = checkedProduct(
        openingCount, sizeof(PlanarPressureRegionOpeningFlow), message);
    const std::size_t regions = checkedProduct(
        regionCount, sizeof(PlanarPressureRegionOpeningBalance), message);
    const std::size_t components = checkedProduct(
        componentCount, sizeof(PlanarPressureRegionOpeningComponent),
        message);
    return checkedSum(checkedSum(openings, regions, message),
                      components, message);
}

std::uint64_t checkedCube(const std::size_t value) {
    const std::uint64_t sample = static_cast<std::uint64_t>(value);
    if (sample != 0
        && sample > std::numeric_limits<std::uint64_t>::max() / sample) {
        throw std::length_error(
            "planar pressure region opening-flow work size overflows");
    }
    const std::uint64_t square = sample * sample;
    if (sample != 0
        && square > std::numeric_limits<std::uint64_t>::max() / sample) {
        throw std::length_error(
            "planar pressure region opening-flow work size overflows");
    }
    return square * sample;
}

double flowTolerance(
    const double reference,
    const PlanarPressureRegionOpeningFlowSettings& settings) {
    const double tolerance = std::max(
        settings.absoluteFlowRateToleranceCubicMetersPerSecond,
        settings.relativeFlowRateTolerance * std::abs(reference));
    if (!std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "planar pressure region opening-flow tolerance is non-finite");
    }
    return tolerance;
}

double volumeTolerance(
    const PlanarPressureRegionOpeningBalance& region,
    const PlanarPressureRegionOpeningFlowSettings& settings) {
    const double reference = std::max({
        std::abs(region.previousVolumeCubicMeters),
        std::abs(region.currentVolumeCubicMeters),
        std::abs(region.geometryVolumeChangeCubicMeters),
        std::abs(region.integratedOutwardRelativeVolumeCubicMeters),
    });
    const double tolerance = std::max(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance * reference);
    if (!std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "planar pressure region opening-flow volume tolerance is non-finite");
    }
    return tolerance;
}

std::vector<double> solveReducedLaplacian(
    std::vector<double> matrix,
    const std::vector<double>& rightHandSide,
    const double relativePivotTolerance) {
    const std::size_t count = rightHandSide.size();
    if (matrix.size() != count * count) {
        throw std::invalid_argument(
            "planar pressure region opening-flow matrix size is invalid");
    }
    double maximumDiagonal = 0.0;
    for (std::size_t row = 0; row < count; ++row) {
        maximumDiagonal = std::max(
            maximumDiagonal, matrix[row * count + row]);
    }
    if (!std::isfinite(maximumDiagonal) || !(maximumDiagonal > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region opening-flow matrix is singular");
    }
    const double pivotTolerance =
        relativePivotTolerance * maximumDiagonal;
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            double value = matrix[row * count + column];
            for (std::size_t inner = 0; inner < column; ++inner) {
                value -= matrix[row * count + inner]
                    * matrix[column * count + inner];
            }
            if (row == column) {
                if (!std::isfinite(value) || value <= pivotTolerance) {
                    throw std::invalid_argument(
                        "planar pressure region opening-flow matrix is ill-conditioned");
                }
                matrix[row * count + column] = std::sqrt(value);
            } else {
                matrix[row * count + column] = value
                    / matrix[column * count + column];
            }
        }
    }

    std::vector<double> solution(rightHandSide);
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column < row; ++column) {
            solution[row] -= matrix[row * count + column]
                * solution[column];
        }
        solution[row] /= matrix[row * count + row];
    }
    for (std::size_t row = count; row-- > 0;) {
        for (std::size_t column = row + 1; column < count; ++column) {
            solution[row] -= matrix[column * count + row]
                * solution[column];
        }
        solution[row] /= matrix[row * count + row];
    }
    if (!std::ranges::all_of(
            solution,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "planar pressure region opening-flow solution is non-finite");
    }
    return solution;
}

std::uint64_t allocationFingerprint(
    const PlanarPressureRegionOpeningFlowAllocation& allocation) {
    Fingerprint fingerprint;
    fingerprint.integer(allocation.version);
    fingerprint.integer(allocation.sourceSweepVersion);
    fingerprint.enumeration(allocation.axis);
    fingerprint.real(allocation.durationSeconds);
    fingerprint.real(
        allocation.settings
            .absoluteFlowRateToleranceCubicMetersPerSecond);
    fingerprint.real(allocation.settings.relativeFlowRateTolerance);
    fingerprint.real(
        allocation.settings.absoluteVolumeToleranceCubicMeters);
    fingerprint.real(allocation.settings.relativeVolumeTolerance);
    fingerprint.real(allocation.settings.relativeCholeskyPivotTolerance);
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.openings.size()));
    for (const auto& opening : allocation.openings) {
        fingerprint.integer(opening.openingStableId);
        fingerprint.integer(opening.negativeSideRegionStableId);
        fingerprint.integer(opening.positiveSideRegionStableId);
        fingerprint.integer(opening.componentStableId);
        fingerprint.real(opening.areaSquareMeters);
        fingerprint.real(
            opening.relativeVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(opening.relativeNormalVelocityMetersPerSecond);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.regions.size()));
    for (const auto& region : allocation.regions) {
        fingerprint.integer(region.regionStableId);
        fingerprint.integer(region.componentStableId);
        fingerprint.real(region.previousVolumeCubicMeters);
        fingerprint.real(region.currentVolumeCubicMeters);
        fingerprint.real(region.geometryVolumeChangeCubicMeters);
        fingerprint.real(
            region.requestedOutwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            region.solvedOutwardRelativeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            region.integratedOutwardRelativeVolumeCubicMeters);
        fingerprint.real(region.continuityResidualCubicMeters);
        fingerprint.real(region.continuityToleranceCubicMeters);
        fingerprint.integer(static_cast<std::uint8_t>(
            region.withinTolerance));
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.components.size()));
    for (const auto& component : allocation.components) {
        fingerprint.integer(component.componentStableId);
        fingerprint.integer(component.correctionRegionStableId);
        fingerprint.integer(
            static_cast<std::uint64_t>(component.regionCount));
        fingerprint.integer(
            static_cast<std::uint64_t>(component.openingCount));
        fingerprint.real(
            component.requestedOutwardFlowRateSumCubicMetersPerSecond);
        fingerprint.real(
            component.appliedCorrectionFlowRateCubicMetersPerSecond);
        fingerprint.real(
            component.compatibilityToleranceCubicMetersPerSecond);
        fingerprint.real(
            component.maximumAbsoluteFlowBalanceResidualCubicMetersPerSecond);
        fingerprint.integer(static_cast<std::uint8_t>(
            component.sourceCompatible));
        fingerprint.integer(static_cast<std::uint8_t>(
            component.linearSolveWithinTolerance));
        fingerprint.integer(static_cast<std::uint8_t>(
            component.feasible));
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.failedComponentCount));
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.failedRegionCount));
    fingerprint.real(
        allocation.maximumAbsoluteOpeningNormalVelocityMetersPerSecond);
    fingerprint.real(
        allocation.maximumAbsoluteContinuityResidualCubicMeters);
    fingerprint.real(allocation.globalGeometryVolumeChangeCubicMeters);
    fingerprint.real(
        allocation.globalIntegratedOutwardRelativeVolumeCubicMeters);
    fingerprint.real(allocation.globalContinuityResidualCubicMeters);
    fingerprint.integer(static_cast<std::uint8_t>(
        allocation.allComponentsFeasible));
    fingerprint.integer(static_cast<std::uint8_t>(
        allocation.allRegionsWithinTolerance));
    fingerprint.integer(
        static_cast<std::uint64_t>(allocation.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionOpeningFlowAllocation buildAllocation(
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition> definitions,
    const PlanarPressureRegionOpeningFlowSettings& settings,
    const PlanarPressureRegionOpeningFlowLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    validatePlanarPressureRegionSweepLedger(
        sweep,
        {limits.maximumIntervals,
         limits.maximumRegions,
         std::numeric_limits<std::size_t>::max()});
    if (definitions.size() > limits.maximumOpenings) {
        throw std::length_error(
            "planar pressure region opening-flow exceeds its opening limit");
    }

    std::vector<PlanarPressureRegionOpeningDefinition> canonical(
        definitions.begin(), definitions.end());
    std::ranges::sort(
        canonical, {},
        &PlanarPressureRegionOpeningDefinition::openingStableId);
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto& opening = canonical[index];
        if (opening.openingStableId == 0
            || opening.negativeSideRegionStableId == 0
            || opening.positiveSideRegionStableId == 0
            || opening.negativeSideRegionStableId
                == opening.positiveSideRegionStableId
            || !std::isfinite(opening.areaSquareMeters)
            || !(opening.areaSquareMeters > 0.0)
            || (index != 0
                && opening.openingStableId
                    == canonical[index - 1].openingStableId)) {
            throw std::invalid_argument(
                "planar pressure region opening definition is invalid");
        }
    }

    PlanarPressureRegionOpeningFlowAllocation result;
    result.sourceSweepVersion = sweep.version;
    result.axis = sweep.axis;
    result.durationSeconds = sweep.durationSeconds;
    result.settings = settings;
    result.regions.reserve(sweep.regions.size());
    std::map<std::uint64_t, std::size_t> regionIndexes;
    for (std::size_t index = 0; index < sweep.regions.size(); ++index) {
        const auto& source = sweep.regions[index];
        PlanarPressureRegionOpeningBalance region;
        region.regionStableId = source.regionStableId;
        region.previousVolumeCubicMeters =
            source.previousVolumeCubicMeters;
        region.currentVolumeCubicMeters =
            source.currentVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters =
            source.geometryVolumeChangeCubicMeters;
        region.requestedOutwardRelativeFlowRateCubicMetersPerSecond =
            -region.geometryVolumeChangeCubicMeters
            / result.durationSeconds;
        if (!std::isfinite(
                region.requestedOutwardRelativeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "planar pressure region opening-flow request is non-finite");
        }
        regionIndexes.emplace(region.regionStableId, index);
        result.regions.push_back(region);
    }

    DisjointSet disjoint(result.regions.size());
    result.openings.reserve(canonical.size());
    std::vector<std::pair<std::size_t, std::size_t>> openingEndpoints;
    openingEndpoints.reserve(canonical.size());
    for (const auto& source : canonical) {
        const auto negative = regionIndexes.find(
            source.negativeSideRegionStableId);
        const auto positive = regionIndexes.find(
            source.positiveSideRegionStableId);
        if (negative == regionIndexes.end()
            || positive == regionIndexes.end()) {
            throw std::invalid_argument(
                "planar pressure region opening references a foreign region");
        }
        disjoint.join(negative->second, positive->second);
        openingEndpoints.emplace_back(
            negative->second, positive->second);
        result.openings.push_back({
            source.openingStableId,
            source.negativeSideRegionStableId,
            source.positiveSideRegionStableId,
            0,
            source.areaSquareMeters,
            0.0,
            0.0,
        });
    }

    std::map<std::size_t, ComponentWork> byRoot;
    for (std::size_t region = 0; region < result.regions.size(); ++region) {
        byRoot[disjoint.root(region)].regions.push_back(region);
    }
    for (std::size_t opening = 0; opening < result.openings.size();
         ++opening) {
        byRoot[disjoint.root(openingEndpoints[opening].first)]
            .openings.push_back(opening);
    }
    std::vector<ComponentWork> work;
    work.reserve(byRoot.size());
    for (auto& entry : byRoot) {
        auto component = std::move(entry.second);
        component.stableId = result.regions[
            component.regions.front()].regionStableId;
        work.push_back(std::move(component));
    }
    std::ranges::sort(work, {}, &ComponentWork::stableId);

    result.components.reserve(work.size());
    std::vector<std::size_t> localByRegion(
        result.regions.size(), std::numeric_limits<std::size_t>::max());
    std::uint64_t factorizationWork = 0;
    for (const auto& source : work) {
        PlanarPressureRegionOpeningComponent component;
        component.componentStableId = source.stableId;
        component.correctionRegionStableId = source.stableId;
        component.regionCount = source.regions.size();
        component.openingCount = source.openings.size();
        double absoluteRequestSum = 0.0;
        for (const std::size_t regionIndex : source.regions) {
            auto& region = result.regions[regionIndex];
            region.componentStableId = component.componentStableId;
            const double request =
                region.requestedOutwardRelativeFlowRateCubicMetersPerSecond;
            component.requestedOutwardFlowRateSumCubicMetersPerSecond +=
                request;
            absoluteRequestSum += std::abs(request);
        }
        component.compatibilityToleranceCubicMetersPerSecond =
            flowTolerance(absoluteRequestSum, settings);
        component.sourceCompatible = std::abs(
            component.requestedOutwardFlowRateSumCubicMetersPerSecond)
            <= component.compatibilityToleranceCubicMetersPerSecond;
        component.appliedCorrectionFlowRateCubicMetersPerSecond =
            component.sourceCompatible
            ? -component.requestedOutwardFlowRateSumCubicMetersPerSecond
            : 0.0;
        for (const std::size_t openingIndex : source.openings) {
            result.openings[openingIndex].componentStableId =
                component.componentStableId;
        }

        std::vector<double> adjustedRequests(source.regions.size(), 0.0);
        for (std::size_t local = 0; local < source.regions.size(); ++local) {
            const std::size_t regionIndex = source.regions[local];
            localByRegion[regionIndex] = local;
            adjustedRequests[local] = result.regions[regionIndex]
                .requestedOutwardRelativeFlowRateCubicMetersPerSecond;
        }
        if (component.sourceCompatible) {
            adjustedRequests.front() +=
                component.appliedCorrectionFlowRateCubicMetersPerSecond;
            std::vector<double> potentials(source.regions.size(), 0.0);
            if (source.regions.size() > 1) {
                const std::size_t dimension = source.regions.size() - 1;
                const std::size_t entryCount = checkedProduct(
                    dimension, dimension,
                    "planar pressure region opening-flow matrix size overflows");
                const std::size_t factorizationBytes = checkedProduct(
                    entryCount, sizeof(double),
                    "planar pressure region opening-flow matrix size overflows");
                if (factorizationBytes > limits.maximumFactorizationBytes) {
                    throw std::length_error(
                        "planar pressure region opening-flow exceeds its factorization byte limit");
                }
                const std::uint64_t nextWork = checkedCube(dimension);
                if (nextWork > limits.maximumFactorizationWork
                    - factorizationWork) {
                    throw std::length_error(
                        "planar pressure region opening-flow exceeds its factorization work limit");
                }
                factorizationWork += nextWork;
                std::vector<double> matrix(entryCount, 0.0);
                for (const std::size_t openingIndex : source.openings) {
                    const auto [negativeIndex, positiveIndex] =
                        openingEndpoints[openingIndex];
                    const std::size_t negative =
                        localByRegion[negativeIndex];
                    const std::size_t positive =
                        localByRegion[positiveIndex];
                    const double area =
                        result.openings[openingIndex].areaSquareMeters;
                    if (negative != 0) {
                        matrix[(negative - 1) * dimension
                               + (negative - 1)] += area;
                    }
                    if (positive != 0) {
                        matrix[(positive - 1) * dimension
                               + (positive - 1)] += area;
                    }
                    if (negative != 0 && positive != 0) {
                        matrix[(negative - 1) * dimension
                               + (positive - 1)] -= area;
                        matrix[(positive - 1) * dimension
                               + (negative - 1)] -= area;
                    }
                }
                std::vector<double> rightHandSide(
                    adjustedRequests.begin() + 1,
                    adjustedRequests.end());
                const auto reduced = solveReducedLaplacian(
                    std::move(matrix), rightHandSide,
                    settings.relativeCholeskyPivotTolerance);
                std::copy(reduced.begin(), reduced.end(),
                          potentials.begin() + 1);
            }
            for (const std::size_t openingIndex : source.openings) {
                const auto [negativeIndex, positiveIndex] =
                    openingEndpoints[openingIndex];
                const double flow =
                    result.openings[openingIndex].areaSquareMeters
                    * (potentials[localByRegion[negativeIndex]]
                       - potentials[localByRegion[positiveIndex]]);
                result.openings[openingIndex]
                    .relativeVolumeFlowRateCubicMetersPerSecond = flow;
                result.openings[openingIndex]
                    .relativeNormalVelocityMetersPerSecond = flow
                    / result.openings[openingIndex].areaSquareMeters;
            }
        }

        for (const std::size_t openingIndex : source.openings) {
            const auto [negativeIndex, positiveIndex] =
                openingEndpoints[openingIndex];
            const double flow = result.openings[openingIndex]
                .relativeVolumeFlowRateCubicMetersPerSecond;
            result.regions[negativeIndex]
                .solvedOutwardRelativeFlowRateCubicMetersPerSecond += flow;
            result.regions[positiveIndex]
                .solvedOutwardRelativeFlowRateCubicMetersPerSecond -= flow;
            result.maximumAbsoluteOpeningNormalVelocityMetersPerSecond =
                std::max(
                    result.maximumAbsoluteOpeningNormalVelocityMetersPerSecond,
                    std::abs(result.openings[openingIndex]
                        .relativeNormalVelocityMetersPerSecond));
        }
        component.linearSolveWithinTolerance = component.sourceCompatible;
        for (std::size_t local = 0; local < source.regions.size(); ++local) {
            const auto& region = result.regions[source.regions[local]];
            const double residual =
                region.solvedOutwardRelativeFlowRateCubicMetersPerSecond
                - adjustedRequests[local];
            component.maximumAbsoluteFlowBalanceResidualCubicMetersPerSecond =
                std::max(
                    component
                        .maximumAbsoluteFlowBalanceResidualCubicMetersPerSecond,
                    std::abs(residual));
            const double tolerance = flowTolerance(
                std::max(
                    std::abs(adjustedRequests[local]),
                    std::abs(region
                        .solvedOutwardRelativeFlowRateCubicMetersPerSecond)),
                settings);
            component.linearSolveWithinTolerance =
                component.linearSolveWithinTolerance
                && std::abs(residual) <= tolerance;
        }
        component.feasible = component.sourceCompatible
            && component.linearSolveWithinTolerance;
        result.failedComponentCount += !component.feasible;
        result.components.push_back(component);
        for (const std::size_t regionIndex : source.regions) {
            localByRegion[regionIndex] =
                std::numeric_limits<std::size_t>::max();
        }
    }

    for (auto& region : result.regions) {
        region.integratedOutwardRelativeVolumeCubicMeters =
            result.durationSeconds
            * region.solvedOutwardRelativeFlowRateCubicMetersPerSecond;
        region.continuityResidualCubicMeters =
            region.geometryVolumeChangeCubicMeters
            + region.integratedOutwardRelativeVolumeCubicMeters;
        region.continuityToleranceCubicMeters = volumeTolerance(
            region, settings);
        region.withinTolerance =
            std::isfinite(region.continuityResidualCubicMeters)
            && std::abs(region.continuityResidualCubicMeters)
                <= region.continuityToleranceCubicMeters;
        result.failedRegionCount += !region.withinTolerance;
        result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
            result.maximumAbsoluteContinuityResidualCubicMeters,
            std::abs(region.continuityResidualCubicMeters));
        result.globalGeometryVolumeChangeCubicMeters +=
            region.geometryVolumeChangeCubicMeters;
        result.globalIntegratedOutwardRelativeVolumeCubicMeters +=
            region.integratedOutwardRelativeVolumeCubicMeters;
    }
    result.globalContinuityResidualCubicMeters =
        result.globalGeometryVolumeChangeCubicMeters
        + result.globalIntegratedOutwardRelativeVolumeCubicMeters;
    result.maximumAbsoluteContinuityResidualCubicMeters = std::max(
        result.maximumAbsoluteContinuityResidualCubicMeters,
        std::abs(result.globalContinuityResidualCubicMeters));
    result.allComponentsFeasible = result.failedComponentCount == 0;
    result.allRegionsWithinTolerance = result.failedRegionCount == 0;
    result.ownedStorageBytes = ownedStorageBytes(
        result.openings.size(), result.regions.size(),
        result.components.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar pressure region opening-flow exceeds its owned byte limit");
    }
    if (!std::isfinite(
            result.maximumAbsoluteOpeningNormalVelocityMetersPerSecond)
        || !std::isfinite(
            result.maximumAbsoluteContinuityResidualCubicMeters)
        || !std::isfinite(result.globalGeometryVolumeChangeCubicMeters)
        || !std::isfinite(
            result.globalIntegratedOutwardRelativeVolumeCubicMeters)
        || !std::isfinite(result.globalContinuityResidualCubicMeters)) {
        throw std::invalid_argument(
            "planar pressure region opening-flow aggregate is non-finite");
    }
    result.fingerprint = allocationFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionOpeningFlowAllocation
solvePlanarPressureRegionOpeningFlow(
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition> openings,
    const PlanarPressureRegionOpeningFlowSettings& settings,
    const PlanarPressureRegionOpeningFlowLimits& limits) {
    auto result = buildAllocation(sweep, openings, settings, limits);
    validatePlanarPressureRegionOpeningFlow(
        result, sweep, openings, limits);
    return result;
}

void validatePlanarPressureRegionOpeningFlow(
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionSweepLedger& sweep,
    const std::span<const PlanarPressureRegionOpeningDefinition> openings,
    const PlanarPressureRegionOpeningFlowLimits& limits) {
    validateLimits(limits);
    if (allocation.openings.size() > limits.maximumOpenings
        || allocation.regions.size() > limits.maximumRegions
        || allocation.components.size() > limits.maximumRegions) {
        throw std::length_error(
            "planar pressure region opening-flow validation limit exceeded");
    }
    const auto expected = buildAllocation(
        sweep, openings, allocation.settings, limits);
    if (expected != allocation) {
        throw std::invalid_argument(
            "planar pressure region opening-flow allocation is invalid");
    }
}

} // namespace simwing::fsi::fluid
