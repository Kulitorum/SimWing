#include "scene_fluid_opening_flux.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
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

struct AxisInterpolation {
    std::size_t first = 0;
    std::size_t second = 0;
    double fraction = 0.0;
};

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

bool finite(const Vec3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
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

double dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

double length(const Vec3& value) {
    return std::hypot(value.x, value.y, value.z);
}

bool closeValue(const double first, const double second) {
    return std::isfinite(first) && std::isfinite(second)
        && std::abs(first - second)
            <= 1.0e-12
                + 1.0e-10 * std::max(std::abs(first), std::abs(second));
}

AxisInterpolation axisInterpolation(const double position,
                                    const double lower,
                                    const double spacing,
                                    const std::size_t count,
                                    const double latticeOffset) {
    const double lattice = (position - lower) / spacing - latticeOffset;
    const double base = std::floor(lattice);
    const double fraction = lattice - base;
    std::size_t first = 0;
    if (base < 0.0) {
        first = count - 1;
    } else if (base >= static_cast<double>(count)) {
        first = 0;
    } else {
        first = static_cast<std::size_t>(base);
    }
    const std::size_t second = first + 1 == count ? 0 : first + 1;
    return {first, second, fraction};
}

double interpolateComponent(
    const fluid::PeriodicCartesianGrid& grid,
    const std::span<const double> values,
    const Vec3& position,
    const std::array<double, 3>& offsets) {
    const auto counts = grid.cellCounts();
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    const AxisInterpolation x = axisInterpolation(
        position.x, lower.x, spacing.x, counts.x, offsets[0]);
    const AxisInterpolation y = axisInterpolation(
        position.y, lower.y, spacing.y, counts.y, offsets[1]);
    const AxisInterpolation z = axisInterpolation(
        position.z, lower.z, spacing.z, counts.z, offsets[2]);
    const std::array<std::size_t, 2> xs{x.first, x.second};
    const std::array<std::size_t, 2> ys{y.first, y.second};
    const std::array<std::size_t, 2> zs{z.first, z.second};
    const std::array<double, 2> xWeights{1.0 - x.fraction, x.fraction};
    const std::array<double, 2> yWeights{1.0 - y.fraction, y.fraction};
    const std::array<double, 2> zWeights{1.0 - z.fraction, z.fraction};
    double result = 0.0;
    for (std::size_t dz = 0; dz < 2; ++dz) {
        for (std::size_t dy = 0; dy < 2; ++dy) {
            for (std::size_t dx = 0; dx < 2; ++dx) {
                result += xWeights[dx] * yWeights[dy] * zWeights[dz]
                    * values[grid.cellIndex(xs[dx], ys[dy], zs[dz])];
            }
        }
    }
    return result;
}

Vec3 interpolateMacVelocity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    const Vec3& position) {
    return {
        interpolateComponent(
            grid, velocity.xFaces(), position, {0.0, 0.5, 0.5}),
        interpolateComponent(
            grid, velocity.yFaces(), position, {0.5, 0.0, 0.5}),
        interpolateComponent(
            grid, velocity.zFaces(), position, {0.5, 0.5, 0.0}),
    };
}

Vec3 barycentricPoint(const Vec3& first,
                      const Vec3& second,
                      const Vec3& third,
                      const std::array<double, 3>& barycentric) {
    return {
        barycentric[0] * first.x + barycentric[1] * second.x
            + barycentric[2] * third.x,
        barycentric[0] * first.y + barycentric[1] * second.y
            + barycentric[2] * third.y,
        barycentric[0] * first.z + barycentric[1] * second.z
            + barycentric[2] * third.z,
    };
}

double faceNormalVelocity(
    const SceneFluidOpeningGridPatch& patch,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    const std::size_t index = grid.cellIndex(
        patch.faceI, patch.faceJ, patch.faceK);
    switch (patch.faceAxis) {
    case SceneFluidOpeningPatchFaceAxis::X:
        return velocity.xFaces()[index]
            * patch.unitNormalNegativeToPositive.x;
    case SceneFluidOpeningPatchFaceAxis::Y:
        return velocity.yFaces()[index]
            * patch.unitNormalNegativeToPositive.y;
    case SceneFluidOpeningPatchFaceAxis::Z:
        return velocity.zFaces()[index]
            * patch.unitNormalNegativeToPositive.z;
    }
    throw std::invalid_argument(
        "scene fluid opening flux has an invalid face axis");
}

double cellPatchFluidFlow(
    const SceneFluidOpeningGridPatch& patch,
    const std::span<const fluid::SceneFluidClippedVertex> vertices,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    std::size_t& evaluationCount,
    const SceneFluidOpeningFluxLimits& limits) {
    constexpr std::array<std::array<double, 3>, 4> barycentrics{{
        {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}},
        {{0.6, 0.2, 0.2}},
        {{0.2, 0.6, 0.2}},
        {{0.2, 0.2, 0.6}},
    }};
    constexpr std::array<double, 4> weights{{
        -27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0,
    }};
    double area = 0.0;
    double flow = 0.0;
    for (std::size_t index = 1; index + 1 < vertices.size(); ++index) {
        const Vec3& first = vertices[0].positionMeters;
        const Vec3& second = vertices[index].positionMeters;
        const Vec3& third = vertices[index + 1].positionMeters;
        const double triangleArea = 0.5 * length(cross(
            subtract(second, first), subtract(third, first)));
        if (!(triangleArea > 0.0)) {
            continue;
        }
        std::size_t nextEvaluationCount = 0;
        if (!checkedAdd(evaluationCount, barycentrics.size(),
                        nextEvaluationCount)
            || nextEvaluationCount > limits.maximumVelocityEvaluations) {
            throw std::length_error(
                "scene fluid opening flux exceeds its velocity-evaluation limit");
        }
        evaluationCount = nextEvaluationCount;
        double meanNormalVelocity = 0.0;
        for (std::size_t sample = 0; sample < barycentrics.size(); ++sample) {
            const Vec3 position = barycentricPoint(
                first, second, third, barycentrics[sample]);
            const Vec3 fluidVelocity = interpolateMacVelocity(
                grid, velocity, position);
            if (!finite(fluidVelocity)) {
                throw std::overflow_error(
                    "scene fluid opening velocity interpolation is not finite");
            }
            meanNormalVelocity += weights[sample]
                * dot(fluidVelocity, patch.unitNormalNegativeToPositive);
        }
        flow += triangleArea * meanNormalVelocity;
        area += triangleArea;
    }
    if (!closeValue(area, patch.areaSquareMeters)
        || !std::isfinite(flow)) {
        throw std::invalid_argument(
            "scene fluid opening flux polygon area or integral is invalid");
    }
    return flow;
}

std::uint64_t velocityFingerprint(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    Fingerprint fingerprint;
    const auto counts = grid.cellCounts();
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    fingerprint.integer(static_cast<std::uint64_t>(counts.x));
    fingerprint.integer(static_cast<std::uint64_t>(counts.y));
    fingerprint.integer(static_cast<std::uint64_t>(counts.z));
    for (const double value : {
             lower.x, lower.y, lower.z, upper.x, upper.y, upper.z}) {
        fingerprint.real(value);
    }
    for (const auto component : {
             velocity.xFaces(), velocity.yFaces(), velocity.zFaces()}) {
        fingerprint.integer(static_cast<std::uint64_t>(component.size()));
        for (const double value : component) {
            fingerprint.real(value);
        }
    }
    return fingerprint.value();
}

std::size_t storageBytesForCounts(const std::size_t regionCount,
                                  const std::size_t openingCount,
                                  const std::size_t sampleCount) {
    std::size_t regionBytes = 0;
    std::size_t openingBytes = 0;
    std::size_t sampleBytes = 0;
    std::size_t firstTotal = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            regionCount, sizeof(SceneFluidOpeningRegionFlux), regionBytes)
        || !checkedMultiply(
            openingCount, sizeof(SceneFluidOpeningFlux), openingBytes)
        || !checkedMultiply(
            sampleCount, sizeof(SceneFluidOpeningFluxSample), sampleBytes)
        || !checkedAdd(regionBytes, openingBytes, firstTotal)
        || !checkedAdd(firstTotal, sampleBytes, total)) {
        throw std::length_error(
            "scene fluid opening-flux storage size overflows");
    }
    return total;
}

std::size_t storageBytes(const SceneFluidOpeningFluxSet& flux) {
    return storageBytesForCounts(
        flux.regions.size(), flux.openings.size(), flux.samples.size());
}

std::uint64_t fluxFingerprint(const SceneFluidOpeningFluxSet& flux) {
    Fingerprint fingerprint;
    fingerprint.integer(flux.version);
    fingerprint.integer(flux.surfaceDefinitionFingerprint);
    fingerprint.integer(flux.surfaceStateFingerprint);
    fingerprint.integer(flux.openingPatchFingerprint);
    fingerprint.integer(flux.velocityFingerprint);
    fingerprint.integer(flux.structureDefinitionFingerprint);
    fingerprint.integer(flux.acceptedStepCount);
    fingerprint.real(flux.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(flux.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(flux.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(flux.cellCounts.z));
    for (const double value : {
             flux.lowerMeters.x, flux.lowerMeters.y, flux.lowerMeters.z,
             flux.upperMeters.x, flux.upperMeters.y, flux.upperMeters.z,
             flux.totalAreaSquareMeters,
             flux.totalFluidVolumeFlowRateCubicMetersPerSecond,
             flux.totalSurfaceSweepRateCubicMetersPerSecond,
             flux.totalRelativeVolumeFlowRateCubicMetersPerSecond,
             flux.globalFluidRegionBalanceResidualCubicMetersPerSecond,
             flux.globalSurfaceRegionBalanceResidualCubicMetersPerSecond,
             flux.globalRelativeRegionBalanceResidualCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        flux.velocityEvaluationCount));
    fingerprint.integer(static_cast<std::uint64_t>(flux.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(flux.regions.size()));
    for (const auto& region : flux.regions) {
        fingerprint.integer(static_cast<std::uint64_t>(region.regionIndex));
        fingerprint.integer(region.regionId);
        fingerprint.enumeration(region.kind);
        fingerprint.real(
            region.outwardFluidVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            region.outwardSurfaceSweepRateCubicMetersPerSecond);
        fingerprint.real(
            region.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(flux.openings.size()));
    for (const auto& opening : flux.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(
            opening.openingIndex));
        fingerprint.integer(opening.openingId);
        fingerprint.integer(opening.negativeSideRegionId);
        fingerprint.integer(opening.positiveSideRegionId);
        fingerprint.enumeration(opening.role);
        fingerprint.integer(static_cast<std::uint64_t>(opening.firstSample));
        fingerprint.integer(static_cast<std::uint64_t>(opening.sampleCount));
        fingerprint.real(opening.areaSquareMeters);
        fingerprint.real(
            opening.fluidVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            opening.surfaceSweepRateCubicMetersPerSecond);
        fingerprint.real(
            opening.relativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(flux.samples.size()));
    for (const auto& sample : flux.samples) {
        fingerprint.integer(sample.patchStableId);
        fingerprint.integer(sample.sourcePointStableId);
        fingerprint.integer(static_cast<std::uint64_t>(sample.openingIndex));
        fingerprint.integer(sample.openingId);
        fingerprint.integer(sample.negativeSideRegionId);
        fingerprint.integer(sample.positiveSideRegionId);
        fingerprint.enumeration(sample.role);
        for (const double value : {
                 sample.areaSquareMeters,
                 sample.fluidNormalVelocityMetersPerSecond,
                 sample.surfaceNormalVelocityMetersPerSecond,
                 sample.relativeNormalVelocityMetersPerSecond,
                 sample.fluidVolumeFlowRateCubicMetersPerSecond,
                 sample.surfaceSweepRateCubicMetersPerSecond,
                 sample.relativeVolumeFlowRateCubicMetersPerSecond}) {
            fingerprint.real(value);
        }
    }
    return fingerprint.value();
}

SceneFluidOpeningFluxSet buildFlux(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    const SceneFluidOpeningFluxLimits& limits) {
    if (!velocity.matches(grid) || !fluid::isFinite(velocity)) {
        throw std::invalid_argument(
            "scene fluid opening flux requires a matching finite MAC field");
    }
    if (surface.regions.size() > limits.maximumRegions
        || quadrature.openings.size() > limits.maximumOpenings
        || patches.patches.size() > limits.maximumPatchSamples) {
        throw std::length_error(
            "scene fluid opening flux exceeds its count limit");
    }
    const std::size_t expectedBytes = storageBytesForCounts(
        surface.regions.size(), quadrature.openings.size(),
        patches.patches.size());
    if (expectedBytes > limits.maximumFluxBytes) {
        throw std::length_error(
            "scene fluid opening flux exceeds its byte limit");
    }

    SceneFluidOpeningFluxSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = patches.surfaceStateFingerprint;
    result.openingPatchFingerprint = patches.fingerprint;
    result.velocityFingerprint = velocityFingerprint(grid, velocity);
    result.structureDefinitionFingerprint =
        patches.structureDefinitionFingerprint;
    result.acceptedStepCount = patches.acceptedStepCount;
    result.simulationTimeSeconds = patches.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.regions.reserve(surface.regions.size());
    for (std::size_t regionIndex = 0;
         regionIndex < surface.regions.size(); ++regionIndex) {
        result.regions.push_back({
            regionIndex,
            surface.regions[regionIndex].id,
            surface.regions[regionIndex].kind,
        });
    }
    result.openings.reserve(quadrature.openings.size());
    result.samples.reserve(patches.patches.size());

    for (const auto& sourceOpening : quadrature.openings) {
        SceneFluidOpeningFlux opening;
        opening.openingIndex = sourceOpening.openingIndex;
        opening.openingId = sourceOpening.openingId;
        opening.negativeSideRegionId =
            sourceOpening.negativeSideRegionId;
        opening.positiveSideRegionId =
            sourceOpening.positiveSideRegionId;
        opening.role = sourceOpening.role;
        opening.firstSample = result.samples.size();
        for (std::size_t localPoint = 0;
             localPoint < sourceOpening.pointCount; ++localPoint) {
            const std::size_t pointIndex =
                sourceOpening.firstPoint + localPoint;
            const auto& range = patches.pointRanges[pointIndex];
            for (const auto& patch : patches.patchesForPoint(range)) {
                double fluidFlow = 0.0;
                if (patch.ownerKind
                    == SceneFluidOpeningPatchOwnerKind::Face) {
                    std::size_t nextEvaluationCount = 0;
                    if (!checkedAdd(result.velocityEvaluationCount, 1,
                                    nextEvaluationCount)
                        || nextEvaluationCount
                            > limits.maximumVelocityEvaluations) {
                        throw std::length_error(
                            "scene fluid opening flux exceeds its velocity-evaluation limit");
                    }
                    result.velocityEvaluationCount = nextEvaluationCount;
                    fluidFlow = patch.areaSquareMeters
                        * faceNormalVelocity(patch, grid, velocity);
                } else {
                    fluidFlow = cellPatchFluidFlow(
                        patch, patches.verticesForPatch(patch), grid,
                        velocity, result.velocityEvaluationCount, limits);
                }
                const double fluidNormal =
                    fluidFlow / patch.areaSquareMeters;
                const double surfaceNormal =
                    patch.surfaceSweepRateCubicMetersPerSecond
                        / patch.areaSquareMeters;
                const double relativeNormal = fluidNormal - surfaceNormal;
                const double relativeFlow = fluidFlow
                    - patch.surfaceSweepRateCubicMetersPerSecond;
                if (!std::isfinite(fluidFlow)
                    || !std::isfinite(fluidNormal)
                    || !std::isfinite(surfaceNormal)
                    || !std::isfinite(relativeNormal)
                    || !std::isfinite(relativeFlow)
                    || !closeValue(
                        relativeNormal * patch.areaSquareMeters,
                        relativeFlow)) {
                    throw std::overflow_error(
                        "scene fluid opening flux sample is not finite or closed");
                }
                result.samples.push_back({
                    patch.stableId,
                    patch.sourcePointStableId,
                    patch.openingIndex,
                    patch.openingId,
                    patch.negativeSideRegionId,
                    patch.positiveSideRegionId,
                    patch.role,
                    patch.areaSquareMeters,
                    fluidNormal,
                    surfaceNormal,
                    relativeNormal,
                    fluidFlow,
                    patch.surfaceSweepRateCubicMetersPerSecond,
                    relativeFlow,
                });
                opening.areaSquareMeters += patch.areaSquareMeters;
                opening.fluidVolumeFlowRateCubicMetersPerSecond +=
                    fluidFlow;
                opening.surfaceSweepRateCubicMetersPerSecond +=
                    patch.surfaceSweepRateCubicMetersPerSecond;
                opening.relativeVolumeFlowRateCubicMetersPerSecond +=
                    relativeFlow;
            }
        }
        opening.sampleCount = result.samples.size() - opening.firstSample;
        if (opening.sampleCount == 0
            || !closeValue(
                opening.areaSquareMeters,
                sourceOpening.areaSquareMeters)
            || !closeValue(
                opening.surfaceSweepRateCubicMetersPerSecond,
                sourceOpening.surfaceSweepRateCubicMetersPerSecond)
            || !closeValue(
                opening.fluidVolumeFlowRateCubicMetersPerSecond
                    - opening.surfaceSweepRateCubicMetersPerSecond,
                opening.relativeVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid opening flux does not close its opening ledger");
        }
        result.totalAreaSquareMeters += opening.areaSquareMeters;
        result.totalFluidVolumeFlowRateCubicMetersPerSecond +=
            opening.fluidVolumeFlowRateCubicMetersPerSecond;
        result.totalSurfaceSweepRateCubicMetersPerSecond +=
            opening.surfaceSweepRateCubicMetersPerSecond;
        result.totalRelativeVolumeFlowRateCubicMetersPerSecond +=
            opening.relativeVolumeFlowRateCubicMetersPerSecond;
        const auto negativeRegion = surface.mappings.regionIndex(
            opening.negativeSideRegionId);
        const auto positiveRegion = surface.mappings.regionIndex(
            opening.positiveSideRegionId);
        if (!negativeRegion || !positiveRegion
            || *negativeRegion == *positiveRegion) {
            throw std::invalid_argument(
                "scene fluid opening flux references invalid side regions");
        }
        auto& negative = result.regions[*negativeRegion];
        auto& positive = result.regions[*positiveRegion];
        negative.outwardFluidVolumeFlowRateCubicMetersPerSecond +=
            opening.fluidVolumeFlowRateCubicMetersPerSecond;
        positive.outwardFluidVolumeFlowRateCubicMetersPerSecond -=
            opening.fluidVolumeFlowRateCubicMetersPerSecond;
        negative.outwardSurfaceSweepRateCubicMetersPerSecond +=
            opening.surfaceSweepRateCubicMetersPerSecond;
        positive.outwardSurfaceSweepRateCubicMetersPerSecond -=
            opening.surfaceSweepRateCubicMetersPerSecond;
        negative.outwardRelativeVolumeFlowRateCubicMetersPerSecond +=
            opening.relativeVolumeFlowRateCubicMetersPerSecond;
        positive.outwardRelativeVolumeFlowRateCubicMetersPerSecond -=
            opening.relativeVolumeFlowRateCubicMetersPerSecond;
        result.openings.push_back(opening);
    }
    if (!closeValue(
            result.totalAreaSquareMeters, patches.totalAreaSquareMeters)
        || !closeValue(
            result.totalSurfaceSweepRateCubicMetersPerSecond,
            patches.totalSurfaceSweepRateCubicMetersPerSecond)
        || !closeValue(
            result.totalFluidVolumeFlowRateCubicMetersPerSecond
                - result.totalSurfaceSweepRateCubicMetersPerSecond,
            result.totalRelativeVolumeFlowRateCubicMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid opening flux does not close its global ledger");
    }
    double absoluteFluidRegionFlow = 0.0;
    double absoluteSurfaceRegionFlow = 0.0;
    double absoluteRelativeRegionFlow = 0.0;
    for (const auto& region : result.regions) {
        result.globalFluidRegionBalanceResidualCubicMetersPerSecond +=
            region.outwardFluidVolumeFlowRateCubicMetersPerSecond;
        result.globalSurfaceRegionBalanceResidualCubicMetersPerSecond +=
            region.outwardSurfaceSweepRateCubicMetersPerSecond;
        result.globalRelativeRegionBalanceResidualCubicMetersPerSecond +=
            region.outwardRelativeVolumeFlowRateCubicMetersPerSecond;
        absoluteFluidRegionFlow += std::abs(
            region.outwardFluidVolumeFlowRateCubicMetersPerSecond);
        absoluteSurfaceRegionFlow += std::abs(
            region.outwardSurfaceSweepRateCubicMetersPerSecond);
        absoluteRelativeRegionFlow += std::abs(
            region.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    const auto balanced = [](const double residual,
                             const double absoluteFlow) {
        return std::isfinite(residual) && std::isfinite(absoluteFlow)
            && std::abs(residual)
                <= 1.0e-12 + 1.0e-10 * absoluteFlow;
    };
    if (!balanced(
            result.globalFluidRegionBalanceResidualCubicMetersPerSecond,
            absoluteFluidRegionFlow)
        || !balanced(
            result.globalSurfaceRegionBalanceResidualCubicMetersPerSecond,
            absoluteSurfaceRegionFlow)
        || !balanced(
            result.globalRelativeRegionBalanceResidualCubicMetersPerSecond,
            absoluteRelativeRegionFlow)) {
        throw std::overflow_error(
            "scene fluid opening region balances do not close globally");
    }
    result.ownedStorageBytes = storageBytes(result);
    result.fingerprint = fluxFingerprint(result);
    return result;
}

bool sameGrid(const SceneFluidOpeningFluxSet& flux,
              const fluid::PeriodicCartesianGrid& grid) {
    return flux.cellCounts == grid.cellCounts()
        && flux.lowerMeters == grid.lowerMeters()
        && flux.upperMeters == grid.upperMeters();
}

} // namespace

std::uint64_t sceneFluidOpeningFluxVelocityFingerprint(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocityMetersPerSecond) {
    if (!velocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(velocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid opening-flux velocity fingerprint input is invalid");
    }
    return velocityFingerprint(grid, velocityMetersPerSecond);
}

SceneFluidOpeningFluxSet evaluateSceneFluidOpeningFlux(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocityMetersPerSecond,
    const SceneFluidOpeningFluxLimits& limits) {
    validateSceneFluidOpeningGridPatches(
        patches, surface, state, caps, quadrature, grid);
    auto result = buildFlux(
        surface, quadrature, patches, grid,
        velocityMetersPerSecond, limits);
    validateSceneFluidOpeningFlux(
        result, surface, state, caps, quadrature, patches,
        grid, velocityMetersPerSecond);
    return result;
}

void validateSceneFluidOpeningFluxIntegrity(
    const SceneFluidOpeningFluxSet& flux) {
    if (flux.version != sceneFluidOpeningFluxVersion
        || flux.fingerprint == 0
        || flux.ownedStorageBytes != storageBytes(flux)
        || flux.fingerprint != fluxFingerprint(flux)) {
        throw std::invalid_argument(
            "scene fluid opening-flux integrity is invalid");
    }
}

void validateSceneFluidOpeningFlux(
    const SceneFluidOpeningFluxSet& flux,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocityMetersPerSecond) {
    validateSceneFluidOpeningGridPatches(
        patches, surface, state, caps, quadrature, grid);
    validateSceneFluidOpeningFluxIntegrity(flux);
    if (!velocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(velocityMetersPerSecond)
        || flux.version != sceneFluidOpeningFluxVersion
        || flux.surfaceDefinitionFingerprint != surface.fingerprint
        || flux.surfaceStateFingerprint != state.fingerprint
        || flux.openingPatchFingerprint != patches.fingerprint
        || flux.velocityFingerprint
            != velocityFingerprint(grid, velocityMetersPerSecond)
        || flux.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || flux.acceptedStepCount != state.acceptedStepCount
        || flux.simulationTimeSeconds != state.simulationTimeSeconds
        || !sameGrid(flux, grid)) {
        throw std::invalid_argument(
            "scene fluid opening-flux identity is invalid");
    }
    const SceneFluidOpeningFluxLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildFlux(
        surface, quadrature, patches, grid,
        velocityMetersPerSecond, unlimited);
    if (flux != expected) {
        throw std::invalid_argument(
            "scene fluid opening-flux payload is invalid");
    }
}

} // namespace simwing::fsi
