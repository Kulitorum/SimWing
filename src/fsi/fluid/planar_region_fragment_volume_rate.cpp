#include "fluid/planar_region_fragment_volume_rate.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
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

struct SurfaceMotion {
    double displacementMeters = 0.0;
    double velocityMetersPerSecond = 0.0;
};

void validateLimits(
    const PlanarPressureRegionFragmentVolumeRateLimits& limits) {
    if (limits.maximumFragments == 0 || limits.maximumCells == 0
        || limits.maximumRegions == 0 || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar regional fragment volume-rate limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional fragment volume-rate storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional fragment volume-rate storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t fragmentCount,
                              const std::size_t cellCount,
                              const std::size_t regionCount,
                              const std::size_t componentCount) {
    std::size_t result = checkedProduct(
        fragmentCount, sizeof(PlanarPressureRegionFragmentVolumeRate));
    result = checkedSum(result, checkedProduct(
        cellCount, sizeof(PlanarPressureRegionFragmentCellVolumeRate)));
    result = checkedSum(result, checkedProduct(
        regionCount, sizeof(PlanarPressureRegionFragmentRegionVolumeRate)));
    return checkedSum(result, checkedProduct(
        componentCount,
        sizeof(PlanarPressureRegionFragmentComponentVolumeRate)));
}

double axisLower(const PeriodicCartesianGrid& grid,
                 const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return grid.lowerMeters().x;
    case GridFaceAxis::Y: return grid.lowerMeters().y;
    case GridFaceAxis::Z: return grid.lowerMeters().z;
    }
    throw std::invalid_argument(
        "planar regional fragment volume-rate axis is invalid");
}

double axisSpacing(const PeriodicCartesianGrid& grid,
                   const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return grid.cellSpacingMeters().x;
    case GridFaceAxis::Y: return grid.cellSpacingMeters().y;
    case GridFaceAxis::Z: return grid.cellSpacingMeters().z;
    }
    throw std::invalid_argument(
        "planar regional fragment volume-rate axis is invalid");
}

std::int64_t segmentOrdinal(const double coordinateMeters,
                            const double lowerMeters,
                            const double spacingMeters) {
    const double ordinal = std::floor(
        (coordinateMeters - lowerMeters) / spacingMeters);
    if (!std::isfinite(ordinal)
        || ordinal < static_cast<double>(
            std::numeric_limits<std::int64_t>::min())
        || ordinal > static_cast<double>(
            std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error(
            "planar regional fragment volume-rate segment overflows");
    }
    return static_cast<std::int64_t>(ordinal);
}

void insertSurfaceMotion(std::map<std::uint64_t, SurfaceMotion>& motions,
                         const std::uint64_t stableId,
                         const double displacementMeters,
                         const double velocityMetersPerSecond) {
    const SurfaceMotion motion{
        displacementMeters, velocityMetersPerSecond};
    const auto [found, inserted] = motions.emplace(stableId, motion);
    if (!inserted && (found->second.displacementMeters
                          != motion.displacementMeters
                      || found->second.velocityMetersPerSecond
                          != motion.velocityMetersPerSecond)) {
        throw std::invalid_argument(
            "planar regional fragment surface motion is inconsistent");
    }
}

std::map<std::uint64_t, SurfaceMotion> surfaceMotions(
    const PlanarPressureRegionSweepLedger& sweep) {
    std::map<std::uint64_t, SurfaceMotion> result;
    for (const auto& interval : sweep.intervals) {
        insertSurfaceMotion(
            result, interval.lowerSurfaceStableId,
            interval.lowerSurfaceDisplacementMeters,
            interval.lowerSurfaceVelocityMetersPerSecond);
        insertSurfaceMotion(
            result, interval.upperSurfaceStableId,
            interval.upperSurfaceDisplacementMeters,
            interval.upperSurfaceVelocityMetersPerSecond);
    }
    return result;
}

SurfaceMotion boundaryMotion(
    const PlanarPressureRegionFragmentBoundary& boundary,
    const std::map<std::uint64_t, SurfaceMotion>& motions) {
    if (boundary.kind
        == PlanarPressureRegionFragmentBoundaryKind::GridFace) {
        return {};
    }
    const auto found = motions.find(boundary.surfaceStableId);
    if (found == motions.end()) {
        throw std::invalid_argument(
            "planar regional fragment boundary motion is missing");
    }
    return found->second;
}

double closureTolerance(const double scale) {
    return 4096.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale);
}

std::uint64_t setFingerprint(
    const PlanarPressureRegionFragmentVolumeRateSet& set) {
    Fingerprint fingerprint;
    fingerprint.integer(set.version);
    fingerprint.integer(set.sourceFragmentFingerprint);
    fingerprint.integer(set.sourceTopologyFingerprint);
    fingerprint.integer(set.sourceSweepVersion);
    fingerprint.enumeration(set.axis);
    fingerprint.real(set.durationSeconds);
    fingerprint.integer(static_cast<std::uint8_t>(set.topologyStable));
    fingerprint.integer(static_cast<std::uint64_t>(set.fragments.size()));
    for (const auto& fragment : set.fragments) {
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.fragmentIndex));
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(fragment.cellIndex));
        for (const double value : {
                 fragment.lowerBoundaryDisplacementMeters,
                 fragment.upperBoundaryDisplacementMeters,
                 fragment.lowerBoundaryVelocityMetersPerSecond,
                 fragment.upperBoundaryVelocityMetersPerSecond,
                 fragment.previousVolumeCubicMeters,
                 fragment.currentVolumeCubicMeters,
                 fragment.geometryVolumeChangeCubicMeters,
                 fragment.geometryVolumeChangeRateCubicMetersPerSecond}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.cells.size()));
    for (const auto& cell : set.cells) {
        fingerprint.integer(static_cast<std::uint64_t>(cell.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(cell.fragmentCount));
        for (const double value : {
                 cell.previousVolumeCubicMeters,
                 cell.currentVolumeCubicMeters,
                 cell.geometryVolumeChangeCubicMeters,
                 cell.geometryVolumeChangeRateCubicMetersPerSecond,
                 cell.fixedCellVolumeClosureResidualCubicMeters}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.regions.size()));
    for (const auto& region : set.regions) {
        fingerprint.integer(region.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(region.fragmentCount));
        for (const double value : {
                 region.previousVolumeCubicMeters,
                 region.currentVolumeCubicMeters,
                 region.geometryVolumeChangeCubicMeters,
                 region.sourceSweepGeometryVolumeChangeCubicMeters,
                 region.geometryVolumeChangeClosureResidualCubicMeters,
                 region.geometryVolumeChangeRateCubicMetersPerSecond}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.components.size()));
    for (const auto& component : set.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.fragmentCount));
        for (const double value : {
                 component.previousVolumeCubicMeters,
                 component.currentVolumeCubicMeters,
                 component.geometryVolumeChangeCubicMeters,
                 component.geometryVolumeChangeRateCubicMetersPerSecond}) {
            fingerprint.real(value);
        }
    }
    for (const double value : {
             set.maximumAbsoluteFragmentVolumeChangeCubicMeters,
             set.maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond,
             set.maximumAbsoluteCellClosureResidualCubicMeters,
             set.maximumAbsoluteRegionClosureResidualCubicMeters,
             set.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
             set.globalGeometryVolumeChangeCubicMeters,
             set.globalGeometryVolumeChangeRateCubicMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(set.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentVolumeRateSet buildSet(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentTopology(
        topology, grid, sweep, fragments, limits.topologyLimits);
    if (fragments.fragments.size() > limits.maximumFragments
        || fragments.cells.size() > limits.maximumCells
        || fragments.regions.size() > limits.maximumRegions
        || topology.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional fragment volume-rate count limit exceeded");
    }

    PlanarPressureRegionFragmentVolumeRateSet result;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceSweepVersion = sweep.version;
    result.axis = sweep.axis;
    result.durationSeconds = sweep.durationSeconds;
    result.topologyStable = true;
    result.ownedStorageBytes = ownedStorageBytes(
        fragments.fragments.size(), fragments.cells.size(),
        fragments.regions.size(), topology.components.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional fragment volume-rate storage limit exceeded");
    }

    const auto motions = surfaceMotions(sweep);
    const double lowerMeters = axisLower(grid, sweep.axis);
    const double spacingMeters = axisSpacing(grid, sweep.axis);
    for (const auto& fragment : fragments.fragments) {
        for (const auto* boundary
             : {&fragment.lowerBoundary, &fragment.upperBoundary}) {
            if (boundary->kind
                != PlanarPressureRegionFragmentBoundaryKind::PressureLayer) {
                continue;
            }
            const SurfaceMotion motion = boundaryMotion(*boundary, motions);
            const double previousCoordinate =
                boundary->unwrappedCoordinateMeters
                - motion.displacementMeters;
            if (segmentOrdinal(
                    boundary->unwrappedCoordinateMeters,
                    lowerMeters, spacingMeters)
                != segmentOrdinal(
                    previousCoordinate, lowerMeters, spacingMeters)) {
                throw std::invalid_argument(
                    "planar regional fragment volume rates require topology-stable layer motion");
            }
        }

        const SurfaceMotion lower = boundaryMotion(
            fragment.lowerBoundary, motions);
        const SurfaceMotion upper = boundaryMotion(
            fragment.upperBoundary, motions);
        PlanarPressureRegionFragmentVolumeRate rate;
        rate.fragmentIndex = result.fragments.size();
        rate.stableId = fragment.stableId;
        rate.regionStableId = fragment.regionStableId;
        rate.componentIndex = topology.fragments[
            rate.fragmentIndex].componentIndex;
        rate.cellIndex = grid.cellIndex(
            fragment.i, fragment.j, fragment.k);
        rate.lowerBoundaryDisplacementMeters = lower.displacementMeters;
        rate.upperBoundaryDisplacementMeters = upper.displacementMeters;
        rate.lowerBoundaryVelocityMetersPerSecond =
            lower.velocityMetersPerSecond;
        rate.upperBoundaryVelocityMetersPerSecond =
            upper.velocityMetersPerSecond;
        rate.currentVolumeCubicMeters = fragment.volumeCubicMeters;
        const double boundaryVolumeChangeCubicMeters =
            fragment.transverseAreaSquareMeters
            * (upper.displacementMeters - lower.displacementMeters);
        rate.previousVolumeCubicMeters =
            rate.currentVolumeCubicMeters
            - boundaryVolumeChangeCubicMeters;
        rate.geometryVolumeChangeCubicMeters =
            rate.currentVolumeCubicMeters
            - rate.previousVolumeCubicMeters;
        rate.geometryVolumeChangeRateCubicMetersPerSecond =
            rate.geometryVolumeChangeCubicMeters / sweep.durationSeconds;
        if (!std::isfinite(rate.previousVolumeCubicMeters)
            || !(rate.previousVolumeCubicMeters > 0.0)
            || !std::isfinite(
                rate.geometryVolumeChangeRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "planar regional fragment reconstructed volume is invalid");
        }
        result.maximumAbsoluteFragmentVolumeChangeCubicMeters = std::max(
            result.maximumAbsoluteFragmentVolumeChangeCubicMeters,
            std::abs(rate.geometryVolumeChangeCubicMeters));
        result.maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond,
                std::abs(
                    rate.geometryVolumeChangeRateCubicMetersPerSecond));
        result.fragments.push_back(rate);
    }

    result.cells.resize(fragments.cells.size());
    for (const auto& source : fragments.cells) {
        auto& cell = result.cells[source.cellIndex];
        cell.cellIndex = source.cellIndex;
    }
    result.regions.reserve(fragments.regions.size());
    std::map<std::uint64_t, std::size_t> regionIndices;
    for (const auto& source : fragments.regions) {
        const std::size_t index = result.regions.size();
        regionIndices.emplace(source.regionStableId, index);
        result.regions.push_back({source.regionStableId});
    }
    result.components.reserve(topology.components.size());
    for (const auto& source : topology.components) {
        result.components.push_back({
            source.componentIndex, source.stableId, source.regionStableId});
    }
    for (const auto& rate : result.fragments) {
        auto& cell = result.cells[rate.cellIndex];
        ++cell.fragmentCount;
        cell.previousVolumeCubicMeters += rate.previousVolumeCubicMeters;
        cell.currentVolumeCubicMeters += rate.currentVolumeCubicMeters;
        cell.geometryVolumeChangeCubicMeters +=
            rate.geometryVolumeChangeCubicMeters;
        cell.geometryVolumeChangeRateCubicMetersPerSecond +=
            rate.geometryVolumeChangeRateCubicMetersPerSecond;

        auto& region = result.regions.at(
            regionIndices.at(rate.regionStableId));
        ++region.fragmentCount;
        region.previousVolumeCubicMeters += rate.previousVolumeCubicMeters;
        region.currentVolumeCubicMeters += rate.currentVolumeCubicMeters;
        region.geometryVolumeChangeCubicMeters +=
            rate.geometryVolumeChangeCubicMeters;
        region.geometryVolumeChangeRateCubicMetersPerSecond +=
            rate.geometryVolumeChangeRateCubicMetersPerSecond;

        auto& component = result.components.at(rate.componentIndex);
        ++component.fragmentCount;
        component.previousVolumeCubicMeters +=
            rate.previousVolumeCubicMeters;
        component.currentVolumeCubicMeters +=
            rate.currentVolumeCubicMeters;
        component.geometryVolumeChangeCubicMeters +=
            rate.geometryVolumeChangeCubicMeters;
        component.geometryVolumeChangeRateCubicMetersPerSecond +=
            rate.geometryVolumeChangeRateCubicMetersPerSecond;
    }

    const double cellVolume = grid.cellVolumeCubicMeters();
    for (auto& cell : result.cells) {
        cell.fixedCellVolumeClosureResidualCubicMeters =
            cell.geometryVolumeChangeCubicMeters;
        result.maximumAbsoluteCellClosureResidualCubicMeters = std::max(
            result.maximumAbsoluteCellClosureResidualCubicMeters,
            std::abs(cell.fixedCellVolumeClosureResidualCubicMeters));
        const double scale = std::max({
            cellVolume, std::abs(cell.previousVolumeCubicMeters),
            std::abs(cell.currentVolumeCubicMeters)});
        if (cell.fragmentCount
                != fragments.cells[cell.cellIndex].fragmentCount
            || std::abs(cell.previousVolumeCubicMeters - cellVolume)
                > closureTolerance(scale)
            || std::abs(cell.currentVolumeCubicMeters - cellVolume)
                > closureTolerance(scale)
            || std::abs(cell.fixedCellVolumeClosureResidualCubicMeters)
                > closureTolerance(scale)) {
            throw std::invalid_argument(
                "planar regional fragment cell volume-rate closure failed");
        }
    }
    for (auto& region : result.regions) {
        const auto source = std::ranges::find(
            sweep.regions, region.regionStableId,
            &PlanarPressureRegionSweepSummary::regionStableId);
        if (source == sweep.regions.end()) {
            throw std::invalid_argument(
                "planar regional fragment volume-rate region is missing");
        }
        region.sourceSweepGeometryVolumeChangeCubicMeters =
            source->geometryVolumeChangeCubicMeters;
        region.geometryVolumeChangeClosureResidualCubicMeters =
            region.geometryVolumeChangeCubicMeters
            - source->geometryVolumeChangeCubicMeters;
        result.maximumAbsoluteRegionClosureResidualCubicMeters = std::max(
            result.maximumAbsoluteRegionClosureResidualCubicMeters,
            std::abs(
                region.geometryVolumeChangeClosureResidualCubicMeters));
        const double scale = std::max({
            std::abs(source->previousVolumeCubicMeters),
            std::abs(source->currentVolumeCubicMeters), 1.0});
        const auto fragmentSource = std::ranges::find(
            fragments.regions, region.regionStableId,
            &PlanarPressureRegionFragmentRegionSummary::regionStableId);
        if (fragmentSource == fragments.regions.end()
            || region.fragmentCount != fragmentSource->fragmentCount
            || std::abs(
                region.previousVolumeCubicMeters
                - source->previousVolumeCubicMeters)
                > closureTolerance(scale)
            || std::abs(
                region.currentVolumeCubicMeters
                - source->currentVolumeCubicMeters)
                > closureTolerance(scale)
            || std::abs(
                region.geometryVolumeChangeClosureResidualCubicMeters)
                > closureTolerance(scale)) {
            throw std::invalid_argument(
                "planar regional fragment region volume-rate closure failed");
        }
    }
    for (const auto& component : result.components) {
        const auto& source = topology.components.at(
            component.componentIndex);
        const double scale = std::max({
            std::abs(component.currentVolumeCubicMeters),
            std::abs(source.volumeCubicMeters), 1.0});
        if (component.stableId != source.stableId
            || component.regionStableId != source.regionStableId
            || component.fragmentCount != source.fragmentCount
            || std::abs(
                component.currentVolumeCubicMeters
                - source.volumeCubicMeters) > closureTolerance(scale)) {
            throw std::invalid_argument(
                "planar regional fragment component volume-rate closure failed");
        }
        result.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
                std::abs(
                    component.geometryVolumeChangeRateCubicMetersPerSecond));
    }
    for (const auto& region : result.regions) {
        result.globalGeometryVolumeChangeCubicMeters +=
            region.geometryVolumeChangeCubicMeters;
    }
    result.globalGeometryVolumeChangeRateCubicMetersPerSecond =
        result.globalGeometryVolumeChangeCubicMeters / sweep.durationSeconds;
    const double domainVolume = grid.cellVolumeCubicMeters()
        * static_cast<double>(grid.cellCount());
    if (std::abs(result.globalGeometryVolumeChangeCubicMeters)
            > closureTolerance(domainVolume)
        || std::abs(
            result.globalGeometryVolumeChangeCubicMeters
            - sweep.globalGeometryVolumeChangeCubicMeters)
            > closureTolerance(domainVolume)) {
        throw std::invalid_argument(
            "planar regional fragment global volume-rate closure failed");
    }
    result.fingerprint = setFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentVolumeRateSet
buildPlanarPressureRegionFragmentVolumeRates(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateLimits& limits) {
    return buildSet(grid, sweep, fragments, topology, limits);
}

void validatePlanarPressureRegionFragmentVolumeRates(
    const PlanarPressureRegionFragmentVolumeRateSet& rates,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateLimits& limits) {
    if (rates != buildSet(grid, sweep, fragments, topology, limits)) {
        throw std::invalid_argument(
            "planar regional fragment volume-rate set is corrupted");
    }
}

} // namespace simwing::fsi::fluid
