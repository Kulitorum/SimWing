#include "scene_fluid_cell_volume.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
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

struct EdgeIncidence {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t negativeRegion = 0;
    std::size_t positiveRegion = 0;
};

struct Contribution {
    std::size_t cellIndex = 0;
    std::size_t regionIndex = 0;
    std::size_t ordinal = 0;
    double volumeCubicMeters = 0.0;
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

double tolerance(const SceneFluidCellVolumeSettings& settings,
                 const double referenceVolume) {
    return std::max(
        settings.absoluteVolumeToleranceCubicMeters,
        settings.relativeVolumeTolerance * std::abs(referenceVolume));
}

void validateSettings(const SceneFluidCellVolumeSettings& settings) {
    if (!std::isfinite(settings.absoluteVolumeToleranceCubicMeters)
        || settings.absoluteVolumeToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeVolumeTolerance)
        || settings.relativeVolumeTolerance < 0.0
        || (settings.absoluteVolumeToleranceCubicMeters == 0.0
            && settings.relativeVolumeTolerance == 0.0)) {
        throw std::invalid_argument(
            "scene fluid cell-volume tolerances are invalid");
    }
}

std::size_t outsideRegionIndex(
    const SceneFluidSurfaceDefinition& surface) {
    std::size_t result = surface.regions.size();
    for (std::size_t index = 0; index < surface.regions.size(); ++index) {
        if (surface.regions[index].kind != RegionKind::Outside) {
            continue;
        }
        if (result != surface.regions.size()) {
            throw std::invalid_argument(
                "scene fluid cell volumes require exactly one Outside region");
        }
        result = index;
    }
    if (result == surface.regions.size()) {
        throw std::invalid_argument(
            "scene fluid cell volumes require one Outside region");
    }
    return result;
}

void validateClosedManifold(
    const SceneFluidSurfaceDefinition& surface) {
    if (!surface.openings.empty()) {
        throw std::invalid_argument(
            "scene fluid cell volumes do not yet support authored openings");
    }
    using Edge = std::pair<std::size_t, std::size_t>;
    std::map<Edge, std::vector<EdgeIncidence>> edges;
    for (const auto& triangle : surface.triangles) {
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const std::size_t from = triangle.vertexIndices[edge];
            const std::size_t to = triangle.vertexIndices[(edge + 1) % 3];
            edges[{std::min(from, to), std::max(from, to)}].push_back({
                from, to,
                triangle.negativeSideRegionIndex,
                triangle.positiveSideRegionIndex,
            });
        }
    }
    for (const auto& [edge, adjacent] : edges) {
        static_cast<void>(edge);
        if (adjacent.size() != 2
            || adjacent[0].from != adjacent[1].to
            || adjacent[0].to != adjacent[1].from
            || adjacent[0].negativeRegion != adjacent[1].negativeRegion
            || adjacent[0].positiveRegion != adjacent[1].positiveRegion) {
            throw std::invalid_argument(
                "scene fluid cell volumes require closed consistently wound two-sided manifolds");
        }
    }
}

void validateResolvedSubset(const SceneFluidGridEpoch& epoch) {
    if (!epoch.ownership.facePatches.empty()
        || epoch.crossings.unpairedContactSegmentCount != 0
        || epoch.facePartitions.unresolvedActiveFaceCount != 0
        || epoch.facePartitions.partitions.size()
            != epoch.faceTopology.activeFaces.size()) {
        throw std::invalid_argument(
            "scene fluid cell volumes require fully resolved transverse face partitions");
    }
    std::vector<bool> resolved(epoch.faceTopology.activeFaces.size(), false);
    for (const auto& partition : epoch.facePartitions.partitions) {
        if (partition.activeFaceIndex >= resolved.size()
            || resolved[partition.activeFaceIndex]) {
            throw std::invalid_argument(
                "scene fluid cell-volume face partition ownership is invalid");
        }
        resolved[partition.activeFaceIndex] = true;
    }
    if (!std::ranges::all_of(resolved, [](const bool value) {
            return value;
        })) {
        throw std::invalid_argument(
            "scene fluid cell-volume face partitions are incomplete");
    }
}

std::size_t regionIndex(const SceneFluidSurfaceDefinition& surface,
                        const StableId id) {
    const auto index = surface.mappings.regionIndex(id);
    if (!index) {
        throw std::invalid_argument(
            "scene fluid cell-volume face references an unknown region");
    }
    return *index;
}

void appendContribution(std::vector<Contribution>& contributions,
                        const SceneFluidCellVolumeLimits& limits,
                        const std::size_t cellIndex,
                        const std::size_t region,
                        const double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "scene fluid cell-volume contribution is non-finite");
    }
    if (contributions.size() == limits.maximumContributionEvents) {
        throw std::length_error(
            "scene fluid cell-volume contributions exceed their limit");
    }
    contributions.push_back({
        cellIndex, region, contributions.size(), value});
}

std::vector<const fluid::SceneFluidFacePartition*>
partitionByActiveFace(const SceneFluidGridEpoch& epoch) {
    std::vector<const fluid::SceneFluidFacePartition*> result(
        epoch.faceTopology.activeFaces.size(), nullptr);
    for (const auto& partition : epoch.facePartitions.partitions) {
        result[partition.activeFaceIndex] = &partition;
    }
    return result;
}

double faceCoordinate(const fluid::SceneFluidActiveFace& face,
                      const fluid::PeriodicCartesianGrid& grid) {
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    if (face.axis == fluid::GridFaceAxis::X) {
        return lower.x + static_cast<double>(face.i) * spacing.x;
    }
    if (face.axis == fluid::GridFaceAxis::Y) {
        return lower.y + static_cast<double>(face.j) * spacing.y;
    }
    return lower.z + static_cast<double>(face.k) * spacing.z;
}

double faceArea(const fluid::GridFaceAxis axis,
                const fluid::PeriodicCartesianGrid& grid) {
    const auto spacing = grid.cellSpacingMeters();
    if (axis == fluid::GridFaceAxis::X) {
        return spacing.y * spacing.z;
    }
    if (axis == fluid::GridFaceAxis::Y) {
        return spacing.z * spacing.x;
    }
    return spacing.x * spacing.y;
}

std::pair<std::size_t, std::size_t> adjacentCells(
    const fluid::SceneFluidActiveFace& face,
    const fluid::PeriodicCartesianGrid& grid) {
    if (face.axis == fluid::GridFaceAxis::X) {
        return {grid.cellIndex(face.i - 1, face.j, face.k),
                grid.cellIndex(face.i, face.j, face.k)};
    }
    if (face.axis == fluid::GridFaceAxis::Y) {
        return {grid.cellIndex(face.i, face.j - 1, face.k),
                grid.cellIndex(face.i, face.j, face.k)};
    }
    return {grid.cellIndex(face.i, face.j, face.k - 1),
            grid.cellIndex(face.i, face.j, face.k)};
}

std::vector<double> wholeSurfaceVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const std::size_t outsideIndex,
    const fluid::PeriodicCartesianGrid& grid) {
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    const double domainVolume = (upper.x - lower.x)
        * (upper.y - lower.y) * (upper.z - lower.z);
    std::vector<double> result(surface.regions.size(), 0.0);
    result[outsideIndex] = domainVolume;
    for (const auto& triangle : surface.triangles) {
        const auto& first =
            state.vertices[triangle.vertexIndices[0]].positionMeters;
        const auto& second =
            state.vertices[triangle.vertexIndices[1]].positionMeters;
        const auto& third =
            state.vertices[triangle.vertexIndices[2]].positionMeters;
        const Vec3 areaVector = cross(
            subtract(second, first), subtract(third, first));
        const Vec3 centroid{
            (first.x + second.x + third.x) / 3.0,
            (first.y + second.y + third.y) / 3.0,
            (first.z + second.z + third.z) / 3.0,
        };
        const double contribution = dot(centroid, areaVector) / 6.0;
        result[triangle.negativeSideRegionIndex] += contribution;
        result[triangle.positiveSideRegionIndex] -= contribution;
    }
    return result;
}

std::uint64_t volumeFingerprint(const SceneFluidCellVolumeSet& volumes) {
    Fingerprint fingerprint;
    fingerprint.integer(volumes.version);
    fingerprint.integer(volumes.surfaceDefinitionFingerprint);
    fingerprint.integer(volumes.surfaceStateFingerprint);
    fingerprint.integer(volumes.gridEpochFingerprint);
    fingerprint.integer(volumes.structureDefinitionFingerprint);
    fingerprint.integer(volumes.acceptedStepCount);
    fingerprint.real(volumes.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cellCounts.z));
    for (const double value : {
             volumes.lowerMeters.x, volumes.lowerMeters.y,
             volumes.lowerMeters.z, volumes.upperMeters.x,
             volumes.upperMeters.y, volumes.upperMeters.z,
             volumes.settings.absoluteVolumeToleranceCubicMeters,
             volumes.settings.relativeVolumeTolerance,
             volumes.cellVolumeCubicMeters,
             volumes.maximumCellVolumeResidualCubicMeters,
             volumes.maximumRegionVolumeResidualCubicMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(volumes.outsideRegionId);
    fingerprint.integer(static_cast<std::uint64_t>(volumes.cells.size()));
    for (const auto& cell : volumes.cells) {
        fingerprint.integer(static_cast<std::uint64_t>(cell.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.i));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.j));
        fingerprint.integer(static_cast<std::uint64_t>(cell.cell.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.firstRegionVolume));
        fingerprint.integer(static_cast<std::uint64_t>(
            cell.regionVolumeCount));
        fingerprint.real(cell.assignedVolumeCubicMeters);
        fingerprint.real(cell.volumeResidualCubicMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.cellRegionVolumes.size()));
    for (const auto& value : volumes.cellRegionVolumes) {
        fingerprint.integer(value.regionId);
        fingerprint.real(value.volumeCubicMeters);
        fingerprint.real(value.volumeFraction);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        volumes.regionVolumes.size()));
    for (const auto& value : volumes.regionVolumes) {
        fingerprint.integer(value.regionId);
        fingerprint.real(value.summedCellVolumeCubicMeters);
        fingerprint.real(value.wholeSurfaceVolumeCubicMeters);
        fingerprint.real(value.volumeResidualCubicMeters);
    }
    return fingerprint.value();
}

SceneFluidCellVolumeSet buildVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidCellVolumeSettings& settings,
    const SceneFluidCellVolumeLimits& limits) {
    validateSettings(settings);
    validateClosedManifold(surface);
    validateResolvedSubset(epoch);
    const std::size_t outsideIndex = outsideRegionIndex(surface);
    if (grid.cellCount() > limits.maximumCells
        || grid.cellCount() > limits.maximumCellRegionVolumes) {
        throw std::length_error(
            "scene fluid cell volumes exceed their cell limits");
    }

    SceneFluidCellVolumeSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.gridEpochFingerprint = epoch.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    result.outsideRegionId = surface.regions[outsideIndex].id;
    result.cellVolumeCubicMeters = grid.cellVolumeCubicMeters();

    std::vector<Contribution> contributions;
    contributions.reserve(std::min(
        limits.maximumContributionEvents, grid.cellCount()));
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        appendContribution(
            contributions, limits, cell, outsideIndex,
            result.cellVolumeCubicMeters);
    }

    for (const auto& owned : epoch.ownership.cellPatches) {
        const auto& patch = epoch.patches.patches[owned.sourcePatchIndex];
        const auto& triangle = surface.triangles[patch.triangleIndex];
        const auto& first =
            state.vertices[triangle.vertexIndices[0]].positionMeters;
        const auto& second =
            state.vertices[triangle.vertexIndices[1]].positionMeters;
        const auto& third =
            state.vertices[triangle.vertexIndices[2]].positionMeters;
        const Vec3 normal = cross(
            subtract(second, first), subtract(third, first));
        const double magnitude = length(normal);
        if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
            throw std::invalid_argument(
                "scene fluid cell-volume triangle normal is invalid");
        }
        const double contribution = dot(patch.centroidMeters, normal)
            * patch.areaSquareMeters / (3.0 * magnitude);
        appendContribution(
            contributions, limits, patch.cellIndex,
            triangle.negativeSideRegionIndex, contribution);
        appendContribution(
            contributions, limits, patch.cellIndex,
            triangle.positiveSideRegionIndex, -contribution);
    }

    const auto partitions = partitionByActiveFace(epoch);
    for (std::size_t faceIndex = 0;
         faceIndex < epoch.faceTopology.activeFaces.size(); ++faceIndex) {
        const auto& face = epoch.faceTopology.activeFaces[faceIndex];
        const auto& partition = *partitions[faceIndex];
        const auto [lowerCell, upperCell] = adjacentCells(face, grid);
        const double coefficient = faceCoordinate(face, grid) / 3.0;
        const double fullArea = faceArea(face.axis, grid);
        for (const auto [cell, sign] : {
                 std::pair{lowerCell, 1.0},
                 std::pair{upperCell, -1.0}}) {
            appendContribution(
                contributions, limits, cell, outsideIndex,
                -sign * coefficient * fullArea);
            for (std::size_t offset = 0;
                 offset < partition.regionAreaCount; ++offset) {
                const auto& area = epoch.facePartitions.regionAreas[
                    partition.firstRegionArea + offset];
                appendContribution(
                    contributions, limits, cell,
                    regionIndex(surface, area.regionId),
                    sign * coefficient * area.areaSquareMeters);
            }
        }
    }

    std::ranges::sort(
        contributions,
        [](const Contribution& first, const Contribution& second) {
            return std::tie(
                       first.cellIndex, first.regionIndex, first.ordinal)
                < std::tie(
                    second.cellIndex, second.regionIndex, second.ordinal);
        });
    const double cellTolerance = tolerance(
        settings, result.cellVolumeCubicMeters);
    std::vector<double> summedRegions(surface.regions.size(), 0.0);
    result.cells.reserve(grid.cellCount());
    std::size_t contributionIndex = 0;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t cellIndex = grid.cellIndex(i, j, k);
                SceneFluidCellVolume cell;
                cell.cellIndex = cellIndex;
                cell.cell = {i, j, k};
                cell.firstRegionVolume = result.cellRegionVolumes.size();
                while (contributionIndex < contributions.size()
                       && contributions[contributionIndex].cellIndex
                           == cellIndex) {
                    const std::size_t region =
                        contributions[contributionIndex].regionIndex;
                    double volume = 0.0;
                    while (contributionIndex < contributions.size()
                           && contributions[contributionIndex].cellIndex
                               == cellIndex
                           && contributions[contributionIndex].regionIndex
                               == region) {
                        volume += contributions[contributionIndex]
                                      .volumeCubicMeters;
                        ++contributionIndex;
                    }
                    if (!std::isfinite(volume) || volume < -cellTolerance) {
                        throw std::invalid_argument(
                            "scene fluid cell region has a negative or non-finite volume");
                    }
                    if (std::abs(volume) <= cellTolerance) {
                        volume = 0.0;
                    }
                    if (volume == 0.0) {
                        continue;
                    }
                    if (result.cellRegionVolumes.size()
                        == limits.maximumCellRegionVolumes) {
                        throw std::length_error(
                            "scene fluid cell-region volumes exceed their limit");
                    }
                    result.cellRegionVolumes.push_back({
                        surface.regions[region].id,
                        volume,
                        volume / result.cellVolumeCubicMeters,
                    });
                    cell.assignedVolumeCubicMeters += volume;
                    summedRegions[region] += volume;
                }
                cell.regionVolumeCount = result.cellRegionVolumes.size()
                    - cell.firstRegionVolume;
                cell.volumeResidualCubicMeters =
                    cell.assignedVolumeCubicMeters
                    - result.cellVolumeCubicMeters;
                result.maximumCellVolumeResidualCubicMeters = std::max(
                    result.maximumCellVolumeResidualCubicMeters,
                    std::abs(cell.volumeResidualCubicMeters));
                if (!std::isfinite(cell.volumeResidualCubicMeters)
                    || std::abs(cell.volumeResidualCubicMeters)
                        > cellTolerance) {
                    throw std::invalid_argument(
                        "scene fluid cell-region volumes do not close the cell");
                }
                result.cells.push_back(cell);
            }
        }
    }
    if (contributionIndex != contributions.size()) {
        throw std::logic_error(
            "scene fluid cell-volume contribution escaped the grid");
    }

    const auto whole = wholeSurfaceVolumes(
        surface, state, outsideIndex, grid);
    const auto lower = grid.lowerMeters();
    const auto upper = grid.upperMeters();
    const double domainVolume = (upper.x - lower.x)
        * (upper.y - lower.y) * (upper.z - lower.z);
    const double regionTolerance = tolerance(settings, domainVolume);
    result.regionVolumes.reserve(surface.regions.size());
    for (std::size_t region = 0; region < surface.regions.size(); ++region) {
        double expected = whole[region];
        if (!std::isfinite(expected) || expected < -regionTolerance) {
            throw std::invalid_argument(
                "scene fluid whole-surface region volume is invalid");
        }
        if (std::abs(expected) <= regionTolerance) {
            expected = 0.0;
        }
        const double residual = summedRegions[region] - expected;
        result.maximumRegionVolumeResidualCubicMeters = std::max(
            result.maximumRegionVolumeResidualCubicMeters,
            std::abs(residual));
        if (!std::isfinite(residual)
            || std::abs(residual) > regionTolerance) {
            throw std::invalid_argument(
                "scene fluid cell volumes do not close their whole-surface region volume");
        }
        result.regionVolumes.push_back({
            surface.regions[region].id,
            summedRegions[region],
            expected,
            residual,
        });
    }

    std::size_t cellBytes = 0;
    std::size_t cellRegionBytes = 0;
    std::size_t regionBytes = 0;
    std::size_t first = 0;
    std::size_t total = 0;
    if (!checkedMultiply(result.cells.size(),
                         sizeof(SceneFluidCellVolume), cellBytes)
        || !checkedMultiply(result.cellRegionVolumes.size(),
                            sizeof(SceneFluidCellRegionVolume),
                            cellRegionBytes)
        || !checkedMultiply(result.regionVolumes.size(),
                            sizeof(SceneFluidRegionVolume), regionBytes)
        || !checkedAdd(cellBytes, cellRegionBytes, first)
        || !checkedAdd(first, regionBytes, total)
        || total > limits.maximumVolumeBytes) {
        throw std::length_error(
            "scene fluid cell-volume result exceeds its byte limit");
    }
    result.fingerprint = volumeFingerprint(result);
    return result;
}

} // namespace

SceneFluidCellVolumeSet buildSceneFluidCellVolumes(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidCellVolumeSettings& settings,
    const SceneFluidCellVolumeLimits& limits) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    auto result = buildVolumes(
        surface, state, grid, epoch, settings, limits);
    validateSceneFluidCellVolumes(
        result, surface, state, grid, transfer, epoch);
    return result;
}

void validateSceneFluidCellVolumes(
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    if (volumes.version != sceneFluidCellVolumeVersion
        || volumes.fingerprint == 0
        || volumes.surfaceDefinitionFingerprint != surface.fingerprint
        || volumes.surfaceStateFingerprint != state.fingerprint
        || volumes.gridEpochFingerprint != epoch.fingerprint
        || volumes.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || volumes.acceptedStepCount != state.acceptedStepCount
        || volumes.simulationTimeSeconds != state.simulationTimeSeconds
        || volumes.cellCounts != grid.cellCounts()
        || volumes.lowerMeters != grid.lowerMeters()
        || volumes.upperMeters != grid.upperMeters()) {
        throw std::invalid_argument(
            "scene fluid cell-volume identity is invalid");
    }
    const SceneFluidCellVolumeLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildVolumes(
        surface, state, grid, epoch, volumes.settings, unlimited);
    if (volumes != expected
        || volumes.fingerprint != volumeFingerprint(volumes)) {
        throw std::invalid_argument(
            "scene fluid cell-volume payload is invalid");
    }
}

} // namespace simwing::fsi
