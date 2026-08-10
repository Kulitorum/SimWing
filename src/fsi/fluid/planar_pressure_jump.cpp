#include "fluid/planar_pressure_jump.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

double domainLength(const PeriodicCartesianGrid& grid,
                    const GridFaceAxis axis) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    switch (axis) {
    case GridFaceAxis::X:
        return upper.x - lower.x;
    case GridFaceAxis::Y:
        return upper.y - lower.y;
    case GridFaceAxis::Z:
        return upper.z - lower.z;
    }
    throw std::invalid_argument(
        "planar pressure-jump axis is invalid");
}

double transverseArea(const PeriodicCartesianGrid& grid,
                      const GridFaceAxis axis) {
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    double area = 0.0;
    switch (axis) {
    case GridFaceAxis::X:
        area = (upper.y - lower.y) * (upper.z - lower.z);
        break;
    case GridFaceAxis::Y:
        area = (upper.x - lower.x) * (upper.z - lower.z);
        break;
    case GridFaceAxis::Z:
        area = (upper.x - lower.x) * (upper.y - lower.y);
        break;
    default:
        throw std::invalid_argument(
            "planar pressure-jump axis is invalid");
    }
    if (!std::isfinite(area) || !(area > 0.0)) {
        throw std::invalid_argument(
            "planar pressure-jump transverse area is invalid");
    }
    return area;
}

std::size_t transverseTileCount(const GridCellCounts counts,
                                const GridFaceAxis axis) {
    std::size_t first = 0;
    std::size_t second = 0;
    switch (axis) {
    case GridFaceAxis::X:
        first = counts.y;
        second = counts.z;
        break;
    case GridFaceAxis::Y:
        first = counts.x;
        second = counts.z;
        break;
    case GridFaceAxis::Z:
        first = counts.x;
        second = counts.y;
        break;
    default:
        throw std::invalid_argument(
            "planar pressure-jump axis is invalid");
    }
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar pressure-jump transverse tile count overflows");
    }
    return first * second;
}

std::vector<PlanarPressureJumpLayerDefinition> canonicalLayers(
    const PeriodicCartesianGrid& grid,
    const std::span<const PlanarPressureJumpLayerDefinition> authored) {
    if (authored.empty()) {
        throw std::invalid_argument(
            "planar pressure-jump layer chain is empty");
    }
    std::vector<PlanarPressureJumpLayerDefinition> result(
        authored.begin(), authored.end());
    std::ranges::sort(
        result,
        [](const auto& first, const auto& second) {
            if (first.physicalPlaneCoordinateMeters
                != second.physicalPlaneCoordinateMeters) {
                return first.physicalPlaneCoordinateMeters
                    < second.physicalPlaneCoordinateMeters;
            }
            return first.surfaceStableId < second.surfaceStableId;
        });

    const GridFaceAxis axis = result.front().topology.axis;
    std::set<std::uint64_t> surfaceIds;
    for (std::size_t index = 0; index < result.size(); ++index) {
        const auto& layer = result[index];
        if (layer.surfaceStableId == 0
            || layer.minusRegionStableId == 0
            || layer.plusRegionStableId == 0) {
            throw std::invalid_argument(
                "planar pressure-jump surface and region IDs must be nonzero");
        }
        if (layer.minusRegionStableId == layer.plusRegionStableId) {
            throw std::invalid_argument(
                "planar pressure-jump layer must separate distinct regions");
        }
        if (layer.topology.axis != axis) {
            throw std::invalid_argument(
                "planar pressure-jump layers must share one axis");
        }
        if (!std::isfinite(layer.physicalPlaneCoordinateMeters)
            || !std::isfinite(layer.pressureJumpPascals)) {
            throw std::invalid_argument(
                "planar pressure-jump layer data must be finite");
        }
        if (!surfaceIds.insert(layer.surfaceStableId).second) {
            throw std::invalid_argument(
                "planar pressure-jump surface IDs must be unique");
        }
        static_cast<void>(movingPlanarCrossingFraction(
            grid, layer.topology,
            layer.physicalPlaneCoordinateMeters));
        if (index != 0) {
            const auto& previous = result[index - 1];
            if (!(previous.physicalPlaneCoordinateMeters
                    < layer.physicalPlaneCoordinateMeters)) {
                throw std::invalid_argument(
                    "planar pressure-jump layers must have distinct positions");
            }
            if (previous.plusRegionStableId
                != layer.minusRegionStableId) {
                throw std::invalid_argument(
                    "planar pressure-jump layers do not form a region chain");
            }
        }
    }
    if (result.back().plusRegionStableId
            != result.front().minusRegionStableId) {
        throw std::invalid_argument(
            "planar pressure-jump region chain is not periodic");
    }
    const double span = result.back().physicalPlaneCoordinateMeters
        - result.front().physicalPlaneCoordinateMeters;
    if (!(span < domainLength(grid, axis))) {
        throw std::invalid_argument(
            "planar pressure-jump layers span more than one period");
    }
    return result;
}

GridFacePressureJump makeCrossing(
    const PlanarPressureJumpLayerDefinition& layer,
    const double crossingFraction,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) {
    return {
        layer.surfaceStableId,
        layer.minusRegionStableId,
        layer.plusRegionStableId,
        layer.topology.axis,
        i, j, k,
        layer.pressureJumpPascals,
        crossingFraction,
    };
}

} // namespace

SharpPressureJumpField makePlanarPressureJumpField(
    const PeriodicCartesianGrid& grid,
    const std::span<const PlanarPressureJumpLayerDefinition> layers) {
    const auto canonical = canonicalLayers(grid, layers);
    const GridCellCounts counts = grid.cellCounts();
    const GridFaceAxis axis = canonical.front().topology.axis;
    const std::size_t tiles = transverseTileCount(counts, axis);
    if (tiles == 0
        || canonical.size()
        > std::numeric_limits<std::size_t>::max() / tiles) {
        throw std::length_error(
            "planar pressure-jump crossing count overflows");
    }
    std::vector<GridFacePressureJump> faces;
    faces.reserve(canonical.size() * tiles);
    for (const auto& layer : canonical) {
        const double fraction = movingPlanarCrossingFraction(
            grid, layer.topology,
            layer.physicalPlaneCoordinateMeters);
        switch (axis) {
        case GridFaceAxis::X:
            for (std::size_t k = 0; k < counts.z; ++k) {
                for (std::size_t j = 0; j < counts.y; ++j) {
                    faces.push_back(makeCrossing(
                        layer, fraction,
                        layer.topology.faceCoordinate, j, k));
                }
            }
            break;
        case GridFaceAxis::Y:
            for (std::size_t k = 0; k < counts.z; ++k) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    faces.push_back(makeCrossing(
                        layer, fraction,
                        i, layer.topology.faceCoordinate, k));
                }
            }
            break;
        case GridFaceAxis::Z:
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    faces.push_back(makeCrossing(
                        layer, fraction,
                        i, j, layer.topology.faceCoordinate));
                }
            }
            break;
        default:
            throw std::invalid_argument(
                "planar pressure-jump axis is invalid");
        }
    }
    return SharpPressureJumpField(grid, std::move(faces));
}

PlanarPressureJumpLayerTranslation translatePlanarPressureJumpLayers(
    const PeriodicCartesianGrid& grid,
    const std::span<const PlanarPressureJumpLayerDefinition> currentLayers,
    const double displacementMeters) {
    if (!std::isfinite(displacementMeters)) {
        throw std::invalid_argument(
            "planar pressure-jump displacement must be finite");
    }
    auto canonical = canonicalLayers(grid, currentLayers);
    static_cast<void>(makePlanarPressureJumpField(grid, canonical));

    PlanarPressureJumpLayerTranslation result;
    result.layers.reserve(canonical.size());
    result.rebaseDirections.reserve(canonical.size());
    for (auto layer : canonical) {
        const double coordinate =
            layer.physicalPlaneCoordinateMeters + displacementMeters;
        if (!std::isfinite(coordinate)) {
            throw std::invalid_argument(
                "planar pressure-jump translated coordinate is non-finite");
        }
        const auto selection = selectMovingPlanarTopology(
            grid, layer.topology, coordinate);
        layer.topology = selection.topology;
        layer.physicalPlaneCoordinateMeters = coordinate;
        result.layers.push_back(layer);
        result.rebaseDirections.push_back(
            selection.rebaseDirection);
    }
    static_cast<void>(makePlanarPressureJumpField(
        grid, result.layers));
    return result;
}

StaticPlanarPressureRegionProfile makeStaticPlanarPressureRegionProfile(
    const PeriodicCartesianGrid& grid,
    const std::span<const PlanarPressureJumpLayerDefinition> layers,
    const double volumeMeanPressurePascals) {
    if (!std::isfinite(volumeMeanPressurePascals)) {
        throw std::invalid_argument(
            "planar pressure region mean pressure must be finite");
    }
    const auto canonical = canonicalLayers(grid, layers);

    StaticPlanarPressureRegionProfile result;
    result.axis = canonical.front().topology.axis;
    result.windowLowerCoordinateMeters =
        canonical.front().physicalPlaneCoordinateMeters;
    const double length = domainLength(grid, result.axis);
    result.windowUpperCoordinateMeters =
        result.windowLowerCoordinateMeters + length;
    const double area = transverseArea(grid, result.axis);
    result.geometricDomainVolumeCubicMeters = area * length;
    result.requestedVolumeMeanPressurePascals =
        volumeMeanPressurePascals;
    if (!std::isfinite(result.windowUpperCoordinateMeters)
        || !std::isfinite(result.geometricDomainVolumeCubicMeters)
        || !(result.geometricDomainVolumeCubicMeters > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region domain geometry is invalid");
    }

    std::map<std::uint64_t, double> pressureOffsets;
    std::map<std::uint64_t, double> regionVolumes;
    pressureOffsets.emplace(
        canonical.front().minusRegionStableId, 0.0);
    result.intervals.reserve(canonical.size());
    double cumulativeJump = 0.0;
    double totalVolume = 0.0;
    double volumeWeightedOffset = 0.0;
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto& layer = canonical[index];
        cumulativeJump += layer.pressureJumpPascals;
        if (!std::isfinite(cumulativeJump)) {
            throw std::invalid_argument(
                "planar pressure region cumulative jump is non-finite");
        }
        if (index + 1 == canonical.size() && cumulativeJump != 0.0) {
            throw std::invalid_argument(
                "planar pressure region jump cycle does not close exactly");
        }
        const auto [offset, inserted] = pressureOffsets.emplace(
            layer.plusRegionStableId, cumulativeJump);
        if (!inserted && offset->second != cumulativeJump) {
            throw std::invalid_argument(
                "planar pressure region has inconsistent repeated pressure potential");
        }

        const auto& next = canonical[(index + 1) % canonical.size()];
        const double upperCoordinate = index + 1 < canonical.size()
            ? next.physicalPlaneCoordinateMeters
            : result.windowUpperCoordinateMeters;
        const double thickness = upperCoordinate
            - layer.physicalPlaneCoordinateMeters;
        const double volume = area * thickness;
        if (!std::isfinite(thickness) || !(thickness > 0.0)
            || !std::isfinite(volume) || !(volume > 0.0)) {
            throw std::invalid_argument(
                "planar pressure region interval geometry is invalid");
        }
        totalVolume += volume;
        volumeWeightedOffset += volume * cumulativeJump;
        regionVolumes[layer.plusRegionStableId] += volume;
        if (!std::isfinite(totalVolume)
            || !std::isfinite(volumeWeightedOffset)
            || !std::isfinite(regionVolumes[layer.plusRegionStableId])) {
            throw std::invalid_argument(
                "planar pressure region accumulation is non-finite");
        }
        result.intervals.push_back({
            layer.surfaceStableId,
            next.surfaceStableId,
            layer.plusRegionStableId,
            layer.physicalPlaneCoordinateMeters,
            upperCoordinate,
            volume,
            cumulativeJump,
        });
    }
    if (!(totalVolume > 0.0)) {
        throw std::invalid_argument(
            "planar pressure region interval volume is empty");
    }

    const double gaugeShift = volumeMeanPressurePascals
        - volumeWeightedOffset / totalVolume;
    if (!std::isfinite(gaugeShift)) {
        throw std::invalid_argument(
            "planar pressure region gauge shift is non-finite");
    }
    double achievedWeightedPressure = 0.0;
    for (auto& interval : result.intervals) {
        interval.pressurePascals += gaugeShift;
        achievedWeightedPressure += interval.volumeCubicMeters
            * interval.pressurePascals;
        if (!std::isfinite(interval.pressurePascals)
            || !std::isfinite(achievedWeightedPressure)) {
            throw std::invalid_argument(
                "planar pressure region pressure accumulation is non-finite");
        }
    }
    result.regions.reserve(pressureOffsets.size());
    for (const auto& [regionStableId, pressureOffset] : pressureOffsets) {
        const auto volume = regionVolumes.find(regionStableId);
        if (volume == regionVolumes.end()
            || !(volume->second > 0.0)) {
            throw std::invalid_argument(
                "planar pressure region has no positive interval volume");
        }
        result.regions.push_back({
            regionStableId,
            volume->second,
            pressureOffset + gaugeShift,
        });
    }
    result.intervalVolumeCubicMeters = totalVolume;
    result.volumeClosureResidualCubicMeters = totalVolume
        - result.geometricDomainVolumeCubicMeters;
    result.achievedVolumeMeanPressurePascals =
        achievedWeightedPressure / totalVolume;
    if (!std::isfinite(result.volumeClosureResidualCubicMeters)
        || !std::isfinite(result.achievedVolumeMeanPressurePascals)) {
        throw std::invalid_argument(
            "planar pressure region diagnostics are non-finite");
    }
    return result;
}

} // namespace simwing::fsi::fluid
