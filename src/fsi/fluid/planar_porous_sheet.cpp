#include "fluid/planar_porous_sheet.h"

#include <cmath>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

void validateDefinition(const PlanarPorousSheetDefinition& definition) {
    if (definition.surfaceStableId == 0
        || definition.minusRegionStableId == 0
        || definition.plusRegionStableId == 0) {
        throw std::invalid_argument(
            "planar porous sheet surface and region IDs must be nonzero");
    }
    if (definition.minusRegionStableId
            == definition.plusRegionStableId) {
        throw std::invalid_argument(
            "planar porous sheet must separate two distinct regions");
    }
    if (!std::isfinite(
            definition.surfaceNormalVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "planar porous sheet normal velocity must be finite");
    }
    const auto& resistance = definition.resistance;
    if (!std::isfinite(resistance.linearPascalSecondsPerMeter)
        || !std::isfinite(
            resistance.quadraticPascalSecondsSquaredPerSquareMeter)
        || resistance.linearPascalSecondsPerMeter < 0.0
        || resistance.quadraticPascalSecondsSquaredPerSquareMeter < 0.0
        || (resistance.linearPascalSecondsPerMeter == 0.0
            && resistance.quadraticPascalSecondsSquaredPerSquareMeter
                == 0.0)) {
        throw std::invalid_argument(
            "planar porous sheet requires a finite nonnegative active resistance");
    }
}

PorousGridFaceCrossing makeCrossing(
    const PlanarPorousSheetDefinition& definition,
    const double crossingFraction,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) {
    return {
        definition.surfaceStableId,
        definition.minusRegionStableId,
        definition.plusRegionStableId,
        definition.topology.axis,
        i,
        j,
        k,
        crossingFraction,
        definition.surfaceNormalVelocityMetersPerSecond,
        definition.resistance,
    };
}

} // namespace

std::vector<PorousGridFaceCrossing> makePlanarPorousSheetCrossings(
    const PeriodicCartesianGrid& grid,
    const PlanarPorousSheetDefinition& definition) {
    validateDefinition(definition);
    const double crossingFraction = movingPorousCrossingFraction(
        grid, definition.topology,
        definition.physicalPlaneCoordinateMeters);
    const GridCellCounts counts = grid.cellCounts();
    std::vector<PorousGridFaceCrossing> result;
    switch (definition.topology.axis) {
    case GridFaceAxis::X:
        result.reserve(counts.y * counts.z);
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                result.push_back(makeCrossing(
                    definition, crossingFraction,
                    definition.topology.faceCoordinate, j, k));
            }
        }
        break;
    case GridFaceAxis::Y:
        result.reserve(counts.x * counts.z);
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                result.push_back(makeCrossing(
                    definition, crossingFraction,
                    i, definition.topology.faceCoordinate, k));
            }
        }
        break;
    case GridFaceAxis::Z:
        result.reserve(counts.x * counts.y);
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                result.push_back(makeCrossing(
                    definition, crossingFraction,
                    i, j, definition.topology.faceCoordinate));
            }
        }
        break;
    default:
        throw std::invalid_argument(
            "planar porous sheet topology axis is invalid");
    }
    return result;
}

} // namespace simwing::fsi::fluid
