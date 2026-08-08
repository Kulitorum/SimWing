#include "fluid/moving_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 add(const Vector3& first, const Vector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

Vector3 subtract(const Vector3& first, const Vector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

double length(const Vector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

Vector3 scale(const Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

std::size_t axisOrdinal(const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return 0;
    case GridFaceAxis::Y:
        return 1;
    case GridFaceAxis::Z:
        return 2;
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

void validateFaceCoordinates(const PeriodicCartesianGrid& grid,
                             const GridFaceMovingInterface& face) {
    const auto counts = grid.cellCounts();
    static_cast<void>(axisOrdinal(face.axis));
    if (face.i >= counts.x || face.j >= counts.y || face.k >= counts.z) {
        throw std::invalid_argument(
            "moving interface face coordinate is out of range");
    }
}

std::size_t faceIndex(const PeriodicCartesianGrid& grid,
                      const GridFaceMovingInterface& face) {
    return grid.cellIndex(face.i, face.j, face.k);
}

double normalFaceVelocity(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocity,
    const GridFaceMovingInterface& face) {
    const std::size_t index = faceIndex(grid, face);
    switch (face.axis) {
    case GridFaceAxis::X:
        return velocity.xFaces()[index];
    case GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

std::pair<std::size_t, std::size_t> adjacentCells(
    const PeriodicCartesianGrid& grid,
    const GridFaceMovingInterface& face) {
    const auto counts = grid.cellCounts();
    const std::size_t plus = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case GridFaceAxis::X:
        return {
            grid.cellIndex(
                (face.i + counts.x - 1) % counts.x, face.j, face.k),
            plus};
    case GridFaceAxis::Y:
        return {
            grid.cellIndex(
                face.i, (face.j + counts.y - 1) % counts.y, face.k),
            plus};
    case GridFaceAxis::Z:
        return {
            grid.cellIndex(
                face.i, face.j, (face.k + counts.z - 1) % counts.z),
            plus};
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

double dot(const CellScalarField& first, const CellScalarField& second) {
    const auto firstValues = first.values();
    const auto secondValues = second.values();
    double result = 0.0;
    for (std::size_t index = 0; index < firstValues.size(); ++index) {
        result += firstValues[index] * secondValues[index];
    }
    return result;
}

void validateSettings(const MovingInterfaceProjectionSettings& settings) {
    const auto& projection = settings.projection;
    if (!std::isfinite(projection.densityKgPerCubicMeter)
        || !(projection.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(projection.timeStepSeconds)
        || !(projection.timeStepSeconds > 0.0)) {
        throw std::invalid_argument(
            "moving-interface projection density and time step must be positive");
    }
    if (!std::isfinite(projection.absoluteResidualTolerance)
        || projection.absoluteResidualTolerance < 0.0
        || !std::isfinite(projection.relativeResidualTolerance)
        || projection.relativeResidualTolerance < 0.0
        || (projection.absoluteResidualTolerance == 0.0
            && projection.relativeResidualTolerance == 0.0)) {
        throw std::invalid_argument(
            "moving-interface projection residual tolerances are invalid");
    }
    if (!std::isfinite(
            settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond)
        || settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond < 0.0) {
        throw std::invalid_argument(
            "moving-interface region volume-rate tolerance is invalid");
    }
}

std::vector<std::size_t> cellRegionIndices(
    const FaceAlignedMovingInterface& interfaces) {
    const auto stableIds = interfaces.regionStableIds();
    std::vector<std::size_t> result;
    result.reserve(interfaces.cellRegionStableIds().size());
    for (const auto stableId : interfaces.cellRegionStableIds()) {
        const auto found = std::lower_bound(
            stableIds.begin(), stableIds.end(), stableId);
        if (found == stableIds.end() || *found != stableId) {
            throw std::logic_error(
                "moving-interface cell region is absent from canonical regions");
        }
        result.push_back(static_cast<std::size_t>(found - stableIds.begin()));
    }
    return result;
}

std::vector<double> regionMeans(
    const CellScalarField& field,
    const std::span<const std::size_t> cellRegions,
    const std::size_t regionCount,
    std::vector<std::size_t>* counts = nullptr) {
    std::vector<double> sums(regionCount, 0.0);
    std::vector<std::size_t> localCounts(regionCount, 0);
    for (std::size_t cell = 0; cell < cellRegions.size(); ++cell) {
        sums[cellRegions[cell]] += field.values()[cell];
        ++localCounts[cellRegions[cell]];
    }
    for (std::size_t region = 0; region < regionCount; ++region) {
        if (localCounts[region] == 0) {
            throw std::logic_error("moving-interface fluid region has no cells");
        }
        sums[region] /= static_cast<double>(localCounts[region]);
    }
    if (counts != nullptr) {
        *counts = std::move(localCounts);
    }
    return sums;
}

void subtractRegionMeans(
    CellScalarField& field,
    const std::span<const std::size_t> cellRegions,
    const std::size_t regionCount) {
    const auto means = regionMeans(field, cellRegions, regionCount);
    for (std::size_t cell = 0; cell < cellRegions.size(); ++cell) {
        field.values()[cell] -= means[cellRegions[cell]];
    }
}

void applyDisconnectedNegativeLaplacian(
    const PeriodicCartesianGrid& grid,
    const FaceAlignedMovingInterface& interfaces,
    const CellScalarField& input,
    CellScalarField& output) {
    std::ranges::fill(output.values(), 0.0);
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    const auto xConstraints = interfaces.xFaceConstraints();
    const auto yConstraints = interfaces.yFaceConstraints();
    const auto zConstraints = interfaces.zFaceConstraints();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto plus = grid.cellIndex(i, j, k);
                const auto xMinus = grid.cellIndex(
                    (i + counts.x - 1) % counts.x, j, k);
                const auto yMinus = grid.cellIndex(
                    i, (j + counts.y - 1) % counts.y, k);
                const auto zMinus = grid.cellIndex(
                    i, j, (k + counts.z - 1) % counts.z);
                if (xConstraints[plus] == 0) {
                    const double contribution =
                        (input.values()[plus] - input.values()[xMinus])
                        / (spacing.x * spacing.x);
                    output.values()[plus] += contribution;
                    output.values()[xMinus] -= contribution;
                }
                if (yConstraints[plus] == 0) {
                    const double contribution =
                        (input.values()[plus] - input.values()[yMinus])
                        / (spacing.y * spacing.y);
                    output.values()[plus] += contribution;
                    output.values()[yMinus] -= contribution;
                }
                if (zConstraints[plus] == 0) {
                    const double contribution =
                        (input.values()[plus] - input.values()[zMinus])
                        / (spacing.z * spacing.z);
                    output.values()[plus] += contribution;
                    output.values()[zMinus] -= contribution;
                }
            }
        }
    }
}

void imposeInterfaceVelocity(
    MacVelocityField& velocity,
    const FaceAlignedMovingInterface& interfaces) {
    const auto xConstraints = interfaces.xFaceConstraints();
    const auto yConstraints = interfaces.yFaceConstraints();
    const auto zConstraints = interfaces.zFaceConstraints();
    const auto xVelocities =
        interfaces.xFaceNormalVelocitiesMetersPerSecond();
    const auto yVelocities =
        interfaces.yFaceNormalVelocitiesMetersPerSecond();
    const auto zVelocities =
        interfaces.zFaceNormalVelocitiesMetersPerSecond();
    for (std::size_t face = 0; face < xConstraints.size(); ++face) {
        if (xConstraints[face] != 0) {
            velocity.xFaces()[face] = xVelocities[face];
        }
        if (yConstraints[face] != 0) {
            velocity.yFaces()[face] = yVelocities[face];
        }
        if (zConstraints[face] != 0) {
            velocity.zFaces()[face] = zVelocities[face];
        }
    }
}

void applyPressureCorrection(
    const PeriodicCartesianGrid& grid,
    const FaceAlignedMovingInterface& interfaces,
    const CellScalarField& pressureCorrection,
    const double correctionScale,
    MacVelocityField& velocity) {
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    const auto xConstraints = interfaces.xFaceConstraints();
    const auto yConstraints = interfaces.yFaceConstraints();
    const auto zConstraints = interfaces.zFaceConstraints();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto plus = grid.cellIndex(i, j, k);
                if (xConstraints[plus] == 0) {
                    const auto minus = grid.cellIndex(
                        (i + counts.x - 1) % counts.x, j, k);
                    velocity.xFaces()[plus] -= correctionScale
                        * (pressureCorrection.values()[plus]
                           - pressureCorrection.values()[minus])
                        / spacing.x;
                }
                if (yConstraints[plus] == 0) {
                    const auto minus = grid.cellIndex(
                        i, (j + counts.y - 1) % counts.y, k);
                    velocity.yFaces()[plus] -= correctionScale
                        * (pressureCorrection.values()[plus]
                           - pressureCorrection.values()[minus])
                        / spacing.y;
                }
                if (zConstraints[plus] == 0) {
                    const auto minus = grid.cellIndex(
                        i, j, (k + counts.z - 1) % counts.z);
                    velocity.zFaces()[plus] -= correctionScale
                        * (pressureCorrection.values()[plus]
                           - pressureCorrection.values()[minus])
                        / spacing.z;
                }
            }
        }
    }
}

double faceArea(const PeriodicCartesianGrid& grid,
                const GridFaceAxis axis) {
    const auto spacing = grid.cellSpacingMeters();
    switch (axis) {
    case GridFaceAxis::X:
        return spacing.y * spacing.z;
    case GridFaceAxis::Y:
        return spacing.x * spacing.z;
    case GridFaceAxis::Z:
        return spacing.x * spacing.y;
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

std::pair<Vector3, Vector3> faceCorners(
    const PeriodicCartesianGrid& grid,
    const GridFaceMovingInterface& face) {
    const Vector3 spacing = grid.cellSpacingMeters();
    Vector3 center;
    switch (face.axis) {
    case GridFaceAxis::X:
        center = grid.xFaceCenterMeters(face.i, face.j, face.k);
        return {
            {center.x, center.y - 0.5 * spacing.y,
             center.z - 0.5 * spacing.z},
            {center.x, center.y + 0.5 * spacing.y,
             center.z + 0.5 * spacing.z},
        };
    case GridFaceAxis::Y:
        center = grid.yFaceCenterMeters(face.i, face.j, face.k);
        return {
            {center.x - 0.5 * spacing.x, center.y,
             center.z - 0.5 * spacing.z},
            {center.x + 0.5 * spacing.x, center.y,
             center.z + 0.5 * spacing.z},
        };
    case GridFaceAxis::Z:
        center = grid.zFaceCenterMeters(face.i, face.j, face.k);
        return {
            {center.x - 0.5 * spacing.x,
             center.y - 0.5 * spacing.y, center.z},
            {center.x + 0.5 * spacing.x,
             center.y + 0.5 * spacing.y, center.z},
        };
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

Vector3 axialVector(const GridFaceAxis axis, const double value) {
    switch (axis) {
    case GridFaceAxis::X:
        return {value, 0.0, 0.0};
    case GridFaceAxis::Y:
        return {0.0, value, 0.0};
    case GridFaceAxis::Z:
        return {0.0, 0.0, value};
    }
    throw std::invalid_argument("moving interface face has an unknown axis");
}

double constrainedVelocityError(
    const MacVelocityField& velocity,
    const FaceAlignedMovingInterface& interfaces) {
    double result = 0.0;
    for (const auto& face : interfaces.faces()) {
        const std::size_t index = face.i
            + interfaces.cellCounts().x
                * (face.j + interfaces.cellCounts().y * face.k);
        double actual = 0.0;
        switch (face.axis) {
        case GridFaceAxis::X:
            actual = velocity.xFaces()[index];
            break;
        case GridFaceAxis::Y:
            actual = velocity.yFaces()[index];
            break;
        case GridFaceAxis::Z:
            actual = velocity.zFaces()[index];
            break;
        }
        result = std::max(
            result,
            std::abs(actual - face.normalVelocityMetersPerSecond));
    }
    return result;
}

} // namespace

FaceAlignedMovingInterface::FaceAlignedMovingInterface(
    const PeriodicCartesianGrid& grid,
    std::vector<GridFaceMovingInterface> faces)
    : cellCounts_(grid.cellCounts()), lowerMeters_(grid.lowerMeters()),
      upperMeters_(grid.upperMeters()), faces_(std::move(faces)),
      xFaceConstraints_(grid.cellCount(), 0),
      yFaceConstraints_(grid.cellCount(), 0),
      zFaceConstraints_(grid.cellCount(), 0),
      xFaceNormalVelocitiesMetersPerSecond_(grid.cellCount(), 0.0),
      yFaceNormalVelocitiesMetersPerSecond_(grid.cellCount(), 0.0),
      zFaceNormalVelocitiesMetersPerSecond_(grid.cellCount(), 0.0) {
    if (faces_.empty()) {
        throw std::invalid_argument(
            "a face-aligned moving interface requires at least one face");
    }
    for (const auto& face : faces_) {
        validateFaceCoordinates(grid, face);
    }
    std::sort(faces_.begin(), faces_.end(), [&](const auto& first, const auto& second) {
        return std::tuple(axisOrdinal(first.axis), faceIndex(grid, first))
            < std::tuple(axisOrdinal(second.axis), faceIndex(grid, second));
    });

    std::map<std::uint64_t, std::pair<std::uint64_t, std::uint64_t>>
        surfaceRegionPairs;
    std::tuple<std::size_t, std::size_t> previousFace{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max()};
    for (const auto& face : faces_) {
        if (face.surfaceStableId == 0
            || face.minusRegionStableId == 0
            || face.plusRegionStableId == 0
            || !std::isfinite(face.normalVelocityMetersPerSecond)) {
            throw std::invalid_argument(
                "moving interface IDs, regions, and velocity must be valid");
        }
        const auto key = std::tuple(
            axisOrdinal(face.axis), faceIndex(grid, face));
        if (key == previousFace) {
            throw std::invalid_argument(
                "multiple moving interfaces on one grid face are unsupported");
        }
        previousFace = key;
        const auto pair = std::pair(
            face.minusRegionStableId, face.plusRegionStableId);
        const auto [found, inserted] = surfaceRegionPairs.emplace(
            face.surfaceStableId, pair);
        if (!inserted && found->second != pair) {
            throw std::invalid_argument(
                "one moving surface ID cannot alias different oriented regions");
        }
        const std::size_t index = faceIndex(grid, face);
        switch (face.axis) {
        case GridFaceAxis::X:
            xFaceConstraints_[index] = 1;
            xFaceNormalVelocitiesMetersPerSecond_[index] =
                face.normalVelocityMetersPerSecond;
            break;
        case GridFaceAxis::Y:
            yFaceConstraints_[index] = 1;
            yFaceNormalVelocitiesMetersPerSecond_[index] =
                face.normalVelocityMetersPerSecond;
            break;
        case GridFaceAxis::Z:
            zFaceConstraints_[index] = 1;
            zFaceNormalVelocitiesMetersPerSecond_[index] =
                face.normalVelocityMetersPerSecond;
            break;
        }
    }

    const auto counts = grid.cellCounts();
    constexpr std::size_t unassigned = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> components(grid.cellCount(), unassigned);
    std::vector<std::size_t> queue;
    std::size_t componentCount = 0;
    for (std::size_t seed = 0; seed < grid.cellCount(); ++seed) {
        if (components[seed] != unassigned) {
            continue;
        }
        queue.clear();
        queue.push_back(seed);
        components[seed] = componentCount;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
            const std::size_t cell = queue[cursor];
            const std::size_t i = cell % counts.x;
            const std::size_t j = (cell / counts.x) % counts.y;
            const std::size_t k = cell / (counts.x * counts.y);
            const auto visit = [&](const std::size_t neighbour,
                                   const bool constrained) {
                if (!constrained && components[neighbour] == unassigned) {
                    components[neighbour] = componentCount;
                    queue.push_back(neighbour);
                }
            };
            const auto xPlus = grid.cellIndex((i + 1) % counts.x, j, k);
            const auto xMinus = grid.cellIndex(
                (i + counts.x - 1) % counts.x, j, k);
            const auto yPlus = grid.cellIndex(i, (j + 1) % counts.y, k);
            const auto yMinus = grid.cellIndex(
                i, (j + counts.y - 1) % counts.y, k);
            const auto zPlus = grid.cellIndex(i, j, (k + 1) % counts.z);
            const auto zMinus = grid.cellIndex(
                i, j, (k + counts.z - 1) % counts.z);
            visit(xPlus, xFaceConstraints_[xPlus] != 0);
            visit(xMinus, xFaceConstraints_[cell] != 0);
            visit(yPlus, yFaceConstraints_[yPlus] != 0);
            visit(yMinus, yFaceConstraints_[cell] != 0);
            visit(zPlus, zFaceConstraints_[zPlus] != 0);
            visit(zMinus, zFaceConstraints_[cell] != 0);
        }
        ++componentCount;
    }

    std::vector<std::uint64_t> componentRegions(componentCount, 0);
    const auto bindRegion = [&](const std::size_t component,
                                const std::uint64_t stableId) {
        if (componentRegions[component] == 0) {
            componentRegions[component] = stableId;
        } else if (componentRegions[component] != stableId) {
            throw std::invalid_argument(
                "moving interface region labels conflict with grid connectivity");
        }
    };
    for (const auto& face : faces_) {
        const auto [minusCell, plusCell] = adjacentCells(grid, face);
        const auto minusComponent = components[minusCell];
        const auto plusComponent = components[plusCell];
        const bool sameRegion = face.minusRegionStableId
            == face.plusRegionStableId;
        if (sameRegion != (minusComponent == plusComponent)) {
            throw std::invalid_argument(
                "moving interface side labels disagree with grid connectivity");
        }
        bindRegion(minusComponent, face.minusRegionStableId);
        bindRegion(plusComponent, face.plusRegionStableId);
    }
    if (std::ranges::find(componentRegions, std::uint64_t{0})
        != componentRegions.end()) {
        throw std::invalid_argument(
            "every moving-interface fluid component requires a stable region ID");
    }
    std::map<std::uint64_t, std::size_t> regionComponents;
    for (std::size_t component = 0;
         component < componentRegions.size(); ++component) {
        if (!regionComponents.emplace(
                componentRegions[component], component).second) {
            throw std::invalid_argument(
                "one moving fluid region cannot alias disconnected components");
        }
    }
    for (const auto& [stableId, component] : regionComponents) {
        static_cast<void>(component);
        regionStableIds_.push_back(stableId);
    }
    cellRegionStableIds_.reserve(grid.cellCount());
    for (const auto component : components) {
        cellRegionStableIds_.push_back(componentRegions[component]);
    }
}

std::uint32_t FaceAlignedMovingInterface::version() const noexcept {
    return version_;
}

GridCellCounts FaceAlignedMovingInterface::cellCounts() const noexcept {
    return cellCounts_;
}

bool FaceAlignedMovingInterface::matches(
    const PeriodicCartesianGrid& grid) const noexcept {
    return cellCounts_ == grid.cellCounts()
        && lowerMeters_ == grid.lowerMeters()
        && upperMeters_ == grid.upperMeters();
}

std::size_t FaceAlignedMovingInterface::faceCount() const noexcept {
    return faces_.size();
}

std::size_t FaceAlignedMovingInterface::regionCount() const noexcept {
    return regionStableIds_.size();
}

std::span<const GridFaceMovingInterface>
FaceAlignedMovingInterface::faces() const noexcept {
    return faces_;
}

std::span<const std::uint64_t>
FaceAlignedMovingInterface::regionStableIds() const noexcept {
    return regionStableIds_;
}

std::span<const std::uint64_t>
FaceAlignedMovingInterface::cellRegionStableIds() const noexcept {
    return cellRegionStableIds_;
}

std::span<const std::uint8_t>
FaceAlignedMovingInterface::xFaceConstraints() const noexcept {
    return xFaceConstraints_;
}

std::span<const std::uint8_t>
FaceAlignedMovingInterface::yFaceConstraints() const noexcept {
    return yFaceConstraints_;
}

std::span<const std::uint8_t>
FaceAlignedMovingInterface::zFaceConstraints() const noexcept {
    return zFaceConstraints_;
}

std::span<const double>
FaceAlignedMovingInterface::xFaceNormalVelocitiesMetersPerSecond() const noexcept {
    return xFaceNormalVelocitiesMetersPerSecond_;
}

std::span<const double>
FaceAlignedMovingInterface::yFaceNormalVelocitiesMetersPerSecond() const noexcept {
    return yFaceNormalVelocitiesMetersPerSecond_;
}

std::span<const double>
FaceAlignedMovingInterface::zFaceNormalVelocitiesMetersPerSecond() const noexcept {
    return zFaceNormalVelocitiesMetersPerSecond_;
}

MovingInterfaceProjectionDiagnostics projectVelocityWithMovingInterfaces(
    const PeriodicCartesianGrid& grid,
    MacVelocityField& predictedVelocityMetersPerSecond,
    CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionSettings& settings) {
    validateSettings(settings);
    if (!interfaces.matches(grid)
        || !predictedVelocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)) {
        throw std::invalid_argument(
            "moving-interface projection inputs do not match their grid");
    }
    if (!isFinite(predictedVelocityMetersPerSecond)
        || !isFinite(pressurePascals)) {
        throw std::invalid_argument(
            "moving-interface projection fields must be finite");
    }

    MovingInterfaceProjectionDiagnostics diagnostics;
    diagnostics.interfaceVersion = interfaces.version();
    diagnostics.interfaceFaceCount = interfaces.faceCount();
    diagnostics.fluidRegionCount = interfaces.regionCount();
    const auto regionIds = interfaces.regionStableIds();
    const auto cellRegions = cellRegionIndices(interfaces);
    std::vector<std::size_t> regionCellCounts;
    const auto pressureMeansBefore = regionMeans(
        pressurePascals, cellRegions, regionIds.size(), &regionCellCounts);
    diagnostics.regions.reserve(regionIds.size());
    for (std::size_t region = 0; region < regionIds.size(); ++region) {
        diagnostics.regions.push_back({
            regionIds[region], regionCellCounts[region], 0.0,
            pressureMeansBefore[region], pressureMeansBefore[region]});
    }

    MacVelocityField candidateVelocity = predictedVelocityMetersPerSecond;
    imposeInterfaceVelocity(candidateVelocity, interfaces);
    CellScalarField divergence(grid);
    computeDivergence(grid, candidateVelocity, divergence);
    auto& projection = diagnostics.projection;
    projection.compatibilityDivergencePerSecond = mean(divergence);
    projection.divergenceL2BeforePerSecond = l2Norm(divergence);
    projection.divergenceMaximumBeforePerSecond =
        maximumAbsoluteValue(divergence);
    projection.kineticEnergyBeforeJoules = kineticEnergyJoules(
        grid, candidateVelocity,
        settings.projection.densityKgPerCubicMeter);
    const double cellVolume = grid.cellVolumeCubicMeters();
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        diagnostics.regions[cellRegions[cell]]
            .compatibilityVolumeRateCubicMetersPerSecond +=
                divergence.values()[cell] * cellVolume;
    }
    for (const auto& region : diagnostics.regions) {
        diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond =
            std::max(
                diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond,
                std::abs(region.compatibilityVolumeRateCubicMetersPerSecond));
    }
    if (diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond
        > settings.absoluteRegionVolumeRateToleranceCubicMetersPerSecond) {
        projection.divergenceL2AfterPerSecond =
            projection.divergenceL2BeforePerSecond;
        projection.divergenceMaximumAfterPerSecond =
            projection.divergenceMaximumBeforePerSecond;
        projection.kineticEnergyAfterJoules =
            projection.kineticEnergyBeforeJoules;
        projection.pressureMeanPascals = mean(pressurePascals);
        return diagnostics;
    }

    CellScalarField rightHandSide(grid);
    const double rightHandSideScale =
        -settings.projection.densityKgPerCubicMeter
        / settings.projection.timeStepSeconds;
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        rightHandSide.values()[cell] =
            rightHandSideScale * divergence.values()[cell];
    }
    // Remove only compatibility roundoff already admitted by the explicit
    // physical volume-rate tolerance. Material incompatibility returned above.
    subtractRegionMeans(rightHandSide, cellRegions, regionIds.size());

    CellScalarField candidateCorrection = pressurePascals;
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        candidateCorrection.values()[cell] -=
            pressureMeansBefore[cellRegions[cell]];
    }
    CellScalarField operatorPressure(grid);
    applyDisconnectedNegativeLaplacian(
        grid, interfaces, candidateCorrection, operatorPressure);
    CellScalarField residual(grid);
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        residual.values()[cell] = rightHandSide.values()[cell]
            - operatorPressure.values()[cell];
    }
    subtractRegionMeans(residual, cellRegions, regionIds.size());
    projection.initialResidualPascalsPerSquareMeter = l2Norm(residual);
    projection.finalResidualPascalsPerSquareMeter =
        projection.initialResidualPascalsPerSquareMeter;
    const double convergenceThreshold = std::max(
        settings.projection.absoluteResidualTolerance,
        settings.projection.relativeResidualTolerance
            * projection.initialResidualPascalsPerSquareMeter);

    CellScalarField direction = residual;
    CellScalarField operatorDirection(grid);
    double residualSquared = dot(residual, residual);
    projection.converged =
        projection.finalResidualPascalsPerSquareMeter <= convergenceThreshold;
    while (!projection.converged
           && projection.iterationCount
               < settings.projection.maximumIterations) {
        applyDisconnectedNegativeLaplacian(
            grid, interfaces, direction, operatorDirection);
        const double denominator = dot(direction, operatorDirection);
        if (!std::isfinite(denominator) || !(denominator > 0.0)) {
            break;
        }
        const double alpha = residualSquared / denominator;
        for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
            candidateCorrection.values()[cell] +=
                alpha * direction.values()[cell];
            residual.values()[cell] -=
                alpha * operatorDirection.values()[cell];
        }
        ++projection.iterationCount;
        subtractRegionMeans(
            candidateCorrection, cellRegions, regionIds.size());
        subtractRegionMeans(residual, cellRegions, regionIds.size());
        const double nextResidualSquared = dot(residual, residual);
        projection.finalResidualPascalsPerSquareMeter = std::sqrt(
            nextResidualSquared / static_cast<double>(grid.cellCount()));
        projection.converged = std::isfinite(
            projection.finalResidualPascalsPerSquareMeter)
            && projection.finalResidualPascalsPerSquareMeter
                <= convergenceThreshold;
        if (projection.converged) {
            residualSquared = nextResidualSquared;
            break;
        }
        if (!std::isfinite(nextResidualSquared) || !(residualSquared > 0.0)) {
            break;
        }
        const double beta = nextResidualSquared / residualSquared;
        for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
            direction.values()[cell] = residual.values()[cell]
                + beta * direction.values()[cell];
        }
        subtractRegionMeans(direction, cellRegions, regionIds.size());
        residualSquared = nextResidualSquared;
    }

    if (projection.converged) {
        applyDisconnectedNegativeLaplacian(
            grid, interfaces, candidateCorrection, operatorPressure);
        for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
            residual.values()[cell] = rightHandSide.values()[cell]
                - operatorPressure.values()[cell];
        }
        subtractRegionMeans(residual, cellRegions, regionIds.size());
        projection.finalResidualPascalsPerSquareMeter = l2Norm(residual);
        projection.converged = std::isfinite(
            projection.finalResidualPascalsPerSquareMeter)
            && projection.finalResidualPascalsPerSquareMeter
                <= convergenceThreshold;
    }
    if (!projection.converged) {
        projection.divergenceL2AfterPerSecond =
            projection.divergenceL2BeforePerSecond;
        projection.divergenceMaximumAfterPerSecond =
            projection.divergenceMaximumBeforePerSecond;
        projection.kineticEnergyAfterJoules =
            projection.kineticEnergyBeforeJoules;
        projection.pressureMeanPascals = mean(pressurePascals);
        return diagnostics;
    }

    applyPressureCorrection(
        grid, interfaces, candidateCorrection,
        settings.projection.timeStepSeconds
            / settings.projection.densityKgPerCubicMeter,
        candidateVelocity);
    CellScalarField candidatePressure = candidateCorrection;
    for (std::size_t cell = 0; cell < grid.cellCount(); ++cell) {
        candidatePressure.values()[cell] +=
            pressureMeansBefore[cellRegions[cell]];
    }
    CellScalarField correctedDivergence(grid);
    computeDivergence(grid, candidateVelocity, correctedDivergence);
    projection.divergenceL2AfterPerSecond = l2Norm(correctedDivergence);
    projection.divergenceMaximumAfterPerSecond =
        maximumAbsoluteValue(correctedDivergence);
    projection.kineticEnergyAfterJoules = kineticEnergyJoules(
        grid, candidateVelocity,
        settings.projection.densityKgPerCubicMeter);
    projection.pressureMeanPascals = mean(candidatePressure);
    diagnostics.maximumNormalVelocityErrorMetersPerSecond =
        constrainedVelocityError(candidateVelocity, interfaces);
    const auto pressureMeansAfter = regionMeans(
        candidatePressure, cellRegions, regionIds.size());
    for (std::size_t region = 0; region < regionIds.size(); ++region) {
        diagnostics.regions[region].pressureMeanAfterPascals =
            pressureMeansAfter[region];
    }

    std::map<std::uint64_t, MovingInterfaceSurfaceDiagnostics> surfaces;
    diagnostics.faces.reserve(interfaces.faceCount());
    for (const auto& face : interfaces.faces()) {
        const auto [minusCell, plusCell] = adjacentCells(grid, face);
        const double area = faceArea(grid, face.axis);
        const double pressureDifferencePascals =
            (candidatePressure.values()[minusCell]
             - candidatePressure.values()[plusCell]);
        const double pressureForceAlongPositiveAxis =
            pressureDifferencePascals * area;
        const Vector3 traction = axialVector(
            face.axis, pressureDifferencePascals);
        const Vector3 force = axialVector(
            face.axis, pressureForceAlongPositiveAxis);
        const double directConstraintForceAlongPositiveAxis =
            -settings.projection.densityKgPerCubicMeter
            * grid.cellVolumeCubicMeters()
            / settings.projection.timeStepSeconds
            * (face.normalVelocityMetersPerSecond
               - normalFaceVelocity(
                   grid, predictedVelocityMetersPerSecond, face));
        const Vector3 directConstraintForce = axialVector(
            face.axis, directConstraintForceAlongPositiveAxis);
        const double constraintReactionForceAlongPositiveAxis =
            pressureForceAlongPositiveAxis
            + directConstraintForceAlongPositiveAxis;
        const Vector3 constraintReactionForce = axialVector(
            face.axis, constraintReactionForceAlongPositiveAxis);
        const Vector3 constraintReactionTraction = scale(
            constraintReactionForce, 1.0 / area);
        const auto [lowerCorner, upperCorner] = faceCorners(grid, face);
        diagnostics.faces.push_back({
            face.surfaceStableId,
            face.minusRegionStableId,
            face.plusRegionStableId,
            face.axis,
            face.i,
            face.j,
            face.k,
            lowerCorner,
            upperCorner,
            area,
            face.normalVelocityMetersPerSecond,
            traction,
            force,
            pressureForceAlongPositiveAxis
                * face.normalVelocityMetersPerSecond,
            directConstraintForce,
            constraintReactionTraction,
            constraintReactionForce,
            constraintReactionForceAlongPositiveAxis
                * face.normalVelocityMetersPerSecond,
        });
        auto& surface = surfaces[face.surfaceStableId];
        surface.stableId = face.surfaceStableId;
        ++surface.faceCount;
        surface.areaSquareMeters += area;
        surface.pressureForceNewtons = add(
            surface.pressureForceNewtons, force);
        surface.pressurePowerWatts +=
            pressureForceAlongPositiveAxis
            * face.normalVelocityMetersPerSecond;
        surface.directConstraintForceNewtons = add(
            surface.directConstraintForceNewtons,
            directConstraintForce);
        surface.constraintReactionForceNewtons = add(
            surface.constraintReactionForceNewtons,
            constraintReactionForce);
        surface.constraintReactionPowerWatts +=
            constraintReactionForceAlongPositiveAxis
            * face.normalVelocityMetersPerSecond;
    }
    for (const auto& face : diagnostics.faces) {
        auto& surface = surfaces.at(face.surfaceStableId);
        const Vector3 meanTraction = scale(
            surface.pressureForceNewtons,
            1.0 / surface.areaSquareMeters);
        surface.maximumPressureTractionDeviationPascals = std::max(
            surface.maximumPressureTractionDeviationPascals,
            length(subtract(face.pressureTractionPascals, meanTraction)));
        const Vector3 meanConstraintReactionTraction = scale(
            surface.constraintReactionForceNewtons,
            1.0 / surface.areaSquareMeters);
        surface.maximumConstraintReactionTractionDeviationPascals = std::max(
            surface.maximumConstraintReactionTractionDeviationPascals,
            length(subtract(
                face.constraintReactionTractionPascals,
                meanConstraintReactionTraction)));
    }
    diagnostics.surfaces.reserve(surfaces.size());
    for (auto& [stableId, surface] : surfaces) {
        static_cast<void>(stableId);
        surface.pressureImpulseNewtonSeconds = scale(
            surface.pressureForceNewtons,
            settings.projection.timeStepSeconds);
        surface.pressureWorkJoules = surface.pressurePowerWatts
            * settings.projection.timeStepSeconds;
        surface.constraintReactionImpulseNewtonSeconds = scale(
            surface.constraintReactionForceNewtons,
            settings.projection.timeStepSeconds);
        surface.constraintReactionWorkJoules =
            surface.constraintReactionPowerWatts
            * settings.projection.timeStepSeconds;
        diagnostics.totalPressureForceNewtons = add(
            diagnostics.totalPressureForceNewtons,
            surface.pressureForceNewtons);
        diagnostics.totalPressureImpulseNewtonSeconds = add(
            diagnostics.totalPressureImpulseNewtonSeconds,
            surface.pressureImpulseNewtonSeconds);
        diagnostics.totalPressurePowerWatts += surface.pressurePowerWatts;
        diagnostics.totalPressureWorkJoules += surface.pressureWorkJoules;
        diagnostics.totalConstraintReactionForceNewtons = add(
            diagnostics.totalConstraintReactionForceNewtons,
            surface.constraintReactionForceNewtons);
        diagnostics.totalConstraintReactionImpulseNewtonSeconds = add(
            diagnostics.totalConstraintReactionImpulseNewtonSeconds,
            surface.constraintReactionImpulseNewtonSeconds);
        diagnostics.totalConstraintReactionPowerWatts +=
            surface.constraintReactionPowerWatts;
        diagnostics.totalConstraintReactionWorkJoules +=
            surface.constraintReactionWorkJoules;
        diagnostics.surfaces.push_back(surface);
    }
    diagnostics.finite = projection.converged
        && std::isfinite(
            diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumNormalVelocityErrorMetersPerSecond)
        && finite(diagnostics.totalPressureForceNewtons)
        && finite(diagnostics.totalPressureImpulseNewtonSeconds)
        && std::isfinite(diagnostics.totalPressurePowerWatts)
        && std::isfinite(diagnostics.totalPressureWorkJoules)
        && finite(diagnostics.totalConstraintReactionForceNewtons)
        && finite(diagnostics.totalConstraintReactionImpulseNewtonSeconds)
        && std::isfinite(
            diagnostics.totalConstraintReactionPowerWatts)
        && std::isfinite(
            diagnostics.totalConstraintReactionWorkJoules)
        && diagnostics.faces.size() == interfaces.faceCount()
        && std::ranges::all_of(
            diagnostics.faces,
            [](const MovingInterfaceFaceDiagnostics& face) {
                return face.surfaceStableId > 0
                    && face.minusRegionStableId > 0
                    && face.plusRegionStableId > 0
                    && finite(face.lowerCornerMeters)
                    && finite(face.upperCornerMeters)
                    && std::isfinite(face.areaSquareMeters)
                    && face.areaSquareMeters > 0.0
                    && std::isfinite(
                        face.normalVelocityMetersPerSecond)
                    && finite(face.pressureTractionPascals)
                    && finite(face.pressureForceNewtons)
                    && std::isfinite(face.pressurePowerWatts)
                    && finite(face.directConstraintForceNewtons)
                    && finite(face.constraintReactionTractionPascals)
                    && finite(face.constraintReactionForceNewtons)
                    && std::isfinite(
                        face.constraintReactionPowerWatts);
            })
        && std::ranges::all_of(
            diagnostics.surfaces,
            [](const MovingInterfaceSurfaceDiagnostics& surface) {
                return surface.faceCount > 0
                    && std::isfinite(surface.areaSquareMeters)
                    && surface.areaSquareMeters > 0.0
                    && finite(surface.pressureForceNewtons)
                    && finite(surface.pressureImpulseNewtonSeconds)
                    && std::isfinite(
                        surface.maximumPressureTractionDeviationPascals)
                    && std::isfinite(surface.pressurePowerWatts)
                    && std::isfinite(surface.pressureWorkJoules)
                    && finite(surface.directConstraintForceNewtons)
                    && finite(surface.constraintReactionForceNewtons)
                    && finite(
                        surface.constraintReactionImpulseNewtonSeconds)
                    && std::isfinite(
                        surface.maximumConstraintReactionTractionDeviationPascals)
                    && std::isfinite(
                        surface.constraintReactionPowerWatts)
                    && std::isfinite(
                        surface.constraintReactionWorkJoules);
            })
        && isFinite(candidateVelocity) && isFinite(candidatePressure);
    if (!diagnostics.finite) {
        projection.converged = false;
        return diagnostics;
    }

    predictedVelocityMetersPerSecond = std::move(candidateVelocity);
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
