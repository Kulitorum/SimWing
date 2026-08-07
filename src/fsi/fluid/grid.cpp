#include "fluid/grid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

std::size_t checkedProduct(std::size_t first,
                           std::size_t second,
                           const char* message) {
    if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::invalid_argument(message);
    }
    return first * second;
}

void requireMatching(const PeriodicCartesianGrid& grid,
                     const CellScalarField& field,
                     const char* message) {
    if (!field.matches(grid)) {
        throw std::invalid_argument(message);
    }
}

void requireMatching(const PeriodicCartesianGrid& grid,
                     const MacVelocityField& field,
                     const char* message) {
    if (!field.matches(grid)) {
        throw std::invalid_argument(message);
    }
}

std::size_t uncheckedIndex(const GridCellCounts counts,
                           const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) {
    return i + counts.x * (j + counts.y * k);
}

} // namespace

PeriodicCartesianGrid::PeriodicCartesianGrid(const GridCellCounts cellCounts,
                                             const Vector3 lowerMeters,
                                             const Vector3 upperMeters)
    : cellCounts_(cellCounts),
      lowerMeters_(lowerMeters),
      upperMeters_(upperMeters) {
    if (cellCounts.x < 2 || cellCounts.y < 2 || cellCounts.z < 2) {
        throw std::invalid_argument(
            "periodic Cartesian grids require at least two cells per axis");
    }
    if (!finiteVector(lowerMeters) || !finiteVector(upperMeters)
        || !(upperMeters.x > lowerMeters.x)
        || !(upperMeters.y > lowerMeters.y)
        || !(upperMeters.z > lowerMeters.z)) {
        throw std::invalid_argument(
            "Cartesian grid bounds must be finite and strictly increasing");
    }
    cellCount_ = checkedProduct(
        checkedProduct(cellCounts.x, cellCounts.y,
                       "Cartesian grid cell count overflows size_t"),
        cellCounts.z,
        "Cartesian grid cell count overflows size_t");
    if (cellCount_ > std::vector<double>().max_size()) {
        throw std::invalid_argument(
            "Cartesian grid has more cells than one field can store");
    }
    cellSpacingMeters_ = {
        (upperMeters.x - lowerMeters.x) / static_cast<double>(cellCounts.x),
        (upperMeters.y - lowerMeters.y) / static_cast<double>(cellCounts.y),
        (upperMeters.z - lowerMeters.z) / static_cast<double>(cellCounts.z),
    };
    if (!finiteVector(cellSpacingMeters_)
        || !(cellSpacingMeters_.x > 0.0)
        || !(cellSpacingMeters_.y > 0.0)
        || !(cellSpacingMeters_.z > 0.0)) {
        throw std::invalid_argument("Cartesian grid spacing is not representable");
    }
    if (!std::isfinite(cellVolumeCubicMeters())
        || !(cellVolumeCubicMeters() > 0.0)) {
        throw std::invalid_argument("Cartesian grid cell volume is not representable");
    }
}

GridCellCounts PeriodicCartesianGrid::cellCounts() const noexcept {
    return cellCounts_;
}

Vector3 PeriodicCartesianGrid::lowerMeters() const noexcept {
    return lowerMeters_;
}

Vector3 PeriodicCartesianGrid::upperMeters() const noexcept {
    return upperMeters_;
}

Vector3 PeriodicCartesianGrid::cellSpacingMeters() const noexcept {
    return cellSpacingMeters_;
}

std::size_t PeriodicCartesianGrid::cellCount() const noexcept {
    return cellCount_;
}

double PeriodicCartesianGrid::cellVolumeCubicMeters() const noexcept {
    return cellSpacingMeters_.x
        * cellSpacingMeters_.y
        * cellSpacingMeters_.z;
}

std::size_t PeriodicCartesianGrid::cellIndex(const std::size_t i,
                                             const std::size_t j,
                                             const std::size_t k) const {
    if (i >= cellCounts_.x || j >= cellCounts_.y || k >= cellCounts_.z) {
        throw std::out_of_range("Cartesian grid cell index is out of range");
    }
    return uncheckedIndex(cellCounts_, i, j, k);
}

Vector3 PeriodicCartesianGrid::cellCenterMeters(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
    static_cast<void>(cellIndex(i, j, k));
    return {
        lowerMeters_.x + (static_cast<double>(i) + 0.5) * cellSpacingMeters_.x,
        lowerMeters_.y + (static_cast<double>(j) + 0.5) * cellSpacingMeters_.y,
        lowerMeters_.z + (static_cast<double>(k) + 0.5) * cellSpacingMeters_.z,
    };
}

Vector3 PeriodicCartesianGrid::xFaceCenterMeters(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
    static_cast<void>(cellIndex(i, j, k));
    return {
        lowerMeters_.x + static_cast<double>(i) * cellSpacingMeters_.x,
        lowerMeters_.y + (static_cast<double>(j) + 0.5) * cellSpacingMeters_.y,
        lowerMeters_.z + (static_cast<double>(k) + 0.5) * cellSpacingMeters_.z,
    };
}

Vector3 PeriodicCartesianGrid::yFaceCenterMeters(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
    static_cast<void>(cellIndex(i, j, k));
    return {
        lowerMeters_.x + (static_cast<double>(i) + 0.5) * cellSpacingMeters_.x,
        lowerMeters_.y + static_cast<double>(j) * cellSpacingMeters_.y,
        lowerMeters_.z + (static_cast<double>(k) + 0.5) * cellSpacingMeters_.z,
    };
}

Vector3 PeriodicCartesianGrid::zFaceCenterMeters(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
    static_cast<void>(cellIndex(i, j, k));
    return {
        lowerMeters_.x + (static_cast<double>(i) + 0.5) * cellSpacingMeters_.x,
        lowerMeters_.y + (static_cast<double>(j) + 0.5) * cellSpacingMeters_.y,
        lowerMeters_.z + static_cast<double>(k) * cellSpacingMeters_.z,
    };
}

CellScalarField::CellScalarField(const PeriodicCartesianGrid& grid,
                                 const double initialValue)
    : cellCounts_(grid.cellCounts()),
      values_(grid.cellCount(), initialValue) {}

GridCellCounts CellScalarField::cellCounts() const noexcept {
    return cellCounts_;
}

bool CellScalarField::matches(const PeriodicCartesianGrid& grid) const noexcept {
    return cellCounts_ == grid.cellCounts() && values_.size() == grid.cellCount();
}

std::span<double> CellScalarField::values() noexcept {
    return values_;
}

std::span<const double> CellScalarField::values() const noexcept {
    return values_;
}

MacVelocityField::MacVelocityField(const PeriodicCartesianGrid& grid,
                                   const double initialValue)
    : cellCounts_(grid.cellCounts()) {
    for (auto& component : components_) {
        component.assign(grid.cellCount(), initialValue);
    }
}

GridCellCounts MacVelocityField::cellCounts() const noexcept {
    return cellCounts_;
}

bool MacVelocityField::matches(const PeriodicCartesianGrid& grid) const noexcept {
    return cellCounts_ == grid.cellCounts()
        && std::ranges::all_of(components_, [&grid](const auto& component) {
            return component.size() == grid.cellCount();
        });
}

std::span<double> MacVelocityField::xFaces() noexcept {
    return components_[0];
}

std::span<const double> MacVelocityField::xFaces() const noexcept {
    return components_[0];
}

std::span<double> MacVelocityField::yFaces() noexcept {
    return components_[1];
}

std::span<const double> MacVelocityField::yFaces() const noexcept {
    return components_[1];
}

std::span<double> MacVelocityField::zFaces() noexcept {
    return components_[2];
}

std::span<const double> MacVelocityField::zFaces() const noexcept {
    return components_[2];
}

bool isFinite(const CellScalarField& field) noexcept {
    return std::ranges::all_of(field.values(), [](const double value) {
        return std::isfinite(value);
    });
}

bool isFinite(const MacVelocityField& field) noexcept {
    const auto finite = [](const double value) { return std::isfinite(value); };
    return std::ranges::all_of(field.xFaces(), finite)
        && std::ranges::all_of(field.yFaces(), finite)
        && std::ranges::all_of(field.zFaces(), finite);
}

void computeDivergence(const PeriodicCartesianGrid& grid,
                       const MacVelocityField& velocity,
                       CellScalarField& divergencePerSecond) {
    requireMatching(grid, velocity, "velocity field does not match its grid");
    requireMatching(grid, divergencePerSecond,
                    "divergence field does not match its grid");
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    const auto xFaces = velocity.xFaces();
    const auto yFaces = velocity.yFaces();
    const auto zFaces = velocity.zFaces();
    auto divergence = divergencePerSecond.values();
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t nextK = k + 1 == counts.z ? 0 : k + 1;
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const auto center = uncheckedIndex(counts, i, j, k);
                divergence[center] =
                    (xFaces[uncheckedIndex(counts, nextI, j, k)]
                     - xFaces[center]) / spacing.x
                    + (yFaces[uncheckedIndex(counts, i, nextJ, k)]
                       - yFaces[center]) / spacing.y
                    + (zFaces[uncheckedIndex(counts, i, j, nextK)]
                       - zFaces[center]) / spacing.z;
            }
        }
    }
}

void computePressureGradient(const PeriodicCartesianGrid& grid,
                             const CellScalarField& pressurePascals,
                             MacVelocityField& gradientPascalsPerMeter) {
    requireMatching(grid, pressurePascals,
                    "pressure field does not match its grid");
    requireMatching(grid, gradientPascalsPerMeter,
                    "pressure-gradient field does not match its grid");
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    const auto pressure = pressurePascals.values();
    auto xGradient = gradientPascalsPerMeter.xFaces();
    auto yGradient = gradientPascalsPerMeter.yFaces();
    auto zGradient = gradientPascalsPerMeter.zFaces();
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t previousK = k == 0 ? counts.z - 1 : k - 1;
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t previousJ = j == 0 ? counts.y - 1 : j - 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t previousI = i == 0 ? counts.x - 1 : i - 1;
                const auto center = uncheckedIndex(counts, i, j, k);
                xGradient[center] =
                    (pressure[center]
                     - pressure[uncheckedIndex(counts, previousI, j, k)])
                    / spacing.x;
                yGradient[center] =
                    (pressure[center]
                     - pressure[uncheckedIndex(counts, i, previousJ, k)])
                    / spacing.y;
                zGradient[center] =
                    (pressure[center]
                     - pressure[uncheckedIndex(counts, i, j, previousK)])
                    / spacing.z;
            }
        }
    }
}

void applyNegativeLaplacian(const PeriodicCartesianGrid& grid,
                            const CellScalarField& input,
                            CellScalarField& output) {
    requireMatching(grid, input, "Laplacian input does not match its grid");
    requireMatching(grid, output, "Laplacian output does not match its grid");
    if (&input == &output) {
        throw std::invalid_argument("Laplacian input and output must not alias");
    }
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    const double inverseDxSquared = 1.0 / (spacing.x * spacing.x);
    const double inverseDySquared = 1.0 / (spacing.y * spacing.y);
    const double inverseDzSquared = 1.0 / (spacing.z * spacing.z);
    const auto values = input.values();
    auto result = output.values();
    for (std::size_t k = 0; k < counts.z; ++k) {
        const std::size_t previousK = k == 0 ? counts.z - 1 : k - 1;
        const std::size_t nextK = k + 1 == counts.z ? 0 : k + 1;
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t previousJ = j == 0 ? counts.y - 1 : j - 1;
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t previousI = i == 0 ? counts.x - 1 : i - 1;
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const auto center = uncheckedIndex(counts, i, j, k);
                result[center] =
                    (2.0 * values[center]
                     - values[uncheckedIndex(counts, previousI, j, k)]
                     - values[uncheckedIndex(counts, nextI, j, k)])
                        * inverseDxSquared
                    + (2.0 * values[center]
                       - values[uncheckedIndex(counts, i, previousJ, k)]
                       - values[uncheckedIndex(counts, i, nextJ, k)])
                        * inverseDySquared
                    + (2.0 * values[center]
                       - values[uncheckedIndex(counts, i, j, previousK)]
                       - values[uncheckedIndex(counts, i, j, nextK)])
                        * inverseDzSquared;
            }
        }
    }
}

double mean(const CellScalarField& field) {
    if (field.values().empty()) {
        throw std::invalid_argument("cannot compute the mean of an empty field");
    }
    return std::accumulate(field.values().begin(), field.values().end(), 0.0)
        / static_cast<double>(field.values().size());
}

double l2Norm(const CellScalarField& field) {
    if (field.values().empty()) {
        throw std::invalid_argument("cannot compute the norm of an empty field");
    }
    double sumSquares = 0.0;
    for (const double value : field.values()) {
        sumSquares += value * value;
    }
    return std::sqrt(sumSquares / static_cast<double>(field.values().size()));
}

double maximumAbsoluteValue(const CellScalarField& field) {
    double maximum = 0.0;
    for (const double value : field.values()) {
        maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

double kineticEnergyJoules(const PeriodicCartesianGrid& grid,
                           const MacVelocityField& velocity,
                           const double densityKgPerCubicMeter) {
    requireMatching(grid, velocity, "velocity field does not match its grid");
    if (!std::isfinite(densityKgPerCubicMeter)
        || !(densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument("fluid density must be finite and positive");
    }
    double sumSquares = 0.0;
    for (const double value : velocity.xFaces()) {
        sumSquares += value * value;
    }
    for (const double value : velocity.yFaces()) {
        sumSquares += value * value;
    }
    for (const double value : velocity.zFaces()) {
        sumSquares += value * value;
    }
    return 0.5 * densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters() * sumSquares;
}

} // namespace simwing::fsi::fluid
