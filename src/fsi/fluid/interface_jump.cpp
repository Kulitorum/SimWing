#include "fluid/interface_jump.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace simwing::fsi::fluid {
namespace {

std::uint8_t axisOrdinal(const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X:
        return 0;
    case GridFaceAxis::Y:
        return 1;
    case GridFaceAxis::Z:
        return 2;
    }
    throw std::invalid_argument("pressure-jump face has an invalid axis");
}

auto canonicalFaceKey(const GridFacePressureJump& face) {
    return std::tuple{
        axisOrdinal(face.axis), face.k, face.j, face.i};
}

std::pair<std::uint64_t, std::uint64_t> unorderedRegionPair(
    const GridFacePressureJump& face) {
    return std::minmax(
        face.minusRegionStableId, face.plusRegionStableId);
}

void requireMatching(const PeriodicCartesianGrid& grid,
                     const SharpPressureJumpField& field) {
    if (!field.matches(grid)) {
        throw std::invalid_argument(
            "sharp pressure-jump field does not match its grid");
    }
}

std::size_t uncheckedIndex(const GridCellCounts counts,
                           const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) {
    return i + counts.x * (j + counts.y * k);
}

} // namespace

SharpPressureJumpField::SharpPressureJumpField(
    const PeriodicCartesianGrid& grid,
    std::vector<GridFacePressureJump> faces)
    : cellCounts_(grid.cellCounts()),
      faces_(std::move(faces)),
      xFaceJumpsPascals_(grid.cellCount(), 0.0),
      yFaceJumpsPascals_(grid.cellCount(), 0.0),
      zFaceJumpsPascals_(grid.cellCount(), 0.0) {
    std::map<std::uint64_t, std::pair<std::uint64_t, std::uint64_t>>
        surfaceRegions;
    for (const auto& face : faces_) {
        static_cast<void>(axisOrdinal(face.axis));
        if (face.surfaceStableId == 0
            || face.minusRegionStableId == 0
            || face.plusRegionStableId == 0) {
            throw std::invalid_argument(
                "pressure-jump surface and region IDs must be nonzero");
        }
        if (face.minusRegionStableId == face.plusRegionStableId) {
            throw std::invalid_argument(
                "a pressure jump must separate two distinct regions");
        }
        if (!std::isfinite(face.pressureJumpPascals)) {
            throw std::invalid_argument("pressure jumps must be finite");
        }
        if (face.i >= cellCounts_.x
            || face.j >= cellCounts_.y
            || face.k >= cellCounts_.z) {
            throw std::invalid_argument(
                "pressure-jump face index is out of range");
        }
        const auto regions = unorderedRegionPair(face);
        const auto [iterator, inserted] = surfaceRegions.emplace(
            face.surfaceStableId, regions);
        if (!inserted && iterator->second != regions) {
            throw std::invalid_argument(
                "one pressure-jump surface ID cannot separate different region pairs");
        }
    }

    std::sort(faces_.begin(), faces_.end(),
              [](const auto& first, const auto& second) {
                  return canonicalFaceKey(first) < canonicalFaceKey(second);
              });
    for (std::size_t index = 1; index < faces_.size(); ++index) {
        if (canonicalFaceKey(faces_[index - 1])
            == canonicalFaceKey(faces_[index])) {
            throw std::invalid_argument(
                "multiple pressure-jump crossings on one grid face are not yet supported");
        }
    }

    for (const auto& face : faces_) {
        const auto index = grid.cellIndex(face.i, face.j, face.k);
        switch (face.axis) {
        case GridFaceAxis::X:
            xFaceJumpsPascals_[index] = face.pressureJumpPascals;
            break;
        case GridFaceAxis::Y:
            yFaceJumpsPascals_[index] = face.pressureJumpPascals;
            break;
        case GridFaceAxis::Z:
            zFaceJumpsPascals_[index] = face.pressureJumpPascals;
            break;
        }
    }
}

GridCellCounts SharpPressureJumpField::cellCounts() const noexcept {
    return cellCounts_;
}

bool SharpPressureJumpField::matches(
    const PeriodicCartesianGrid& grid) const noexcept {
    return cellCounts_ == grid.cellCounts()
        && xFaceJumpsPascals_.size() == grid.cellCount()
        && yFaceJumpsPascals_.size() == grid.cellCount()
        && zFaceJumpsPascals_.size() == grid.cellCount();
}

bool SharpPressureJumpField::empty() const noexcept {
    return faces_.empty();
}

std::size_t SharpPressureJumpField::faceCount() const noexcept {
    return faces_.size();
}

std::span<const GridFacePressureJump>
SharpPressureJumpField::faces() const noexcept {
    return faces_;
}

std::span<const double>
SharpPressureJumpField::xFaceJumpsPascals() const noexcept {
    return xFaceJumpsPascals_;
}

std::span<const double>
SharpPressureJumpField::yFaceJumpsPascals() const noexcept {
    return yFaceJumpsPascals_;
}

std::span<const double>
SharpPressureJumpField::zFaceJumpsPascals() const noexcept {
    return zFaceJumpsPascals_;
}

void computePressureGradientWithJumps(
    const PeriodicCartesianGrid& grid,
    const CellScalarField& pressurePascals,
    const SharpPressureJumpField& pressureJumps,
    MacVelocityField& gradientPascalsPerMeter) {
    requireMatching(grid, pressureJumps);
    computePressureGradient(grid, pressurePascals, gradientPascalsPerMeter);
    const auto spacing = grid.cellSpacingMeters();
    const auto xJumps = pressureJumps.xFaceJumpsPascals();
    const auto yJumps = pressureJumps.yFaceJumpsPascals();
    const auto zJumps = pressureJumps.zFaceJumpsPascals();
    auto xGradient = gradientPascalsPerMeter.xFaces();
    auto yGradient = gradientPascalsPerMeter.yFaces();
    auto zGradient = gradientPascalsPerMeter.zFaces();
    const auto counts = grid.cellCounts();
    for (const auto& face : pressureJumps.faces()) {
        const auto index = uncheckedIndex(counts, face.i, face.j, face.k);
        switch (face.axis) {
        case GridFaceAxis::X:
            xGradient[index] -= xJumps[index] / spacing.x;
            break;
        case GridFaceAxis::Y:
            yGradient[index] -= yJumps[index] / spacing.y;
            break;
        case GridFaceAxis::Z:
            zGradient[index] -= zJumps[index] / spacing.z;
            break;
        }
    }
}

void computePressureJumpSource(
    const PeriodicCartesianGrid& grid,
    const SharpPressureJumpField& pressureJumps,
    CellScalarField& sourcePascalsPerSquareMeter) {
    requireMatching(grid, pressureJumps);
    if (!sourcePascalsPerSquareMeter.matches(grid)) {
        throw std::invalid_argument(
            "pressure-jump source field does not match its grid");
    }
    const auto spacing = grid.cellSpacingMeters();
    const double inverseDxSquared = 1.0 / (spacing.x * spacing.x);
    const double inverseDySquared = 1.0 / (spacing.y * spacing.y);
    const double inverseDzSquared = 1.0 / (spacing.z * spacing.z);
    auto source = sourcePascalsPerSquareMeter.values();
    std::fill(source.begin(), source.end(), 0.0);
    const auto counts = grid.cellCounts();
    for (const auto& face : pressureJumps.faces()) {
        const auto plusCell = uncheckedIndex(
            counts, face.i, face.j, face.k);
        std::size_t minusCell = 0;
        double inverseSpacingSquared = 0.0;
        switch (face.axis) {
        case GridFaceAxis::X:
            minusCell = uncheckedIndex(
                counts,
                face.i == 0 ? counts.x - 1 : face.i - 1,
                face.j, face.k);
            inverseSpacingSquared = inverseDxSquared;
            break;
        case GridFaceAxis::Y:
            minusCell = uncheckedIndex(
                counts, face.i,
                face.j == 0 ? counts.y - 1 : face.j - 1,
                face.k);
            inverseSpacingSquared = inverseDySquared;
            break;
        case GridFaceAxis::Z:
            minusCell = uncheckedIndex(
                counts, face.i, face.j,
                face.k == 0 ? counts.z - 1 : face.k - 1);
            inverseSpacingSquared = inverseDzSquared;
            break;
        }
        const double contribution =
            face.pressureJumpPascals * inverseSpacingSquared;
        source[minusCell] += contribution;
        source[plusCell] -= contribution;
    }
}

} // namespace simwing::fsi::fluid
