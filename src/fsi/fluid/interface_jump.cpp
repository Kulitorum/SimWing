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

auto canonicalCrossingKey(const GridFacePressureJump& face) {
    return std::tuple{
        axisOrdinal(face.axis), face.k, face.j, face.i,
        face.crossingFraction, face.surfaceStableId,
        face.minusRegionStableId, face.plusRegionStableId};
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

template<typename Callback>
void forEachUniqueFace(const SharpPressureJumpField& field,
                       Callback&& callback) {
    const auto faces = field.faces();
    std::size_t first = 0;
    while (first < faces.size()) {
        callback(faces[first]);
        std::size_t next = first + 1;
        while (next < faces.size()
               && canonicalFaceKey(faces[next])
                   == canonicalFaceKey(faces[first])) {
            ++next;
        }
        first = next;
    }
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
        if (!std::isfinite(face.pressureJumpPascals)
            || !std::isfinite(face.crossingFraction)) {
            throw std::invalid_argument(
                "pressure jumps and crossing fractions must be finite");
        }
        if (!(face.crossingFraction > 0.0)
            || !(face.crossingFraction < 1.0)) {
            throw std::invalid_argument(
                "pressure-jump crossing fractions must lie inside the face-normal cell segment");
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
                  return canonicalCrossingKey(first)
                      < canonicalCrossingKey(second);
              });
    for (std::size_t index = 1; index < faces_.size(); ++index) {
        if (canonicalFaceKey(faces_[index - 1])
            == canonicalFaceKey(faces_[index])) {
            if (!(faces_[index - 1].crossingFraction
                    < faces_[index].crossingFraction)) {
                throw std::invalid_argument(
                    "pressure-jump crossings on one grid face must have distinct ordered positions");
            }
            if (faces_[index - 1].plusRegionStableId
                != faces_[index].minusRegionStableId) {
                throw std::invalid_argument(
                    "pressure-jump crossings on one grid face do not form a continuous region chain");
            }
        }
    }

    for (const auto& face : faces_) {
        const auto index = grid.cellIndex(face.i, face.j, face.k);
        double* aggregate = nullptr;
        switch (face.axis) {
        case GridFaceAxis::X:
            aggregate = &xFaceJumpsPascals_[index];
            break;
        case GridFaceAxis::Y:
            aggregate = &yFaceJumpsPascals_[index];
            break;
        case GridFaceAxis::Z:
            aggregate = &zFaceJumpsPascals_[index];
            break;
        }
        *aggregate += face.pressureJumpPascals;
        if (!std::isfinite(*aggregate)) {
            throw std::invalid_argument(
                "aggregate pressure jump on one grid face is not finite");
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
    forEachUniqueFace(pressureJumps, [&](const auto& face) {
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
    });
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
    const auto xJumps = pressureJumps.xFaceJumpsPascals();
    const auto yJumps = pressureJumps.yFaceJumpsPascals();
    const auto zJumps = pressureJumps.zFaceJumpsPascals();
    forEachUniqueFace(pressureJumps, [&](const auto& face) {
        const auto plusCell = uncheckedIndex(
            counts, face.i, face.j, face.k);
        std::size_t minusCell = 0;
        double inverseSpacingSquared = 0.0;
        double aggregateJumpPascals = 0.0;
        switch (face.axis) {
        case GridFaceAxis::X:
            minusCell = uncheckedIndex(
                counts,
                face.i == 0 ? counts.x - 1 : face.i - 1,
                face.j, face.k);
            inverseSpacingSquared = inverseDxSquared;
            aggregateJumpPascals = xJumps[plusCell];
            break;
        case GridFaceAxis::Y:
            minusCell = uncheckedIndex(
                counts, face.i,
                face.j == 0 ? counts.y - 1 : face.j - 1,
                face.k);
            inverseSpacingSquared = inverseDySquared;
            aggregateJumpPascals = yJumps[plusCell];
            break;
        case GridFaceAxis::Z:
            minusCell = uncheckedIndex(
                counts, face.i, face.j,
                face.k == 0 ? counts.z - 1 : face.k - 1);
            inverseSpacingSquared = inverseDzSquared;
            aggregateJumpPascals = zJumps[plusCell];
            break;
        }
        const double contribution =
            aggregateJumpPascals * inverseSpacingSquared;
        source[minusCell] += contribution;
        source[plusCell] -= contribution;
    });
}

} // namespace simwing::fsi::fluid
