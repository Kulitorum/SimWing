#include "cfd_slice.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <span>

namespace simwing::viewer {
namespace {

constexpr std::size_t axisIndex(const CfdSliceAxis axis) noexcept {
    return static_cast<std::size_t>(axis);
}

double component(const Vec3d& value, const std::size_t axis) noexcept {
    if (axis == 0) {
        return value.x;
    }
    if (axis == 1) {
        return value.y;
    }
    return value.z;
}

void setComponent(Vec3d& value,
                  const std::size_t axis,
                  const double componentValue) noexcept {
    if (axis == 0) {
        value.x = componentValue;
    } else if (axis == 1) {
        value.y = componentValue;
    } else {
        value.z = componentValue;
    }
}

bool metadataSize(const DiagnosticFrame& frame,
                  const std::string_view name,
                  std::size_t& result) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [&](const ScalarField& field) {
            return field.name == name
                && field.association == FieldAssociation::Global;
        });
    if (found == frame.scalarFields.end() || found->values.size() != 1) {
        return false;
    }
    const double value = found->values.front();
    if (!std::isfinite(value) || value < 0.0
        || std::floor(value) != value
        || value > static_cast<double>(
            std::numeric_limits<std::size_t>::max())) {
        return false;
    }
    result = static_cast<std::size_t>(value);
    return true;
}

bool multiplyBounded(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) noexcept {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

bool uniformIncreasingCoordinates(std::span<const double> coordinates,
                                  double& spacing) noexcept {
    if (coordinates.empty()) {
        return false;
    }
    if (coordinates.size() == 1) {
        spacing = 0.0;
        return std::isfinite(coordinates.front());
    }
    spacing = coordinates[1] - coordinates[0];
    if (!std::isfinite(spacing) || !(spacing > 0.0)) {
        return false;
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(coordinates.front()),
                    std::abs(coordinates.back()), std::abs(spacing)});
    for (std::size_t index = 1; index < coordinates.size(); ++index) {
        const double delta = coordinates[index] - coordinates[index - 1];
        if (!std::isfinite(coordinates[index])
            || std::abs(delta - spacing) > tolerance) {
            return false;
        }
    }
    return true;
}

std::size_t gridIndex(const std::array<std::size_t, 3>& counts,
                      const std::size_t i,
                      const std::size_t j,
                      const std::size_t k) noexcept {
    return i + counts[0] * (j + counts[1] * k);
}

} // namespace

std::size_t CfdGridDescriptor::cellCount() const noexcept {
    return cellCounts[0] * cellCounts[1] * cellCounts[2];
}

bool CfdSliceQuad::operator==(const CfdSliceQuad& other) const noexcept {
    if (sourceVertexIndex != other.sourceVertexIndex) {
        return false;
    }
    for (std::size_t index = 0; index < cornersMetres.size(); ++index) {
        const Vec3d& first = cornersMetres[index];
        const Vec3d& second = other.cornersMetres[index];
        if (first.x != second.x || first.y != second.y
            || first.z != second.z) {
            return false;
        }
    }
    return true;
}

std::optional<CfdGridDescriptor> describeCfdGrid(
    const DiagnosticFrame& frame) {
    CfdGridDescriptor result;
    if (!metadataSize(
            frame, cfdGridVertexBeginFieldName, result.vertexBegin)
        || !metadataSize(
            frame, cfdGridCellCountXFieldName, result.cellCounts[0])
        || !metadataSize(
            frame, cfdGridCellCountYFieldName, result.cellCounts[1])
        || !metadataSize(
            frame, cfdGridCellCountZFieldName, result.cellCounts[2])
        || std::ranges::any_of(
            result.cellCounts,
            [](const std::size_t count) { return count == 0; })) {
        return std::nullopt;
    }

    std::size_t xy = 0;
    std::size_t count = 0;
    if (!multiplyBounded(result.cellCounts[0], result.cellCounts[1], xy)
        || !multiplyBounded(xy, result.cellCounts[2], count)
        || result.vertexBegin > frame.vertices.size()
        || count > frame.vertices.size() - result.vertexBegin) {
        return std::nullopt;
    }

    std::vector<bool> referenced(frame.vertices.size(), false);
    for (const DiagnosticTriangle& triangle : frame.triangles) {
        if (triangle.vertex0 >= referenced.size()
            || triangle.vertex1 >= referenced.size()
            || triangle.vertex2 >= referenced.size()) {
            return std::nullopt;
        }
        referenced[triangle.vertex0] = true;
        referenced[triangle.vertex1] = true;
        referenced[triangle.vertex2] = true;
    }
    for (const DiagnosticLine& line : frame.lines) {
        if (line.vertex0 >= referenced.size()
            || line.vertex1 >= referenced.size()) {
            return std::nullopt;
        }
        referenced[line.vertex0] = true;
        referenced[line.vertex1] = true;
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (referenced[result.vertexBegin + index]) {
            return std::nullopt;
        }
    }

    for (std::size_t axis = 0; axis < 3; ++axis) {
        auto& coordinates = result.cellCentreCoordinatesMetres[axis];
        coordinates.reserve(result.cellCounts[axis]);
        for (std::size_t index = 0; index < result.cellCounts[axis]; ++index) {
            std::array<std::size_t, 3> coordinateIndex{};
            coordinateIndex[axis] = index;
            const std::size_t vertexIndex = result.vertexBegin + gridIndex(
                result.cellCounts, coordinateIndex[0], coordinateIndex[1],
                coordinateIndex[2]);
            coordinates.push_back(component(
                frame.vertices[vertexIndex].positionMetres, axis));
        }
        if (!uniformIncreasingCoordinates(
                coordinates, result.cellSpacingMetres[axis])) {
            return std::nullopt;
        }
    }

    for (std::size_t k = 0; k < result.cellCounts[2]; ++k) {
        for (std::size_t j = 0; j < result.cellCounts[1]; ++j) {
            for (std::size_t i = 0; i < result.cellCounts[0]; ++i) {
                const Vec3d& position = frame.vertices[
                    result.vertexBegin
                    + gridIndex(result.cellCounts, i, j, k)].positionMetres;
                if (position.x
                        != result.cellCentreCoordinatesMetres[0][i]
                    || position.y
                        != result.cellCentreCoordinatesMetres[1][j]
                    || position.z
                        != result.cellCentreCoordinatesMetres[2][k]) {
                    return std::nullopt;
                }
            }
        }
    }
    return result;
}

std::optional<CfdSliceGeometry> buildCfdSliceGeometry(
    const DiagnosticFrame& frame,
    const CfdSliceAxis axis,
    const std::size_t sliceIndex) {
    const auto descriptor = describeCfdGrid(frame);
    const std::size_t normalAxis = axisIndex(axis);
    if (!descriptor || normalAxis >= 3
        || sliceIndex >= descriptor->cellCounts[normalAxis]) {
        return std::nullopt;
    }
    const std::size_t firstTangent = (normalAxis + 1) % 3;
    const std::size_t secondTangent = (normalAxis + 2) % 3;
    if (!(descriptor->cellSpacingMetres[firstTangent] > 0.0)
        || !(descriptor->cellSpacingMetres[secondTangent] > 0.0)) {
        return std::nullopt;
    }

    CfdSliceGeometry result;
    result.axis = axis;
    result.sliceIndex = sliceIndex;
    result.coordinateMetres =
        descriptor->cellCentreCoordinatesMetres[normalAxis][sliceIndex];
    result.quads.reserve(
        descriptor->cellCounts[firstTangent]
        * descriptor->cellCounts[secondTangent]);
    const double firstHalf =
        0.5 * descriptor->cellSpacingMetres[firstTangent];
    const double secondHalf =
        0.5 * descriptor->cellSpacingMetres[secondTangent];

    for (std::size_t second = 0;
         second < descriptor->cellCounts[secondTangent]; ++second) {
        for (std::size_t first = 0;
             first < descriptor->cellCounts[firstTangent]; ++first) {
            std::array<std::size_t, 3> coordinateIndex{};
            coordinateIndex[normalAxis] = sliceIndex;
            coordinateIndex[firstTangent] = first;
            coordinateIndex[secondTangent] = second;
            const std::size_t source = descriptor->vertexBegin + gridIndex(
                descriptor->cellCounts, coordinateIndex[0],
                coordinateIndex[1], coordinateIndex[2]);
            const Vec3d centre = frame.vertices[source].positionMetres;
            CfdSliceQuad quad;
            quad.sourceVertexIndex = static_cast<std::uint32_t>(source);
            quad.cornersMetres.fill(centre);
            constexpr std::array<double, 4> firstSigns{-1.0, 1.0, 1.0, -1.0};
            constexpr std::array<double, 4> secondSigns{-1.0, -1.0, 1.0, 1.0};
            for (std::size_t corner = 0; corner < 4; ++corner) {
                setComponent(
                    quad.cornersMetres[corner], firstTangent,
                    component(centre, firstTangent)
                        + firstSigns[corner] * firstHalf);
                setComponent(
                    quad.cornersMetres[corner], secondTangent,
                    component(centre, secondTangent)
                        + secondSigns[corner] * secondHalf);
            }
            result.quads.push_back(std::move(quad));
        }
    }
    return result;
}

} // namespace simwing::viewer
