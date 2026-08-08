#include "vector_glyphs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::viewer {
namespace {

bool isFinite(const Vec3d& value) noexcept {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vec3d add(const Vec3d& left, const Vec3d& right) noexcept {
    return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Vec3d subtract(const Vec3d& left, const Vec3d& right) noexcept {
    return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Vec3d multiply(const Vec3d& value, const double scale) noexcept {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vec3d cross(const Vec3d& left, const Vec3d& right) noexcept {
    return {
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    };
}

double magnitude(const Vec3d& value) noexcept {
    return std::hypot(value.x, value.y, value.z);
}

Vec3d normalized(const Vec3d& value, const double length) noexcept {
    return multiply(value, 1.0 / length);
}

void validateSettings(const VectorGlyphSettings& settings) {
    if (!std::isfinite(settings.maximumLengthMetres)
        || settings.maximumLengthMetres < 0.0
        || !std::isfinite(settings.automaticLengthFactor)
        || settings.automaticLengthFactor <= 0.0
        || !std::isfinite(settings.arrowheadLengthFraction)
        || settings.arrowheadLengthFraction <= 0.0
        || settings.arrowheadLengthFraction > 1.0
        || !std::isfinite(settings.arrowheadWidthFraction)
        || settings.arrowheadWidthFraction <= 0.0
        || settings.arrowheadWidthFraction > 1.0
        || !std::isfinite(settings.minimumRelativeMagnitude)
        || settings.minimumRelativeMagnitude < 0.0
        || settings.minimumRelativeMagnitude > 1.0
        || settings.maximumGlyphCount == 0) {
        throw std::invalid_argument("invalid vector glyph settings");
    }
}

double characteristicSpacing(
    const std::vector<DiagnosticVertex>& vertices) {
    for (const DiagnosticVertex& vertex : vertices) {
        if (!isFinite(vertex.positionMetres)) {
            throw std::invalid_argument(
                "vector glyph vertices must have finite positions");
        }
    }
    if (vertices.size() < 2) {
        return 0.0;
    }
    Vec3d minimum = vertices.front().positionMetres;
    Vec3d maximum = minimum;
    for (const DiagnosticVertex& vertex : vertices) {
        minimum.x = std::min(minimum.x, vertex.positionMetres.x);
        minimum.y = std::min(minimum.y, vertex.positionMetres.y);
        minimum.z = std::min(minimum.z, vertex.positionMetres.z);
        maximum.x = std::max(maximum.x, vertex.positionMetres.x);
        maximum.y = std::max(maximum.y, vertex.positionMetres.y);
        maximum.z = std::max(maximum.z, vertex.positionMetres.z);
    }

    std::array<double, 3> extents{
        maximum.x - minimum.x,
        maximum.y - minimum.y,
        maximum.z - minimum.z,
    };
    double logarithmicMeasure = 0.0;
    std::size_t dimension = 0;
    for (const double extent : extents) {
        if (!std::isfinite(extent)) {
            throw std::invalid_argument(
                "vector glyph position extent must be finite");
        }
        if (extent > 0.0) {
            logarithmicMeasure += std::log(extent);
            ++dimension;
        }
    }
    if (dimension == 0) {
        return 0.0;
    }
    const double spacing = std::exp(
        (logarithmicMeasure
         - std::log(static_cast<double>(vertices.size())))
        / static_cast<double>(dimension));
    if (!std::isfinite(spacing)) {
        throw std::invalid_argument(
            "vector glyph characteristic spacing must be finite");
    }
    return spacing;
}

} // namespace

VectorGlyphGeometry buildVertexVectorGlyphs(
    const std::vector<DiagnosticVertex>& vertices,
    const VectorField& field,
    const VectorGlyphSettings& settings) {
    validateSettings(settings);
    if (field.association != FieldAssociation::Vertex
        || field.values.size() != vertices.size()) {
        throw std::invalid_argument(
            "vector glyph field must contain one vector per vertex");
    }
    if (vertices.size()
        > static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max())) {
        throw std::invalid_argument(
            "vector glyph vertex count exceeds index capacity");
    }

    VectorGlyphGeometry result;
    result.characteristicSpacingMetres = characteristicSpacing(vertices);
    for (const Vec3d& value : field.values) {
        if (!isFinite(value)) {
            throw std::invalid_argument(
                "vector glyph field values must be finite");
        }
        const double valueMagnitude = magnitude(value);
        if (!std::isfinite(valueMagnitude)) {
            throw std::invalid_argument(
                "vector glyph magnitude must be finite");
        }
        result.maximumVectorMagnitude = std::max(
            result.maximumVectorMagnitude, valueMagnitude);
    }
    result.maximumGlyphLengthMetres = settings.maximumLengthMetres > 0.0
        ? settings.maximumLengthMetres
        : result.characteristicSpacingMetres
            * settings.automaticLengthFactor;
    if (!std::isfinite(result.maximumGlyphLengthMetres)) {
        throw std::invalid_argument(
            "vector glyph maximum length must be finite");
    }
    if (vertices.empty()
        || result.maximumVectorMagnitude == 0.0
        || result.maximumGlyphLengthMetres == 0.0) {
        return result;
    }

    const std::size_t stride = std::max<std::size_t>(
        1,
        (vertices.size()
         + static_cast<std::size_t>(settings.maximumGlyphCount) - 1)
            / static_cast<std::size_t>(settings.maximumGlyphCount));
    result.segments.reserve(
        std::min(vertices.size(),
                 static_cast<std::size_t>(settings.maximumGlyphCount)) * 3);
    const double minimumMagnitude = result.maximumVectorMagnitude
        * settings.minimumRelativeMagnitude;
    for (std::size_t index = 0; index < vertices.size(); index += stride) {
        const Vec3d value = field.values[index];
        const double vectorMagnitude = magnitude(value);
        if (vectorMagnitude <= minimumMagnitude) {
            continue;
        }
        const Vec3d direction = normalized(value, vectorMagnitude);
        const double glyphLength = result.maximumGlyphLengthMetres
            * vectorMagnitude / result.maximumVectorMagnitude;
        const Vec3d start = vertices[index].positionMetres;
        const Vec3d tip = add(start, multiply(direction, glyphLength));
        const Vec3d reference = std::abs(direction.z) < 0.9
            ? Vec3d{0.0, 0.0, 1.0}
            : Vec3d{0.0, 1.0, 0.0};
        const Vec3d sideValue = cross(direction, reference);
        const Vec3d side = normalized(sideValue, magnitude(sideValue));
        const Vec3d arrowBase = subtract(
            tip,
            multiply(direction,
                     glyphLength * settings.arrowheadLengthFraction));
        const Vec3d arrowWidth = multiply(
            side, glyphLength * settings.arrowheadWidthFraction);
        const Vec3d firstWing = add(arrowBase, arrowWidth);
        const Vec3d secondWing = subtract(arrowBase, arrowWidth);
        if (!isFinite(tip) || !isFinite(firstWing) || !isFinite(secondWing)) {
            throw std::invalid_argument(
                "vector glyph derived geometry must be finite");
        }
        const std::uint32_t source = static_cast<std::uint32_t>(index);
        result.segments.push_back({start, tip, source});
        result.segments.push_back({tip, firstWing, source});
        result.segments.push_back({tip, secondWing, source});
        ++result.glyphCount;
    }
    return result;
}

} // namespace simwing::viewer
