#pragma once

#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::viewer {

struct VectorGlyphSettings {
    // Zero selects automatic scaling from the dimensional point spacing.
    double maximumLengthMetres = 0.0;
    double automaticLengthFactor = 0.8;
    double arrowheadLengthFraction = 0.28;
    double arrowheadWidthFraction = 0.16;
    double minimumRelativeMagnitude = 1.0e-6;
    std::uint32_t maximumGlyphCount = 100'000;
};

struct VectorGlyphSegment {
    Vec3d startMetres;
    Vec3d endMetres;
    std::uint32_t sourceVertexIndex = 0;

    [[nodiscard]] bool operator==(
        const VectorGlyphSegment& other) const noexcept {
        return startMetres.x == other.startMetres.x
            && startMetres.y == other.startMetres.y
            && startMetres.z == other.startMetres.z
            && endMetres.x == other.endMetres.x
            && endMetres.y == other.endMetres.y
            && endMetres.z == other.endMetres.z
            && sourceVertexIndex == other.sourceVertexIndex;
    }
};

struct VectorGlyphGeometry {
    std::vector<VectorGlyphSegment> segments;
    std::uint32_t glyphCount = 0;
    double maximumVectorMagnitude = 0.0;
    double characteristicSpacingMetres = 0.0;
    double maximumGlyphLengthMetres = 0.0;

    [[nodiscard]] bool operator==(
        const VectorGlyphGeometry&) const noexcept = default;
};

// Builds one shaft and two arrowhead segments for each retained nonzero
// vertex vector. Length is normalized by the field maximum, and large fields
// are sampled by a deterministic integer stride. The result owns its geometry
// and never aliases a frame or renderer buffer.
[[nodiscard]] VectorGlyphGeometry buildVertexVectorGlyphs(
    const std::vector<DiagnosticVertex>& vertices,
    const VectorField& field,
    const VectorGlyphSettings& settings = VectorGlyphSettings{});

} // namespace simwing::viewer
