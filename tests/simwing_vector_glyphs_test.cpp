#include "vector_glyphs.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::viewer;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

bool sameVector(const Vec3d& actual, const Vec3d& expected) {
    return actual.x == expected.x
        && actual.y == expected.y
        && actual.z == expected.z;
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

std::vector<DiagnosticVertex> squareVertices() {
    return {
        {1, {0.0, 0.0, 0.0}},
        {2, {2.0, 0.0, 0.0}},
        {3, {0.0, 2.0, 0.0}},
        {4, {2.0, 2.0, 0.0}},
    };
}

void testDeterministicArrowGeometry() {
    auto vertices = squareVertices();
    VectorField field{
        "velocity", "m/s", FieldAssociation::Vertex,
        {{2.0, 0.0, 0.0},
         {0.0, 1.0, 0.0},
         {},
         {-1.0, 0.0, 0.0}}};
    VectorGlyphSettings settings;
    settings.maximumLengthMetres = 2.0;
    settings.minimumRelativeMagnitude = 0.0;
    const VectorGlyphGeometry first =
        buildVertexVectorGlyphs(vertices, field, settings);
    const VectorGlyphGeometry second =
        buildVertexVectorGlyphs(vertices, field, settings);

    check(first == second,
          "glyphs: repeated geometry construction is deterministic");
    check(first.glyphCount == 3 && first.segments.size() == 9,
          "glyphs: each retained nonzero vector has three segments");
    check(first.segments[0].sourceVertexIndex == 0
              && first.segments[3].sourceVertexIndex == 1
              && first.segments[6].sourceVertexIndex == 3,
          "glyphs: source vertex identity survives deterministic ordering");
    check(sameVector(first.segments[0].startMetres,
                     vertices[0].positionMetres)
              && sameVector(first.segments[0].endMetres,
                            {2.0, 0.0, 0.0})
              && sameVector(first.segments[3].endMetres,
                            {2.0, 1.0, 0.0})
              && sameVector(first.segments[6].endMetres,
                            {1.0, 2.0, 0.0}),
          "glyphs: shaft direction and relative magnitude are exact");
    checkNear(first.maximumVectorMagnitude, 2.0, 0.0,
              "glyphs: maximum vector magnitude is retained");
    checkNear(first.characteristicSpacingMetres, 1.0, 0.0,
              "glyphs: planar characteristic spacing uses area per point");
    checkNear(first.maximumGlyphLengthMetres, 2.0, 0.0,
              "glyphs: explicit maximum length is retained");

    const VectorGlyphGeometry owned = first;
    vertices[0].positionMetres.x = 99.0;
    field.values[0].x = 99.0;
    check(owned == first,
          "glyphs: output owns geometry independent of later input mutation");
}

void testAutomaticScaleAndBoundedSampling() {
    std::vector<DiagnosticVertex> vertices;
    VectorField field{
        "velocity", "m/s", FieldAssociation::Vertex, {}};
    for (std::uint32_t index = 0; index < 10; ++index) {
        vertices.push_back({index + 1,
                            {static_cast<double>(index), 0.0, 0.0}});
        field.values.push_back({1.0, 0.0, 0.0});
    }
    VectorGlyphSettings settings;
    settings.automaticLengthFactor = 0.5;
    settings.maximumGlyphCount = 3;
    const VectorGlyphGeometry geometry =
        buildVertexVectorGlyphs(vertices, field, settings);
    checkNear(geometry.characteristicSpacingMetres, 0.9, 1.0e-15,
              "glyphs: line spacing uses extent per point");
    checkNear(geometry.maximumGlyphLengthMetres, 0.45, 1.0e-15,
              "glyphs: automatic length scales from point spacing");
    check(geometry.glyphCount == 3 && geometry.segments.size() == 9
              && geometry.segments[0].sourceVertexIndex == 0
              && geometry.segments[3].sourceVertexIndex == 4
              && geometry.segments[6].sourceVertexIndex == 8,
          "glyphs: large fields use a bounded deterministic integer stride");
}

void testZeroAndInvalidInputs() {
    auto vertices = squareVertices();
    VectorField field{
        "velocity", "m/s", FieldAssociation::Vertex,
        std::vector<Vec3d>(vertices.size())};
    const VectorGlyphGeometry zero =
        buildVertexVectorGlyphs(vertices, field);
    check(zero.glyphCount == 0 && zero.segments.empty()
              && zero.maximumVectorMagnitude == 0.0,
          "glyphs: a zero field produces no arbitrary directions");

    VectorField invalid = field;
    invalid.association = FieldAssociation::Global;
    expectRejected(
        [&] { static_cast<void>(buildVertexVectorGlyphs(vertices, invalid)); },
        "glyphs: non-vertex fields are rejected");
    invalid = field;
    invalid.values.pop_back();
    expectRejected(
        [&] { static_cast<void>(buildVertexVectorGlyphs(vertices, invalid)); },
        "glyphs: mismatched field cardinality is rejected");
    invalid = field;
    invalid.values[0].x = std::numeric_limits<double>::infinity();
    expectRejected(
        [&] { static_cast<void>(buildVertexVectorGlyphs(vertices, invalid)); },
        "glyphs: non-finite vectors are rejected");
    invalid = field;
    invalid.values[0] = {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()};
    expectRejected(
        [&] { static_cast<void>(buildVertexVectorGlyphs(vertices, invalid)); },
        "glyphs: overflowing derived magnitudes are rejected");
    auto invalidVertices = vertices;
    invalidVertices[0].positionMetres.x =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(
            buildVertexVectorGlyphs(invalidVertices, field)); },
        "glyphs: non-finite positions are rejected");
    VectorGlyphSettings invalidSettings;
    invalidSettings.maximumGlyphCount = 0;
    expectRejected(
        [&] { static_cast<void>(
            buildVertexVectorGlyphs(vertices, field, invalidSettings)); },
        "glyphs: an empty glyph budget is rejected");
}

} // namespace

int main() {
    testDeterministicArrowGeometry();
    testAutomaticScaleAndBoundedSampling();
    testZeroAndInvalidInputs();
    if (failures != 0) {
        std::fprintf(stderr, "%d vector glyph check(s) failed\n", failures);
        return 1;
    }
    std::puts("all vector glyph checks passed");
    return 0;
}
