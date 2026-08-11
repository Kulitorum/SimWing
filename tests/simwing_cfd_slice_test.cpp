#include "cfd_slice.h"

#include <cstdio>
#include <string>

namespace {

using namespace simwing::viewer;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

bool same(const Vec3d& first, const Vec3d& second) {
    return first.x == second.x && first.y == second.y
        && first.z == second.z;
}

void addMetadata(DiagnosticFrame& frame,
                 const std::size_t begin,
                 const std::size_t x,
                 const std::size_t y,
                 const std::size_t z) {
    const auto add = [&](const std::string_view name, const std::size_t value) {
        frame.scalarFields.push_back({
            std::string(name), "1", FieldAssociation::Global,
            {static_cast<double>(value)},
        });
    };
    add(cfdGridVertexBeginFieldName, begin);
    add(cfdGridCellCountXFieldName, x);
    add(cfdGridCellCountYFieldName, y);
    add(cfdGridCellCountZFieldName, z);
}

DiagnosticFrame mixedFrame() {
    DiagnosticFrame frame;
    frame.vertices = {
        {1, {-2.0, 0.0, 0.0}},
        {2, {-1.0, 1.0, 0.0}},
        {3, {-1.0, 0.0, 1.0}},
    };
    frame.triangles = {{10, 0, 1, 2, 1, 2}};
    constexpr std::size_t xCount = 2;
    constexpr std::size_t yCount = 3;
    constexpr std::size_t zCount = 2;
    for (std::size_t k = 0; k < zCount; ++k) {
        for (std::size_t j = 0; j < yCount; ++j) {
            for (std::size_t i = 0; i < xCount; ++i) {
                frame.vertices.push_back({
                    100 + i + xCount * (j + yCount * k),
                    {0.5 + static_cast<double>(i),
                     -1.0 + 2.0 * static_cast<double>(j),
                     10.0 + 4.0 * static_cast<double>(k)},
                });
            }
        }
    }
    addMetadata(frame, 3, xCount, yCount, zCount);
    return frame;
}

void testDescriptorAndAxisSlices() {
    const DiagnosticFrame frame = mixedFrame();
    const auto descriptor = describeCfdGrid(frame);
    check(descriptor && descriptor->vertexBegin == 3
              && descriptor->cellCounts
                  == std::array<std::size_t, 3>{2, 3, 2}
              && descriptor->cellCount() == 12
              && descriptor->cellSpacingMetres
                  == std::array<double, 3>{1.0, 2.0, 4.0},
          "CFD slice descriptor retains the explicit X-fast grid");

    const auto ySlice = buildCfdSliceGeometry(frame, CfdSliceAxis::Y, 1);
    check(ySlice && ySlice->coordinateMetres == 1.0
              && ySlice->quads.size() == 4
              && ySlice->quads.front().sourceVertexIndex == 5
              && same(ySlice->quads.front().cornersMetres[0],
                      Vec3d{0.0, 1.0, 8.0})
              && same(ySlice->quads.front().cornersMetres[2],
                      Vec3d{1.0, 1.0, 12.0}),
          "Y slice emits one correctly sized quad per X/Z cell");

    const auto zSlice = buildCfdSliceGeometry(frame, CfdSliceAxis::Z, 1);
    check(zSlice && zSlice->coordinateMetres == 14.0
              && zSlice->quads.size() == 6
              && zSlice->quads.front().sourceVertexIndex == 9,
          "Z slice selects the requested X/Y cell layer deterministically");
}

void testMalformedMetadataIsIgnored() {
    DiagnosticFrame frame = mixedFrame();
    frame.scalarFields.pop_back();
    check(!describeCfdGrid(frame),
          "CFD slice ignores incomplete grid metadata");

    frame = mixedFrame();
    frame.scalarFields.back().values[0] = 200.0;
    check(!describeCfdGrid(frame),
          "CFD slice rejects metadata that escapes the frame vertex range");

    frame = mixedFrame();
    frame.lines.push_back({20, 3, 4, 0});
    check(!describeCfdGrid(frame),
          "CFD slice rejects a grid sample referenced by structural geometry");

    frame = mixedFrame();
    check(!buildCfdSliceGeometry(frame, CfdSliceAxis::X, 2),
          "CFD slice rejects an out-of-range plane index");
}

} // namespace

int main() {
    testDescriptorAndAxisSlices();
    testMalformedMetadataIsIgnored();
    if (failures != 0) {
        std::fprintf(stderr, "%d CFD slice check(s) failed\n", failures);
        return 1;
    }
    std::puts("all CFD slice checks passed");
    return 0;
}
