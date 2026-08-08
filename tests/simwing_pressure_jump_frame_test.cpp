#include "pressure_jump_frame.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing;
namespace fluid = simwing::fsi::fluid;

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
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [&](const auto& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

const viewer::VectorField* vectorField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.vectorFields,
        [&](const auto& field) { return field.name == name; });
    return found == frame.vectorFields.end() ? nullptr : &*found;
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "pressure-jump frame serializes");
    return bytes;
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

fluid::ProjectionSettings strictSettings() {
    fluid::ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.absoluteResidualTolerance = 1.0e-11;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    return settings;
}

viewer::PressureJumpFrameContext context() {
    return {
        "sha256:pressure-jump-frame-test",
        "simwing-pressure-jump-frame-test-v1",
        1,
        0.01,
        0.01,
        1.2,
    };
}

std::vector<fluid::GridFacePressureJump> splitSlab(
    const fluid::PeriodicCartesianGrid& grid) {
    std::vector<fluid::GridFacePressureJump> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                10, 1, 2, fluid::GridFaceAxis::X,
                1, j, k, 40.0, 0.25});
            result.push_back({
                20, 2, 3, fluid::GridFaceAxis::X,
                1, j, k, 60.0, 0.75});
            result.push_back({
                30, 3, 2, fluid::GridFaceAxis::X,
                3, j, k, -60.0, 0.25});
            result.push_back({
                40, 2, 1, fluid::GridFaceAxis::X,
                3, j, k, -40.0, 0.75});
        }
    }
    return result;
}

void testOwningLayerGeometryAndDeterminism() {
    const fluid::PeriodicCartesianGrid grid(
        {4, 2, 2}, {}, {4.0, 2.0, 2.0});
    auto authored = splitSlab(grid);
    auto reversed = authored;
    std::reverse(reversed.begin(), reversed.end());
    const fluid::SharpPressureJumpField firstJumps(grid, authored);
    const fluid::SharpPressureJumpField secondJumps(grid, reversed);
    fluid::MacVelocityField firstVelocity(grid);
    fluid::MacVelocityField secondVelocity(grid);
    fluid::CellScalarField firstPressure(grid);
    fluid::CellScalarField secondPressure(grid);
    const auto settings = strictSettings();
    const auto firstDiagnostics = fluid::projectVelocityWithPressureJumps(
        grid, firstVelocity, firstPressure, firstJumps, settings);
    const auto secondDiagnostics = fluid::projectVelocityWithPressureJumps(
        grid, secondVelocity, secondPressure, secondJumps, settings);
    const auto first = viewer::buildPressureJumpFrame(
        grid, firstVelocity, firstPressure,
        firstJumps, firstDiagnostics, context());
    const auto second = viewer::buildPressureJumpFrame(
        grid, secondVelocity, secondPressure,
        secondJumps, secondDiagnostics, context());

    check(first.vertices.size() == grid.cellCount() + 4 * 16
              && first.triangles.size() == 2 * 16,
          "pressure-jump frame owns cell samples and one quad per crossing");
    check(serialized(first) == serialized(second),
          "pressure-jump frame is byte deterministic after canonicalization");
    const std::size_t firstSurfaceVertex = grid.cellCount();
    checkNear(first.vertices[firstSurfaceVertex].positionMetres.x,
              0.75, 0.0,
              "pressure-jump frame places the first layer at its fraction");
    checkNear(first.vertices[firstSurfaceVertex + 4].positionMetres.x,
              1.25, 0.0,
              "pressure-jump frame separates a second same-face layer");
    const auto& triangle = first.triangles.front();
    const auto& a = first.vertices[triangle.vertex0].positionMetres;
    const auto& b = first.vertices[triangle.vertex1].positionMetres;
    const auto& c = first.vertices[triangle.vertex2].positionMetres;
    const double orientedX = (b.y - a.y) * (c.z - a.z)
        - (b.z - a.z) * (c.y - a.y);
    check(orientedX > 0.0
              && triangle.negativeRegionId == 1
              && triangle.positiveRegionId == 2,
          "pressure-jump quad orientation follows minus-to-plus X traversal");

    const auto* pressure = scalarField(first, "pressure sample");
    const auto* jump = scalarField(first, "pressure jump");
    const auto* fraction = scalarField(first, "crossing fraction");
    const auto* normal = vectorField(first, "interface normal");
    check(pressure != nullptr
              && pressure->values.size() == first.vertices.size()
              && jump != nullptr
              && jump->values.size() == first.triangles.size()
              && fraction != nullptr
              && fraction->values.size() == first.triangles.size()
              && normal != nullptr
              && normal->values.size() == first.triangles.size(),
          "pressure-jump frame fields match their entity associations");
    checkNear(
        pressure->values[firstSurfaceVertex + 4]
            - pressure->values[firstSurfaceVertex],
        50.0, 2.0e-13,
        "pressure-jump frame reconstructs ordered intermediate pressure");
    check(jump->values[0] == 40.0
              && jump->values[1] == 40.0
              && jump->values[2] == 60.0
              && fraction->values[0] == 0.25
              && fraction->values[2] == 0.75
              && normal->values[0].x == 1.0,
          "pressure-jump triangle fields preserve each authored crossing");
}

void testPeriodicWrapAndAllAxisOrientation() {
    const fluid::PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    const fluid::SharpPressureJumpField jumps(grid, {
        {10, 1, 2, fluid::GridFaceAxis::X,
         0, 0, 0, 1.0, 0.25},
        {20, 3, 4, fluid::GridFaceAxis::Y,
         2, 1, 0, 2.0, 0.5},
        {30, 5, 6, fluid::GridFaceAxis::Z,
         3, 2, 1, 3.0, 0.75},
    });
    fluid::MacVelocityField velocity(grid);
    fluid::CellScalarField pressure(grid);
    const auto diagnostics = fluid::projectVelocityWithPressureJumps(
        grid, velocity, pressure, jumps, strictSettings());
    const auto frame = viewer::buildPressureJumpFrame(
        grid, velocity, pressure, jumps, diagnostics, context());
    const std::size_t firstSurfaceVertex = grid.cellCount();
    checkNear(frame.vertices[firstSurfaceVertex].positionMetres.x,
              3.75, 0.0,
              "pressure-jump frame wraps a periodic pre-face crossing");
    const auto* normals = vectorField(frame, "interface normal");
    check(normals != nullptr
              && normals->values[0].x == 1.0
              && normals->values[2].y == 1.0
              && normals->values[4].z == 1.0,
          "pressure-jump frame preserves positive X/Y/Z crossing normals");
}

void testRejectedInputs() {
    const fluid::PeriodicCartesianGrid grid(
        {4, 2, 2}, {}, {4.0, 2.0, 2.0});
    const fluid::SharpPressureJumpField jumps(grid, splitSlab(grid));
    fluid::MacVelocityField velocity(grid);
    fluid::CellScalarField pressure(grid);
    const auto diagnostics = fluid::projectVelocityWithPressureJumps(
        grid, velocity, pressure, jumps, strictSettings());
    auto invalidDiagnostics = diagnostics;
    invalidDiagnostics.pressureJumpFaceCount = 0;
    expectRejected(
        [&] { static_cast<void>(viewer::buildPressureJumpFrame(
            grid, velocity, pressure, jumps,
            invalidDiagnostics, context())); },
        "pressure-jump frame rejects inconsistent accepted diagnostics");
    auto invalidContext = context();
    invalidContext.step = 0;
    expectRejected(
        [&] { static_cast<void>(viewer::buildPressureJumpFrame(
            grid, velocity, pressure, jumps,
            diagnostics, invalidContext)); },
        "pressure-jump frame rejects an uncommitted context");
    auto invalidPressure = pressure;
    invalidPressure.values()[0] =
        std::numeric_limits<double>::infinity();
    expectRejected(
        [&] { static_cast<void>(viewer::buildPressureJumpFrame(
            grid, velocity, invalidPressure, jumps,
            diagnostics, context())); },
        "pressure-jump frame rejects non-finite solver state");
    const fluid::PeriodicCartesianGrid otherGrid(
        {3, 2, 2}, {}, {3.0, 2.0, 2.0});
    const fluid::SharpPressureJumpField foreign(otherGrid);
    expectRejected(
        [&] { static_cast<void>(viewer::buildPressureJumpFrame(
            grid, velocity, pressure, foreign,
            diagnostics, context())); },
        "pressure-jump frame rejects foreign crossing topology");
}

} // namespace

int main() {
    testOwningLayerGeometryAndDeterminism();
    testPeriodicWrapAndAllAxisOrientation();
    testRejectedInputs();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing pressure-jump frame check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing pressure-jump frame checks passed");
    return 0;
}
