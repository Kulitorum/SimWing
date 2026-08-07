#include "fluid/interface_jump.h"
#include "fluid/projection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFacePressureJump;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::ProjectionSettings;
using simwing::fsi::fluid::SharpPressureJumpField;
using simwing::fsi::fluid::computePressureGradient;
using simwing::fsi::fluid::computePressureGradientWithJumps;
using simwing::fsi::fluid::computePressureJumpSource;
using simwing::fsi::fluid::maximumAbsoluteValue;
using simwing::fsi::fluid::mean;
using simwing::fsi::fluid::projectVelocity;
using simwing::fsi::fluid::projectVelocityWithPressureJumps;

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
                     "FAIL: %s (actual %.17g, expected %.17g, tolerance %.3g)\n",
                     message, actual, expected, tolerance);
        ++failures;
    }
}

PeriodicCartesianGrid testGrid() {
    return PeriodicCartesianGrid({16, 4, 3}, {}, {1.0, 1.0, 1.0});
}

ProjectionSettings strictSettings() {
    ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.025;
    settings.absoluteResidualTolerance = 1.0e-9;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    return settings;
}

std::vector<GridFacePressureJump> slabFaces(
    const PeriodicCartesianGrid& grid,
    const std::size_t firstInsideCell,
    const std::size_t firstOutsideCell,
    const double insideMinusOutsidePascals) {
    std::vector<GridFacePressureJump> faces;
    const auto counts = grid.cellCounts();
    faces.reserve(2 * counts.y * counts.z);
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                100, 1, 2, GridFaceAxis::X,
                firstInsideCell, j, k, insideMinusOutsidePascals});
            faces.push_back({
                101, 2, 1, GridFaceAxis::X,
                firstOutsideCell, j, k, -insideMinusOutsidePascals});
        }
    }
    return faces;
}

SharpPressureJumpField staticSlabJumps(
    const PeriodicCartesianGrid& grid,
    const double jumpPascals = 250.0) {
    return SharpPressureJumpField(grid, slabFaces(grid, 4, 11, jumpPascals));
}

void expectRejected(const PeriodicCartesianGrid& grid,
                    std::vector<GridFacePressureJump> faces,
                    const char* message) {
    bool rejected = false;
    try {
        static_cast<void>(SharpPressureJumpField(grid, std::move(faces)));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void testValidationAndCanonicalization() {
    const PeriodicCartesianGrid grid({4, 3, 2}, {}, {1.0, 1.0, 1.0});
    std::vector<GridFacePressureJump> authored = {
        {20, 2, 1, GridFaceAxis::Z, 3, 2, 1, -4.0},
        {10, 1, 2, GridFaceAxis::X, 1, 0, 0, 4.0},
        {10, 2, 1, GridFaceAxis::Y, 2, 1, 0, -4.0},
    };
    auto reversed = authored;
    std::reverse(reversed.begin(), reversed.end());
    const SharpPressureJumpField first(grid, authored);
    const SharpPressureJumpField second(grid, reversed);
    check(first == second,
          "validation: authored order canonicalizes to an identical jump field");
    check(first.faceCount() == 3 && !first.empty(),
          "validation: every unique crossing is retained");
    check(first.faces()[0].axis == GridFaceAxis::X
              && first.faces()[1].axis == GridFaceAxis::Y
              && first.faces()[2].axis == GridFaceAxis::Z,
          "validation: canonical order is axis then face index");
    checkNear(first.xFaceJumpsPascals()[grid.cellIndex(1, 0, 0)],
              4.0, 0.0,
              "validation: signed x jump is stored on its unique face");

    auto invalid = authored;
    invalid[0].surfaceStableId = 0;
    expectRejected(grid, invalid,
                   "validation: zero surface stable ID is rejected");
    invalid = authored;
    invalid[0].plusRegionStableId = invalid[0].minusRegionStableId;
    expectRejected(grid, invalid,
                   "validation: equal side regions are rejected");
    invalid = authored;
    invalid[0].pressureJumpPascals =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(grid, invalid,
                   "validation: non-finite pressure jump is rejected");
    invalid = authored;
    invalid[0].i = grid.cellCounts().x;
    expectRejected(grid, invalid,
                   "validation: out-of-range face coordinate is rejected");
    invalid = authored;
    invalid[0].axis = static_cast<GridFaceAxis>(99);
    expectRejected(grid, invalid,
                   "validation: unknown face axis is rejected");
    invalid = authored;
    invalid.push_back({30, 1, 2, GridFaceAxis::X, 1, 0, 0, 7.0});
    expectRejected(grid, invalid,
                   "validation: multiple crossings on one face are explicit unsupported topology");
    invalid = authored;
    invalid.push_back({10, 3, 4, GridFaceAxis::X, 2, 0, 0, 1.0});
    expectRejected(grid, invalid,
                   "validation: one surface ID cannot alias different region pairs");
}

void testSharpGradientAndSourcePairing() {
    const auto grid = testGrid();
    constexpr double jumpPascals = 250.0;
    const auto jumps = staticSlabJumps(grid, jumpPascals);
    CellScalarField pressure(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                pressure.values()[grid.cellIndex(i, j, k)] =
                    i >= 4 && i < 11 ? jumpPascals : 0.0;
            }
        }
    }
    MacVelocityField ordinaryGradient(grid);
    MacVelocityField sharpGradient(grid);
    computePressureGradient(grid, pressure, ordinaryGradient);
    computePressureGradientWithJumps(
        grid, pressure, jumps, sharpGradient);
    double ordinaryMaximum = 0.0;
    double sharpMaximum = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        ordinaryMaximum = std::max({
            ordinaryMaximum,
            std::abs(ordinaryGradient.xFaces()[index]),
            std::abs(ordinaryGradient.yFaces()[index]),
            std::abs(ordinaryGradient.zFaces()[index]),
        });
        sharpMaximum = std::max({
            sharpMaximum,
            std::abs(sharpGradient.xFaces()[index]),
            std::abs(sharpGradient.yFaces()[index]),
            std::abs(sharpGradient.zFaces()[index]),
        });
    }
    checkNear(ordinaryMaximum,
              jumpPascals / grid.cellSpacingMeters().x,
              0.0,
              "jump operator: ordinary gradient sees the discontinuity");
    checkNear(sharpMaximum, 0.0, 0.0,
              "jump operator: prescribed discontinuity has zero sharp gradient");

    CellScalarField source(grid);
    computePressureJumpSource(grid, jumps, source);
    checkNear(mean(source), 0.0, 0.0,
              "jump operator: periodic pressure source is exactly compatible");
    check(maximumAbsoluteValue(source) > 0.0,
          "jump operator: membrane discontinuity contributes a Poisson source");
}

void testAllAxisJumpStencilsAndPeriodicMinusCells() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {2.0, 3.0, 4.0});
    const std::array<GridFacePressureJump, 3> cases = {{
        {10, 1, 2, GridFaceAxis::X, 0, 1, 1, 6.0},
        {20, 1, 2, GridFaceAxis::Y, 2, 0, 1, 6.0},
        {30, 1, 2, GridFaceAxis::Z, 3, 2, 0, 6.0},
    }};
    const auto counts = grid.cellCounts();
    const auto spacing = grid.cellSpacingMeters();
    for (const auto& face : cases) {
        const SharpPressureJumpField jumps(grid, {face});
        CellScalarField zeroPressure(grid);
        MacVelocityField sharpGradient(grid);
        computePressureGradientWithJumps(
            grid, zeroPressure, jumps, sharpGradient);
        const auto plusCell = grid.cellIndex(face.i, face.j, face.k);
        std::size_t minusCell = 0;
        double spacingMeters = 0.0;
        std::span<const double> selectedGradient;
        switch (face.axis) {
        case GridFaceAxis::X:
            minusCell = grid.cellIndex(counts.x - 1, face.j, face.k);
            spacingMeters = spacing.x;
            selectedGradient = sharpGradient.xFaces();
            break;
        case GridFaceAxis::Y:
            minusCell = grid.cellIndex(face.i, counts.y - 1, face.k);
            spacingMeters = spacing.y;
            selectedGradient = sharpGradient.yFaces();
            break;
        case GridFaceAxis::Z:
            minusCell = grid.cellIndex(face.i, face.j, counts.z - 1);
            spacingMeters = spacing.z;
            selectedGradient = sharpGradient.zFaces();
            break;
        }
        checkNear(selectedGradient[plusCell],
                  -face.pressureJumpPascals / spacingMeters,
                  0.0,
                  "jump stencil: sharp gradient uses the signed face jump");

        CellScalarField source(grid, 123.0);
        computePressureJumpSource(grid, jumps, source);
        const double sourceMagnitude = face.pressureJumpPascals
            / (spacingMeters * spacingMeters);
        checkNear(source.values()[minusCell], sourceMagnitude, 0.0,
                  "jump stencil: periodic minus cell receives positive source");
        checkNear(source.values()[plusCell], -sourceMagnitude, 0.0,
                  "jump stencil: plus cell receives negative source");
        double absoluteSourceSum = 0.0;
        for (const double value : source.values()) {
            absoluteSourceSum += std::abs(value);
        }
        checkNear(absoluteSourceSum, 2.0 * sourceMagnitude, 0.0,
                  "jump stencil: source overwrites output and touches only adjacent cells");
    }
}

void testStaticPressureJumpProjection() {
    const auto grid = testGrid();
    constexpr double jumpPascals = 250.0;
    const auto jumps = staticSlabJumps(grid, jumpPascals);
    MacVelocityField firstVelocity(grid);
    MacVelocityField secondVelocity(grid);
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    const auto settings = strictSettings();
    const auto firstDiagnostics = projectVelocityWithPressureJumps(
        grid, firstVelocity, firstPressure, jumps, settings);
    const auto secondDiagnostics = projectVelocityWithPressureJumps(
        grid, secondVelocity, secondPressure, jumps, settings);

    check(firstDiagnostics.converged,
          "static jump: sharp pressure projection converges");
    check(firstDiagnostics.pressureJumpFaceCount == 24,
          "static jump: diagnostics report every membrane face crossing");
    checkNear(firstDiagnostics.pressureJumpSourceCompatibilityPascalsPerSquareMeter,
              0.0, 0.0,
              "static jump: the closed periodic source is exactly compatible");
    check(firstDiagnostics.divergenceL2AfterPerSecond < 1.0e-11,
          "static jump: projected velocity remains divergence-free");
    check(firstVelocity == secondVelocity
              && firstPressure == secondPressure
              && firstDiagnostics == secondDiagnostics,
          "static jump: identical sharp projections replay bit-for-bit");

    const double insideFraction = 7.0 / 16.0;
    const double expectedOutside = -jumpPascals * insideFraction;
    const double expectedInside = jumpPascals + expectedOutside;
    double maximumPressureError = 0.0;
    double maximumVelocity = 0.0;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const double expected = i >= 4 && i < 11
                    ? expectedInside : expectedOutside;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(firstPressure.values()[index] - expected));
                maximumVelocity = std::max({
                    maximumVelocity,
                    std::abs(firstVelocity.xFaces()[index]),
                    std::abs(firstVelocity.yFaces()[index]),
                    std::abs(firstVelocity.zFaces()[index]),
                });
            }
        }
    }
    check(maximumPressureError < 2.0e-12,
          "static jump: cell pressure remains a sharp two-level field");
    check(maximumVelocity < 2.0e-13,
          "static jump: balanced pressure discontinuity creates no spurious flow");
    checkNear(mean(firstPressure), 0.0, 2.0e-14,
              "static jump: pressure retains the periodic zero-mean gauge");
}

void testPressureJumpAcrossPeriodicBoundary() {
    const auto grid = testGrid();
    constexpr double jumpPascals = 80.0;
    const SharpPressureJumpField jumps(
        grid, slabFaces(grid, 0, 4, jumpPascals));
    checkNear(jumps.xFaceJumpsPascals()[grid.cellIndex(0, 0, 0)],
              jumpPascals, 0.0,
              "periodic jump: face zero separates the last and first cells");

    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const auto diagnostics = projectVelocityWithPressureJumps(
        grid, velocity, pressure, jumps, strictSettings());
    check(diagnostics.converged,
          "periodic jump: boundary-crossing membrane projection converges");
    const double expectedOutside = -0.25 * jumpPascals;
    const double expectedInside = 0.75 * jumpPascals;
    double maximumPressureError = 0.0;
    double maximumVelocity = 0.0;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto index = grid.cellIndex(i, j, k);
                const double expected = i < 4
                    ? expectedInside : expectedOutside;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(pressure.values()[index] - expected));
                maximumVelocity = std::max({
                    maximumVelocity,
                    std::abs(velocity.xFaces()[index]),
                    std::abs(velocity.yFaces()[index]),
                    std::abs(velocity.zFaces()[index]),
                });
            }
        }
    }
    check(maximumPressureError < 1.0e-12,
          "periodic jump: pressure discontinuity is sharp across the domain boundary");
    check(maximumVelocity < 1.0e-13,
          "periodic jump: boundary placement creates no spurious flow");
}

void fillDeterministicVelocity(MacVelocityField& velocity) {
    for (std::size_t index = 0; index < velocity.xFaces().size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        velocity.xFaces()[index] = std::sin(0.31 * sample);
        velocity.yFaces()[index] = std::cos(0.23 * sample);
        velocity.zFaces()[index] = std::sin(0.17 * sample + 0.2);
    }
}

void testEmptyEquivalenceMismatchAndRollback() {
    const auto grid = testGrid();
    SharpPressureJumpField emptyJumps(grid);
    MacVelocityField directVelocity(grid);
    fillDeterministicVelocity(directVelocity);
    MacVelocityField jumpVelocity = directVelocity;
    CellScalarField directPressure(grid);
    CellScalarField jumpPressure(grid);
    const auto settings = strictSettings();
    const auto directDiagnostics = projectVelocity(
        grid, directVelocity, directPressure, settings);
    const auto jumpDiagnostics = projectVelocityWithPressureJumps(
        grid, jumpVelocity, jumpPressure, emptyJumps, settings);
    check(directVelocity == jumpVelocity
              && directPressure == jumpPressure
              && directDiagnostics == jumpDiagnostics,
          "empty jump: overload takes the exact no-interface arithmetic path");

    const PeriodicCartesianGrid otherGrid({8, 4, 3}, {}, {1.0, 1.0, 1.0});
    SharpPressureJumpField foreignJumps(otherGrid);
    bool rejected = false;
    try {
        static_cast<void>(projectVelocityWithPressureJumps(
            grid, jumpVelocity, jumpPressure, foreignJumps, settings));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected,
          "validation: even an empty jump field is bound to its grid shape");

    const auto jumps = staticSlabJumps(grid);
    MacVelocityField failedVelocity(grid);
    CellScalarField failedPressure(grid, 3.0);
    const auto originalVelocity = failedVelocity;
    const auto originalPressure = failedPressure;
    auto failingSettings = settings;
    failingSettings.absoluteResidualTolerance = 1.0e-30;
    failingSettings.relativeResidualTolerance = 0.0;
    failingSettings.maximumIterations = 1;
    const auto failedDiagnostics = projectVelocityWithPressureJumps(
        grid, failedVelocity, failedPressure, jumps, failingSettings);
    check(!failedDiagnostics.converged,
          "rollback: intentionally truncated sharp projection does not converge");
    check(failedVelocity == originalVelocity && failedPressure == originalPressure,
          "rollback: failed sharp projection commits neither pressure nor velocity");
}

} // namespace

int main() {
    testValidationAndCanonicalization();
    testSharpGradientAndSourcePairing();
    testAllAxisJumpStencilsAndPeriodicMinusCells();
    testStaticPressureJumpProjection();
    testPressureJumpAcrossPeriodicBoundary();
    testEmptyEquivalenceMismatchAndRollback();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid interface-jump check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid interface-jump checks passed");
    return 0;
}
