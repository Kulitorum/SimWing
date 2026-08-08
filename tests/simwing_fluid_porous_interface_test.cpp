#include "fluid/porous_interface.h"
#include "fluid/projection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::fsi::fluid;

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

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void testConstitutiveLawAndInverse() {
    const DarcyForchheimerResistance combined{100.0, 25.0};
    checkNear(porousPressureJumpPascals(combined, 2.0),
              -300.0, 0.0,
              "porous law opposes positive normal flow");
    checkNear(porousPressureJumpPascals(combined, -2.0),
              300.0, 0.0,
              "porous law opposes negative normal flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  combined, -300.0),
              2.0, 5.0e-16,
              "porous inverse recovers positive combined-law flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  combined, 300.0),
              -2.0, 5.0e-16,
              "porous inverse recovers negative combined-law flow");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {50.0, 0.0}, -100.0),
              2.0, 0.0,
              "porous inverse supports a pure Darcy fit");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {0.0, 25.0}, -100.0),
              2.0, 0.0,
              "porous inverse supports a pure quadratic fit");
    check(porousPressureJumpPascals(combined, 0.0) == 0.0
              && porousRelativeNormalVelocityMetersPerSecond(
                     combined, 0.0) == 0.0,
          "porous law preserves the exact zero state");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {1.0e308, 1.0e308}, -1.0e308),
              0.5 * (std::sqrt(5.0) - 1.0), 2.0e-16,
              "porous inverse remains stable at large matched scales");
    checkNear(porousRelativeNormalVelocityMetersPerSecond(
                  {0.0, 1.0e-308}, -1.0e308),
              1.0e308, 2.0e292,
              "porous quadratic inverse avoids an intermediate overflow");

    for (const double velocity : {-20.0, -1.0, -1.0e-6,
                                  1.0e-6, 1.0, 20.0}) {
        const double jump = porousPressureJumpPascals(combined, velocity);
        checkNear(porousRelativeNormalVelocityMetersPerSecond(
                      combined, jump),
                  velocity, 3.0e-15 * std::max(1.0, std::abs(velocity)),
                  "porous combined law round-trips across flow scales");
        check(-jump * velocity >= 0.0,
              "porous pressure work is dissipative");
    }
}

std::vector<PorousGridFaceCrossing> porousPlane(
    const PeriodicCartesianGrid& grid) {
    std::vector<PorousGridFaceCrossing> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back({
                10, 1, 2, GridFaceAxis::X,
                1, j, k, 0.35, 0.5, {100.0, 20.0}});
        }
    }
    return result;
}

void testFluxDrivenProjectionCanonical() {
    const PeriodicCartesianGrid grid(
        {4, 3, 2}, {}, {4.0, 3.0, 2.0});
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 2.0);
    auto definitions = porousPlane(grid);
    auto reversedDefinitions = definitions;
    std::reverse(reversedDefinitions.begin(), reversedDefinitions.end());
    const PorousPressureJumpField porous(
        grid, velocity, definitions);
    const PorousPressureJumpField reversed(
        grid, velocity, reversedDefinitions);

    check(porous.pressureJumps() == reversed.pressureJumps()
              && std::ranges::equal(porous.samples(), reversed.samples()),
          "porous field canonicalizes reversed authored crossings");
    check(porous.samples().size() == 6
              && porous.pressureJumps().faceCount() == 6,
          "porous field retains every resolved plane tile");
    for (const auto& sample : porous.samples()) {
        checkNear(sample.fluidNormalVelocityMetersPerSecond,
                  2.0, 0.0,
                  "porous sample reads the owning MAC normal velocity");
        checkNear(sample.relativeNormalVelocityMetersPerSecond,
                  1.5, 0.0,
                  "porous sample subtracts sheet normal velocity");
        checkNear(sample.pressureJump.pressureJumpPascals,
                  -195.0, 0.0,
                  "porous sample applies the calibrated signed jump");
        checkNear(sample.volumeFlowRateCubicMetersPerSecond,
                  1.5, 0.0,
                  "porous sample retains resolved tile volume flux");
        checkNear(sample.dissipationWatts,
                  292.5, 0.0,
                  "porous sample closes positive pressure dissipation");
    }
    checkNear(porous.totalDissipationWatts(),
              1755.0, 0.0,
              "porous field aggregates tile dissipation");

    std::vector<GridFacePressureJump> drivenJumps(
        porous.pressureJumps().faces().begin(),
        porous.pressureJumps().faces().end());
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            drivenJumps.push_back({
                20, 2, 1, GridFaceAxis::X,
                3, j, k, 195.0, 0.5});
        }
    }
    const SharpPressureJumpField completeJumps(
        grid, std::move(drivenJumps));
    CellScalarField pressure(grid);
    ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.absoluteResidualTolerance = 1.0e-11;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    const auto diagnostics = projectVelocityWithPressureJumps(
        grid, velocity, pressure, completeJumps, settings);
    check(diagnostics.converged
              && diagnostics.pressureJumpFaceCount == 12
              && diagnostics.divergenceL2AfterPerSecond < 1.0e-13,
          "porous plane plus pressure source projects a compatible flow");
    check(std::ranges::all_of(
              velocity.xFaces(),
              [](const double sample) { return sample == 2.0; })
              && std::ranges::all_of(
                  velocity.yFaces(),
                  [](const double sample) { return sample == 0.0; })
              && std::ranges::all_of(
                  velocity.zFaces(),
                  [](const double sample) { return sample == 0.0; }),
          "flux-driven porous canonical creates no spurious velocity");
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i == 0 || i == 3
                    ? 97.5 : -97.5;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(pressure.values()[grid.cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 2.0e-12,
          "flux-driven porous canonical preserves its analytic pressure loss");
}

void testAxisAreaAndValidation() {
    const PeriodicCartesianGrid grid(
        {2, 2, 2}, {}, {4.0, 6.0, 8.0});
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 3.0);
    std::ranges::fill(velocity.yFaces(), 4.0);
    std::ranges::fill(velocity.zFaces(), 5.0);
    const DarcyForchheimerResistance linear{10.0, 0.0};
    const PorousPressureJumpField axes(grid, velocity, {
        {10, 1, 2, GridFaceAxis::X, 0, 0, 0, 0.2, 0.0, linear},
        {20, 3, 4, GridFaceAxis::Y, 1, 0, 0, 0.4, 0.0, linear},
        {30, 5, 6, GridFaceAxis::Z, 1, 1, 0, 0.6, 0.0, linear},
    });
    check(axes.samples().size() == 3
              && axes.samples()[0].faceAreaSquareMeters == 12.0
              && axes.samples()[1].faceAreaSquareMeters == 8.0
              && axes.samples()[2].faceAreaSquareMeters == 6.0
              && axes.samples()[0].volumeFlowRateCubicMetersPerSecond == 36.0
              && axes.samples()[1].volumeFlowRateCubicMetersPerSecond == 32.0
              && axes.samples()[2].volumeFlowRateCubicMetersPerSecond == 30.0,
          "porous sampling uses the correct X/Y/Z MAC tile areas");

    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals({}, 1.0)); },
        "porous law rejects zero resistance");
    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals(
            {-1.0, 1.0}, 1.0)); },
        "porous law rejects negative resistance");
    expectInvalid(
        [&] { static_cast<void>(porousPressureJumpPascals(
            linear, std::numeric_limits<double>::infinity())); },
        "porous law rejects non-finite velocity");
    expectInvalid(
        [&] { static_cast<void>(
            porousRelativeNormalVelocityMetersPerSecond(
                linear, std::numeric_limits<double>::quiet_NaN())); },
        "porous inverse rejects non-finite pressure");

    const PeriodicCartesianGrid otherGrid(
        {3, 2, 2}, {}, {6.0, 6.0, 8.0});
    const MacVelocityField foreignVelocity(otherGrid);
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, foreignVelocity, {})); },
        "porous field rejects foreign MAC topology");
    auto invalid = PorousGridFaceCrossing{
        10, 1, 2, GridFaceAxis::X, 2, 0, 0,
        0.5, 0.0, linear};
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field rejects an out-of-range face");
    invalid.i = 0;
    invalid.surfaceNormalVelocityMetersPerSecond =
        std::numeric_limits<double>::infinity();
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field rejects non-finite sheet velocity");
    invalid.surfaceNormalVelocityMetersPerSecond = 0.0;
    invalid.surfaceStableId = 0;
    expectInvalid(
        [&] { static_cast<void>(PorousPressureJumpField(
            grid, velocity, {invalid})); },
        "porous field delegates stable topology validation");
}

} // namespace

int main() {
    testConstitutiveLawAndInverse();
    testFluxDrivenProjectionCanonical();
    testAxisAreaAndValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing porous-interface check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing porous-interface checks passed");
    return 0;
}
