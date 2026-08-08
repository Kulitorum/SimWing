#include "fluid/projected_advection.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace {

using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::ProjectedMacAdvectionFailureStage;
using simwing::fsi::fluid::ProjectedMacAdvectionSspRk2Settings;
using simwing::fsi::fluid::ProjectionSettings;
using simwing::fsi::fluid::VariableMacAdvectionSettings;
using simwing::fsi::fluid::VariableMacReconstruction;
using simwing::fsi::fluid::advectVelocityByMacFlow;
using simwing::fsi::fluid::advectVelocityProjectedSspRk2;
using simwing::fsi::fluid::projectVelocity;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
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

double streamFunction(const double x, const double y) {
    return std::sin(x) * std::sin(y)
        + 0.18 * std::sin(2.0 * x + y)
        - 0.11 * std::cos(x - 2.0 * y);
}

MacVelocityField vorticalVelocity(
    const PeriodicCartesianGrid& grid) {
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    const auto lower = grid.lowerMeters();
    const auto spacing = grid.cellSpacingMeters();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            const std::size_t nextJ = j + 1 == counts.y ? 0 : j + 1;
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t nextI = i + 1 == counts.x ? 0 : i + 1;
                const double x = lower.x
                    + static_cast<double>(i) * spacing.x;
                const double nextX = lower.x
                    + static_cast<double>(nextI) * spacing.x;
                const double y = lower.y
                    + static_cast<double>(j) * spacing.y;
                const double nextY = lower.y
                    + static_cast<double>(nextJ) * spacing.y;
                const std::size_t index = grid.cellIndex(i, j, k);
                result.xFaces()[index] =
                    (streamFunction(x, nextY) - streamFunction(x, y))
                    / spacing.y;
                result.yFaces()[index] =
                    -(streamFunction(nextX, y) - streamFunction(x, y))
                    / spacing.x;
            }
        }
    }
    return result;
}

MacVelocityField translatingTaylorGreenVelocity(
    const PeriodicCartesianGrid& grid,
    const double timeSeconds) {
    constexpr double backgroundX = 0.35;
    constexpr double backgroundY = -0.2;
    MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const auto xFace = grid.xFaceCenterMeters(i, j, k);
                const auto yFace = grid.yFaceCenterMeters(i, j, k);
                result.xFaces()[index] =
                    backgroundX
                    + std::sin(xFace.x - backgroundX * timeSeconds)
                        * std::cos(xFace.y - backgroundY * timeSeconds);
                result.yFaces()[index] =
                    backgroundY
                    - std::cos(yFace.x - backgroundX * timeSeconds)
                        * std::sin(yFace.y - backgroundY * timeSeconds);
            }
        }
    }
    return result;
}

ProjectedMacAdvectionSspRk2Settings settings() {
    ProjectedMacAdvectionSspRk2Settings result;
    result.densityKgPerCubicMeter = 1.2;
    result.timeStepSeconds = 0.01;
    result.absoluteDivergenceTolerancePerSecond = 2.0e-10;
    result.relativeDivergenceTolerance = 1.0e-12;
    result.projectionAbsoluteResidualTolerance = 1.0e-12;
    result.projectionRelativeResidualTolerance = 1.0e-12;
    result.projectionMaximumIterations = 1000;
    result.absoluteMomentumToleranceNewtonSeconds = 4.0e-12;
    result.relativeMomentumTolerance = 1.0e-12;
    result.absoluteEnergyToleranceJoules = 4.0e-12;
    result.relativeEnergyTolerance = 1.0e-12;
    return result;
}

VariableMacAdvectionSettings transportSettings(
    const ProjectedMacAdvectionSspRk2Settings& source) {
    VariableMacAdvectionSettings result;
    result.densityKgPerCubicMeter = source.densityKgPerCubicMeter;
    result.timeStepSeconds = source.timeStepSeconds;
    result.reconstruction = source.reconstruction;
    result.maximumLocalOutgoingCourantNumber =
        source.maximumLocalOutgoingCourantNumber;
    result.absoluteDivergenceTolerancePerSecond =
        source.absoluteDivergenceTolerancePerSecond;
    result.relativeDivergenceTolerance = source.relativeDivergenceTolerance;
    result.absoluteMomentumToleranceNewtonSeconds =
        source.absoluteMomentumToleranceNewtonSeconds;
    result.relativeMomentumTolerance = source.relativeMomentumTolerance;
    result.absoluteEnergyToleranceJoules =
        source.absoluteEnergyToleranceJoules;
    result.relativeEnergyTolerance = source.relativeEnergyTolerance;
    if (source.reconstruction
        == VariableMacReconstruction::MonotonizedCentral) {
        result.enforceEulerEnergyNonIncrease = false;
    }
    return result;
}

ProjectionSettings pressureSettings(
    const ProjectedMacAdvectionSspRk2Settings& source) {
    ProjectionSettings result;
    result.densityKgPerCubicMeter = source.densityKgPerCubicMeter;
    result.timeStepSeconds = source.timeStepSeconds;
    result.absoluteResidualTolerance =
        source.projectionAbsoluteResidualTolerance;
    result.relativeResidualTolerance =
        source.projectionRelativeResidualTolerance;
    result.maximumIterations = source.projectionMaximumIterations;
    return result;
}

void average(const std::span<const double> original,
             const std::span<const double> twiceAdvanced,
             const std::span<double> destination) {
    for (std::size_t index = 0; index < original.size(); ++index) {
        destination[index] = original[index]
            + 0.5 * (twiceAdvanced[index] - original[index]);
    }
}

void testExactCompositionAndDeterminism() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {14, 12, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid);
    const auto integrationSettings = settings();
    const auto advection = transportSettings(integrationSettings);
    const auto projection = pressureSettings(integrationSettings);

    auto firstStage = originalVelocity;
    auto firstPressure = originalPressure;
    const auto expectedFirstAdvection = advectVelocityByMacFlow(
        grid, firstStage, firstStage, advection);
    const auto expectedFirstProjection = projectVelocity(
        grid, firstStage, firstPressure, projection);
    auto twiceAdvanced = firstStage;
    auto expectedPressure = firstPressure;
    const auto expectedSecondAdvection = advectVelocityByMacFlow(
        grid, twiceAdvanced, twiceAdvanced, advection);
    auto expectedVelocity = originalVelocity;
    average(originalVelocity.xFaces(), twiceAdvanced.xFaces(),
            expectedVelocity.xFaces());
    average(originalVelocity.yFaces(), twiceAdvanced.yFaces(),
            expectedVelocity.yFaces());
    average(originalVelocity.zFaces(), twiceAdvanced.zFaces(),
            expectedVelocity.zFaces());
    const auto expectedSecondProjection = projectVelocity(
        grid, expectedVelocity, expectedPressure, projection);

    auto firstVelocity = originalVelocity;
    auto secondVelocity = originalVelocity;
    auto actualFirstPressure = originalPressure;
    auto actualSecondPressure = originalPressure;
    const auto first = advectVelocityProjectedSspRk2(
        grid, firstVelocity, actualFirstPressure, integrationSettings);
    const auto second = advectVelocityProjectedSspRk2(
        grid, secondVelocity, actualSecondPressure, integrationSettings);
    check(expectedFirstAdvection.accepted
              && expectedFirstProjection.converged
              && expectedSecondAdvection.accepted
              && expectedSecondProjection.converged
              && first.accepted
              && first.failureStage
                  == ProjectedMacAdvectionFailureStage::None
              && firstVelocity == expectedVelocity
              && actualFirstPressure == expectedPressure,
          "projected SSPRK2: result equals its four manual stages exactly");
    check(first.firstAdvection == expectedFirstAdvection
              && first.firstProjection == expectedFirstProjection
              && first.secondAdvection == expectedSecondAdvection
              && first.secondProjection == expectedSecondProjection,
          "projected SSPRK2: every stage diagnostic retains its oracle contract");
    check(first == second
              && firstVelocity == secondVelocity
              && actualFirstPressure == actualSecondPressure,
          "projected SSPRK2: identical runs replay bit-for-bit");
    check(first.momentumResidualNormNewtonSeconds < 4.0e-12
              && first.kineticEnergyAfterJoules
                  < first.kineticEnergyBeforeJoules
              && first.totalKineticEnergyLossJoules > 0.0
              && first.finalDivergenceL2PerSecond < 2.0e-10,
          "projected SSPRK2: momentum, energy, and continuity ledgers close");
}

void testUniformNoOpAndRepeatedEligibility() {
    const PeriodicCartesianGrid uniformGrid(
        {8, 7, 3}, {}, {2.0, 3.0, 4.0});
    MacVelocityField uniform(uniformGrid);
    std::ranges::fill(uniform.xFaces(), 0.5);
    std::ranges::fill(uniform.yFaces(), -0.25);
    std::ranges::fill(uniform.zFaces(), 0.125);
    const auto uniformBefore = uniform;
    CellScalarField uniformPressure(uniformGrid);
    const auto noOp = advectVelocityProjectedSspRk2(
        uniformGrid, uniform, uniformPressure, settings());
    check(noOp.accepted && uniform == uniformBefore
              && noOp.maximumVelocityChangeMetersPerSecond == 0.0
              && noOp.totalKineticEnergyLossJoules == 0.0,
          "projected SSPRK2: uniform self-advection is a bit-exact no-op");

    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    auto firstVelocity = vorticalVelocity(grid);
    auto secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    for (std::size_t step = 0; step < 12; ++step) {
        const auto first = advectVelocityProjectedSspRk2(
            grid, firstVelocity, firstPressure, settings());
        const auto second = advectVelocityProjectedSspRk2(
            grid, secondVelocity, secondPressure, settings());
        check(first.accepted && second.accepted,
              "projected SSPRK2: every repeated nonlinear step is accepted");
        check(first == second && firstVelocity == secondVelocity
                  && firstPressure == secondPressure,
              "projected SSPRK2: repeated nonlinear steps replay bit-for-bit");
        check(first.secondAdvection.divergenceCompatible
                  && first.finalDivergenceL2PerSecond < 2.0e-10,
              "projected SSPRK2: each accepted result is eligible for the next step");
    }
}

void testLimitedReconstructionNonlinearPath() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {16, 14, 2}, {}, {twoPi, twoPi, 1.0});
    auto firstVelocity = vorticalVelocity(grid);
    auto secondVelocity = firstVelocity;
    CellScalarField firstPressure(grid);
    CellScalarField secondPressure(grid);
    auto integrationSettings = settings();
    integrationSettings.reconstruction =
        VariableMacReconstruction::MonotonizedCentral;
    integrationSettings.timeStepSeconds = 0.005;
    const auto first = advectVelocityProjectedSspRk2(
        grid, firstVelocity, firstPressure, integrationSettings);
    const auto second = advectVelocityProjectedSspRk2(
        grid, secondVelocity, secondPressure, integrationSettings);
    check(first.accepted && second.accepted
              && first.reconstruction
                  == VariableMacReconstruction::MonotonizedCentral
              && first.firstAdvection.reconstruction
                  == VariableMacReconstruction::MonotonizedCentral
              && !first.firstAdvection.energyCriterionEnabled
              && !first.secondAdvection.energyCriterionEnabled,
          "limited nonlinear path: projected SSPRK2 encloses both MUSCL Euler stages");
    check(first == second && firstVelocity == secondVelocity
              && firstPressure == secondPressure,
          "limited nonlinear path: reconstructed projected transport replays exactly");
    check(first.momentumResidualNormNewtonSeconds < 4.0e-12
              && first.kineticEnergyAfterJoules
                  <= first.kineticEnergyBeforeJoules
              && first.finalDivergenceL2PerSecond < 2.0e-10,
          "limited nonlinear path: aggregate momentum, energy, and continuity close");
}

struct FlowState {
    MacVelocityField velocity;
    CellScalarField pressure;

    explicit FlowState(const PeriodicCartesianGrid& grid)
        : velocity(vorticalVelocity(grid)), pressure(grid) {}
};

FlowState integrate(const PeriodicCartesianGrid& grid,
                    const std::size_t steps,
                    const double finalTime) {
    FlowState result(grid);
    auto integrationSettings = settings();
    integrationSettings.timeStepSeconds =
        finalTime / static_cast<double>(steps);
    integrationSettings.projectionAbsoluteResidualTolerance = 1.0e-13;
    integrationSettings.projectionRelativeResidualTolerance = 1.0e-13;
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advectVelocityProjectedSspRk2(
            grid, result.velocity, result.pressure, integrationSettings);
        check(diagnostics.accepted,
              "temporal refinement: every projected SSPRK2 step is accepted");
    }
    return result;
}

double velocityError(const MacVelocityField& actual,
                     const MacVelocityField& expected) {
    double squaredError = 0.0;
    const auto accumulate = [&squaredError](
                                const std::span<const double> first,
                                const std::span<const double> second) {
        for (std::size_t index = 0; index < first.size(); ++index) {
            const double error = first[index] - second[index];
            squaredError += error * error;
        }
    };
    accumulate(actual.xFaces(), expected.xFaces());
    accumulate(actual.yFaces(), expected.yFaces());
    accumulate(actual.zFaces(), expected.zFaces());
    return std::sqrt(
        squaredError
        / static_cast<double>(3 * actual.xFaces().size()));
}

void testObservedSecondOrderTemporalRefinement() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    constexpr double finalTime = 0.08;
    const auto reference = integrate(grid, 256, finalTime);
    const auto coarse = integrate(grid, 8, finalTime);
    const auto medium = integrate(grid, 16, finalTime);
    const auto fine = integrate(grid, 32, finalTime);
    const double coarseError = velocityError(
        coarse.velocity, reference.velocity);
    const double mediumError = velocityError(
        medium.velocity, reference.velocity);
    const double fineError = velocityError(
        fine.velocity, reference.velocity);
    const double coarseRatio = coarseError / mediumError;
    const double fineRatio = mediumError / fineError;
    check(coarseRatio > 3.6 && coarseRatio < 4.5,
          "temporal refinement: first nonlinear SSPRK2 ratio is second order");
    check(fineRatio > 3.4 && fineRatio < 4.6,
          "temporal refinement: refined nonlinear SSPRK2 ratio remains second order");
}

double translatingTaylorGreenError(
    const std::size_t resolution,
    const VariableMacReconstruction reconstruction) {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {resolution, resolution, 2}, {}, {twoPi, twoPi, 1.0});
    constexpr double finalTime = 0.08;
    const auto expected = translatingTaylorGreenVelocity(grid, finalTime);
    auto velocity = translatingTaylorGreenVelocity(grid, 0.0);
    CellScalarField pressure(grid);
    const double spacing = grid.cellSpacingMeters().x;
    // dt scales with h^2 so the SSPRK2 temporal error is fourth order in h
    // and cannot masquerade as the measured spatial convergence.
    const double requestedTimeStep = 0.12 * spacing * spacing;
    const std::size_t steps = static_cast<std::size_t>(
        std::ceil(finalTime / requestedTimeStep));
    auto integrationSettings = settings();
    integrationSettings.reconstruction = reconstruction;
    integrationSettings.timeStepSeconds =
        finalTime / static_cast<double>(steps);
    integrationSettings.projectionAbsoluteResidualTolerance = 1.0e-10;
    integrationSettings.projectionRelativeResidualTolerance = 1.0e-12;
    for (std::size_t step = 0; step < steps; ++step) {
        const auto diagnostics = advectVelocityProjectedSspRk2(
            grid, velocity, pressure, integrationSettings);
        check(diagnostics.accepted,
              "spatial refinement: every projected Taylor-Green step is accepted");
        if (!diagnostics.accepted) {
            std::fprintf(
                stderr,
                "Taylor-Green rejection n=%zu step=%zu stage=%u "
                "adv=%d/%d proj=%d/%d energy=%.17g -> %.17g "
                "momentum=%.17g div=%.17g\n",
                resolution, step,
                static_cast<unsigned int>(diagnostics.failureStage),
                diagnostics.firstAdvection.accepted ? 1 : 0,
                diagnostics.secondAdvection.accepted ? 1 : 0,
                diagnostics.firstProjection.converged ? 1 : 0,
                diagnostics.secondProjection.converged ? 1 : 0,
                diagnostics.kineticEnergyBeforeJoules,
                diagnostics.kineticEnergyAfterJoules,
                diagnostics.momentumResidualNormNewtonSeconds,
                diagnostics.finalDivergenceL2PerSecond);
            break;
        }
    }
    double absoluteError = 0.0;
    for (std::size_t index = 0; index < grid.cellCount(); ++index) {
        absoluteError += std::abs(
            velocity.xFaces()[index] - expected.xFaces()[index]);
        absoluteError += std::abs(
            velocity.yFaces()[index] - expected.yFaces()[index]);
    }
    return absoluteError
        / static_cast<double>(2 * grid.cellCount());
}

void testObservedNonlinearSpatialRefinement() {
    const double donorCoarse = translatingTaylorGreenError(
        16, VariableMacReconstruction::DonorCell);
    const double donorMedium = translatingTaylorGreenError(
        32, VariableMacReconstruction::DonorCell);
    const double donorFine = translatingTaylorGreenError(
        64, VariableMacReconstruction::DonorCell);
    const double musclCoarse = translatingTaylorGreenError(
        16, VariableMacReconstruction::MonotonizedCentral);
    const double musclMedium = translatingTaylorGreenError(
        32, VariableMacReconstruction::MonotonizedCentral);
    const double musclFine = translatingTaylorGreenError(
        64, VariableMacReconstruction::MonotonizedCentral);
    const double donorCoarseRatio = donorCoarse / donorMedium;
    const double donorFineRatio = donorMedium / donorFine;
    const double musclCoarseRatio = musclCoarse / musclMedium;
    const double musclFineRatio = musclMedium / musclFine;
    if (!(donorCoarseRatio > 1.7 && donorCoarseRatio < 2.2)
        || !(donorFineRatio > 1.7 && donorFineRatio < 2.2)
        || !(musclCoarseRatio > 3.0 && musclCoarseRatio < 5.0)
        || !(musclFineRatio > 3.0 && musclFineRatio < 5.0)
        || !(musclFine < donorFine)) {
        std::fprintf(
            stderr,
            "translating Taylor-Green spatial: donor %.17g %.17g "
            "(%.17g %.17g %.17g), MC %.17g %.17g "
            "(%.17g %.17g %.17g)\n",
            donorCoarseRatio, donorFineRatio,
            donorCoarse, donorMedium, donorFine,
            musclCoarseRatio, musclFineRatio,
            musclCoarse, musclMedium, musclFine);
    }
    check(donorCoarseRatio > 1.7 && donorCoarseRatio < 2.2,
          "spatial refinement: translating-vortex donor transport is first order");
    check(donorFineRatio > 1.7 && donorFineRatio < 2.2,
          "spatial refinement: refined translating-vortex donor transport "
          "remains first order");
    check(musclCoarseRatio > 3.0 && musclCoarseRatio < 5.0,
          "spatial refinement: first nonlinear MC ratio approaches second order");
    check(musclFineRatio > 3.0 && musclFineRatio < 5.0,
          "spatial refinement: refined nonlinear MC ratio remains near second order");
    check(musclFine < donorFine,
          "spatial refinement: limited MC improves the fine nonlinear solution");
}

void testTransactionalRejection() {
    const double twoPi = 2.0 * std::numbers::pi;
    const PeriodicCartesianGrid grid(
        {12, 12, 2}, {}, {twoPi, twoPi, 1.0});
    const auto originalVelocity = vorticalVelocity(grid);
    const CellScalarField originalPressure(grid, 0.125);

    auto unstableVelocity = originalVelocity;
    auto unstablePressure = originalPressure;
    auto unstableSettings = settings();
    unstableSettings.timeStepSeconds = 10.0;
    const auto unstable = advectVelocityProjectedSspRk2(
        grid, unstableVelocity, unstablePressure, unstableSettings);
    check(!unstable.accepted
              && unstable.failureStage
                  == ProjectedMacAdvectionFailureStage::FirstAdvection
              && unstableVelocity == originalVelocity
              && unstablePressure == originalPressure,
          "rollback: unstable first transport commits neither field");

    auto failedVelocity = originalVelocity;
    auto failedPressure = originalPressure;
    auto failedSettings = settings();
    failedSettings.projectionAbsoluteResidualTolerance = 1.0e-30;
    failedSettings.projectionRelativeResidualTolerance = 0.0;
    failedSettings.projectionMaximumIterations = 1;
    const auto failedProjection = advectVelocityProjectedSspRk2(
        grid, failedVelocity, failedPressure, failedSettings);
    check(!failedProjection.accepted
              && failedProjection.failureStage
                  == ProjectedMacAdvectionFailureStage::FirstProjection
              && failedVelocity == originalVelocity
              && failedPressure == originalPressure,
          "rollback: failed intermediate projection discards every candidate");

    auto invalidVelocity = originalVelocity;
    auto invalidPressure = originalPressure;
    auto invalidSettings = settings();
    invalidSettings.maximumLocalOutgoingCourantNumber = 1.01;
    expectRejected(
        [&] { static_cast<void>(advectVelocityProjectedSspRk2(
            grid, invalidVelocity, invalidPressure, invalidSettings)); },
        "validation: an unsafe projected SSPRK2 CFL ceiling is rejected");
    check(invalidVelocity == originalVelocity
              && invalidPressure == originalPressure,
          "validation: invalid settings mutate neither field");

    invalidVelocity = originalVelocity;
    invalidPressure = originalPressure;
    invalidSettings = settings();
    invalidSettings.reconstruction =
        static_cast<VariableMacReconstruction>(255);
    expectRejected(
        [&] { static_cast<void>(advectVelocityProjectedSspRk2(
            grid, invalidVelocity, invalidPressure, invalidSettings)); },
        "validation: an unknown projected reconstruction is rejected");
    check(invalidVelocity == originalVelocity
              && invalidPressure == originalPressure,
          "validation: rejected reconstruction mutates neither field");

    auto nonfiniteVelocity = originalVelocity;
    auto nonfinitePressure = originalPressure;
    nonfiniteVelocity.xFaces().front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(advectVelocityProjectedSspRk2(
            grid, nonfiniteVelocity, nonfinitePressure, settings())); },
        "validation: non-finite projected SSPRK2 input is rejected");
    check(std::isnan(nonfiniteVelocity.xFaces().front())
              && nonfinitePressure == originalPressure,
          "validation: non-finite rejection does not rewrite either field");
}

} // namespace

int main() {
    testExactCompositionAndDeterminism();
    testUniformNoOpAndRepeatedEligibility();
    testLimitedReconstructionNonlinearPath();
    testObservedSecondOrderTemporalRefinement();
    testObservedNonlinearSpatialRefinement();
    testTransactionalRejection();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d projected nonlinear advection check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all projected nonlinear advection checks passed");
    return 0;
}
