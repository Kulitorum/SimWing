#include "amr_external_flow_projection.h"

#include <cstdio>
#include <exception>

namespace {

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        simwing::fsi::amr::Runtime runtime(argc, argv);
        const auto projected =
            simwing::fsi::amr::evaluateWindTunnelProjectedCoarseGrid();
        const auto& diagnostics = projected.diagnostics;
        std::printf(
            "AMR projection: active=%zu iterations=%zu residual=%.17g "
            "divergence %.17g -> %.17g (ratio %.17g), pressure=%.17g Pa, "
            "inlet-error=%.17g m/s, outlet-change=%.17g m/s\n",
            diagnostics.activeCompositeCellCount,
            diagnostics.solverIterations,
            diagnostics.solverFinalResidual,
            diagnostics.initialMaximumDivergencePerSecond,
            diagnostics.projectedMaximumDivergencePerSecond,
            diagnostics.maximumDivergenceReductionRatio,
            diagnostics.maximumPressureCorrectionPascals,
            diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond,
            diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond);
        check(diagnostics.accepted && diagnostics.finite,
              "two-level open-boundary projection is accepted");
        check(diagnostics.hierarchy.accepted
                  && diagnostics.activeCompositeCellCount == 70'776,
              "projection retains the accepted two-level hierarchy and composite cell ownership");
        check(projected.cellCounts
                      == simwing::fsi::fluid::GridCellCounts{32, 48, 16}
                  && projected.velocityMetersPerSecond.size() == 24'576
                  && projected.pressurePascals.size() == 24'576
                  && projected.divergencePerSecond.size() == 24'576,
              "projection publishes one owning X-fast coarse diagnostic sample per cell");
        check(diagnostics.initialMaximumDivergencePerSecond > 0.0
                  && diagnostics.maximumDivergenceReductionRatio < 1.0e-7,
              "multigrid pressure correction strongly contracts composite divergence");
        check(diagnostics.maximumPressureCorrectionPascals > 0.0
                  && diagnostics.pressureOutletReferenceOwned,
              "projection owns a nonzero pressure correction and the +Y pressure reference");
        check(diagnostics.lowerYInflowNormalVelocityErrorMetersPerSecond
                      < 1.0e-11
                  && diagnostics.upperYOutletNormalVelocityChangeMetersPerSecond
                      > 0.0,
              "projection preserves prescribed -Y normal inflow and adjusts +Y outlet flux");

        const auto momentum =
            simwing::fsi::amr::evaluateWindTunnelMomentumAdvance();
        std::printf(
            "AMR momentum: CFL=%.17g velocity-change=%.17g m/s, "
            "energy %.17g -> %.17g J, predictor-divergence=%.17g 1/s, "
            "projected-divergence=%.17g 1/s\n",
            momentum.diagnostics.maximumOutgoingCourantNumber,
            momentum.diagnostics.maximumCellVelocityChangeMetersPerSecond,
            momentum.diagnostics.kineticEnergyBeforeJoules,
            momentum.diagnostics.kineticEnergyAfterJoules,
            momentum.diagnostics.correctedProjection
                .initialMaximumDivergencePerSecond,
            momentum.diagnostics.correctedProjection
                .projectedMaximumDivergencePerSecond);
        check(momentum.diagnostics.accepted
                  && momentum.diagnostics.maximumOutgoingCourantNumber < 1.0
                  && momentum.diagnostics
                         .maximumCellVelocityChangeMetersPerSecond > 0.0,
              "two-level donor-cell momentum predictor advances a nonuniform velocity within CFL");
        check(momentum.diagnostics.correctedProjection.accepted
                  && momentum.diagnostics.correctedProjection
                         .maximumDivergenceReductionRatio < 1.0e-7,
              "momentum predictor receives a fresh accepted open-boundary projection");

        simwing::fsi::amr::WindTunnelMomentumState state;
        simwing::fsi::fluid::Vector3 previousVelocity =
            state.projectedCoarseGrid().velocityMetersPerSecond[12'288];
        bool evolvedAcrossSteps = false;
        for (std::size_t step = 0; step < 12; ++step) {
            const auto advanced = state.advance();
            check(advanced.diagnostics.accepted
                      && advanced.diagnostics.correctedProjection
                             .projectedMaximumDivergencePerSecond < 1.0e-9,
                  "persistent AMR momentum state accepts each projected timestep");
            const auto currentVelocity =
                state.projectedCoarseGrid()
                    .velocityMetersPerSecond[12'288];
            evolvedAcrossSteps = evolvedAcrossSteps
                || currentVelocity != previousVelocity;
            previousVelocity = currentVelocity;
        }
        check(evolvedAcrossSteps,
              "persistent AMR momentum state evolves rather than replaying one static projection");
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr, "%d AMR projection check(s) failed\n", failures);
        return 1;
    }
    std::puts("all AMR external-flow projection checks passed");
    return 0;
}
