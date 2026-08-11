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
        const auto diagnostics =
            simwing::fsi::amr::evaluateWindTunnelPressureProjection();
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
                  && diagnostics.activeCompositeCellCount == 46'080,
              "projection retains the accepted two-level hierarchy and composite cell ownership");
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
