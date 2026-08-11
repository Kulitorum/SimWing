#include "amr_external_flow.h"

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
            simwing::fsi::amr::evaluateWindTunnelBoundaryInitialization();
        check(diagnostics.accepted && diagnostics.finite,
              "AMR wind-tunnel boundary initialization is accepted");
        check(diagnostics.coarseCellCounts
                  == simwing::fsi::fluid::GridCellCounts{32, 48, 16}
                  && diagnostics.refinedDomainCellCounts
                      == simwing::fsi::fluid::GridCellCounts{64, 96, 32}
                  && diagnostics.coarseCellSpacingMeters
                      == simwing::fsi::fluid::Vector3{0.375, 0.25, 0.375},
              "AMR wind-tunnel hierarchy retains explicit physical spacing");
        check(diagnostics.physicalBoundaries
                  == simwing::fsi::amr::positiveYWindTunnelBoundaries
                  && diagnostics.allPhysicalDirectionsNonPeriodic,
              "positive-Y wind tunnel owns inflow, outlet, and far field without periodic faces");
        check(diagnostics.coarseBlockCount == 6
                  && diagnostics.refinedBlockCount > 0
                  && diagnostics.coarseValidCellCount == 24'576
                  && diagnostics.refinedValidCellCount == 24'576
                  && diagnostics.allocatedVelocityBytes > 0,
              "AMR hierarchy allocates bounded coarse and refined face fields");
        check(diagnostics.maximumInteriorWakePerturbationMetersPerSecond > 0.0
                  && diagnostics.lowerYInflowMaximumErrorMetersPerSecond
                      == 0.0
                  && diagnostics
                         .upperYOutflowGradientMaximumErrorMetersPerSecond
                      == 0.0
                  && diagnostics.farFieldMaximumErrorMetersPerSecond == 0.0,
              "inflow, pressure-outlet velocity gradient, and far-field ghost contracts close exactly");
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr, "%d AMR external-flow check(s) failed\n", failures);
        return 1;
    }
    std::puts("all AMR external-flow checks passed");
    return 0;
}
