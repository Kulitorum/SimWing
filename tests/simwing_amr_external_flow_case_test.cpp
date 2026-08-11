#include "amr_external_flow_case.h"
#include "cfd_slice.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace {

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

const simwing::viewer::ScalarField* field(
    const simwing::viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [&](const simwing::viewer::ScalarField& candidate) {
            return candidate.name == name;
        });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

simwing::fsi::Scene staticPlateScene() {
    using namespace simwing::fsi;
    Scene scene;
    scene.metadata.designChecksum = "sha256:amr-static-plate";
    scene.metadata.exporterVersion = "amr-static-plate-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "inside"},
    };
    scene.vertices = {
        {0x1000000000000010ULL, {-1.0, 0.93, -1.0}},
        {0x1000000000000011ULL, {1.0, 0.93, -1.0}},
        {0x1000000000000012ULL, {1.0, 0.93, 1.0}},
        {0x1000000000000013ULL, {-1.0, 0.93, 1.0}},
    };
    scene.fabricMaterials = {
        {100, "ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {0x2000000000000500ULL,
         {0x1000000000000010ULL, 0x1000000000000011ULL,
          0x1000000000000012ULL},
         {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {0x2000000000000501ULL,
         {0x1000000000000010ULL, 0x1000000000000012ULL,
          0x1000000000000013ULL},
         {{{0.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    return scene;
}

std::vector<std::uint8_t> serialized(
    const simwing::viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    simwing::viewer::ProtocolError error;
    check(simwing::viewer::serializeFrame(frame, bytes, &error),
          "AMR external-flow frame serializes");
    return bytes;
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        simwing::fsi::amr::Runtime runtime(argc, argv);
        simwing::fsi::amr::ExternalFlowTransportCase first;
        simwing::fsi::amr::ExternalFlowTransportCase second;
        simwing::viewer::DiagnosticFrame firstFrame;
        simwing::viewer::DiagnosticFrame secondFrame;
        for (std::size_t step = 0; step < 8; ++step) {
            firstFrame = first.advance();
            secondFrame = second.advance();
            check(serialized(firstFrame) == serialized(secondFrame),
                  "AMR external-flow transport replays deterministically");
        }

        simwing::viewer::ProtocolError error;
        check(simwing::viewer::validateFrame(firstFrame, &error)
                  && firstFrame.step == 8
                  && firstFrame.vertices.size() == 24'576,
              "AMR external-flow case publishes an accepted coarse-grid frame");
        const auto descriptor = simwing::viewer::describeCfdGrid(firstFrame);
        check(descriptor
                  && descriptor->cellCounts
                      == std::array<std::size_t, 3>{32, 48, 16}
                  && descriptor->cellSpacingMetres
                      == std::array<double, 3>{0.375, 0.25, 0.375},
              "AMR external-flow frame drives the axis-aligned CFD slice controls");
        const auto* marker = field(firstFrame, "external_flow.marker");
        bool markerMovedDownstream = false;
        if (marker != nullptr) {
            for (std::size_t k = 0; k < 16; ++k) {
                for (std::size_t j = 1; j < 48; ++j) {
                    for (std::size_t i = 0; i < 32; ++i) {
                        if (marker->values[i + 32 * (j + 48 * k)] > 0.0) {
                            markerMovedDownstream = true;
                        }
                    }
                }
            }
        }
        check(marker != nullptr && markerMovedDownstream,
              "passive flow marker visibly advances in positive Y");
        check(first.diagnostics().accepted
                  && first.diagnostics().momentum.accepted
                  && first.diagnostics().momentum
                         .maximumCellVelocityChangeMetersPerSecond > 0.0
                  && first.diagnostics().maximumOutgoingCourantNumber < 1.0
                  && first.diagnostics().minimumMarker >= 0.0
                  && first.diagnostics().maximumMarker <= 1.0,
              "projected momentum and passive marker transport remain within their donor-cell bounds");
        check(first.diagnostics()
                      .projection.maximumDivergenceReductionRatio
                      < 1.0e-7,
              "visible transport retains the accepted two-level pressure projection");

        simwing::fsi::amr::ExternalFlowTransportCase staticWing(
            staticPlateScene());
        const auto staticWingFrame = staticWing.advance();
        const auto& staticWingDiagnostics =
            staticWing.diagnostics().momentum.staticWing;
        std::printf(
            "AMR static wing: triangles=%zu cut-cells=%zu normal=%.17g->%.17g->%.17g m/s divergence=%.17g 1/s slice=%d\n",
            staticWingDiagnostics.binding.triangleCount,
            staticWingDiagnostics.binding.activeCompositeCutCellCount,
            staticWingDiagnostics
                .maximumSurfaceNormalSpeedBeforeMetersPerSecond,
            staticWingDiagnostics
                .maximumSurfaceNormalSpeedAfterForcingMetersPerSecond,
            staticWingDiagnostics
                .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond,
            staticWing.diagnostics().projection
                .projectedMaximumDivergencePerSecond,
            simwing::viewer::describeCfdGrid(staticWingFrame).has_value()
                ? 1 : 0);
        check(staticWingDiagnostics.active
                  && staticWingDiagnostics.accepted
                  && staticWingDiagnostics.binding.accepted
                  && staticWingDiagnostics.binding.triangleCount == 2
                  && staticWingDiagnostics.binding
                         .activeCompositeCutCellCount > 0,
              "authoritative static scene triangles bind to positive-area cells on the AMR hierarchy");
        check(staticWingDiagnostics
                      .maximumSurfaceNormalSpeedBeforeMetersPerSecond > 0.0
                  && staticWingDiagnostics
                         .maximumSurfaceNormalSpeedAfterForcingMetersPerSecond
                      < staticWingDiagnostics
                            .maximumSurfaceNormalSpeedBeforeMetersPerSecond
                  && std::isfinite(
                      staticWingDiagnostics
                          .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond)
                  && staticWingDiagnostics.forcingProjectionIterations == 12
                  && staticWingDiagnostics
                         .maximumSurfaceNormalSpeedAfterProjectionMetersPerSecond
                      < 0.05 * staticWingDiagnostics
                                   .maximumSurfaceNormalSpeedBeforeMetersPerSecond,
              "alternating static-wing forcing and projection strongly reduce incident normal velocity");
        check(staticWingFrame.sceneChecksum == "sha256:amr-static-plate"
                  && staticWingFrame.solverCommit
                      == simwing::fsi::amr::staticWingExternalFlowSolverId
                  && staticWingFrame.vertices.size() == 24'580
                  && staticWingFrame.triangles.size() == 2
                  && field(staticWingFrame, "static wing cut cell") != nullptr,
              "static-wing frames retain the scene surface beside the CFD slice");
        check(simwing::viewer::describeCfdGrid(staticWingFrame).has_value()
                  && staticWing.diagnostics().projection
                         .projectedMaximumDivergencePerSecond < 1.0e-9,
              "surface geometry does not break the coarse CFD slice contract or composite projection");

        const auto slabSettings =
            simwing::fsi::amr::makeThreeSliceSpanwiseSlabSettings();
        simwing::fsi::amr::ExternalFlowTransportCase slab(
            staticPlateScene(), slabSettings);
        const auto slabFrame = slab.advance();
        const auto slabGrid = simwing::viewer::describeCfdGrid(slabFrame);
        const auto& slabInterface =
            slab.diagnostics().momentum.staticWing.binding;
        check(slabGrid
                  && slabGrid->cellCounts
                      == std::array<std::size_t, 3>{3, 48, 16}
                  && slabFrame.vertices.size() == 2'308,
              "three-slice spanwise slab publishes only its compact CFD volume and retained local surface");
        check(slabInterface.accepted
                  && slabInterface.clippedToWindTunnel
                  && !slabInterface.entirelyInsideWindTunnel
                  && slabInterface.triangleCount == 2
                  && std::abs(slabInterface.surfaceAreaSquareMeters - 2.25)
                      < 1.0e-12
                  && slabInterface.activeCompositeCutCellCount > 0,
              "spanwise slab clips authoritative material area without periodic images or invented closing faces");
        check(slabFrame.solverCommit
                      == simwing::fsi::amr::staticWingSlabExternalFlowSolverId
                  && slab.diagnostics().accepted
                  && slab.diagnostics().projection
                         .projectedMaximumDivergencePerSecond < 1.0e-9,
              "three-slice slab retains separate provenance and an accepted projected advance");
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr, "%d AMR external-flow case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all AMR external-flow case checks passed");
    return 0;
}
