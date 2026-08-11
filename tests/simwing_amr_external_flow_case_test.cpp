#include "amr_external_flow_case.h"
#include "cfd_slice.h"

#include <algorithm>
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
