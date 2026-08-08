#include "ram_air_cell_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "ram-air cell frame serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

const viewer::VectorField* vectorField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.vectorFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.vectorFields.end() ? nullptr : &*found;
}

double length(const fsi::StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

void testTopologyAndOpenMouthClamp() {
    fsi::RamAirCellCase simulation;
    const auto& definition = simulation.structure().definition();
    constexpr std::size_t nodeCount = 89;
    constexpr std::size_t triangleCount =
        fsi::ramAirCellPanelCount
        * 2 * fsi::ramAirCellTilesPerEdge
        * fsi::ramAirCellTilesPerEdge;
    constexpr std::size_t dihedralCount =
        fsi::ramAirCellPanelCount * 40;
    check(definition.nodes.size() == nodeCount
              && definition.triangles.size() == triangleCount
              && definition.membranes.size() == triangleCount
              && definition.dihedrals.size() == dihedralCount
              && definition.constraints.empty(),
          "ram-air cell contains five shared-node flexible panels");

    std::size_t fixedCount = 0;
    bool clampOnlyAtMouth = true;
    for (const auto& node : definition.nodes) {
        if (!node.fixed) {
            continue;
        }
        ++fixedCount;
        const bool firstTwoRows = node.positionMeters.x <= 1.25;
        const bool perimeter = node.positionMeters.y == -0.5
            || node.positionMeters.y == 0.5
            || node.positionMeters.z == 1.0
            || node.positionMeters.z == 2.0;
        clampOnlyAtMouth = clampOnlyAtMouth
            && firstTwoRows && perimeter;
    }
    check(fixedCount == 32 && clampOnlyAtMouth,
          "two perimeter rows clamp the open mouth's position and slope");
}

void testDeterministicFivePanelDeformation() {
    fsi::RamAirCellCase first;
    fsi::RamAirCellCase second;
    const auto initial = first.structure().checkpoint();
    viewer::DiagnosticFrame frame;
    double peakDisplacement = 0.0;
    double peakOutwardInflation = 0.0;
    double peakReaction = 0.0;
    double peakDirectReactionDifference = 0.0;
    double maximumOpeningMeanVelocity = 0.0;
    double maximumOpeningRmsVelocity = 0.0;
    for (std::size_t step = 0; step < 240; ++step) {
        frame = first.advance();
        check(serialized(frame) == serialized(second.advance()),
              "ram-air cell accepted frames are deterministic");
        peakDisplacement = std::max(
            peakDisplacement, first.maximumDisplacementMeters());
        peakOutwardInflation = std::max(
            peakOutwardInflation, first.maximumOutwardInflationMeters());

        const auto& diagnostics = first.diagnostics();
        peakReaction = std::max(
            peakReaction,
            length(diagnostics.fluidReactionForceNewtons));
        peakDirectReactionDifference = std::max(
            peakDirectReactionDifference,
            length({
                diagnostics.fluidReactionForceNewtons.x
                    - diagnostics.fluidPressureForceNewtons.x,
                diagnostics.fluidReactionForceNewtons.y
                    - diagnostics.fluidPressureForceNewtons.y,
                diagnostics.fluidReactionForceNewtons.z
                    - diagnostics.fluidPressureForceNewtons.z,
            }));
        check(first.fluidDiagnostics().projection.converged
                  && first.fluidDiagnostics().finite
                  && diagnostics.finite
                  && diagnostics.panels.size()
                         == fsi::ramAirCellPanelCount
                  && diagnostics.fluidDivergenceL2PerSecond < 1.0e-8
                  && diagnostics.maximumPanelForceResidualNewtons < 1.0e-8
                  && diagnostics.maximumPanelMomentResidualNewtonMeters
                         < 1.0e-8
                  && length(diagnostics.forceResidualNewtons) < 1.0e-8
                  && length(diagnostics.momentResidualNewtonMeters) < 1.0e-8,
              "ram-air cell closes every projection and conservative panel transfer");
        for (const auto& panel : diagnostics.panels) {
            check(panel.loadKind
                      == fsi::PlanarFaceResolvedLoadKind::
                             CompleteConstraintReaction,
                  "every ram-air cell panel receives the complete CFD reaction");
        }

        const auto* opening = scalarField(
            frame, "cell.opening_mean_velocity");
        if (opening != nullptr && opening->values.size() == 1) {
            maximumOpeningMeanVelocity = std::max(
                maximumOpeningMeanVelocity,
                std::abs(opening->values.front()));
        }
        const auto* openingRms = scalarField(
            frame, "cell.opening_rms_velocity");
        if (openingRms != nullptr && openingRms->values.size() == 1) {
            maximumOpeningRmsVelocity = std::max(
                maximumOpeningRmsVelocity,
                openingRms->values.front());
        }
    }

    const auto final = first.structure().checkpoint();
    bool anchorsUnchanged = true;
    for (std::size_t index = 0; index < initial.nodes.size(); ++index) {
        if (!first.structure().definition().nodes[index].fixed) {
            continue;
        }
        anchorsUnchanged = anchorsUnchanged
            && final.nodes[index].positionMeters
                   == initial.nodes[index].positionMeters
            && final.nodes[index].velocityMetersPerSecond
                   == fsi::StructureVector3{};
    }
    check(anchorsUnchanged,
          "ram-air flow preserves every mouth clamp node exactly");
    check(peakDisplacement > 1.0e-4
              && peakOutwardInflation > 1.0e-5
              && first.structure().diagnostics().finite,
          "complete CFD wall reaction visibly deforms the flexible cell");
    check(peakReaction > 0.01
              && peakDirectReactionDifference > 0.01
              && maximumOpeningMeanVelocity < 1.0e-8
              && maximumOpeningRmsVelocity > 0.01,
          "oscillating cross-flow conserves net mouth flux while exercising local exchange and direct reaction");
    check(frame.step == 240
              && frame.vertices.size() == initial.nodes.size()
              && frame.triangles.size()
                     == first.structure().definition().triangles.size()
              && frame.lines.empty()
              && frame.sceneChecksum == fsi::ramAirCellCaseChecksum
              && frame.solverCommit == fsi::ramAirCellCaseSolverId,
          "ram-air frame publishes the complete accepted five-panel shell");

    const auto* displacement = scalarField(frame, "cell.displacement");
    const auto* panel = scalarField(frame, "cell.panel");
    const auto* pressure = scalarField(
        frame, "cell.cfd_pressure_traction");
    const auto* reaction = scalarField(
        frame, "cell.cfd_complete_reaction_traction");
    const auto* gust = scalarField(frame, "cell.gust_speed");
    const auto* openingRms = scalarField(
        frame, "cell.opening_rms_velocity");
    const auto* divergence = scalarField(
        frame, "cell.fluid_divergence_l2");
    const auto* nodalForce = vectorField(frame, "cell.cfd_nodal_force");
    const auto* totalReaction = vectorField(
        frame, "cell.total_cfd_reaction");
    check(displacement != nullptr
              && displacement->association
                     == viewer::FieldAssociation::Vertex
              && displacement->values.size() == frame.vertices.size()
              && panel != nullptr
              && panel->association
                     == viewer::FieldAssociation::Triangle
              && panel->values.size() == frame.triangles.size()
              && pressure != nullptr
              && pressure->values.size() == frame.triangles.size()
              && reaction != nullptr
              && reaction->values.size() == frame.triangles.size()
              && gust != nullptr && gust->values.size() == 1
              && openingRms != nullptr
              && openingRms->values.size() == 1
              && divergence != nullptr
              && divergence->association
                     == viewer::FieldAssociation::Global
              && divergence->values.size() == 1
              && nodalForce != nullptr
              && nodalForce->association
                     == viewer::FieldAssociation::Vertex
              && nodalForce->values.size() == frame.vertices.size()
              && totalReaction != nullptr
              && totalReaction->association
                     == viewer::FieldAssociation::Global
              && totalReaction->values.size() == 1,
          "ram-air frames expose panel identity, deformation, flow, and mapped loads");
}

} // namespace

int main() {
    try {
        testTopologyAndOpenMouthClamp();
        testDeterministicFivePanelDeformation();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing ram-air cell check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing ram-air cell checks passed");
    return 0;
}
