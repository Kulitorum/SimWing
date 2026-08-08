#include "projected_flag_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
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
          "projected flag frame serializes");
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

void testTopologyAndFixedReferenceContract() {
    fsi::ProjectedGustFlagCase simulation;
    const auto& definition = simulation.structure().definition();
    constexpr std::size_t nodeCount =
        fsi::projectedFlagNodesPerSide
        * fsi::projectedFlagNodesPerSide;
    constexpr std::size_t triangleCount =
        2 * fsi::projectedFlagTilesPerSide
        * fsi::projectedFlagTilesPerSide;
    constexpr std::size_t interiorEdgeCount = 40;
    check(definition.nodes.size() == nodeCount
              && definition.triangles.size() == triangleCount
              && definition.membranes.size() == triangleCount
              && definition.dihedrals.size() == interiorEdgeCount
              && definition.constraints.empty(),
          "projected flag contains the complete flexible fabric grid");

    std::size_t fixedCount = 0;
    for (std::size_t index = 0; index < definition.nodes.size(); ++index) {
        const auto& node = definition.nodes[index];
        if (node.fixed) {
            ++fixedCount;
            check(index % fsi::projectedFlagNodesPerSide <= 1,
                  "only the authored leading clamp strip is fixed");
        }
    }
    check(fixedCount == 2 * fsi::projectedFlagNodesPerSide,
          "two node rows encode clamp position and slope");
}

void testDeterministicProjectedGustMotion() {
    fsi::ProjectedGustFlagCase first;
    fsi::ProjectedGustFlagCase second;
    const auto initial = first.structure().checkpoint();
    viewer::DiagnosticFrame frame;
    double peakNormalMotion = 0.0;
    double peakReactionDifference = 0.0;
    double minimumReaction = std::numeric_limits<double>::infinity();
    double maximumReaction = -std::numeric_limits<double>::infinity();
    for (std::size_t step = 0; step < 240; ++step) {
        frame = first.advance();
        check(serialized(frame) == serialized(second.advance()),
              "projected flag accepted frames are deterministic");
        peakNormalMotion = std::max(
            peakNormalMotion,
            first.maximumNormalDisplacementMeters());
        const auto& fluid = first.fluidDiagnostics();
        const auto surface = std::ranges::find_if(
            fluid.surfaces,
            [](const auto& candidate) {
                return candidate.stableId
                    == fsi::projectedFlagSurfaceStableId;
            });
        if (surface != fluid.surfaces.end()) {
            peakReactionDifference = std::max(
                peakReactionDifference,
                std::abs(surface->constraintReactionForceNewtons.x
                         - surface->pressureForceNewtons.x));
            minimumReaction = std::min(
                minimumReaction,
                surface->constraintReactionForceNewtons.x);
            maximumReaction = std::max(
                maximumReaction,
                surface->constraintReactionForceNewtons.x);
        }
        const auto& transfer = first.transferDiagnostics();
        check(fluid.projection.converged && fluid.finite
                  && fluid.projection.divergenceL2AfterPerSecond < 1.0e-8
                  && transfer.loadKind
                      == fsi::PlanarFaceResolvedLoadKind::
                             CompleteConstraintReaction
                  && transfer.forceResidualNormNewtons < 1.0e-8
                  && transfer.momentResidualNormNewtonMeters < 1.0e-8,
              "projected flag closes each CFD projection and conservative load transfer");
    }

    const auto final = first.structure().checkpoint();
    check(frame.step == 240
              && frame.vertices.size() == initial.nodes.size()
              && frame.triangles.size()
                     == first.structure().definition().triangles.size()
              && frame.lines.empty()
              && frame.sceneChecksum
                     == fsi::projectedGustFlagCaseChecksum
              && frame.solverCommit
                     == fsi::projectedGustFlagCaseSolverId,
          "projected flag frame publishes the complete accepted fabric surface");

    bool anchorsUnchanged = true;
    for (std::size_t z = 0;
         z < fsi::projectedFlagNodesPerSide; ++z) {
        for (std::size_t y = 0; y < 2; ++y) {
            const std::size_t node =
                z * fsi::projectedFlagNodesPerSide + y;
            anchorsUnchanged = anchorsUnchanged
                && final.nodes[node].positionMeters
                    == initial.nodes[node].positionMeters
                && final.nodes[node].velocityMetersPerSecond
                    == fsi::StructureVector3{};
        }
    }
    check(anchorsUnchanged,
          "projected gust preserves the clamped edge exactly");
    check(peakNormalMotion > 0.02
              && first.maximumFreeEdgeDisplacementMeters() > 0.01
              && first.structure().diagnostics().finite,
          "CFD reaction produces visible non-rigid free-edge fabric motion");
    check(peakReactionDifference > 0.1
              && minimumReaction < -0.1
              && maximumReaction > 0.1,
          "accelerating cross-flow exercises pressure plus signed direct constraint impulse");

    const auto* displacement = scalarField(
        frame, "flag.normal_displacement");
    const auto* pressure = scalarField(
        frame, "flag.cfd_pressure_traction");
    const auto* reaction = scalarField(
        frame, "flag.cfd_complete_reaction_traction");
    const auto* divergence = scalarField(
        frame, "flag.fluid_divergence_l2");
    const auto* nodalForce = vectorField(
        frame, "flag.cfd_nodal_force");
    check(displacement != nullptr
              && displacement->association
                     == viewer::FieldAssociation::Vertex
              && displacement->values.size() == frame.vertices.size()
              && pressure != nullptr
              && pressure->association
                     == viewer::FieldAssociation::Triangle
              && pressure->values.size() == frame.triangles.size()
              && reaction != nullptr
              && reaction->values.size() == frame.triangles.size()
              && divergence != nullptr
              && divergence->association
                     == viewer::FieldAssociation::Global
              && divergence->values.size() == 1
              && nodalForce != nullptr
              && nodalForce->association
                     == viewer::FieldAssociation::Vertex
              && nodalForce->values.size() == frame.vertices.size(),
          "projected flag frames expose deformation, CFD loads, divergence, and mapped nodal force");
}

} // namespace

int main() {
    try {
        testTopologyAndFixedReferenceContract();
        testDeterministicProjectedGustMotion();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing projected flag check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing projected flag checks passed");
    return 0;
}
