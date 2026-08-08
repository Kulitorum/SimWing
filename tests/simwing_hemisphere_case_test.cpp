#include "hemisphere_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(bool condition, const char* message) {
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
          "hemisphere frame serializes");
    return bytes;
}

void testTopologyAndBoundary() {
    fsi::AnchoredHemisphereCase simulation;
    const fsi::StructureDefinition& definition =
        simulation.structure().definition();
    constexpr std::size_t expectedNodes = 1
        + fsi::hemisphereLatitudeSegments
              * fsi::hemisphereRadialSegments;
    constexpr std::size_t expectedTriangles =
        fsi::hemisphereRadialSegments
        * (2 * fsi::hemisphereLatitudeSegments - 1);
    constexpr std::size_t expectedInteriorEdges =
        expectedNodes + expectedTriangles - 1
        - fsi::hemisphereRadialSegments;
    check(definition.nodes.size() == expectedNodes
              && definition.triangles.size() == expectedTriangles
              && definition.membranes.size() == expectedTriangles
              && definition.dihedrals.size() == expectedInteriorEdges
              && definition.constraints.size()
                     == fsi::hemisphereRadialSegments,
          "hemisphere topology contains the complete membrane disk");

    std::size_t fixedNodes = 0;
    for (std::size_t nodeIndex = 0;
         nodeIndex < definition.nodes.size();
         ++nodeIndex) {
        const fsi::StructureNodeDefinition& node = definition.nodes[nodeIndex];
        if (node.fixed) {
            ++fixedNodes;
            check(std::abs(node.positionMeters.z) < 1.0e-14,
                  "only equatorial nodes are fixed");
            const std::size_t boundaryBegin = 1
                + (fsi::hemisphereLatitudeSegments - 1)
                      * fsi::hemisphereRadialSegments;
            check(nodeIndex >= boundaryBegin
                      && fsi::isHemisphereAnchorLongitude(
                          nodeIndex - boundaryBegin),
                  "only the retained original rim points are fixed");
        }
    }
    check(fixedNodes == fsi::hemisphereAnchorCount,
          "exactly two original equatorial points are anchored");
    const bool compliantRim = std::all_of(
        definition.constraints.begin(), definition.constraints.end(),
        [](const fsi::StructureConstraintDefinition& constraint) {
            return constraint.kind == fsi::StructureConstraintKind::Distance
                && constraint.complianceMetersPerNewton
                    == fsi::hemisphereRimComplianceMetersPerNewton;
        });
    check(compliantRim,
          "every visible rim segment is compliant rather than rigid");
}

void testDeterministicPressureMotion() {
    fsi::AnchoredHemisphereCase first;
    fsi::AnchoredHemisphereCase second;
    const auto initial = first.structure().checkpoint();
    viewer::DiagnosticFrame frame;
    for (std::size_t step = 0; step < 180; ++step) {
        frame = first.advance();
        check(serialized(frame) == serialized(second.advance()),
              "hemisphere accepted frames are deterministic");
    }

    const auto final = first.structure().checkpoint();
    check(frame.step == 180
              && frame.vertices.size() == initial.nodes.size()
              && frame.triangles.size()
                     == first.structure().definition().triangles.size()
              && frame.lines.size()
                     == fsi::hemisphereRadialSegments
              && frame.sceneChecksum
                     == fsi::anchoredHemisphereCaseChecksum
              && frame.solverCommit
                     == fsi::anchoredHemisphereCaseSolverId,
          "hemisphere frame publishes the complete accepted surface");
    check(std::abs(first.apexRadialDisplacementMeters()) > 1.0e-3
              && first.pressurePascals() > 0.0
              && first.structure().diagnostics().finite,
          "pulsed follower pressure produces finite visible dome motion");

    const std::size_t boundaryBegin = 1
        + (fsi::hemisphereLatitudeSegments - 1)
              * fsi::hemisphereRadialSegments;
    bool anchorsUnchanged = true;
    double maximumFreeBoundaryMotion = 0.0;
    for (std::size_t node = boundaryBegin; node < final.nodes.size(); ++node) {
        const std::size_t longitude = node - boundaryBegin;
        const bool anchored =
            fsi::isHemisphereAnchorLongitude(longitude);
        if (anchored) {
            anchorsUnchanged = anchorsUnchanged
                && final.nodes[node].positionMeters
                       == initial.nodes[node].positionMeters
                && final.nodes[node].velocityMetersPerSecond
                       == fsi::StructureVector3{};
        } else {
            const fsi::StructureVector3 delta{
                final.nodes[node].positionMeters.x
                    - initial.nodes[node].positionMeters.x,
                final.nodes[node].positionMeters.y
                    - initial.nodes[node].positionMeters.y,
                final.nodes[node].positionMeters.z
                    - initial.nodes[node].positionMeters.z};
            maximumFreeBoundaryMotion = std::max(
                maximumFreeBoundaryMotion,
                std::sqrt(delta.x * delta.x + delta.y * delta.y
                          + delta.z * delta.z));
        }
    }
    check(anchorsUnchanged,
          "pressure motion preserves both anchors exactly");
    check(maximumFreeBoundaryMotion > 0.5,
          "pressure motion produces a large free-rim flap");
    check(first.maximumFreeRimDisplacementMeters()
              == maximumFreeBoundaryMotion,
          "reported free-rim motion matches the accepted structure state");

    const auto radialField = std::find_if(
        frame.scalarFields.begin(), frame.scalarFields.end(),
        [](const viewer::ScalarField& field) {
            return field.name == "dome.radial_displacement";
        });
    check(radialField != frame.scalarFields.end()
              && radialField->association == viewer::FieldAssociation::Vertex
              && radialField->values.size() == frame.vertices.size(),
          "hemisphere frames expose radial deformation for colouring");
}

} // namespace

int main() {
    testTopologyAndBoundary();
    testDeterministicPressureMotion();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing hemisphere check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing hemisphere checks passed");
    return 0;
}
