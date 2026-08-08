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
    fsi::ClampedHemisphereCase simulation;
    const fsi::StructureDefinition& definition =
        simulation.structure().definition();
    constexpr std::size_t expectedNodes = 1
        + fsi::clampedHemisphereLatitudeSegments
              * fsi::clampedHemisphereRadialSegments;
    constexpr std::size_t expectedTriangles =
        fsi::clampedHemisphereRadialSegments
        * (2 * fsi::clampedHemisphereLatitudeSegments - 1);
    constexpr std::size_t expectedInteriorEdges =
        expectedNodes + expectedTriangles - 1
        - fsi::clampedHemisphereRadialSegments;
    check(definition.nodes.size() == expectedNodes
              && definition.triangles.size() == expectedTriangles
              && definition.membranes.size() == expectedTriangles
              && definition.dihedrals.size() == expectedInteriorEdges
              && definition.constraints.size()
                     == fsi::clampedHemisphereRadialSegments,
          "hemisphere topology contains the complete membrane disk");

    std::size_t fixedNodes = 0;
    for (const fsi::StructureNodeDefinition& node : definition.nodes) {
        if (node.fixed) {
            ++fixedNodes;
            check(std::abs(node.positionMeters.z) < 1.0e-14,
                  "only equatorial nodes are fixed");
        }
    }
    check(fixedNodes == fsi::clampedHemisphereRadialSegments,
          "the complete equatorial ring is clamped");
}

void testDeterministicPressureMotion() {
    fsi::ClampedHemisphereCase first;
    fsi::ClampedHemisphereCase second;
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
                     == fsi::clampedHemisphereRadialSegments
              && frame.sceneChecksum
                     == fsi::clampedHemisphereCaseChecksum
              && frame.solverCommit
                     == fsi::clampedHemisphereCaseSolverId,
          "hemisphere frame publishes the complete accepted surface");
    check(std::abs(first.apexRadialDisplacementMeters()) > 1.0e-3
              && first.pressurePascals() > 0.0
              && first.structure().diagnostics().finite,
          "pulsed follower pressure produces finite visible dome motion");

    const std::size_t boundaryBegin = 1
        + (fsi::clampedHemisphereLatitudeSegments - 1)
              * fsi::clampedHemisphereRadialSegments;
    bool boundaryUnchanged = true;
    for (std::size_t node = boundaryBegin; node < final.nodes.size(); ++node) {
        boundaryUnchanged = boundaryUnchanged
            && final.nodes[node].positionMeters
                   == initial.nodes[node].positionMeters
            && final.nodes[node].velocityMetersPerSecond
                   == fsi::StructureVector3{};
    }
    check(boundaryUnchanged,
          "pressure motion preserves the clamped equatorial ring exactly");

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
