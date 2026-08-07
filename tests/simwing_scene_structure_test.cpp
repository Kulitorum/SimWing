#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(double actual,
               double expected,
               double tolerance,
               const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message,
                     actual,
                     expected);
        ++failures;
    }
}

bool contains(const SceneStructureAssembly& assembly,
              SceneStructureDiagnosticCode code) {
    return std::ranges::any_of(
        assembly.diagnostics,
        [code](const SceneStructureDiagnostic& diagnostic) {
            return diagnostic.code == code;
        });
}

Scene surfaceScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-structure-test";
    scene.metadata.exporterVersion = "scene-structure-test/1";
    scene.regions = {
        {2, RegionKind::Cell, "cell-1"},
        {1, RegionKind::Outside, "outside"},
    };
    scene.vertices = {
        {13, {0.0, 1.0, 0.0}},
        {11, {3.0, 0.0, 0.0}},
        {10, {0.0, 0.0, 0.0}},
        {12, {3.0, 1.0, 0.0}},
    };
    scene.fabricMaterials = {
        {100, "test fabric", 900.0, 650.0, 220.0, 0.015,
         0.05, 0.02, 0.0, 0.0},
    };
    scene.triangles = {
        {501, {10, 12, 13}, {{{0.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {500, {10, 11, 12}, {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    scene.lineMaterials = {
        {110, "test line", 0.001, 0.0008, 1000.0, 1.0},
    };
    scene.attachments = {
        {301, AttachmentKind::SurfaceVertex, 12, 0, {}},
        {300, AttachmentKind::SurfaceVertex, 10, 0, {}},
    };
    scene.suspensionLines = {
        {400, 300, 301, 110, 2.5, SuspensionLineRole::Suspension},
    };
    return scene;
}

void reverseCollections(Scene& scene) {
    std::ranges::reverse(scene.regions);
    std::ranges::reverse(scene.vertices);
    std::ranges::reverse(scene.fabricMaterials);
    std::ranges::reverse(scene.seamMaterials);
    std::ranges::reverse(scene.triangles);
    std::ranges::reverse(scene.seams);
    std::ranges::reverse(scene.lineMaterials);
    std::ranges::reverse(scene.suspensionJunctions);
    std::ranges::reverse(scene.attachments);
    std::ranges::reverse(scene.suspensionLines);
}

void testDeterministicMappingAndStableIds() {
    const Scene scene = surfaceScene();
    const SceneStructureAssembly first = assembleSceneStructure(scene);
    check(first.ok(), "canonical surface scene assembles");

    Scene reordered = scene;
    reverseCollections(reordered);
    const SceneStructureAssembly second = assembleSceneStructure(reordered);
    check(second.ok(), "reordered surface scene assembles");
    check(first.mappings.nodeVertexIds == second.mappings.nodeVertexIds
              && first.mappings.triangleIds == second.mappings.triangleIds
              && first.mappings.membraneTriangleIds
                     == second.mappings.membraneTriangleIds
              && first.mappings.dihedralTriangleIds
                     == second.mappings.dihedralTriangleIds
              && first.mappings.nodeSuspensionJunctionIds
                     == second.mappings.nodeSuspensionJunctionIds
              && first.mappings.constraintSuspensionLineIds
                     == second.mappings.constraintSuspensionLineIds,
          "mapping arrays are deterministic under collection reordering");
    check(first.mappings.nodeVertexIds
              == std::vector<StableId>({10, 11, 12, 13}),
          "node mapping uses canonical stable vertex IDs");
    check(first.mappings.triangleIds
              == std::vector<StableId>({500, 501}),
          "triangle mapping uses canonical stable triangle IDs");
    check(first.mappings.membraneTriangleIds == first.mappings.triangleIds,
          "each scene triangle has a corresponding membrane mapping");
    check(first.mappings.constraintSuspensionLineIds
              == std::vector<StableId>({400}),
          "cable mapping retains its stable suspension-line ID");
    check(first.mappings.nodeIndex(12) == std::optional<std::size_t>(2)
              && !first.mappings.nodeIndex(999).has_value()
              && first.mappings.triangleIndex(501)
                     == std::optional<std::size_t>(1)
              && first.mappings.membraneIndex(500)
                     == std::optional<std::size_t>(0)
              && first.mappings.constraintIndex(400)
                     == std::optional<std::size_t>(0),
          "stable-ID reverse lookups agree with definition indices");

    check(first.definition.nodes.size() == second.definition.nodes.size()
              && first.definition.triangles.size()
                     == second.definition.triangles.size()
              && first.definition.membranes.size()
                     == second.definition.membranes.size()
              && first.definition.dihedrals.size()
                     == second.definition.dihedrals.size()
              && first.definition.constraints.size()
                     == second.definition.constraints.size(),
          "reordered assembly has the same structural cardinalities");
    for (std::size_t i = 0; i < first.definition.nodes.size(); ++i) {
        check(first.definition.nodes[i].positionMeters
                  == second.definition.nodes[i].positionMeters
                  && first.definition.nodes[i].massKg
                         == second.definition.nodes[i].massKg,
              "reordered assembly has identical node definitions");
    }
}

void testConservativeFabricMassAndMembranes() {
    const SceneStructureAssembly assembly = assembleSceneStructure(
        surfaceScene());
    check(assembly.ok(), "mass test scene assembles");
    const double nodeMass = std::accumulate(
        assembly.definition.nodes.begin(),
        assembly.definition.nodes.end(),
        0.0,
        [](double total, const StructureNodeDefinition& node) {
            return total + node.massKg;
        });
    checkNear(assembly.totalFabricMassKg, 0.1, 1.0e-15,
              "two material square metres at 0.05 kg/m2 has expected mass");
    checkNear(nodeMass, assembly.totalFabricMassKg, 1.0e-15,
              "all triangle fabric mass is conservatively lumped to nodes");
    check(assembly.totalFabricMassKg != 0.15,
          "fabric mass uses undeformed material-chart area, not spatial area");
    checkNear(assembly.definition.nodes[*assembly.mappings.nodeIndex(10)].massKg,
              0.1 / 3.0, 1.0e-15,
              "shared diagonal vertex receives two triangle shares");
    checkNear(assembly.definition.nodes[*assembly.mappings.nodeIndex(11)].massKg,
              0.1 / 6.0, 1.0e-15,
              "corner vertex receives one triangle share");

    check(assembly.definition.triangles.size() == 2
              && assembly.definition.membranes.size() == 2,
          "one membrane is assigned to every structural triangle");
    const StructureMembraneDefinition& first =
        assembly.definition.membranes.front();
    check(first.triangle == 0
              && first.materialCoordinates[1].x == 2.0
              && first.materialCoordinates[2].y == 1.0,
          "membrane retains the scene material chart");
    check(first.material.warpStiffnessNewtonsPerMeter == 900.0
              && first.material.weftStiffnessNewtonsPerMeter == 650.0
              && first.material.shearStiffnessNewtonsPerMeter == 220.0
              && first.material.dampingSeconds == 0.02
              && first.material.couplingStiffnessNewtonsPerMeter == 0.0
              && first.material.compressionStiffnessRatio == 1.0
              && first.role == StructureMaterialRole::Bulk,
          "scene fabric maps to a compatible Bulk orthotropic membrane");

    check(assembly.definition.dihedrals.size() == 1
              && assembly.mappings.dihedralTriangleIds
                     == std::vector<std::array<StableId, 2>>({{500, 501}}),
          "flat two-triangle strip creates one stable manifold hinge");
    const StructureDihedralDefinition& hinge =
        assembly.definition.dihedrals.front();
    check(hinge.nodes[0] == *assembly.mappings.nodeIndex(12)
              && hinge.nodes[1] == *assembly.mappings.nodeIndex(10)
              && hinge.nodes[2] == *assembly.mappings.nodeIndex(11)
              && hinge.nodes[3] == *assembly.mappings.nodeIndex(13),
          "hinge node order follows the lower stable-ID triangle orientation");
    checkNear(hinge.restAngleRadians, 0.0, 1.0e-15,
              "flat analytic strip has zero rest angle");
    checkNear(hinge.complianceRadiansPerNewtonMeter, 80.0 / 3.0, 1.0e-12,
              "strip-width bending rule gives analytic angular compliance");

    const StructureConstraintDefinition& cable =
        assembly.definition.constraints.front();
    check(cable.kind == StructureConstraintKind::Cable
              && cable.firstNode == *assembly.mappings.nodeIndex(10)
              && cable.secondNode == *assembly.mappings.nodeIndex(12),
          "surface-to-surface suspension line maps to its stable endpoint nodes");
    checkNear(cable.complianceMetersPerNewton, 0.0025, 0.0,
              "line cable compliance is rest length divided by axial stiffness");

    Structure constructed(assembly.definition);
    check(constructed.definition().nodes.size() == 4
              && constructed.definition().membranes.size() == 2,
          "successful bridge output is accepted by the XPBD Structure boundary");
}

void testFoldedStripRestAngle() {
    Scene scene = surfaceScene();
    for (Vertex& vertex : scene.vertices) {
        if (vertex.id == 10) {
            vertex.positionMeters = {0.0, 0.0, 0.0};
        } else if (vertex.id == 11) {
            vertex.positionMeters = {1.0, 0.0, 0.0};
        } else if (vertex.id == 12) {
            vertex.positionMeters = {0.0, 1.0, 0.0};
        } else if (vertex.id == 13) {
            vertex.positionMeters = {0.0, 0.0, -1.0};
        }
    }
    for (Triangle& triangle : scene.triangles) {
        if (triangle.id == 500) {
            triangle.vertexIds = {10, 11, 12};
        } else {
            triangle.vertexIds = {11, 10, 13};
        }
        triangle.materialCoordinates = {
            Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}};
    }
    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    check(assembly.ok() && assembly.definition.dihedrals.size() == 1,
          "folded strip: consistently oriented right-angle strip assembles");
    const StructureDihedralDefinition& hinge =
        assembly.definition.dihedrals.front();
    checkNear(hinge.restAngleRadians, std::acos(-1.0) / 2.0, 1.0e-15,
              "folded strip: signed rest angle comes from scene geometry");
    checkNear(hinge.complianceRadiansPerNewtonMeter, 200.0 / 3.0,
              1.0e-12,
              "folded strip: unit material charts give analytic compliance");
}

void testWeldedDifferentSheetsDoNotCreateHinge() {
    Scene scene = surfaceScene();
    scene.triangles[1].sheetId = 901;

    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    check(assembly.ok(),
          "sheet identity: welded triangles from different sheets assemble");
    check(assembly.definition.dihedrals.empty()
              && assembly.mappings.dihedralTriangleIds.empty(),
          "sheet identity: a welded inter-sheet edge is not a fabric hinge");
}

void testUnsupportedPilotHarnessIsExplicitAndDeterministic() {
    Scene scene = surfaceScene();
    scene.pilots = {
        {200, "pilot", 90.0, {0.0, 0.0, -1.0}, {}, {},
         {8.0, 10.0, 6.0}},
    };
    scene.attachments[1] =
        {302, AttachmentKind::PilotHarness, 0, 200, {0.0, 0.0, 0.2}};
    scene.suspensionLines[0].startAttachmentId = 301;
    scene.suspensionLines[0].endAttachmentId = 302;

    const SceneStructureAssembly first = assembleSceneStructure(scene);
    check(!first.ok() && !first.assembled,
          "pilot-harness scene is rejected transactionally");
    check(first.definition.nodes.empty()
              && first.mappings.nodeVertexIds.empty(),
          "unsupported assembly exposes no partial structure or mapping");
    check(contains(first, SceneStructureDiagnosticCode::UnsupportedPilot)
              && contains(first,
                          SceneStructureDiagnosticCode::UnsupportedPilotHarness)
              && contains(first,
                          SceneStructureDiagnosticCode::UnsupportedSuspensionTopology),
          "pilot, harness, and full suspension limitations are explicit");

    Scene reordered = scene;
    reverseCollections(reordered);
    const SceneStructureAssembly second = assembleSceneStructure(reordered);
    check(first.diagnostics == second.diagnostics,
          "unsupported diagnostics are deterministic under reordering");
}

void testSuspensionJunctionSegmentGraph() {
    Scene scene = surfaceScene();
    scene.suspensionJunctions = {
        {50, {1.5, 0.5, -0.5}, 0.02, false},
    };
    scene.attachments.push_back(
        {302, AttachmentKind::SuspensionJunction, 0, 0, {}, 50});
    scene.suspensionLines[0].endAttachmentId = 302;
    scene.suspensionLines.push_back(
        {401, 302, 301, 110, 1.5, SuspensionLineRole::Riser});

    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    check(assembly.ok(),
          "junction graph: explicit surface-junction segments assemble");
    check(assembly.definition.nodes.size() == 5
              && assembly.mappings.nodeSuspensionJunctionIds
                     == std::vector<StableId>({50})
              && assembly.mappings.junctionNodeIndex(50)
                     == std::optional<std::size_t>(4),
          "junction graph: stable junction maps to its own XPBD node");
    const std::size_t junction =
        *assembly.mappings.junctionNodeIndex(50);
    checkNear(assembly.definition.nodes[junction].massKg, 0.02, 0.0,
              "junction graph: explicit dynamic junction mass is preserved");
    check(assembly.definition.constraints.size() == 2
              && assembly.mappings.constraintSuspensionLineIds
                     == std::vector<StableId>({400, 401}),
          "junction graph: line segments remain separate stable constraints");
    check(assembly.definition.constraints[0].firstNode
                  == *assembly.mappings.nodeIndex(10)
              && assembly.definition.constraints[0].secondNode == junction
              && assembly.definition.constraints[1].firstNode == junction
              && assembly.definition.constraints[1].secondNode
                     == *assembly.mappings.nodeIndex(12),
          "junction graph: connectivity is not collapsed through the junction");

    Scene reordered = scene;
    reverseCollections(reordered);
    const SceneStructureAssembly reorderedAssembly =
        assembleSceneStructure(reordered);
    check(reorderedAssembly.ok()
              && reorderedAssembly.mappings.nodeSuspensionJunctionIds
                     == assembly.mappings.nodeSuspensionJunctionIds
              && reorderedAssembly.mappings.constraintSuspensionLineIds
                     == assembly.mappings.constraintSuspensionLineIds,
          "junction graph: mapping is deterministic under scene reordering");

    Scene collision = scene;
    collision.suspensionJunctions.front().id = 10;
    collision.attachments.back().suspensionJunctionId = 10;
    const SceneStructureAssembly collisionAssembly =
        assembleSceneStructure(collision);
    check(!collisionAssembly.ok()
              && contains(collisionAssembly,
                          SceneStructureDiagnosticCode::InvalidScene)
              && collisionAssembly.definition.nodes.empty(),
          "junction graph: vertex/junction stable-ID collision is rejected transactionally");
}

void testSeamAndBendingTopologyRejections() {
    Scene seamScene = surfaceScene();
    seamScene.seamMaterials = {
        {120, "test seam", 0.001, 3500.0},
    };
    seamScene.seams = {
        {610, 120, {10, 11}, {13, 12}},
    };
    check(validateScene(seamScene).ok(),
          "seam bridge: paired seam topology is valid scene data");
    const SceneStructureAssembly seamAssembly =
        assembleSceneStructure(seamScene);
    check(!seamAssembly.ok()
              && contains(seamAssembly,
                          SceneStructureDiagnosticCode::UnsupportedSeam)
              && seamAssembly.definition.nodes.empty(),
          "seam bridge: unsupported stitch physics rejects transactionally");

    Scene inconsistent = surfaceScene();
    inconsistent.triangles[0].vertexIds = {12, 10, 13};
    inconsistent.triangles[0].materialCoordinates = {
        Vec2{0.0, 0.0}, Vec2{2.0, 0.0}, Vec2{0.0, 1.0}};
    check(validateScene(inconsistent).ok(),
          "bending topology: oppositely oriented face remains a valid scene");
    const SceneStructureAssembly inconsistentAssembly =
        assembleSceneStructure(inconsistent);
    check(!inconsistentAssembly.ok()
              && contains(
                  inconsistentAssembly,
                  SceneStructureDiagnosticCode::UnsupportedBendingTopology),
          "bending topology: same-direction shared edge is diagnosed");

    Scene nonManifold = surfaceScene();
    nonManifold.vertices.push_back({14, {1.5, 0.5, 1.0}});
    nonManifold.triangles.push_back(
        {502,
         {10, 12, 14},
         {{{0.0, 0.0}, {2.0, 0.0}, {1.0, 1.0}}},
         1,
         2,
         100,
         900,
         SurfaceRole::Skin});
    check(validateScene(nonManifold).ok(),
          "bending topology: three-face edge remains serializable scene data");
    const SceneStructureAssembly nonManifoldAssembly =
        assembleSceneStructure(nonManifold);
    check(!nonManifoldAssembly.ok()
              && contains(
                  nonManifoldAssembly,
                  SceneStructureDiagnosticCode::UnsupportedBendingTopology),
          "bending topology: non-manifold shared edge is diagnosed");
}

void testExplicitValidationAndAssemblyRejections() {
    Scene invalid = surfaceScene();
    invalid.triangles.front().materialId = 999;
    const SceneStructureAssembly invalidResult = assembleSceneStructure(invalid);
    check(!invalidResult.ok()
              && contains(invalidResult,
                          SceneStructureDiagnosticCode::InvalidScene),
          "invalid scene is rejected before mapping");
    check(std::ranges::any_of(
              invalidResult.diagnostics,
              [](const SceneStructureDiagnostic& diagnostic) {
                  return diagnostic.sceneValidationCode
                      == ValidationCode::MissingMaterialReference;
              }),
          "wrapped invalid-scene diagnostic retains its validation code");

    Scene isolated = surfaceScene();
    isolated.vertices.push_back({99, {4.0, 4.0, 4.0}});
    const SceneStructureAssembly isolatedResult =
        assembleSceneStructure(isolated);
    check(!isolatedResult.ok()
              && contains(isolatedResult,
                          SceneStructureDiagnosticCode::ZeroMassDynamicNode),
          "dynamic scene vertex without fabric mass is rejected");

    Scene badChart = surfaceScene();
    badChart.triangles.front().materialCoordinates = {
        Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{2.0, 0.0}};
    const SceneStructureAssembly chartResult =
        assembleSceneStructure(badChart);
    check(!chartResult.ok()
              && contains(chartResult,
                          SceneStructureDiagnosticCode::DegenerateMaterialCoordinates),
          "material chart incompatible with XPBD is rejected explicitly");

    Scene reversedChart = surfaceScene();
    std::swap(reversedChart.triangles.front().materialCoordinates[1],
              reversedChart.triangles.front().materialCoordinates[2]);
    const SceneStructureAssembly reversedChartResult =
        assembleSceneStructure(reversedChart);
    check(!reversedChartResult.ok()
              && contains(reversedChartResult,
                          SceneStructureDiagnosticCode::DegenerateMaterialCoordinates),
          "negative-area material chart rejected before Structure construction");

    Scene illConditioned = surfaceScene();
    illConditioned.fabricMaterials.front().warpStiffnessNewtonsPerMeter =
        1.0e100;
    const SceneStructureAssembly materialResult =
        assembleSceneStructure(illConditioned);
    check(!materialResult.ok()
              && contains(materialResult,
                          SceneStructureDiagnosticCode::MaterialIncompatible),
          "ill-conditioned scene material is rejected before Structure construction");

    SceneStructureLimits limits;
    limits.maximumNodes = 3;
    const SceneStructureAssembly bounded = assembleSceneStructure(
        surfaceScene(), limits);
    check(!bounded.ok()
              && contains(bounded,
                          SceneStructureDiagnosticCode::MappingOverflow),
          "configured node mapping bound rejects an oversized mapping");

    Scene sameEndpoint = surfaceScene();
    sameEndpoint.attachments[0].vertexId = 10;
    const SceneStructureAssembly topology = assembleSceneStructure(
        sameEndpoint);
    check(!topology.ok()
              && contains(topology,
                          SceneStructureDiagnosticCode::UnsupportedSuspensionTopology),
          "surface cable with identical endpoint nodes is rejected");
}

} // namespace

int main() {
    testDeterministicMappingAndStableIds();
    testConservativeFabricMassAndMembranes();
    testFoldedStripRestAngle();
    testWeldedDifferentSheetsDoNotCreateHinge();
    testUnsupportedPilotHarnessIsExplicitAndDeterministic();
    testSuspensionJunctionSegmentGraph();
    testSeamAndBendingTopologyRejections();
    testExplicitValidationAndAssemblyRejections();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing scene-structure check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing scene-structure checks passed");
    return 0;
}
