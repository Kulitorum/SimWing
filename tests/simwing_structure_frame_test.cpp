#include "structure_frame.h"

#include "scene_structure.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using simwing::fsi::Structure;
using simwing::fsi::StructureConstraintKind;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureNodeDefinition;
using simwing::fsi::StructureStepSettings;
using namespace simwing::viewer;

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

const ScalarField* scalarField(
    const DiagnosticFrame& frame,
    std::string_view name) {
    for (const ScalarField& field : frame.scalarFields) {
        if (field.name == name) {
            return &field;
        }
    }
    return nullptr;
}

const VectorField* vectorField(
    const DiagnosticFrame& frame,
    std::string_view name) {
    for (const VectorField& field : frame.vectorFields) {
        if (field.name == name) {
            return &field;
        }
    }
    return nullptr;
}

StructureDefinition sampleDefinition() {
    StructureDefinition definition;
    definition.nodes = {
        StructureNodeDefinition{{0.0, 0.0, 0.0}, 0.0, true},
        StructureNodeDefinition{{1.0, 0.0, 0.0}, 2.0, false},
        StructureNodeDefinition{{0.0, 1.0, 0.0}, 1.0, false},
    };
    definition.triangles.push_back({{0, 1, 2}});
    definition.constraints.push_back(
        {StructureConstraintKind::Distance, 0, 1, 1.0, 0.0});
    definition.constraints.push_back(
        {StructureConstraintKind::Cable, 1, 2, 0.5, 0.0});

    simwing::fsi::StructureMembraneDefinition membrane;
    membrane.triangle = 0;
    membrane.materialCoordinates = {
        simwing::fsi::StructureVector2{0.0, 0.0},
        simwing::fsi::StructureVector2{1.0, 0.0},
        simwing::fsi::StructureVector2{0.0, 1.0},
    };
    membrane.material.warpStiffnessNewtonsPerMeter = 800.0;
    membrane.material.weftStiffnessNewtonsPerMeter = 500.0;
    membrane.material.couplingStiffnessNewtonsPerMeter = 100.0;
    membrane.material.shearStiffnessNewtonsPerMeter = 180.0;
    definition.membranes.push_back(membrane);
    return definition;
}

StructureFrameMappingDefinition sampleMappingDefinition() {
    StructureFrameMappingDefinition mapping;
    mapping.vertexStableIds = {101, 102, 103};
    mapping.triangles = {{201, 11, 12}};
    mapping.lines = {{301, 7}, {302, 8}};
    return mapping;
}

StructureFrameContext sampleContext() {
    StructureFrameContext context;
    context.sceneChecksum = "sha256:structure-frame-fixture";
    context.solverCommit = "abcdef123456";
    context.timeStepSeconds = 0.25;
    context.couplingIteration = 4;
    context.couplingResiduals = {1.0e-6, 2.0e-5, 3.0e-7, 4.0e-8,
                                 -5.0e-4};
    context.conservation.fluidMassKilograms = 3.5;
    context.conservation.totalEnergyJoules = 9.0;
    return context;
}

simwing::fsi::Scene sampleAssemblyScene() {
    using namespace simwing::fsi;
    Scene scene;
    scene.metadata.designChecksum = "sha256:structure-frame-scene";
    scene.metadata.exporterVersion = "structure-frame-test/1";
    scene.regions = {
        {91, RegionKind::Cell, "cell"},
        {90, RegionKind::Outside, "outside"},
    };
    scene.vertices = {
        {13, {0.0, 1.0, 0.0}},
        {11, {2.0, 0.0, 0.0}},
        {10, {0.0, 0.0, 0.0}},
        {12, {2.0, 1.0, 0.0}},
    };
    scene.fabricMaterials = {
        {100, "viewer fabric", 900.0, 650.0, 220.0, 0.015,
         0.05, 0.02, 0.0, 0.0},
    };
    scene.triangles = {
        {501, {10, 12, 13}, {{{0.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}}},
         91, 90, 100, 901, SurfaceRole::Rib},
        {500, {10, 11, 12}, {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}}},
         90, 91, 100, 900, SurfaceRole::Skin},
    };
    scene.lineMaterials = {
        {110, "viewer line", 0.001, 0.0008, 1000.0, 1.0},
    };
    scene.suspensionJunctions = {
        {50, {1.0, 0.5, -0.5}, 0.02, false},
    };
    scene.attachments = {
        {314, AttachmentKind::SuspensionJunction, 0, 0, {}, 50},
        {313, AttachmentKind::SurfaceVertex, 13, 0, {}},
        {311, AttachmentKind::SurfaceVertex, 12, 0, {}},
        {310, AttachmentKind::SurfaceVertex, 10, 0, {}},
        {312, AttachmentKind::SurfaceVertex, 11, 0, {}},
    };
    scene.suspensionLines = {
        {702, 314, 313, 110, 1.0, SuspensionLineRole::Riser},
        {701, 312, 314, 110, 1.0, SuspensionLineRole::Riser},
        {700, 310, 311, 110, 2.5, SuspensionLineRole::Brake},
    };
    return scene;
}

simwing::fsi::Scene samplePilotAssemblyScene() {
    using namespace simwing::fsi;
    Scene scene = sampleAssemblyScene();
    scene.pilots = {
        {800, "viewer pilot", 90.0, {1.0, 0.5, -2.0}, {}, {},
         {8.0, 10.0, 6.0}},
    };
    scene.attachments.push_back(
        {315, AttachmentKind::PilotHarness, 0, 800,
         {0.0, 0.0, 0.5}});
    scene.suspensionLines = {
        {700, 310, 314, 110, 1.25, SuspensionLineRole::Suspension},
        {701, 312, 314, 110, 1.25, SuspensionLineRole::Suspension},
        {702, 314, 315, 110, 1.0, SuspensionLineRole::Riser},
    };
    return scene;
}

void reverseSceneCollections(simwing::fsi::Scene& scene) {
    std::ranges::reverse(scene.regions);
    std::ranges::reverse(scene.vertices);
    std::ranges::reverse(scene.fabricMaterials);
    std::ranges::reverse(scene.triangles);
    std::ranges::reverse(scene.lineMaterials);
    std::ranges::reverse(scene.suspensionJunctions);
    std::ranges::reverse(scene.attachments);
    std::ranges::reverse(scene.suspensionLines);
}

void testExactInitialFrame() {
    Structure structure(sampleDefinition());
    structure.addExternalForce(1, {4.0, -2.0, 1.0});
    const StructureFrameMapping mapping(
        structure, sampleMappingDefinition());
    const StructureFrameContext context = sampleContext();
    const DiagnosticFrame frame =
        buildStructureFrame(structure, mapping, context);

    check(frame.sceneChecksum == context.sceneChecksum
              && frame.solverCommit == context.solverCommit
              && frame.step == 0
              && frame.simulationTimeSeconds == 0.0
              && frame.timeStepSeconds == 0.25
              && frame.couplingIteration == 4,
          "content: provenance and committed time are exact");
    check(frame.vertices.size() == 3
              && frame.vertices[0].stableId == 101
              && frame.vertices[1].stableId == 102
              && frame.vertices[2].stableId == 103
              && frame.vertices[1].positionMetres.x == 1.0,
          "content: positions retain structure order and stable IDs");
    check(frame.triangles.size() == 1
              && frame.triangles[0].stableId == 201
              && frame.triangles[0].vertex0 == 0
              && frame.triangles[0].vertex1 == 1
              && frame.triangles[0].vertex2 == 2
              && frame.triangles[0].negativeRegionId == 11
              && frame.triangles[0].positiveRegionId == 12,
          "content: triangle topology and two-sided regions are exact");
    check(frame.lines.size() == 2
              && frame.lines[0].stableId == 301
              && frame.lines[0].vertex0 == 0
              && frame.lines[0].vertex1 == 1
              && frame.lines[0].role == 7
              && frame.lines[1].stableId == 302
              && frame.lines[1].vertex0 == 1
              && frame.lines[1].vertex1 == 2
              && frame.lines[1].role == 8,
          "content: constraints become stable diagnostic lines");

    const VectorField* velocity =
        vectorField(frame, "structure.velocity");
    const VectorField* pending =
        vectorField(frame, "structure.pending_external_force");
    const VectorField* totalPending =
        vectorField(frame, "structure.total_pending_external_force");
    check(velocity != nullptr
              && velocity->association == FieldAssociation::Vertex
              && velocity->values.size() == 3
              && velocity->values[1].x == 0.0,
          "content: per-node velocity is published");
    check(pending != nullptr
              && pending->association == FieldAssociation::Vertex
              && pending->values.size() == 3
              && pending->values[0].x == 0.0
              && pending->values[1].x == 4.0
              && pending->values[1].y == -2.0
              && pending->values[1].z == 1.0,
          "content: exact per-node pending force is published");
    check(totalPending != nullptr
              && totalPending->association == FieldAssociation::Global
              && totalPending->values.size() == 1
              && totalPending->values[0].x == 4.0
              && totalPending->values[0].y == -2.0
              && totalPending->values[0].z == 1.0,
          "content: aggregate pending force remains available");

    const ScalarField* lengths =
        scalarField(frame, "structure.constraint_length");
    const ScalarField* violations =
        scalarField(frame, "structure.constraint_violation");
    check(lengths != nullptr && lengths->values.size() == 2,
          "content: each line receives a geometric length");
    checkNear(lengths->values[0], 1.0, 0.0,
              "content: bilateral line length is exact");
    checkNear(lengths->values[1], std::sqrt(2.0), 1.0e-15,
              "content: cable line length is exact");
    check(violations != nullptr && violations->values.size() == 2,
          "content: each line receives a constraint violation");
    checkNear(violations->values[0], 0.0, 0.0,
              "content: satisfied bilateral line has no violation");
    checkNear(violations->values[1], std::sqrt(2.0) - 0.5, 1.0e-15,
              "content: cable violation is extension only");

    const ScalarField* globalStrain = scalarField(
        frame, "structure.maximum_absolute_membrane_strain");
    check(globalStrain != nullptr
              && globalStrain->association == FieldAssociation::Global
              && globalStrain->values == std::vector<double>{0.0},
          "content: only the public global membrane strain is emitted");
    bool hasTriangleScalar = false;
    for (const ScalarField& field : frame.scalarFields) {
        hasTriangleScalar = hasTriangleScalar
            || field.association == FieldAssociation::Triangle;
    }
    check(!hasTriangleScalar,
          "content: unavailable per-triangle strain is not invented");

    ProtocolError error;
    check(validateFrame(frame, &error),
          "content: constructed frame satisfies the viewer protocol");
}

void testSerializationCompatibility() {
    Structure structure(sampleDefinition());
    const StructureFrameMapping mapping(
        structure, sampleMappingDefinition());
    const DiagnosticFrame frame =
        buildStructureFrame(structure, mapping, sampleContext());
    ProtocolError error;
    std::vector<std::uint8_t> bytes;
    check(serializeFrame(frame, bytes, &error),
          "serialization: adapter frame serializes with protocol v1");
    DiagnosticFrame decoded;
    check(deserializeFrame(bytes, decoded, &error),
          "serialization: adapter frame deserializes with protocol v1");
    check(decoded.sceneChecksum == frame.sceneChecksum
              && decoded.solverCommit == frame.solverCommit
              && decoded.vertices.size() == frame.vertices.size()
              && decoded.vertices[2].stableId == 103
              && decoded.triangles[0].negativeRegionId == 11
              && decoded.triangles[0].positiveRegionId == 12
              && decoded.scalarFields.size() == frame.scalarFields.size()
              && decoded.vectorFields.size() == frame.vectorFields.size(),
          "serialization: topology, regions, provenance, and fields survive");
    std::vector<std::uint8_t> reencoded;
    check(serializeFrame(decoded, reencoded, &error) && reencoded == bytes,
          "serialization: adapter frame round trip is byte deterministic");
}

void testEvolvingStructureFrames() {
    Structure structure(sampleDefinition());
    const StructureFrameMapping mapping(
        structure, sampleMappingDefinition());
    StructureFrameContext context = sampleContext();
    const DiagnosticFrame before =
        buildStructureFrame(structure, mapping, context);

    structure.addExternalForce(1, {4.0, 0.0, 0.0});
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.25;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    static_cast<void>(structure.step(settings));
    const DiagnosticFrame after =
        buildStructureFrame(structure, mapping, context);

    check(before.step == 0 && after.step == 1
              && after.simulationTimeSeconds == 0.25,
          "evolution: frame time follows accepted structure steps");
    checkNear(after.vertices[1].positionMetres.x, 1.125, 1.0e-15,
              "evolution: frame geometry follows the committed node state");
    const VectorField* velocity =
        vectorField(after, "structure.velocity");
    const VectorField* pending =
        vectorField(after, "structure.pending_external_force");
    const VectorField* applied =
        vectorField(after, "structure.last_applied_external_force");
    checkNear(velocity->values[1].x, 0.5, 1.0e-15,
              "evolution: node velocity follows the accepted step");
    check(pending->values[1].x == 0.0,
          "evolution: accepted loads disappear from pending-node field");
    check(applied->values[0].x == 4.0,
          "evolution: aggregate accepted load remains visible");
}

template <class Callable>
void checkInvalid(Callable&& callable, const char* message) {
    bool rejected = false;
    try {
        callable();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void testInvalidMappings() {
    Structure structure(sampleDefinition());

    StructureFrameMappingDefinition invalid = sampleMappingDefinition();
    invalid.vertexStableIds.pop_back();
    checkInvalid(
        [&] { StructureFrameMapping mapping(structure, invalid); },
        "mapping: incomplete vertex mapping is rejected");

    invalid = sampleMappingDefinition();
    invalid.vertexStableIds[1] = 0;
    checkInvalid(
        [&] { StructureFrameMapping mapping(structure, invalid); },
        "mapping: zero vertex stable ID is rejected");

    invalid = sampleMappingDefinition();
    invalid.vertexStableIds[1] = invalid.vertexStableIds[0];
    checkInvalid(
        [&] { StructureFrameMapping mapping(structure, invalid); },
        "mapping: duplicate vertex stable ID is rejected");

    invalid = sampleMappingDefinition();
    invalid.triangles[0].positiveRegionId =
        invalid.triangles[0].negativeRegionId;
    const StructureFrameMapping sameRegionMapping(structure, invalid);
    const DiagnosticFrame sameRegionFrame = buildStructureFrame(
        structure, sameRegionMapping, sampleContext());
    check(sameRegionFrame.triangles[0].negativeRegionId
              == sameRegionFrame.triangles[0].positiveRegionId,
          "mapping: internal sheets preserve one connected region on both sides");

    invalid = sampleMappingDefinition();
    invalid.lines[1].stableId = 0;
    checkInvalid(
        [&] { StructureFrameMapping mapping(structure, invalid); },
        "mapping: zero line stable ID is rejected");

    const StructureFrameMapping mapping(
        structure, sampleMappingDefinition());
    StructureDefinition foreignDefinition = sampleDefinition();
    foreignDefinition.nodes[1].massKg = 3.0;
    Structure foreign(std::move(foreignDefinition));
    checkInvalid(
        [&] {
            static_cast<void>(
                buildStructureFrame(foreign, mapping, sampleContext()));
        },
        "mapping: mapping cannot be used with a foreign structure");

    StructureFrameContext invalidContext = sampleContext();
    invalidContext.sceneChecksum.clear();
    checkInvalid(
        [&] {
            static_cast<void>(
                buildStructureFrame(structure, mapping, invalidContext));
        },
        "mapping: invalid provenance is rejected before publication");
}

void testSceneAssemblyMappingAndReordering() {
    using namespace simwing::fsi;
    const Scene scene = sampleAssemblyScene();
    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    check(assembly.ok(), "scene mapping: fixture assembly succeeds");
    Structure structure(assembly.definition);
    const StructureFrameMapping mapping =
        makeStructureFrameMapping(scene, assembly, structure);

    check(mapping.vertexStableIds()
              == std::vector<StableId>({10, 11, 12, 13, 50}),
          "scene mapping: surface and junction IDs follow canonical assembly order");
    check(mapping.triangles().size() == 2
              && mapping.triangles()[0].stableId == 500
              && mapping.triangles()[0].negativeRegionId == 90
              && mapping.triangles()[0].positiveRegionId == 91
              && mapping.triangles()[1].stableId == 501
              && mapping.triangles()[1].negativeRegionId == 91
              && mapping.triangles()[1].positiveRegionId == 90,
          "scene mapping: triangle IDs and oriented side regions are exact");
    check(mapping.lines().size() == 3
              && mapping.lines()[0].stableId == 700
              && mapping.lines()[0].role
                     == static_cast<std::uint32_t>(
                         SuspensionLineRole::Brake)
               && mapping.lines()[1].stableId == 701
               && mapping.lines()[1].role
                      == static_cast<std::uint32_t>(
                          SuspensionLineRole::Riser)
              && mapping.lines()[2].stableId == 702,
          "scene mapping: junction graph line IDs and roles are exact");

    Scene reordered = scene;
    reverseSceneCollections(reordered);
    const SceneStructureAssembly reorderedAssembly =
        assembleSceneStructure(reordered);
    check(reorderedAssembly.ok(),
          "scene mapping: reordered fixture assembly succeeds");
    Structure reorderedStructure(reorderedAssembly.definition);
    const StructureFrameMapping reorderedMapping =
        makeStructureFrameMapping(
            reordered, reorderedAssembly, reorderedStructure);

    const StructureFrameContext context = sampleContext();
    const DiagnosticFrame first =
        buildStructureFrame(structure, mapping, context);
    const DiagnosticFrame second = buildStructureFrame(
        reorderedStructure, reorderedMapping, context);
    check(first.vertices.size() == 5
              && first.vertices[4].stableId == 50
              && first.lines.size() == 3
              && first.lines[1].vertex0 == 1
              && first.lines[1].vertex1 == 4
              && first.lines[2].vertex0 == 4
              && first.lines[2].vertex1 == 3,
          "scene mapping: junction node and segmented line endpoints are visible");
    ProtocolError error;
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> secondBytes;
    check(serializeFrame(first, firstBytes, &error)
              && serializeFrame(second, secondBytes, &error)
              && firstBytes == secondBytes,
          "scene mapping: collection reordering produces identical frames");
}

void testSceneAssemblyMappingRejections() {
    using namespace simwing::fsi;
    const Scene scene = sampleAssemblyScene();
    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    Structure structure(assembly.definition);

    Scene staleScene = scene;
    staleScene.vertices[0].positionMeters.z = 0.25;
    checkInvalid(
        [&] {
            static_cast<void>(makeStructureFrameMapping(
                staleScene, assembly, structure));
        },
        "scene mapping: assembly from an older scene is rejected");

    SceneStructureAssembly alteredMapping = assembly;
    std::ranges::reverse(alteredMapping.mappings.triangleIds);
    checkInvalid(
        [&] {
            static_cast<void>(makeStructureFrameMapping(
                scene, alteredMapping, structure));
        },
        "scene mapping: hand-edited stable-ID mapping is rejected");

    SceneStructureAssembly alteredDefinition = assembly;
    alteredDefinition.definition.nodes[0].massKg *= 2.0;
    checkInvalid(
        [&] {
            static_cast<void>(makeStructureFrameMapping(
                scene, alteredDefinition, structure));
        },
        "scene mapping: hand-edited assembly definition is rejected");

    StructureDefinition foreignDefinition = assembly.definition;
    foreignDefinition.nodes[0].massKg *= 2.0;
    Structure foreign(std::move(foreignDefinition));
    checkInvalid(
        [&] {
            static_cast<void>(makeStructureFrameMapping(
                scene, assembly, foreign));
        },
        "scene mapping: Structure from another assembly is rejected");

    Scene invalidScene = scene;
    invalidScene.triangles[0].materialId = 999;
    const SceneStructureAssembly failed =
        assembleSceneStructure(invalidScene);
    check(!failed.ok(),
          "scene mapping: invalid fixture produces a failed assembly");
    checkInvalid(
        [&] {
            static_cast<void>(makeStructureFrameMapping(
                invalidScene, failed, structure));
        },
        "scene mapping: failed assembly is rejected transactionally");
}

void testPilotSuspensionFrameAndReplay() {
    using namespace simwing::fsi;
    const Scene scene = samplePilotAssemblyScene();
    const SceneStructureAssembly assembly = assembleSceneStructure(scene);
    check(assembly.ok(), "pilot frame: scene assembly succeeds");
    Structure structure(assembly.definition);
    const StructureFrameMapping mapping = makeStructureFrameMapping(
        scene, assembly, structure);
    check(mapping.vertexStableIds()
              == std::vector<StableId>({10, 11, 12, 13, 50, 315}),
          "pilot frame: harness is an immutable diagnostic vertex");
    check(mapping.lines().size() == 3
              && mapping.lines()[0].vertex0 == 0
              && mapping.lines()[0].vertex1 == 4
              && mapping.lines()[1].vertex0 == 1
              && mapping.lines()[1].vertex1 == 4
              && mapping.lines()[2].vertex0 == 4
              && mapping.lines()[2].vertex1 == 5,
          "pilot frame: suspension lines retain canopy, junction, and harness endpoints");

    StructureStepSettings settings;
    settings.gravityMetersPerSecondSquared = {};
    settings.constraintIterations =
        assembly.settings.suspensionSolverIterations;
    const auto saved = structure.checkpoint();
    static_cast<void>(structure.step(settings));
    const DiagnosticFrame first = buildStructureFrame(
        structure, mapping, sampleContext());
    structure.restore(saved);
    static_cast<void>(structure.step(settings));
    const DiagnosticFrame replay = buildStructureFrame(
        structure, mapping, sampleContext());
    ProtocolError error;
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> replayBytes;
    check(first.vertices.size() == 6 && first.lines.size() == 3
              && serializeFrame(first, firstBytes, &error)
              && serializeFrame(replay, replayBytes, &error)
              && firstBytes == replayBytes,
          "pilot frame: composite checkpoint continuation serializes bit-identically");
}

} // namespace

int main() {
    testExactInitialFrame();
    testSerializationCompatibility();
    testEvolvingStructureFrames();
    testInvalidMappings();
    testSceneAssemblyMappingAndReordering();
    testSceneAssemblyMappingRejections();
    testPilotSuspensionFrameAndReplay();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing structure-frame check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing structure-frame checks passed");
    return 0;
}
