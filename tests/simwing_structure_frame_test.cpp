#include "structure_frame.h"

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
    checkInvalid(
        [&] { StructureFrameMapping mapping(structure, invalid); },
        "mapping: equal membrane side regions are rejected");

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

} // namespace

int main() {
    testExactInitialFrame();
    testSerializationCompatibility();
    testEvolvingStructureFrames();
    testInvalidMappings();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing structure-frame check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing structure-frame checks passed");
    return 0;
}
