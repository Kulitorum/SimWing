#include "scene.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

bool contains(const ValidationReport& report, ValidationCode code) {
    return std::any_of(
        report.diagnostics.begin(), report.diagnostics.end(),
        [code](const ValidationDiagnostic& diagnostic) {
            return diagnostic.code == code;
        });
}

Scene validScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:canonical-test-design";
    scene.metadata.exporterVersion = "simwing-scene-test/1";
    scene.metadata.sourceLengthToMeters = 0.001;

    scene.regions = {
        {2, RegionKind::Cell, "cell-1"},
        {1, RegionKind::Outside, "outside"},
    };
    scene.vertices = {
        {13, {0.0, 1.0, 0.0}},
        {11, {1.0, 0.0, 0.0}},
        {10, {0.0, 0.0, 0.0}},
        {12, {1.0, 1.0, 0.0}},
    };
    scene.fabricMaterials = {
        {100, "synthetic-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.01, 2.5e-12},
    };
    scene.triangles = {
        {501, {10, 12, 13}, {{{0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}},
         1, 2, 100, SurfaceRole::Skin},
        {500, {10, 11, 12}, {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}}},
         1, 2, 100, SurfaceRole::Skin},
    };
    scene.openings = {
        {600, {10, 11, 12, 13}, 1, 2, OpeningRole::Intake},
    };
    scene.lineMaterials = {
        {110, "aramid-1mm", 0.001, 0.0008, 18500.0, 1.0},
    };
    scene.pilots = {
        {200, "pilot", 90.0, {0.5, 0.5, -1.0}, {}, {},
         {8.0, 10.0, 6.0}},
    };
    scene.attachments = {
        {301, AttachmentKind::PilotHarness, 0, 200, {0.0, 0.0, 0.25}},
        {300, AttachmentKind::SurfaceVertex, 10, 0, {}},
    };
    scene.suspensionLines = {
        {400, 300, 301, 110, 1.25, SuspensionLineRole::Suspension},
    };
    return scene;
}

std::string encode(const Scene& scene, bool* success = nullptr) {
    std::ostringstream output(std::ios::out | std::ios::binary);
    std::string error;
    const bool result = writeScene(scene, output, &error);
    if (success) {
        *success = result;
    }
    if (!result) {
        std::fprintf(stderr, "encode error: %s\n", error.c_str());
    }
    return output.str();
}

void reverseTopLevelCollections(Scene& scene) {
    std::reverse(scene.regions.begin(), scene.regions.end());
    std::reverse(scene.vertices.begin(), scene.vertices.end());
    std::reverse(scene.fabricMaterials.begin(), scene.fabricMaterials.end());
    std::reverse(scene.triangles.begin(), scene.triangles.end());
    std::reverse(scene.openings.begin(), scene.openings.end());
    std::reverse(scene.lineMaterials.begin(), scene.lineMaterials.end());
    std::reverse(scene.pilots.begin(), scene.pilots.end());
    std::reverse(scene.attachments.begin(), scene.attachments.end());
    std::reverse(scene.suspensionLines.begin(), scene.suspensionLines.end());
}

void testValidRoundTripAndDeterminism() {
    const Scene original = validScene();
    const ValidationReport report = validateScene(original);
    check(report.ok(), "valid scene passes validation");
    check(report.diagnostics.empty(), "valid scene has no diagnostics");

    bool wrote = false;
    const std::string first = encode(original, &wrote);
    check(wrote, "valid scene serializes");
    check(!first.empty(), "valid scene payload is not empty");

    Scene reordered = original;
    reverseTopLevelCollections(reordered);
    bool wroteReordered = false;
    const std::string second = encode(reordered, &wroteReordered);
    check(wroteReordered, "reordered valid scene serializes");
    check(first == second,
          "binary payload is independent of top-level insertion order");

    Scene decoded;
    std::istringstream input(first, std::ios::in | std::ios::binary);
    std::string error;
    check(readScene(input, decoded, &error), "valid payload deserializes");
    check(error.empty(), "valid payload has no read error");
    check(validateScene(decoded).ok(), "round-tripped scene validates");
    check(encode(decoded) == first,
          "deserialize-serialize round trip is byte deterministic");
    check(decoded.metadata.lengthUnit == LengthUnit::Metre
              && decoded.metadata.handedness
                     == CoordinateHandedness::RightHanded
              && decoded.metadata.upAxis == UpAxis::PositiveZ,
          "round trip retains SI right-handed Z-up metadata");
    check(decoded.triangles.front().id == 500,
          "decoded entity collections use canonical stable-ID order");
    check(decoded.triangles.front().negativeSideRegionId == 1
              && decoded.triangles.front().positiveSideRegionId == 2,
          "round trip preserves oriented face side regions");
}

void testDuplicateIdsAndDeterministicDiagnostics() {
    Scene scene = validScene();
    scene.vertices.push_back({10, {2.0, 2.0, 0.0}});
    scene.triangles.push_back(scene.triangles.front());

    const ValidationReport first = validateScene(scene);
    check(!first.ok(), "duplicate IDs reject a scene");
    check(contains(first, ValidationCode::DuplicateId),
          "duplicate IDs produce the expected diagnostic");

    Scene reordered = scene;
    reverseTopLevelCollections(reordered);
    const ValidationReport second = validateScene(reordered);
    const auto key = [](const ValidationReport& report) {
        std::vector<std::tuple<EntityKind, StableId, ValidationCode,
                               std::string>> result;
        for (const ValidationDiagnostic& diagnostic : report.diagnostics) {
            result.emplace_back(diagnostic.entityKind, diagnostic.entityId,
                                diagnostic.code, diagnostic.message);
        }
        return result;
    };
    check(key(first) == key(second),
          "validation diagnostics are deterministic under reordering");
}

void testNonFiniteAndDegenerateFaces() {
    Scene nonFinite = validScene();
    nonFinite.vertices.front().positionMeters.x =
        std::numeric_limits<double>::infinity();
    check(contains(validateScene(nonFinite), ValidationCode::NonFiniteValue),
          "non-finite vertices are rejected");

    Scene degenerate = validScene();
    const auto vertex12 = std::find_if(
        degenerate.vertices.begin(), degenerate.vertices.end(),
        [](const Vertex& vertex) { return vertex.id == 12; });
    vertex12->positionMeters = {2.0, 0.0, 0.0};
    const ValidationReport report = validateScene(degenerate);
    check(contains(report, ValidationCode::DegenerateTriangle),
          "collinear face geometry is rejected");
}

void testMissingReferencesAndSideRegions() {
    Scene scene = validScene();
    scene.triangles.front().vertexIds[0] = 9999;
    scene.triangles.front().materialId = 9998;
    scene.triangles.front().negativeSideRegionId = 9997;
    const ValidationReport missing = validateScene(scene);
    check(contains(missing, ValidationCode::MissingVertexReference),
          "missing face vertices are rejected");
    check(contains(missing, ValidationCode::MissingMaterialReference),
          "missing face materials are rejected");
    check(contains(missing, ValidationCode::MissingRegionReference),
          "missing face side regions are rejected");

    Scene sameSide = validScene();
    sameSide.triangles.front().negativeSideRegionId = 2;
    check(contains(validateScene(sameSide),
                   ValidationCode::InvalidSideRegions),
          "identical face side regions are rejected");
}

void testOpeningsAttachmentsAndLines() {
    Scene invalidOpening = validScene();
    invalidOpening.openings.front().orderedVertexIds = {10, 11, 10};
    check(contains(validateScene(invalidOpening),
                   ValidationCode::InvalidOpening),
          "opening loops require unique boundary vertices");

    Scene missingPilot = validScene();
    missingPilot.attachments.front().pilotId = 9999;
    check(contains(validateScene(missingPilot),
                   ValidationCode::MissingPilotReference),
          "pilot harness attachment requires an existing pilot");

    Scene dangling = validScene();
    dangling.suspensionLines.front().endAttachmentId = 9999;
    check(contains(validateScene(dangling),
                   ValidationCode::DanglingLineAttachment),
          "dangling line attachments are rejected");

    Scene missingLineMaterial = validScene();
    missingLineMaterial.suspensionLines.front().materialId = 9999;
    check(contains(validateScene(missingLineMaterial),
                   ValidationCode::MissingMaterialReference),
          "suspension lines require an existing line material");
}

void testInvalidWriteAndMalformedRead() {
    Scene invalid = validScene();
    invalid.triangles.front().materialId = 9999;
    std::ostringstream rejected(std::ios::out | std::ios::binary);
    std::string error;
    check(!writeScene(invalid, rejected, &error),
          "writer rejects an invalid scene");
    check(rejected.str().empty(),
          "writer emits nothing before validation succeeds");
    check(!error.empty(), "invalid scene write explains its failure");

    const std::string complete = encode(validScene());
    const std::string truncated = complete.substr(0, complete.size() - 7);
    std::istringstream input(truncated, std::ios::in | std::ios::binary);
    Scene destination = validScene();
    destination.metadata.designChecksum = "destination-sentinel";
    check(!readScene(input, destination, &error),
          "reader rejects a truncated payload");
    check(destination.metadata.designChecksum == "destination-sentinel",
          "failed read does not partially replace the destination scene");

    Scene oversized = validScene();
    oversized.metadata.exporterVersion.assign(1'048'577, 'x');
    std::ostringstream oversizedOutput(std::ios::out | std::ios::binary);
    check(!writeScene(oversized, oversizedOutput, &error),
          "writer rejects strings beyond the binary safety limit");
    check(oversizedOutput.str().empty(),
          "binary safety preflight fails before writing partial output");
}

} // namespace

int main() {
    testValidRoundTripAndDeterminism();
    testDuplicateIdsAndDeterministicDiagnostics();
    testNonFiniteAndDegenerateFaces();
    testMissingReferencesAndSideRegions();
    testOpeningsAttachmentsAndLines();
    testInvalidWriteAndMalformedRead();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing scene check(s) failed\n", failures);
        return 1;
    }
    std::printf("all SimWing scene checks passed\n");
    return 0;
}
