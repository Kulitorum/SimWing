#include "engine_paths.h"
#include "input_migration.h"
#include "nurbs_model.h"
#include "scene_structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <system_error>
#include <vector>

extern "C" int MAIN__();
extern "C" void f_exit();

namespace {

int failures = 0;

void check(bool condition, const char *message)
{
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

std::string pathToUtf8(const std::filesystem::path &path)
{
    const auto encoded = path.u8string();
    return {reinterpret_cast<const char *>(encoded.data()), encoded.size()};
}

class CurrentPathGuard
{
public:
    CurrentPathGuard()
        : original_(std::filesystem::current_path())
    {
    }

    ~CurrentPathGuard()
    {
        std::error_code error;
        std::filesystem::current_path(original_, error);
    }

private:
    std::filesystem::path original_;
};

lep::SimWingFabricExportSettings fabric(
    simwing::fsi::SurfaceRole role, const std::string &name)
{
    lep::SimWingFabricExportSettings value;
    value.role = role;
    value.name = name;
    value.warpStiffnessNewtonsPerMeter = 1000.0;
    value.weftStiffnessNewtonsPerMeter = 900.0;
    value.shearStiffnessNewtonsPerMeter = 100.0;
    value.bendingStiffnessNewtonMeters = 0.001;
    value.arealDensityKgPerSquareMeter = 0.04;
    value.dampingSeconds = 0.02;
    value.porosityFraction = 0.01;
    value.permeabilitySquareMeters = 1.0e-12;
    return value;
}

lep::SimWingLineExportSettings line(
    const std::string &captureTypeName,
    double diameterMeters)
{
    lep::SimWingLineExportSettings value;
    value.captureTypeName = captureTypeName;
    value.name = "synthetic fixture " + captureTypeName;
    value.diameterMeters = diameterMeters;
    value.linearDensityKgPerMeter = 0.002;
    value.axialStiffnessNewtons = 20000.0;
    value.dragCoefficient = 1.0;
    return value;
}

lep::SimWingSceneExportSettings settings()
{
    using simwing::fsi::SurfaceRole;
    lep::SimWingSceneExportSettings value;
    value.designChecksum = "fixture:3.28/gnuC2-27";
    value.exporterVersion = "real-capture-test/1";
    value.fabricMaterials = {
        fabric(SurfaceRole::Skin, "synthetic fixture skin"),
        fabric(SurfaceRole::Rib, "synthetic fixture rib"),
        fabric(SurfaceRole::Diagonal, "synthetic fixture diagonal"),
        fabric(SurfaceRole::MiniRib, "synthetic fixture mini-rib"),
    };
    value.lineMaterials = {
        line("Riser", 0.005),
        line("Line275", 0.0019),
        line("Line160", 0.0014),
        line("Line120", 0.00115),
        line("Line200B", 0.002),
    };
    value.pilot.name = "synthetic fixture pilot";
    value.pilot.massKg = 70.0;
    value.pilot.centerOfMassPositionMeters = {0.0, 0.842486, -7.99559};
    value.pilot.bodyToWorld = {1.0, 0.0, 0.0, 0.0};
    value.pilot.principalInertiaKgSquareMeters = {8.0, 10.0, 6.0};
    value.pilot.harnessPoints = {
        {{0.210, 0.842486, -7.79559},
         {0.210, 0.0, 0.2}},
        {{-0.210, 0.842486, -7.79559},
         {-0.210, 0.0, 0.2}},
    };
    value.pilot.endpointMatchToleranceMeters = 0.002;
    value.suspensionJunctionMassKg = 0.02;
    // This diagnostic fixture uses the exporter's existing vertex-only
    // attachment contract. The coarse captured skin is up to about 12 mm
    // from a few authored brake-anchor endpoints; a production design source
    // must eventually author exact attachment vertices.
    value.surfaceEndpointMatchToleranceMeters = 0.015;
    return value;
}

bool runRealCapture(const std::filesystem::path &input,
                    const std::filesystem::path &output)
{
    std::filesystem::create_directories(output);
    PreparedInput prepared = PreparedInput::forVersion328(input, output);
    const std::string inputUtf8 = pathToUtf8(prepared.path());
    const std::string outputUtf8 = pathToUtf8(output);
    lep_configure_paths(inputUtf8.c_str(), outputUtf8.c_str());

    CurrentPathGuard currentPath;
    std::filesystem::current_path(input.parent_path());
    lep::resetNurbsModel();
    const int result = MAIN__();
    f_exit();
    return result == 0;
}

void testRealDesignCapture(const std::filesystem::path &input,
                           const std::filesystem::path &output)
{
    check(runRealCapture(input, output),
          "the real 3.28 fixture calculation succeeds");
    if (failures != 0) {
        return;
    }

    const lep::SimWingSceneExportResult result =
        lep::buildSimWingScene(settings());
    if (!result.success) {
        for (const std::string &error : result.errors) {
            std::fprintf(stderr, "export: %s\n", error.c_str());
        }
    }
    check(result.success,
          "the real calculation capture exports a valid scene-v2.1 scene");
    if (!result.success) {
        return;
    }

    const auto roleCount = [&result](simwing::fsi::SurfaceRole role) {
        return std::count_if(
            result.scene.triangles.begin(), result.scene.triangles.end(),
            [role](const simwing::fsi::Triangle &triangle) {
                return triangle.role == role;
            });
    };
    const auto attachmentCount =
        [&result](simwing::fsi::AttachmentKind kind) {
            return std::count_if(
                result.scene.attachments.begin(),
                result.scene.attachments.end(),
                [kind](const simwing::fsi::Attachment &attachment) {
                    return attachment.kind == kind;
                });
        };

    check(result.validation.ok(), "the real exported scene validates");
    check(roleCount(simwing::fsi::SurfaceRole::Skin) > 0,
          "real export contains skin triangles");
    check(roleCount(simwing::fsi::SurfaceRole::Rib) > 0,
          "real export contains triangulated ribs");
    check(roleCount(simwing::fsi::SurfaceRole::Diagonal) > 0,
          "real export contains captured internal sheets");
    check(roleCount(simwing::fsi::SurfaceRole::MiniRib) > 0,
          "real export zipper-triangulates captured mini-ribs");
    check(!result.scene.openings.empty(),
          "real export contains authored intake and crossport openings");
    std::set<std::pair<simwing::fsi::StableId, simwing::fsi::StableId>>
        triangleEdges;
    for (const auto &triangle : result.scene.triangles) {
        for (std::size_t corner = 0; corner < 3; ++corner) {
            auto edge = std::pair{
                triangle.vertexIds[corner],
                triangle.vertexIds[(corner + 1) % 3]};
            if (edge.second < edge.first) {
                std::swap(edge.first, edge.second);
            }
            triangleEdges.insert(edge);
        }
    }
    const bool crossportsUseMeshEdges = std::ranges::all_of(
        result.scene.openings, [&](const auto &opening) {
            if (opening.role != simwing::fsi::OpeningRole::Crossport) {
                return true;
            }
            for (std::size_t index = 0;
                 index < opening.orderedVertexIds.size(); ++index) {
                auto edge = std::pair{
                    opening.orderedVertexIds[index],
                    opening.orderedVertexIds[
                        (index + 1) % opening.orderedVertexIds.size()]};
                if (edge.second < edge.first) {
                    std::swap(edge.first, edge.second);
                }
                if (!triangleEdges.contains(edge)) {
                    return false;
                }
            }
            return true;
        });
    check(crossportsUseMeshEdges,
          "real crossport loops reuse exact rib-triangulation boundary edges");
    check(result.scene.suspensionLines.size() == 190,
          "real export preserves every captured suspension segment");
    check(attachmentCount(simwing::fsi::AttachmentKind::PilotHarness) == 2,
          "real export connects both captured riser roots to the pilot");
    check(attachmentCount(simwing::fsi::AttachmentKind::SurfaceVertex) > 0,
          "real export connects terminal lines to authoritative surface vertices");
    check(!result.scene.suspensionJunctions.empty(),
          "real export preserves intermediate suspension junctions");

    const simwing::fsi::SceneStructureAssembly assembly =
        simwing::fsi::assembleSceneStructure(result.scene);
    if (!assembly.ok()) {
        std::set<simwing::fsi::StableId> triangleVertices;
        std::set<simwing::fsi::StableId> openingVertices;
        std::set<simwing::fsi::StableId> attachmentVertices;
        for (const auto &triangle : result.scene.triangles) {
            triangleVertices.insert(triangle.vertexIds.begin(),
                                    triangle.vertexIds.end());
        }
        for (const auto &opening : result.scene.openings) {
            openingVertices.insert(opening.orderedVertexIds.begin(),
                                   opening.orderedVertexIds.end());
        }
        for (const auto &attachment : result.scene.attachments) {
            if (attachment.kind
                == simwing::fsi::AttachmentKind::SurfaceVertex) {
                attachmentVertices.insert(attachment.vertexId);
            }
        }
        std::size_t orphanCount = 0;
        std::size_t orphanOpenings = 0;
        std::size_t orphanAttachments = 0;
        for (const auto &vertex : result.scene.vertices) {
            if (!triangleVertices.contains(vertex.id)) {
                ++orphanCount;
                orphanOpenings += openingVertices.contains(vertex.id);
                orphanAttachments += attachmentVertices.contains(vertex.id);
            }
        }
        std::fprintf(stderr,
                     "scene vertices %zu, triangle-used %zu, orphan %zu, orphan-opening %zu, orphan-attachment %zu\n",
                     result.scene.vertices.size(), triangleVertices.size(),
                     orphanCount, orphanOpenings, orphanAttachments);
        std::map<simwing::fsi::SceneStructureDiagnosticCode, std::size_t>
            diagnosticCounts;
        for (const auto &diagnostic : assembly.diagnostics) {
            ++diagnosticCounts[diagnostic.code];
        }
        for (const auto &[code, count] : diagnosticCounts) {
            std::fprintf(stderr, "assembly diagnostic %u: %zu\n",
                         static_cast<unsigned>(code), count);
        }
        for (std::size_t index = 0;
             index < std::min<std::size_t>(20, assembly.diagnostics.size());
             ++index) {
            const auto &diagnostic = assembly.diagnostics[index];
            std::fprintf(stderr, "assembly: %s\n",
                         diagnostic.message.c_str());
        }
    }
    check(assembly.ok(),
          "real scene assembles through the composite SimWing structure boundary");
    if (!assembly.ok()) {
        return;
    }
    check(assembly.definition.suspension.has_value()
              && assembly.definition.suspension->segments.size() == 190
              && assembly.mappings.suspensionSegmentLineIds.size() == 190
              && assembly.mappings.pilotHarnessAttachmentIds.size() == 2,
          "real assembly retains all suspension segments and both harness roots");

    simwing::fsi::Structure structure(assembly.definition);
    const simwing::viewer::StructureFrameMapping mapping =
        simwing::viewer::makeStructureFrameMapping(
            result.scene, assembly, structure);
    simwing::fsi::StructureStepSettings step;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    step.constraintIterations =
        assembly.settings.suspensionSolverIterations;
    const simwing::fsi::StructureCheckpoint saved = structure.checkpoint();
    const simwing::fsi::StructureDiagnostics firstDiagnostics =
        structure.step(step);
    simwing::viewer::StructureFrameContext context;
    context.sceneChecksum = result.scene.metadata.designChecksum;
    context.solverCommit = "real-structural-worker-test/1";
    context.timeStepSeconds = step.timeStepSeconds;
    const simwing::viewer::DiagnosticFrame firstFrame =
        simwing::viewer::buildStructureFrame(structure, mapping, context);
    structure.restore(saved);
    const simwing::fsi::StructureDiagnostics replayDiagnostics =
        structure.step(step);
    const simwing::viewer::DiagnosticFrame replayFrame =
        simwing::viewer::buildStructureFrame(structure, mapping, context);
    simwing::viewer::ProtocolError protocolError;
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> replayBytes;
    check(firstDiagnostics.finite
              && firstDiagnostics == replayDiagnostics
              && simwing::viewer::serializeFrame(
                     firstFrame, firstBytes, &protocolError)
              && simwing::viewer::serializeFrame(
                     replayFrame, replayBytes, &protocolError)
              && firstBytes == replayBytes,
          "real structural step and composite checkpoint replay are bit-identical");
    check(firstFrame.lines.size() == 190
              && firstFrame.vertices.size()
                     == assembly.definition.nodes.size() + 2,
          "real diagnostic frame contains all lines and rigid harness vertices");

    std::stringstream trace(
        std::ios::in | std::ios::out | std::ios::binary);
    simwing::viewer::TraceWriter writer(trace);
    const simwing::viewer::TraceHeader header{
        result.scene.metadata.designChecksum,
        context.solverCommit};
    check(writer.writeHeader(header)
              && writer.writeFrame(firstFrame)
              && writer.finish(),
          "real accepted structure state writes a completed viewer trace");
    trace.seekg(0);
    simwing::viewer::TraceReader reader(trace);
    simwing::viewer::TraceHeader decodedHeader;
    simwing::viewer::DiagnosticFrame decodedFrame;
    check(reader.readHeader(decodedHeader)
              && decodedHeader.sceneChecksum == header.sceneChecksum
              && decodedHeader.solverCommit == header.solverCommit
              && reader.readNext(decodedFrame)
                     == simwing::viewer::TraceReadStatus::Frame
              && decodedFrame.lines.size() == 190
              && reader.readNext(decodedFrame)
                     == simwing::viewer::TraceReadStatus::End,
          "real viewer trace replays one complete accepted frame");
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <design> <output-directory>\n", argv[0]);
        return 2;
    }
    try {
        testRealDesignCapture(
            std::filesystem::absolute(std::filesystem::u8path(argv[1]))
                .lexically_normal(),
            std::filesystem::absolute(std::filesystem::u8path(argv[2]))
                .lexically_normal());
    } catch (const std::exception &exception) {
        std::fprintf(stderr, "FAIL: %s\n", exception.what());
        ++failures;
    }
    lep::resetNurbsModel();
    if (failures == 0) {
        std::printf("simwing real model scene export tests passed\n");
    }
    return failures == 0 ? 0 : 1;
}
