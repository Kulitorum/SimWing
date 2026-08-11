#include "engine_paths.h"
#include "input_migration.h"
#include "nurbs_model.h"
#include "scene_fluid_cell_volume.h"
#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_mimetic_condensed_trace_system.h"
#include "scene_fluid_mimetic_control_cell.h"
#include "scene_fluid_mimetic_pressure_audit.h"
#include "scene_fluid_mimetic_pressure_sampling.h"
#include "scene_fluid_mimetic_pressure_solve.h"
#include "scene_fluid_mimetic_pressure_state.h"
#include "scene_fluid_mimetic_pressure_state_persistence.h"
#include "scene_fluid_pressure_shadow_comparison.h"
#include "scene_fluid_pressure_owner_transition.h"
#include "scene_fluid_mimetic_trace_flow.h"
#include "scene_fluid_mimetic_trace_solve.h"
#include "scene_fluid_mimetic_trace_system.h"
#include "scene_fluid_opening_cap.h"
#include "scene_fluid_opening_face_crossing.h"
#include "scene_fluid_pressure_control_volume.h"
#include "scene_fluid_pressure_face_link.h"
#include "scene_fluid_surface.h"
#include "scene_structure.h"
#include "structure_frame.h"
#include "structure_rest_audit.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
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
    value.collapsedBoundarySeam = {
        "synthetic fixture collapsed wingtip seam", 0.001, 5000.0};
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
          "the real calculation capture exports a valid scene-v2.2 scene");
    if (!result.success) {
        return;
    }
    const std::filesystem::path scenePath =
        output / "simwing-scene-v2.bin";
    std::ofstream sceneOutput(
        scenePath, std::ios::binary | std::ios::trunc);
    std::string sceneWriteError;
    check(sceneOutput
              && simwing::fsi::writeScene(
                  result.scene, sceneOutput, &sceneWriteError),
          "real scene-v2 export writes a reusable worker input payload");
    sceneOutput.close();

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
    const std::size_t intakeCount = std::count_if(
        result.scene.openings.begin(), result.scene.openings.end(),
        [](const auto &opening) {
            return opening.role == simwing::fsi::OpeningRole::Intake;
        });
    const bool intakesHaveBoundaryDisks = std::ranges::all_of(
        result.scene.openings, [](const auto &opening) {
            return opening.role != simwing::fsi::OpeningRole::Intake
                || (opening.orderedVertexIds.size() >= 3
                    && opening.capTriangleVertexIds.size() + 2
                        == opening.orderedVertexIds.size());
        });
    check(intakeCount > 0 && intakesHaveBoundaryDisks,
          "real intakes carry explicit boundary-vertex cap disks");
    const auto vertexPosition = [&result](simwing::fsi::StableId id) {
        const auto vertex = std::find_if(
            result.scene.vertices.begin(), result.scene.vertices.end(),
            [id](const auto &candidate) { return candidate.id == id; });
        return vertex->positionMeters;
    };
    const bool collapsedSeamsArePaired =
        result.scene.seams.size() == 2
        && result.scene.seamMaterials.size() == 1
        && std::ranges::all_of(
            result.scene.seams,
            [&](const simwing::fsi::Seam &seam) {
                if (seam.firstOrderedVertexIds.size() != 33
                    || seam.secondOrderedVertexIds.size() != 33
                    || seam.firstOrderedVertexIds.front()
                        != seam.secondOrderedVertexIds.front()
                    || seam.firstOrderedVertexIds.back()
                        != seam.secondOrderedVertexIds.back()) {
                    return false;
                }
                for (std::size_t index = 0;
                     index < seam.firstOrderedVertexIds.size(); ++index) {
                    const auto first = vertexPosition(
                        seam.firstOrderedVertexIds[index]);
                    const auto second = vertexPosition(
                        seam.secondOrderedVertexIds[index]);
                    if (first.x != second.x || first.y != second.y
                        || first.z != second.z
                        || (index != 0
                            && index + 1
                                != seam.firstOrderedVertexIds.size()
                            && seam.firstOrderedVertexIds[index]
                                == seam.secondOrderedVertexIds[index])) {
                        return false;
                    }
                }
                return true;
            });
    check(collapsedSeamsArePaired,
          "real collapsed wingtips become two exact 33-pair sewn chains with shared endpoints");
    const auto subtract = [](const auto &first, const auto &second) {
        return simwing::fsi::Vec3{
            first.x - second.x,
            first.y - second.y,
            first.z - second.z};
    };
    const auto cross = [](const auto &first, const auto &second) {
        return simwing::fsi::Vec3{
            first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
    };
    const auto dot = [](const auto &first, const auto &second) {
        return first.x * second.x + first.y * second.y
            + first.z * second.z;
    };
    bool intakeCapsAreOriented = true;
    bool hasNonPlanarIntake = false;
    for (const auto &opening : result.scene.openings) {
        if (opening.role != simwing::fsi::OpeningRole::Intake) {
            continue;
        }
        simwing::fsi::Vec3 newell;
        for (std::size_t index = 0;
             index < opening.orderedVertexIds.size(); ++index) {
            const auto current = vertexPosition(
                opening.orderedVertexIds[index]);
            const auto next = vertexPosition(
                opening.orderedVertexIds[
                    (index + 1) % opening.orderedVertexIds.size()]);
            newell.x += (current.y - next.y) * (current.z + next.z);
            newell.y += (current.z - next.z) * (current.x + next.x);
            newell.z += (current.x - next.x) * (current.y + next.y);
        }
        const double normalLength =
            std::hypot(newell.x, newell.y, newell.z);
        intakeCapsAreOriented = intakeCapsAreOriented
            && normalLength > 0.0;
        if (!(normalLength > 0.0)) {
            continue;
        }
        const simwing::fsi::Vec3 normal{
            newell.x / normalLength,
            newell.y / normalLength,
            newell.z / normalLength};
        const auto planePoint = vertexPosition(
            opening.orderedVertexIds.front());
        for (const auto id : opening.orderedVertexIds) {
            hasNonPlanarIntake = hasNonPlanarIntake
                || std::abs(dot(
                       subtract(vertexPosition(id), planePoint), normal))
                    > 1.0e-10;
        }
        for (const auto &triangle : opening.capTriangleVertexIds) {
            const auto first = vertexPosition(triangle[0]);
            const auto second = vertexPosition(triangle[1]);
            const auto third = vertexPosition(triangle[2]);
            const auto areaVector = cross(
                subtract(second, first), subtract(third, first));
            intakeCapsAreOriented = intakeCapsAreOriented
                && dot(areaVector, normal) > 1.0e-18;
        }
    }
    check(intakeCapsAreOriented && hasNonPlanarIntake,
          "real intake disks are nonplanar, nondegenerate, and consistently oriented");
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
    const bool openingsUseMeshEdges = std::ranges::all_of(
        result.scene.openings, [&](const auto &opening) {
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
    check(openingsUseMeshEdges,
          "real intake and crossport loops reuse exact fabric-mesh boundary edges");
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
              && assembly.mappings.pilotHarnessAttachmentIds.size() == 2
              && assembly.definition.constraints.size() == 190
              && assembly.mappings.constraintSeamRanges.size() == 2
              && assembly.mappings.constraintSeamRanges[0].constraintCount
                  == 95
              && assembly.mappings.constraintSeamRanges[1].constraintCount
                  == 95,
          "real assembly retains suspension and both 95-constraint wingtip seams");

    simwing::fsi::Structure structure(assembly.definition);
    const simwing::fsi::SceneFluidSurfaceAssembly fluidSurface =
        simwing::fsi::assembleSceneFluidSurface(result.scene);
    if (!fluidSurface.ok()) {
        for (const auto &diagnostic : fluidSurface.diagnostics) {
            std::fprintf(stderr, "fluid surface: %s\n",
                         diagnostic.message.c_str());
        }
    }
    check(fluidSurface.ok(),
          "real scene assembles through the fluid-surface boundary");
    if (!fluidSurface.ok()) {
        return;
    }
    const simwing::fsi::SceneFluidSurfaceState fluidState =
        simwing::fsi::captureSceneFluidSurfaceState(
            fluidSurface.definition, assembly.mappings, structure);
    check(fluidState.vertices.size()
              == fluidSurface.definition.vertices.size()
              && fluidState.fingerprint != 0,
          "every real opening vertex has live Structure-owned fluid motion");
    const simwing::fsi::SceneFluidOpeningCapSet openingCaps =
        simwing::fsi::buildSceneFluidOpeningCaps(
            fluidSurface.definition, fluidState);
    check(openingCaps.caps.size() == result.scene.openings.size()
              && openingCaps.triangles.size()
                  >= result.scene.openings.size()
              && openingCaps.totalAreaSquareMeters > 0.0,
          "real intake and crossport loops produce accepted fluid caps through their junctions");
    const simwing::fsi::SceneFluidSurfaceTransfer fluidTransfer(
        fluidSurface.definition, assembly.mappings, structure);
    simwing::fsi::Vec3 minimum = fluidState.vertices.front().positionMeters;
    simwing::fsi::Vec3 maximum = minimum;
    for (const auto &vertex : fluidState.vertices) {
        minimum.x = std::min(minimum.x, vertex.positionMeters.x);
        minimum.y = std::min(minimum.y, vertex.positionMeters.y);
        minimum.z = std::min(minimum.z, vertex.positionMeters.z);
        maximum.x = std::max(maximum.x, vertex.positionMeters.x);
        maximum.y = std::max(maximum.y, vertex.positionMeters.y);
        maximum.z = std::max(maximum.z, vertex.positionMeters.z);
    }
    constexpr double fluidDomainPaddingMeters = 0.5;
    const simwing::fsi::fluid::PeriodicCartesianGrid fluidGrid(
        {2, 2, 2},
        {minimum.x - fluidDomainPaddingMeters,
         minimum.y - fluidDomainPaddingMeters,
         minimum.z - fluidDomainPaddingMeters},
        {maximum.x + fluidDomainPaddingMeters,
         maximum.y + fluidDomainPaddingMeters,
         maximum.z + fluidDomainPaddingMeters});
    const simwing::fsi::SceneFluidGridEpoch fluidEpoch =
        simwing::fsi::buildSceneFluidGridEpoch(
            fluidSurface.definition, fluidState, fluidGrid, fluidTransfer);
    const simwing::fsi::SceneFluidCellVolumeSet fluidVolumes =
        simwing::fsi::buildSceneFluidCellVolumes(
            fluidSurface.definition, fluidState, fluidGrid, fluidTransfer,
            fluidEpoch);
    const simwing::fsi::SceneFluidOpeningQuadratureSet openingQuadrature =
        simwing::fsi::buildSceneFluidOpeningQuadrature(
            fluidSurface.definition, fluidState, openingCaps);
    const simwing::fsi::SceneFluidOpeningGridPatchSet openingPatches =
        simwing::fsi::buildSceneFluidOpeningGridPatches(
            fluidSurface.definition, fluidState, openingCaps,
            openingQuadrature, fluidGrid);
    const simwing::fsi::SceneFluidOpeningFaceCrossingSet
        openingFaceCrossings =
            simwing::fsi::buildSceneFluidOpeningFaceCrossings(
                fluidSurface.definition, fluidState, openingCaps,
                openingQuadrature, openingPatches, fluidGrid);
    const simwing::fsi::SceneFluidCappedFacePartitionSet
        cappedFacePartitions =
            simwing::fsi::buildSceneFluidCappedFacePartitions(
                fluidSurface.definition, fluidState, fluidGrid,
                fluidTransfer, fluidEpoch, openingCaps,
                openingQuadrature, openingPatches,
                openingFaceCrossings);
    const simwing::fsi::SceneFluidRegionConnectivity fluidConnectivity =
        simwing::fsi::buildSceneFluidRegionConnectivity(
            fluidSurface.definition);
    const simwing::fsi::SceneFluidPressureControlVolumeSet pressureVolumes =
        simwing::fsi::buildSceneFluidPressureControlVolumes(
            fluidSurface.definition, fluidVolumes, fluidConnectivity);
    const simwing::fsi::SceneFluidPressureFaceLinkSet pressureFaceLinks =
        simwing::fsi::buildSceneFluidPressureFaceLinks(
            fluidSurface.definition, fluidState, fluidGrid, fluidTransfer,
            fluidEpoch, openingCaps, openingQuadrature, openingPatches,
            openingFaceCrossings, cappedFacePartitions, fluidVolumes,
            fluidConnectivity, pressureVolumes);
    const simwing::fsi::SceneFluidMimeticControlCellSet mimeticControlCells =
        simwing::fsi::buildSceneFluidMimeticControlCells(
            fluidSurface.definition, fluidState, fluidGrid, fluidEpoch,
            openingCaps, openingQuadrature, openingPatches, pressureVolumes,
            pressureFaceLinks);
    const std::size_t cellOwnedOpeningPatchCount = std::ranges::count(
        openingPatches.patches,
        simwing::fsi::SceneFluidOpeningPatchOwnerKind::Cell,
        &simwing::fsi::SceneFluidOpeningGridPatch::ownerKind);
    const simwing::fsi::SceneFluidMimeticTraceSystem mimeticTraceSystem =
        simwing::fsi::buildSceneFluidMimeticTraceSystem(
            mimeticControlCells);
    simwing::fsi::fluid::MacVelocityField mimeticPredictedVelocity(
        fluidGrid);
    for (std::size_t index = 0;
         index < fluidGrid.cellCount(); ++index) {
        const double sample = static_cast<double>(index + 1);
        mimeticPredictedVelocity.xFaces()[index] = -0.85 + 0.01 * sample;
        mimeticPredictedVelocity.yFaces()[index] = 0.02 * sample;
        mimeticPredictedVelocity.zFaces()[index] = -0.015 * sample;
    }
    const auto mimeticOpeningFlux =
        simwing::fsi::evaluateSceneFluidOpeningFlux(
            fluidSurface.definition, fluidState, openingCaps,
            openingQuadrature, openingPatches, fluidGrid,
            mimeticPredictedVelocity);
    simwing::fsi::SceneFluidMimeticPressureAuditSettings
        realAuditSettings;
    realAuditSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    realAuditSettings.pressureSolve.relativeResidualTolerance = 1.0e-5;
    realAuditSettings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-8;
    realAuditSettings.pressureSolve.maximumIterations = 1000;
    const auto realPressureAudit =
        simwing::fsi::buildSceneFluidMimeticPressureAuditEndpoint(
            fluidSurface.definition, fluidState, fluidGrid, fluidEpoch,
            openingCaps, openingQuadrature, openingPatches,
            pressureVolumes, pressureFaceLinks, mimeticOpeningFlux,
            mimeticPredictedVelocity, realAuditSettings);
    const auto mimeticTraceFlow =
        simwing::fsi::sampleSceneFluidMimeticTraceFlows(
            mimeticControlCells, mimeticTraceSystem, pressureFaceLinks,
            mimeticOpeningFlux, fluidGrid, mimeticPredictedVelocity);
    const auto mimeticSampledSources =
        simwing::fsi::buildSceneFluidMimeticPressureSources(
            mimeticControlCells, mimeticTraceSystem, mimeticTraceFlow);
    check(mimeticTraceFlow.traces.size()
                  == mimeticTraceSystem.sharedTraceCount
              && mimeticTraceFlow.authoredOpeningTraceCount
                  == cellOwnedOpeningPatchCount
              && mimeticTraceFlow.cartesianTraceCount
                      + mimeticTraceFlow.authoredOpeningTraceCount
                  == mimeticTraceSystem.sharedTraceCount
              && mimeticTraceFlow.openingFluxFingerprint
                  == mimeticOpeningFlux.fingerprint
              && mimeticTraceFlow
                    .maximumAbsoluteComponentBalanceResidualCubicMetersPerSecond
                  == 0.0
              && mimeticSampledSources.mimeticTraceFlowFingerprint
                  == mimeticTraceFlow.fingerprint
              && mimeticSampledSources.pressureVolumeRateFingerprint == 0
              && realPressureAudit.controlCells == mimeticControlCells
              && realPressureAudit.fullTraceSystem == mimeticTraceSystem
              && realPressureAudit.predictedTraceFlows == mimeticTraceFlow
              && realPressureAudit.pressureSources
                  == mimeticSampledSources
              && realPressureAudit.pressureEpoch.diagnostics.accepted
              && !realPressureAudit.usesConsecutiveWarmStart
              && !realPressureAudit.usesRegionWallPrediction
              && mimeticSampledSources
                    .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond
                  < 1.0e-12,
          "real mimetic predictor samples every Cartesian and previously rejected opening trace into compatible physical sources");
    const auto condensedTraceSystem =
        simwing::fsi::buildSceneFluidMimeticCondensedTraceSystem(
            mimeticTraceSystem);
    std::vector<double> condensedComponentConstants(
        condensedTraceSystem.traces.size(), 0.0);
    for (const auto& trace : condensedTraceSystem.traces) {
        condensedComponentConstants[trace.traceIndex] =
            1.0 + static_cast<double>(trace.componentIndex);
    }
    const auto condensedConstantAction =
        simwing::fsi::applySceneFluidMimeticCondensedTraceOperator(
            condensedTraceSystem, mimeticTraceSystem,
            condensedComponentConstants);
    double maximumCondensedConstantAction = 0.0;
    for (const double value : condensedConstantAction) {
        maximumCondensedConstantAction = std::max(
            maximumCondensedConstantAction, std::abs(value));
    }
    std::vector<double> componentConstantTraces(
        mimeticTraceSystem.traces.size(), 0.0);
    for (const auto& trace : mimeticTraceSystem.traces) {
        componentConstantTraces[trace.traceIndex] =
            1.0 + static_cast<double>(trace.componentIndex);
    }
    const auto constantTraceAction =
        simwing::fsi::applySceneFluidMimeticTraceOperator(
            mimeticTraceSystem, componentConstantTraces);
    double maximumConstantTraceAction = 0.0;
    for (const double value : constantTraceAction) {
        maximumConstantTraceAction = std::max(
            maximumConstantTraceAction, std::abs(value));
    }
    std::vector<double> manufacturedTraceValues(
        mimeticTraceSystem.traces.size(), 0.0);
    for (const auto& trace : mimeticTraceSystem.traces) {
        manufacturedTraceValues[trace.traceIndex] = std::sin(
            0.001 * static_cast<double>(trace.traceIndex + 1));
    }
    std::vector<double> manufacturedGaugeValues(
        mimeticTraceSystem.componentCount, 0.0);
    for (std::size_t component = 0;
         component < mimeticTraceSystem.componentCount; ++component) {
        manufacturedGaugeValues[component] = manufacturedTraceValues[
            mimeticTraceSystem.componentGaugeTraceIndices[component]];
    }
    for (const auto& trace : mimeticTraceSystem.traces) {
        manufacturedTraceValues[trace.traceIndex] -=
            manufacturedGaugeValues[trace.componentIndex];
    }
    for (const std::size_t gauge :
         mimeticTraceSystem.componentGaugeTraceIndices) {
        manufacturedTraceValues[gauge] = 0.0;
    }
    const auto manufacturedTraceRightHandSide =
        simwing::fsi::applySceneFluidMimeticTraceOperator(
            mimeticTraceSystem, manufacturedTraceValues);
    std::vector<double> truncatedTraceWarmStart(
        mimeticTraceSystem.traces.size(), 0.0);
    const auto originalTruncatedTraceWarmStart = truncatedTraceWarmStart;
    simwing::fsi::SceneFluidMimeticTraceSolveSettings traceSolveSettings;
    traceSolveSettings.absoluteResidualTolerancePascalsMeters = 1.0e-30;
    traceSolveSettings.relativeResidualTolerance = 0.0;
    traceSolveSettings
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    traceSolveSettings.maximumIterations = 1;
    const auto traceSolveDiagnostics =
        simwing::fsi::solveSceneFluidMimeticTraceSystem(
            mimeticTraceSystem, manufacturedTraceRightHandSide,
            truncatedTraceWarmStart, traceSolveSettings);
    std::vector<double> manufacturedReducedTraceValues(
        condensedTraceSystem.traces.size(), 0.0);
    for (const auto& trace : condensedTraceSystem.traces) {
        manufacturedReducedTraceValues[trace.traceIndex] = std::sin(
            0.001 * static_cast<double>(trace.traceIndex + 1));
    }
    std::vector<double> manufacturedReducedGaugeValues(
        condensedTraceSystem.componentCount, 0.0);
    for (std::size_t component = 0;
         component < condensedTraceSystem.componentCount; ++component) {
        manufacturedReducedGaugeValues[component] =
            manufacturedReducedTraceValues[
                condensedTraceSystem
                    .componentGaugeTraceIndices[component]];
    }
    for (const auto& trace : condensedTraceSystem.traces) {
        manufacturedReducedTraceValues[trace.traceIndex] -=
            manufacturedReducedGaugeValues[trace.componentIndex];
    }
    for (const std::size_t gauge :
         condensedTraceSystem.componentGaugeTraceIndices) {
        manufacturedReducedTraceValues[gauge] = 0.0;
    }
    const auto manufacturedReducedTraceRightHandSide =
        simwing::fsi::applySceneFluidMimeticCondensedTraceOperator(
            condensedTraceSystem, mimeticTraceSystem,
            manufacturedReducedTraceValues);
    std::vector<double> truncatedReducedTraceWarmStart(
        condensedTraceSystem.traces.size(), 0.0);
    const auto originalTruncatedReducedTraceWarmStart =
        truncatedReducedTraceWarmStart;
    auto reducedTraceSolveSettings = traceSolveSettings;
    reducedTraceSolveSettings
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-8;
    const auto reducedTraceSolveDiagnostics =
        simwing::fsi::solveSceneFluidMimeticCondensedTraceSystem(
            condensedTraceSystem, mimeticTraceSystem,
            manufacturedReducedTraceRightHandSide,
            truncatedReducedTraceWarmStart, reducedTraceSolveSettings);
    auto manufacturedCondensedFullTraceValues =
        manufacturedTraceValues;
    std::vector<double> condensedFullGaugeValues(
        condensedTraceSystem.componentCount, 0.0);
    for (std::size_t component = 0;
         component < condensedTraceSystem.componentCount; ++component) {
        const auto& reducedGauge = condensedTraceSystem.traces[
            condensedTraceSystem.componentGaugeTraceIndices[component]];
        condensedFullGaugeValues[component] =
            manufacturedCondensedFullTraceValues[
                reducedGauge.fullTraceIndex];
    }
    for (const auto& trace : mimeticTraceSystem.traces) {
        manufacturedCondensedFullTraceValues[trace.traceIndex] -=
            condensedFullGaugeValues[trace.componentIndex];
    }
    const auto manufacturedCondensedFullRightHandSide =
        simwing::fsi::applySceneFluidMimeticTraceOperator(
            mimeticTraceSystem,
            manufacturedCondensedFullTraceValues);
    const auto manufacturedCondensedRightHandSide =
        simwing::fsi::condenseSceneFluidMimeticTraceRightHandSide(
            condensedTraceSystem, mimeticTraceSystem,
            manufacturedCondensedFullRightHandSide);
    std::vector<double> convergedReducedTraceValues(
        condensedTraceSystem.traces.size(), 0.0);
    auto convergedReducedTraceSolveSettings = reducedTraceSolveSettings;
    convergedReducedTraceSolveSettings.relativeResidualTolerance = 1.0e-5;
    convergedReducedTraceSolveSettings.maximumIterations = 300;
    const auto convergedReducedTraceSolveDiagnostics =
        simwing::fsi::solveSceneFluidMimeticCondensedTraceSystem(
            condensedTraceSystem, mimeticTraceSystem,
            manufacturedCondensedRightHandSide,
            convergedReducedTraceValues,
            convergedReducedTraceSolveSettings);
    const auto reconstructedConvergedFullTraceValues =
        simwing::fsi::reconstructSceneFluidMimeticFullTraces(
            condensedTraceSystem, mimeticTraceSystem,
            manufacturedCondensedFullRightHandSide,
            convergedReducedTraceValues);
    const auto reconstructedConvergedFullAction =
        simwing::fsi::applySceneFluidMimeticTraceOperator(
            mimeticTraceSystem,
            reconstructedConvergedFullTraceValues);
    double maximumReconstructedFullResidual = 0.0;
    for (std::size_t trace = 0;
         trace < reconstructedConvergedFullAction.size(); ++trace) {
        maximumReconstructedFullResidual = std::max(
            maximumReconstructedFullResidual,
            std::abs(reconstructedConvergedFullAction[trace]
                     - manufacturedCondensedFullRightHandSide[trace]));
    }
    std::vector<double> realPredictedContinuityRates(
        mimeticTraceSystem.localOperators.size(), 0.0);
    std::size_t sourceReceiver = 1;
    while (sourceReceiver < mimeticControlCells.controlCells.size()
           && mimeticControlCells.controlCells[sourceReceiver]
                   .componentIndex
               != mimeticControlCells.controlCells[0].componentIndex) {
        ++sourceReceiver;
    }
    simwing::fsi::SceneFluidMimeticPressureSolveResult
        realSourcePressure;
    simwing::fsi::SceneFluidMimeticPressureSourceSet realPhysicalSources;
    simwing::fsi::SceneFluidMimeticPressureState realPressureState;
    simwing::fsi::SceneFluidMimeticPressureSampleSet realPressureSamples;
    bool realPressureTransferCloses = false;
    bool realPressureStatePersists = false;
    bool realPressureComparisonCloses = false;
    if (sourceReceiver < realPredictedContinuityRates.size()) {
        simwing::fsi::SceneFluidMimeticPressureSourceSettings
            realSourceSettings;
        const double sourceContinuityRate = -0.02
            * realSourceSettings.timeStepSeconds
            / realSourceSettings.densityKgPerCubicMeter;
        realPredictedContinuityRates[0] = sourceContinuityRate;
        realPredictedContinuityRates[sourceReceiver] =
            -sourceContinuityRate;
        realPhysicalSources =
            simwing::fsi::buildSceneFluidMimeticPressureSources(
                mimeticControlCells, realPredictedContinuityRates,
                realSourceSettings);
        auto realSourceSolveSettings =
            convergedReducedTraceSolveSettings;
        realSourceSolveSettings.maximumIterations = 1000;
        realSourcePressure =
            simwing::fsi::solveSceneFluidMimeticPressureSystem(
                condensedTraceSystem, mimeticTraceSystem,
                realPhysicalSources,
                std::vector<double>(
                    condensedTraceSystem.traces.size(), 0.0),
                realSourceSolveSettings);
        realPressureState =
            simwing::fsi::captureSceneFluidMimeticPressureState(
                mimeticControlCells, mimeticTraceSystem,
                condensedTraceSystem, realPhysicalSources,
                realSourcePressure);
        realPressureSamples =
            simwing::fsi::sampleSceneFluidMimeticPressure(
                fluidEpoch.quadrature, pressureVolumes,
                mimeticControlCells, mimeticTraceSystem,
                condensedTraceSystem, realPressureState);
        const auto realPressureTransfer =
            simwing::fsi::evaluateSceneFluidMimeticPressureQuadrature(
                fluidSurface.definition, fluidState, fluidTransfer,
                fluidEpoch.quadrature, realPressureSamples);
        realPressureTransferCloses =
            realPressureTransfer.diagnostics().finite
            && realPressureTransfer.diagnostics()
                    .forceResidualNormNewtons < 1.0e-9
            && realPressureTransfer.diagnostics()
                    .momentResidualNormNewtonMeters < 1.0e-9;
        const auto realPressureComparison =
            simwing::fsi::compareSceneFluidPressureShadow(
                fluidSurface.definition, fluidState, fluidTransfer,
                fluidEpoch.quadrature, realPressureSamples,
                realPressureAudit.pressureEpoch.acceptedPressureSamples);
        const auto realPressureZeroComparison =
            simwing::fsi::compareSceneFluidPressureShadow(
                fluidSurface.definition, fluidState, fluidTransfer,
                fluidEpoch.quadrature, realPressureSamples,
                realPressureSamples);
        simwing::fsi::validateSceneFluidPressureShadowComparisonIntegrity(
            realPressureComparison);
        simwing::fsi::validateSceneFluidPressureShadowComparisonIntegrity(
            realPressureZeroComparison);
        simwing::fsi::SceneFluidPressureOwnerTransitionPolicy
            exactComparisonPolicy;
        exactComparisonPolicy.requireSourceComparison = false;
        const auto exactComparisonTransition =
            simwing::fsi::decideSceneFluidPressureOwnerTransition(
                realPressureZeroComparison, exactComparisonPolicy);
        const auto missingSourceTransition =
            simwing::fsi::decideSceneFluidPressureOwnerTransition(
                realPressureZeroComparison);
        simwing::fsi::
            validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
                exactComparisonTransition, realPressureZeroComparison);
        simwing::fsi::
            validateSceneFluidPressureOwnerTransitionDecisionIntegrity(
                missingSourceTransition, realPressureZeroComparison);
        realPressureComparisonCloses =
            realPressureComparison.diagnostics.finite
            && !realPressureComparison.includesSourceComparison
            && !realPressureZeroComparison.includesSourceComparison
            && realPressureComparison.samples.size()
                == fluidEpoch.quadrature.points.size()
            && realPressureComparison.diagnostics
                    .maximumAbsolutePressureDifferenceDeltaPascals > 0.0
            && realPressureZeroComparison.diagnostics
                    .pressureDifferenceDeltaL2Pascals == 0.0
            && realPressureZeroComparison.diagnostics
                    .forceDeltaNormNewtons == 0.0
            && realPressureZeroComparison.diagnostics
                    .momentDeltaNormNewtonMeters == 0.0
            && realPressureZeroComparison.diagnostics
                    .maximumNodalForceDeltaNewtons == 0.0
            && realPressureZeroComparison.diagnostics
                    .bestFitShadowPressureScale > 1.0 - 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .bestFitShadowPressureScale < 1.0 + 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .pressureDifferenceCosineSimilarity > 1.0 - 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .pressureDifferenceCosineSimilarity < 1.0 + 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .bestFitPressureShapeResidualL2Pascals == 0.0
            && realPressureZeroComparison.diagnostics
                    .bestFitShadowNodalForceScale > 1.0 - 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .bestFitShadowNodalForceScale < 1.0 + 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .nodalForceCosineSimilarity > 1.0 - 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .nodalForceCosineSimilarity < 1.0 + 1.0e-15
            && realPressureZeroComparison.diagnostics
                    .bestFitNodalForceShapeResidualL2Newtons == 0.0
            && exactComparisonTransition.selectedOwner
                == simwing::fsi::SceneFluidPressureOwner::ShadowMimetic
            && exactComparisonTransition.rejectionMask == 0
            && missingSourceTransition.selectedOwner
                == simwing::fsi::SceneFluidPressureOwner::ReferenceGraph
            && simwing::fsi::sceneFluidPressureOwnerTransitionRejectedFor(
                missingSourceTransition,
                simwing::fsi::SceneFluidPressureOwnerTransitionRejection::
                    MissingSourceComparison)
            && realPressureComparison.diagnostics.referenceTransfer
                    .forceResidualNormNewtons < 1.0e-8
            && realPressureComparison.diagnostics.shadowTransfer
                    .forceResidualNormNewtons < 1.0e-8
            && realPressureComparison.diagnostics.referenceTransfer
                    .momentResidualNormNewtonMeters < 1.0e-8
            && realPressureComparison.diagnostics.shadowTransfer
                    .momentResidualNormNewtonMeters < 1.0e-8;
        if (!realPressureComparisonCloses) {
            std::fprintf(
                stderr,
                "real pressure comparison: samples=%zu max-delta=%.17g reference-force-residual=%.17g shadow-force-residual=%.17g reference-moment-residual=%.17g shadow-moment-residual=%.17g\n",
                realPressureComparison.samples.size(),
                realPressureComparison.diagnostics
                    .maximumAbsolutePressureDifferenceDeltaPascals,
                realPressureComparison.diagnostics.referenceTransfer
                    .forceResidualNormNewtons,
                realPressureComparison.diagnostics.shadowTransfer
                    .forceResidualNormNewtons,
                realPressureComparison.diagnostics.referenceTransfer
                    .momentResidualNormNewtonMeters,
                realPressureComparison.diagnostics.shadowTransfer
                    .momentResidualNormNewtonMeters);
        }
        std::vector<std::uint8_t> encodedRealPressureState;
        simwing::fsi::SceneFluidMimeticPressureState decodedRealPressureState;
        simwing::fsi::SceneFluidMimeticPressureStatePersistenceError
            pressurePersistenceError;
        realPressureStatePersists =
            simwing::fsi::serializeSceneFluidMimeticPressureState(
                realPressureState, mimeticControlCells,
                mimeticTraceSystem, condensedTraceSystem,
                encodedRealPressureState, &pressurePersistenceError)
            && simwing::fsi::deserializeSceneFluidMimeticPressureState(
                encodedRealPressureState, mimeticControlCells,
                mimeticTraceSystem, condensedTraceSystem,
                decodedRealPressureState, &pressurePersistenceError)
            && decodedRealPressureState == realPressureState;
    }
    const bool mimeticAuditMatches =
        realPressureAudit.condensedTraceSystem == condensedTraceSystem
        && realPressureAudit.pressureEpoch.acceptedPressureSamples
               .bindings.size()
            == fluidEpoch.quadrature.points.size()
        && mimeticControlCells.readyControlCellCount
            == mimeticControlCells.controlCells.size()
        && mimeticControlCells.incompleteTopologyControlCellCount == 0
        && mimeticControlCells.nonclosingControlCellCount == 0
        && mimeticControlCells.unresolvedCartesianFaceCount == 0
        && pressureFaceLinks.unresolvedAmbiguousFaceCount == 0
        && pressureFaceLinks.surfaceClassifiedFullFaceCount == 10
        && mimeticControlCells.omittedZeroVolumeMaterialSideCount == 0
        && mimeticControlCells.missingOpeningControlSideCount == 0
        && mimeticControlCells.openingHalfFaceCount
            == 2 * cellOwnedOpeningPatchCount
        && mimeticTraceSystem.localOperators.size()
            == mimeticControlCells.controlCells.size()
        && mimeticTraceSystem.incidences.size()
            == mimeticControlCells.halfFaces.size()
        && 2 * mimeticTraceSystem.sharedTraceCount
                + mimeticTraceSystem.materialWallTraceCount
            == mimeticTraceSystem.incidences.size()
        && mimeticTraceSystem.localOperatorStorageBytes
            == 7 * mimeticControlCells.halfFaces.size() * sizeof(double)
        && mimeticTraceSystem.sharedTraceCount == 42'927
        && mimeticTraceSystem.materialWallTraceCount == 148'776
        && condensedTraceSystem.traces.size()
            == mimeticTraceSystem.sharedTraceCount
        && condensedTraceSystem.eliminatedMaterialWallTraceCount
            == mimeticTraceSystem.materialWallTraceCount
        && condensedTraceSystem.localCondensations.size()
            == mimeticTraceSystem.localOperators.size()
        && condensedTraceSystem.localCondensationStorageBytes
            == mimeticControlCells.halfFaces.size()
                * (sizeof(std::uint8_t) + 2 * sizeof(double))
        && condensedTraceSystem.minimumPositiveOperatorDiagonal > 0.0
        && maximumCondensedConstantAction < 1.0e-8
        && maximumConstantTraceAction < 1.0e-9
        && traceSolveDiagnostics.compatible
        && !traceSolveDiagnostics.converged
        && traceSolveDiagnostics.finite
        && traceSolveDiagnostics.iterationCount == 1
        && traceSolveDiagnostics.finalResidualL2PascalsMeters
            < traceSolveDiagnostics.initialResidualL2PascalsMeters
        && truncatedTraceWarmStart == originalTruncatedTraceWarmStart
        && reducedTraceSolveDiagnostics.compatible
        && !reducedTraceSolveDiagnostics.converged
        && reducedTraceSolveDiagnostics.finite
        && reducedTraceSolveDiagnostics.iterationCount == 1
        && reducedTraceSolveDiagnostics.finalResidualL2PascalsMeters
            < reducedTraceSolveDiagnostics.initialResidualL2PascalsMeters
        && truncatedReducedTraceWarmStart
            == originalTruncatedReducedTraceWarmStart
        && convergedReducedTraceSolveDiagnostics.compatible
        && convergedReducedTraceSolveDiagnostics.converged
        && convergedReducedTraceSolveDiagnostics.finite
        && convergedReducedTraceSolveDiagnostics.iterationCount > 1
        && convergedReducedTraceSolveDiagnostics.iterationCount <= 300
        && convergedReducedTraceSolveDiagnostics
                .finalResidualL2PascalsMeters
            <= 1.0e-5
                * convergedReducedTraceSolveDiagnostics
                    .initialResidualL2PascalsMeters
        && maximumReconstructedFullResidual < 2.0e-4
        && sourceReceiver < realPredictedContinuityRates.size()
        && realPhysicalSources.fingerprint != 0
        && realPhysicalSources
                .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond
            < 1.0e-18
        && realSourcePressure.diagnostics.accepted
        && realSourcePressure.pressureSourceFingerprint
            == realPhysicalSources.fingerprint
        && realSourcePressure.diagnostics
            .reconstructedFullResidualConverged
        && realSourcePressure.diagnostics.reducedTraceSolve.converged
        && realSourcePressure.diagnostics
                .reconstructedFullResidualMaximumPascalsMeters
            < 2.0e-4
        && realSourcePressure.diagnostics.maximumCellConservationResidual
            < 1.0e-9
        && realPressureState.fingerprint != 0
        && realPressureState.fullTraceSystemFingerprint
            == mimeticTraceSystem.fingerprint
        && realPressureState.condensedTraceSystemFingerprint
            == condensedTraceSystem.fingerprint
        && realPressureState.pressureSourceFingerprint
            == realPhysicalSources.fingerprint
        && realPressureState.controls.size()
            == mimeticControlCells.controlCells.size()
        && realPressureState.traces.size()
            == condensedTraceSystem.traces.size()
        && realPressureSamples.pressureStateFingerprint
            == realPressureState.fingerprint
        && realPressureSamples.bindings.size()
            == fluidEpoch.quadrature.points.size()
        && realPressureSamples.pressures.size()
            == fluidEpoch.quadrature.points.size()
        && realPressureTransferCloses
        && realPressureComparisonCloses
        && realPressureStatePersists;
    if (!mimeticAuditMatches) {
        std::fprintf(
            stderr,
            "real mimetic shell audit: ready=%zu/%zu incomplete=%zu nonclosing=%zu unresolved-faces=%zu [active=%zu ambiguous=%zu classified=%zu opening=%zu capped=%zu] omitted-wall-sides=%zu missing-opening-sides=%zu opening-halves=%zu max-halves/control=%zu traces=%zu shared=%zu walls=%zu compact-bytes=%zu wall-condensation-bytes=%zu condensed-traces=%zu condensed-locals=%zu max-null=%.17g max-condensed-null=%.17g solve-compatible=%d solve-converged=%d solve-finite=%d solve-iterations=%zu solve-compatibility=%.17g solve-initial=%.17g solve-final=%.17g solve-final-max=%.17g reduced-solve-compatible=%d reduced-solve-converged=%d reduced-solve-finite=%d reduced-solve-iterations=%zu reduced-solve-compatibility=%.17g reduced-solve-initial=%.17g reduced-solve-final=%.17g reduced-solve-final-max=%.17g converged-reduced-compatible=%d converged-reduced-converged=%d converged-reduced-finite=%d converged-reduced-iterations=%zu converged-reduced-compatibility=%.17g converged-reduced-initial=%.17g converged-reduced-final=%.17g converged-reduced-final-max=%.17g reconstructed-full-max=%.17g source-pair=%d source-accepted=%d source-full-converged=%d source-reduced-compatible=%d source-reduced-converged=%d source-reduced-iterations=%zu source-reduced-initial=%.17g source-reduced-final=%.17g source-full-tolerance=%.17g source-full-rms=%.17g source-full-max=%.17g source-cell-max=%.17g max-area=%.17g max-moment=%.17g\n",
            mimeticControlCells.readyControlCellCount,
            mimeticControlCells.controlCells.size(),
            mimeticControlCells.incompleteTopologyControlCellCount,
            mimeticControlCells.nonclosingControlCellCount,
            mimeticControlCells.unresolvedCartesianFaceCount,
            pressureFaceLinks.unresolvedActiveFaceCount,
            pressureFaceLinks.unresolvedAmbiguousFaceCount,
            pressureFaceLinks.surfaceClassifiedFullFaceCount,
            pressureFaceLinks.unresolvedOpeningFaceCount,
            pressureFaceLinks.unresolvedCappedFaceCount,
            mimeticControlCells.omittedZeroVolumeMaterialSideCount,
            mimeticControlCells.missingOpeningControlSideCount,
            mimeticControlCells.openingHalfFaceCount,
            mimeticControlCells.maximumHalfFaceCountPerControl,
            mimeticTraceSystem.traces.size(),
            mimeticTraceSystem.sharedTraceCount,
            mimeticTraceSystem.materialWallTraceCount,
            mimeticTraceSystem.localOperatorStorageBytes,
            condensedTraceSystem.localCondensationStorageBytes,
            condensedTraceSystem.traces.size(),
            condensedTraceSystem.localCondensations.size(),
            maximumConstantTraceAction,
            maximumCondensedConstantAction,
            traceSolveDiagnostics.compatible ? 1 : 0,
            traceSolveDiagnostics.converged ? 1 : 0,
            traceSolveDiagnostics.finite ? 1 : 0,
            traceSolveDiagnostics.iterationCount,
            traceSolveDiagnostics
                .maximumAbsoluteComponentCompatibilityPascalsMeters,
            traceSolveDiagnostics.initialResidualL2PascalsMeters,
            traceSolveDiagnostics.finalResidualL2PascalsMeters,
            traceSolveDiagnostics.finalResidualMaximumPascalsMeters,
            reducedTraceSolveDiagnostics.compatible ? 1 : 0,
            reducedTraceSolveDiagnostics.converged ? 1 : 0,
            reducedTraceSolveDiagnostics.finite ? 1 : 0,
            reducedTraceSolveDiagnostics.iterationCount,
            reducedTraceSolveDiagnostics
                .maximumAbsoluteComponentCompatibilityPascalsMeters,
            reducedTraceSolveDiagnostics
                .initialResidualL2PascalsMeters,
            reducedTraceSolveDiagnostics
                .finalResidualL2PascalsMeters,
            reducedTraceSolveDiagnostics
                .finalResidualMaximumPascalsMeters,
            convergedReducedTraceSolveDiagnostics.compatible ? 1 : 0,
            convergedReducedTraceSolveDiagnostics.converged ? 1 : 0,
            convergedReducedTraceSolveDiagnostics.finite ? 1 : 0,
            convergedReducedTraceSolveDiagnostics.iterationCount,
            convergedReducedTraceSolveDiagnostics
                .maximumAbsoluteComponentCompatibilityPascalsMeters,
            convergedReducedTraceSolveDiagnostics
                .initialResidualL2PascalsMeters,
            convergedReducedTraceSolveDiagnostics
                .finalResidualL2PascalsMeters,
            convergedReducedTraceSolveDiagnostics
                .finalResidualMaximumPascalsMeters,
            maximumReconstructedFullResidual,
            sourceReceiver < realPredictedContinuityRates.size() ? 1 : 0,
            realSourcePressure.diagnostics.accepted ? 1 : 0,
            realSourcePressure.diagnostics
                    .reconstructedFullResidualConverged
                ? 1 : 0,
            realSourcePressure.diagnostics.reducedTraceSolve.compatible
                ? 1 : 0,
            realSourcePressure.diagnostics.reducedTraceSolve.converged
                ? 1 : 0,
            realSourcePressure.diagnostics.reducedTraceSolve.iterationCount,
            realSourcePressure.diagnostics.reducedTraceSolve
                .initialResidualL2PascalsMeters,
            realSourcePressure.diagnostics.reducedTraceSolve
                .finalResidualL2PascalsMeters,
            realSourcePressure.diagnostics
                .reconstructedFullResidualTolerancePascalsMeters,
            realSourcePressure.diagnostics
                .reconstructedFullResidualL2PascalsMeters,
            realSourcePressure.diagnostics
                .reconstructedFullResidualMaximumPascalsMeters,
            realSourcePressure.diagnostics
                .maximumCellConservationResidual,
            mimeticControlCells.maximumAreaClosureErrorSquareMeters,
            mimeticControlCells.maximumDivergenceTheoremErrorCubicMeters);
        for (const auto& cell : mimeticControlCells.controlCells) {
            if (!cell.areaVectorClosed || !cell.divergenceTheoremClosed) {
                std::fprintf(
                    stderr,
                    "  nonclosing control=%zu grid-cell=%zu region=%llu halves=%zu cart=%zu wall=%zu opening=%zu incidents=%zu area=%.17g moment=%.17g\n",
                    cell.controlVolumeIndex, cell.cellIndex,
                    static_cast<unsigned long long>(cell.regionId),
                    cell.halfFaceCount, cell.cartesianHalfFaceCount,
                    cell.materialWallHalfFaceCount,
                    cell.openingHalfFaceCount,
                    cell.unresolvedCartesianIncidentCount,
                    cell.maximumAreaClosureErrorSquareMeters,
                    cell.maximumDivergenceTheoremErrorCubicMeters);
            }
        }
    }
    const bool hasBoundaryChainPartition = std::ranges::any_of(
        fluidEpoch.facePartitions.partitions,
        [](const simwing::fsi::fluid::SceneFluidFacePartition &partition) {
            return partition.kind
                == simwing::fsi::fluid::SceneFluidFacePartitionKind::
                    BoundaryOpenChain;
        });
    const auto unresolvedCappedFaceRecordCount = std::ranges::count(
        cappedFacePartitions.faces,
        simwing::fsi::invalidSceneFluidCappedFacePartitionIndex,
        &simwing::fsi::SceneFluidCappedFace::partitionIndex);
    const auto unpairedMaterialFaceCount = std::ranges::count(
        cappedFacePartitions.faces,
        simwing::fsi::SceneFluidCappedFaceStatus::
            UnpairedMaterialEndpoint,
        &simwing::fsi::SceneFluidCappedFace::status);
    std::set<simwing::fsi::StableId> cappedFailureSourceIds;
    for (const auto& face : cappedFacePartitions.faces) {
        if (face.failureSourceStableId != simwing::fsi::invalidStableId) {
            cappedFailureSourceIds.insert(face.failureSourceStableId);
        }
    }
    std::set<simwing::fsi::StableId> rejectedEmbeddedOpeningIds;
    using EmbeddedRejectionStatus =
        simwing::fsi::SceneFluidEmbeddedOpeningRejectionStatus;
    const bool embeddedRejectionsAreNonAdmissible =
        std::ranges::all_of(
            pressureFaceLinks.embeddedOpeningRejections,
            [&](const simwing::fsi::SceneFluidEmbeddedOpeningRejection&
                    rejection) {
                rejectedEmbeddedOpeningIds.insert(rejection.openingId);
                return rejection.openingPatchStableId != 0
                    && rejection.status
                        == EmbeddedRejectionStatus::
                            NonPositiveProjectedDistance
                    && rejection.projectedCenterDistanceMeters < 0.0
                    && rejection.negativeCentroidSignedDistanceMeters > 0.0
                    && rejection.positiveCentroidSignedDistanceMeters > 0.0;
            });
    const bool embeddedOneRingSupportIsTwoSided =
        std::ranges::all_of(
            pressureFaceLinks.embeddedOpeningRejections,
            [&](const simwing::fsi::SceneFluidEmbeddedOpeningRejection&
                    rejection) {
                if (rejection.oneRingStatus
                        != simwing::fsi::
                            SceneFluidEmbeddedOpeningOneRingStatus::BothSides
                    || rejection.negativeAdmissibleOneRingSupportCount == 0
                    || rejection.positiveAdmissibleOneRingSupportCount == 0
                    || rejection.firstOneRingSupport
                        > pressureFaceLinks
                            .embeddedOpeningOneRingSupports.size()
                    || rejection.oneRingSupportCount
                        > pressureFaceLinks
                                .embeddedOpeningOneRingSupports.size()
                            - rejection.firstOneRingSupport) {
                    return false;
                }
                for (std::size_t offset = 0;
                     offset < rejection.oneRingSupportCount; ++offset) {
                    const auto& support = pressureFaceLinks
                        .embeddedOpeningOneRingSupports[
                            rejection.firstOneRingSupport + offset];
                    if (support.rejectionIndex != rejection.rejectionIndex
                        || support.cartesianFaceLinkIndex
                            >= pressureFaceLinks.links.size()) {
                        return false;
                    }
                    const auto& link = pressureFaceLinks.links[
                        support.cartesianFaceLinkIndex];
                    if (link.stableId
                            != support.cartesianFaceLinkStableId
                        || link.faceIndex
                            >= pressureFaceLinks.faces.size()
                        || link.kind
                            != simwing::fsi::
                                SceneFluidPressureFaceLinkKind::SameRegion
                        || link.geometryKind
                            != simwing::fsi::
                                SceneFluidPressureLinkGeometryKind::
                                    CartesianFace) {
                        return false;
                    }
                    const bool rootIsMinus =
                        link.minusControlVolumeIndex
                            == support.rootControlVolumeIndex;
                    if ((!rootIsMinus
                            && link.plusControlVolumeIndex
                                != support.rootControlVolumeIndex)
                        || (rootIsMinus
                                ? link.plusControlVolumeIndex
                                : link.minusControlVolumeIndex)
                            != support.donorControlVolumeIndex) {
                        return false;
                    }
                }
                return true;
            });
    check(fluidVolumes.regionVolumes.size() == result.scene.regions.size()
              && fluidVolumes.openingCapCount
                  == result.scene.openings.size()
              && fluidEpoch.faceGraph.higherDegreeNodeCount > 0
              && fluidEpoch.facePartitions.unresolvedActiveFaceCount > 0
              && hasBoundaryChainPartition
              && fluidEpoch.facePartitions.partitions.size() == 1
              && fluidEpoch.facePartitions.unresolvedActiveFaceCount
                      + fluidEpoch.facePartitions.partitions.size()
                  == fluidEpoch.faceTopology.activeFaces.size()
              && fluidVolumes.maximumCellVolumeResidualCubicMeters
                  < 1.0e-10
              && fluidVolumes.maximumRegionVolumeResidualCubicMeters
                  < 1.0e-10,
          "real multi-region wing retains junction faces and closes its coarse-grid volume ledger");
    check(pressureVolumes.regions.size() == result.scene.regions.size()
              && !pressureVolumes.controlVolumes.empty()
              && openingPatches.patches.size()
                  >= openingQuadrature.points.size()
              && !openingFaceCrossings.crossings.empty()
              && openingFaceCrossings.candidateSegmentCount
                  == 2 * openingFaceCrossings.crossings.size()
              && openingFaceCrossings.unpairedContactSegmentCount == 0
              && openingFaceCrossings.crossingLengthMeters > 0.0
              && cappedFacePartitions.touchedFaceCount == 9
              && cappedFacePartitions.faces.size() == 9
              && cappedFacePartitions.partitions.size() == 9
              && cappedFacePartitions.unresolvedTouchedFaceCount == 0
              && unresolvedCappedFaceRecordCount == 0
              && unpairedMaterialFaceCount == 0
              && cappedFailureSourceIds.empty()
              && cappedFacePartitions.segmentPairTestCount > 0
              && std::abs(openingPatches.totalAreaSquareMeters
                          - openingCaps.totalAreaSquareMeters)
                  < 1.0e-10
              && pressureFaceLinks.unresolvedActiveFaceCount == 0
              && pressureFaceLinks.unresolvedCappedFaceCount == 0
              && pressureFaceLinks.resolvedPartitionFaceCount == 10
              && pressureFaceLinks.unresolvedEmbeddedOpeningPatchCount == 24
              && pressureFaceLinks.embeddedOpeningRejections.size() == 24
              && pressureFaceLinks.embeddedOpeningBothSideOneRingCount == 24
              && pressureFaceLinks.embeddedOpeningSingleSideOneRingCount == 0
              && pressureFaceLinks.embeddedOpeningNoSideOneRingCount == 0
              && !pressureFaceLinks.embeddedOpeningOneRingSupports.empty()
              && rejectedEmbeddedOpeningIds.size() == 2
              && embeddedRejectionsAreNonAdmissible
              && embeddedOneRingSupportIsTwoSided
              && pressureFaceLinks.unresolvedEmbeddedOpeningAreaSquareMeters
                  > 0.0
              && mimeticAuditMatches,
          "real wing closes cap-crossed pressure faces and retains explicit embedded-opening limits");
    const simwing::viewer::StructureFrameMapping mapping =
        simwing::viewer::makeStructureFrameMapping(
            result.scene, assembly, structure);
    simwing::fsi::StructureStepSettings step;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    step.constraintIterations =
        assembly.settings.suspensionSolverIterations;
    const auto restAudit = simwing::fsi::auditStructureRestState(
        structure, step);
    check(restAudit.stationary
              && restAudit.maximumInitialSuspensionExtensionMeters == 0.0,
          "real structural export is stationary under a transactional zero-load probe");
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
    check(firstFrame.lines.size()
                  == assembly.definition.constraints.size() + 190
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
              && decodedFrame.lines.size()
                  == assembly.definition.constraints.size() + 190
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
