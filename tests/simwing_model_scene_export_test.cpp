#include "../src/model/nurbs_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string &message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

std::size_t sourceIndex(int rib, int point, int field)
{
    return static_cast<std::size_t>(
        rib + (point + field * 500) * 101 - 50601);
}

std::size_t shapingIndex(int panel, int point)
{
    return static_cast<std::size_t>(panel + point * 101 - 101);
}

std::size_t tessellationIndex(int panel, int point, int segment)
{
    return static_cast<std::size_t>(
        (panel + (point + segment * 500) * 101) * 3 - 151803);
}

std::size_t holeIndex(int rib, int hole, int field)
{
    return static_cast<std::size_t>(
        rib + (hole + field * 200) * 101 - 20301);
}

struct CaptureArrays
{
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> shaping;
    std::vector<double> tessellation;
    std::vector<double> holes;
};

CaptureArrays makeCapture()
{
    constexpr int pointCount = 8;
    const std::array<double, pointCount> chordCm{
        100.0, 70.0, 30.0, 0.0, 5.0, 20.0, 70.0, 100.0};
    const std::array<double, pointCount> heightCm{
        0.0, 15.0, 20.0, 0.0, -3.0, -10.0, -10.0, 0.0};

    CaptureArrays arrays;
    const std::size_t sourceSize = sourceIndex(1, pointCount, 49) + 1;
    arrays.u.assign(sourceSize, 0.0);
    arrays.v.assign(sourceSize, 0.0);
    arrays.w.assign(sourceSize, 0.0);
    arrays.shaping.assign(shapingIndex(1, pointCount) + 1, 0.0);
    arrays.tessellation.assign(
        tessellationIndex(1, pointCount, 2) + 3, 0.0);
    arrays.holes.assign(holeIndex(1, 1, 9) + 1, 0.0);

    for (int rib = 0; rib <= 1; ++rib) {
        const double spanCm = rib == 0 ? 0.0 : 100.0;
        for (int point = 1; point <= pointCount; ++point) {
            const std::size_t planar = sourceIndex(rib, point, 3);
            arrays.u[planar] = chordCm[point - 1];
            arrays.v[planar] = heightCm[point - 1];

            const std::size_t spatial = sourceIndex(rib, point, 47);
            arrays.u[spatial] = spanCm;
            arrays.v[spatial] = chordCm[point - 1];
            arrays.w[spatial] = -heightCm[point - 1];
        }
    }
    for (int point = 1; point <= pointCount; ++point) {
        for (int segment = 1; segment <= 2; ++segment) {
            const int rib = segment == 1 ? 1 : 0;
            const std::size_t source = sourceIndex(rib, point, 47);
            const std::size_t tessellation =
                tessellationIndex(1, point, segment);
            arrays.tessellation[tessellation] = arrays.u[source];
            arrays.tessellation[tessellation + 1] = arrays.v[source];
            arrays.tessellation[tessellation + 2] = arrays.w[source];
        }
    }

    arrays.holes[0] = 1.0;
    arrays.holes[holeIndex(0, 1, 2)] = 50.0;
    arrays.holes[holeIndex(0, 1, 3)] = 0.0;
    arrays.holes[holeIndex(0, 1, 4)] = 10.0;
    arrays.holes[holeIndex(0, 1, 5)] = 5.0;
    arrays.holes[holeIndex(0, 1, 6)] = 0.0;
    arrays.holes[holeIndex(0, 1, 7)] = 0.0;
    arrays.holes[holeIndex(0, 1, 9)] = 1.0;
    return arrays;
}

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

lep::SimWingSceneExportSettings settings()
{
    lep::SimWingSceneExportSettings value;
    value.designChecksum = "fixture-checksum";
    value.exporterVersion = "fixture-exporter";
    value.fabricMaterials = {
        fabric(simwing::fsi::SurfaceRole::Skin, "skin"),
        fabric(simwing::fsi::SurfaceRole::Rib, "rib")};
    value.lineMaterials.push_back(
        {"Dyneema", "Dyneema fixture", 0.001, 0.002, 20000.0, 1.0});
    value.pilot.name = "fixture pilot";
    value.pilot.massKg = 80.0;
    value.pilot.centerOfMassPositionMeters = {0.0, 0.0, -1.0};
    value.pilot.bodyToWorld = {1.0, 0.0, 0.0, 0.0};
    value.pilot.principalInertiaKgSquareMeters = {8.0, 9.0, 4.0};
    value.pilot.harnessPoints.push_back(
        {{0.0, 0.0, -1.0}, {0.0, 0.0, 0.0}});
    value.pilot.endpointMatchToleranceMeters = 1.0e-4;
    value.suspensionJunctionMassKg = 0.02;
    value.surfaceEndpointMatchToleranceMeters = 1.0e-4;
    return value;
}

void captureFixture()
{
    constexpr int pointCount = 8;
    CaptureArrays arrays = makeCapture();
    lep::resetNurbsModel();
    lep_nurbs_capture_panel(arrays.u.data(),
                            arrays.v.data(),
                            arrays.w.data(),
                            arrays.shaping.data(),
                            arrays.tessellation.data(),
                            1,
                            pointCount,
                            4,
                            2,
                            1,
                            0,
                            0);
    for (int rib = 0; rib <= 1; ++rib) {
        lep_nurbs_capture_rib(arrays.u.data(),
                              arrays.v.data(),
                              arrays.w.data(),
                              arrays.holes.data(),
                              100.0,
                              rib,
                              pointCount);
    }

    lep_nurbs_set_line_capture(1);
    const std::string type = "Dyneema";
    lep_nurbs_set_line_tag("A1",
                           2,
                           1,
                           0,
                           type.data(),
                           static_cast<int>(type.size()),
                           1.0);
    // A surface vertex -> explicit line junction -> pilot harness graph.
    lep_nurbs_capture_line(100.0, 70.0, -15.0,
                           50.0, 50.0, 50.0, 1);
    lep_nurbs_capture_line(50.0, 50.0, 50.0,
                           0.0, 0.0, 100.0, 1);
    lep_nurbs_set_line_capture(0);
}

void testAuthoritativeSceneExport()
{
    captureFixture();
    const lep::SimWingSceneExportResult first =
        lep::buildSimWingScene(settings());
    if (!first.success) {
        for (const std::string &error : first.errors) {
            std::cerr << "export: " << error << '\n';
        }
    }
    check(first.success, "captured analytical model exports scene-v2");
    if (!first.success) {
        return;
    }
    check(first.validation.ok(), "exported scene validates");
    check(first.scene.metadata.schemaMajor == simwing::fsi::sceneSchemaMajor
              && first.scene.metadata.schemaMinor
                     == simwing::fsi::sceneSchemaMinor,
          "export uses the current scene-v2 schema");
    check(std::abs(first.scene.metadata.sourceLengthToMeters - 0.01)
              < 1.0e-15,
          "metadata records centimetre source units");

    const auto roleCount = [&first](simwing::fsi::SurfaceRole role) {
        return std::count_if(
            first.scene.triangles.begin(), first.scene.triangles.end(),
            [role](const simwing::fsi::Triangle &triangle) {
                return triangle.role == role;
            });
    };
    check(roleCount(simwing::fsi::SurfaceRole::Skin) > 0,
          "analytical skins are triangulated");
    check(roleCount(simwing::fsi::SurfaceRole::Rib) > 0,
          "exact rib faces are triangulated");
    check(std::count_if(
              first.scene.openings.begin(), first.scene.openings.end(),
              [](const simwing::fsi::Opening &opening) {
                  return opening.role == simwing::fsi::OpeningRole::Intake;
              }) == 2,
          "the authored open intake is preserved on both halves");
    check(std::count_if(
              first.scene.openings.begin(), first.scene.openings.end(),
              [](const simwing::fsi::Opening &opening) {
                  return opening.role == simwing::fsi::OpeningRole::Crossport;
              }) == 1,
          "the exact centre-rib hole becomes one crossport");

    double maximumCoordinate = 0.0;
    for (const simwing::fsi::Vertex &vertex : first.scene.vertices) {
        maximumCoordinate = std::max(
            {maximumCoordinate,
             std::abs(vertex.positionMeters.x),
             std::abs(vertex.positionMeters.y),
             std::abs(vertex.positionMeters.z)});
    }
    check(maximumCoordinate > 0.99 && maximumCoordinate < 1.01,
          "geometry crosses the capture boundary in SI metres");
    check(first.scene.suspensionLines.size() == 2,
          "line graph segments are preserved");
    check(first.scene.suspensionJunctions.size() == 1,
          "unmatched shared endpoint is an explicit junction");
    check(first.scene.attachments.size() == 3,
          "surface, junction and pilot endpoints are distinct attachments");

    const lep::SimWingSceneExportResult second =
        lep::buildSimWingScene(settings());
    std::ostringstream firstBytes(std::ios::binary);
    std::ostringstream secondBytes(std::ios::binary);
    std::string error;
    check(second.success, "repeat export succeeds");
    check(simwing::fsi::writeScene(first.scene, firstBytes, &error),
          "first scene serializes");
    check(simwing::fsi::writeScene(second.scene, secondBytes, &error),
          "second scene serializes");
    check(firstBytes.str() == secondBytes.str(),
          "stable IDs and binary scene bytes are deterministic");
}

void testMissingPhysicalAssignmentFailsTransactionally()
{
    captureFixture();
    lep::SimWingSceneExportSettings incomplete = settings();
    incomplete.lineMaterials.clear();
    const lep::SimWingSceneExportResult result =
        lep::buildSimWingScene(incomplete);
    check(!result.success, "missing captured line material is rejected");
    check(result.scene.vertices.empty() && result.scene.triangles.empty(),
          "failed export returns no partial scene");
}

} // namespace

int main()
{
    testAuthoritativeSceneExport();
    testMissingPhysicalAssignmentFailsTransactionally();
    lep::resetNurbsModel();
    if (failures == 0) {
        std::cout << "simwing model scene export tests passed\n";
    }
    return failures == 0 ? 0 : 1;
}
