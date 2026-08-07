#include "engine_paths.h"
#include "input_migration.h"
#include "nurbs_model.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
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
    value.surfaceEndpointMatchToleranceMeters = 0.002;
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
    check(result.scene.suspensionLines.size() == 190,
          "real export preserves every captured suspension segment");
    check(attachmentCount(simwing::fsi::AttachmentKind::PilotHarness) == 2,
          "real export connects both captured riser roots to the pilot");
    check(attachmentCount(simwing::fsi::AttachmentKind::SurfaceVertex) > 0,
          "real export connects terminal lines to authoritative surface vertices");
    check(!result.scene.suspensionJunctions.empty(),
          "real export preserves intermediate suspension junctions");
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
