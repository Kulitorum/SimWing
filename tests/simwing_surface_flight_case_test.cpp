#include "surface_flight_case.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <string_view>

namespace {

using namespace simwing::fsi;
namespace viewer = simwing::viewer;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

Scene flightScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:surface-flight-case";
    scene.metadata.exporterVersion = "surface-flight-case-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {0.0, 0.0, 0.0}},
        {11, {1.0, 0.0, 0.0}},
        {12, {0.0, 1.0, 0.0}},
        {13, {0.0, 0.0, 1.0}},
    };
    scene.fabricMaterials = {
        {100, "test-fabric", 800.0, 600.0, 200.0, 0.01,
         0.04, 0.01, 0.0, 0.0},
    };
    const std::array<Vec2, 3> chart{
        Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{0.0, 1.0}};
    scene.triangles = {
        {500, {11, 12, 13}, chart, 2, 1, 100, 900, SurfaceRole::Skin},
        {501, {10, 13, 12}, chart, 2, 1, 100, 901, SurfaceRole::Skin},
        {502, {10, 11, 13}, chart, 2, 1, 100, 902, SurfaceRole::Skin},
        {503, {10, 12, 11}, chart, 2, 1, 100, 903, SurfaceRole::Skin},
    };
    return scene;
}

bool hasScalar(const viewer::DiagnosticFrame& frame,
               const std::string_view name,
               const viewer::FieldAssociation association) {
    return std::ranges::any_of(frame.scalarFields, [&](const auto& field) {
        return field.name == name && field.association == association;
    });
}

bool hasVector(const viewer::DiagnosticFrame& frame,
               const std::string_view name,
               const viewer::FieldAssociation association) {
    return std::ranges::any_of(frame.vectorFields, [&](const auto& field) {
        return field.name == name && field.association == association;
    });
}

void testTransactionalSurfaceFlightFrame() {
    SurfaceFlightCaseSettings settings;
    settings.structureStep.gravityMetersPerSecondSquared = {};
    settings.structureStep.workerThreads = 0;
    settings.aerodynamics.windRampSeconds = 0.0;
    SurfaceFlightCase first(flightScene(), settings);
    SurfaceFlightCase second(flightScene(), settings);

    const viewer::DiagnosticFrame firstFrame = first.advance();
    const viewer::DiagnosticFrame replayFrame = second.advance();
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> replayBytes;
    const bool serialized = viewer::serializeFrame(firstFrame, firstBytes)
        && viewer::serializeFrame(replayFrame, replayBytes);
    check(firstFrame.step == 1 && firstFrame.simulationTimeSeconds > 0.0
              && first.acceptedStepCount() == 1,
          "surface flight: one accepted transaction advances gas and structure");
    check(serialized && firstBytes == replayBytes,
          "surface flight: equal scenes and settings replay deterministically");
    check(hasScalar(firstFrame, "surface_aero.pressure_jump",
                    viewer::FieldAssociation::Triangle)
              && hasScalar(firstFrame, "surface_aero.external_traction",
                           viewer::FieldAssociation::Triangle)
              && hasVector(firstFrame, "surface_aero.applied_force",
                           viewer::FieldAssociation::Vertex),
          "surface flight: live pressure, traction, and force fields reach viewer");
    check(first.diagnostics().finite
              && first.diagnostics().transfer.forceResidualNormNewtons < 1.0e-10
              && first.diagnostics().aerodynamics.liftNewtons > 0.0,
          "surface flight: accepted diagnostic ledgers are finite and conservative");
}

} // namespace

int main() {
    testTransactionalSurfaceFlightFrame();
    if (failures != 0) {
        std::fprintf(stderr, "%d surface-flight test(s) failed\n", failures);
        return 1;
    }
    std::puts("surface-flight tests passed");
    return 0;
}
