#include "scene_fluid_pressure_operator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

void addTetra(Scene& scene,
              const StableId vertexBase,
              const StableId triangleBase,
              const StableId insideRegion,
              const StableId outsideRegion,
              const double leftX,
              const double rightX,
              const double yRadius,
              const double zRadius,
              const StableId sheet,
              const bool openMouth = false) {
    scene.vertices.insert(scene.vertices.end(), {
        {vertexBase + 0, {leftX, 1.5, 1.45}},
        {vertexBase + 1, {rightX, 1.5 - yRadius, 1.45 - zRadius}},
        {vertexBase + 2, {rightX, 1.5 + yRadius, 1.45 - zRadius}},
        {vertexBase + 3, {rightX, 1.5, 1.45 + zRadius}},
    });
    const std::array<Vec2, 3> chart{{{0.0, 0.0},
                                      {1.0, 0.0},
                                      {0.0, 1.0}}};
    scene.triangles.insert(scene.triangles.end(), {
        {triangleBase + 0,
         {vertexBase + 0, vertexBase + 2, vertexBase + 1},
         chart, insideRegion, outsideRegion, 100, sheet, SurfaceRole::Skin},
        {triangleBase + 1,
         {vertexBase + 0, vertexBase + 1, vertexBase + 3},
         chart, insideRegion, outsideRegion, 100, sheet, SurfaceRole::Skin},
        {triangleBase + 2,
         {vertexBase + 0, vertexBase + 3, vertexBase + 2},
         chart, insideRegion, outsideRegion, 100, sheet, SurfaceRole::Skin},
    });
    if (openMouth) {
        scene.openings.push_back({
            700, {vertexBase + 1, vertexBase + 2, vertexBase + 3},
            insideRegion, outsideRegion, OpeningRole::Intake,
        });
    } else {
        scene.triangles.push_back({
            triangleBase + 3,
            {vertexBase + 1, vertexBase + 2, vertexBase + 3},
            chart, insideRegion, outsideRegion, 100, sheet,
            SurfaceRole::Skin,
        });
    }
}

Scene nestedScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-solve";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-solve-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "outer"},
        {3, RegionKind::Cell, "inner"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, 2.8, 0.3, 0.3, 900);
    addTetra(scene, 20, 600, 3, 2, 1.6, 2.4, 0.1, 0.1, 901);
    return scene;
}

Scene openScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-solve-open";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-solve-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, 2.0, 0.3, 0.3, 900, true);
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;
    SceneFluidSurfaceState state;
    SceneFluidGridEpoch epoch;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet openingQuadrature;
    SceneFluidOpeningGridPatchSet openingPatches;
    SceneFluidOpeningFaceCrossingSet openingFaceCrossings;
    SceneFluidCappedFacePartitionSet cappedFacePartitions;
    SceneFluidCellVolumeSet volumes;
    SceneFluidRegionConnectivity connectivity;
    SceneFluidPressureControlVolumeSet pressureVolumes;
    SceneFluidPressureFaceLinkSet faceLinks;
    SceneFluidPressureOperator pressureOperator;

    explicit Fixture(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(surface.definition,
                   structureAssembly.mappings,
                   structure),
          state(captureSceneFluidSurfaceState(
              surface.definition, structureAssembly.mappings, structure)),
          epoch(buildSceneFluidGridEpoch(
              surface.definition, state, grid(), transfer)),
          caps(buildSceneFluidOpeningCaps(surface.definition, state)),
          openingQuadrature(buildSceneFluidOpeningQuadrature(
              surface.definition, state, caps)),
          openingPatches(buildSceneFluidOpeningGridPatches(
              surface.definition, state, caps, openingQuadrature, grid())),
          openingFaceCrossings(buildSceneFluidOpeningFaceCrossings(
              surface.definition, state, caps, openingQuadrature,
              openingPatches, grid())),
          cappedFacePartitions(buildSceneFluidCappedFacePartitions(
              surface.definition, state, grid(), transfer, epoch, caps,
              openingQuadrature, openingPatches, openingFaceCrossings)),
          volumes(buildSceneFluidCellVolumes(
              surface.definition, state, grid(), transfer, epoch)),
          connectivity(buildSceneFluidRegionConnectivity(
              surface.definition)),
          pressureVolumes(buildSceneFluidPressureControlVolumes(
              surface.definition, volumes, connectivity)),
          faceLinks(buildSceneFluidPressureFaceLinks(
              surface.definition, state, grid(), transfer, epoch, caps,
              openingQuadrature, openingPatches, openingFaceCrossings,
              cappedFacePartitions, volumes, connectivity,
              pressureVolumes)),
          pressureOperator(buildSceneFluidPressureOperator(
              surface.definition, state, grid(), transfer, epoch, caps,
              openingQuadrature, openingPatches, volumes, connectivity,
              pressureVolumes, faceLinks)) {}
};

SceneFluidPressureSolveSettings strictSettings() {
    SceneFluidPressureSolveSettings settings;
    settings.absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-11;
    settings.maximumIterations = 4000;
    return settings;
}

void normalizeGauges(const SceneFluidPressureOperator& pressureOperator,
                     std::vector<double>& pressure) {
    for (const auto& component : pressureOperator.components) {
        const double gauge = pressure[component.gaugeControlVolumeIndex];
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            pressure[pressureOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]] -= gauge;
        }
        pressure[component.gaugeControlVolumeIndex] = 0.0;
    }
}

double maximumError(const std::vector<double>& first,
                    const std::vector<double>& second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result = std::max(result, std::abs(first[index] - second[index]));
    }
    return result;
}

void testManufacturedMultiComponentRecovery() {
    Fixture fixture(nestedScene());
    const auto& pressureOperator = fixture.pressureOperator;
    std::vector<double> expected(pressureOperator.rows.size(), 0.0);
    for (std::size_t row = 0; row < expected.size(); ++row) {
        const double sample = static_cast<double>(row + 1);
        expected[row] = 0.7 * std::sin(0.17 * sample)
            + 0.2 * std::cos(0.31 * sample)
            + 0.03 * sample;
    }
    normalizeGauges(pressureOperator, expected);
    const auto rightHandSide = applySceneFluidPressureOperator(
        pressureOperator, expected);
    std::vector<double> firstPressure(expected.size(), 0.0);
    std::vector<double> secondPressure(expected.size(), 0.0);
    const auto settings = strictSettings();
    const auto first = solveSceneFluidPressureSystem(
        pressureOperator, rightHandSide, firstPressure, settings);
    const auto second = solveSceneFluidPressureSystem(
        pressureOperator, rightHandSide, secondPressure, settings);
    check(first == second
              && firstPressure == secondPressure
              && first.compatible && first.converged && first.finite
              && first.iterationCount > 0
              && first.componentCount == 3,
          "manufactured pressure solve converges deterministically in all components");
    check(maximumError(firstPressure, expected) < 2.0e-10,
          "manufactured pressure solve recovers the gauge-normalized field");
    check(first.finalResidualL2PascalsMeters < 1.0e-11
              && first.finalResidualMaximumPascalsMeters < 1.0e-10,
          "manufactured pressure solve closes its explicitly recomputed residual");
    for (const auto& component : first.components) {
        checkNear(component.rightHandSideSumPascalsMeters,
                  0.0, 2.0e-12,
                  "manufactured component RHS satisfies compatibility");
        check(component.pressureGaugeAfterPascals == 0.0
                  && firstPressure[
                      component.gaugeControlVolumeIndex] == 0.0,
              "manufactured pressure solve commits exact component gauges");
    }
}

void testFaceAlignedIntakeRecovery() {
    Fixture fixture(openScene());
    const auto& pressureOperator = fixture.pressureOperator;
    std::vector<double> expected(pressureOperator.rows.size(), 0.0);
    for (const auto& control : fixture.pressureVolumes.controlVolumes) {
        expected[control.controlVolumeIndex] =
            control.regionId == 2 ? 125.0 : -25.0;
    }
    normalizeGauges(pressureOperator, expected);
    const auto rightHandSide = applySceneFluidPressureOperator(
        pressureOperator, expected);
    std::vector<double> pressure(expected.size(), 7.0);
    const auto diagnostics = solveSceneFluidPressureSystem(
        pressureOperator, rightHandSide, pressure, strictSettings());
    check(diagnostics.compatible && diagnostics.converged
              && diagnostics.finite
              && diagnostics.componentCount == 1,
          "face-aligned intake pressure system solves as one component");
    check(maximumError(pressure, expected) < 2.0e-10,
          "face-aligned intake solve recovers its cross-region pressure jump");
}

void testCompatibilityRollbackAndValidation() {
    Fixture fixture(nestedScene());
    const auto& pressureOperator = fixture.pressureOperator;
    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t row = 0; row < manufactured.size(); ++row) {
        manufactured[row] = std::sin(0.23 * static_cast<double>(row + 1));
    }
    const auto compatible = applySceneFluidPressureOperator(
        pressureOperator, manufactured);

    auto incompatible = compatible;
    incompatible.front() += 1.0e-4;
    std::vector<double> pressure(pressureOperator.rows.size(), 3.5);
    const auto originalPressure = pressure;
    const auto incompatibleDiagnostics = solveSceneFluidPressureSystem(
        pressureOperator, incompatible, pressure, strictSettings());
    check(!incompatibleDiagnostics.compatible
              && !incompatibleDiagnostics.converged
              && incompatibleDiagnostics.finite
              && incompatibleDiagnostics
                    .maximumAbsoluteComponentCompatibilityPascalsMeters
                    > 9.0e-5
              && pressure == originalPressure,
          "incompatible component source is reported without pressure mutation");

    auto truncatedSettings = strictSettings();
    truncatedSettings.absoluteResidualTolerancePascalsMeters = 1.0e-30;
    truncatedSettings.relativeResidualTolerance = 0.0;
    truncatedSettings.maximumIterations = 1;
    const auto truncated = solveSceneFluidPressureSystem(
        pressureOperator, compatible, pressure, truncatedSettings);
    check(truncated.compatible && !truncated.converged
              && truncated.finite && truncated.iterationCount == 1
              && pressure == originalPressure,
          "truncated pressure solve rolls its warm start back exactly");

    auto invalidSettings = strictSettings();
    invalidSettings.absoluteResidualTolerancePascalsMeters = -1.0;
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidPressureSystem(
            pressureOperator, compatible, pressure, invalidSettings)); },
        "pressure solve rejects invalid residual settings");
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidPressureSystem(
            pressureOperator,
            std::span<const double>(compatible).first(
                compatible.size() - 1),
            pressure, strictSettings())); },
        "pressure solve rejects a short right-hand side");
    auto nonFinitePressure = originalPressure;
    nonFinitePressure.front() =
        std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidPressureSystem(
            pressureOperator, compatible, nonFinitePressure,
            strictSettings())); },
        "pressure solve rejects a non-finite warm start");
}

} // namespace

int main() {
    try {
        testManufacturedMultiComponentRecovery();
        testFaceAlignedIntakeRecovery();
        testCompatibilityRollbackAndValidation();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid pressure-solve check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-solve checks passed");
    return 0;
}
