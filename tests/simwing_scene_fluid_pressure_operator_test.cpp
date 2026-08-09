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

template<typename Callback>
void expectLimited(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::length_error&) {
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
        "sha256:scene-fluid-pressure-operator";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-operator-test/1";
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

Scene openScene(const double mouthX = 2.0) {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-operator-open";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-operator-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, mouthX, 0.3, 0.3, 900, true);
    return scene;
}

Scene tiltedOpenScene() {
    auto scene = openScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-operator-tilted-open";
    scene.vertices[1].positionMeters.x = 2.6;
    scene.vertices[2].positionMeters.x = 2.8;
    scene.vertices[3].positionMeters.x = 2.7;
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
    SceneFluidCellVolumeSet volumes;
    SceneFluidRegionConnectivity connectivity;
    SceneFluidPressureControlVolumeSet pressureVolumes;
    SceneFluidPressureFaceLinkSet faceLinks;

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
          volumes(buildSceneFluidCellVolumes(
              surface.definition, state, grid(), transfer, epoch)),
          connectivity(buildSceneFluidRegionConnectivity(
              surface.definition)),
          pressureVolumes(buildSceneFluidPressureControlVolumes(
              surface.definition, volumes, connectivity)),
          faceLinks(buildSceneFluidPressureFaceLinks(
              surface.definition, state, grid(), transfer, epoch, caps,
              openingQuadrature, openingPatches, volumes, connectivity,
              pressureVolumes)) {}

    SceneFluidPressureOperator pressureOperator(
        const SceneFluidPressureOperatorLimits& limits = {}) const {
        return buildSceneFluidPressureOperator(
            surface.definition, state, grid(), transfer, epoch, caps,
            openingQuadrature, openingPatches, volumes, connectivity,
            pressureVolumes, faceLinks, limits);
    }

    void validate(const SceneFluidPressureOperator& value) const {
        validateSceneFluidPressureOperator(
            value, surface.definition, state, grid(), transfer, epoch,
            caps, openingQuadrature, openingPatches, volumes,
            connectivity, pressureVolumes, faceLinks);
    }
};

double dot(const std::vector<double>& first,
           const std::vector<double>& second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result += first[index] * second[index];
    }
    return result;
}

void testSymmetricConservativeOperator() {
    Fixture fixture(nestedScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "pressure-operator nested fixture assembles");
    const auto first = fixture.pressureOperator();
    const auto repeated = fixture.pressureOperator();
    check(first == repeated
              && first.version == sceneFluidPressureOperatorVersion
              && first.fingerprint != 0
              && first.rows.size()
                  == fixture.pressureVolumes.controlVolumes.size()
              && first.entries.size() == 2 * fixture.faceLinks.links.size()
              && first.components.size()
                  == fixture.pressureVolumes.components.size(),
          "pressure operator is deterministic and retains its exact graph");
    std::vector<double> firstPressure(first.rows.size());
    std::vector<double> secondPressure(first.rows.size());
    for (std::size_t row = 0; row < first.rows.size(); ++row) {
        const double sample = static_cast<double>(row + 1);
        firstPressure[row] = std::sin(0.37 * sample) + 0.1 * sample;
        secondPressure[row] = std::cos(0.19 * sample) - 0.03 * sample;
    }
    const auto firstApplied = applySceneFluidPressureOperator(
        first, firstPressure);
    const auto secondApplied = applySceneFluidPressureOperator(
        first, secondPressure);
    checkNear(dot(firstPressure, secondApplied),
              dot(secondPressure, firstApplied), 2.0e-12,
              "pressure operator is symmetric");

    double expectedEnergy = 0.0;
    for (const auto& link : fixture.faceLinks.links) {
        const double difference =
            firstPressure[link.minusControlVolumeIndex]
            - firstPressure[link.plusControlVolumeIndex];
        expectedEnergy += link.geometryWeightMeters
            * difference * difference;
    }
    checkNear(dot(firstPressure, firstApplied), expectedEnergy, 2.0e-12,
              "pressure operator energy equals its conservative link sum");
    check(expectedEnergy > 0.0,
          "nonconstant pressure has positive graph energy");

    for (const auto& component : first.components) {
        double componentSum = 0.0;
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            componentSum += firstApplied[
                first.componentControlVolumeIndices[
                    component.firstControlVolumeMember + offset]];
        }
        checkNear(componentSum, 0.0, 2.0e-12,
                  "each pressure component conserves its integrated row sum");
        check(first.rows[component.gaugeControlVolumeIndex].isGauge
                  && first.rows[
                      component.gaugeControlVolumeIndex].componentIndex
                      == component.componentIndex,
              "each operator component retains its canonical gauge owner");
    }

    std::vector<double> constant(first.rows.size(), 1.0);
    const auto nullResult = applySceneFluidPressureOperator(first, constant);
    check(std::ranges::all_of(nullResult, [](const double value) {
              return value == 0.0;
          }),
          "constant pressure lies exactly in every component nullspace");
    fixture.validate(first);
}

void testOpeningTopologyBoundaries() {
    Fixture fixture(openScene());
    const auto faceAligned = fixture.pressureOperator();
    check(fixture.faceLinks.resolvedOpeningFaceCount == 1
              && fixture.faceLinks.unresolvedOpeningFaceCount == 0
              && faceAligned.components.size() == 1,
          "face-aligned intake closes one connected pressure component");
    const auto constantResult = applySceneFluidPressureOperator(
        faceAligned, std::vector<double>(faceAligned.rows.size(), 1.0));
    check(std::ranges::all_of(constantResult, [](const double value) {
              return value == 0.0;
          }),
          "connected intake operator retains its exact constant null mode");
    std::vector<double> regionJump(faceAligned.rows.size(), 0.0);
    for (const auto& control : fixture.pressureVolumes.controlVolumes) {
        regionJump[control.controlVolumeIndex] =
            control.regionId == 2 ? 1.0 : 0.0;
    }
    const auto jumpResult = applySceneFluidPressureOperator(
        faceAligned, regionJump);
    checkNear(dot(regionJump, jumpResult), 0.18, 2.0e-15,
              "face-aligned intake contributes its exact aperture graph energy");

    Fixture offFace(tiltedOpenScene());
    const std::size_t unresolvedFaceCount =
        offFace.faceLinks.unresolvedActiveFaceCount
        + offFace.faceLinks.unresolvedAmbiguousFaceCount
        + offFace.faceLinks.unresolvedOpeningFaceCount;
    const auto embeddedOperator = offFace.pressureOperator();
    const auto embeddedLink = std::ranges::find(
        offFace.faceLinks.links,
        SceneFluidPressureLinkGeometryKind::EmbeddedOpening,
        &SceneFluidPressureFaceLink::geometryKind);
    check(unresolvedFaceCount == 0
              && offFace.pressureVolumes.components.size() == 1
              && offFace.faceLinks.embeddedOpeningLinkCount == 1
              && embeddedLink != offFace.faceLinks.links.end()
              && embeddedOperator.components.size() == 1,
          "tilted off-face intake closes its authored pressure component");
    const auto embeddedConstant = applySceneFluidPressureOperator(
        embeddedOperator,
        std::vector<double>(embeddedOperator.rows.size(), 1.0));
    check(std::ranges::all_of(
              embeddedConstant,
              [](const double value) { return value == 0.0; }),
          "embedded intake retains the exact constant pressure null mode");
    std::vector<double> embeddedJump(embeddedOperator.rows.size(), 0.0);
    for (const auto& control : offFace.pressureVolumes.controlVolumes) {
        embeddedJump[control.controlVolumeIndex] =
            control.regionId == 2 ? 1.0 : 0.0;
    }
    const auto embeddedResult = applySceneFluidPressureOperator(
        embeddedOperator, embeddedJump);
    if (embeddedLink != offFace.faceLinks.links.end()) {
        checkNear(dot(embeddedJump, embeddedResult),
                  embeddedLink->geometryWeightMeters, 2.0e-15,
                  "embedded intake contributes its exact centroid-weighted graph energy");
    }
}

void testCorruptionInputsAndLimits() {
    Fixture fixture(nestedScene());
    const auto accepted = fixture.pressureOperator();
    auto corrupt = accepted;
    corrupt.entries.front().geometryWeightMeters += 0.01;
    expectInvalid(
        [&] { fixture.validate(corrupt); },
        "pressure-operator validation rejects entry corruption");
    expectInvalid(
        [&] { static_cast<void>(applySceneFluidPressureOperator(
            corrupt, std::vector<double>(corrupt.rows.size()))); },
        "pressure-operator application rejects corrupt topology");
    expectInvalid(
        [&] { static_cast<void>(applySceneFluidPressureOperator(
            accepted, std::vector<double>(accepted.rows.size() - 1))); },
        "pressure-operator application rejects a short pressure field");
    auto nonFinite = std::vector<double>(accepted.rows.size(), 0.0);
    nonFinite.front() = std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(applySceneFluidPressureOperator(
            accepted, nonFinite)); },
        "pressure-operator application rejects non-finite pressure");

    SceneFluidPressureOperatorLimits limits;
    limits.maximumRows = accepted.rows.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.pressureOperator(limits)); },
        "pressure-operator assembly bounds row count");
    limits = {};
    limits.maximumEntries = accepted.entries.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.pressureOperator(limits)); },
        "pressure-operator assembly bounds entry count");
    limits = {};
    limits.maximumComponents = accepted.components.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.pressureOperator(limits)); },
        "pressure-operator assembly bounds component count");
    limits = {};
    limits.maximumOperatorBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.pressureOperator(limits)); },
        "pressure-operator assembly bounds owned storage");
}

} // namespace

int main() {
    try {
        testSymmetricConservativeOperator();
        testOpeningTopologyBoundaries();
        testCorruptionInputsAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid pressure-operator check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-operator checks passed");
    return 0;
}
