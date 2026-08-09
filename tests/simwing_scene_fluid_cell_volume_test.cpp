#include "scene_fluid_cell_volume.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <utility>

namespace {

using namespace simwing::fsi;
using namespace simwing::fsi::fluid;

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

void addTetraPositions(Scene& scene,
                       const StableId vertexBase,
                       const StableId triangleBase,
                       const StableId insideRegion,
                       const StableId outsideRegion,
                       const StableId sheet,
                       const std::array<Vec3, 4>& positions) {
    scene.vertices.insert(scene.vertices.end(), {
        {vertexBase + 0, positions[0]},
        {vertexBase + 1, positions[1]},
        {vertexBase + 2, positions[2]},
        {vertexBase + 3, positions[3]},
    });
    const auto chart = [&](const std::size_t firstIndex,
                           const std::size_t secondIndex,
                           const std::size_t thirdIndex) {
        const Vec3 first = positions[firstIndex];
        const Vec3 second = positions[secondIndex];
        const Vec3 third = positions[thirdIndex];
        const Vec3 edge{second.x - first.x,
                        second.y - first.y,
                        second.z - first.z};
        const Vec3 diagonal{third.x - first.x,
                            third.y - first.y,
                            third.z - first.z};
        const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
        const double thirdX = (diagonal.x * edge.x
                               + diagonal.y * edge.y
                               + diagonal.z * edge.z) / edgeLength;
        const double diagonalLengthSquared =
            diagonal.x * diagonal.x
            + diagonal.y * diagonal.y
            + diagonal.z * diagonal.z;
        return std::array<Vec2, 3>{{
            {0.0, 0.0},
            {edgeLength, 0.0},
            {thirdX, std::sqrt(std::max(
                0.0, diagonalLengthSquared - thirdX * thirdX))},
        }};
    };
    scene.triangles.insert(scene.triangles.end(), {
        {triangleBase + 0,
         {vertexBase + 0, vertexBase + 2, vertexBase + 1},
         chart(0, 2, 1), insideRegion, outsideRegion, 100, sheet,
         SurfaceRole::Skin},
        {triangleBase + 1,
         {vertexBase + 0, vertexBase + 1, vertexBase + 3},
         chart(0, 1, 3), insideRegion, outsideRegion, 100, sheet,
         SurfaceRole::Skin},
        {triangleBase + 2,
         {vertexBase + 0, vertexBase + 3, vertexBase + 2},
         chart(0, 3, 2), insideRegion, outsideRegion, 100, sheet,
         SurfaceRole::Skin},
        {triangleBase + 3,
         {vertexBase + 1, vertexBase + 2, vertexBase + 3},
         chart(1, 2, 3), insideRegion, outsideRegion, 100, sheet,
         SurfaceRole::Skin},
    });
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
              const StableId sheet) {
    addTetraPositions(
        scene, vertexBase, triangleBase, insideRegion, outsideRegion, sheet,
        {{{leftX, 1.5, 1.45},
          {rightX, 1.5 - yRadius, 1.45 - zRadius},
          {rightX, 1.5 + yRadius, 1.45 - zRadius},
          {rightX, 1.5, 1.45 + zRadius}}});
}

Scene nestedScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-cell-volume";
    scene.metadata.exporterVersion = "scene-fluid-cell-volume-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "outer-cell"},
        {3, RegionKind::Cell, "inner-cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, 2.8, 0.3, 0.3, 900);
    addTetra(scene, 20, 600, 3, 2, 1.6, 2.4, 0.1, 0.1, 901);
    return scene;
}

Scene largeClosedScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-full-cell-volume";
    scene.metadata.exporterVersion = "scene-fluid-cell-volume-test/2";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "large-cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetraPositions(
        scene, 10, 500, 2, 1, 900,
        {{{0.11, 0.18, 0.27},
          {3.83, 0.31, 3.57},
          {3.69, 3.77, 0.22},
          {0.27, 3.63, 3.84}}});
    return scene;
}

Scene openTetraScene() {
    Scene scene = nestedScene();
    scene.metadata.designChecksum = "sha256:scene-fluid-open-tetra-volume";
    scene.metadata.exporterVersion = "scene-fluid-cell-volume-test/3";
    scene.regions.resize(2);
    scene.vertices.resize(4);
    scene.triangles.resize(3);
    scene.openings = {
        {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

Scene translatedLargeClosedScene() {
    Scene scene = largeClosedScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-translated-full-cell-volume";
    for (auto& vertex : scene.vertices) {
        vertex.positionMeters.x -= 2.0;
        vertex.positionMeters.y -= 3.0;
        vertex.positionMeters.z += 0.5;
    }
    return scene;
}

Scene openScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-open-volume";
    scene.metadata.exporterVersion = "scene-fluid-cell-volume-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.2, 1.5}},
        {11, {1.8, 1.2, 1.5}},
        {12, {1.8, 1.8, 1.5}},
        {13, {1.2, 1.8, 1.5}},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {0.6, 0.0}, {0.6, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 12, 13},
         {{{0.0, 0.0}, {0.6, 0.6}, {0.0, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    scene.openings = {
        {700, {10, 11, 12, 13}, 1, 2, OpeningRole::Intake},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

PeriodicCartesianGrid fineGrid() {
    return {{8, 8, 8}, {}, {4.0, 4.0, 4.0}};
}

PeriodicCartesianGrid translatedFineGrid() {
    return {{8, 8, 8}, {-2.0, -3.0, 0.5}, {2.0, 1.0, 4.5}};
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;

    explicit Fixture(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(surface.definition,
                   structureAssembly.mappings,
                   structure) {}
};

std::map<StableId, SceneFluidRegionVolume> regionVolumes(
    const SceneFluidCellVolumeSet& volumes) {
    std::map<StableId, SceneFluidRegionVolume> result;
    for (const auto& region : volumes.regionVolumes) {
        result.emplace(region.regionId, region);
    }
    return result;
}

void testNestedAnalyticVolumesAndAcceptedMotion() {
    Fixture fixture(nestedScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "cell-volume nested fixture assembles");
    const auto initialState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto initialEpoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, initialState, grid(), fixture.transfer);
    const auto first = buildSceneFluidCellVolumes(
        fixture.surface.definition, initialState, grid(), fixture.transfer,
        initialEpoch);
    const auto repeated = buildSceneFluidCellVolumes(
        fixture.surface.definition, initialState, grid(), fixture.transfer,
        initialEpoch);
    check(first == repeated
              && first.version == sceneFluidCellVolumeVersion
              && first.fingerprint != 0
              && first.cells.size() == grid().cellCount()
              && first.maximumTetrahedronVolumeResidualCubicMeters < 1.0e-12
              && first.maximumCellVolumeResidualCubicMeters < 1.0e-12
              && first.maximumRegionVolumeResidualCubicMeters < 1.0e-12,
          "nested cell volumes are deterministic and close every cell and region");

    const auto regions = regionVolumes(first);
    constexpr double outerVolume = 0.096;
    constexpr double innerVolume = 0.005333333333333333;
    checkNear(regions.at(1).summedCellVolumeCubicMeters,
              64.0 - outerVolume, 2.0e-12,
              "outside volume excludes the outer tetrahedron");
    checkNear(regions.at(2).summedCellVolumeCubicMeters,
              outerVolume - innerVolume, 2.0e-12,
              "nested annular region subtracts the inner tetrahedron");
    checkNear(regions.at(3).summedCellVolumeCubicMeters,
              innerVolume, 2.0e-12,
              "inner region recovers analytic tetrahedron volume");
    std::size_t mixedCellCount = 0;
    for (const auto& cell : first.cells) {
        mixedCellCount += cell.regionVolumeCount > 1;
        double fraction = 0.0;
        for (std::size_t offset = 0;
             offset < cell.regionVolumeCount; ++offset) {
            fraction += first.cellRegionVolumes[
                cell.firstRegionVolume + offset].volumeFraction;
        }
        checkNear(fraction, 1.0, 1.0e-10,
                  "each sparse cell fraction set sums to one");
    }
    check(mixedCellCount == 2,
          "the resolved nested tetrahedra cut exactly two Cartesian cells");
    validateSceneFluidCellVolumes(
        first, fixture.surface.definition, initialState, grid(),
        fixture.transfer, initialEpoch);

    for (std::size_t node = 0;
         node < fixture.structure.definition().nodes.size(); ++node) {
        fixture.structure.addExternalForce(
            node,
            {fixture.structure.definition().nodes[node].massKg,
             0.0, 0.0});
    }
    StructureStepSettings step;
    step.timeStepSeconds = 0.25;
    step.substeps = 1;
    step.constraintIterations = 4;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(step);
    const auto movedState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto movedEpoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, movedState, grid(), fixture.transfer);
    const auto moved = buildSceneFluidCellVolumes(
        fixture.surface.definition, movedState, grid(), fixture.transfer,
        movedEpoch);
    const auto movedRegions = regionVolumes(moved);
    check(diagnostics.finite
              && moved.fingerprint != first.fingerprint
              && moved.gridEpochFingerprint != first.gridEpochFingerprint,
          "accepted structural motion creates a distinct bound cell-volume epoch");
    for (const StableId region : {StableId{1}, StableId{2}, StableId{3}}) {
        checkNear(movedRegions.at(region).summedCellVolumeCubicMeters,
                  regions.at(region).summedCellVolumeCubicMeters,
                  3.0e-12,
                  "rigid accepted translation preserves each physical region volume");
    }
}

void testFullInteriorCellClassification() {
    Fixture fixture(largeClosedScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "large closed volume fixture assembles");
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto epoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, fineGrid(), fixture.transfer);
    check(epoch.facePartitions.unresolvedActiveFaceCount > 0,
          "large tetrahedron exercises boundary-crossing open face chains");
    const auto volumes = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, fineGrid(), fixture.transfer,
        epoch);
    const auto regions = regionVolumes(volumes);
    constexpr double insideVolume = 85.518598 / 6.0;
    checkNear(regions.at(2).summedCellVolumeCubicMeters,
              insideVolume, 2.0e-11,
              "large tetrahedron recovers analytic volume across full cells");
    checkNear(regions.at(1).summedCellVolumeCubicMeters,
              64.0 - insideVolume, 2.0e-11,
              "large tetrahedron leaves the analytic outside volume");
    std::size_t fullInsideCellCount = 0;
    for (const auto& cell : volumes.cells) {
        if (cell.regionVolumeCount != 1) {
            continue;
        }
        const auto& region =
            volumes.cellRegionVolumes[cell.firstRegionVolume];
        if (region.regionId == 2) {
            ++fullInsideCellCount;
            checkNear(region.volumeFraction, 1.0, 1.0e-12,
                      "classified full interior cell has unit volume fraction");
        }
    }
    check(fullInsideCellCount == 24
              && volumes.tetrahedronCellClipCount > 0
              && volumes.nonzeroTetrahedronCellClipCount > 0
              && volumes.nonzeroTetrahedronCellClipCount
                  <= volumes.tetrahedronCellClipCount
              && volumes.maximumTetrahedronVolumeResidualCubicMeters
                  < 2.0e-12,
          "signed tetrahedron clipping reaches full interior cells");
    validateSceneFluidCellVolumes(
        volumes, fixture.surface.definition, state, fineGrid(),
        fixture.transfer, epoch);

    Fixture translated(translatedLargeClosedScene());
    const auto translatedState = captureSceneFluidSurfaceState(
        translated.surface.definition,
        translated.structureAssembly.mappings,
        translated.structure);
    const auto translatedEpoch = buildSceneFluidGridEpoch(
        translated.surface.definition, translatedState, translatedFineGrid(),
        translated.transfer);
    const auto translatedVolumes = buildSceneFluidCellVolumes(
        translated.surface.definition, translatedState, translatedFineGrid(),
        translated.transfer, translatedEpoch);
    const auto translatedRegions = regionVolumes(translatedVolumes);
    checkNear(translatedRegions.at(2).summedCellVolumeCubicMeters,
              insideVolume, 3.0e-11,
              "signed clipping is invariant to a nonzero grid origin");
    checkNear(translatedRegions.at(1).summedCellVolumeCubicMeters,
              64.0 - insideVolume, 3.0e-11,
              "translated grid preserves the outside region volume");
}

void testPlanarOpeningCapVolume() {
    Fixture fixture(openTetraScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "open tetrahedron volume fixture assembles");
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto epoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer);
    const auto volumes = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, grid(), fixture.transfer, epoch);
    const auto regions = regionVolumes(volumes);
    check(volumes.openingCapFingerprint != 0
              && volumes.openingCapCount == 1,
          "open-cell volume binds one explicit virtual cap");
    checkNear(volumes.openingCapAreaSquareMeters, 0.18, 1.0e-14,
              "open-cell volume retains the analytic mouth area");
    checkNear(regions.at(2).summedCellVolumeCubicMeters,
              0.096, 2.0e-12,
              "virtual mouth cap recovers the open tetrahedron volume");
    checkNear(regions.at(1).summedCellVolumeCubicMeters,
              64.0 - 0.096, 2.0e-12,
              "virtual mouth cap preserves the outside volume");
    validateSceneFluidCellVolumes(
        volumes, fixture.surface.definition, state, grid(),
        fixture.transfer, epoch);
}

void testUnsupportedTopologyCorruptionAndLimits() {
    Fixture open(openScene());
    check(open.surface.ok() && open.structureAssembly.ok(),
          "open cell-volume rejection fixture assembles");
    const auto openState = captureSceneFluidSurfaceState(
        open.surface.definition,
        open.structureAssembly.mappings,
        open.structure);
    const auto openEpoch = buildSceneFluidGridEpoch(
        open.surface.definition, openState, grid(), open.transfer);
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            open.surface.definition, openState, grid(), open.transfer,
            openEpoch)); },
        "zero-thickness capped region is rejected as nonphysical");

    auto nonManifoldScene = openScene();
    nonManifoldScene.openings.clear();
    Fixture nonManifold(std::move(nonManifoldScene));
    check(nonManifold.surface.ok() && nonManifold.structureAssembly.ok(),
          "open surface without authored intake assembles for volume rejection");
    const auto nonManifoldState = captureSceneFluidSurfaceState(
        nonManifold.surface.definition,
        nonManifold.structureAssembly.mappings,
        nonManifold.structure);
    const auto nonManifoldEpoch = buildSceneFluidGridEpoch(
        nonManifold.surface.definition, nonManifoldState, grid(),
        nonManifold.transfer);
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            nonManifold.surface.definition, nonManifoldState, grid(),
            nonManifold.transfer, nonManifoldEpoch)); },
        "open non-manifold surface rejects without relying on opening metadata");

    Fixture fixture(nestedScene());
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto epoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer);
    const auto accepted = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, grid(), fixture.transfer, epoch);
    auto corrupt = accepted;
    corrupt.cellRegionVolumes.front().volumeCubicMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidCellVolumes(
            corrupt, fixture.surface.definition, state, grid(),
            fixture.transfer, epoch); },
        "cell-volume validation rejects nested payload corruption");

    SceneFluidCellVolumeSettings invalidSettings;
    invalidSettings.absoluteVolumeToleranceCubicMeters = 0.0;
    invalidSettings.relativeVolumeTolerance = 0.0;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            fixture.surface.definition, state, grid(), fixture.transfer,
            epoch, invalidSettings)); },
        "cell-volume construction rejects zero closure authority");

    SceneFluidCellVolumeLimits limits;
    limits.openingCaps.maximumCaps = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            open.surface.definition, openState, grid(), open.transfer,
            openEpoch, {}, limits)); },
        "cell-volume opening-cap work is bounded");
    limits = {};
    limits.maximumTetrahedronCellClips = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            fixture.surface.definition, state, grid(), fixture.transfer,
            epoch, {}, limits)); },
        "cell-volume tetrahedron clipping work is bounded");
    limits = {};
    limits.maximumContributionEvents = grid().cellCount() - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            fixture.surface.definition, state, grid(), fixture.transfer,
            epoch, {}, limits)); },
        "cell-volume contribution work is bounded before publication");
    limits = {};
    limits.maximumCellRegionVolumes = grid().cellCount() - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            fixture.surface.definition, state, grid(), fixture.transfer,
            epoch, {}, limits)); },
        "cell-volume sparse output count is bounded");
    limits = {};
    limits.maximumVolumeBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidCellVolumes(
            fixture.surface.definition, state, grid(), fixture.transfer,
            epoch, {}, limits)); },
        "cell-volume result bytes are bounded");
    validateSceneFluidCellVolumes(
        accepted, fixture.surface.definition, state, grid(),
        fixture.transfer, epoch);
}

} // namespace

int main() {
    try {
        testNestedAnalyticVolumesAndAcceptedMotion();
        testFullInteriorCellClassification();
        testPlanarOpeningCapVolume();
        testUnsupportedTopologyCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid cell-volume check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid cell-volume checks passed");
    return 0;
}
