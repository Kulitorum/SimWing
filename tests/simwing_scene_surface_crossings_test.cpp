#include "fluid/scene_surface_crossings.h"

#include <cmath>
#include <cstdio>
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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
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

Scene transverseScene(const bool reverse = false) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-surface-crossings";
    scene.metadata.exporterVersion = "scene-surface-crossings-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.vertices = {
        {10, {1.2, 1.2, 1.2}},
        {11, {2.8, 1.2, 1.2}},
        {12, {1.2, 1.8, 1.8}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {500, reverse ? std::array<StableId, 3>{10, 12, 11}
                      : std::array<StableId, 3>{10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {0.0, 0.6}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    return scene;
}

Scene tangentialScene() {
    Scene scene = transverseScene();
    scene.vertices = {
        {10, {2.0, 1.2, 1.2}},
        {11, {2.8, 1.2, 1.2}},
        {12, {2.0, 1.8, 1.8}},
    };
    return scene;
}

Scene coplanarScene() {
    Scene scene = transverseScene();
    scene.vertices = {
        {10, {1.2, 1.2, 2.0}},
        {11, {1.8, 1.2, 2.0}},
        {12, {1.2, 1.8, 2.0}},
    };
    return scene;
}

PeriodicCartesianGrid grid() {
    return PeriodicCartesianGrid({4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

PeriodicCartesianGrid nonAssociativePlaneGrid() {
    constexpr double plane = 2.8991341562581114;
    return PeriodicCartesianGrid(
        {4, 4, 4}, {-2.0 * plane, 0.0, 0.0},
        {2.0 * plane, 4.0, 4.0});
}

Scene nonAssociativePlaneScene() {
    constexpr double plane = 2.8991341562581114;
    Scene scene = transverseScene();
    scene.metadata.designChecksum =
        "sha256:scene-surface-crossings-non-associative-plane";
    scene.vertices = {
        {10, {plane - 0.005, 1.2, 1.2}},
        {11, {plane + 0.005, 1.2, 1.2}},
        {12, {plane - 0.005, 1.8, 1.8}},
    };
    return scene;
}

struct Pipeline {
    Scene scene;
    PeriodicCartesianGrid fluidGrid;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceState state;
    SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections;
    SceneFluidGridPatchSet patches;
    SceneFluidPatchOwnership ownership;

    explicit Pipeline(Scene source,
                      PeriodicCartesianGrid sourceGrid = grid())
        : scene(std::move(source)),
          fluidGrid(std::move(sourceGrid)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          state(captureSceneFluidSurfaceState(
              surface.definition, structureAssembly.mappings, structure)),
          candidates(buildSceneFluidGridCandidates(
              surface.definition, state, fluidGrid)),
          intersections(intersectSceneFluidSurfaceWithGrid(
              surface.definition, state, fluidGrid, candidates)),
          patches(clipSceneFluidSurfaceToCells(
              surface.definition,
              state,
              fluidGrid,
              candidates,
              intersections)),
          ownership(ownSceneFluidSurfacePatches(
              surface.definition,
              state,
              fluidGrid,
              candidates,
              intersections,
              patches)) {}
};

SceneFluidFaceCrossingSet crossings(const Pipeline& pipeline,
                                    const SceneFluidFaceCrossingLimits& limits = {}) {
    return buildSceneFluidFaceCrossings(
        pipeline.surface.definition,
        pipeline.state,
        pipeline.fluidGrid,
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership,
        limits);
}

void testCanonicalAdjacentPlaneCoordinate() {
    const auto offsetGrid = nonAssociativePlaneGrid();
    const auto lower = offsetGrid.lowerMeters();
    const auto spacing = offsetGrid.cellSpacingMeters();
    const double directPlane = lower.x + 3.0 * spacing.x;
    const double steppedPlane =
        (lower.x + 2.0 * spacing.x) + spacing.x;
    check(directPlane != steppedPlane,
          "scene face crossings: fixture exercises non-associative grid-plane arithmetic");

    Pipeline pipeline(nonAssociativePlaneScene(), offsetGrid);
    const auto result = crossings(pipeline);
    check(result.candidateSegmentCount == 2
              && result.unpairedContactSegmentCount == 0
              && result.crossings.size() == 1,
          "scene face crossings: adjacent cells share one canonical grid plane");
    if (result.crossings.size() == 1) {
        const auto& crossing = result.crossings.front();
        check(crossing.axis == GridFaceAxis::X
                  && crossing.i == 3
                  && crossing.first.positionMeters.x == directPlane
                  && crossing.second.positionMeters.x == directPlane,
              "scene face crossings: canonical plane survives both cell clips exactly");
    }
}

void testTransverseCrossing() {
    Pipeline pipeline(transverseScene());
    const auto first = crossings(pipeline);
    const auto second = crossings(pipeline);
    check(first == second && first.fingerprint != 0
              && first.candidateSegmentCount == 2
              && first.unpairedContactSegmentCount == 0
              && first.coplanarAreaPatchCount == 0
              && first.crossings.size() == 1,
          "scene face crossings: exact adjacent segments pair deterministically");
    const auto& crossing = first.crossings.front();
    check(crossing.stableId != 0
              && crossing.axis == GridFaceAxis::X
              && crossing.i == 2 && crossing.j == 1 && crossing.k == 1
              && crossing.lowerCellSourcePatchIndex
                  < crossing.upperCellSourcePatchIndex
              && crossing.triangleId == 500
              && crossing.negativeSideRegionId == 1
              && crossing.positiveSideRegionId == 2
              && crossing.materialId == 100
              && crossing.sheetId == 900
              && crossing.role == SurfaceRole::Skin,
          "scene face crossings: face and authored physical identity are retained");
    checkNear(crossing.first.positionMeters.x, 2.0, 0.0,
              "scene face crossings: first endpoint lies exactly on the face");
    checkNear(crossing.second.positionMeters.x, 2.0, 0.0,
              "scene face crossings: second endpoint lies exactly on the face");
    checkNear(crossing.lengthMeters, std::sqrt(0.18), 2.0e-15,
              "scene face crossings: transverse segment length is analytic");
    checkNear(first.crossingLengthMeters, crossing.lengthMeters, 0.0,
              "scene face crossings: total length counts the crossing once");
    checkNear(crossing.midpointMeters.y, 1.35, 2.0e-15,
              "scene face crossings: midpoint is exact");
    checkNear(crossing.negativeToPositiveDirectionInFace.y,
              -std::sqrt(0.5), 2.0e-15,
              "scene face crossings: in-face direction follows winding in Y");
    checkNear(crossing.negativeToPositiveDirectionInFace.z,
              std::sqrt(0.5), 2.0e-15,
              "scene face crossings: in-face direction follows winding in Z");
    validateSceneFluidFaceCrossings(
        first,
        pipeline.surface.definition,
        pipeline.state,
        pipeline.fluidGrid,
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership);
}

void testTranslatedTransverseCrossing() {
    const Vec3 translation{256.0, -512.0, 1024.0};
    Scene scene = transverseScene();
    for (auto& vertex : scene.vertices) {
        vertex.positionMeters.x += translation.x;
        vertex.positionMeters.y += translation.y;
        vertex.positionMeters.z += translation.z;
    }
    const PeriodicCartesianGrid translatedGrid(
        {4, 4, 4},
        {translation.x, translation.y, translation.z},
        {translation.x + 4.0,
         translation.y + 4.0,
         translation.z + 4.0});
    Pipeline pipeline(std::move(scene), translatedGrid);
    const auto result = crossings(pipeline);
    check(result.candidateSegmentCount == 2
              && result.unpairedContactSegmentCount == 0
              && result.crossings.size() == 1,
          "scene face crossings: distant common translation preserves adjacent pairing");
    if (result.crossings.size() == 1) {
        const auto& crossing = result.crossings.front();
        check(crossing.axis == GridFaceAxis::X
                  && crossing.i == 2 && crossing.j == 1
                  && crossing.k == 1,
              "scene face crossings: distant translation preserves face identity");
        checkNear(crossing.lengthMeters, std::sqrt(0.18), 2.0e-10,
                  "scene face crossings: distant translation preserves length");
        checkNear(crossing.midpointMeters.y - translation.y, 1.35,
                  2.0e-10,
                  "scene face crossings: distant translation preserves local midpoint");
    }
}

void testWindingAndContactClassification() {
    Pipeline reversed(transverseScene(true));
    const auto reversedCrossings = crossings(reversed);
    check(reversedCrossings.crossings.size() == 1,
          "scene face crossings: reversed triangle still has one crossing");
    checkNear(reversedCrossings.crossings.front()
                  .negativeToPositiveDirectionInFace.y,
              std::sqrt(0.5), 2.0e-15,
              "scene face crossings: reversed winding reverses in-face side direction");
    check(reversedCrossings.crossings.front().negativeSideRegionId == 1
              && reversedCrossings.crossings.front().positiveSideRegionId == 2,
          "scene face crossings: reversing winding does not relabel authored sides");

    Pipeline tangential(tangentialScene());
    const auto tangentialCrossings = crossings(tangential);
    check(tangentialCrossings.crossings.empty()
              && tangentialCrossings.candidateSegmentCount == 1
              && tangentialCrossings.unpairedContactSegmentCount == 1,
          "scene face crossings: an unpaired triangle edge remains contact only");

    Pipeline coplanar(coplanarScene());
    const auto coplanarCrossings = crossings(coplanar);
    check(coplanarCrossings.crossings.empty()
              && coplanarCrossings.candidateSegmentCount == 0
              && coplanarCrossings.unpairedContactSegmentCount == 0
              && coplanarCrossings.coplanarAreaPatchCount == 1,
          "scene face crossings: face-owned triangle area is not converted to a line");
}

void testLimitsAndTransactionalValidation() {
    Pipeline pipeline(transverseScene());
    SceneFluidFaceCrossingLimits limits;
    limits.maximumCandidateSegments = 1;
    expectLimited(
        [&] { static_cast<void>(crossings(pipeline, limits)); },
        "scene face crossings: candidate count is bounded");
    limits = {};
    limits.maximumCrossings = 0;
    expectLimited(
        [&] { static_cast<void>(crossings(pipeline, limits)); },
        "scene face crossings: crossing count is bounded");
    limits = {};
    limits.maximumCrossingBytes = sizeof(SceneFluidFaceCrossing) - 1;
    expectLimited(
        [&] { static_cast<void>(crossings(pipeline, limits)); },
        "scene face crossings: crossing storage is byte-bounded");

    const auto accepted = crossings(pipeline);
    auto corrupt = accepted;
    corrupt.crossings.front().lengthMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidFaceCrossings(
            corrupt,
            pipeline.surface.definition,
            pipeline.state,
            pipeline.fluidGrid,
            pipeline.candidates,
            pipeline.intersections,
            pipeline.patches,
            pipeline.ownership); },
        "scene face crossings: corrupted derived geometry is rejected");
    validateSceneFluidFaceCrossings(
        accepted,
        pipeline.surface.definition,
        pipeline.state,
        pipeline.fluidGrid,
        pipeline.candidates,
        pipeline.intersections,
        pipeline.patches,
        pipeline.ownership);
}

} // namespace

int main() {
    testTransverseCrossing();
    testTranslatedTransverseCrossing();
    testCanonicalAdjacentPlaneCoordinate();
    testWindingAndContactClassification();
    testLimitsAndTransactionalValidation();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene surface-crossing test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene surface-crossing tests passed\n");
    return 0;
}
