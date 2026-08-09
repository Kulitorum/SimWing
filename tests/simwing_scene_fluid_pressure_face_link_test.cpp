#include "scene_fluid_pressure_face_link.h"
#include "scene_fluid_pressure_operator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <utility>

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
        "sha256:scene-fluid-pressure-face-links";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-face-link-test/1";
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
        "sha256:scene-fluid-pressure-face-open";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-face-link-test/1";
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

Scene tiltedOpenScene() {
    auto scene = openScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-face-tilted-open";
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
              surface.definition, volumes, connectivity)) {}

    SceneFluidPressureFaceLinkSet links(
        const SceneFluidPressureFaceLinkSettings& settings = {},
        const SceneFluidPressureFaceLinkLimits& limits = {}) const {
        return buildSceneFluidPressureFaceLinks(
            surface.definition, state, grid(), transfer, epoch, caps,
            openingQuadrature, openingPatches, volumes, connectivity,
            pressureVolumes, settings, limits);
    }
};

void testExactNestedFaceLinks() {
    Fixture fixture(nestedScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "pressure-face nested fixture assembles");
    const auto first = fixture.links();
    const auto repeated = fixture.links();
    check(first == repeated
              && first.version == sceneFluidPressureFaceLinkVersion
              && first.fingerprint != 0
              && first.faces.size() == 3 * grid().cellCount()
              && first.resolvedPartitionFaceCount == 1
              && first.unresolvedActiveFaceCount == 0
              && first.unresolvedAmbiguousFaceCount == 0,
          "nested pressure face links resolve deterministically without fallback");
    const auto partitionFace = std::ranges::find(
        first.faces,
        SceneFluidPressureFaceStatus::ResolvedPartition,
        &SceneFluidPressureFace::status);
    check(partitionFace != first.faces.end()
              && partitionFace->linkCount == 3,
          "nested interface face owns one exact link per region area");
    std::map<StableId, double> areas;
    for (std::size_t offset = 0;
         offset < partitionFace->linkCount; ++offset) {
        const auto& link = first.links[partitionFace->firstLink + offset];
        areas.emplace(link.minusRegionId, link.areaSquareMeters);
        const auto& minus = fixture.pressureVolumes.controlVolumes[
            link.minusControlVolumeIndex];
        const auto& plus = fixture.pressureVolumes.controlVolumes[
            link.plusControlVolumeIndex];
        check(link.kind == SceneFluidPressureFaceLinkKind::SameRegion
                  && link.minusRegionId == link.plusRegionId
                  && minus.regionId == link.minusRegionId
                  && plus.regionId == link.plusRegionId
                  && minus.componentIndex == link.componentIndex
                  && plus.componentIndex == link.componentIndex,
              "pressure face link joins matching same-region control volumes");
        checkNear(link.geometryWeightMeters,
                  link.areaSquareMeters / link.centerDistanceMeters,
                  0.0,
                  "pressure face link retains its exact geometric weight");
    }
    checkNear(areas.at(1), 0.955, 3.0e-15,
              "partitioned pressure face retains exterior area");
    checkNear(areas.at(2), 0.04, 3.0e-15,
              "partitioned pressure face retains annular-cell area");
    checkNear(areas.at(3), 0.005, 3.0e-15,
              "partitioned pressure face retains inner-cell area");
    check(first.resolvedFullFaceCount + first.resolvedPartitionFaceCount
              + first.resolvedOpeningFaceCount
              == first.faces.size()
              && first.maximumResolvedAreaResidualSquareMeters < 4.0e-15,
          "every nested-scene periodic face closes its exact linked area");
    validateSceneFluidPressureFaceLinks(
        first, fixture.surface.definition, fixture.state, grid(),
        fixture.transfer, fixture.epoch, fixture.caps,
        fixture.openingQuadrature, fixture.openingPatches,
        fixture.volumes, fixture.connectivity, fixture.pressureVolumes);
}

void testFaceAndEmbeddedOpeningLinks() {
    Fixture fixture(openScene());
    const auto links = fixture.links();
    const auto faceOwnedPatchCount = std::ranges::count(
        fixture.openingPatches.patches,
        SceneFluidOpeningPatchOwnerKind::Face,
        &SceneFluidOpeningGridPatch::ownerKind);
    check(faceOwnedPatchCount == 1
              && links.resolvedOpeningFaceCount == 1
              && links.unresolvedOpeningFaceCount == 0,
          "face-aligned authored opening owns one resolved pressure face");
    const auto openingFace = std::ranges::find(
        links.faces,
        SceneFluidPressureFaceStatus::ResolvedOpening,
        &SceneFluidPressureFace::status);
    check(openingFace != links.faces.end()
              && openingFace->linkCount == 2,
          "resolved opening face owns aperture and complementary links");
    const SceneFluidPressureFaceLink* aperture = nullptr;
    const SceneFluidPressureFaceLink* complement = nullptr;
    for (std::size_t offset = 0;
         offset < openingFace->linkCount; ++offset) {
        const auto& link = links.links[openingFace->firstLink + offset];
        if (link.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            aperture = &link;
        } else {
            complement = &link;
        }
    }
    check(aperture != nullptr
              && aperture->minusRegionId == 2
              && aperture->plusRegionId == 1
              && aperture->openingId == 700
              && aperture->openingPatchStableId != 0,
          "opening pressure link follows authored negative-to-positive orientation");
    check(complement != nullptr
              && complement->minusRegionId == 1
              && complement->plusRegionId == 1,
          "opening pressure face retains its same-region complement");
    if (aperture != nullptr && complement != nullptr) {
        checkNear(aperture->areaSquareMeters, 0.18, 2.0e-15,
                  "opening pressure link retains exact aperture area");
        checkNear(complement->areaSquareMeters, 0.82, 2.0e-15,
                  "opening pressure link retains exact complementary area");
    }
    checkNear(openingFace->areaResidualSquareMeters, 0.0, 2.0e-15,
              "opening and complementary links close the Cartesian face");

    Fixture tilted(tiltedOpenScene());
    const auto embeddedLinks = tilted.links();
    const auto cellOwnedPatchCount = std::ranges::count(
        tilted.openingPatches.patches,
        SceneFluidOpeningPatchOwnerKind::Cell,
        &SceneFluidOpeningGridPatch::ownerKind);
    const auto embedded = std::ranges::find(
        embeddedLinks.links,
        SceneFluidPressureLinkGeometryKind::EmbeddedOpening,
        &SceneFluidPressureFaceLink::geometryKind);
    check(cellOwnedPatchCount == 1
              && embeddedLinks.embeddedOpeningLinkCount == 1
              && embeddedLinks.unresolvedEmbeddedOpeningPatchCount == 0
              && embeddedLinks.unresolvedEmbeddedOpeningAreaSquareMeters
                  == 0.0
              && embedded != embeddedLinks.links.end()
              && embeddedLinks.unresolvedOpeningFaceCount == 0,
          "tilted off-face intake owns one embedded pressure link");
    if (embedded != embeddedLinks.links.end()) {
        const auto& patch = tilted.openingPatches.patches.front();
        const auto& minus = tilted.pressureVolumes.controlVolumes[
            embedded->minusControlVolumeIndex];
        const auto& plus = tilted.pressureVolumes.controlVolumes[
            embedded->plusControlVolumeIndex];
        const double projectedDistance =
            (plus.centroidMeters.x - minus.centroidMeters.x)
                * patch.unitNormalNegativeToPositive.x
            + (plus.centroidMeters.y - minus.centroidMeters.y)
                * patch.unitNormalNegativeToPositive.y
            + (plus.centroidMeters.z - minus.centroidMeters.z)
                * patch.unitNormalNegativeToPositive.z;
        check(embedded->kind
                      == SceneFluidPressureFaceLinkKind::AuthoredOpening
                  && embedded->faceIndex
                      == invalidSceneFluidPressureFaceIndex
                  && std::abs(patch.unitNormalNegativeToPositive.y) > 0.1
                  && minus.cellIndex == patch.cellIndex
                  && plus.cellIndex == patch.cellIndex
                  && minus.regionId == patch.negativeSideRegionId
                  && plus.regionId == patch.positiveSideRegionId,
              "embedded intake joins its exact same-cell side controls");
        checkNear(embedded->areaSquareMeters, patch.areaSquareMeters, 0.0,
                  "embedded intake preserves exact clipped patch area");
        checkNear(embedded->centerDistanceMeters, projectedDistance, 0.0,
                  "embedded intake uses projected pressure-centroid separation");
        checkNear(embedded->geometryWeightMeters,
                  patch.areaSquareMeters / projectedDistance, 0.0,
                  "embedded intake publishes exact centroid conductance weight");
        checkNear(embedded->unitNormalMinusToPlus.x,
                  patch.unitNormalNegativeToPositive.x, 0.0,
                  "embedded intake preserves authored normal x");
        checkNear(embedded->unitNormalMinusToPlus.y,
                  patch.unitNormalNegativeToPositive.y, 0.0,
                  "embedded intake preserves authored normal y");
        checkNear(embedded->unitNormalMinusToPlus.z,
                  patch.unitNormalNegativeToPositive.z, 0.0,
                  "embedded intake preserves authored normal z");
    }
    checkNear(embeddedLinks.totalEmbeddedOpeningAreaSquareMeters,
              tilted.openingPatches.totalAreaSquareMeters, 0.0,
              "embedded pressure links consume every off-face opening patch");
}

void testCorruptionSettingsAndLimits() {
    Fixture fixture(nestedScene());
    const auto accepted = fixture.links();
    auto corrupt = accepted;
    corrupt.links.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidPressureFaceLinks(
            corrupt, fixture.surface.definition, fixture.state, grid(),
            fixture.transfer, fixture.epoch, fixture.caps,
            fixture.openingQuadrature, fixture.openingPatches,
            fixture.volumes, fixture.connectivity,
            fixture.pressureVolumes); },
        "pressure-face-link validation rejects nested corruption");

    SceneFluidPressureFaceLinkSettings invalidSettings;
    invalidSettings.areaToleranceSquareMeters = 0.0;
    expectInvalid(
        [&] { static_cast<void>(fixture.links(invalidSettings)); },
        "pressure-face-link assembly rejects zero area tolerance");
    invalidSettings = {};
    invalidSettings.minimumCenterDistanceMeters = 0.0;
    expectInvalid(
        [&] { static_cast<void>(fixture.links(invalidSettings)); },
        "pressure-face-link assembly rejects zero centroid distance tolerance");

    Fixture tilted(tiltedOpenScene());
    SceneFluidPressureFaceLinkSettings excessiveDistance;
    excessiveDistance.minimumCenterDistanceMeters = 10.0;
    const auto unresolvedEmbedded = tilted.links(excessiveDistance);
    check(unresolvedEmbedded.embeddedOpeningLinkCount == 0
              && unresolvedEmbedded.unresolvedEmbeddedOpeningPatchCount == 1
              && unresolvedEmbedded.unresolvedEmbeddedOpeningAreaSquareMeters
                  == tilted.openingPatches.totalAreaSquareMeters,
          "embedded pressure link retains degenerate centroid separation as unresolved");
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidPressureOperator(
            tilted.surface.definition, tilted.state, grid(),
            tilted.transfer, tilted.epoch, tilted.caps,
            tilted.openingQuadrature, tilted.openingPatches,
            tilted.volumes, tilted.connectivity,
            tilted.pressureVolumes, unresolvedEmbedded)); },
        "pressure operator rejects unresolved embedded openings");

    SceneFluidPressureFaceLinkLimits limits;
    limits.maximumFaces = 3 * grid().cellCount() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.links({}, limits)); },
        "pressure-face-link assembly bounds face count");
    limits = {};
    limits.maximumLinks = accepted.links.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.links({}, limits)); },
        "pressure-face-link assembly bounds link count");
    limits = {};
    limits.maximumLinkBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.links({}, limits)); },
        "pressure-face-link assembly bounds owned storage");
}

} // namespace

int main() {
    try {
        testExactNestedFaceLinks();
        testFaceAndEmbeddedOpeningLinks();
        testCorruptionSettingsAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid pressure-face-link check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-face-link checks passed");
    return 0;
}
