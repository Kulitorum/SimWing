#include "fluid/mimetic_local_cell.h"
#include "scene_fluid_mimetic_control_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
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
                     "FAIL: %s (actual %.17g expected %.17g tolerance %.3g)\n",
                     message, actual, expected, tolerance);
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
    const std::array<Vec2, 3> chart{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
    }};
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
        "sha256:scene-fluid-mimetic-control-cell";
    scene.metadata.exporterVersion =
        "scene-fluid-mimetic-control-cell-test/1";
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

Scene openScene(const bool tilted) {
    Scene scene;
    scene.metadata.designChecksum = tilted
        ? "sha256:scene-fluid-mimetic-tilted-open"
        : "sha256:scene-fluid-mimetic-face-open";
    scene.metadata.exporterVersion =
        "scene-fluid-mimetic-control-cell-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, 2.0, 0.3, 0.3, 900, true);
    if (tilted) {
        scene.vertices[1].positionMeters.x = 2.6;
        scene.vertices[2].positionMeters.x = 2.8;
        scene.vertices[3].positionMeters.x = 2.7;
    }
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
              pressureVolumes)) {}

    SceneFluidMimeticControlCellSet mimetic(
        const SceneFluidMimeticControlCellSettings& settings = {},
        const SceneFluidMimeticControlCellLimits& limits = {}) const {
        return mimetic(faceLinks, settings, limits);
    }

    SceneFluidMimeticControlCellSet mimetic(
        const SceneFluidPressureFaceLinkSet& selectedFaceLinks,
        const SceneFluidMimeticControlCellSettings& settings = {},
        const SceneFluidMimeticControlCellLimits& limits = {}) const {
        return buildSceneFluidMimeticControlCells(
            surface.definition, state, grid(), epoch, caps,
            openingQuadrature, openingPatches, pressureVolumes,
            selectedFaceLinks, settings, limits);
    }

    void validate(const SceneFluidMimeticControlCellSet& value) const {
        validateSceneFluidMimeticControlCells(
            value, surface.definition, state, grid(), epoch, caps,
            openingQuadrature, openingPatches, pressureVolumes,
            faceLinks);
    }
};

void checkLocalOperators(
    const SceneFluidMimeticControlCellSet& set,
    const char* message) {
    bool allBuilt = true;
    for (const auto& cell : set.controlCells) {
        if (!cell.readyForMimeticOperator) {
            allBuilt = false;
            continue;
        }
        fluid::MimeticLocalCellGeometry geometry;
        geometry.volumeCubicMeters = cell.volumeCubicMeters;
        geometry.centroidMeters = cell.centroidMeters;
        for (std::size_t offset = 0; offset < cell.halfFaceCount; ++offset) {
            const auto& face = set.halfFaces[cell.firstHalfFace + offset];
            geometry.halfFaces.push_back({
                face.areaSquareMeters,
                face.centroidMeters,
                face.outwardUnitNormal,
            });
        }
        try {
            static_cast<void>(
                fluid::buildMimeticLocalCellOperator(geometry));
        } catch (const std::exception&) {
            allBuilt = false;
        }
    }
    check(allBuilt, message);
}

void checkPairedHalfFaces(
    const SceneFluidMimeticControlCellSet& set,
    const char* message) {
    using Key = std::pair<SceneFluidMimeticHalfFaceKind, std::uint64_t>;
    std::map<Key, std::vector<const SceneFluidMimeticHalfFace*>> traces;
    std::map<std::uint64_t,
             std::vector<const SceneFluidMimeticHalfFace*>> walls;
    for (const auto& face : set.halfFaces) {
        if (face.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            walls[face.sourceStableId].push_back(&face);
        } else {
            traces[{face.kind, face.traceStableId}].push_back(&face);
        }
    }
    bool paired = std::ranges::all_of(
        traces,
        [](const auto& item) {
            if (item.second.size() != 2) {
                return false;
            }
            const auto& first = *item.second[0];
            const auto& second = *item.second[1];
            return first.areaSquareMeters == second.areaSquareMeters
                && first.otherControlVolumeIndex
                    == second.controlVolumeIndex
                && second.otherControlVolumeIndex
                    == first.controlVolumeIndex
                && first.outwardUnitNormal.x
                    == -second.outwardUnitNormal.x
                && first.outwardUnitNormal.y
                    == -second.outwardUnitNormal.y
                && first.outwardUnitNormal.z
                    == -second.outwardUnitNormal.z;
        });
    paired = paired && std::ranges::all_of(
        walls,
        [](const auto& item) { return item.second.size() == 2; });
    check(paired, message);
}

void testNestedClosedShells() {
    Fixture fixture(nestedScene());
    const auto first = fixture.mimetic();
    const auto repeated = fixture.mimetic();
    check(first == repeated
              && first.version == sceneFluidMimeticControlCellVersion
              && first.fingerprint != 0
              && first.controlCells.size()
                  == fixture.pressureVolumes.controlVolumes.size()
              && first.readyControlCellCount == first.controlCells.size()
              && first.incompleteTopologyControlCellCount == 0
              && first.nonclosingControlCellCount == 0
              && first.unresolvedCartesianFaceCount == 0
              && first.omittedZeroVolumeMaterialSideCount == 0
              && first.missingOpeningControlSideCount == 0
              && first.cartesianHalfFaceCount > 0
              && first.materialWallHalfFaceCount > 0
              && first.openingHalfFaceCount == 0,
          "nested regions assemble deterministic complete mimetic shells");
    check(first.maximumAreaClosureErrorSquareMeters < 2.0e-12
              && first.maximumDivergenceTheoremErrorCubicMeters < 2.0e-12,
          "nested region shells close area vectors and first normal moments");
    checkPairedHalfFaces(
        first, "nested material and trace half-faces retain exact pairs");
    checkLocalOperators(
        first, "every nested scene shell builds the generic SPD local kernel");
    fixture.validate(first);

    bool foundPeriodicImage = false;
    for (const auto& face : first.halfFaces) {
        foundPeriodicImage = foundPeriodicImage
            || face.periodicImageOffsetMeters.x != 0.0
            || face.periodicImageOffsetMeters.y != 0.0
            || face.periodicImageOffsetMeters.z != 0.0;
    }
    check(foundPeriodicImage,
          "periodic boundary half-faces retain an explicit local image offset");
}

void testFaceAndEmbeddedOpeningShells() {
    Fixture faceAligned(openScene(false));
    const auto faceSet = faceAligned.mimetic();
    check(faceAligned.faceLinks.resolvedOpeningFaceCount == 1
              && faceAligned.faceLinks.unresolvedEmbeddedOpeningPatchCount == 0
              && faceSet.readyControlCellCount == faceSet.controlCells.size()
              && faceSet.openingHalfFaceCount == 0
              && faceSet.omittedZeroVolumeMaterialSideCount == 0
              && faceSet.missingOpeningControlSideCount == 0,
          "face-aligned opening is retained through its Cartesian trace pair");
    checkPairedHalfFaces(
        faceSet, "face-aligned opening shell retains paired trace ownership");
    checkLocalOperators(
        faceSet, "face-aligned opening shells build local SPD kernels");

    Fixture embedded(openScene(true));
    auto rejectionSettings = SceneFluidPressureFaceLinkSettings{};
    rejectionSettings.minimumCenterDistanceMeters = 10.0;
    const auto rejectedLinks = buildSceneFluidPressureFaceLinks(
        embedded.surface.definition, embedded.state, grid(),
        embedded.transfer, embedded.epoch, embedded.caps,
        embedded.openingQuadrature, embedded.openingPatches,
        embedded.openingFaceCrossings, embedded.cappedFacePartitions,
        embedded.volumes, embedded.connectivity, embedded.pressureVolumes,
        rejectionSettings);
    const auto embeddedSet = embedded.mimetic(rejectedLinks);
    const bool embeddedComplete =
        rejectedLinks.unresolvedEmbeddedOpeningPatchCount == 1
              && embeddedSet.openingHalfFaceCount == 2
              && embeddedSet.readyControlCellCount
                  == embeddedSet.controlCells.size()
              && embeddedSet.incompleteTopologyControlCellCount == 0
              && embeddedSet.nonclosingControlCellCount == 0
              && embeddedSet.omittedZeroVolumeMaterialSideCount == 0
              && embeddedSet.missingOpeningControlSideCount == 0;
    if (!embeddedComplete) {
        std::fprintf(
            stderr,
            "embedded diagnostics: rejected=%zu opening-halves=%zu "
            "ready=%zu/%zu incomplete=%zu nonclosing=%zu "
            "omitted-wall-sides=%zu missing-opening-sides=%zu "
            "max-area=%.17g max-moment=%.17g\n",
            rejectedLinks.unresolvedEmbeddedOpeningPatchCount,
            embeddedSet.openingHalfFaceCount,
            embeddedSet.readyControlCellCount,
            embeddedSet.controlCells.size(),
            embeddedSet.incompleteTopologyControlCellCount,
            embeddedSet.nonclosingControlCellCount,
            embeddedSet.omittedZeroVolumeMaterialSideCount,
            embeddedSet.missingOpeningControlSideCount,
            embeddedSet.maximumAreaClosureErrorSquareMeters,
            embeddedSet.maximumDivergenceTheoremErrorCubicMeters);
    }
    check(embeddedComplete,
          "non-admissible two-point opening still completes geometric shells");
    checkPairedHalfFaces(
        embeddedSet, "embedded opening half-faces retain exact trace pairs");
    checkLocalOperators(
        embeddedSet, "embedded opening shells build local SPD kernels");
}

void testLimitsAndCorruption() {
    Fixture fixture(nestedScene());
    const auto accepted = fixture.mimetic();
    auto corrupt = accepted;
    corrupt.halfFaces.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidMimeticControlCellIntegrity(corrupt); },
        "mimetic shell validation rejects fingerprinted corruption");

    SceneFluidMimeticControlCellLimits limits;
    limits.maximumControlCells = accepted.controlCells.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.mimetic({}, limits)); },
        "mimetic shell assembly bounds control-cell count");
    limits = {};
    limits.maximumHalfFaces = accepted.halfFaces.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.mimetic({}, limits)); },
        "mimetic shell assembly bounds half-face count");
    limits = {};
    limits.maximumOwnedBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.mimetic({}, limits)); },
        "mimetic shell assembly bounds owned storage");

    auto badSettings = SceneFluidMimeticControlCellSettings{};
    badSettings.unitNormalTolerance = -1.0;
    expectInvalid(
        [&] { static_cast<void>(fixture.mimetic(badSettings)); },
        "mimetic shell assembly rejects invalid tolerances");
}

} // namespace

int main() {
    try {
        testNestedClosedShells();
        testFaceAndEmbeddedOpeningShells();
        testLimitsAndCorruption();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid mimetic control-cell check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid mimetic control-cell checks passed");
    return 0;
}
