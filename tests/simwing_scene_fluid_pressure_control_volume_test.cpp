#include "scene_fluid_pressure_control_volume.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
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

std::array<Vec2, 3> intrinsicChart(
    const std::array<Vec3, 4>& positions,
    const std::array<std::size_t, 3>& vertices) {
    const Vec3& first = positions[vertices[0]];
    const Vec3& second = positions[vertices[1]];
    const Vec3& third = positions[vertices[2]];
    const Vec3 edge{second.x - first.x,
                    second.y - first.y,
                    second.z - first.z};
    const Vec3 diagonal{third.x - first.x,
                        third.y - first.y,
                        third.z - first.z};
    const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
    const double projected = (diagonal.x * edge.x
                              + diagonal.y * edge.y
                              + diagonal.z * edge.z) / edgeLength;
    const double diagonalSquared = diagonal.x * diagonal.x
        + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
    return {{{0.0, 0.0},
             {edgeLength, 0.0},
             {projected, std::sqrt(std::max(
                 0.0, diagonalSquared - projected * projected))}}};
}

Scene tetraScene(const std::array<Vec3, 4>& positions,
                 const bool openMouth) {
    Scene scene;
    scene.metadata.designChecksum = openMouth
        ? "sha256:scene-fluid-pressure-open-tetra"
        : "sha256:scene-fluid-pressure-large-tetra";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-control-volume-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}},
    }};
    const std::size_t faceCount = openMouth ? 3 : 4;
    for (std::size_t face = 0; face < faceCount; ++face) {
        scene.triangles.push_back({
            500 + face,
            {10 + faces[face][0], 10 + faces[face][1],
             10 + faces[face][2]},
            intrinsicChart(positions, faces[face]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    if (openMouth) {
        scene.openings = {
            {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
        };
    }
    return scene;
}

Scene largeClosedScene() {
    return tetraScene(
        {{{0.11, 0.18, 0.27},
          {3.83, 0.31, 3.57},
          {3.69, 3.77, 0.22},
          {0.27, 3.63, 3.84}}},
        false);
}

Scene openTetraScene() {
    return tetraScene(
        {{{1.2, 1.5, 1.45},
          {2.8, 1.2, 1.15},
          {2.8, 1.8, 1.15},
          {2.8, 1.5, 1.75}}},
        true);
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

fluid::PeriodicCartesianGrid fineGrid() {
    return {{8, 8, 8}, {}, {4.0, 4.0, 4.0}};
}

struct Fixture {
    Scene scene;
    fluid::PeriodicCartesianGrid grid;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;
    SceneFluidRegionConnectivity connectivity;

    Fixture(Scene source, fluid::PeriodicCartesianGrid sourceGrid)
        : scene(std::move(source)),
          grid(std::move(sourceGrid)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(surface.definition,
                   structureAssembly.mappings,
                   structure),
          connectivity(buildSceneFluidRegionConnectivity(
              surface.definition)) {}

    SceneFluidCellVolumeSet volumes() const {
        const auto state = captureSceneFluidSurfaceState(
            surface.definition, structureAssembly.mappings, structure);
        const auto epoch = buildSceneFluidGridEpoch(
            surface.definition, state, grid, transfer);
        return buildSceneFluidCellVolumes(
            surface.definition, state, grid, transfer, epoch);
    }
};

void testSparsePressureVolumesAndGauges() {
    Fixture fixture(largeClosedScene(), fineGrid());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "pressure-volume large tetrahedron fixture assembles");
    const auto volumes = fixture.volumes();
    const auto first = buildSceneFluidPressureControlVolumes(
        fixture.surface.definition, volumes, fixture.connectivity);
    const auto repeated = buildSceneFluidPressureControlVolumes(
        fixture.surface.definition, volumes, fixture.connectivity);
    check(first == repeated
              && first.version == sceneFluidPressureControlVolumeVersion
              && first.fingerprint != 0
              && first.cells.size() == fineGrid().cellCount()
              && first.controlVolumes.size()
                  == volumes.cellRegionVolumes.size()
              && first.regions.size() == 2
              && first.components.size() == 2,
          "sparse pressure volumes are deterministic and preserve all owners");
    std::size_t mixedCellCount = 0;
    std::size_t fullInsideCellCount = 0;
    for (const auto& cell : first.cells) {
        mixedCellCount += cell.controlVolumeCount > 1;
        double volume = 0.0;
        for (std::size_t offset = 0;
             offset < cell.controlVolumeCount; ++offset) {
            const auto& control = first.controlVolumes[
                cell.firstControlVolume + offset];
            volume += control.volumeCubicMeters;
            fullInsideCellCount += cell.controlVolumeCount == 1
                && control.regionId == 2;
        }
        checkNear(volume, first.cellVolumeCubicMeters, 2.0e-12,
                  "pressure unknowns close each Cartesian cell");
    }
    check(mixedCellCount > 0 && fullInsideCellCount == 24,
          "pressure topology retains mixed and full interior control volumes");
    std::set<std::uint64_t> stableIds;
    for (const auto& control : first.controlVolumes) {
        check(control.stableId != 0
                  && stableIds.insert(control.stableId).second,
              "pressure control-volume IDs are nonzero and unique");
    }
    for (const auto& component : first.components) {
        const auto& gauge =
            first.controlVolumes[component.gaugeControlVolumeIndex];
        check(gauge.componentIndex == component.componentIndex
                  && gauge.regionId == component.gaugeRegionId
                  && gauge.belongsToGaugeRegion,
              "each disconnected pressure component owns a canonical gauge volume");
    }
    checkNear(first.domainVolumeCubicMeters, 64.0, 0.0,
              "pressure volumes retain the exact grid-domain volume");
    checkNear(first.domainVolumeResidualCubicMeters, 0.0, 3.0e-11,
              "pressure volumes close globally");
    validateSceneFluidPressureControlVolumes(
        first, fixture.surface.definition, volumes, fixture.connectivity);
}

void testIntakeConnectivityAndAcceptedMotionBinding() {
    Fixture open(openTetraScene(), grid());
    const auto openVolumes = open.volumes();
    const auto openPressure = buildSceneFluidPressureControlVolumes(
        open.surface.definition, openVolumes, open.connectivity);
    check(open.connectivity.components.size() == 1
              && openPressure.components.size() == 1
              && openPressure.components.front().gaugeRegionId == 1,
          "intake connectivity places cell and Outside pressure volumes under one gauge");
    check(std::ranges::all_of(
              openPressure.controlVolumes,
              [](const SceneFluidPressureControlVolume& control) {
                  return control.componentIndex == 0;
              }),
          "every intake-linked pressure volume maps to the joined component");

    Fixture moving(largeClosedScene(), fineGrid());
    const auto previousVolumes = moving.volumes();
    const auto previous = buildSceneFluidPressureControlVolumes(
        moving.surface.definition, previousVolumes, moving.connectivity);
    for (std::size_t node = 0;
         node < moving.structure.definition().nodes.size(); ++node) {
        const double mass = moving.structure.definition().nodes[node].massKg;
        moving.structure.addExternalForce(node, {mass * 0.1, 0.0, 0.0});
    }
    StructureStepSettings step;
    step.timeStepSeconds = 0.25;
    step.substeps = 1;
    step.constraintIterations = 0;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto diagnostics = moving.structure.step(step);
    const auto currentVolumes = moving.volumes();
    const auto current = buildSceneFluidPressureControlVolumes(
        moving.surface.definition, currentVolumes, moving.connectivity);
    check(diagnostics.finite
              && current.fingerprint != previous.fingerprint
              && current.surfaceStateFingerprint
                  != previous.surfaceStateFingerprint
              && current.acceptedStepCount
                  == previous.acceptedStepCount + 1,
          "accepted Structure motion creates a distinct bound pressure-volume epoch");
    std::map<std::uint64_t, std::pair<std::size_t, StableId>> previousKeys;
    for (const auto& control : previous.controlVolumes) {
        previousKeys.emplace(
            control.stableId,
            std::pair{control.cellIndex, control.regionId});
    }
    std::size_t persistentCellControlCount = 0;
    for (const auto& control : current.controlVolumes) {
        const auto found = previousKeys.find(control.stableId);
        if (found == previousKeys.end()) continue;
        ++persistentCellControlCount;
        check(found->second
                  == std::pair{control.cellIndex, control.regionId},
              "persistent pressure ID retains its exact cell-region key");
    }
    check(persistentCellControlCount >= fineGrid().cellCount(),
          "accepted motion preserves stable IDs for persistent cell-region pairs");
}

void testCorruptionAndLimits() {
    Fixture fixture(openTetraScene(), grid());
    const auto volumes = fixture.volumes();
    const auto accepted = buildSceneFluidPressureControlVolumes(
        fixture.surface.definition, volumes, fixture.connectivity);
    auto corrupt = accepted;
    corrupt.controlVolumes.front().volumeCubicMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidPressureControlVolumes(
            corrupt, fixture.surface.definition, volumes,
            fixture.connectivity); },
        "pressure-control-volume validation rejects nested corruption");

    auto corruptSource = volumes;
    corruptSource.cellRegionVolumes.front().volumeCubicMeters += 0.01;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, corruptSource,
            fixture.connectivity)); },
        "pressure-control-volume assembly rejects a corrupted source epoch");

    SceneFluidPressureControlVolumeLimits limits;
    limits.maximumCells = grid().cellCount() - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds cell count");
    limits = {};
    limits.maximumControlVolumes = volumes.cellRegionVolumes.size() - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds sparse unknown count");
    limits = {};
    limits.maximumRegions = 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds regions");
    limits = {};
    limits.maximumOpenings = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds authored openings");
    limits = {};
    limits.maximumComponents = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds pressure components");
    limits = {};
    limits.maximumControlVolumeBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureControlVolumes(
            fixture.surface.definition, volumes,
            fixture.connectivity, limits)); },
        "pressure-control-volume assembly bounds owned storage");
}

} // namespace

int main() {
    try {
        testSparsePressureVolumesAndGauges();
        testIntakeConnectivityAndAcceptedMotionBinding();
        testCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid pressure-control-volume check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-control-volume checks passed");
    return 0;
}
