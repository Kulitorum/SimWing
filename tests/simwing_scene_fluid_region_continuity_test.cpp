#include "scene_fluid_region_continuity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <ranges>
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

template<std::size_t VertexCount>
std::array<Vec2, 3> intrinsicChart(
    const std::array<Vec3, VertexCount>& positions,
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

Scene openSquarePyramidScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-region-continuity";
    scene.metadata.exporterVersion =
        "scene-fluid-region-continuity-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 5> positions{{
        {1.0, 1.0, 1.0},
        {2.0, 0.5, 0.5},
        {2.0, 1.5, 0.5},
        {2.0, 1.5, 1.5},
        {2.0, 0.5, 1.5},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        scene.triangles.push_back({
            500 + face,
            {10 + faces[face][0], 10 + faces[face][1],
             10 + faces[face][2]},
            intrinsicChart(positions, faces[face]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {0.0, 0.0, 0.0}, {4.0, 4.0, 4.0}};
}

SceneStructureAssembly fixedMouthAssembly(const Scene& scene) {
    auto assembly = assembleSceneStructure(scene);
    for (std::size_t node = 0;
         node < assembly.mappings.nodeVertexIds.size(); ++node) {
        if (assembly.mappings.nodeVertexIds[node] != 10) {
            assembly.definition.nodes[node].fixed = true;
        }
    }
    return assembly;
}

struct Fixture {
    Scene scene = openSquarePyramidScene();
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structureAssembly = fixedMouthAssembly(scene);
    Structure structure{structureAssembly.definition};
    SceneFluidSurfaceTransfer transfer{
        surface.definition, structureAssembly.mappings, structure};
};

struct Endpoint {
    SceneFluidSurfaceState state;
    SceneFluidGridEpoch gridEpoch;
    SceneFluidCellVolumeSet volumes;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;
};

Endpoint captureEndpoint(const Fixture& fixture) {
    auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    auto gridEpoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer);
    auto volumes = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, grid(), fixture.transfer,
        gridEpoch);
    auto caps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    auto quadrature = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, state, caps);
    auto patches = buildSceneFluidOpeningGridPatches(
        fixture.surface.definition, state, caps, quadrature, grid());
    return {
        std::move(state), std::move(gridEpoch), std::move(volumes),
        std::move(caps), std::move(quadrature), std::move(patches),
    };
}

fluid::MacVelocityField constantXVelocity(const double value) {
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), value);
    return velocity;
}

SceneFluidOpeningFluxSet evaluateFlux(
    const Fixture& fixture,
    const Endpoint& endpoint,
    const fluid::MacVelocityField& velocity) {
    return evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, endpoint.state, endpoint.caps,
        endpoint.quadrature, endpoint.patches, grid(), velocity);
}

const SceneFluidRegionVolume& volumeFor(
    const SceneFluidCellVolumeSet& volumes,
    const StableId regionId) {
    const auto found = std::ranges::find(
        volumes.regionVolumes, regionId, &SceneFluidRegionVolume::regionId);
    if (found == volumes.regionVolumes.end()) {
        throw std::runtime_error("continuity test volume region is missing");
    }
    return *found;
}

const SceneFluidRegionContinuity& continuityFor(
    const SceneFluidRegionContinuitySet& continuity,
    const StableId regionId) {
    const auto found = std::ranges::find(
        continuity.regions, regionId,
        &SceneFluidRegionContinuity::regionId);
    if (found == continuity.regions.end()) {
        throw std::runtime_error(
            "continuity test audit region is missing");
    }
    return *found;
}

void advanceApex(Fixture& fixture, const double accelerationX) {
    const auto found = std::ranges::find(
        fixture.structureAssembly.mappings.nodeVertexIds,
        StableId{10});
    if (found == fixture.structureAssembly.mappings.nodeVertexIds.end()) {
        throw std::runtime_error("continuity test apex node is missing");
    }
    const std::size_t node = static_cast<std::size_t>(
        found - fixture.structureAssembly.mappings.nodeVertexIds.begin());
    const double mass = fixture.structure.definition().nodes[node].massKg;
    fixture.structure.addExternalForce(
        node, {mass * accelerationX, 0.0, 0.0});
    StructureStepSettings step;
    step.timeStepSeconds = 0.25;
    step.substeps = 1;
    step.constraintIterations = 0;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(step);
    check(diagnostics.finite, "continuity test structural step is finite");
}

struct MovingEpochPair {
    Fixture fixture;
    Endpoint previous;
    Endpoint current;
    fluid::MacVelocityField previousVelocity;
    fluid::MacVelocityField currentVelocity;
    SceneFluidOpeningFluxSet previousFlux;
    SceneFluidOpeningFluxSet currentFlux;

    MovingEpochPair()
        : previous(captureEndpoint(fixture)),
          current([&] {
              advanceApex(fixture, -0.8);
              return captureEndpoint(fixture);
          }()),
          previousVelocity(constantXVelocity(0.0)),
          currentVelocity([&] {
              const double change =
                  volumeFor(current.volumes, 2)
                      .summedCellVolumeCubicMeters
                  - volumeFor(previous.volumes, 2)
                      .summedCellVolumeCubicMeters;
              const double duration = current.state.simulationTimeSeconds
                  - previous.state.simulationTimeSeconds;
              return constantXVelocity(
                  -2.0 * change
                  / (duration * current.caps.totalAreaSquareMeters));
          }()),
          previousFlux(evaluateFlux(
              fixture, previous, previousVelocity)),
          currentFlux(evaluateFlux(
              fixture, current, currentVelocity)) {}
};

void testCompatibleMovingCellContinuity() {
    MovingEpochPair pair;
    const double cellVolumeChange =
        volumeFor(pair.current.volumes, 2).summedCellVolumeCubicMeters
        - volumeFor(pair.previous.volumes, 2).summedCellVolumeCubicMeters;
    check(cellVolumeChange > 1.0e-5,
          "moving continuity fixture has a nonzero expanding cell");
    const auto first = auditSceneFluidRegionContinuity(
        pair.previous.volumes, pair.current.volumes,
        pair.previousFlux, pair.currentFlux);
    const auto repeated = auditSceneFluidRegionContinuity(
        pair.previous.volumes, pair.current.volumes,
        pair.previousFlux, pair.currentFlux);
    check(first == repeated
              && first.version == sceneFluidRegionContinuityVersion
              && first.fingerprint != 0
              && first.regions.size() == 2
              && first.allRegionsWithinTolerance
              && first.failedRegionCount == 0,
          "compatible moving-cell continuity is deterministic and accepted");
    const auto& cell = continuityFor(first, 2);
    const auto& outside = continuityFor(first, 1);
    checkNear(cell.geometryVolumeChangeCubicMeters,
              cellVolumeChange, 0.0,
              "cell audit retains measured geometry volume change");
    checkNear(cell.integratedOutwardRelativeVolumeCubicMeters,
              -cellVolumeChange, 2.0e-14,
              "endpoint-trapezoidal intake transport supplies expanding cell volume");
    checkNear(cell.continuityResidualCubicMeters,
              0.0, 2.0e-14,
              "expanding cell volume and intake transport close");
    checkNear(outside.geometryVolumeChangeCubicMeters,
              -cellVolumeChange, 3.0e-12,
              "Outside geometry change opposes cell expansion");
    checkNear(outside.integratedOutwardRelativeVolumeCubicMeters,
              cellVolumeChange, 2.0e-14,
              "Outside transport opposes cell intake");
    checkNear(first.globalGeometryVolumeChangeCubicMeters,
              0.0, 3.0e-12,
              "all authored geometry volumes conserve the fixed domain");
    checkNear(first.globalIntegratedOutwardRelativeVolumeCubicMeters,
              0.0, 0.0,
              "all authored opening transports cancel globally");
    checkNear(first.globalContinuityResidualCubicMeters,
              0.0, 3.0e-12,
              "global two-epoch continuity closes");
    validateSceneFluidRegionContinuity(
        first, pair.previous.volumes, pair.current.volumes,
        pair.previousFlux, pair.currentFlux);
}

void testStationaryAndMismatchedFlowDiagnostics() {
    Fixture stationaryFixture;
    const auto previous = captureEndpoint(stationaryFixture);
    advanceApex(stationaryFixture, 0.0);
    const auto current = captureEndpoint(stationaryFixture);
    const auto zeroVelocity = constantXVelocity(0.0);
    const auto previousFlux = evaluateFlux(
        stationaryFixture, previous, zeroVelocity);
    const auto currentFlux = evaluateFlux(
        stationaryFixture, current, zeroVelocity);
    const auto stationary = auditSceneFluidRegionContinuity(
        previous.volumes, current.volumes, previousFlux, currentFlux);
    check(stationary.allRegionsWithinTolerance
              && stationary.maximumAbsoluteContinuityResidualCubicMeters
                  == 0.0,
          "stationary consecutive epochs close continuity exactly");

    MovingEpochPair moving;
    const auto missingPreviousFlux = evaluateFlux(
        moving.fixture, moving.previous, zeroVelocity);
    const auto missingCurrentFlux = evaluateFlux(
        moving.fixture, moving.current, zeroVelocity);
    const auto mismatch = auditSceneFluidRegionContinuity(
        moving.previous.volumes, moving.current.volumes,
        missingPreviousFlux, missingCurrentFlux);
    check(!mismatch.allRegionsWithinTolerance
              && mismatch.failedRegionCount == 2
              && mismatch.maximumAbsoluteContinuityResidualCubicMeters
                  > 1.0e-5,
          "missing intake flow is reported as a per-region continuity mismatch");
    checkNear(mismatch.globalContinuityResidualCubicMeters,
              0.0, 3.0e-12,
              "equal-and-opposite local mismatches still cancel globally");
}

void testBindingCorruptionAndLimits() {
    MovingEpochPair pair;
    const auto accepted = auditSceneFluidRegionContinuity(
        pair.previous.volumes, pair.current.volumes,
        pair.previousFlux, pair.currentFlux);

    auto corrupt = accepted;
    corrupt.regions.front().continuityResidualCubicMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidRegionContinuity(
            corrupt, pair.previous.volumes, pair.current.volumes,
            pair.previousFlux, pair.currentFlux); },
        "continuity validation rejects nested ledger corruption");

    auto corruptFlux = pair.currentFlux;
    corruptFlux.regions.front()
        .outwardRelativeVolumeFlowRateCubicMetersPerSecond += 0.01;
    expectInvalid(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, pair.current.volumes,
            pair.previousFlux, corruptFlux)); },
        "continuity audit rejects a corrupted source flux epoch");

    auto corruptVolumes = pair.current.volumes;
    corruptVolumes.regionVolumes.front().summedCellVolumeCubicMeters += 0.01;
    expectInvalid(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, corruptVolumes,
            pair.previousFlux, pair.currentFlux)); },
        "continuity audit rejects a corrupted source volume epoch");

    SceneFluidRegionContinuityLimits limits;
    limits.maximumRegions = 1;
    expectLimited(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, pair.current.volumes,
            pair.previousFlux, pair.currentFlux, {}, limits)); },
        "continuity audit bounds region count");
    limits = {};
    limits.maximumContinuityBytes = 0;
    expectLimited(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, pair.current.volumes,
            pair.previousFlux, pair.currentFlux, {}, limits)); },
        "continuity audit bounds owned storage");

    SceneFluidRegionContinuitySettings invalidSettings;
    invalidSettings.absoluteVolumeToleranceCubicMeters = 0.0;
    invalidSettings.relativeVolumeTolerance = 0.0;
    expectInvalid(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, pair.current.volumes,
            pair.previousFlux, pair.currentFlux, invalidSettings)); },
        "continuity audit rejects zero tolerance authority");

    advanceApex(pair.fixture, 0.0);
    const auto skippedCurrent = captureEndpoint(pair.fixture);
    const auto skippedFlux = evaluateFlux(
        pair.fixture, skippedCurrent, constantXVelocity(0.0));
    expectInvalid(
        [&] { static_cast<void>(auditSceneFluidRegionContinuity(
            pair.previous.volumes, skippedCurrent.volumes,
            pair.previousFlux, skippedFlux)); },
        "continuity audit rejects skipped accepted epochs");
}

} // namespace

int main() {
    try {
        testCompatibleMovingCellContinuity();
        testStationaryAndMismatchedFlowDiagnostics();
        testBindingCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid region-continuity check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid region-continuity checks passed");
    return 0;
}
