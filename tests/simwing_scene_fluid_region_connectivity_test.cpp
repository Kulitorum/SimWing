#include "scene_fluid_region_connectivity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <ranges>
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

template<std::size_t VertexCount>
void addTriangle(Scene& scene,
                 const StableId id,
                 const std::array<std::size_t, 3>& vertices,
                 const std::array<Vec3, VertexCount>& positions,
                 const StableId negativeRegion,
                 const StableId positiveRegion,
                 const StableId sheet) {
    scene.triangles.push_back({
        id,
        {10 + vertices[0], 10 + vertices[1], 10 + vertices[2]},
        intrinsicChart(positions, vertices),
        negativeRegion, positiveRegion, 100, sheet, SurfaceRole::Skin,
    });
}

Scene crossportBipyramidScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-crossport-connectivity";
    scene.metadata.exporterVersion =
        "scene-fluid-region-connectivity-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "left-cell"},
        {3, RegionKind::Cell, "right-cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 6> positions{{
        {1.0, 1.0, 1.0},
        {2.0, 0.5, 0.5},
        {2.0, 1.5, 0.5},
        {2.0, 1.5, 1.5},
        {2.0, 0.5, 1.5},
        {3.0, 1.0, 1.0},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 4> leftFaces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
    }};
    const std::array<std::array<std::size_t, 3>, 4> rightFaces{{
        {{5, 1, 2}}, {{5, 2, 3}}, {{5, 3, 4}}, {{5, 4, 1}},
    }};
    for (std::size_t face = 0; face < leftFaces.size(); ++face) {
        addTriangle(
            scene, 500 + face, leftFaces[face], positions, 2, 1, 900);
        addTriangle(
            scene, 510 + face, rightFaces[face], positions, 3, 1, 901);
    }
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 3, OpeningRole::Crossport},
    };
    return scene;
}

Scene intakePyramidScene() {
    Scene scene = crossportBipyramidScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-intake-connectivity";
    scene.vertices.resize(5);
    std::erase_if(
        scene.triangles,
        [](const Triangle& triangle) { return triangle.id >= 510; });
    scene.regions.resize(2);
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

Scene separateSealedPyramidsScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-separated-sealed-connectivity";
    scene.metadata.exporterVersion =
        "scene-fluid-region-connectivity-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "left-cell"},
        {3, RegionKind::Cell, "right-cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 10> positions{{
        {0.5, 1.0, 1.0},
        {1.5, 0.5, 0.5},
        {1.5, 1.5, 0.5},
        {1.5, 1.5, 1.5},
        {1.5, 0.5, 1.5},
        {3.5, 1.0, 1.0},
        {2.5, 0.5, 0.5},
        {2.5, 1.5, 0.5},
        {2.5, 1.5, 1.5},
        {2.5, 0.5, 1.5},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 6> leftFaces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
        {{1, 2, 3}}, {{1, 3, 4}},
    }};
    const std::array<std::array<std::size_t, 3>, 6> rightFaces{{
        {{5, 6, 7}}, {{5, 7, 8}}, {{5, 8, 9}}, {{5, 9, 6}},
        {{6, 8, 7}}, {{6, 9, 8}},
    }};
    for (std::size_t face = 0; face < leftFaces.size(); ++face) {
        addTriangle(
            scene, 500 + face, leftFaces[face], positions, 2, 1, 900);
        addTriangle(
            scene, 510 + face, rightFaces[face], positions, 3, 1, 901);
    }
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {0.0, 0.0, 0.0}, {4.0, 4.0, 4.0}};
}

SceneStructureAssembly assemblyWithFixedVertices(
    const Scene& scene,
    const std::set<StableId>& fixedVertices) {
    auto assembly = assembleSceneStructure(scene);
    for (std::size_t node = 0;
         node < assembly.mappings.nodeVertexIds.size(); ++node) {
        assembly.definition.nodes[node].fixed = fixedVertices.contains(
            assembly.mappings.nodeVertexIds[node]);
    }
    return assembly;
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;

    explicit Fixture(Scene source,
                     const std::set<StableId>& fixedVertices = {})
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(
              assemblyWithFixedVertices(scene, fixedVertices)),
          structure(structureAssembly.definition),
          transfer(surface.definition,
                   structureAssembly.mappings,
                   structure) {}
};

struct Endpoint {
    SceneFluidSurfaceState state;
    SceneFluidGridEpoch gridEpoch;
    SceneFluidCellVolumeSet volumes;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;
    SceneFluidOpeningFluxSet flux;
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
    fluid::MacVelocityField zeroVelocity(grid());
    auto flux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, state, caps, quadrature, patches,
        grid(), zeroVelocity);
    return {
        std::move(state), std::move(gridEpoch), std::move(volumes),
        std::move(caps), std::move(quadrature), std::move(patches),
        std::move(flux),
    };
}

void advance(Fixture& fixture,
             const std::set<StableId>& acceleratedVertices,
             const double accelerationX) {
    for (std::size_t node = 0;
         node < fixture.structureAssembly.mappings.nodeVertexIds.size();
         ++node) {
        if (!acceleratedVertices.contains(
                fixture.structureAssembly.mappings.nodeVertexIds[node])) {
            continue;
        }
        const double mass = fixture.structure.definition().nodes[node].massKg;
        fixture.structure.addExternalForce(
            node, {mass * accelerationX, 0.0, 0.0});
    }
    StructureStepSettings step;
    step.timeStepSeconds = 0.25;
    step.substeps = 1;
    step.constraintIterations = 0;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(step);
    check(diagnostics.finite, "connectivity structural step is finite");
}

const SceneFluidRegionComponentContinuity& componentForGauge(
    const SceneFluidRegionComponentContinuitySet& components,
    const StableId gaugeRegionId) {
    const auto found = std::ranges::find(
        components.components, gaugeRegionId,
        &SceneFluidRegionComponentContinuity::gaugeRegionId);
    if (found == components.components.end()) {
        throw std::runtime_error(
            "connectivity test gauge component is missing");
    }
    return *found;
}

void testCanonicalIntakeAndCrossportConnectivity() {
    const auto intakeSurface = assembleSceneFluidSurface(
        intakePyramidScene());
    check(intakeSurface.ok(), "intake connectivity fixture assembles");
    const auto intake = buildSceneFluidRegionConnectivity(
        intakeSurface.definition);
    const auto intakeRepeated = buildSceneFluidRegionConnectivity(
        intakeSurface.definition);
    check(intake == intakeRepeated
              && intake.version == sceneFluidRegionConnectivityVersion
              && intake.fingerprint != 0
              && intake.components.size() == 1
              && intake.regions.size() == 2
              && intake.openings.size() == 1,
          "intake connectivity is deterministic and joins both regions");
    check(intake.components.front().gaugeRegionId == 1
              && intake.components.front().containsOutside
              && intake.components.front().regionCount == 2
              && intake.components.front().openingCount == 1
              && intake.openings.front().role == OpeningRole::Intake,
          "intake component owns the canonical Outside pressure gauge");
    validateSceneFluidRegionConnectivity(intake, intakeSurface.definition);

    const auto crossportSurface = assembleSceneFluidSurface(
        crossportBipyramidScene());
    check(crossportSurface.ok(), "crossport connectivity fixture assembles");
    const auto crossport = buildSceneFluidRegionConnectivity(
        crossportSurface.definition);
    check(crossport.components.size() == 2
              && crossport.components[0].gaugeRegionId == 1
              && crossport.components[0].containsOutside
              && crossport.components[0].regionCount == 1
              && crossport.components[0].openingCount == 0
              && crossport.components[1].gaugeRegionId == 2
              && !crossport.components[1].containsOutside
              && crossport.components[1].regionCount == 2
              && crossport.components[1].openingCount == 1,
          "crossport joins its cells but leaves Outside pressure-disconnected");
    check(crossport.regions[1].componentIndex == 1
              && crossport.regions[2].componentIndex == 1
              && crossport.openings.front().componentIndex == 1
              && crossport.openings.front().role == OpeningRole::Crossport,
          "crossport members retain canonical component identity");
}

void testConnectedAndSealedComponentCompatibility() {
    Fixture connected(intakePyramidScene());
    check(connected.surface.ok() && connected.structureAssembly.ok(),
          "connected component-continuity fixture assembles");
    const auto connectedPrevious = captureEndpoint(connected);
    advance(connected, {}, 0.0);
    const auto connectedCurrent = captureEndpoint(connected);
    const auto connectedContinuity = auditSceneFluidRegionContinuity(
        connectedPrevious.volumes, connectedCurrent.volumes,
        connectedPrevious.flux, connectedCurrent.flux);
    const auto connectivity = buildSceneFluidRegionConnectivity(
        connected.surface.definition);
    const auto compatible = auditSceneFluidRegionComponentContinuity(
        connectivity, connectedContinuity, connected.surface.definition,
        connectedPrevious.volumes, connectedCurrent.volumes,
        connectedPrevious.flux, connectedCurrent.flux);
    check(compatible.components.size() == 1
              && compatible.allComponentsCompatible
              && compatible.failedComponentCount == 0
              && compatible.fingerprint != 0,
          "stationary intake component is pressure compatible");
    checkNear(
        componentForGauge(compatible, 1)
            .previousOpeningSourceResidualCubicMetersPerSecond,
        0.0, 0.0,
        "intake endpoint source cancels within its connected component");
    validateSceneFluidRegionComponentContinuity(
        compatible, connectivity, connectedContinuity,
        connected.surface.definition,
        connectedPrevious.volumes, connectedCurrent.volumes,
        connectedPrevious.flux, connectedCurrent.flux);

    Fixture sealed(
        separateSealedPyramidsScene(),
        {StableId{11}, StableId{12}, StableId{13}, StableId{14},
         StableId{16}, StableId{17}, StableId{18}, StableId{19}});
    check(sealed.surface.ok() && sealed.structureAssembly.ok(),
          "sealed component-continuity fixture assembles");
    const auto sealedPrevious = captureEndpoint(sealed);
    advance(
        sealed, {StableId{10}, StableId{15}}, -0.4);
    const auto sealedCurrent = captureEndpoint(sealed);
    const auto sealedContinuity = auditSceneFluidRegionContinuity(
        sealedPrevious.volumes, sealedCurrent.volumes,
        sealedPrevious.flux, sealedCurrent.flux);
    const auto sealedConnectivity = buildSceneFluidRegionConnectivity(
        sealed.surface.definition);
    const auto incompatible = auditSceneFluidRegionComponentContinuity(
        sealedConnectivity, sealedContinuity, sealed.surface.definition,
        sealedPrevious.volumes, sealedCurrent.volumes,
        sealedPrevious.flux, sealedCurrent.flux);
    check(sealedConnectivity.components.size() == 3
              && sealedConnectivity.openings.empty()
              && !incompatible.allComponentsCompatible
              && incompatible.failedComponentCount == 2,
          "moving sealed divider exposes two incompatible pressure components");
    const auto& outside = componentForGauge(incompatible, 1);
    const auto& left = componentForGauge(incompatible, 2);
    const auto& right = componentForGauge(incompatible, 3);
    check(outside.compatible
              && !left.compatible
              && !right.compatible
              && std::abs(left.continuityResidualCubicMeters) > 1.0e-5,
          "only the two changing sealed cells fail component continuity");
    checkNear(left.continuityResidualCubicMeters,
              -right.continuityResidualCubicMeters, 3.0e-12,
              "sealed-cell incompatibilities are equal and opposite");
    checkNear(incompatible.globalContinuityResidualCubicMeters,
              0.0, 3.0e-12,
              "sealed local incompatibility still preserves global volume");

    auto corruptComponents = incompatible;
    corruptComponents.components[1].continuityResidualCubicMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidRegionComponentContinuity(
            corruptComponents, sealedConnectivity, sealedContinuity,
            sealed.surface.definition,
            sealedPrevious.volumes, sealedCurrent.volumes,
            sealedPrevious.flux, sealedCurrent.flux); },
        "component-continuity validation rejects nested corruption");

    SceneFluidRegionConnectivityLimits componentLimits;
    componentLimits.maximumComponents = 2;
    expectLimited(
        [&] { static_cast<void>(auditSceneFluidRegionComponentContinuity(
            sealedConnectivity, sealedContinuity, sealed.surface.definition,
            sealedPrevious.volumes, sealedCurrent.volumes,
            sealedPrevious.flux, sealedCurrent.flux, componentLimits)); },
        "component continuity bounds pressure-component count");
    componentLimits = {};
    componentLimits.maximumConnectivityBytes = 0;
    expectLimited(
        [&] { static_cast<void>(auditSceneFluidRegionComponentContinuity(
            sealedConnectivity, sealedContinuity, sealed.surface.definition,
            sealedPrevious.volumes, sealedCurrent.volumes,
            sealedPrevious.flux, sealedCurrent.flux, componentLimits)); },
        "component continuity bounds owned storage");
}

void testCorruptionAndLimits() {
    const auto surface = assembleSceneFluidSurface(
        crossportBipyramidScene());
    const auto accepted = buildSceneFluidRegionConnectivity(
        surface.definition);
    auto corrupt = accepted;
    corrupt.components.front().gaugeRegionId = 99;
    expectInvalid(
        [&] { validateSceneFluidRegionConnectivity(
            corrupt, surface.definition); },
        "region-connectivity validation rejects gauge corruption");

    SceneFluidRegionConnectivityLimits limits;
    limits.maximumRegions = 2;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidRegionConnectivity(
            surface.definition, limits)); },
        "region connectivity bounds authored regions");
    limits = {};
    limits.maximumOpenings = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidRegionConnectivity(
            surface.definition, limits)); },
        "region connectivity bounds authored openings");
    limits = {};
    limits.maximumComponents = 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidRegionConnectivity(
            surface.definition, limits)); },
        "region connectivity bounds connected components");
    limits = {};
    limits.maximumConnectivityBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidRegionConnectivity(
            surface.definition, limits)); },
        "region connectivity bounds owned storage");
}

} // namespace

int main() {
    try {
        testCanonicalIntakeAndCrossportConnectivity();
        testConnectedAndSealedComponentCompatibility();
        testCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid region-connectivity check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid region-connectivity checks passed");
    return 0;
}
