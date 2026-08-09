#include "scene_fluid_surface.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

Scene surfaceScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-surface";
    scene.metadata.exporterVersion = "scene-fluid-surface-test/1";
    scene.regions = {
        {2, RegionKind::Cell, "cell"},
        {1, RegionKind::Outside, "outside"},
    };
    scene.vertices = {
        {12, {1.0, 1.0, 0.0}},
        {10, {0.0, 0.0, 0.0}},
        {13, {0.0, 1.0, 0.0}},
        {11, {1.0, 0.0, 0.0}},
    };
    scene.fabricMaterials = {
        {100, "porous-ripstop", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.triangles = {
        {501, {10, 12, 13}, {{{0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {500, {10, 11, 12}, {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    scene.openings = {
        {600, {10, 11, 12, 13}, 1, 2, OpeningRole::Intake,
         {{{10, 11, 12}}, {{10, 12, 13}}}},
    };
    return scene;
}

void reverseRelevantCollections(Scene& scene) {
    std::reverse(scene.regions.begin(), scene.regions.end());
    std::reverse(scene.vertices.begin(), scene.vertices.end());
    std::reverse(scene.fabricMaterials.begin(), scene.fabricMaterials.end());
    std::reverse(scene.triangles.begin(), scene.triangles.end());
    std::reverse(scene.openings.begin(), scene.openings.end());
}

void testDeterministicSurfaceAssembly() {
    const Scene scene = surfaceScene();
    const auto first = assembleSceneFluidSurface(scene);
    Scene reordered = scene;
    reverseRelevantCollections(reordered);
    const auto second = assembleSceneFluidSurface(reordered);

    check(first.ok() && second.ok()
              && first.definition == second.definition
              && first.definition.fingerprint != 0,
          "scene fluid surface: validated input assembles deterministically");
    const auto& definition = first.definition;
    check(definition.version == sceneFluidSurfaceDefinitionVersion
              && definition.regions.size() == 2
              && definition.materials.size() == 1
              && definition.vertices.size() == 4
              && definition.triangles.size() == 2
              && definition.openings.size() == 1,
          "scene fluid surface: compact production counts are complete");
    check(definition.mappings.regionIds == std::vector<StableId>({1, 2})
              && definition.mappings.materialIds
                  == std::vector<StableId>({100})
              && definition.mappings.vertexIds
                  == std::vector<StableId>({10, 11, 12, 13})
              && definition.mappings.triangleIds
                  == std::vector<StableId>({500, 501})
              && definition.mappings.openingIds
                  == std::vector<StableId>({600}),
          "scene fluid surface: every entity uses canonical stable-ID order");
    check(definition.materials.front().porosityFraction == 0.0125
              && definition.materials.front().permeabilitySquareMeters
                  == 2.5e-12
              && definition.triangles.front().vertexIndices
                  == std::array<std::size_t, 3>({0, 1, 2})
              && definition.triangles.front().negativeSideRegionIndex == 0
              && definition.triangles.front().positiveSideRegionIndex == 1
              && definition.openings.front().orderedVertexIndices
                  == std::vector<std::size_t>({0, 1, 2, 3})
              && definition.openings.front().capTriangleVertexIndices
                  == std::vector<std::array<std::size_t, 3>>({
                      {{0, 1, 2}}, {{0, 2, 3}}}),
          "scene fluid surface: winding, side ownership, permeability, and opening order survive");
    check(definition.mappings.vertexIndex(12) == 2
              && definition.mappings.triangleIndex(501) == 1
              && !definition.mappings.openingIndex(999),
          "scene fluid surface: deterministic reverse mappings are bounded");

    Scene internal = scene;
    internal.triangles.front().role = SurfaceRole::Diagonal;
    internal.triangles.front().positiveSideRegionId =
        internal.triangles.front().negativeSideRegionId;
    const auto internalAssembly = assembleSceneFluidSurface(internal);
    check(internalAssembly.ok()
              && internalAssembly.definition.triangles[1]
                     .negativeSideRegionIndex
                  == internalAssembly.definition.triangles[1]
                         .positiveSideRegionIndex
              && internalAssembly.definition.triangles[1].role
                  == SurfaceRole::Diagonal,
          "scene fluid surface: authored same-region internal sheets remain explicit");
}

void testAcceptedStructureCapture() {
    const Scene scene = surfaceScene();
    const auto fluidAssembly = assembleSceneFluidSurface(scene);
    const auto structureAssembly = assembleSceneStructure(scene);
    check(fluidAssembly.ok() && structureAssembly.ok(),
          "scene fluid state: both authoritative adapters accept the scene");
    Structure structure(structureAssembly.definition);

    const auto initial = captureSceneFluidSurfaceState(
        fluidAssembly.definition, structureAssembly.mappings, structure);
    check(initial.version == sceneFluidSurfaceStateVersion
              && initial.fingerprint != 0
              && initial.definitionFingerprint
                  == fluidAssembly.definition.fingerprint
              && initial.structureDefinitionFingerprint
                  == structure.checkpoint().definitionFingerprint
              && initial.acceptedStepCount == 0
              && initial.simulationTimeSeconds == 0.0
              && initial.vertices.size() == 4,
          "scene fluid state: initial accepted epoch binds both definitions");
    validateSceneFluidSurfaceState(fluidAssembly.definition, initial);

    for (std::size_t index = 0;
         index < structureAssembly.definition.nodes.size(); ++index) {
        structure.addExternalForce(
            index,
            {structureAssembly.definition.nodes[index].massKg, 0.0, 0.0});
    }
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.01;
    settings.gravityMetersPerSecondSquared = {};
    settings.constraintIterations = 4;
    const auto diagnostics = structure.step(settings);
    const auto captured = captureSceneFluidSurfaceState(
        fluidAssembly.definition, structureAssembly.mappings, structure);
    const auto replay = captureSceneFluidSurfaceState(
        fluidAssembly.definition, structureAssembly.mappings, structure);
    const auto nodeStates = structure.nodeStates();
    bool matchesNodes = true;
    for (std::size_t index = 0; index < captured.vertices.size(); ++index) {
        const auto nodeIndex = structureAssembly.mappings.nodeIndex(
            captured.vertices[index].id);
        if (!nodeIndex) {
            matchesNodes = false;
            continue;
        }
        const auto& node = nodeStates[*nodeIndex];
        matchesNodes = matchesNodes
            && captured.vertices[index].positionMeters.x
                == node.positionMeters.x
            && captured.vertices[index].positionMeters.y
                == node.positionMeters.y
            && captured.vertices[index].positionMeters.z
                == node.positionMeters.z
            && captured.vertices[index].velocityMetersPerSecond.x
                == node.velocityMetersPerSecond.x
            && captured.vertices[index].velocityMetersPerSecond.y
                == node.velocityMetersPerSecond.y
            && captured.vertices[index].velocityMetersPerSecond.z
                == node.velocityMetersPerSecond.z;
    }
    check(diagnostics.finite,
          "scene fluid state: fixture Structure step remains finite");
    check(captured == replay,
          "scene fluid state: repeated accepted-state capture is exact");
    check(captured.acceptedStepCount == 1
              && captured.simulationTimeSeconds == 0.01,
          "scene fluid state: accepted epoch follows Structure time");
    check(matchesNodes,
          "scene fluid state: accepted Structure motion is captured exactly in stable order");
    auto corruptState = captured;
    corruptState.vertices.front().positionMeters.x += 0.01;
    expectRejected(
        [&] { validateSceneFluidSurfaceState(
            fluidAssembly.definition, corruptState); },
        "scene fluid state: payload corruption invalidates its epoch fingerprint");

    auto corruptDefinition = fluidAssembly.definition;
    ++corruptDefinition.fingerprint;
    expectRejected(
        [&] { static_cast<void>(captureSceneFluidSurfaceState(
            corruptDefinition, structureAssembly.mappings, structure)); },
        "scene fluid state: corrupt surface identity is rejected");
    auto corruptMappings = structureAssembly.mappings;
    corruptMappings.nodeVertexIds.erase(
        corruptMappings.nodeVertexIds.begin());
    expectRejected(
        [&] { static_cast<void>(captureSceneFluidSurfaceState(
            fluidAssembly.definition, corruptMappings, structure)); },
        "scene fluid state: incomplete Structure mapping is rejected");
    corruptMappings = structureAssembly.mappings;
    std::reverse(corruptMappings.nodeVertexIds.begin(),
                 corruptMappings.nodeVertexIds.end());
    expectRejected(
        [&] { static_cast<void>(captureSceneFluidSurfaceState(
            fluidAssembly.definition, corruptMappings, structure)); },
        "scene fluid state: reordered Structure mapping is rejected");

    Scene foreignScene = surfaceScene();
    foreignScene.vertices.front().positionMeters.x += 0.125;
    const auto foreignFluid = assembleSceneFluidSurface(foreignScene);
    check(foreignFluid.ok(),
          "scene fluid state: foreign geometry fixture remains valid");
    expectRejected(
        [&] { static_cast<void>(captureSceneFluidSurfaceState(
            foreignFluid.definition,
            structureAssembly.mappings,
            structure)); },
        "scene fluid state: foreign reference geometry is rejected");
}

void testTransactionalRejectionAndLimits() {
    Scene invalid = surfaceScene();
    invalid.triangles.front().materialId = 999;
    const auto invalidAssembly = assembleSceneFluidSurface(invalid);
    check(!invalidAssembly.ok()
              && !invalidAssembly.assembled
              && invalidAssembly.definition.triangles.empty()
              && !invalidAssembly.diagnostics.empty()
              && invalidAssembly.diagnostics.front().code
                  == SceneFluidSurfaceDiagnosticCode::InvalidScene
              && invalidAssembly.diagnostics.front().sceneValidationCode
                  == ValidationCode::MissingMaterialReference,
          "scene fluid surface: scene validation failure publishes no production data");

    SceneFluidSurfaceLimits limits;
    limits.maximumTriangles = 1;
    const auto limited = assembleSceneFluidSurface(surfaceScene(), limits);
    check(!limited.ok() && limited.definition.vertices.empty()
              && limited.diagnostics.front().code
                  == SceneFluidSurfaceDiagnosticCode::LimitExceeded,
          "scene fluid surface: entity limits reject transactionally");
    limits = {};
    limits.maximumOpeningCapTriangles = 1;
    const auto capLimited =
        assembleSceneFluidSurface(surfaceScene(), limits);
    check(!capLimited.ok()
              && capLimited.definition.vertices.empty()
              && capLimited.diagnostics.front().code
                  == SceneFluidSurfaceDiagnosticCode::LimitExceeded,
          "scene fluid surface: authored cap triangles are bounded");
    limits = {};
    limits.maximumMappingBytes = 1;
    const auto mappingLimited =
        assembleSceneFluidSurface(surfaceScene(), limits);
    check(!mappingLimited.ok()
              && mappingLimited.definition.vertices.empty()
              && mappingLimited.diagnostics.front().code
                  == SceneFluidSurfaceDiagnosticCode::MappingOverflow,
          "scene fluid surface: mapping byte limit rejects before publication");

    Scene unsupportedOpening = surfaceScene();
    unsupportedOpening.vertices.push_back({14, {0.5, 0.5, 0.0}});
    unsupportedOpening.openings.front().orderedVertexIds = {10, 11, 14};
    unsupportedOpening.openings.front().capTriangleVertexIds.clear();
    check(validateScene(unsupportedOpening).ok(),
          "scene fluid surface: opening-only source geometry remains scene-valid");
    const auto unsupported =
        assembleSceneFluidSurface(unsupportedOpening);
    check(!unsupported.ok() && unsupported.definition.openings.empty()
              && unsupported.diagnostics.front().code
                  == SceneFluidSurfaceDiagnosticCode::
                      OpeningVertexWithoutFabricMotion,
          "scene fluid surface: opening-only motion is rejected instead of frozen");
}

} // namespace

int main() {
    testDeterministicSurfaceAssembly();
    testAcceptedStructureCapture();
    testTransactionalRejectionAndLimits();
    if (failures != 0) {
        std::fprintf(stderr, "%d scene fluid-surface test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("scene fluid-surface tests passed\n");
    return 0;
}
