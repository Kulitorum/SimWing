#include "scene_fluid_mimetic_pressure_source.h"
#include "scene_fluid_mimetic_trace_flow.h"
#include "scene_fluid_pressure_projection.h"
#include "scene_fluid_pressure_link_flow.h"
#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_region_momentum.h"
#include "scene_fluid_region_link_flow.h"
#include "scene_fluid_region_transport.h"
#include "scene_fluid_region_wall.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <ranges>
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
        "sha256:scene-fluid-pressure-projection";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-projection-test/1";
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
        "sha256:scene-fluid-pressure-projection-open";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-projection-test/1";
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
        "sha256:scene-fluid-pressure-projection-tilted-open";
    scene.vertices[1].positionMeters.x = 2.6;
    scene.vertices[2].positionMeters.x = 2.8;
    scene.vertices[3].positionMeters.x = 2.7;
    return scene;
}

Scene concaveOpenScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-projection-concave-open";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-projection-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.vertices = {
        {10, {1.2, 1.5, 1.45}},
        {11, {2.4, 1.2, 1.2}},
        {12, {2.4, 1.8, 1.2}},
        {13, {2.4, 1.5, 1.38}},
        {14, {2.4, 1.2, 1.8}},
    };
    const std::array<Vec2, 3> chart{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}};
    const std::array<std::array<StableId, 3>, 4> faces{{
        {{10, 12, 11}}, {{10, 13, 12}},
        {{10, 14, 13}}, {{10, 11, 14}},
    }};
    for (std::size_t index = 0; index < faces.size(); ++index) {
        scene.triangles.push_back({
            500 + index, faces[index], chart,
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

Scene authoredNonPlanarOpenScene() {
    auto scene = concaveOpenScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-projection-authored-nonplanar";
    scene.vertices[3].positionMeters.x = 2.5;
    scene.openings.front().capTriangleVertexIds = {
        {{11, 12, 13}}, {{11, 13, 14}},
    };
    return scene;
}

Scene reversedOpenScene() {
    auto scene = openScene();
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-projection-open-reversed";
    scene.vertices.front().positionMeters.x = 2.8;
    for (auto& triangle : scene.triangles) {
        std::swap(triangle.vertexIds[1], triangle.vertexIds[2]);
    }
    std::swap(scene.openings.front().orderedVertexIds[1],
              scene.openings.front().orderedVertexIds[2]);
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

SceneStructureAssembly structureAssembly(
    const Scene& scene, const bool fixed, const bool fixedMouth) {
    auto result = assembleSceneStructure(scene);
    for (std::size_t index = 0;
         index < result.definition.nodes.size(); ++index) {
        if (fixed
            || (fixedMouth
                && result.mappings.nodeVertexIds[index] != StableId{10})) {
            result.definition.nodes[index].fixed = true;
        }
    }
    return result;
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

    explicit Fixture(Scene source,
                     const bool fixed = false,
                     const bool fixedMouth = false)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(::structureAssembly(
              scene, fixed, fixedMouth)),
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

    SceneFluidOpeningFluxSet flux(
        const fluid::MacVelocityField& velocity) const {
        return evaluateSceneFluidOpeningFlux(
            surface.definition, state, caps, openingQuadrature,
            openingPatches, grid(), velocity);
    }

    SceneFluidPressureProjection project(
        const fluid::MacVelocityField& velocity,
        const SceneFluidOpeningFluxSet& openingFlux,
        const std::vector<double>& warmPressure,
        const SceneFluidPressureProjectionSettings& settings = {},
        const SceneFluidPressureProjectionLimits& limits = {}) const {
        return projectSceneFluidPressureLinkFlows(
            surface.definition, state, grid(), transfer, epoch, caps,
            openingQuadrature, openingPatches, openingFlux, velocity,
            volumes, connectivity, pressureVolumes, faceLinks,
            pressureOperator, warmPressure, settings, limits);
    }
};

fluid::MacVelocityField manufacturedVelocity() {
    fluid::MacVelocityField velocity(grid());
    for (std::size_t index = 0; index < grid().cellCount(); ++index) {
        const double sample = static_cast<double>(index + 1);
        velocity.xFaces()[index] = 0.4 * std::sin(0.17 * sample);
        velocity.yFaces()[index] = -0.3 * std::cos(0.11 * sample);
        velocity.zFaces()[index] = 0.2 * std::sin(0.07 * sample + 0.2);
    }
    return velocity;
}

SceneFluidPressureProjectionSettings strictSettings() {
    SceneFluidPressureProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.225;
    settings.timeStepSeconds = 0.02;
    settings.absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond =
        2.0e-11;
    settings.relativeCorrectedVolumeRateTolerance = 1.0e-11;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-13;
    settings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-11;
    settings.pressureSolve.maximumIterations = 4000;
    return settings;
}

double faceValue(const SceneFluidPressureFace& face,
                 const fluid::MacVelocityField& velocity) {
    const auto index = grid().cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case fluid::GridFaceAxis::X:
        return velocity.xFaces()[index];
    case fluid::GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case fluid::GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    throw std::runtime_error("invalid pressure-projection test face axis");
}

void testClosedLinkProjection() {
    Fixture fixture(nestedScene());
    const auto velocity = manufacturedVelocity();
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 0.0);
    const auto settings = strictSettings();
    const auto first = fixture.project(
        velocity, openingFlux, warm, settings);
    const auto repeated = fixture.project(
        velocity, openingFlux, warm, settings);
    check(first == repeated
              && first.version == sceneFluidPressureProjectionVersion
              && first.fingerprint != 0
              && first.diagnostics.accepted
              && first.diagnostics.finite
              && first.diagnostics.pressureSolve.compatible
              && first.diagnostics.pressureSolve.converged
              && first.diagnostics.authoredOpeningLinkCount == 0
              && first.pressurePascals.size()
                  == fixture.pressureVolumes.controlVolumes.size(),
          "closed link projection is deterministic, finite, and accepted");
    check(first.diagnostics
              .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond
              > 0.01
              && first.diagnostics
                  .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  < 2.0e-11,
          "pressure correction removes manufactured control-volume divergence");
    check(first.diagnostics
              .maximumPredictedComponentBalanceResidualCubicMetersPerSecond
              < 2.0e-15
              && first.diagnostics
                  .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond
                  < 2.0e-15,
          "link incidence preserves every closed component balance");

    for (const auto& projected : first.links) {
        const auto& source = fixture.faceLinks.links[projected.linkIndex];
        const auto& face = fixture.faceLinks.faces[source.faceIndex];
        checkNear(
            projected.predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            source.areaSquareMeters * faceValue(face, velocity),
            0.0,
            "same-region projection reads the exact owning MAC face");
        const double expectedCorrection = settings.timeStepSeconds
            / settings.densityKgPerCubicMeter
            * source.geometryWeightMeters
            * (first.pressurePascals[source.minusControlVolumeIndex]
               - first.pressurePascals[source.plusControlVolumeIndex]);
        checkNear(
            projected.pressureCorrectionVolumeFlowRateCubicMetersPerSecond,
            expectedCorrection, 0.0,
            "projected link retains its exact pressure-flow correction");
        checkNear(
            projected.correctedRelativeVolumeFlowRateCubicMetersPerSecond,
            projected.predictedRelativeVolumeFlowRateCubicMetersPerSecond
                + expectedCorrection,
            0.0,
            "projected link closes predicted and pressure-correction flow");
    }
    validateSceneFluidPressureProjectionIntegrity(first);
}

void testFaceAlignedOpeningProjection() {
    Fixture fixture(openScene());
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 2.0);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 17.0);
    const auto projected = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    check(projected.diagnostics.accepted
              && projected.diagnostics.authoredOpeningLinkCount == 1
              && openingFlux.samples.size() == 1,
          "face-aligned intake participates in an accepted link projection");
    const auto aperture = std::ranges::find(
        projected.links,
        SceneFluidPressureFaceLinkKind::AuthoredOpening,
        &SceneFluidPressureProjectedLink::kind);
    check(aperture != projected.links.end()
              && aperture->openingPatchStableId
                  == openingFlux.samples.front().patchStableId,
          "aperture projection binds the exact opening-patch identity");
    if (aperture != projected.links.end()) {
        checkNear(
            aperture->predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            openingFlux.samples.front()
                .relativeVolumeFlowRateCubicMetersPerSecond,
            0.0,
            "aperture prediction reuses accepted relative opening flow");
        checkNear(
            aperture->predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            0.36, 5.0e-15,
            "0.18-square-metre intake carries its exact two-metre-per-second flow");
    }
    const auto complement = aperture == projected.links.end()
        ? projected.links.end()
        : std::ranges::find_if(
            projected.links,
            [&](const auto& link) {
                return link.faceIndex == aperture->faceIndex
                    && link.kind
                        == SceneFluidPressureFaceLinkKind::SameRegion;
            });
    check(complement != projected.links.end(),
          "opening face retains a separately projected complement");
    if (complement != projected.links.end()) {
        checkNear(
            complement->predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            1.64, 5.0e-15,
            "opening complement retains its separate MAC volume flow");
    }
    check(projected.diagnostics
              .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
              < 2.0e-11,
          "cross-region intake projection closes control-volume continuity");

    Fixture reversedFixture(reversedOpenScene());
    const auto reversedFlux = reversedFixture.flux(velocity);
    std::vector<double> reversedWarm(
        reversedFixture.pressureOperator.rows.size(), 0.0);
    const auto reversed = reversedFixture.project(
        velocity, reversedFlux, reversedWarm, strictSettings());
    const auto reversedAperture = std::ranges::find(
        reversed.links,
        SceneFluidPressureFaceLinkKind::AuthoredOpening,
        &SceneFluidPressureProjectedLink::kind);
    check(reversed.diagnostics.accepted
              && reversedAperture != reversed.links.end()
              && reversedFlux.samples.size() == 1
              && reversedFixture.faceLinks.links[
                     reversedAperture->linkIndex].minusRegionId == 1
              && reversedFixture.faceLinks.links[
                     reversedAperture->linkIndex].plusRegionId == 2,
          "reversed intake resolves from spatial Outside to Cell");
    if (reversedAperture != reversed.links.end()) {
        checkNear(
            reversedFlux.samples.front()
                .relativeVolumeFlowRateCubicMetersPerSecond,
            -0.36, 5.0e-15,
            "reversed intake retains negative authored relative flow");
        checkNear(
            reversedAperture
                ->predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            0.36, 5.0e-15,
            "projection flips authored flow into positive spatial face orientation");
    }
}

void testEmbeddedOpeningProjection() {
    Fixture fixture(tiltedOpenScene(), true);
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 2.0);
    std::ranges::fill(velocity.yFaces(), -0.5);
    std::ranges::fill(velocity.zFaces(), 0.25);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 0.0);
    const auto projected = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    const auto aperture = std::ranges::find_if(
        projected.links,
        [](const auto& link) {
            return link.kind
                    == SceneFluidPressureFaceLinkKind::AuthoredOpening
                && link.faceIndex == invalidSceneFluidPressureFaceIndex;
        });
    check(projected.diagnostics.accepted
              && projected.diagnostics.authoredOpeningLinkCount == 1
              && openingFlux.samples.size() == 1
              && aperture != projected.links.end(),
          "tilted off-face intake participates in an accepted projection");
    if (aperture != projected.links.end()) {
        checkNear(
            aperture->predictedRelativeVolumeFlowRateCubicMetersPerSecond,
            openingFlux.samples.front()
                .relativeVolumeFlowRateCubicMetersPerSecond,
            0.0,
            "embedded aperture reuses exact off-face staggered flux");
    }
    check(projected.diagnostics
              .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
              < 2.0e-11,
          "embedded aperture projection closes control-volume continuity");

    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projected, velocity);
    const auto repeatedMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projected, velocity);
    check(momentum == repeatedMomentum
              && momentum.diagnostics.embeddedOpeningLinkCount == 1
              && momentum.diagnostics.normalEquationControlCount == 2,
          "tilted aperture deterministically reconstructs both incident controls through normal equations");
    check(openingFlux.samples.front()
              .surfaceSweepRateCubicMetersPerSecond == 0.0,
          "fixed tilted aperture has exact zero cap sweep for reconstruction audit");
    std::vector<bool> embeddedControls(momentum.controlVolumes.size(), false);
    for (const auto& link : fixture.faceLinks.links) {
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening) {
            embeddedControls[link.minusControlVolumeIndex] = true;
            embeddedControls[link.plusControlVolumeIndex] = true;
        }
    }
    bool stationaryNormalEquations = true;
    for (std::size_t controlIndex = 0;
         controlIndex < momentum.controlVolumes.size(); ++controlIndex) {
        if (!embeddedControls[controlIndex]) {
            continue;
        }
        fluid::Vector3 gradient;
        for (std::size_t linkIndex = 0;
             linkIndex < fixture.faceLinks.links.size(); ++linkIndex) {
            const auto& link = fixture.faceLinks.links[linkIndex];
            if (link.minusControlVolumeIndex != controlIndex
                && link.plusControlVolumeIndex != controlIndex) {
                continue;
            }
            const auto& normal = link.unitNormalMinusToPlus;
            const auto& reconstructed =
                momentum.controlVolumes[controlIndex]
                    .velocityMetersPerSecond;
            const double observed = projected.links[linkIndex]
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond
                / link.areaSquareMeters;
            const double residual =
                reconstructed.x * normal.x
                + reconstructed.y * normal.y
                + reconstructed.z * normal.z
                - observed;
            gradient.x += link.areaSquareMeters * normal.x * residual;
            gradient.y += link.areaSquareMeters * normal.y * residual;
            gradient.z += link.areaSquareMeters * normal.z * residual;
        }
        stationaryNormalEquations = stationaryNormalEquations
            && std::sqrt(
                   gradient.x * gradient.x
                   + gradient.y * gradient.y
                   + gradient.z * gradient.z) < 1.0e-11;
    }
    check(stationaryNormalEquations,
          "tilted aperture reconstruction satisfies its vector normal equations");

    SceneFluidRegionTransportSettings transportSettings;
    transportSettings.timeStepSeconds = strictSettings().timeStepSeconds;
    const auto transport = advanceSceneFluidRegionMomentum(
        momentum, fixture.faceLinks, projected, transportSettings);
    check(transport.diagnostics.accepted,
          "tilted aperture momentum advances through conservative transport");
    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = transportSettings.timeStepSeconds;
    stepSettings.substeps = 1;
    stepSettings.constraintIterations = 1;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    check(fixture.structure.step(stepSettings).finite,
          "fixed tilted aperture advances to a consecutive epoch");
    const auto currentState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structureAssembly.mappings,
        fixture.structure);
    const auto currentEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto currentFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, currentState,
        currentEpoch.openingCaps, currentEpoch.openingQuadrature,
        currentEpoch.openingPatches, grid(), velocity);
    const auto prediction = predictSceneFluidRegionLinkFlows(
        transport, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    const auto embeddedPrediction = std::ranges::find_if(
        prediction.links,
        [](const auto& link) {
            return link.faceIndex == invalidSceneFluidPressureFaceIndex;
        });
    bool exactEmbeddedPrediction = embeddedPrediction
        != prediction.links.end();
    if (exactEmbeddedPrediction) {
        const auto& source = currentEpoch.pressureFaceLinks.links[
            embeddedPrediction->linkIndex];
        const auto& minus = transport.controlVolumes[
            source.minusControlVolumeIndex].velocityMetersPerSecond;
        const auto& plus = transport.controlVolumes[
            source.plusControlVolumeIndex].velocityMetersPerSecond;
        const auto& normal = source.unitNormalMinusToPlus;
        const double minusNormal = minus.x * normal.x
            + minus.y * normal.y + minus.z * normal.z;
        const double plusNormal = plus.x * normal.x
            + plus.y * normal.y + plus.z * normal.z;
        const double expected = 0.5 * (minusNormal + plusNormal);
        exactEmbeddedPrediction =
            embeddedPrediction->predictedAbsoluteVelocityMetersPerSecond
            == expected;
    }
    check(prediction.diagnostics.embeddedOpeningLinkCount == 1
              && exactEmbeddedPrediction,
          "transport predictor projects endpoint momentum onto the tilted aperture normal");

    SceneFluidRegionWallSettings wallSettings;
    wallSettings.timeStepSeconds = transportSettings.timeStepSeconds;
    wallSettings.kinematicViscositySquareMetersPerSecond = 0.0;
    const auto wallExchange = exchangeSceneFluidRegionWallMomentum(
        transport, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, wallSettings);
    const auto wallPrediction = predictSceneFluidRegionLinkFlows(
        wallExchange, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    check(wallPrediction.diagnostics.embeddedOpeningLinkCount == 1
              && wallPrediction.links == prediction.links,
          "zero wall exchange preserves tilted aperture prediction exactly");
    SceneFluidPressureFaceLinkSettings rejectedMimeticLinkSettings;
    rejectedMimeticLinkSettings.minimumCenterDistanceMeters = 10.0;
    const auto rejectedMimeticFaceLinks =
        buildSceneFluidPressureFaceLinks(
            fixture.surface.definition, currentState, grid(),
            fixture.transfer, currentEpoch.gridEpoch,
            currentEpoch.openingCaps, currentEpoch.openingQuadrature,
            currentEpoch.openingPatches,
            currentEpoch.openingFaceCrossings,
            currentEpoch.cappedFacePartitions,
            currentEpoch.cellVolumes, fixture.connectivity,
            currentEpoch.pressureControlVolumes,
            rejectedMimeticLinkSettings);
    const auto rejectedMimeticControls =
        buildSceneFluidMimeticControlCells(
            fixture.surface.definition, currentState, grid(),
            currentEpoch.gridEpoch, currentEpoch.openingCaps,
            currentEpoch.openingQuadrature, currentEpoch.openingPatches,
            currentEpoch.pressureControlVolumes,
            rejectedMimeticFaceLinks);
    const auto rejectedMimeticSystem =
        buildSceneFluidMimeticTraceSystem(rejectedMimeticControls);
    const auto rejectedWallMimeticFlows =
        sampleSceneFluidMimeticTraceFlows(
            rejectedMimeticControls, rejectedMimeticSystem,
            rejectedMimeticFaceLinks, currentFlux, wallExchange);
    const auto rejectedOpeningTrace = std::ranges::find(
        rejectedWallMimeticFlows.traces,
        SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace,
        &SceneFluidMimeticPredictedTraceFlow::kind);
    check(rejectedMimeticFaceLinks.unresolvedEmbeddedOpeningPatchCount == 1
              && rejectedOpeningTrace
                  != rejectedWallMimeticFlows.traces.end()
              && std::isfinite(rejectedOpeningTrace
                    ->predictedRelativeVolumeFlowRateCubicMetersPerSecond)
              && rejectedWallMimeticFlows.regionWallExchangeFingerprint
                  == wallExchange.fingerprint,
          "material-wall-adjusted prediction directly covers an authored opening after its two-point graph link is rejected");
}

void testConcaveEmbeddedOpeningProjection() {
    Fixture fixture(concaveOpenScene(), true);
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 1.25);
    std::ranges::fill(velocity.yFaces(), -0.1);
    std::ranges::fill(velocity.zFaces(), 0.2);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 0.0);
    const auto projection = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
    check(fixture.caps.caps.size() == 1
              && fixture.caps.caps.front().triangleCount == 2
              && fixture.openingQuadrature.points.size() == 2
              && fixture.openingPatches.patches.size() == 2
              && fixture.faceLinks.embeddedOpeningLinkCount == 2,
          "concave intake triangulation reaches two exact off-face pressure links");
    check(projection.diagnostics.accepted
              && projection.diagnostics.authoredOpeningLinkCount == 2
              && projection.diagnostics
                     .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  < 2.0e-11
              && momentum.diagnostics.embeddedOpeningLinkCount == 2
              && momentum.diagnostics.normalEquationControlCount == 2,
          "concave off-face intake closes pressure and explicit-normal region momentum");
    validateSceneFluidRegionMomentumState(
        momentum, grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
}

void testAuthoredNonPlanarEmbeddedOpeningProjection() {
    Fixture fixture(authoredNonPlanarOpenScene(), true);
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 1.25);
    std::ranges::fill(velocity.yFaces(), -0.1);
    std::ranges::fill(velocity.zFaces(), 0.2);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 0.0);
    const auto projection = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
    check(fixture.caps.triangles.size() == 2
              && fixture.volumes.openingCapCount == 1
              && fixture.volumes.openingCapAreaSquareMeters
                  == fixture.caps.totalAreaSquareMeters
              && fixture.openingPatches.patches.size() == 2
              && fixture.faceLinks.embeddedOpeningLinkCount == 2,
          "authored nonplanar intake reaches exact cut-cell volume and embedded links");
    check(fixture.caps.triangles[0].unitNormalNegativeToPositive
                  .x
              != fixture.caps.triangles[1]
                     .unitNormalNegativeToPositive
                     .x
              || fixture.caps.triangles[0]
                     .unitNormalNegativeToPositive
                     .y
                  != fixture.caps.triangles[1]
                         .unitNormalNegativeToPositive
                         .y
              || fixture.caps.triangles[0]
                     .unitNormalNegativeToPositive
                     .z
                  != fixture.caps.triangles[1]
                         .unitNormalNegativeToPositive
                         .z,
          "authored nonplanar intake keeps distinct per-triangle normals");
    check(projection.diagnostics.accepted
              && projection.diagnostics.authoredOpeningLinkCount == 2
              && projection.diagnostics
                     .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  < 2.0e-11
              && momentum.diagnostics.embeddedOpeningLinkCount == 2
              && momentum.diagnostics.normalEquationControlCount == 2,
          "authored nonplanar intake closes pressure and explicit-normal momentum");
    validateSceneFluidRegionMomentumState(
        momentum, grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
}

void testZeroFlowAndRejectedAttempt() {
    Fixture fixture(nestedScene());
    fluid::MacVelocityField zeroVelocity(grid());
    const auto zeroFlux = fixture.flux(zeroVelocity);
    std::vector<double> zeroWarm(fixture.pressureOperator.rows.size(), 0.0);
    const auto zero = fixture.project(
        zeroVelocity, zeroFlux, zeroWarm, strictSettings());
    check(zero.diagnostics.accepted
              && zero.diagnostics.pressureSolve.iterationCount == 0
              && zero.diagnostics
                  .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  == 0.0
              && zero.diagnostics
                  .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  == 0.0,
          "zero predicted flow is an exact zero-iteration projection");

    const auto velocity = manufacturedVelocity();
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(
        fixture.pressureOperator.rows.size(), 4.25);
    auto settings = strictSettings();
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-30;
    settings.pressureSolve.relativeResidualTolerance = 0.0;
    settings.pressureSolve.maximumIterations = 0;
    const auto rejected = fixture.project(
        velocity, openingFlux, warm, settings);
    check(!rejected.diagnostics.accepted
              && rejected.diagnostics.finite
              && rejected.diagnostics.pressureSolve.compatible
              && !rejected.diagnostics.pressureSolve.converged
              && rejected.pressurePascals.empty()
              && std::ranges::all_of(
                  rejected.links,
                  [](const auto& link) {
                      return link
                                  .pressureCorrectionVolumeFlowRateCubicMetersPerSecond
                                  == 0.0
                          && link
                                  .correctedRelativeVolumeFlowRateCubicMetersPerSecond
                                  == 0.0;
                  }),
          "non-converged attempt exposes diagnostics but no corrected state");
    validateSceneFluidPressureProjectionIntegrity(rejected);
}

void testLinkResolvedContinuation() {
    Fixture fixture(openScene(), true);
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 2.0);
    const auto previousFlux = fixture.flux(velocity);
    std::vector<double> previousWarm(
        fixture.pressureOperator.rows.size(), 0.0);
    const auto previousProjection = fixture.project(
        velocity, previousFlux, previousWarm, strictSettings());

    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = strictSettings().timeStepSeconds;
    stepSettings.substeps = 1;
    stepSettings.constraintIterations = 1;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    check(fixture.structure.step(stepSettings).finite,
          "link continuation current Structure state advances");
    const auto currentState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structureAssembly.mappings,
        fixture.structure);
    const auto currentEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto currentFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, currentState,
        currentEpoch.openingCaps, currentEpoch.openingQuadrature,
        currentEpoch.openingPatches, grid(), velocity);
    const auto continuation = continueSceneFluidPressureLinkFlows(
        grid(), fixture.faceLinks, fixture.openingPatches,
        previousProjection, currentEpoch.pressureFaceLinks, currentFlux,
        velocity);
    const auto repeated = continueSceneFluidPressureLinkFlows(
        grid(), fixture.faceLinks, fixture.openingPatches,
        previousProjection, currentEpoch.pressureFaceLinks, currentFlux,
        velocity);
    check(continuation == repeated
              && continuation.diagnostics.finite
              && continuation.diagnostics.multiLinkFaceCount > 0
              && continuation.diagnostics.openingLinkCount == 1
              && continuation.diagnostics
                     .maximumCarriedAbsoluteVelocityDeviationMetersPerSecond
                  > 0.0
              && continuation.diagnostics
                     .maximumAbsoluteFaceFlowClosureCubicMetersPerSecond
                  < 1.0e-15,
          "link continuation deterministically preserves subface velocity with exact face-total closure");

    std::vector<double> currentWarm(
        currentEpoch.pressureOperator.rows.size(), 0.0);
    const auto continuedProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        currentFlux, velocity, continuation, currentEpoch.cellVolumes,
        fixture.connectivity, currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.pressureOperator,
        currentWarm, strictSettings());
    const auto directProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        currentFlux, velocity, currentEpoch.cellVolumes,
        fixture.connectivity, currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.pressureOperator,
        currentWarm, strictSettings());
    bool exactPrediction = continuedProjection.diagnostics.accepted
        && continuedProjection.linkFlowContinuationFingerprint
            == continuation.fingerprint
        && continuedProjection.links.size() == continuation.links.size();
    for (std::size_t index = 0;
         exactPrediction && index < continuation.links.size(); ++index) {
        exactPrediction = continuedProjection.links[index]
            .predictedRelativeVolumeFlowRateCubicMetersPerSecond
            == continuation.links[index]
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond;
    }
    check(exactPrediction,
          "pressure projection consumes the exact topology-bound continued link flows");
    const auto maximumPressure = [](const SceneFluidPressureProjection& value) {
        double result = 0.0;
        for (const double pressure : value.pressurePascals) {
            result = std::max(result, std::abs(pressure));
        }
        return result;
    };
    check(maximumPressure(directProjection) > 1.0e-6
              && maximumPressure(continuedProjection)
                  < 0.2 * maximumPressure(directProjection),
          "static subface continuation cannot replace region-resolved momentum transport without strongly reducing steady pressure load");

    auto corrupt = continuation;
    corrupt.links.front()
        .predictedRelativeVolumeFlowRateCubicMetersPerSecond += 0.01;
    expectInvalid(
        [&] {
            validateSceneFluidPressureLinkFlowContinuationIntegrity(corrupt);
        },
        "link continuation integrity rejects predicted-flow corruption");

    auto foreignVelocity = velocity;
    foreignVelocity.xFaces().front() += 0.01;
    expectInvalid(
        [&] {
            static_cast<void>(continueSceneFluidPressureLinkFlows(
                grid(), fixture.faceLinks, fixture.openingPatches,
                previousProjection, currentEpoch.pressureFaceLinks,
                currentFlux, foreignVelocity));
        },
        "link continuation rejects a bulk field foreign to opening flux");

    SceneFluidPressureLinkFlowContinuationLimits limits;
    limits.maximumLinks = currentEpoch.pressureFaceLinks.links.size() - 1;
    expectLimited(
        [&] {
            static_cast<void>(continueSceneFluidPressureLinkFlows(
                grid(), fixture.faceLinks, fixture.openingPatches,
                previousProjection, currentEpoch.pressureFaceLinks,
                currentFlux, velocity, limits));
        },
        "link continuation bounds its output link count");
}

void testAreaChangingLinkContinuation() {
    Fixture fixture(openScene(), false, true);
    const auto velocity = manufacturedVelocity();
    const auto previousFlux = fixture.flux(velocity);
    std::vector<double> warm(
        fixture.pressureOperator.rows.size(), 0.0);
    const auto previousProjection = fixture.project(
        velocity, previousFlux, warm, strictSettings());
    const auto previousMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, previousProjection, velocity);
    SceneFluidRegionTransportSettings transportSettings;
    transportSettings.timeStepSeconds = strictSettings().timeStepSeconds;
    const auto transport = advanceSceneFluidRegionMomentum(
        previousMomentum, fixture.faceLinks, previousProjection,
        transportSettings);
    check(transport.diagnostics.accepted,
          "area-changing link continuation advances region momentum");

    const auto apex = std::ranges::find(
        fixture.structureAssembly.mappings.nodeVertexIds, StableId{10});
    check(apex != fixture.structureAssembly.mappings.nodeVertexIds.end(),
          "area-changing link continuation finds its free apex");
    if (apex == fixture.structureAssembly.mappings.nodeVertexIds.end()) {
        return;
    }
    fixture.structure.addExternalForce(
        static_cast<std::size_t>(
            apex - fixture.structureAssembly.mappings.nodeVertexIds.begin()),
        {-0.2, 0.0, 0.0});
    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = strictSettings().timeStepSeconds;
    stepSettings.substeps = 4;
    stepSettings.constraintIterations = 10;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    check(fixture.structure.step(stepSettings).finite,
          "area-changing link continuation advances its free apex");
    const auto currentState = captureSceneFluidSurfaceState(
        fixture.surface.definition, fixture.structureAssembly.mappings,
        fixture.structure);
    const auto currentEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto currentFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, currentState,
        currentEpoch.openingCaps, currentEpoch.openingQuadrature,
        currentEpoch.openingPatches, grid(), velocity);
    const auto continuation = continueSceneFluidPressureLinkFlows(
        grid(), fixture.faceLinks, fixture.openingPatches,
        previousProjection, currentEpoch.pressureFaceLinks, currentFlux,
        velocity);
    check(continuation.diagnostics.finite
              && continuation.diagnostics
                     .maximumAreaRenormalizationMetersPerSecond > 0.0
              && continuation.diagnostics
                     .maximumAbsoluteFaceFlowClosureCubicMetersPerSecond
                  < 1.0e-15,
          "area-changing link continuation recenters carried deviations without changing face-total flow");

    SceneFluidRegionWallSettings wallSettings;
    wallSettings.timeStepSeconds = transportSettings.timeStepSeconds;
    const auto wallExchange = exchangeSceneFluidRegionWallMomentum(
        transport, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, wallSettings);
    const auto repeatedWallExchange = exchangeSceneFluidRegionWallMomentum(
        transport, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, wallSettings);
    check(wallExchange == repeatedWallExchange
              && wallExchange.diagnostics.accepted
              && wallExchange.diagnostics.finite
              && wallExchange.diagnostics.quadraturePointCount
                  == currentEpoch.gridEpoch.quadrature.points.size()
              && wallExchange.diagnostics
                     .maximumRelativeTangentialSpeedMetersPerSecond > 0.0
              && wallExchange.diagnostics.viscousDissipationJoules >= 0.0
              && wallExchange.diagnostics
                     .momentumResidualNormKilogramMetersPerSecond < 1.0e-12,
          "two-sided material-wall exchange is deterministic, dissipative, and momentum conservative");
    validateSceneFluidRegionWallExchange(
        wallExchange, transport, grid(),
        currentEpoch.pressureControlVolumes, fixture.surface.definition,
        currentState, currentEpoch.gridEpoch.quadrature);
    bool tangentialTraction = true;
    for (const auto& sample : wallExchange.samples) {
        tangentialTraction = tangentialTraction
            && std::abs(
                sample.structureTraction.tractionPascals.x
                    * sample.unitNormalNegativeToPositive.x
                + sample.structureTraction.tractionPascals.y
                    * sample.unitNormalNegativeToPositive.y
                + sample.structureTraction.tractionPascals.z
                    * sample.unitNormalNegativeToPositive.z) < 1.0e-12;
    }
    check(tangentialTraction,
          "material-wall exchange leaves normal traction exclusively to pressure");
    const auto acceptedWall = captureSceneFluidAcceptedWallTractions(
        wallExchange);
    validateSceneFluidAcceptedWallTractions(
        acceptedWall, currentEpoch.gridEpoch.quadrature,
        wallExchange.fingerprint);
    auto corruptAcceptedWall = acceptedWall;
    corruptAcceptedWall.tractions.front().tractionPascals.x += 1.0;
    expectInvalid(
        [&] {
            validateSceneFluidAcceptedWallTractionSetIntegrity(
                corruptAcceptedWall);
        },
        "accepted material-wall endpoint rejects traction corruption");

    auto zeroWallSettings = wallSettings;
    zeroWallSettings.kinematicViscositySquareMetersPerSecond = 0.0;
    const auto zeroWallExchange = exchangeSceneFluidRegionWallMomentum(
        transport, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, zeroWallSettings);
    check(zeroWallExchange.diagnostics.accepted
              && zeroWallExchange.diagnostics.fluidImpulseKilogramMetersPerSecond
                  == fluid::Vector3{}
              && zeroWallExchange.diagnostics.structureImpulseKilogramMetersPerSecond
                  == fluid::Vector3{}
              && zeroWallExchange.diagnostics.viscousDissipationJoules == 0.0,
          "zero material-wall viscosity is an exact no-exchange fixed point");

    auto limitedWallSettings = wallSettings;
    limitedWallSettings.kinematicViscositySquareMetersPerSecond = 1.0;
    limitedWallSettings.maximumSubsteps = 1;
    const auto limitedWallExchange = exchangeSceneFluidRegionWallMomentum(
        transport, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, limitedWallSettings);
    check(!limitedWallExchange.diagnostics.accepted
              && limitedWallExchange.diagnostics.failureStage
                  == SceneFluidRegionWallFailureStage::SubstepLimit
              && limitedWallExchange.controlVolumes.empty()
              && limitedWallExchange.samples.empty(),
          "material-wall exchange rejects an unsafe explicit step without publishing partial state");

    auto corruptWallExchange = wallExchange;
    corruptWallExchange.samples.front().structureTraction.tractionPascals.x +=
        1.0;
    expectInvalid(
        [&] {
            validateSceneFluidRegionWallExchangeIntegrity(
                corruptWallExchange);
        },
        "material-wall exchange integrity rejects traction corruption");

    const auto regionPrediction = predictSceneFluidRegionLinkFlows(
        transport, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    const auto wallPrediction = predictSceneFluidRegionLinkFlows(
        wallExchange, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    const auto zeroWallPrediction = predictSceneFluidRegionLinkFlows(
        zeroWallExchange, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    const auto repeatedPrediction = predictSceneFluidRegionLinkFlows(
        transport, grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    check(regionPrediction == repeatedPrediction
              && regionPrediction.diagnostics.finite
              && regionPrediction.diagnostics.openingLinkCount == 1
              && regionPrediction.diagnostics
                     .maximumAbsoluteVolumeChangeCubicMeters > 0.0
              && regionPrediction.links.size()
                  == currentEpoch.pressureFaceLinks.links.size(),
          "transported region momentum deterministically predicts the moving epoch link flow and GCL remap");
    check(wallPrediction.sourceTransportFingerprint == 0
              && wallPrediction.sourceWallExchangeFingerprint
                  == wallExchange.fingerprint
              && wallPrediction.currentPressureControlVolumeFingerprint
                  == currentEpoch.pressureControlVolumes.fingerprint
              && wallPrediction.links.size()
                  == currentEpoch.pressureFaceLinks.links.size(),
          "material-wall-adjusted region momentum predicts the same current pressure-link topology");
    bool zeroWallPreservesPrediction =
        zeroWallPrediction.links.size() == regionPrediction.links.size();
    for (std::size_t index = 0;
         zeroWallPreservesPrediction
             && index < zeroWallPrediction.links.size(); ++index) {
        zeroWallPreservesPrediction =
            zeroWallPrediction.links[index]
                .predictedAbsoluteVelocityMetersPerSecond
                == regionPrediction.links[index]
                    .predictedAbsoluteVelocityMetersPerSecond
            && zeroWallPrediction.links[index]
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond
                == regionPrediction.links[index]
                    .predictedRelativeVolumeFlowRateCubicMetersPerSecond;
    }
    check(zeroWallPreservesPrediction,
          "zero wall viscosity preserves the transported link-flow predictor exactly");
    validateSceneFluidRegionLinkFlowPrediction(
        regionPrediction, transport, grid(),
        currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);
    validateSceneFluidRegionLinkFlowPrediction(
        wallPrediction, wallExchange, grid(),
        currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentFlux);

    const auto volumeRates = buildSceneFluidPressureVolumeRates(
        fixture.volumes, currentEpoch.cellVolumes,
        currentEpoch.pressureControlVolumes);
    const auto mimeticControls = buildSceneFluidMimeticControlCells(
        fixture.surface.definition, currentState, grid(),
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks);
    const auto mimeticSystem = buildSceneFluidMimeticTraceSystem(
        mimeticControls);
    const auto mimeticTraceFlows = sampleSceneFluidMimeticTraceFlows(
        mimeticControls, mimeticSystem,
        currentEpoch.pressureFaceLinks, currentFlux, grid(), velocity);
    const auto wallMimeticTraceFlows = sampleSceneFluidMimeticTraceFlows(
        mimeticControls, mimeticSystem,
        currentEpoch.pressureFaceLinks, currentFlux, wallExchange);
    SceneFluidMimeticPressureSourceSettings mimeticSourceSettings;
    mimeticSourceSettings.densityKgPerCubicMeter =
        strictSettings().densityKgPerCubicMeter;
    mimeticSourceSettings.timeStepSeconds =
        strictSettings().timeStepSeconds;
    const auto mimeticSources = buildSceneFluidMimeticPressureSources(
        mimeticControls, mimeticSystem, mimeticTraceFlows, volumeRates,
        mimeticSourceSettings);
    const auto wallMimeticSources = buildSceneFluidMimeticPressureSources(
        mimeticControls, mimeticSystem, wallMimeticTraceFlows, volumeRates,
        mimeticSourceSettings);
    bool exactMimeticVolumeRates = mimeticSources.controls.size()
        == mimeticControls.controlCells.size();
    for (std::size_t index = 0;
         exactMimeticVolumeRates && index < mimeticSources.controls.size();
         ++index) {
        const auto& control = mimeticControls.controlCells[index];
        const auto& source = mimeticSources.controls[index];
        const auto& rate = volumeRates.controlVolumes[
            control.controlVolumeIndex];
        exactMimeticVolumeRates =
            source.geometryVolumeChangeRateCubicMetersPerSecond
                == rate.volumeChangeRateCubicMetersPerSecond
            && source.predictedContinuityResidualCubicMetersPerSecond
                == source.geometryVolumeChangeRateCubicMetersPerSecond
                    + source
                        .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond
            && source.integratedSourcePascalsMeters
                == -mimeticSourceSettings.densityKgPerCubicMeter
                    / mimeticSourceSettings.timeStepSeconds
                    * source
                        .predictedContinuityResidualCubicMetersPerSecond;
    }
    check(mimeticSources.mimeticTraceFlowFingerprint
              == mimeticTraceFlows.fingerprint
              && mimeticSources.pressureVolumeRateFingerprint
                  == volumeRates.fingerprint
              && exactMimeticVolumeRates
              && mimeticSources
                    .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond
                  < 2.0e-12,
          "mimetic pressure sources bind exact shared-trace flow and accepted consecutive-epoch GCL rates");
    bool exactWallGraphOverlap = true;
    for (const auto& trace : wallMimeticTraceFlows.traces) {
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            const auto graph = std::ranges::find(
                wallPrediction.links, trace.sourceStableId,
                &SceneFluidRegionPredictedLinkFlow::stableId);
            exactWallGraphOverlap = graph != wallPrediction.links.end()
                && graph
                    ->predictedRelativeVolumeFlowRateCubicMetersPerSecond
                    == trace
                        .predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        } else {
            const auto graph = std::ranges::find(
                wallPrediction.links, trace.sourceStableId,
                &SceneFluidRegionPredictedLinkFlow::openingPatchStableId);
            if (graph != wallPrediction.links.end()) {
                exactWallGraphOverlap = graph
                    ->predictedRelativeVolumeFlowRateCubicMetersPerSecond
                    == trace
                        .predictedRelativeVolumeFlowRateCubicMetersPerSecond;
            }
        }
        if (!exactWallGraphOverlap) break;
    }
    check(wallMimeticTraceFlows.regionWallExchangeFingerprint
              == wallExchange.fingerprint
              && wallMimeticTraceFlows.sourceDensityKgPerCubicMeter
                  == wallExchange.densityKgPerCubicMeter
              && wallMimeticTraceFlows.traces.size()
                  == mimeticSystem.sharedTraceCount
              && exactWallGraphOverlap
              && wallMimeticSources.mimeticTraceFlowFingerprint
                  == wallMimeticTraceFlows.fingerprint
              && wallMimeticSources.pressureVolumeRateFingerprint
                  == volumeRates.fingerprint,
          "material-wall-adjusted region velocities reproduce every existing graph predictor and extend over the complete mimetic trace topology");
    auto wrongMimeticDensitySettings = mimeticSourceSettings;
    wrongMimeticDensitySettings.densityKgPerCubicMeter *= 0.5;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidMimeticPressureSources(
            mimeticControls, mimeticSystem, wallMimeticTraceFlows,
            volumeRates, wrongMimeticDensitySettings)); },
        "mimetic pressure sources reject a wall-predictor density mismatch");
    auto wrongMimeticSourceSettings = mimeticSourceSettings;
    wrongMimeticSourceSettings.timeStepSeconds *= 0.5;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidMimeticPressureSources(
            mimeticControls, mimeticSystem, mimeticTraceFlows, volumeRates,
            wrongMimeticSourceSettings)); },
        "mimetic pressure sources reject a volume-rate duration mismatch");
    std::vector<double> currentWarm(
        currentEpoch.pressureOperator.rows.size(), 0.0);
    const auto regionProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        currentFlux, velocity, regionPrediction, currentEpoch.cellVolumes,
        fixture.connectivity, currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.pressureOperator,
        volumeRates, currentWarm, strictSettings());
    const auto wallProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        currentFlux, velocity, wallPrediction, currentEpoch.cellVolumes,
        fixture.connectivity, currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.pressureOperator,
        volumeRates, currentWarm, strictSettings());
    bool exactPredictedFlow = regionProjection.links.size()
        == regionPrediction.links.size();
    for (std::size_t index = 0;
         exactPredictedFlow && index < regionProjection.links.size();
         ++index) {
        exactPredictedFlow = regionProjection.links[index]
            .predictedRelativeVolumeFlowRateCubicMetersPerSecond
            == regionPrediction.links[index]
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond;
    }
    check(regionProjection.diagnostics.accepted
              && regionProjection.regionLinkFlowPredictionFingerprint
                  == regionPrediction.fingerprint
              && regionProjection.linkFlowContinuationFingerprint == 0
              && regionProjection.pressureVolumeRateFingerprint
                  == volumeRates.fingerprint
              && exactPredictedFlow,
          "moving pressure projection consumes the transported region link predictor exactly");
    check(wallProjection.diagnostics.accepted
              && wallProjection.regionLinkFlowPredictionFingerprint
                  == wallPrediction.fingerprint
              && wallProjection.regionWallExchangeFingerprint
                  == wallExchange.fingerprint,
          "moving pressure projection retains exact material-wall predictor provenance");

    const auto currentMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.openingPatches,
        regionProjection, velocity);
    const auto movingTransport = advanceSceneFluidRegionMomentum(
        currentMomentum, currentEpoch.pressureFaceLinks, regionProjection,
        transportSettings);
    bool exactGclVolumes = movingTransport.controlVolumes.size()
        == currentMomentum.controlVolumes.size();
    for (std::size_t index = 0;
         exactGclVolumes && index < movingTransport.controlVolumes.size();
         ++index) {
        exactGclVolumes = movingTransport.controlVolumes[index]
            .volumeCubicMeters
            == currentMomentum.controlVolumes[index].volumeCubicMeters
                + transportSettings.timeStepSeconds
                    * regionProjection.controlVolumes[index]
                        .geometryVolumeChangeRateCubicMetersPerSecond;
    }
    check(movingTransport.diagnostics.accepted
              && movingTransport.diagnostics.usesMovingVolumeRates
              && movingTransport.pressureVolumeRateFingerprint
                  == volumeRates.fingerprint
              && movingTransport.diagnostics
                     .maximumAbsoluteGeometryVolumeChangeCubicMeters > 0.0
              && exactGclVolumes,
          "region transport advances moving control volumes through the accepted discrete GCL rate");
    validateSceneFluidRegionTransport(
        movingTransport, currentMomentum,
        currentEpoch.pressureFaceLinks, regionProjection);

    auto corruptPrediction = regionPrediction;
    corruptPrediction.links.front()
        .predictedRelativeVolumeFlowRateCubicMetersPerSecond += 0.01;
    expectInvalid(
        [&] {
            validateSceneFluidRegionLinkFlowPredictionIntegrity(
                corruptPrediction);
        },
        "region link-flow integrity rejects predicted-flow corruption");

    SceneFluidRegionLinkFlowLimits predictionLimits;
    predictionLimits.maximumLinks =
        currentEpoch.pressureFaceLinks.links.size() - 1;
    expectLimited(
        [&] {
            static_cast<void>(predictSceneFluidRegionLinkFlows(
                transport, grid(), currentEpoch.pressureControlVolumes,
                currentEpoch.pressureFaceLinks, currentFlux,
                predictionLimits));
        },
        "region link-flow prediction bounds its moving link state");
}

void testRegionMomentumReconstruction() {
    Fixture fixture(openScene(), true);
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 2.0);
    std::ranges::fill(velocity.yFaces(), -0.25);
    std::ranges::fill(velocity.zFaces(), 0.5);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(
        fixture.pressureOperator.rows.size(), 0.0);
    const auto projection = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
    const auto repeated = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
    check(momentum == repeated
              && momentum.fingerprint != 0
              && momentum.diagnostics.finite
              && momentum.diagnostics.controlVolumeCount
                  == fixture.pressureVolumes.controlVolumes.size()
              && momentum.diagnostics.linkCount
                  == fixture.faceLinks.links.size()
              && momentum.diagnostics.openingLinkCount == 1
              && momentum.diagnostics.sampledComponentCount
                     + momentum.diagnostics.fallbackComponentCount
                  == 3 * momentum.controlVolumes.size()
              && momentum.diagnostics.kineticEnergyJoules > 0.0
              && momentum.diagnostics
                     .maximumAbsoluteVelocityMetersPerSecond > 0.0,
          "region momentum deterministically reconstructs finite cell-region vectors from corrected links");
    validateSceneFluidRegionMomentumState(
        momentum, grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);

    auto corrupt = momentum;
    corrupt.controlVolumes.front().velocityMetersPerSecond.x += 0.01;
    expectInvalid(
        [&] { validateSceneFluidRegionMomentumStateIntegrity(corrupt); },
        "region momentum integrity rejects velocity corruption");

    auto foreignVelocity = velocity;
    foreignVelocity.xFaces().front() += 0.01;
    expectInvalid(
        [&] {
            static_cast<void>(reconstructSceneFluidRegionMomentumState(
                grid(), fixture.pressureVolumes, fixture.faceLinks,
                fixture.openingPatches, projection, foreignVelocity));
        },
        "region momentum rejects fallback velocity foreign to the pressure projection");

    SceneFluidRegionMomentumLimits limits;
    limits.maximumControlVolumes =
        fixture.pressureVolumes.controlVolumes.size() - 1;
    expectLimited(
        [&] {
            static_cast<void>(reconstructSceneFluidRegionMomentumState(
                grid(), fixture.pressureVolumes, fixture.faceLinks,
                fixture.openingPatches, projection, velocity, limits));
        },
        "region momentum bounds its control-volume state");

    Fixture zeroFixture(nestedScene(), true);
    fluid::MacVelocityField zeroVelocity(grid());
    const auto zeroFlux = zeroFixture.flux(zeroVelocity);
    std::vector<double> zeroWarm(
        zeroFixture.pressureOperator.rows.size(), 0.0);
    const auto zeroProjection = zeroFixture.project(
        zeroVelocity, zeroFlux, zeroWarm, strictSettings());
    const auto zeroMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), zeroFixture.pressureVolumes, zeroFixture.faceLinks,
        zeroFixture.openingPatches, zeroProjection, zeroVelocity);
    check(zeroMomentum.diagnostics.finite
              && zeroMomentum.diagnostics
                     .totalMomentumKilogramMetersPerSecond
                  == fluid::Vector3{}
              && zeroMomentum.diagnostics.kineticEnergyJoules == 0.0
              && zeroMomentum.diagnostics
                     .maximumAbsoluteVelocityMetersPerSecond == 0.0
              && zeroMomentum.diagnostics
                     .maximumLinkNormalVelocityResidualMetersPerSecond
                  == 0.0,
          "zero corrected flow reconstructs exact zero region momentum");
}

void testRegionMomentumTransport() {
    Fixture fixture(nestedScene(), true);
    const auto velocity = manufacturedVelocity();
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(
        fixture.pressureOperator.rows.size(), 0.0);
    const auto projection = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    const auto momentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, projection, velocity);
    SceneFluidRegionTransportSettings settings;
    settings.timeStepSeconds = 0.01;
    const auto transport = advanceSceneFluidRegionMomentum(
        momentum, fixture.faceLinks, projection, settings);
    const auto repeated = advanceSceneFluidRegionMomentum(
        momentum, fixture.faceLinks, projection, settings);
    const double energyTolerance = settings.absoluteEnergyToleranceJoules
        + settings.relativeEnergyTolerance
            * std::max(1.0, transport.diagnostics.kineticEnergyBeforeJoules);
    const double momentumTolerance =
        settings.absoluteMomentumToleranceKilogramMetersPerSecond
        + settings.relativeMomentumTolerance
            * std::max(
                1.0,
                std::sqrt(
                    transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.x
                        * transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.x
                    + transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.y
                        * transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.y
                    + transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.z
                        * transport.diagnostics
                            .momentumBeforeKilogramMetersPerSecond.z));
    check(transport == repeated
              && transport.fingerprint != 0
              && transport.diagnostics.finite
              && transport.diagnostics.accepted
              && transport.diagnostics.failureStage
                  == SceneFluidRegionTransportFailureStage::None
              && transport.diagnostics.substepCount > 0
              && transport.diagnostics
                     .maximumAcceptedSubstepOutgoingCourantNumber
                  <= settings.maximumOutgoingCourantNumber
              && transport.diagnostics
                     .maximumAcceptedSubstepViscousNumber
                  <= settings.maximumViscousNumber
              && transport.controlVolumes.size()
                  == momentum.controlVolumes.size()
              && transport.diagnostics.kineticEnergyAfterJoules
                  <= transport.diagnostics.kineticEnergyBeforeJoules
                      + energyTolerance
              && transport.diagnostics.advectiveEnergyLossJoules
                  >= -energyTolerance
              && transport.diagnostics.viscousEnergyLossJoules
                  >= -energyTolerance
              && transport.diagnostics
                     .momentumResidualNormKilogramMetersPerSecond
                  <= momentumTolerance,
          "region transport deterministically advances conservative donor-cell momentum with dissipative graph viscosity");
    checkNear(
        transport.diagnostics.advectiveEnergyLossJoules
            + transport.diagnostics.viscousEnergyLossJoules,
        transport.diagnostics.kineticEnergyBeforeJoules
            - transport.diagnostics.kineticEnergyAfterJoules,
        energyTolerance,
        "region transport closes its split energy-loss ledger");
    validateSceneFluidRegionTransport(
        transport, momentum, fixture.faceLinks, projection);

    auto acceleratedVelocity = velocity;
    for (double& value : acceleratedVelocity.xFaces()) {
        value += 0.125;
    }
    const auto accelerated = advanceSceneFluidRegionMomentum(
        momentum, fixture.faceLinks, projection, grid(), velocity,
        acceleratedVelocity, settings);
    const double expectedImpulse = projection.settings.densityKgPerCubicMeter
        * (grid().upperMeters().x - grid().lowerMeters().x)
        * (grid().upperMeters().y - grid().lowerMeters().y)
        * (grid().upperMeters().z - grid().lowerMeters().z) * 0.125;
    check(accelerated.diagnostics.accepted,
          "bulk-increment region transport remains accepted");
    check(accelerated.diagnostics.usesBulkVelocityIncrement
              && accelerated.previousBulkVelocityFingerprint
                  == sceneFluidOpeningFluxVelocityFingerprint(
                      grid(), velocity)
              && accelerated.currentBulkVelocityFingerprint
                  == sceneFluidOpeningFluxVelocityFingerprint(
                      grid(), acceleratedVelocity),
          "region transport explicitly binds the split bulk-MAC velocity increment");
    checkNear(
        accelerated.diagnostics
            .maximumBulkVelocityIncrementMetersPerSecond,
        0.125, 1.0e-15,
          "region transport explicitly binds the split bulk-MAC velocity increment");
    checkNear(
        accelerated.diagnostics
            .bulkVelocityIncrementImpulseKilogramMetersPerSecond.x,
        expectedImpulse, 1.0e-12,
        "uniform bulk-MAC increment delivers its exact domain impulse to all regions");

    auto corrupt = transport;
    corrupt.controlVolumes.front().velocityMetersPerSecond.x += 0.01;
    expectInvalid(
        [&] { validateSceneFluidRegionTransportIntegrity(corrupt); },
        "region transport integrity rejects velocity corruption");

    auto limitedSettings = settings;
    limitedSettings.timeStepSeconds = 100.0;
    limitedSettings.maximumSubsteps = 1;
    const auto limited = advanceSceneFluidRegionMomentum(
        momentum, fixture.faceLinks, projection, limitedSettings);
    check(!limited.diagnostics.accepted
              && limited.diagnostics.finite
              && limited.diagnostics.failureStage
                  == SceneFluidRegionTransportFailureStage::SubstepLimit
              && limited.controlVolumes.empty(),
          "region transport rejects an excessive substep demand without publishing momentum");

    SceneFluidRegionTransportLimits limits;
    limits.maximumLinks = fixture.faceLinks.links.size() - 1;
    expectLimited(
        [&] {
            static_cast<void>(advanceSceneFluidRegionMomentum(
                momentum, fixture.faceLinks, projection, settings, limits));
        },
        "region transport bounds its pressure-link input");
    limits.maximumLinks = fixture.faceLinks.links.size();
    limits.maximumTransportBytes = 0;
    expectLimited(
        [&] {
            static_cast<void>(advanceSceneFluidRegionMomentum(
                momentum, fixture.faceLinks, projection, settings, limits));
        },
        "region transport bounds its complete working storage");

    fluid::MacVelocityField zeroVelocity(grid());
    const auto zeroFlux = fixture.flux(zeroVelocity);
    std::vector<double> zeroWarm(
        fixture.pressureOperator.rows.size(), 0.0);
    const auto zeroProjection = fixture.project(
        zeroVelocity, zeroFlux, zeroWarm, strictSettings());
    const auto zeroMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), fixture.pressureVolumes, fixture.faceLinks,
        fixture.openingPatches, zeroProjection, zeroVelocity);
    const auto zeroTransport = advanceSceneFluidRegionMomentum(
        zeroMomentum, fixture.faceLinks, zeroProjection, settings);
    check(zeroTransport.diagnostics.accepted
              && zeroTransport.diagnostics
                     .momentumBeforeKilogramMetersPerSecond
                  == fluid::Vector3{}
              && zeroTransport.diagnostics
                     .momentumAfterKilogramMetersPerSecond
                  == fluid::Vector3{}
              && zeroTransport.diagnostics.kineticEnergyBeforeJoules == 0.0
              && zeroTransport.diagnostics.kineticEnergyAfterJoules == 0.0
              && zeroTransport.diagnostics
                     .maximumVelocityChangeMetersPerSecond == 0.0,
          "zero region momentum is an exact fixed point of transport and viscosity");
}

void testValidationAndLimits() {
    Fixture fixture(openScene());
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.xFaces(), 2.0);
    const auto openingFlux = fixture.flux(velocity);
    std::vector<double> warm(fixture.pressureOperator.rows.size(), 0.0);
    const auto accepted = fixture.project(
        velocity, openingFlux, warm, strictSettings());
    auto corrupt = accepted;
    corrupt.links.front()
        .predictedRelativeVolumeFlowRateCubicMetersPerSecond += 0.01;
    expectInvalid(
        [&] { validateSceneFluidPressureProjectionIntegrity(corrupt); },
        "pressure-projection integrity rejects link corruption");

    auto changedVelocity = velocity;
    changedVelocity.yFaces().front() = 0.5;
    expectInvalid(
        [&] { static_cast<void>(fixture.project(
            changedVelocity, openingFlux, warm, strictSettings())); },
        "pressure projection rejects an opening-flux velocity mismatch");

    auto invalidSettings = strictSettings();
    invalidSettings.densityKgPerCubicMeter = 0.0;
    expectInvalid(
        [&] { static_cast<void>(fixture.project(
            velocity, openingFlux, warm, invalidSettings)); },
        "pressure projection rejects non-positive density");

    SceneFluidPressureProjectionLimits limits;
    limits.maximumLinks = fixture.faceLinks.links.size() - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.project(
            velocity, openingFlux, warm, strictSettings(), limits)); },
        "pressure projection bounds link count");
    limits = {};
    limits.maximumProjectionBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.project(
            velocity, openingFlux, warm, strictSettings(), limits)); },
        "pressure projection bounds owned storage");
}

} // namespace

int main() {
    try {
        testClosedLinkProjection();
        testFaceAlignedOpeningProjection();
        testEmbeddedOpeningProjection();
        testConcaveEmbeddedOpeningProjection();
        testAuthoredNonPlanarEmbeddedOpeningProjection();
        testZeroFlowAndRejectedAttempt();
        testLinkResolvedContinuation();
        testAreaChangingLinkContinuation();
        testRegionMomentumReconstruction();
        testRegionMomentumTransport();
        testValidationAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid pressure-projection check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-projection checks passed");
    return 0;
}
