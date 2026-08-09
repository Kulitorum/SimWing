#include "scene_fluid_pressure_projection.h"
#include "scene_fluid_pressure_link_flow.h"
#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_region_momentum.h"

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
          volumes(buildSceneFluidCellVolumes(
              surface.definition, state, grid(), transfer, epoch)),
          connectivity(buildSceneFluidRegionConnectivity(
              surface.definition)),
          pressureVolumes(buildSceneFluidPressureControlVolumes(
              surface.definition, volumes, connectivity)),
          faceLinks(buildSceneFluidPressureFaceLinks(
              surface.definition, state, grid(), transfer, epoch, caps,
              openingQuadrature, openingPatches, volumes, connectivity,
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
        testZeroFlowAndRejectedAttempt();
        testLinkResolvedContinuation();
        testAreaChangingLinkContinuation();
        testRegionMomentumReconstruction();
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
