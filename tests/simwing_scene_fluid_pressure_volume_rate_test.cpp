#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_pressure_sampling.h"
#include "scene_fluid_mimetic_pressure_sampling.h"
#include "scene_fluid_mimetic_pressure_warm_start.h"
#include "scene_fluid_region_rebase.h"
#include "scene_fluid_region_wall.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string_view>
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

Scene tetraScene(const bool open = true) {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-pressure-volume-rate";
    scene.metadata.exporterVersion =
        "scene-fluid-pressure-volume-rate-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const double mouthX = open ? 2.0 : 2.2;
    const std::array<Vec3, 4> positions{{
        {1.2, 1.5, 1.45},
        {mouthX, 1.2, 1.15},
        {mouthX, 1.8, 1.15},
        {mouthX, 1.5, 1.75},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 3> faces{{
        {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}},
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
    if (open) {
        scene.openings = {
            {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
        };
    } else {
        const std::array<std::size_t, 3> mouth{{1, 2, 3}};
        scene.triangles.push_back({
            503, {11, 12, 13}, intrinsicChart(positions, mouth),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
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
    explicit Fixture(const bool open = true) : scene(tetraScene(open)) {}

    Scene scene;
    SceneFluidSurfaceAssembly surface = assembleSceneFluidSurface(scene);
    SceneStructureAssembly structureAssembly = fixedMouthAssembly(scene);
    Structure structure{structureAssembly.definition};
    SceneFluidSurfaceTransfer transfer{
        surface.definition, structureAssembly.mappings, structure};
    SceneFluidRegionConnectivity connectivity =
        buildSceneFluidRegionConnectivity(surface.definition);
};

struct Endpoint {
    SceneFluidSurfaceState state;
    SceneFluidGridEpoch epoch;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet openingQuadrature;
    SceneFluidOpeningGridPatchSet openingPatches;
    SceneFluidOpeningFaceCrossingSet openingFaceCrossings;
    SceneFluidCappedFacePartitionSet cappedFacePartitions;
    SceneFluidCellVolumeSet volumes;
    SceneFluidPressureControlVolumeSet pressureVolumes;
};

Endpoint captureEndpoint(const Fixture& fixture) {
    auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    auto epoch = buildSceneFluidGridEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer);
    auto caps = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    auto openingQuadrature = buildSceneFluidOpeningQuadrature(
        fixture.surface.definition, state, caps);
    auto openingPatches = buildSceneFluidOpeningGridPatches(
        fixture.surface.definition, state, caps, openingQuadrature, grid());
    auto openingFaceCrossings = buildSceneFluidOpeningFaceCrossings(
        fixture.surface.definition, state, caps, openingQuadrature,
        openingPatches, grid());
    auto cappedFacePartitions = buildSceneFluidCappedFacePartitions(
        fixture.surface.definition, state, grid(), fixture.transfer, epoch,
        caps, openingQuadrature, openingPatches, openingFaceCrossings);
    auto volumes = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, grid(), fixture.transfer, epoch);
    auto pressureVolumes = buildSceneFluidPressureControlVolumes(
        fixture.surface.definition, volumes, fixture.connectivity);
    return {std::move(state), std::move(epoch), std::move(caps),
            std::move(openingQuadrature), std::move(openingPatches),
            std::move(openingFaceCrossings),
            std::move(cappedFacePartitions), std::move(volumes),
            std::move(pressureVolumes)};
}

void advanceApex(Fixture& fixture, const double accelerationX) {
    const auto found = std::ranges::find(
        fixture.structureAssembly.mappings.nodeVertexIds,
        StableId{10});
    if (found == fixture.structureAssembly.mappings.nodeVertexIds.end()) {
        throw std::runtime_error("pressure-volume-rate apex is missing");
    }
    const std::size_t node = static_cast<std::size_t>(
        found - fixture.structureAssembly.mappings.nodeVertexIds.begin());
    const double mass = fixture.structure.definition().nodes[node].massKg;
    fixture.structure.addExternalForce(
        node, {mass * accelerationX, 0.0, 0.0});
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.25;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(settings);
    check(diagnostics.finite,
          "pressure-volume-rate structural step is finite");
}

void testTopologyStableMovingRates() {
    Fixture fixture;
    const auto previous = captureEndpoint(fixture);
    advanceApex(fixture, -0.8);
    const auto current = captureEndpoint(fixture);
    const auto first = buildSceneFluidPressureVolumeRates(
        previous.volumes, current.volumes, current.pressureVolumes);
    const auto repeated = buildSceneFluidPressureVolumeRates(
        previous.volumes, current.volumes, current.pressureVolumes);
    check(first == repeated
              && first.version == sceneFluidPressureVolumeRateVersion
              && first.fingerprint != 0
              && first.durationSeconds == 0.25
              && first.controlVolumes.size()
                  == current.pressureVolumes.controlVolumes.size()
              && first.components.size() == 1
              && first.maximumAbsoluteControlVolumeRateCubicMetersPerSecond
                  > 1.0e-4,
          "topology-stable moving pressure volumes publish deterministic rates");
    for (const auto& rate : first.controlVolumes) {
        checkNear(rate.volumeChangeCubicMeters,
                  rate.currentVolumeCubicMeters
                      - rate.previousVolumeCubicMeters,
                  0.0,
                  "pressure control-volume change retains exact endpoints");
        checkNear(rate.volumeChangeRateCubicMetersPerSecond,
                  rate.volumeChangeCubicMeters / first.durationSeconds,
                  0.0,
                  "pressure control-volume rate uses the exact epoch duration");
    }
    checkNear(first.globalVolumeChangeCubicMeters, 0.0, 3.0e-12,
              "moving sparse pressure volumes preserve fixed-domain volume");
    checkNear(first.components.front().volumeChangeCubicMeters,
              first.globalVolumeChangeCubicMeters, 3.0e-15,
              "open Cell and Outside share one globally closed component");
    validateSceneFluidPressureVolumeRates(
        first, previous.volumes, current.volumes,
        current.pressureVolumes);
}

void testStationaryAndTopologyRebase() {
    Fixture stationaryFixture;
    const auto stationaryPrevious = captureEndpoint(stationaryFixture);
    advanceApex(stationaryFixture, 0.0);
    const auto stationaryCurrent = captureEndpoint(stationaryFixture);
    const auto stationary = buildSceneFluidPressureVolumeRates(
        stationaryPrevious.volumes, stationaryCurrent.volumes,
        stationaryCurrent.pressureVolumes);
    check(stationary.maximumAbsoluteControlVolumeChangeCubicMeters == 0.0
              && stationary
                  .maximumAbsoluteControlVolumeRateCubicMetersPerSecond
                  == 0.0
              && stationary.globalVolumeChangeCubicMeters == 0.0,
          "stationary consecutive pressure volumes have exact zero rates");

    Fixture crossingFixture;
    const auto crossingPrevious = captureEndpoint(crossingFixture);
    advanceApex(crossingFixture, -8.0);
    const auto crossingCurrent = captureEndpoint(crossingFixture);
    check(crossingPrevious.volumes.cellRegionVolumes.size()
              != crossingCurrent.volumes.cellRegionVolumes.size(),
          "crossing fixture changes sparse cell-region topology");
    const auto crossing = buildSceneFluidPressureVolumeRates(
        crossingPrevious.volumes, crossingCurrent.volumes,
        crossingCurrent.pressureVolumes);
    check(crossing.previousControlVolumeCount == 65
              && crossing.retainedControlVolumeCount == 65
              && crossing.appearedControlVolumeCount == 1
              && crossing.controlVolumes.size() == 66,
          "pressure-volume rate accepts one explicitly marked control appearance");
    const auto appeared = std::ranges::find_if(
        crossing.controlVolumes,
        [](const auto& control) { return control.appearedThisEpoch; });
    check(appeared != crossing.controlVolumes.end()
              && appeared->cellIndex == 20
              && appeared->regionId == 2
              && appeared->previousVolumeCubicMeters == 0.0
              && appeared->currentVolumeCubicMeters > 0.0
              && appeared->volumeChangeCubicMeters
                  == appeared->currentVolumeCubicMeters,
          "new pressure control receives an exact zero-volume previous endpoint");
    checkNear(crossing.globalVolumeChangeCubicMeters, 0.0, 3.0e-12,
              "control appearance preserves the fixed-domain volume ledger");
    validateSceneFluidPressureVolumeRates(
        crossing, crossingPrevious.volumes, crossingCurrent.volumes,
        crossingCurrent.pressureVolumes);

    advanceApex(crossingFixture, 16.0);
    const auto retreatCurrent = captureEndpoint(crossingFixture);
    check(retreatCurrent.pressureVolumes.controlVolumes.size()
              < crossingCurrent.pressureVolumes.controlVolumes.size(),
          "reverse crossing fixture removes one sparse pressure control");
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidPressureVolumeRates(
            crossingCurrent.volumes, retreatCurrent.volumes,
            retreatCurrent.pressureVolumes)); },
        "pressure-volume rate still rejects a disappearing control");
}

void testMovingVolumeProjection() {
    Fixture fixture;
    const auto previous = captureEndpoint(fixture);
    advanceApex(fixture, -0.8);
    const auto current = captureEndpoint(fixture);
    const auto rates = buildSceneFluidPressureVolumeRates(
        previous.volumes, current.volumes, current.pressureVolumes);
    const auto faceLinks = buildSceneFluidPressureFaceLinks(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, current.openingFaceCrossings,
        current.cappedFacePartitions, current.volumes,
        fixture.connectivity, current.pressureVolumes);
    const auto pressureOperator = buildSceneFluidPressureOperator(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, current.volumes, fixture.connectivity,
        current.pressureVolumes, faceLinks);
    fluid::MacVelocityField predictedVelocity(grid());
    const auto openingFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, current.state, current.caps,
        current.openingQuadrature, current.openingPatches, grid(),
        predictedVelocity);
    SceneFluidPressureProjectionSettings settings;
    settings.timeStepSeconds = rates.durationSeconds;
    settings.absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond =
        2.0e-11;
    settings.relativeCorrectedVolumeRateTolerance = 1.0e-11;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-13;
    settings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-10;
    std::vector<double> warm(pressureOperator.rows.size(), 0.0);
    const auto projected = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, openingFlux, predictedVelocity,
        current.volumes, fixture.connectivity, current.pressureVolumes,
        faceLinks, pressureOperator, rates, warm, settings);
    check(projected.diagnostics.accepted
              && projected.diagnostics.usesMovingVolumeRates
              && projected.pressureVolumeRateFingerprint == rates.fingerprint
              && projected.diagnostics.pressureSolve.compatible
              && projected.diagnostics.pressureSolve.converged
              && projected.diagnostics
                  .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond
                  > 1.0e-4
              && projected.diagnostics
                  .predictedContinuityResidualMaximumCubicMetersPerSecond
                  > 1.0e-4,
          "moving-volume pressure projection consumes the nonzero dV/dt source");
    check(projected.diagnostics
              .correctedContinuityResidualMaximumCubicMetersPerSecond
              < 2.0e-11
              && projected.diagnostics
                  .maximumCorrectedComponentContinuityResidualCubicMetersPerSecond
                  < 2.0e-11,
          "moving-volume pressure correction closes local and component continuity");
    for (const auto& control : projected.controlVolumes) {
        checkNear(
            control.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond,
            -control.geometryVolumeChangeRateCubicMetersPerSecond,
            2.0e-11,
            "corrected link flow exactly opposes local geometry-volume rate");
    }
    double correctedApertureFlow = 0.0;
    for (const auto& link : projected.links) {
        if (link.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            correctedApertureFlow +=
                link.correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        }
    }
    check(correctedApertureFlow < -1.0e-4,
          "expanding cell pressure correction draws flow inward through its intake");

    const auto samples = sampleSceneFluidProjectedPressure(
        current.epoch.quadrature, current.pressureVolumes, projected);
    const auto repeatedSamples = sampleSceneFluidProjectedPressure(
        current.epoch.quadrature, current.pressureVolumes, projected);
    check(samples == repeatedSamples
              && samples.fingerprint != 0
              && samples.pressures.size()
                  == current.epoch.quadrature.points.size()
              && samples.maximumAbsolutePressureDifferencePascals > 1.0e-4,
          "accepted moving pressure samples deterministically on every material patch");
    for (std::size_t index = 0; index < samples.bindings.size(); ++index) {
        const auto& binding = samples.bindings[index];
        const auto& pressure = samples.pressures[index];
        check(binding.sampleIndex == index
                  && binding.stableId == pressure.stableId
                  && binding.componentIndex
                      == current.pressureVolumes.controlVolumes[
                             binding.negativeSideControlVolumeIndex]
                             .componentIndex
                  && binding.componentIndex
                      == current.pressureVolumes.controlVolumes[
                             binding.positiveSideControlVolumeIndex]
                             .componentIndex,
              "surface pressure sample binds both sides in one gauge component");
        checkNear(binding.pressureDifferencePascals,
                  pressure.negativeSidePressurePascals
                      - pressure.positiveSidePressurePascals,
                  0.0,
                  "surface pressure sample retains its exact one-sided difference");
    }
    const auto transferred = evaluateSceneFluidProjectedPressureQuadrature(
        fixture.surface.definition, current.state, fixture.transfer,
        current.epoch.quadrature, samples);
    check(transferred.diagnostics().finite
              && transferred.diagnostics().forceResidualNormNewtons
                  < 1.0e-12
              && transferred.diagnostics().momentResidualNormNewtonMeters
                  < 1.0e-12
              && std::hypot(
                  transferred.diagnostics().integratedSurfaceForceNewtons.x,
                  transferred.diagnostics().integratedSurfaceForceNewtons.y,
                  transferred.diagnostics().integratedSurfaceForceNewtons.z)
                  > 1.0e-5,
          "projected pressure reaches conservative nonzero structural traction");
    std::vector<double> shiftedWarm(pressureOperator.rows.size(), 73.0);
    const auto shiftedProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, openingFlux, predictedVelocity,
        current.volumes, fixture.connectivity, current.pressureVolumes,
        faceLinks, pressureOperator, rates, shiftedWarm, settings);
    const auto shiftedSamples = sampleSceneFluidProjectedPressure(
        current.epoch.quadrature, current.pressureVolumes,
        shiftedProjection);
    const auto shiftedTransfer =
        evaluateSceneFluidProjectedPressureQuadrature(
            fixture.surface.definition, current.state, fixture.transfer,
            current.epoch.quadrature, shiftedSamples);
    check(shiftedProjection.pressurePascals == projected.pressurePascals
              && shiftedSamples.pressures == samples.pressures
              && shiftedSamples.bindings == samples.bindings
              && shiftedTransfer == transferred,
          "uniform pressure warm-start gauge shifts leave sampled traction exact");
    fixture.transfer.addLoadsTo(fixture.structure, transferred);
    check(fixture.structure.diagnostics().pendingExternalForceNewtons
              == transferred.diagnostics().transferredNodalForceNewtons,
          "projected pressure load reaches the live Structure load accumulator");
    validateSceneFluidProjectedPressureSamples(
        samples, current.epoch.quadrature, current.pressureVolumes,
        projected);
    auto corruptSamples = samples;
    corruptSamples.bindings.front().pressureDifferencePascals += 1.0;
    expectInvalid(
        [&] { validateSceneFluidPressureSampleIntegrity(corruptSamples); },
        "projected-pressure sampling rejects nested ledger corruption");
    SceneFluidPressureSamplingLimits samplingLimits;
    samplingLimits.maximumSamples = samples.bindings.size() - 1;
    expectLimited(
        [&] { static_cast<void>(sampleSceneFluidProjectedPressure(
            current.epoch.quadrature, current.pressureVolumes, projected,
            samplingLimits)); },
        "projected-pressure sampling bounds sample count");
    samplingLimits = {};
    samplingLimits.maximumSamplingBytes = samples.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(sampleSceneFluidProjectedPressure(
            current.epoch.quadrature, current.pressureVolumes, projected,
            samplingLimits)); },
        "projected-pressure sampling bounds owned storage");

    auto truncatedSettings = settings;
    truncatedSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-30;
    truncatedSettings.pressureSolve.relativeResidualTolerance = 0.0;
    truncatedSettings.pressureSolve.maximumIterations = 0;
    const auto rejectedProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, openingFlux, predictedVelocity,
        current.volumes, fixture.connectivity, current.pressureVolumes,
        faceLinks, pressureOperator, rates, warm, truncatedSettings);
    expectInvalid(
        [&] { static_cast<void>(sampleSceneFluidProjectedPressure(
            current.epoch.quadrature, current.pressureVolumes,
            rejectedProjection)); },
        "pressure sampling rejects a non-converged projection attempt");

    const auto frozen = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, openingFlux, predictedVelocity,
        current.volumes, fixture.connectivity, current.pressureVolumes,
        faceLinks, pressureOperator, warm, settings);
    check(frozen.diagnostics.accepted
              && !frozen.diagnostics.usesMovingVolumeRates
              && frozen.diagnostics
                  .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                  == 0.0
              && std::ranges::all_of(
                  frozen.pressurePascals,
                  [](const double pressure) { return pressure == 0.0; }),
          "frozen projection remains an exact zero-flow baseline without dV/dt");

    auto mismatchedSettings = settings;
    mismatchedSettings.timeStepSeconds *= 0.5;
    expectInvalid(
        [&] { static_cast<void>(projectSceneFluidPressureLinkFlows(
            fixture.surface.definition, current.state, grid(),
            fixture.transfer, current.epoch, current.caps,
            current.openingQuadrature, current.openingPatches, openingFlux,
            predictedVelocity, current.volumes, fixture.connectivity,
            current.pressureVolumes, faceLinks, pressureOperator, rates,
            warm, mismatchedSettings)); },
        "moving-volume projection rejects a pressure-step duration mismatch");
}

void testAppearedControlRegionRebase() {
    Fixture fixture;
    const auto previousState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto previousEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, previousState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto previousMimeticControls =
        buildSceneFluidMimeticControlCells(
            fixture.surface.definition, previousState, grid(),
            previousEpoch.gridEpoch, previousEpoch.openingCaps,
            previousEpoch.openingQuadrature,
            previousEpoch.openingPatches,
            previousEpoch.pressureControlVolumes,
            previousEpoch.pressureFaceLinks);
    const auto previousMimeticFull =
        buildSceneFluidMimeticTraceSystem(previousMimeticControls);
    const auto previousMimeticCondensed =
        buildSceneFluidMimeticCondensedTraceSystem(previousMimeticFull);
    std::vector<double> previousPredictedVolumeRates(
        previousMimeticControls.controlCells.size(), 0.0);
    std::size_t mimeticReceiver = 1;
    while (mimeticReceiver < previousMimeticControls.controlCells.size()
           && previousMimeticControls.controlCells[mimeticReceiver]
                  .componentIndex
               != previousMimeticControls.controlCells.front()
                  .componentIndex) {
        ++mimeticReceiver;
    }
    check(mimeticReceiver < previousPredictedVolumeRates.size(),
          "mimetic warm-start fixture finds a balanced source pair");
    if (mimeticReceiver < previousPredictedVolumeRates.size()) {
        previousPredictedVolumeRates.front() = 0.001;
        previousPredictedVolumeRates[mimeticReceiver] = -0.001;
    }
    SceneFluidMimeticPressureSourceSettings mimeticSourceSettings;
    mimeticSourceSettings.densityKgPerCubicMeter = 1.0;
    mimeticSourceSettings.timeStepSeconds = 0.25;
    const auto previousMimeticSources =
        buildSceneFluidMimeticPressureSources(
            previousMimeticControls, previousPredictedVolumeRates,
            mimeticSourceSettings);
    SceneFluidMimeticTraceSolveSettings mimeticSolveSettings;
    mimeticSolveSettings.absoluteResidualTolerancePascalsMeters = 1.0e-12;
    mimeticSolveSettings.relativeResidualTolerance = 1.0e-13;
    mimeticSolveSettings
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-11;
    mimeticSolveSettings.maximumIterations = 4000;
    const std::vector<double> zeroPreviousMimeticWarm(
        previousMimeticCondensed.traces.size(), 0.0);
    const auto previousMimeticPressure =
        solveSceneFluidMimeticPressureSystem(
            previousMimeticCondensed, previousMimeticFull,
            previousMimeticSources, zeroPreviousMimeticWarm,
            mimeticSolveSettings);
    check(previousMimeticPressure.diagnostics.accepted,
          "mimetic warm-start source pressure solve is accepted");
    const auto previousMimeticState =
        captureSceneFluidMimeticPressureState(
            previousMimeticControls, previousMimeticFull,
            previousMimeticCondensed, previousMimeticSources,
            previousMimeticPressure);
    const auto previousMimeticSamples =
        sampleSceneFluidMimeticPressure(
            previousEpoch.gridEpoch.quadrature,
            previousEpoch.pressureControlVolumes,
            previousMimeticControls, previousMimeticFull,
            previousMimeticCondensed, previousMimeticState);
    const auto repeatedPreviousMimeticSamples =
        sampleSceneFluidMimeticPressure(
            previousEpoch.gridEpoch.quadrature,
            previousEpoch.pressureControlVolumes,
            previousMimeticControls, previousMimeticFull,
            previousMimeticCondensed, previousMimeticState);
    bool exactMimeticSamples = previousMimeticSamples.bindings.size()
        == previousEpoch.gridEpoch.quadrature.points.size();
    for (const auto& binding : previousMimeticSamples.bindings) {
        const auto& pressure = previousMimeticSamples.pressures[
            binding.sampleIndex];
        exactMimeticSamples = exactMimeticSamples
            && pressure.negativeSidePressurePascals
                == previousMimeticState.controls[
                    binding.negativeSideControlVolumeIndex]
                      .pressurePascals
            && pressure.positiveSidePressurePascals
                == previousMimeticState.controls[
                    binding.positiveSideControlVolumeIndex]
                      .pressurePascals
            && binding.pressureDifferencePascals
                == pressure.negativeSidePressurePascals
                    - pressure.positiveSidePressurePascals;
    }
    check(previousMimeticSamples
                  == repeatedPreviousMimeticSamples
              && previousMimeticSamples.fingerprint != 0
              && previousMimeticSamples.pressureStateFingerprint
                  == previousMimeticState.fingerprint
              && previousMimeticSamples.mimeticControlCellFingerprint
                  == previousMimeticControls.fingerprint
              && exactMimeticSamples,
          "accepted mimetic control pressures sample both material sides with exact gauge-safe differences");
    const auto previousMimeticTransfer =
        evaluateSceneFluidMimeticPressureQuadrature(
            fixture.surface.definition, previousState, fixture.transfer,
            previousEpoch.gridEpoch.quadrature,
            previousMimeticSamples);
    check(previousMimeticTransfer.diagnostics().finite
              && previousMimeticTransfer.diagnostics()
                     .forceResidualNormNewtons < 1.0e-12
              && previousMimeticTransfer.diagnostics()
                     .momentResidualNormNewtonMeters < 1.0e-12,
          "accepted mimetic pressure uses the conservative Structure traction path");
    validateSceneFluidMimeticPressureSamples(
        previousMimeticSamples,
        previousEpoch.gridEpoch.quadrature,
        previousEpoch.pressureControlVolumes,
        previousMimeticControls, previousMimeticFull,
        previousMimeticCondensed, previousMimeticState);
    auto corruptPreviousMimeticSamples = previousMimeticSamples;
    corruptPreviousMimeticSamples.bindings.front()
        .pressureDifferencePascals += 0.01;
    expectInvalid(
        [&] { validateSceneFluidMimeticPressureSampleIntegrity(
            corruptPreviousMimeticSamples); },
        "mimetic pressure sampling rejects nested ledger corruption");
    SceneFluidPressureSamplingLimits mimeticSamplingLimits;
    mimeticSamplingLimits.maximumSamples =
        previousMimeticSamples.bindings.size() - 1;
    expectLimited(
        [&] { static_cast<void>(sampleSceneFluidMimeticPressure(
            previousEpoch.gridEpoch.quadrature,
            previousEpoch.pressureControlVolumes,
            previousMimeticControls, previousMimeticFull,
            previousMimeticCondensed, previousMimeticState,
            mimeticSamplingLimits)); },
        "mimetic pressure sampling bounds sample count");
    mimeticSamplingLimits = {};
    mimeticSamplingLimits.maximumSamplingBytes =
        previousMimeticSamples.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(sampleSceneFluidMimeticPressure(
            previousEpoch.gridEpoch.quadrature,
            previousEpoch.pressureControlVolumes,
            previousMimeticControls, previousMimeticFull,
            previousMimeticCondensed, previousMimeticState,
            mimeticSamplingLimits)); },
        "mimetic pressure sampling bounds owned storage");
    fluid::MacVelocityField velocity(grid());
    std::ranges::fill(velocity.yFaces(), 0.4);
    std::ranges::fill(velocity.zFaces(), -0.2);
    const auto previousFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, previousState,
        previousEpoch.openingCaps, previousEpoch.openingQuadrature,
        previousEpoch.openingPatches, grid(), velocity);
    SceneFluidPressureProjectionSettings projectionSettings;
    projectionSettings.timeStepSeconds = 0.25;
    std::vector<double> warm(
        previousEpoch.pressureOperator.rows.size(), 0.0);
    const auto previousProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, previousState, grid(), fixture.transfer,
        previousEpoch.gridEpoch, previousEpoch.openingCaps,
        previousEpoch.openingQuadrature, previousEpoch.openingPatches,
        previousFlux, velocity, previousEpoch.cellVolumes,
        fixture.connectivity, previousEpoch.pressureControlVolumes,
        previousEpoch.pressureFaceLinks, previousEpoch.pressureOperator,
        warm, projectionSettings);
    check(previousProjection.diagnostics.accepted,
          "region-rebase source pressure projection is accepted");
    const auto previousMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), previousEpoch.pressureControlVolumes,
        previousEpoch.pressureFaceLinks, previousEpoch.openingPatches,
        previousProjection, velocity);
    SceneFluidRegionTransportSettings transportSettings;
    transportSettings.timeStepSeconds = 0.25;
    const auto transport = advanceSceneFluidRegionMomentum(
        previousMomentum, previousEpoch.pressureFaceLinks,
        previousProjection, transportSettings);
    check(transport.diagnostics.accepted,
          "region-rebase source momentum transport is accepted");

    advanceApex(fixture, -8.0);
    const auto currentState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto currentEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto topologyTransition =
        buildSceneFluidPressureTopologyTransition(
            previousEpoch.pressureControlVolumes,
            previousEpoch.pressureFaceLinks,
            currentEpoch.pressureControlVolumes,
            currentEpoch.pressureFaceLinks);
    const auto repeatedTransition =
        buildSceneFluidPressureTopologyTransition(
            previousEpoch.pressureControlVolumes,
            previousEpoch.pressureFaceLinks,
            currentEpoch.pressureControlVolumes,
            currentEpoch.pressureFaceLinks);
    check(topologyTransition == repeatedTransition
              && topologyTransition.retainedControlVolumeCount == 65
              && topologyTransition.appearedControlVolumeCount == 1
              && topologyTransition.disappearedControlVolumeCount == 0
              && !topologyTransition.appearanceDonors.empty()
              && topologyTransition.retirementRecipients.empty()
              && topologyTransition.maximumAppearanceDonorCount > 0,
          "pressure topology transition deterministically owns the appeared row and its donors");
    const auto currentMimeticControls =
        buildSceneFluidMimeticControlCells(
            fixture.surface.definition, currentState, grid(),
            currentEpoch.gridEpoch, currentEpoch.openingCaps,
            currentEpoch.openingQuadrature, currentEpoch.openingPatches,
            currentEpoch.pressureControlVolumes,
            currentEpoch.pressureFaceLinks);
    const auto currentMimeticFull =
        buildSceneFluidMimeticTraceSystem(currentMimeticControls);
    const auto currentMimeticCondensed =
        buildSceneFluidMimeticCondensedTraceSystem(currentMimeticFull);
    const auto mimeticWarm = buildSceneFluidMimeticPressureWarmStart(
        previousMimeticState, previousMimeticControls,
        previousMimeticFull, previousMimeticCondensed,
        currentMimeticControls, currentMimeticFull,
        currentMimeticCondensed, topologyTransition);
    const auto repeatedMimeticWarm =
        buildSceneFluidMimeticPressureWarmStart(
            previousMimeticState, previousMimeticControls,
            previousMimeticFull, previousMimeticCondensed,
            currentMimeticControls, currentMimeticFull,
            currentMimeticCondensed, topologyTransition);
    std::vector<double> previousMimeticControlPressures;
    previousMimeticControlPressures.reserve(
        previousMimeticState.controls.size());
    for (const auto& control : previousMimeticState.controls) {
        previousMimeticControlPressures.push_back(
            control.pressurePascals);
    }
    const auto rebasedMimeticControlPressures =
        rebaseSceneFluidPressureWarmStart(
            previousEpoch.pressureControlVolumes,
            currentEpoch.pressureControlVolumes, topologyTransition,
            previousMimeticControlPressures);
    bool exactMimeticTracePolicy = true;
    std::size_t observedRetainedMimeticTraces = 0;
    std::size_t observedAppearedMimeticTraces = 0;
    for (const auto& trace : currentMimeticCondensed.traces) {
        const auto previousTrace = std::ranges::find(
            previousMimeticState.traces, trace.stableId,
            &SceneFluidMimeticAcceptedTracePressure::stableId);
        double expectedPressure = 0.0;
        if (previousTrace != previousMimeticState.traces.end()) {
            expectedPressure = previousTrace->pressurePascals;
            ++observedRetainedMimeticTraces;
        } else {
            const auto& fullTrace =
                currentMimeticFull.traces[trace.fullTraceIndex];
            for (std::size_t offset = 0;
                 offset < fullTrace.incidenceCount; ++offset) {
                const auto& incidence = currentMimeticFull.incidences[
                    fullTrace.firstIncidence + offset];
                expectedPressure += rebasedMimeticControlPressures[
                    incidence.controlCellIndex];
            }
            expectedPressure /= static_cast<double>(
                fullTrace.incidenceCount);
            ++observedAppearedMimeticTraces;
        }
        expectedPressure -= mimeticWarm.componentGaugeShiftsPascals[
            trace.componentIndex];
        if (trace.isGauge) {
            expectedPressure = 0.0;
        }
        exactMimeticTracePolicy = exactMimeticTracePolicy
            && mimeticWarm.reducedTracePascals[trace.traceIndex]
                == expectedPressure;
    }
    check(mimeticWarm == repeatedMimeticWarm
              && mimeticWarm.fingerprint != 0
              && mimeticWarm.sourcePressureStateFingerprint
                  == previousMimeticState.fingerprint
              && mimeticWarm.sourceTopologyTransitionFingerprint
                  == topologyTransition.fingerprint
              && mimeticWarm.retainedTraceCount
                  == observedRetainedMimeticTraces
              && mimeticWarm.appearedTraceCount
                  == observedAppearedMimeticTraces
              && mimeticWarm.appearedTraceCount > 0
              && mimeticWarm.disappearedTraceCount
                  == previousMimeticState.traces.size()
                      - observedRetainedMimeticTraces
              && exactMimeticTracePolicy,
          "mimetic pressure warm start preserves retained traces, initializes appearances from rebased endpoints, and normalizes current gauges");
    bool exactMimeticGauges = true;
    for (const std::size_t gauge :
         currentMimeticCondensed.componentGaugeTraceIndices) {
        exactMimeticGauges = exactMimeticGauges
            && mimeticWarm.reducedTracePascals[gauge] == 0.0;
    }
    check(exactMimeticGauges,
          "mimetic pressure warm start fixes every current gauge exactly");
    validateSceneFluidMimeticPressureWarmStart(
        mimeticWarm, previousMimeticState, previousMimeticControls,
        previousMimeticFull, previousMimeticCondensed,
        currentMimeticControls, currentMimeticFull,
        currentMimeticCondensed, topologyTransition);
    std::vector<double> zeroCurrentPredictedVolumeRates(
        currentMimeticControls.controlCells.size(), 0.0);
    const auto zeroCurrentMimeticSources =
        buildSceneFluidMimeticPressureSources(
            currentMimeticControls, zeroCurrentPredictedVolumeRates,
            mimeticSourceSettings);
    const auto warmStartedCurrentPressure =
        solveSceneFluidMimeticPressureSystem(
            currentMimeticCondensed, currentMimeticFull,
            zeroCurrentMimeticSources, mimeticWarm.reducedTracePascals,
            mimeticSolveSettings);
    check(warmStartedCurrentPressure.diagnostics.accepted,
          "atomic mimetic pressure solve consumes the consecutive-epoch warm start");
    auto corruptMimeticWarm = mimeticWarm;
    corruptMimeticWarm.reducedTracePascals.back() += 0.01;
    expectInvalid(
        [&] { validateSceneFluidMimeticPressureWarmStartIntegrity(
            corruptMimeticWarm); },
        "mimetic pressure warm-start integrity rejects trace corruption");
    SceneFluidMimeticPressureWarmStartLimits mimeticWarmLimits;
    mimeticWarmLimits.maximumControlCells =
        currentMimeticControls.controlCells.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticPressureWarmStart(
                previousMimeticState, previousMimeticControls,
                previousMimeticFull, previousMimeticCondensed,
                currentMimeticControls, currentMimeticFull,
                currentMimeticCondensed, topologyTransition,
                mimeticWarmLimits)); },
        "mimetic pressure warm start bounds control count");
    mimeticWarmLimits = {};
    mimeticWarmLimits.maximumOwnedBytes =
        mimeticWarm.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticPressureWarmStart(
                previousMimeticState, previousMimeticControls,
                previousMimeticFull, previousMimeticCondensed,
                currentMimeticControls, currentMimeticFull,
                currentMimeticCondensed, topologyTransition,
                mimeticWarmLimits)); },
        "mimetic pressure warm start bounds owned storage");
    mimeticWarmLimits = {};
    mimeticWarmLimits.maximumWorkingBytes =
        mimeticWarm.workingStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticPressureWarmStart(
                previousMimeticState, previousMimeticControls,
                previousMimeticFull, previousMimeticCondensed,
                currentMimeticControls, currentMimeticFull,
                currentMimeticCondensed, topologyTransition,
                mimeticWarmLimits)); },
        "mimetic pressure warm start bounds complete working storage");
    validateSceneFluidPressureTopologyTransition(
        topologyTransition,
        previousEpoch.pressureControlVolumes,
        previousEpoch.pressureFaceLinks,
        currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks);
    auto corruptTransition = topologyTransition;
    corruptTransition.appearanceDonors.front().normalizedWeight += 0.01;
    expectInvalid(
        [&] {
            validateSceneFluidPressureTopologyTransitionIntegrity(
                corruptTransition);
        },
        "pressure topology-transition integrity rejects donor corruption");
    SceneFluidPressureTopologyTransitionLimits transitionLimits;
    transitionLimits.maximumMappings =
        topologyTransition.appearanceDonors.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidPressureTopologyTransition(
                previousEpoch.pressureControlVolumes,
                previousEpoch.pressureFaceLinks,
                currentEpoch.pressureControlVolumes,
                currentEpoch.pressureFaceLinks,
                transitionLimits)); },
        "pressure topology transition bounds its mapping count");
    transitionLimits = {};
    transitionLimits.maximumTransitionBytes =
        topologyTransition.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidPressureTopologyTransition(
                previousEpoch.pressureControlVolumes,
                previousEpoch.pressureFaceLinks,
                currentEpoch.pressureControlVolumes,
                currentEpoch.pressureFaceLinks,
                transitionLimits)); },
        "pressure topology transition bounds its owned storage");
    const auto rebase = rebaseSceneFluidRegionTransport(
        transport, previousEpoch.pressureControlVolumes,
        currentEpoch.pressureControlVolumes,
        topologyTransition);
    const auto repeated = rebaseSceneFluidRegionTransport(
        transport, previousEpoch.pressureControlVolumes,
        currentEpoch.pressureControlVolumes,
        topologyTransition);
    check(rebase == repeated
              && rebase.diagnostics.previousControlVolumeCount == 65
              && rebase.diagnostics.currentControlVolumeCount == 66
              && rebase.diagnostics.retainedControlVolumeCount == 65
              && rebase.diagnostics.appearedControlVolumeCount == 1
              && rebase.controlVolumes.size() == 66,
          "region transport deterministically rebases one appeared control");
    const auto appeared = std::ranges::find_if(
        rebase.controlVolumes,
        [](const auto& control) { return control.appearedThisEpoch; });
    check(appeared != rebase.controlVolumes.end()
              && appeared->cellIndex == 20
              && appeared->regionId == 2
              && appeared->sourceControlVolumeIndex
                  == std::numeric_limits<std::size_t>::max()
              && appeared->donorControlVolumeCount > 0
              && appeared->donorLinkAreaSquareMeters > 0.0
              && appeared->velocityMetersPerSecond
                  == fluid::Vector3{0.0, 0.4, -0.2},
          "appeared control inherits only retained same-region one-ring velocity");
    bool retainedVelocityExact = true;
    for (const auto& control : rebase.controlVolumes) {
        if (control.appearedThisEpoch) {
            continue;
        }
        retainedVelocityExact = retainedVelocityExact
            && control.sourceControlVolumeIndex
                < transport.controlVolumes.size()
            && control.stableId
                == transport.controlVolumes[control.sourceControlVolumeIndex]
                    .stableId
            && control.velocityMetersPerSecond
                == transport.controlVolumes[control.sourceControlVolumeIndex]
                    .velocityMetersPerSecond;
    }
    check(retainedVelocityExact,
          "region rebase preserves every retained stable-ID velocity across index insertion");
    validateSceneFluidRegionRebase(
        rebase, transport, previousEpoch.pressureControlVolumes,
        currentEpoch.pressureControlVolumes,
        topologyTransition);
    std::vector<double> previousWarm(
        previousEpoch.pressureControlVolumes.controlVolumes.size(), 73.0);
    const auto currentWarm = rebaseSceneFluidPressureWarmStart(
        previousEpoch.pressureControlVolumes,
        currentEpoch.pressureControlVolumes,
        topologyTransition, previousWarm);
    check(currentWarm.size()
              == currentEpoch.pressureControlVolumes.controlVolumes.size()
              && std::ranges::all_of(
                  currentWarm,
                  [](const double pressure) { return pressure == 73.0; }),
          "pressure warm-start rebase preserves a uniform gauge through control appearance");
    std::ranges::fill(previousWarm, 0.0);
    const auto donorPressure = std::ranges::find_if(
        previousEpoch.pressureControlVolumes.controlVolumes,
        [](const auto& control) {
            return control.cellIndex == 21 && control.regionId == 2;
        });
    check(donorPressure
              != previousEpoch.pressureControlVolumes.controlVolumes.end(),
          "pressure warm-start fixture finds the retained same-region donor");
    if (donorPressure
        != previousEpoch.pressureControlVolumes.controlVolumes.end()) {
        previousWarm[donorPressure->controlVolumeIndex] = 17.0;
        const auto nonuniformWarm = rebaseSceneFluidPressureWarmStart(
            previousEpoch.pressureControlVolumes,
            currentEpoch.pressureControlVolumes,
            topologyTransition, previousWarm);
        check(nonuniformWarm[appeared->controlVolumeIndex] == 17.0,
              "appeared pressure warm start uses its exact retained same-region donor");
    }
    SceneFluidRegionWallSettings wallSettings;
    wallSettings.timeStepSeconds = 0.25;
    const auto wall = exchangeSceneFluidRegionWallMomentum(
        rebase, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature, wallSettings);
    check(wall.diagnostics.accepted
              && wall.sourceTransportFingerprint == transport.fingerprint
              && wall.sourceRegionRebaseFingerprint == rebase.fingerprint
              && wall.controlVolumes.size()
                  == currentEpoch.pressureControlVolumes.controlVolumes.size(),
          "material-wall exchange consumes the rebased current topology with explicit provenance");
    validateSceneFluidRegionWallExchange(
        wall, rebase, grid(), currentEpoch.pressureControlVolumes,
        fixture.surface.definition, currentState,
        currentEpoch.gridEpoch.quadrature);

    auto corrupt = rebase;
    corrupt.controlVolumes.back().momentumKilogramMetersPerSecond.x += 1.0;
    expectInvalid(
        [&] { validateSceneFluidRegionRebaseIntegrity(corrupt); },
        "region-rebase integrity rejects momentum corruption");
    auto corruptLinks = currentEpoch.pressureFaceLinks;
    corruptLinks.links.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { static_cast<void>(
            buildSceneFluidPressureTopologyTransition(
                previousEpoch.pressureControlVolumes,
                previousEpoch.pressureFaceLinks,
                currentEpoch.pressureControlVolumes, corruptLinks)); },
        "topology transition rejects nested face-link corruption");
    SceneFluidRegionRebaseLimits limits;
    limits.maximumControlVolumes = rebase.controlVolumes.size() - 1;
    expectLimited(
        [&] { static_cast<void>(rebaseSceneFluidRegionTransport(
            transport, previousEpoch.pressureControlVolumes,
            currentEpoch.pressureControlVolumes,
            topologyTransition, limits)); },
        "region rebase bounds its current control count");

    const auto crossingFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, currentState,
        currentEpoch.openingCaps, currentEpoch.openingQuadrature,
        currentEpoch.openingPatches, grid(), velocity);
    std::vector<double> crossingWarm(
        currentEpoch.pressureOperator.rows.size(), 0.0);
    const auto crossingProjection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, currentState, grid(), fixture.transfer,
        currentEpoch.gridEpoch, currentEpoch.openingCaps,
        currentEpoch.openingQuadrature, currentEpoch.openingPatches,
        crossingFlux, velocity, currentEpoch.cellVolumes,
        fixture.connectivity, currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.pressureOperator,
        crossingWarm, projectionSettings);
    check(crossingProjection.diagnostics.accepted,
          "region-retirement source pressure projection is accepted");
    const auto crossingMomentum = reconstructSceneFluidRegionMomentumState(
        grid(), currentEpoch.pressureControlVolumes,
        currentEpoch.pressureFaceLinks, currentEpoch.openingPatches,
        crossingProjection, velocity);
    const auto crossingTransport = advanceSceneFluidRegionMomentum(
        crossingMomentum, currentEpoch.pressureFaceLinks,
        crossingProjection, transportSettings);
    check(crossingTransport.diagnostics.accepted,
          "region-retirement source momentum transport is accepted");
    advanceApex(fixture, 16.0);
    const auto retreatState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto retreatEpoch = buildSceneFluidPressureEpoch(
        fixture.surface.definition, retreatState, grid(), fixture.transfer,
        fixture.connectivity);
    const auto retirementTransition =
        buildSceneFluidPressureTopologyTransition(
            currentEpoch.pressureControlVolumes,
            currentEpoch.pressureFaceLinks,
            retreatEpoch.pressureControlVolumes,
            retreatEpoch.pressureFaceLinks);
    check(retirementTransition.retainedControlVolumeCount == 65
              && retirementTransition.appearedControlVolumeCount == 0
              && retirementTransition.disappearedControlVolumeCount == 1
              && retirementTransition.appearanceDonors.empty()
              && retirementTransition.retirementRecipients.size() == 1
              && retirementTransition.retirementRecipients.front()
                     .normalizedWeight == 1.0,
          "pressure topology transition owns the unique disappeared-row retirement");
    const auto retirementRates = buildSceneFluidPressureVolumeRates(
        currentEpoch.cellVolumes, retreatEpoch.cellVolumes,
        retreatEpoch.pressureControlVolumes, retirementTransition);
    check(retirementRates.previousControlVolumeCount == 66
              && retirementRates.retainedControlVolumeCount == 65
              && retirementRates.appearedControlVolumeCount == 0
              && retirementRates.disappearedControlVolumeCount == 1
              && retirementRates.retiredPreviousVolumeCubicMeters > 0.0
              && retirementRates.pressureTopologyTransitionFingerprint
                  == retirementTransition.fingerprint
              && std::abs(retirementRates.globalVolumeChangeCubicMeters)
                  < 3.0e-12,
          "pressure-volume rates conservatively retire one disappeared row");
    const auto retirementRateRecipient = std::ranges::find_if(
        retirementRates.controlVolumes,
        [](const auto& control) {
            return control.cellIndex == 21 && control.regionId == 2;
        });
    check(retirementRateRecipient
                  != retirementRates.controlVolumes.end()
              && retirementRateRecipient
                     ->retiredPreviousVolumeCubicMeters > 0.0
              && retirementRateRecipient->previousVolumeCubicMeters
                  > retirementRateRecipient
                        ->retiredPreviousVolumeCubicMeters,
          "retained pressure row owns the retired previous volume endpoint");
    validateSceneFluidPressureVolumeRates(
        retirementRates, currentEpoch.cellVolumes,
        retreatEpoch.cellVolumes,
        retreatEpoch.pressureControlVolumes, retirementTransition);
    const auto retirement = rebaseSceneFluidRegionTransport(
        crossingTransport, currentEpoch.pressureControlVolumes,
        retreatEpoch.pressureControlVolumes,
        retirementTransition);
    check(retirement.diagnostics.previousControlVolumeCount == 66
              && retirement.diagnostics.currentControlVolumeCount == 65
              && retirement.diagnostics.retainedControlVolumeCount == 65
              && retirement.diagnostics.appearedControlVolumeCount == 0
              && retirement.diagnostics.disappearedControlVolumeCount == 1
              && retirement.sourceTopologyTransitionFingerprint
                  == retirementTransition.fingerprint
              && retirement.diagnostics
                     .maximumRetirementRecipientCount > 0
              && std::abs(retirement.diagnostics
                     .sourceVolumeMappingResidualCubicMeters) < 1.0e-12
              && retirement.diagnostics
                     .sourceMomentumMappingResidualNormKilogramMetersPerSecond
                  < 1.0e-12,
          "region transport conservatively retires one disappeared control "
          "onto retained same-region recipients");
    const auto recipient = std::ranges::find_if(
        retirement.controlVolumes,
        [](const auto& control) {
            return control.cellIndex == 21 && control.regionId == 2;
        });
    check(recipient != retirement.controlVolumes.end()
              && recipient->retiredSourceControlVolumeCount == 1
              && recipient->retiredSourceLinkAreaSquareMeters > 0.0
              && recipient->mappedSourceVolumeCubicMeters
                  > crossingTransport.controlVolumes[
                        recipient->sourceControlVolumeIndex]
                        .volumeCubicMeters,
          "retained same-region neighbour owns the retired source volume and momentum");
    validateSceneFluidRegionRebase(
        retirement, crossingTransport,
        currentEpoch.pressureControlVolumes,
        retreatEpoch.pressureControlVolumes,
        retirementTransition);
    std::vector<double> retirementWarm(
        currentEpoch.pressureControlVolumes.controlVolumes.size(), 0.0);
    const auto retainedPressure = std::ranges::find_if(
        currentEpoch.pressureControlVolumes.controlVolumes,
        [](const auto& control) {
            return control.cellIndex == 21 && control.regionId == 2;
        });
    const auto disappearedPressure = std::ranges::find_if(
        currentEpoch.pressureControlVolumes.controlVolumes,
        [](const auto& control) {
            return control.cellIndex == 20 && control.regionId == 2;
        });
    check(retainedPressure
                  != currentEpoch.pressureControlVolumes.controlVolumes.end()
              && disappearedPressure
                  != currentEpoch.pressureControlVolumes.controlVolumes.end(),
          "pressure-retirement fixture finds retained and disappeared rows");
    if (retainedPressure
            != currentEpoch.pressureControlVolumes.controlVolumes.end()
        && disappearedPressure
            != currentEpoch.pressureControlVolumes.controlVolumes.end()) {
        retirementWarm[retainedPressure->controlVolumeIndex] = 17.0;
        retirementWarm[disappearedPressure->controlVolumeIndex] = 99.0;
        const auto retiredWarm = rebaseSceneFluidPressureWarmStart(
            currentEpoch.pressureControlVolumes,
            retreatEpoch.pressureControlVolumes,
            retirementTransition, retirementWarm);
        check(retiredWarm[recipient->controlVolumeIndex] == 17.0,
              "pressure warm-start retirement preserves the retained row "
              "and drops the retired value");
    }
}

void testIndependentGaugeSamplingRejection() {
    Fixture fixture(false);
    const auto current = captureEndpoint(fixture);
    const auto faceLinks = buildSceneFluidPressureFaceLinks(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, current.openingFaceCrossings,
        current.cappedFacePartitions, current.volumes,
        fixture.connectivity, current.pressureVolumes);
    const auto pressureOperator = buildSceneFluidPressureOperator(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, current.volumes, fixture.connectivity,
        current.pressureVolumes, faceLinks);
    fluid::MacVelocityField velocity(grid());
    const auto openingFlux = evaluateSceneFluidOpeningFlux(
        fixture.surface.definition, current.state, current.caps,
        current.openingQuadrature, current.openingPatches, grid(), velocity);
    std::vector<double> warm(pressureOperator.rows.size(), 0.0);
    const auto projection = projectSceneFluidPressureLinkFlows(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, openingFlux, velocity, current.volumes,
        fixture.connectivity, current.pressureVolumes, faceLinks,
        pressureOperator, warm);
    check(projection.diagnostics.accepted
              && pressureOperator.components.size() == 2,
          "sealed tetra has an accepted zero-flow solve with independent gauges");
    expectInvalid(
        [&] { static_cast<void>(sampleSceneFluidProjectedPressure(
            current.epoch.quadrature, current.pressureVolumes,
            projection)); },
        "surface sampling rejects pressure differences across independent gauges");
}

void testComposedPressureEpoch() {
    Fixture fixture;
    const auto state = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    const auto first = buildSceneFluidPressureEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer,
        fixture.connectivity);
    const auto repeated = buildSceneFluidPressureEpoch(
        fixture.surface.definition, state, grid(), fixture.transfer,
        fixture.connectivity);
    check(first == repeated
              && first.version == sceneFluidPressureEpochVersion
              && first.fingerprint != 0
              && first.surfaceStateFingerprint == state.fingerprint
              && first.gridEpoch.fingerprint != 0
              && first.openingCaps.caps.size() == 1
              && first.openingPatches.patches.size() == 1
              && first.openingFaceCrossings.fingerprint != 0
              && first.openingFaceCrossings.crossings.empty()
              && first.openingFaceCrossings.faceOwnedPatchCount == 1
              && first.cappedFacePartitions.fingerprint != 0
              && first.cappedFacePartitions.touchedFaceCount == 0
              && first.cappedFacePartitions.partitions.empty()
              && first.cellVolumes.cellRegionVolumes.size()
                  == first.pressureControlVolumes.controlVolumes.size()
              && first.pressureFaceLinks.unresolvedActiveFaceCount == 0
              && first.pressureFaceLinks.unresolvedCappedFaceCount == 0
              && first.pressureFaceLinks.unresolvedAmbiguousFaceCount == 0
              && first.pressureFaceLinks.unresolvedOpeningFaceCount == 0
              && first.pressureFaceLinks.unresolvedEmbeddedOpeningPatchCount
                  == 0
              && first.pressureOperator.rows.size()
                  == first.pressureControlVolumes.controlVolumes.size(),
          "pressure epoch atomically composes one fully resolved accepted topology");
    validateSceneFluidPressureEpoch(
        first, fixture.surface.definition, state, grid(), fixture.transfer,
        fixture.connectivity);

    auto corrupt = first;
    ++corrupt.pressureOperator.entries.front().columnControlVolumeIndex;
    expectInvalid(
        [&] { validateSceneFluidPressureEpoch(
            corrupt, fixture.surface.definition, state, grid(),
            fixture.transfer, fixture.connectivity); },
        "pressure epoch rejects nested operator corruption");

    SceneFluidPressureEpochLimits limits;
    limits.maximumEpochBytes = first.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureEpoch(
            fixture.surface.definition, state, grid(), fixture.transfer,
            fixture.connectivity, {}, limits)); },
        "pressure epoch enforces its aggregate owned-storage limit");

    advanceApex(fixture, -0.8);
    const auto movedState = captureSceneFluidSurfaceState(
        fixture.surface.definition,
        fixture.structureAssembly.mappings,
        fixture.structure);
    expectInvalid(
        [&] { validateSceneFluidPressureEpoch(
            first, fixture.surface.definition, movedState, grid(),
            fixture.transfer, fixture.connectivity); },
        "pressure epoch rejects a foreign accepted Structure state");
}

void testCorruptionAndLimits() {
    Fixture fixture;
    const auto previous = captureEndpoint(fixture);
    advanceApex(fixture, -0.8);
    const auto current = captureEndpoint(fixture);
    const auto accepted = buildSceneFluidPressureVolumeRates(
        previous.volumes, current.volumes, current.pressureVolumes);
    auto corrupt = accepted;
    corrupt.controlVolumes.front().volumeChangeCubicMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidPressureVolumeRateIntegrity(corrupt); },
        "pressure-volume-rate integrity rejects nested corruption");

    auto corruptSource = current.pressureVolumes;
    corruptSource.controlVolumes.front().volumeCubicMeters += 0.01;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidPressureVolumeRates(
            previous.volumes, current.volumes, corruptSource)); },
        "pressure-volume rate rejects a corrupted pressure topology");

    SceneFluidPressureVolumeRateLimits limits;
    limits.maximumControlVolumes = accepted.controlVolumes.size() - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureVolumeRates(
            previous.volumes, current.volumes,
            current.pressureVolumes, limits)); },
        "pressure-volume rate bounds control-volume count");
    limits = {};
    limits.maximumVolumeRateBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidPressureVolumeRates(
            previous.volumes, current.volumes,
            current.pressureVolumes, limits)); },
        "pressure-volume rate bounds owned storage");
}

} // namespace

int main(const int argc, const char* const argv[]) {
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--sampling") {
            testMovingVolumeProjection();
            testIndependentGaugeSamplingRejection();
        } else if (argc == 2
                   && std::string_view(argv[1]) == "--epoch") {
            testComposedPressureEpoch();
        } else if (argc == 1) {
            testTopologyStableMovingRates();
            testStationaryAndTopologyRebase();
            testAppearedControlRegionRebase();
            testCorruptionAndLimits();
        } else {
            std::fprintf(stderr, "unexpected test argument\n");
            return 2;
        }
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d pressure-volume-rate check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid pressure-volume-rate/sampling checks passed");
    return 0;
}
