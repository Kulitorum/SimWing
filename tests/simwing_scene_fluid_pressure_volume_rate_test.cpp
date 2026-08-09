#include "scene_fluid_pressure_epoch.h"
#include "scene_fluid_pressure_sampling.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
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
    auto volumes = buildSceneFluidCellVolumes(
        fixture.surface.definition, state, grid(), fixture.transfer, epoch);
    auto pressureVolumes = buildSceneFluidPressureControlVolumes(
        fixture.surface.definition, volumes, fixture.connectivity);
    return {std::move(state), std::move(epoch), std::move(caps),
            std::move(openingQuadrature), std::move(openingPatches),
            std::move(volumes), std::move(pressureVolumes)};
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
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidPressureVolumeRates(
            crossingPrevious.volumes, crossingCurrent.volumes,
            crossingCurrent.pressureVolumes)); },
        "pressure-volume rate rejects a cell-crossing topology rebase");
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
        current.openingPatches, current.volumes, fixture.connectivity,
        current.pressureVolumes);
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

void testIndependentGaugeSamplingRejection() {
    Fixture fixture(false);
    const auto current = captureEndpoint(fixture);
    const auto faceLinks = buildSceneFluidPressureFaceLinks(
        fixture.surface.definition, current.state, grid(), fixture.transfer,
        current.epoch, current.caps, current.openingQuadrature,
        current.openingPatches, current.volumes, fixture.connectivity,
        current.pressureVolumes);
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
              && first.cellVolumes.cellRegionVolumes.size()
                  == first.pressureControlVolumes.controlVolumes.size()
              && first.pressureFaceLinks.unresolvedActiveFaceCount == 0
              && first.pressureFaceLinks.unresolvedAmbiguousFaceCount == 0
              && first.pressureFaceLinks.unresolvedOpeningFaceCount == 0
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
