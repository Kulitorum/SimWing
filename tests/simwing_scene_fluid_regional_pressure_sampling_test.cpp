#include "fluid/planar_region_fragment_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_pressure_state.h"
#include "fluid/planar_region_fragment_projection_energy.h"
#include "fluid/planar_region_fragment_pressure_jump_energy.h"
#include "fluid/planar_region_fragment_surface_load.h"
#include "fluid/planar_region_fragment_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"
#include "scene_fluid_regional_pressure_sampling.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::fsi;
using namespace simwing::fsi::fluid;

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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected, message);
}

PeriodicCartesianGrid grid() {
    return {{4, 2, 2}, {-2.0, -1.0, -1.0}, {2.0, 1.0, 1.0}};
}

std::vector<PlanarPressureJumpLayerDefinition> pocketLayers() {
    return {
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -0.8, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         -0.2, -70.0},
    };
}

std::vector<double> metricVelocity(
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::vector<double>& linkVelocity,
    const PlanarPressureRegionFragmentVolumeRateSet* volumeRates) {
    std::vector<double> result(metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            result[dof.dofIndex] = linkVelocity[dof.sourceFaceLinkIndex];
        } else if (volumeRates != nullptr) {
            const auto& rate = volumeRates->fragments[dof.ownerFragmentIndex];
            result[dof.dofIndex] = dof.kind
                    == PlanarPressureRegionFragmentVelocityDofKind::
                        PressureLayerMinusTrace
                ? rate.upperBoundaryVelocityMetersPerSecond
                : rate.lowerBoundaryVelocityMetersPerSecond;
        }
    }
    return result;
}

struct RegionalEndpoint {
    PeriodicCartesianGrid geometry = grid();
    std::vector<PlanarPressureJumpLayerDefinition> previousLayers =
        pocketLayers();
    std::vector<PlanarPressureJumpLayerDefinition> currentLayers;
    PlanarPressureRegionSweepLedger sweep;
    PlanarPressureRegionFragmentSet fragments;
    PlanarPressureRegionFragmentTopology topology;
    PlanarPressureRegionFragmentVelocityMetric metric;
    PlanarPressureRegionFragmentPressureOperator pressureOperator;
    PlanarPressureRegionFragmentVolumeRateSet volumeRates;
    PlanarPressureRegionFragmentVelocityState before;
    PlanarPressureRegionFragmentVelocityState after;
    PlanarPressureRegionFragmentProjectionEnergyAudit projectionEnergy;
    PlanarPressureRegionFragmentPressureJumpEnergyAudit jumpEnergy;
    PlanarPressureRegionFragmentPressureState pressure;
    PlanarPressureRegionFragmentSurfaceLoadLedger surfaceLoads;
    PlanarPressureRegionFragmentAcceptedState accepted;
    bool moving = false;

    explicit RegionalEndpoint(const bool movingGeometry)
        : moving(movingGeometry) {
        currentLayers = moving
            ? translatePlanarPressureJumpLayers(
                  geometry, previousLayers, 0.1).layers
            : previousLayers;
        const double durationSeconds = moving ? 1.0 : 0.01;
        sweep = makePlanarPressureRegionSweepLedger(
            geometry, previousLayers, currentLayers, durationSeconds);
        fragments = buildPlanarPressureRegionFragments(geometry, sweep);
        topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        metric = buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, sweep, fragments, topology);
        pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        if (moving) {
            volumeRates = buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology);
        }

        PlanarPressureRegionFragmentPressureProjectionSettings settings;
        settings.densityKgPerCubicMeter = 1.2;
        settings.timeStepSeconds = durationSeconds;
        settings.absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
        settings.relativeContinuityTolerance = 1.0e-10;
        settings.pressureSolve.absoluteResidualTolerancePascalsMeters =
            1.0e-13;
        settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
        settings.pressureSolve.maximumIterations = 200;

        std::vector<double> linkVelocity(topology.links.size(), 0.0);
        const auto beforeLinkVelocity = linkVelocity;
        std::vector<double> pressureCorrection(
            pressureOperator.rows.size(), 0.0);
        if (moving) {
            const auto diagnostics =
                projectMovingPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    volumeRates, linkVelocity, pressureCorrection, settings);
            if (!diagnostics.accepted) {
                throw std::runtime_error(
                    "moving regional test projection did not converge");
            }
        } else {
            const auto diagnostics =
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    linkVelocity, pressureCorrection, settings);
            if (!diagnostics.accepted) {
                throw std::runtime_error(
                    "static regional test projection did not converge");
            }
        }

        const auto* rates = moving ? &volumeRates : nullptr;
        before = buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            metricVelocity(metric, beforeLinkVelocity, rates),
            settings.densityKgPerCubicMeter);
        after = buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            metricVelocity(metric, linkVelocity, rates),
            settings.densityKgPerCubicMeter);

        PlanarPressureRegionFragmentProjectionEnergySettings energySettings;
        energySettings.densityKgPerCubicMeter =
            settings.densityKgPerCubicMeter;
        energySettings.timeStepSeconds = settings.timeStepSeconds;
        energySettings.absoluteContinuityToleranceCubicMetersPerSecond =
            settings.absoluteContinuityToleranceCubicMetersPerSecond;
        energySettings.relativeContinuityTolerance =
            settings.relativeContinuityTolerance;
        PlanarPressureRegionFragmentPressureJumpEnergySettings jumpSettings;
        jumpSettings.timeStepSeconds = settings.timeStepSeconds;
        if (moving) {
            projectionEnergy =
                auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, before, after, pressureCorrection,
                    energySettings);
            jumpEnergy =
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, after, jumpSettings);
            pressure =
                composeMovingPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, before, after, projectionEnergy, jumpEnergy);
        } else {
            projectionEnergy =
                auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                    geometry, sweep, fragments, topology, metric, before,
                    after, pressureCorrection, energySettings);
            jumpEnergy =
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, after,
                    jumpSettings);
            pressure =
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, projectionEnergy, jumpEnergy);
        }
        surfaceLoads =
            capturePlanarPressureRegionFragmentSurfaceLoads(pressure);
        accepted = capturePlanarPressureRegionFragmentAcceptedState(
            geometry, sweep, fragments, topology, metric, after, pressure,
            surfaceLoads);
    }
};

Scene pressureScene(const bool incomplete = false,
                    const bool reverseFirstSheet = false) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:regional-pressure-sampling";
    scene.metadata.exporterVersion = "regional-pressure-sampling-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "pocket"},
    };
    scene.vertices = {
        {100, {-0.8, -1.0, -1.0}},
        {101, {-0.8, 1.0, -1.0}},
        {102, {-0.8, 1.0, 1.0}},
        {103, {-0.8, -1.0, 1.0}},
        {200, {-0.2, -1.0, -1.0}},
        {201, {-0.2, 1.0, -1.0}},
        {202, {-0.2, 1.0, 1.0}},
        {203, {-0.2, -1.0, 1.0}},
    };
    scene.fabricMaterials = {
        {300, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<StableId, 3> firstTriangle = reverseFirstSheet
        ? std::array<StableId, 3>{100, 102, 101}
        : std::array<StableId, 3>{100, 101, 102};
    const std::array<StableId, 3> secondTriangle = reverseFirstSheet
        ? std::array<StableId, 3>{100, 103, 102}
        : std::array<StableId, 3>{100, 102, 103};
    scene.triangles = {
        {1000, firstTriangle,
         {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}}},
         1, 2, 300, 10, SurfaceRole::Skin},
        {1001, secondTriangle,
         {{{0.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}},
         1, 2, 300, 10, SurfaceRole::Skin},
        {2000, {200, 201, 202},
         {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}}},
         2, 1, 300, 20, SurfaceRole::Skin},
        {2001, {200, 202, 203},
         {{{0.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}},
         2, 1, 300, 20, SurfaceRole::Skin},
    };
    if (incomplete) {
        scene.triangles.pop_back();
    }
    return scene;
}

SceneFluidSurfaceState prepareState(
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& mappings,
    Structure& structure,
    const bool moving) {
    if (moving) {
        StructureStepSettings settings;
        settings.timeStepSeconds = 1.0;
        settings.substeps = 1;
        settings.constraintIterations = 8;
        settings.gravityMetersPerSecondSquared = {0.1, 0.0, 0.0};
        settings.velocityDampingPerSecond = 0.0;
        const auto diagnostics = structure.step(settings);
        if (!diagnostics.finite) {
            throw std::runtime_error(
                "moving scene pressure fixture did not advance");
        }
    }
    return captureSceneFluidSurfaceState(surface, mappings, structure);
}

struct SceneFixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;
    SceneFluidSurfaceState state;
    SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections;
    SceneFluidGridPatchSet patches;
    SceneFluidPatchOwnership ownership;
    SceneFluidQuadratureDefinition quadrature;

    SceneFixture(const bool moving,
                 const bool incomplete = false,
                 const bool reverseFirstSheet = false)
        : scene(pressureScene(incomplete, reverseFirstSheet)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(
              surface.definition, structureAssembly.mappings, structure),
          state(prepareState(
              surface.definition, structureAssembly.mappings, structure,
              moving)),
          candidates(buildSceneFluidGridCandidates(
              surface.definition, state, grid())),
          intersections(intersectSceneFluidSurfaceWithGrid(
              surface.definition, state, grid(), candidates)),
          patches(clipSceneFluidSurfaceToCells(
              surface.definition, state, grid(), candidates,
              intersections)),
          ownership(ownSceneFluidSurfacePatches(
              surface.definition, state, grid(), candidates,
              intersections, patches)),
          quadrature(buildSceneFluidQuadrature(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership, transfer)) {}
};

SceneFluidRegionalPressureSampleSet sample(
    const RegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalPressureSamplingLimits& limits = {}) {
    return sampleSceneFluidRegionalAcceptedPressure(
        endpoint.accepted, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.metric,
        scene.surface.definition, scene.state, scene.quadrature, limits);
}

void testStaticSamplingAndTransfer() {
    const RegionalEndpoint endpoint(false);
    SceneFixture fixture(false);
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "regional sampling: authoritative scene adapters accept both sheets");
    const auto samples = sample(endpoint, fixture);
    const auto repeated = sample(endpoint, fixture);
    check(samples == repeated
              && samples.version
                  == sceneFluidRegionalPressureSamplingVersion
              && samples.fingerprint != 0
              && samples.regionalAcceptedStateFingerprint
                  == endpoint.accepted.fingerprint
              && samples.regionalPressureStateFingerprint
                  == endpoint.pressure.fingerprint
              && samples.regionalSurfaceLoadFingerprint
                  == endpoint.surfaceLoads.fingerprint
              && samples.regionalTopologyFingerprint
                  == endpoint.topology.fingerprint
              && samples.quadratureFingerprint
                  == fixture.quadrature.fingerprint
              && samples.staticGeometry
              && !samples.usesMovingVolumeRates
              && samples.tiles.size() == 8
              && samples.bindings.size() == fixture.quadrature.points.size()
              && samples.bindings.size() > samples.tiles.size()
              && samples.ownedStorageBytes > 0
              && samples.workingStorageBytes > 0,
          "regional sampling: capture is deterministic and fully source-bound");
    checkNear(samples.sampledAreaSquareMeters, 8.0, 2.0e-13,
              "regional sampling: scene quadrature covers both full planes");
    check(samples.maximumAbsoluteTileAreaResidualSquareMeters < 2.0e-14
              && std::max({
                     std::abs(samples.sourceForceResidualNewtons.x),
                     std::abs(samples.sourceForceResidualNewtons.y),
                     std::abs(samples.sourceForceResidualNewtons.z)})
                  < 2.0e-13
              && std::max({
                     std::abs(samples.sourceMomentResidualNewtonMeters.x),
                     std::abs(samples.sourceMomentResidualNewtonMeters.y),
                     std::abs(samples.sourceMomentResidualNewtonMeters.z)})
                  < 2.0e-13
              && samples.sourcePowerResidualWatts == 0.0,
          "regional sampling: tile area, wrench and static power close");
    for (const auto& binding : samples.bindings) {
        checkNear(binding.pressureDifferencePascals,
                  binding.surfaceStableId == 10 ? -70.0 : 70.0,
                  3.0e-14,
                  "regional sampling: each scene patch receives its exact one-sided pressure");
    }
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    validateSceneFluidRegionalAcceptedPressureSamples(
        samples, endpoint.accepted, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.metric,
        fixture.surface.definition, fixture.state, fixture.quadrature);

    const auto transferResult =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            fixture.surface.definition, fixture.state, fixture.transfer,
            fixture.quadrature, samples);
    checkNear(transferResult.diagnostics().surfaceAreaSquareMeters,
              8.0, 2.0e-13,
              "regional sampling: conservative transfer retains sampled area");
    checkNear(transferResult.diagnostics().integratedSurfaceForceNewtons.x,
              0.0, 3.0e-13,
              "regional sampling: opposite plane resultants close globally");
    double firstSheetForce = 0.0;
    double secondSheetForce = 0.0;
    for (const auto& load : transferResult.nodeLoads()) {
        if (load.stableId < 200) {
            firstSheetForce += load.forceNewtons.x;
        } else {
            secondSheetForce += load.forceNewtons.x;
        }
    }
    checkNear(firstSheetForce, -280.0, 4.0e-13,
              "regional sampling: first authored sheet receives its pressure force");
    checkNear(secondSheetForce, 280.0, 4.0e-13,
              "regional sampling: second authored sheet receives its pressure force");
    check(fixture.structure.diagnostics().pendingExternalForceNewtons
              == StructureVector3{},
          "regional sampling: evaluation does not mutate Structure loads");
}

void testMovingSamplingAndPower() {
    const RegionalEndpoint endpoint(true);
    SceneFixture fixture(true);
    const auto samples = sample(endpoint, fixture);
    check(!samples.staticGeometry && samples.usesMovingVolumeRates
              && samples.acceptedStepCount == 1
              && samples.sampledPressurePowerToSheetWatts != 0.0,
          "regional sampling: moving scene retains epoch and nonzero sheet power");
    checkNear(samples.sampledPressurePowerToSheetWatts,
              endpoint.accepted.pressureWorkToSheetJoules
                  / endpoint.accepted.timeStepSeconds,
              2.0e-12,
              "regional sampling: moving quadrature power closes accepted sheet work");
    checkNear(samples.sourcePowerResidualWatts, 0.0, 2.0e-12,
              "regional sampling: moving source-power residual closes");
    const auto transferResult =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            fixture.surface.definition, fixture.state, fixture.transfer,
            fixture.quadrature, samples);
    checkNear(transferResult.diagnostics().integratedSurfacePowerWatts,
              samples.sampledPressurePowerToSheetWatts, 2.0e-12,
              "regional sampling: conservative Structure transfer preserves moving power");
}

void testRejectionAndLimits() {
    const RegionalEndpoint endpoint(false);
    SceneFixture fixture(false);
    const auto samples = sample(endpoint, fixture);
    auto corrupt = samples;
    corrupt.pressures[0].negativeSidePressurePascals += 1.0;
    expectRejected(
        [&] { validateSceneFluidRegionalPressureSampleIntegrity(corrupt); },
        "regional sampling: nested pressure corruption rejects");
    corrupt = samples;
    corrupt.tiles[0].sampledAreaSquareMeters += 0.1;
    expectRejected(
        [&] { validateSceneFluidRegionalPressureSampleIntegrity(corrupt); },
        "regional sampling: nested tile-coverage corruption rejects");

    auto limits = SceneFluidRegionalPressureSamplingLimits{};
    limits.maximumSamples = samples.bindings.size() - 1;
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, fixture, limits)); },
        "regional sampling: sample-count limit rejects transactionally");
    limits = {};
    limits.maximumTiles = samples.tiles.size() - 1;
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, fixture, limits)); },
        "regional sampling: tile-count limit rejects transactionally");
    limits = {};
    limits.maximumOwnedBytes = samples.ownedStorageBytes - 1;
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, fixture, limits)); },
        "regional sampling: owned-byte limit rejects transactionally");
    limits = {};
    limits.maximumWorkingBytes = samples.workingStorageBytes - 1;
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, fixture, limits)); },
        "regional sampling: working-byte limit rejects transactionally");

    SceneFixture incomplete(false, true);
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, incomplete)); },
        "regional sampling: incomplete authored-sheet coverage rejects");
    expectRejected(
        [&] {
            SceneFixture reversed(false, false, true);
            static_cast<void>(sample(endpoint, reversed));
        },
        "regional sampling: reversed scene normal rejects");

    auto foreignAccepted = endpoint.accepted;
    foreignAccepted.surfaceLoads.tiles[0].surfaceStableId = 999;
    expectRejected(
        [&] {
            static_cast<void>(sampleSceneFluidRegionalAcceptedPressure(
                foreignAccepted, endpoint.geometry, endpoint.sweep,
                endpoint.fragments, endpoint.topology, endpoint.metric,
                fixture.surface.definition, fixture.state,
                fixture.quadrature));
        },
        "regional sampling: foreign accepted endpoint rejects");
}

} // namespace

int main() {
    try {
        testStaticSamplingAndTransfer();
        testMovingSamplingAndPower();
        testRejectionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d regional pressure-sampling test(s) failed\n",
                     failures);
        return 1;
    }
    std::printf("regional pressure-sampling tests passed\n");
    return 0;
}
