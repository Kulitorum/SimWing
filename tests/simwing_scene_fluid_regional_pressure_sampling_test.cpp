#include "fluid/planar_region_fragment_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_pressure_state.h"
#include "fluid/planar_region_fragment_projection_energy.h"
#include "fluid/planar_region_fragment_pressure_jump_energy.h"
#include "fluid/planar_region_fragment_surface_load.h"
#include "fluid/planar_region_fragment_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"
#include "scene_fluid_regional_pressure_sampling.h"
#include "scene_fluid_regional_opening_load_epoch.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>
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

bool samePublicCheckpoint(const StructureCheckpoint& first,
                          const StructureCheckpoint& second) {
    return first.version == second.version
        && first.definitionFingerprint == second.definitionFingerprint
        && first.acceptedStepCount == second.acceptedStepCount
        && first.simulationTimeSeconds == second.simulationTimeSeconds
        && first.nodes == second.nodes
        && first.pendingExternalForcesNewtons
            == second.pendingExternalForcesNewtons
        && first.lastAppliedExternalForceNewtons
            == second.lastAppliedExternalForceNewtons;
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

struct OpeningRegionalEndpoint {
    PeriodicCartesianGrid geometry = grid();
    std::vector<PlanarPressureJumpLayerDefinition> layers = pocketLayers();
    PlanarPressureRegionSweepLedger sweep;
    PlanarPressureRegionFragmentSet fragments;
    PlanarPressureRegionFragmentTopology topology;
    PlanarPressureRegionFragmentPressureOperator basePressureOperator;
    PlanarPressureRegionFragmentVolumeRateSet volumeRates;
    std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions;
    PlanarPressureRegionFragmentOpeningSet openings;
    PlanarPressureRegionFragmentOpeningPressureOperator pressureOperator;
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedFlow;
    PlanarPressureRegionFragmentOpeningPressureState pressure;
    PlanarPressureRegionFragmentSurfaceLoadLedger surfaceLoads;
    PlanarPressureRegionFragmentOpeningSurfaceLoadLedger openingSurfaceLoads;
    PlanarPressureRegionFragmentOpeningLoadState loadState;

    explicit OpeningRegionalEndpoint(const bool fullyOpenFirstSurface = true) {
        constexpr double durationSeconds = 0.01;
        sweep = makePlanarPressureRegionSweepLedger(
            geometry, layers, layers, durationSeconds);
        fragments = buildPlanarPressureRegionFragments(geometry, sweep);
        topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        basePressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        volumeRates = buildPlanarPressureRegionFragmentVolumeRates(
            geometry, sweep, fragments, topology);
        for (const auto& link : topology.links) {
            if (link.kind
                    != PlanarPressureRegionFragmentFaceKind::PressureLayerWall
                || link.surfaceStableId != 10) {
                continue;
            }
            if (!fullyOpenFirstSurface
                && (link.j != 0 || link.k != 0)) {
                continue;
            }
            const std::uint64_t patchId =
                100 + openingDefinitions.size();
            PlanarPressureRegionFragmentOpeningPatchDefinition definition{
                patchId,
                1000,
                link.surfaceStableId,
                link.axis,
                link.i,
                link.j,
                link.k,
                link.minusRegionStableId,
                link.plusRegionStableId,
                fullyOpenFirstSurface
                    ? link.areaSquareMeters
                    : 0.5 * link.areaSquareMeters,
            };
            if (!fullyOpenFirstSurface) {
                definition.authoredWrappedCentroidMeters =
                    Vector3{-0.8, -0.75, -0.5};
            }
            openingDefinitions.push_back(definition);
            resistanceDefinitions.push_back({patchId, {0.0, 0.0}});
        }
        const std::size_t expectedDefinitionCount =
            fullyOpenFirstSurface ? 4 : 1;
        if (openingDefinitions.size() != expectedDefinitionCount) {
            throw std::runtime_error(
                "opening regional fixture found the wrong wall-tile count");
        }
        openings = buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, openingDefinitions);
        pressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                basePressureOperator, geometry, sweep, fragments, topology,
                openingDefinitions, openings);

        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            samples;
        samples.reserve(openingDefinitions.size());
        for (const auto& definition : openingDefinitions)
            samples.push_back({definition.patchStableId, 0.0});
        auto flux = buildPlanarPressureRegionFragmentOpeningFluxState(
            geometry, sweep, fragments, topology, openingDefinitions,
            openings, samples);
        std::vector<double> velocity(topology.links.size(), 0.0);
        std::vector<double> pressureCorrection(
            basePressureOperator.rows.size(), 0.0);
        PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
        settings.projection.densityKgPerCubicMeter = 1.2;
        settings.projection.timeStepSeconds = durationSeconds;
        settings.projection
            .absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
        settings.projection.relativeContinuityTolerance = 1.0e-10;
        settings.projection
            .absoluteMomentumResidualToleranceKilogramMetersPerSecond =
            1.0e-12;
        settings.projection.relativeMomentumResidualTolerance = 1.0e-10;
        settings.projection.absoluteEnergyResidualToleranceJoules = 1.0e-12;
        settings.projection.relativeEnergyResidualTolerance = 1.0e-10;
        settings.projection.pressureSolve
            .absoluteResidualTolerancePascalsMeters = 1.0e-13;
        settings.projection.pressureSolve.relativeResidualTolerance = 0.0;
        settings.projection.pressureSolve.maximumIterations = 200;
        const auto diagnostics =
            advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
                pressureOperator, basePressureOperator, geometry, sweep,
                fragments, topology, volumeRates, openingDefinitions,
                openings, resistanceDefinitions, velocity, samples, flux,
                pressureCorrection, settings);
        if (!diagnostics.accepted) {
            throw std::runtime_error(
                "opening regional fixture pressure step did not converge");
        }
        acceptedFlow =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                pressureOperator, basePressureOperator, geometry, sweep,
                fragments, topology, volumeRates, openingDefinitions,
                openings, resistanceDefinitions, diagnostics, velocity,
                samples, flux, pressureCorrection, settings);
        pressure = composePlanarPressureRegionFragmentOpeningPressureState(
            acceptedFlow, pressureOperator, basePressureOperator, geometry,
            sweep, fragments, topology, volumeRates, openingDefinitions,
            openings, resistanceDefinitions);
        surfaceLoads = capturePlanarPressureRegionFragmentSurfaceLoads(
            pressure);
        openingSurfaceLoads =
            capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                surfaceLoads, pressure, geometry, sweep, fragments, topology,
                openingDefinitions, openings);
        loadState = capturePlanarPressureRegionFragmentOpeningLoadState(
            acceptedFlow, pressure, surfaceLoads, openingSurfaceLoads,
            pressureOperator, basePressureOperator, geometry, sweep,
            fragments, topology, volumeRates, openingDefinitions, openings,
            resistanceDefinitions);
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

Scene openingPressureScene() {
    auto scene = pressureScene();
    std::erase_if(
        scene.triangles,
        [](const Triangle& triangle) { return triangle.sheetId == 10; });
    std::erase_if(
        scene.vertices,
        [](const Vertex& vertex) { return vertex.id < 200; });
    scene.metadata.designChecksum =
        "sha256:regional-opening-pressure-sampling";
    return scene;
}

Scene partialOpeningPressureScene() {
    auto scene = pressureScene();
    std::erase_if(
        scene.triangles,
        [](const Triangle& triangle) { return triangle.sheetId == 10; });
    std::erase_if(
        scene.vertices,
        [](const Vertex& vertex) { return vertex.id < 200; });
    const std::vector<Vertex> retainedVertices{
        {100, {-0.8, -0.5, -1.0}},
        {101, {-0.8, 0.0, -1.0}},
        {102, {-0.8, 0.0, 0.0}},
        {103, {-0.8, -0.5, 0.0}},
        {104, {-0.8, -1.0, 0.0}},
        {105, {-0.8, 0.0, 1.0}},
        {106, {-0.8, -1.0, 1.0}},
        {107, {-0.8, 1.0, -1.0}},
        {108, {-0.8, 1.0, 0.0}},
        {109, {-0.8, 1.0, 1.0}},
    };
    scene.vertices.insert(
        scene.vertices.begin(),
        retainedVertices.begin(), retainedVertices.end());
    const auto materialPoint = [](const double y, const double z) {
        return Vec2{y + 1.0, z + 1.0};
    };
    const auto addRectangle = [&](const StableId firstTriangleId,
                                  const std::array<StableId, 4>& vertices,
                                  const double lowerY,
                                  const double upperY,
                                  const double lowerZ,
                                  const double upperZ) {
        scene.triangles.push_back({
            firstTriangleId,
            {vertices[0], vertices[1], vertices[2]},
            {{materialPoint(lowerY, lowerZ),
              materialPoint(upperY, lowerZ),
              materialPoint(upperY, upperZ)}},
            1, 2, 300, 10, SurfaceRole::Skin,
        });
        scene.triangles.push_back({
            firstTriangleId + 1,
            {vertices[0], vertices[2], vertices[3]},
            {{materialPoint(lowerY, lowerZ),
              materialPoint(upperY, upperZ),
              materialPoint(lowerY, upperZ)}},
            1, 2, 300, 10, SurfaceRole::Skin,
        });
    };
    addRectangle(1000, {100, 101, 102, 103}, -0.5, 0.0, -1.0, 0.0);
    addRectangle(1010, {104, 102, 105, 106}, -1.0, 0.0, 0.0, 1.0);
    addRectangle(1020, {101, 107, 108, 102}, 0.0, 1.0, -1.0, 0.0);
    addRectangle(1030, {102, 108, 109, 105}, 0.0, 1.0, 0.0, 1.0);
    std::ranges::sort(
        scene.triangles, {}, &Triangle::id);
    scene.metadata.designChecksum =
        "sha256:regional-partial-opening-pressure-sampling";
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

    SceneFixture(Scene inputScene, const bool moving)
        : scene(std::move(inputScene)),
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

    SceneFixture(const bool moving,
                 const bool incomplete = false,
                 const bool reverseFirstSheet = false)
        : SceneFixture(
              pressureScene(incomplete, reverseFirstSheet), moving) {}
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

SceneFluidRegionalPressureSampleSet sample(
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalPressureSamplingLimits& limits = {}) {
    return sampleSceneFluidRegionalOpeningPressure(
        endpoint.loadState, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, scene.surface.definition,
        scene.state, scene.quadrature, {}, limits);
}

SceneFluidRegionalOpeningLoadEpoch applyOpeningLoadEpoch(
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningLoadEpochLimits& limits = {}) {
    return applySceneFluidRegionalOpeningLoadEpoch(
        endpoint.acceptedFlow, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, scene.surface.definition,
        scene.state, scene.transfer, scene.quadrature, scene.structure,
        settings, limits);
}

void validateOpeningLoadEpoch(
    const SceneFluidRegionalOpeningLoadEpoch& result,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningLoadEpochLimits& limits = {}) {
    validateSceneFluidRegionalOpeningLoadEpoch(
        result, endpoint.acceptedFlow, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, scene.surface.definition,
        scene.state, scene.transfer, scene.quadrature, settings, limits);
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

void testOpeningAwareSamplingAndApplication() {
    const OpeningRegionalEndpoint endpoint;
    SceneFixture fixture(openingPressureScene(), false);
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "regional opening sampling: retained sheet assembles");
    const auto samples = sample(endpoint, fixture);
    const auto repeated = sample(endpoint, fixture);
    check(samples == repeated && samples.fingerprint != 0
              && samples.openingAware
              && samples.regionalAcceptedStateFingerprint == 0
              && samples.regionalOpeningLoadStateFingerprint
                  == endpoint.loadState.fingerprint
              && samples.regionalPressureStateFingerprint
                  == endpoint.pressure.fingerprint
              && samples.regionalSurfaceLoadFingerprint
                  == endpoint.surfaceLoads.fingerprint
              && samples.staticGeometry
              && samples.usesMovingVolumeRates
              && samples.tiles.size() == 8
              && samples.bindings.size() == fixture.quadrature.points.size(),
          "regional opening sampling: capture is deterministic and bound to the atomic endpoint");
    std::size_t zeroSolidTiles = 0;
    for (const auto& tile : samples.tiles) {
        if (tile.surfaceStableId == 10) {
            check(tile.sourceAreaSquareMeters == 0.0
                      && tile.sampledAreaSquareMeters == 0.0
                      && tile.sampleCount == 0
                      && tile.areaResidualSquareMeters == 0.0,
                  "regional opening sampling: fully open tile receives no fabric sample");
            ++zeroSolidTiles;
        } else {
            check(tile.surfaceStableId == 20
                      && tile.sourceAreaSquareMeters > 0.0
                      && tile.sampleCount > 0,
                  "regional opening sampling: retained tile remains covered");
        }
    }
    check(zeroSolidTiles == 4,
          "regional opening sampling: all removed wall tiles remain explicit");
    checkNear(samples.sampledAreaSquareMeters,
              endpoint.loadState.solidAreaSquareMeters, 2.0e-13,
              "regional opening sampling: quadrature covers retained solid area only");
    checkNear(samples.sampledPressureForceOnSheetNewtons.x,
              endpoint.loadState.solidPressureForceOnSheetNewtons.x,
              4.0e-13,
              "regional opening sampling: retained pressure force closes");
    check(std::max({
              std::abs(samples.sourceForceResidualNewtons.x),
              std::abs(samples.sourceForceResidualNewtons.y),
              std::abs(samples.sourceForceResidualNewtons.z),
              std::abs(samples.sourceMomentResidualNewtonMeters.x),
              std::abs(samples.sourceMomentResidualNewtonMeters.y),
              std::abs(samples.sourceMomentResidualNewtonMeters.z),
              std::abs(samples.sourcePowerResidualWatts)}) < 4.0e-13,
          "regional opening sampling: retained wrench and power close");
    for (const auto& binding : samples.bindings) {
        check(binding.surfaceStableId == 20,
              "regional opening sampling: no aperture patch receives traction");
    }
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    validateSceneFluidRegionalOpeningPressureSamples(
        samples, endpoint.loadState, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, fixture.surface.definition,
        fixture.state, fixture.quadrature);

    const auto transferResult =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            fixture.surface.definition, fixture.state, fixture.transfer,
            fixture.quadrature, samples);
    checkNear(transferResult.diagnostics().surfaceAreaSquareMeters,
              endpoint.loadState.solidAreaSquareMeters, 2.0e-13,
              "regional opening sampling: read-only transfer retains solid area");
    checkNear(transferResult.diagnostics().integratedSurfaceForceNewtons.x,
              endpoint.loadState.solidPressureForceOnSheetNewtons.x,
              4.0e-13,
              "regional opening sampling: read-only transfer retains solid force");
    const auto application = applySceneFluidRegionalAcceptedPressureLoads(
        fixture.surface.definition, fixture.state, fixture.transfer,
        fixture.quadrature, samples, fixture.structure);
    check(application.applied
              && application.sourceSamplingFingerprint
                  == samples.fingerprint
              && application.appliedPressureForceNewtons
                  == samples.sampledPressureForceOnSheetNewtons
              && application.applicationResidualNewtons
                  == StructureVector3{},
          "regional opening sampling: transactional application retains the solid load");
    check(fixture.structure.diagnostics().pendingExternalForceNewtons
              == samples.sampledPressureForceOnSheetNewtons,
          "regional opening sampling: Structure receives no aperture traction");
    validateSceneFluidRegionalPressureLoadApplicationIntegrity(application);
    validateSceneFluidRegionalPressureLoadApplication(
        application, fixture.surface.definition, fixture.state,
        fixture.transfer, fixture.quadrature, samples);

    auto corrupt = samples;
    corrupt.regionalOpeningLoadStateFingerprint ^= 1U;
    expectRejected(
        [&] { validateSceneFluidRegionalPressureSampleIntegrity(corrupt); },
        "regional opening sampling: endpoint fingerprint corruption rejects");
    SceneFluidRegionalPressureSamplingLimits limits;
    limits.maximumSamples = samples.bindings.size() - 1;
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, fixture, limits)); },
        "regional opening sampling: sample limit rejects");
    SceneFixture foreignFullSheet(false);
    expectRejected(
        [&] { static_cast<void>(sample(endpoint, foreignFullSheet)); },
        "regional opening sampling: traction over removed aperture area rejects");
}

void testPartialOpeningSamplingAndApplication() {
    const OpeningRegionalEndpoint endpoint(false);
    SceneFixture fixture(partialOpeningPressureScene(), false);
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "regional partial opening: retained fabric assembles");
    const auto touchedLoad = std::ranges::find(
        endpoint.openingSurfaceLoads.tiles, true,
        &PlanarPressureRegionFragmentOpeningSurfaceLoadTile::
            touchedByOpening);
    check(touchedLoad != endpoint.openingSurfaceLoads.tiles.end()
              && touchedLoad->hasExactSubtileCentroids
              && touchedLoad->openingAreaSquareMeters == 0.5
              && touchedLoad->solidAreaSquareMeters == 0.5,
          "regional partial opening: atomic endpoint retains exact half-tile geometry");
    if (touchedLoad == endpoint.openingSurfaceLoads.tiles.end()) return;

    const auto samples = sample(endpoint, fixture);
    const auto repeated = sample(endpoint, fixture);
    check(samples == repeated && samples.openingAware
              && samples.regionalOpeningLoadStateFingerprint
                  == endpoint.loadState.fingerprint
              && samples.sampledAreaSquareMeters == 7.5,
          "regional partial opening: sampling is deterministic and source-bound");
    const auto touchedCoverage = std::ranges::find(
        samples.tiles, touchedLoad->sourceFaceLinkStableId,
        &SceneFluidRegionalPressureTileCoverage::
            sourceFaceLinkStableId);
    check(touchedCoverage != samples.tiles.end()
              && touchedCoverage->sourceAreaSquareMeters == 0.5
              && touchedCoverage->sampledAreaSquareMeters == 0.5
              && touchedCoverage->sampleCount > 0
              && touchedCoverage->areaResidualSquareMeters == 0.0,
          "regional partial opening: quadrature covers the retained half tile exactly");
    checkNear(
        samples.sampledPressureForceOnSheetNewtons.x,
        endpoint.loadState.solidPressureForceOnSheetNewtons.x,
        5.0e-13,
        "regional partial opening: retained force closes to the atomic endpoint");
    check(std::max({
              std::abs(samples.sourceForceResidualNewtons.x),
              std::abs(samples.sourceForceResidualNewtons.y),
              std::abs(samples.sourceForceResidualNewtons.z),
              std::abs(samples.sourceMomentResidualNewtonMeters.x),
              std::abs(samples.sourceMomentResidualNewtonMeters.y),
              std::abs(samples.sourceMomentResidualNewtonMeters.z),
              std::abs(samples.sourcePowerResidualWatts)}) < 5.0e-13,
          "regional partial opening: exact retained wrench and power close");
    validateSceneFluidRegionalOpeningPressureSamples(
        samples, endpoint.loadState, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, fixture.surface.definition,
        fixture.state, fixture.quadrature);

    const auto transferResult =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            fixture.surface.definition, fixture.state, fixture.transfer,
            fixture.quadrature, samples);
    checkNear(
        transferResult.diagnostics().surfaceAreaSquareMeters,
        endpoint.loadState.solidAreaSquareMeters, 3.0e-13,
        "regional partial opening: conservative transfer retains exact solid area");
    checkNear(
        transferResult.diagnostics().integratedSurfaceMomentNewtonMeters.z,
        endpoint.loadState.solidPressureMomentOnSheetNewtonMeters.z,
        5.0e-13,
        "regional partial opening: conservative transfer retains exact solid moment");
    check(fixture.structure.diagnostics().pendingExternalForceNewtons
              == StructureVector3{},
          "regional partial opening: read-only transfer does not mutate Structure");
    const auto beforeApplication = fixture.structure.checkpoint();
    SceneFluidRegionalPressureLoadApplicationLimits applicationLimits;
    applicationLimits.maximumNodeLoads = fixture.transfer.nodes().size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(applySceneFluidRegionalAcceptedPressureLoads(
                fixture.surface.definition, fixture.state, fixture.transfer,
                fixture.quadrature, samples, fixture.structure, {},
                applicationLimits));
        },
        "regional partial opening: application limit rejects before mutation");
    check(samePublicCheckpoint(
              fixture.structure.checkpoint(), beforeApplication),
          "regional partial opening: rejected application preserves the full checkpoint");
    const auto application = applySceneFluidRegionalAcceptedPressureLoads(
        fixture.surface.definition, fixture.state, fixture.transfer,
        fixture.quadrature, samples, fixture.structure);
    check(application.applied
              && application.applicationResidualNewtons
                  == StructureVector3{},
          "regional partial opening: transactional application closes the retained load");
    for (const auto& component : {
             std::pair{application.appliedPressureForceNewtons.x,
                       samples.sampledPressureForceOnSheetNewtons.x},
             std::pair{application.appliedPressureForceNewtons.y,
                       samples.sampledPressureForceOnSheetNewtons.y},
             std::pair{application.appliedPressureForceNewtons.z,
                       samples.sampledPressureForceOnSheetNewtons.z},
             std::pair{application.resultingPendingForceNewtons.x,
                       samples.sampledPressureForceOnSheetNewtons.x},
             std::pair{application.resultingPendingForceNewtons.y,
                       samples.sampledPressureForceOnSheetNewtons.y},
             std::pair{application.resultingPendingForceNewtons.z,
                       samples.sampledPressureForceOnSheetNewtons.z},
             std::pair{
                 fixture.structure.diagnostics()
                     .pendingExternalForceNewtons.x,
                 samples.sampledPressureForceOnSheetNewtons.x},
             std::pair{
                 fixture.structure.diagnostics()
                     .pendingExternalForceNewtons.y,
                 samples.sampledPressureForceOnSheetNewtons.y},
             std::pair{
                 fixture.structure.diagnostics()
                     .pendingExternalForceNewtons.z,
                 samples.sampledPressureForceOnSheetNewtons.z}}) {
        checkNear(
            component.first, component.second, 8.0e-13,
            "regional partial opening: Structure receives the retained resultant");
    }
    validateSceneFluidRegionalPressureLoadApplicationIntegrity(application);
    validateSceneFluidRegionalPressureLoadApplication(
        application, fixture.surface.definition, fixture.state,
        fixture.transfer, fixture.quadrature, samples);
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

void testTransactionalLoadApplication() {
    const RegionalEndpoint endpoint(false);
    SceneFixture fixture(false);
    const auto samples = sample(endpoint, fixture);
    std::vector<StructureVector3> priorLoads(
        fixture.structure.definition().nodes.size());
    for (std::size_t index = 0; index < priorLoads.size(); ++index) {
        const double value = 0.01 * static_cast<double>(index + 1);
        priorLoads[index] = {value, -2.0 * value, 0.5 * value};
    }
    fixture.structure.setExternalForces(priorLoads);
    const auto before = fixture.structure.checkpoint();
    const auto application = applySceneFluidRegionalAcceptedPressureLoads(
        fixture.surface.definition, fixture.state, fixture.transfer,
        fixture.quadrature, samples, fixture.structure);
    const auto after = fixture.structure.checkpoint();
    check(application.version
                  == sceneFluidRegionalPressureLoadApplicationVersion
              && application.fingerprint != 0
              && application.applied
              && application.sourceSamplingFingerprint
                  == samples.fingerprint
              && application.sourceSurfaceStateFingerprint
                  == fixture.state.fingerprint
              && application.couplingSurfaceFingerprint
                  == fixture.transfer.couplingSurfaceFingerprint()
              && application.targetDefinitionFingerprint
                  == fixture.structure.definitionFingerprint()
              && application.acceptedStepCount == 0
              && application.simulationTimeSeconds == 0.0
              && application.structureNodeCount == priorLoads.size()
              && application.nodeLoads.size()
                  == fixture.transfer.nodes().size()
              && application.ownedStorageBytes > 0
              && application.workingStorageBytes > 0,
          "regional application: immutable receipt binds sampling and Structure epoch");
    check(after.acceptedStepCount == before.acceptedStepCount
              && after.simulationTimeSeconds == before.simulationTimeSeconds
              && after.nodes == before.nodes
              && after.lastAppliedExternalForceNewtons
                  == before.lastAppliedExternalForceNewtons,
          "regional application: pending-load mutation changes no committed state");
    std::vector<StructureVector3> expected =
        before.pendingExternalForcesNewtons;
    for (const auto& load : application.nodeLoads) {
        expected[load.structureNode].x +=
            load.appliedPressureForceNewtons.x;
        expected[load.structureNode].y +=
            load.appliedPressureForceNewtons.y;
        expected[load.structureNode].z +=
            load.appliedPressureForceNewtons.z;
        check(load.priorPendingForceNewtons
                      == before.pendingExternalForcesNewtons[
                          load.structureNode]
                  && load.resultingPendingForceNewtons
                      == expected[load.structureNode]
                  && load.applicationResidualNewtons
                      == StructureVector3{},
              "regional application: each pressure load preserves its prior node load");
    }
    check(after.pendingExternalForcesNewtons == expected,
          "regional application: all resulting pending loads match the receipt");
    validateSceneFluidRegionalPressureLoadApplicationIntegrity(application);
    validateSceneFluidRegionalPressureLoadApplication(
        application, fixture.surface.definition, fixture.state,
        fixture.transfer, fixture.quadrature, samples);
    auto validationLimits =
        SceneFluidRegionalPressureLoadApplicationLimits{};
    validationLimits.maximumNodeLoads = application.nodeLoads.size() - 1;
    expectRejected(
        [&] {
            validateSceneFluidRegionalPressureLoadApplication(
                application, fixture.surface.definition, fixture.state,
                fixture.transfer, fixture.quadrature, samples, {},
                validationLimits);
        },
        "regional application: source-aware validation enforces its node limit");

    const RegionalEndpoint movingEndpoint(true);
    SceneFixture movingFixture(true);
    const auto movingSamples = sample(movingEndpoint, movingFixture);
    expectRejected(
        [&] {
            validateSceneFluidRegionalPressureLoadApplication(
                application, movingFixture.surface.definition,
                movingFixture.state, movingFixture.transfer,
                movingFixture.quadrature, movingSamples);
        },
        "regional application: source-aware validation rejects a foreign accepted source");
    auto corrupt = application;
    corrupt.nodeLoads[0].resultingPendingForceNewtons.x += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalPressureLoadApplicationIntegrity(
                corrupt);
        },
        "regional application: nested receipt corruption rejects");
}

void testMovingLoadApplication() {
    const RegionalEndpoint endpoint(true);
    SceneFixture fixture(true);
    const auto samples = sample(endpoint, fixture);
    const auto before = fixture.structure.checkpoint();
    const auto application = applySceneFluidRegionalAcceptedPressureLoads(
        fixture.surface.definition, fixture.state, fixture.transfer,
        fixture.quadrature, samples, fixture.structure);
    const auto after = fixture.structure.checkpoint();
    check(application.applied && application.acceptedStepCount == 1
              && application.simulationTimeSeconds == 1.0
              && after.nodes == before.nodes
              && after.acceptedStepCount == before.acceptedStepCount
              && after.simulationTimeSeconds == before.simulationTimeSeconds
              && after.pendingExternalForcesNewtons
                  != before.pendingExternalForcesNewtons,
          "regional application: moving accepted pressure reaches pending XPBD loads without stepping");
    validateSceneFluidRegionalPressureLoadApplication(
        application, fixture.surface.definition, fixture.state,
        fixture.transfer, fixture.quadrature, samples);
}

void testApplicationRollbackAndLimits() {
    const RegionalEndpoint endpoint(false);
    SceneFixture fixture(false);
    const auto samples = sample(endpoint, fixture);
    fixture.structure.addExternalForce(0, {1.0, 2.0, 3.0});
    const auto before = fixture.structure.checkpoint();

    auto limits = SceneFluidRegionalPressureLoadApplicationLimits{};
    limits.maximumNodeLoads = fixture.transfer.nodes().size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, samples,
                    fixture.structure, {}, limits));
        },
        "regional application: node-load limit rejects before mutation");
    check(samePublicCheckpoint(before, fixture.structure.checkpoint()),
          "regional application: node-load limit preserves exact target state");
    limits = {};
    limits.maximumStructureNodes =
        fixture.structure.definition().nodes.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, samples,
                    fixture.structure, {}, limits));
        },
        "regional application: Structure-node limit rejects before mutation");
    check(samePublicCheckpoint(before, fixture.structure.checkpoint()),
          "regional application: Structure-node limit preserves exact target state");

    limits = {};
    limits.maximumOwnedBytes =
        fixture.transfer.nodes().size()
            * sizeof(SceneFluidRegionalPressureAppliedNodeLoad) - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, samples,
                    fixture.structure, {}, limits));
        },
        "regional application: owned-byte limit rejects before mutation");
    check(samePublicCheckpoint(before, fixture.structure.checkpoint()),
          "regional application: owned-byte limit preserves exact target state");

    limits = {};
    limits.maximumWorkingBytes =
        fixture.structure.definition().nodes.size()
            * sizeof(StructureVector3) - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, samples,
                    fixture.structure, {}, limits));
        },
        "regional application: working-byte limit rejects before mutation");
    check(samePublicCheckpoint(before, fixture.structure.checkpoint()),
          "regional application: working-byte limit preserves exact target state");

    auto corruptSamples = samples;
    corruptSamples.fingerprint = 0;
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, corruptSamples,
                    fixture.structure));
        },
        "regional application: corrupt samples reject before mutation");
    check(samePublicCheckpoint(before, fixture.structure.checkpoint()),
          "regional application: corrupt samples preserve exact target state");

    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = 0.01;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    static_cast<void>(fixture.structure.step(stepSettings));
    const auto staleBefore = fixture.structure.checkpoint();
    expectRejected(
        [&] {
            static_cast<void>(
                applySceneFluidRegionalAcceptedPressureLoads(
                    fixture.surface.definition, fixture.state,
                    fixture.transfer, fixture.quadrature, samples,
                    fixture.structure));
        },
        "regional application: stale Structure epoch rejects before mutation");
    check(samePublicCheckpoint(
              staleBefore, fixture.structure.checkpoint()),
          "regional application: stale epoch preserves exact current target state");
}

void testAtomicOpeningLoadEpoch() {
    const OpeningRegionalEndpoint endpoint(false);
    SceneFixture fixture(partialOpeningPressureScene(), false);
    SceneFixture repeatedFixture(partialOpeningPressureScene(), false);
    const auto before = fixture.structure.checkpoint();
    const auto result = applyOpeningLoadEpoch(endpoint, fixture);
    const auto repeated = applyOpeningLoadEpoch(
        endpoint, repeatedFixture);
    check(result == repeated
              && result.version
                  == sceneFluidRegionalOpeningLoadEpochVersion
              && result.fingerprint != 0 && result.applied
              && result.loadState == endpoint.loadState
              && result.samples.openingAware
              && result.samples.regionalOpeningLoadStateFingerprint
                  == result.loadState.fingerprint
              && result.application.sourceSamplingFingerprint
                  == result.samples.fingerprint
              && result.application.applied
              && result.ownedStorageBytes > 0
              && result.workingStorageBytes > 0,
          "regional opening epoch: complete load transaction is deterministic and source-bound");
    check(fixture.structure.checkpoint().pendingExternalForcesNewtons
              != before.pendingExternalForcesNewtons
              && fixture.structure.acceptedStepCount()
                  == before.acceptedStepCount
              && fixture.structure.simulationTimeSeconds()
                  == before.simulationTimeSeconds,
          "regional opening epoch: only pending Structure loads change");
    validateSceneFluidRegionalOpeningLoadEpochIntegrity(result);
    validateOpeningLoadEpoch(result, endpoint, fixture);

    auto corrupt = result;
    corrupt.samples.pressures[0].negativeSidePressurePascals += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningLoadEpochIntegrity(corrupt);
        },
        "regional opening epoch: nested sample corruption rejects");

    auto foreignSettings =
        SceneFluidRegionalOpeningLoadEpochSettings{};
    foreignSettings.transfer.momentReferenceMeters.x = 0.25;
    expectRejected(
        [&] {
            validateOpeningLoadEpoch(
                result, endpoint, fixture, foreignSettings);
        },
        "regional opening epoch: transfer-settings provenance rejects");

    const OpeningRegionalEndpoint foreignEndpoint(true);
    expectRejected(
        [&] {
            validateOpeningLoadEpoch(
                result, foreignEndpoint, fixture);
        },
        "regional opening epoch: foreign accepted-flow source rejects");

    SceneFixture limitedFixture(partialOpeningPressureScene(), false);
    const auto limitedBefore = limitedFixture.structure.checkpoint();
    auto limits = SceneFluidRegionalOpeningLoadEpochLimits{};
    limits.maximumOwnedBytes = result.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                applyOpeningLoadEpoch(endpoint, limitedFixture, {}, limits));
        },
        "regional opening epoch: late aggregate limit rejects");
    check(samePublicCheckpoint(
              limitedBefore, limitedFixture.structure.checkpoint()),
          "regional opening epoch: late rejection restores the exact Structure checkpoint");

    SceneFixture applicationLimited(
        partialOpeningPressureScene(), false);
    const auto applicationBefore =
        applicationLimited.structure.checkpoint();
    limits = {};
    limits.application.maximumNodeLoads =
        applicationLimited.transfer.nodes().size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningLoadEpoch(
                endpoint, applicationLimited, {}, limits));
        },
        "regional opening epoch: nested application limit rejects");
    check(samePublicCheckpoint(
              applicationBefore,
              applicationLimited.structure.checkpoint()),
          "regional opening epoch: nested rejection preserves Structure");
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
        testOpeningAwareSamplingAndApplication();
        testPartialOpeningSamplingAndApplication();
        testMovingSamplingAndPower();
        testTransactionalLoadApplication();
        testMovingLoadApplication();
        testApplicationRollbackAndLimits();
        testAtomicOpeningLoadEpoch();
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
