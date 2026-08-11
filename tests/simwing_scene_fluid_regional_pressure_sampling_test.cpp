#include "fluid/planar_region_fragment_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_pressure_state.h"
#include "fluid/planar_region_fragment_projection_energy.h"
#include "fluid/planar_region_fragment_pressure_jump_energy.h"
#include "fluid/planar_region_fragment_surface_load.h"
#include "fluid/planar_region_fragment_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state_persistence.h"
#include "fluid/planar_region_fragment_opening_momentum_prediction.h"
#include "fluid/planar_region_fragment_opening_velocity_metric.h"
#include "fluid/planar_region_fragment_opening_velocity_state.h"
#include "scene_fluid_regional_pressure_sampling.h"
#include "scene_fluid_mimetic_geometry_epoch_transition.h"
#include "scene_fluid_mimetic_pressure_audit.h"
#include "scene_fluid_mimetic_pressure_flow.h"
#include "scene_fluid_opening_flux.h"
#include "scene_fluid_pressure_volume_rate.h"
#include "scene_fluid_region_momentum.h"
#include "scene_fluid_region_rebase.h"
#include "scene_fluid_region_transport.h"
#include "scene_fluid_region_wall.h"
#include "scene_fluid_regional_opening_load_epoch.h"
#include "scene_fluid_regional_opening_momentum_load_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_exchange.h"
#include "scene_fluid_regional_opening_momentum_wall_cycle_owner.h"
#include "scene_fluid_regional_opening_momentum_wall_coupled_state.h"
#include "scene_fluid_regional_opening_momentum_wall_coupled_state_persistence.h"
#include "scene_fluid_regional_opening_momentum_wall_cycle_state_persistence.h"
#include "scene_fluid_regional_opening_momentum_wall_input.h"
#include "scene_fluid_regional_opening_momentum_wall_load_application.h"
#include "scene_fluid_regional_opening_momentum_wall_load_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_post_step_geometry.h"
#include "scene_fluid_regional_opening_momentum_wall_pressure_epoch.h"
#include "scene_fluid_regional_opening_momentum_wall_structure_step_epoch.h"
#include "scene_pressure_cell_geometry.h"

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

void refreshOpeningMomentumWallCycleStatePersistenceChecksum(
    std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 24;
    constexpr std::size_t checksumOffset = 16;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t value = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        value ^= bytes[index];
        value *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            value >> (8U * byte));
    }
}

void refreshOpeningMomentumWallCoupledStatePersistenceChecksum(
    std::vector<std::uint8_t>& bytes) {
    refreshOpeningMomentumWallCycleStatePersistenceChecksum(bytes);
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
    PlanarPressureRegionFragmentVelocityMetric baseMetric;
    PlanarPressureRegionFragmentOpeningVelocityMetric openingMetric;
    PlanarPressureRegionFragmentOpeningVelocityState openingVelocity;
    PlanarPressureRegionFragmentOpeningMomentumCycleResult momentumCycle;
    PlanarPressureRegionFragmentOpeningMomentumCycleState momentumCycleState;
    PlanarPressureRegionFragmentOpeningPressureState pressure;
    PlanarPressureRegionFragmentSurfaceLoadLedger surfaceLoads;
    PlanarPressureRegionFragmentOpeningSurfaceLoadLedger openingSurfaceLoads;
    PlanarPressureRegionFragmentOpeningLoadState loadState;

    explicit OpeningRegionalEndpoint(
        const bool fullyOpenFirstSurface = true,
        const Vector3 fluidVelocityMetersPerSecond = {}) {
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
        const auto component = [&](const GridFaceAxis axis) {
            switch (axis) {
            case GridFaceAxis::X:
                return fluidVelocityMetersPerSecond.x;
            case GridFaceAxis::Y:
                return fluidVelocityMetersPerSecond.y;
            case GridFaceAxis::Z:
                return fluidVelocityMetersPerSecond.z;
            }
            throw std::logic_error(
                "opening regional fixture axis is invalid");
        };
        for (const auto& link : topology.links) {
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
                velocity[link.linkIndex] = component(link.axis);
            }
        }
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
        baseMetric = buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, sweep, fragments, topology);
        openingMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, sweep, fragments, topology, baseMetric,
                openingDefinitions, openings);
        openingVelocity =
            capturePlanarPressureRegionFragmentOpeningVelocityState(
                acceptedFlow, pressureOperator, basePressureOperator,
                geometry, sweep, fragments, topology, volumeRates,
                openingDefinitions, openings, resistanceDefinitions,
                baseMetric, openingMetric);
        momentumCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                openingVelocity, openingMetric, acceptedFlow,
                pressureOperator, basePressureOperator, geometry, sweep,
                fragments, topology, volumeRates, openingDefinitions,
                openings, resistanceDefinitions, baseMetric, openingMetric,
                pressureOperator, basePressureOperator, sweep, fragments,
                topology, volumeRates, openingDefinitions, openings,
                resistanceDefinitions, baseMetric, openingMetric, {},
                settings);
        if (!momentumCycle.diagnostics.accepted) {
            throw std::runtime_error(
                "opening regional fixture momentum cycle did not converge");
        }
        momentumCycleState =
            capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
                momentumCycle, openingMetric, openingMetric);
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

SceneFluidRegionalOpeningMomentumLoadEpoch applyOpeningMomentumLoadEpoch(
    const PlanarPressureRegionFragmentOpeningMomentumCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits = {}) {
    return applySceneFluidRegionalOpeningMomentumLoadEpoch(
        cycleState, endpoint.volumeRates,
        endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.openingMetric,
        scene.surface.definition, scene.state, scene.transfer,
        scene.quadrature, scene.structure, settings, limits);
}

SceneFluidRegionalOpeningMomentumLoadEpoch applyOpeningMomentumLoadEpoch(
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits = {}) {
    return applyOpeningMomentumLoadEpoch(
        endpoint.momentumCycleState, endpoint, scene, settings, limits);
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

void validateOpeningMomentumLoadEpoch(
    const SceneFluidRegionalOpeningMomentumLoadEpoch& result,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumLoadEpochLimits& limits = {}) {
    validateSceneFluidRegionalOpeningMomentumLoadEpoch(
        result, endpoint.momentumCycleState, endpoint.volumeRates,
        endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.openingMetric,
        scene.surface.definition, scene.state, scene.transfer,
        scene.quadrature, settings, limits);
}

SceneFluidRegionalOpeningMomentumWallInput captureOpeningMomentumWallInput(
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits = {}) {
    return captureSceneFluidRegionalOpeningMomentumWallInput(
        transport, endpoint.acceptedFlow, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.quadrature, settings, limits);
}

SceneFluidRegionalOpeningMomentumWallLoadApplication
applyOpeningMomentumWallLoads(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings = {},
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits =
        {}) {
    return applySceneFluidRegionalOpeningMomentumWallLoads(
        cycleState, endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.state, scene.transfer,
        scene.quadrature, scene.structure, settings, limits);
}

void validateOpeningMomentumWallLoads(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings = {},
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits =
        {}) {
    validateSceneFluidRegionalOpeningMomentumWallLoadApplication(
        application, cycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.state, scene.transfer,
        scene.quadrature, settings, limits);
}

SceneFluidRegionalOpeningMomentumWallLoadEpoch
applyOpeningMomentumWallLoadEpoch(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings =
        {},
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits = {}) {
    return applySceneFluidRegionalOpeningMomentumWallLoadEpoch(
        cycleState, endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.transfer, scene.quadrature, scene.structure, settings, limits);
}

void validateOpeningMomentumWallLoadEpoch(
    const SceneFluidRegionalOpeningMomentumWallLoadEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallLoadEpochSettings& settings =
        {},
    const SceneFluidRegionalOpeningMomentumWallLoadEpochLimits& limits = {}) {
    validateSceneFluidRegionalOpeningMomentumWallLoadEpoch(
        epoch, cycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.transfer, scene.quadrature, settings, limits);
}

SceneFluidRegionalOpeningMomentumWallStructureStepEpoch
advanceOpeningMomentumWallStructureStepEpoch(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits = {}) {
    return advanceSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
        cycleState, endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.transfer, scene.quadrature, scene.structure, settings, limits);
}

void validateOpeningMomentumWallStructureStepEpoch(
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpoch& epoch,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits&
        limits = {}) {
    validateSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
        epoch, cycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.transfer, scene.quadrature, scene.structure, settings, limits);
}

void validateOpeningMomentumWallInput(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const PlanarPressureRegionFragmentOpeningMomentumTransport& transport,
    const OpeningRegionalEndpoint& endpoint,
    const SceneFixture& scene,
    const SceneFluidRegionalOpeningMomentumWallInputSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallInputLimits& limits = {}) {
    validateSceneFluidRegionalOpeningMomentumWallInput(
        input, transport, endpoint.acceptedFlow, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, scene.surface.definition, scene.state,
        scene.quadrature, settings, limits);
}

void testMimeticGeometryEpochTransition() {
    const auto scene = makeScenePressureCellGeometry();
    const auto surface = assembleSceneFluidSurface(scene);
    const auto assembly = makeScenePressureCellAssembly(scene);
    check(surface.ok() && assembly.ok(),
          "mimetic geometry transition: open-cell scene assembles");
    Structure structure(assembly.definition);
    SceneFluidSurfaceTransfer transfer(
        surface.definition, assembly.mappings, structure);
    const auto geometryGrid = makeScenePressureCellGrid({4, 4, 4});
    const auto connectivity = buildSceneFluidRegionConnectivity(
        surface.definition);
    const auto previousState = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    const auto previousEpoch = buildSceneFluidMimeticGeometryEpoch(
        surface.definition, previousState, geometryGrid, transfer,
        connectivity);
    MacVelocityField predictedVelocity(geometryGrid);
    std::ranges::fill(predictedVelocity.xFaces(), 0.1);
    const auto previousOpeningFlux = evaluateSceneFluidOpeningFlux(
        surface.definition, previousState, previousEpoch.openingCaps,
        previousEpoch.openingQuadrature, previousEpoch.openingPatches,
        geometryGrid, predictedVelocity);
    SceneFluidMimeticPressureAuditSettings pressureSettings;
    pressureSettings.timeStepSeconds = 0.01;
    pressureSettings.pressureSolve.maximumIterations = 4000;
    pressureSettings.pressureSolve.relativeResidualTolerance = 1.0e-5;
    pressureSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    pressureSettings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-8;
    const auto previousPressure =
        buildSceneFluidMimeticPressureAuditEndpoint(
            surface.definition, previousState, geometryGrid,
            previousEpoch.gridEpoch, previousEpoch.openingCaps,
            previousEpoch.openingQuadrature,
            previousEpoch.openingPatches,
            previousEpoch.pressureControlVolumes,
            previousEpoch.pressureFaceLinks, previousOpeningFlux,
            predictedVelocity, pressureSettings);
    const auto previousCorrectedFlow =
        correctSceneFluidMimeticTraceFlows(previousPressure);
    const auto previousCollapsedVelocity =
        collapseSceneFluidMimeticCorrectedMacVelocity(
            previousCorrectedFlow, previousEpoch.pressureFaceLinks,
            previousEpoch.openingPatches, geometryGrid);
    const auto previousMomentum = reconstructSceneFluidRegionMomentumState(
        geometryGrid, previousEpoch.pressureControlVolumes,
        previousEpoch.pressureFaceLinks, previousEpoch.openingPatches,
        previousPressure.controlCells, previousCorrectedFlow,
        previousCollapsedVelocity.velocityMetersPerSecond);
    SceneFluidRegionTransportSettings transportSettings;
    transportSettings.timeStepSeconds = 0.01;
    const auto transportedMomentum = advanceSceneFluidRegionMomentum(
        previousMomentum, previousEpoch.pressureFaceLinks,
        previousPressure.controlCells, previousCorrectedFlow,
        transportSettings);

    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = 0.01;
    stepSettings.substeps = 1;
    stepSettings.constraintIterations = 8;
    stepSettings.gravityMetersPerSecondSquared = {0.1, 0.0, 0.0};
    stepSettings.velocityDampingPerSecond = 0.0;
    const auto structureDiagnostics = structure.step(stepSettings);
    const auto currentState = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    const auto acceptedGridEpoch = buildSceneFluidGridEpoch(
        surface.definition, currentState, geometryGrid, transfer);
    const auto transition =
        buildSceneFluidMimeticGeometryEpochTransition(
            previousEpoch, surface.definition, previousState,
            currentState, geometryGrid, transfer, connectivity,
            acceptedGridEpoch);
    const auto volumeRates = buildSceneFluidPressureVolumeRates(
        previousEpoch.cellVolumes,
        transition.currentGeometryEpoch.cellVolumes,
        transition.currentGeometryEpoch.pressureControlVolumes,
        transition.topologyTransition);
    const auto regionalRebase = rebaseSceneFluidRegionTransport(
        transportedMomentum, previousEpoch.pressureControlVolumes,
        transition.currentGeometryEpoch.pressureControlVolumes,
        transition.topologyTransition);
    const auto currentOpeningFlux = evaluateSceneFluidOpeningFlux(
        surface.definition, currentState,
        transition.currentGeometryEpoch.openingCaps,
        transition.currentGeometryEpoch.openingQuadrature,
        transition.currentGeometryEpoch.openingPatches, geometryGrid,
        predictedVelocity);
    SceneFluidRegionWallSettings wallSettings;
    wallSettings.timeStepSeconds = volumeRates.durationSeconds;
    const auto wallExchange = exchangeSceneFluidRegionWallMomentum(
        regionalRebase, geometryGrid,
        transition.currentGeometryEpoch.pressureControlVolumes,
        surface.definition, currentState,
        transition.currentGeometryEpoch.gridEpoch.quadrature,
        wallSettings);
    const auto currentPressure =
        buildSceneFluidMimeticPressureAuditEndpoint(
            surface.definition, currentState, geometryGrid,
            transition.currentGeometryEpoch, currentOpeningFlux,
            wallExchange, volumeRates, transition.topologyTransition,
            previousPressure, pressureSettings);
    validateSceneFluidMimeticGeometryEpochTransition(
        transition, previousEpoch, surface.definition, previousState,
        currentState, geometryGrid, transfer, connectivity,
        acceptedGridEpoch);
    check(structureDiagnostics.finite
              && currentState != previousState
              && transition.version
                  == sceneFluidMimeticGeometryEpochTransitionVersion
              && transition.fingerprint != 0
              && transition.previousGeometryEpochFingerprint
                  == previousEpoch.fingerprint
              && transition.currentGeometryEpoch.gridEpoch
                  == acceptedGridEpoch
              && transition.controlVolumeTopologyStable
              && transition.topologyTransition
                     .appearedControlVolumeCount == 0
              && transition.topologyTransition
                     .disappearedControlVolumeCount == 0,
          "mimetic geometry transition: moving authored intake retains a consecutive graph-free topology");
    check(previousPressure.pressureEpoch.diagnostics.accepted
              && previousCorrectedFlow.accepted
              && transportedMomentum.diagnostics.accepted
              && regionalRebase.diagnostics.finite
              && wallExchange.diagnostics.accepted
              && currentPressure.pressureEpoch.diagnostics.accepted
              && currentPressure.usesRegionWallPrediction
              && currentPressure.usesConsecutiveWarmStart
              && currentPressure.scenePressureEpochFingerprint
                  == transition.currentGeometryEpoch.fingerprint
              && currentPressure.pressureTopologyTransitionFingerprint
                  == transition.topologyTransition.fingerprint,
          "mimetic geometry transition: regional rebase and wall-adjusted flow accept the next moving pressure endpoint");
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

    SceneFixture momentumFixture(partialOpeningPressureScene(), false);
    SceneFixture repeatedMomentumFixture(
        partialOpeningPressureScene(), false);
    const auto momentumBefore = momentumFixture.structure.checkpoint();
    const auto momentumResult =
        applyOpeningMomentumLoadEpoch(endpoint, momentumFixture);
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError
        momentumPersistenceError;
    std::vector<std::uint8_t> momentumStateBytes;
    check(serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
              endpoint.momentumCycleState, endpoint.volumeRates,
              endpoint.openingMetric, endpoint.pressureOperator,
              endpoint.basePressureOperator, endpoint.geometry,
              endpoint.sweep, endpoint.fragments, endpoint.topology,
              endpoint.volumeRates, endpoint.openingDefinitions,
              endpoint.openings, endpoint.resistanceDefinitions,
              endpoint.openingMetric, momentumStateBytes,
              &momentumPersistenceError),
          "regional opening momentum epoch: source cycle serializes before restart replay");
    PlanarPressureRegionFragmentOpeningMomentumCycleState
        restoredMomentumCycleState;
    check(deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
              momentumStateBytes, endpoint.volumeRates,
              endpoint.openingMetric, endpoint.pressureOperator,
              endpoint.basePressureOperator, endpoint.geometry,
              endpoint.sweep, endpoint.fragments, endpoint.topology,
              endpoint.volumeRates, endpoint.openingDefinitions,
              endpoint.openings, endpoint.resistanceDefinitions,
              endpoint.openingMetric, restoredMomentumCycleState,
              &momentumPersistenceError)
              && restoredMomentumCycleState
                  == endpoint.momentumCycleState,
          "regional opening momentum epoch: source cycle restores bit-exactly");
    const auto repeatedMomentumResult =
        applyOpeningMomentumLoadEpoch(
            restoredMomentumCycleState, endpoint,
            repeatedMomentumFixture);
    check(momentumResult == repeatedMomentumResult
              && momentumResult.version
                  == sceneFluidRegionalOpeningMomentumLoadEpochVersion
              && momentumResult.fingerprint != 0
              && momentumResult.applied
              && momentumResult.sourceCycleStateFingerprint
                  == endpoint.momentumCycleState.fingerprint
              && momentumResult.sourceTransportFingerprint
                  == endpoint.momentumCycleState.transport.fingerprint
              && momentumResult.sourceAcceptedStateFingerprint
                  == endpoint.momentumCycleState.acceptedState.fingerprint
              && momentumResult.loadEpoch.applied
              && momentumResult.loadEpoch.loadState.acceptedFlow
                  == endpoint.momentumCycleState.acceptedState
              && momentumResult.workingStorageBytes
                  == endpoint.momentumCycleState.ownedStorageBytes
                      + momentumResult.loadEpoch.workingStorageBytes,
          "regional opening momentum epoch: transported provenance and Structure application commit together deterministically");
    check(momentumFixture.structure.checkpoint()
                  .pendingExternalForcesNewtons
              != momentumBefore.pendingExternalForcesNewtons,
          "regional opening momentum epoch: transported pressure changes pending Structure loads");
    validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(
        momentumResult);
    validateOpeningMomentumLoadEpoch(
        momentumResult, endpoint, momentumFixture);

    auto corruptMomentumResult = momentumResult;
    ++corruptMomentumResult.sourceTransportFingerprint;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumLoadEpochIntegrity(
                corruptMomentumResult);
        },
        "regional opening momentum epoch: transport provenance corruption rejects");

    SceneFixture momentumLimited(
        partialOpeningPressureScene(), false);
    const auto momentumLimitedBefore =
        momentumLimited.structure.checkpoint();
    auto momentumLimits =
        SceneFluidRegionalOpeningMomentumLoadEpochLimits{};
    momentumLimits.maximumOwnedBytes =
        momentumResult.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningMomentumLoadEpoch(
                endpoint, momentumLimited, {}, momentumLimits));
        },
        "regional opening momentum epoch: late outer aggregate limit rejects");
    check(samePublicCheckpoint(
              momentumLimitedBefore,
              momentumLimited.structure.checkpoint()),
          "regional opening momentum epoch: late outer rejection restores exact Structure state");

    SceneFixture momentumSourceLimited(
        partialOpeningPressureScene(), false);
    const auto momentumSourceLimitedBefore =
        momentumSourceLimited.structure.checkpoint();
    momentumLimits = {};
    momentumLimits.cycleState.transport.maximumFragments =
        endpoint.momentumCycleState.transport.controls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningMomentumLoadEpoch(
                endpoint, momentumSourceLimited, {}, momentumLimits));
        },
        "regional opening momentum epoch: source-cycle limit rejects");
    check(samePublicCheckpoint(
              momentumSourceLimitedBefore,
              momentumSourceLimited.structure.checkpoint()),
          "regional opening momentum epoch: source rejection preserves Structure state");

    const OpeningRegionalEndpoint foreignMomentumEndpoint(true);
    expectRejected(
        [&] {
            validateOpeningMomentumLoadEpoch(
                momentumResult, foreignMomentumEndpoint,
                momentumFixture);
        },
        "regional opening momentum epoch: foreign transported cycle rejects");

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

void testOpeningMomentumWallInput() {
    const OpeningRegionalEndpoint endpoint(
        false, Vector3{0.0, 0.4, -0.2});
    SceneFixture fixture(partialOpeningPressureScene(), false);
    SceneFixture repeatedFixture(partialOpeningPressureScene(), false);

    const auto input = captureOpeningMomentumWallInput(
        endpoint.momentumCycle.transport, endpoint, fixture);
    PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError
        persistenceError;
    std::vector<std::uint8_t> stateBytes;
    check(serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
              endpoint.momentumCycleState, endpoint.volumeRates,
              endpoint.openingMetric, endpoint.pressureOperator,
              endpoint.basePressureOperator, endpoint.geometry,
              endpoint.sweep, endpoint.fragments, endpoint.topology,
              endpoint.volumeRates, endpoint.openingDefinitions,
              endpoint.openings, endpoint.resistanceDefinitions,
              endpoint.openingMetric, stateBytes, &persistenceError),
          "regional opening wall input: source staggered state serializes before same-epoch replay");
    PlanarPressureRegionFragmentOpeningMomentumCycleState restoredState;
    check(deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
              stateBytes, endpoint.volumeRates, endpoint.openingMetric,
              endpoint.pressureOperator, endpoint.basePressureOperator,
              endpoint.geometry, endpoint.sweep, endpoint.fragments,
              endpoint.topology, endpoint.volumeRates,
              endpoint.openingDefinitions, endpoint.openings,
              endpoint.resistanceDefinitions, endpoint.openingMetric,
              restoredState, &persistenceError),
          "regional opening wall input: source staggered state restores for same-epoch replay");
    const auto repeated = captureOpeningMomentumWallInput(
        restoredState.transport, endpoint, repeatedFixture);

    bool zeroExchange = true;
    for (const auto& sample : input.samples) {
        zeroExchange = zeroExchange
            && sample.negativeSideFluidImpulseKilogramMetersPerSecond
                == Vector3{}
            && sample.positiveSideFluidImpulseKilogramMetersPerSecond
                == Vector3{}
            && sample.structureTraction.tractionPascals
                == StructureVector3{};
    }
    double maximumTangentialSpeed = 0.0;
    for (const auto& control : input.controlVolumes) {
        maximumTangentialSpeed = std::max(
            maximumTangentialSpeed,
            std::hypot(
                control.velocityMetersPerSecond.y,
                control.velocityMetersPerSecond.z));
    }
    check(input == repeated
              && input.version
                  == sceneFluidRegionalOpeningMomentumWallInputVersion
              && input.fingerprint != 0
              && input.sourceTransportFingerprint
                  == endpoint.momentumCycle.transport.fingerprint
              && input.sourceAcceptedStateFingerprint
                  == endpoint.acceptedFlow.fingerprint
              && input.sourceTransportMetricFingerprint
                  == endpoint.openingMetric.fingerprint
              && input.sourceSamplingFingerprint != 0
              && input.quadratureFingerprint
                  == fixture.quadrature.fingerprint
              && input.controlVolumes.size()
                  == endpoint.momentumCycle.transport.controls.size()
              && input.samples.size() == fixture.quadrature.points.size()
              && input.activeControlVolumeCount > 0
              && input.maximumIncidentWallAreaSquareMeters > 0.0
              && maximumTangentialSpeed > 0.4
              && zeroExchange,
          "regional opening wall input: restored current transport maps deterministically to zero-exchange kernel descriptors");
    checkNear(input.wallSampleAreaSquareMeters, 7.5, 2.0e-13,
              "regional opening wall input: retained-solid scene area is complete");
    checkNear(input.controlIncidentWallAreaSquareMeters, 15.0, 4.0e-13,
              "regional opening wall input: both fluid sides own every retained-solid sample");
    validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(input);
    validateOpeningMomentumWallInput(
        input, endpoint.momentumCycle.transport, endpoint, fixture);

    auto corrupt = input;
    corrupt.samples.front().negativeSideControlVolumeIndex =
        corrupt.controlVolumes.size();
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(
                corrupt);
        },
        "regional opening wall input: corrupt fragment ownership rejects");

    auto limits = SceneFluidRegionalOpeningMomentumWallInputLimits{};
    limits.maximumSamples = input.samples.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(captureOpeningMomentumWallInput(
                endpoint.momentumCycle.transport, endpoint, fixture, {},
                limits));
        },
        "regional opening wall input: sample limit rejects before publication");

    const OpeningRegionalEndpoint foreignEndpoint(false);
    expectRejected(
        [&] {
            static_cast<void>(captureOpeningMomentumWallInput(
                foreignEndpoint.momentumCycle.transport, endpoint,
                fixture));
        },
        "regional opening wall input: foreign current transported flow rejects");
}

void testOpeningMomentumWallExchange() {
    const OpeningRegionalEndpoint endpoint(
        false, Vector3{0.0, 0.4, -0.2});
    SceneFixture fixture(partialOpeningPressureScene(), false);
    const auto input = captureOpeningMomentumWallInput(
        endpoint.momentumCycle.transport, endpoint, fixture);
    SceneFluidRegionWallSettings settings;
    settings.timeStepSeconds = input.timeStepSeconds;

    const auto exchange = exchangeSceneFluidRegionalOpeningMomentumWall(
        input, settings);
    const auto repeated = exchangeSceneFluidRegionalOpeningMomentumWall(
        input, settings);
    bool tangentialTractions = true;
    bool nonzeroTraction = false;
    for (const auto& sample : exchange.samples) {
        const auto& traction = sample.structureTraction.tractionPascals;
        const auto& normal = sample.unitNormalNegativeToPositive;
        tangentialTractions = tangentialTractions
            && std::abs(
                traction.x * normal.x + traction.y * normal.y
                + traction.z * normal.z) < 1.0e-12;
        nonzeroTraction = nonzeroTraction
            || traction != StructureVector3{};
    }
    check(exchange == repeated
              && exchange.version
                  == sceneFluidRegionalOpeningMomentumWallExchangeVersion
              && exchange.fingerprint != 0
              && exchange.sourceWallInputFingerprint == input.fingerprint
              && exchange.sourceInput == input
              && exchange.diagnostics.accepted
              && exchange.diagnostics.finite
              && exchange.diagnostics.controlVolumeCount
                  == input.controlVolumes.size()
              && exchange.diagnostics.quadraturePointCount
                  == input.samples.size()
              && exchange.diagnostics
                     .maximumRelativeTangentialSpeedMetersPerSecond > 0.4
              && exchange.diagnostics.viscousDissipationJoules > 0.0
              && exchange.diagnostics.kineticEnergyAfterJoules
                  < exchange.diagnostics.kineticEnergyBeforeJoules
              && exchange.diagnostics
                     .momentumResidualNormKilogramMetersPerSecond < 1.0e-12
              && exchange.controlVolumes != input.controlVolumes
              && tangentialTractions && nonzeroTraction,
          "regional opening wall exchange: shared kernel publishes deterministic dissipative tangential action/reaction");
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(exchange);
    validateSceneFluidRegionalOpeningMomentumWallExchange(
        exchange, input, settings);

    auto zeroSettings = settings;
    zeroSettings.kinematicViscositySquareMetersPerSecond = 0.0;
    const auto zeroExchange = exchangeSceneFluidRegionalOpeningMomentumWall(
        input, zeroSettings);
    check(zeroExchange.diagnostics.accepted
              && zeroExchange.controlVolumes == input.controlVolumes
              && zeroExchange.diagnostics
                     .fluidImpulseKilogramMetersPerSecond == Vector3{}
              && zeroExchange.diagnostics
                     .structureImpulseKilogramMetersPerSecond == Vector3{}
              && zeroExchange.diagnostics.viscousDissipationJoules == 0.0,
          "regional opening wall exchange: zero viscosity preserves the mapped transport exactly");

    const auto adjustment =
        captureSceneFluidRegionalOpeningMomentumAdjustmentState(
            exchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric);
    validateSceneFluidRegionalOpeningMomentumAdjustmentState(
        adjustment, exchange, endpoint.momentumCycle.transport,
        endpoint.openingMetric);
    check(adjustment.sourceTransportFingerprint
              == endpoint.momentumCycle.transport.fingerprint
              && adjustment.sourceAdjustmentFingerprint
                  == exchange.fingerprint
              && adjustment.sourceMetricFingerprint
                  == endpoint.openingMetric.fingerprint
              && adjustment.controls.size()
                  == endpoint.momentumCycle.transport.controls.size()
              && adjustment.adjustedMomentumKilogramMetersPerSecond
                  == exchange.diagnostics
                         .fluidMomentumAfterKilogramMetersPerSecond
              && adjustment.adjustmentImpulseKilogramMetersPerSecond
                  == exchange.diagnostics
                         .fluidImpulseKilogramMetersPerSecond
              && adjustment.adjustedKineticEnergyJoules
                  == exchange.diagnostics.kineticEnergyAfterJoules
              && adjustment.kineticEnergyChangeJoules < 0.0,
          "regional opening wall exchange: accepted output captures an exact fluid-owned adjusted momentum endpoint");

    const auto baselinePrediction =
        predictPlanarPressureRegionFragmentOpeningMomentum(
            endpoint.momentumCycle.transport, endpoint.openingMetric,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.baseMetric, endpoint.openingMetric);
    const auto adjustedPrediction =
        predictPlanarPressureRegionFragmentOpeningMomentum(
            adjustment, endpoint.openingMetric, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.baseMetric,
            endpoint.openingMetric);
    validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
        adjustedPrediction, adjustment, endpoint.openingMetric,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.baseMetric, endpoint.openingMetric);
    check(adjustedPrediction.sourceTransportFingerprint
              == baselinePrediction.sourceTransportFingerprint
              && adjustedPrediction.sourceAdjustmentStateFingerprint
                  == adjustment.fingerprint
              && adjustedPrediction.predictedVelocityState
                  != baselinePrediction.predictedVelocityState
              && adjustedPrediction.diagnostics.sourceMomentumKilogramMetersPerSecond
                  == adjustment.adjustedMomentumKilogramMetersPerSecond
              && adjustedPrediction.diagnostics.sourceKineticEnergyJoules
                  == adjustment.adjustedKineticEnergyJoules,
          "regional opening wall exchange: opt-in prediction consumes adjusted tangential momentum with exact provenance");

    const auto zeroAdjustment =
        captureSceneFluidRegionalOpeningMomentumAdjustmentState(
            zeroExchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric);
    const auto zeroPrediction =
        predictPlanarPressureRegionFragmentOpeningMomentum(
            zeroAdjustment, endpoint.openingMetric, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.baseMetric,
            endpoint.openingMetric);
    check(zeroAdjustment.controls
              == endpoint.momentumCycle.transport.controls
              && zeroPrediction.sourceAdjustmentStateFingerprint
                  == zeroAdjustment.fingerprint
              && zeroPrediction.predictedVelocityState
                  == baselinePrediction.predictedVelocityState
              && zeroPrediction.diagnostics
                  == baselinePrediction.diagnostics
              && zeroPrediction.fingerprint != baselinePrediction.fingerprint,
          "regional opening wall exchange: zero viscosity preserves the exact numeric prediction while retaining adjustment lineage");

    const auto baselinePressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
            baselinePrediction, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            endpoint.momentumCycle.pressureSettings);
    const auto adjustedPressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
            adjustedPrediction, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            endpoint.momentumCycle.pressureSettings);
    const auto zeroPressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
            zeroPrediction, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            endpoint.momentumCycle.pressureSettings);
    validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
        adjustedPressureEpoch, adjustedPrediction,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, endpoint.momentumCycle.pressureSettings);
    check(baselinePressureEpoch.diagnostics.accepted
              && adjustedPressureEpoch.diagnostics.accepted
              && zeroPressureEpoch.diagnostics.accepted
              && adjustedPressureEpoch.sourcePredictionFingerprint
                  == adjustedPrediction.fingerprint
              && adjustedPressureEpoch.acceptedState
                  != baselinePressureEpoch.acceptedState
              && zeroPressureEpoch.acceptedState
                  == baselinePressureEpoch.acceptedState
              && zeroPressureEpoch.diagnostics
                  == baselinePressureEpoch.diagnostics
              && zeroPressureEpoch.sourcePredictionFingerprint
                  != baselinePressureEpoch.sourcePredictionFingerprint,
          "regional opening wall exchange: opt-in adjusted prediction reaches pressure acceptance while zero viscosity preserves the exact numeric endpoint");

    const auto wallPressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            exchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric, {},
            endpoint.momentumCycle.pressureSettings);
    const auto zeroWallPressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            zeroExchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric, {},
            endpoint.momentumCycle.pressureSettings);
    validateSceneFluidRegionalOpeningMomentumWallPressureEpoch(
        wallPressureEpoch, exchange, endpoint.momentumCycle.transport,
        endpoint.openingMetric, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric);
    check(wallPressureEpoch.fingerprint != 0
              && wallPressureEpoch.sourceWallExchangeFingerprint
                  == exchange.fingerprint
              && wallPressureEpoch.adjustmentState == adjustment
              && wallPressureEpoch.prediction == adjustedPrediction
              && wallPressureEpoch.pressureEpoch == adjustedPressureEpoch
              && wallPressureEpoch.pressureEpoch.diagnostics.accepted
              && zeroWallPressureEpoch.adjustmentState == zeroAdjustment
              && zeroWallPressureEpoch.prediction == zeroPrediction
              && zeroWallPressureEpoch.pressureEpoch
                  == zeroPressureEpoch
              && zeroWallPressureEpoch.pressureEpoch.acceptedState
                  == baselinePressureEpoch.acceptedState,
          "regional opening wall pressure epoch: one immutable receipt retains the exact adjustment-prediction-pressure chain");

    const auto pressureWarmStart =
        buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
            baselinePressureEpoch.acceptedState,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.sweep,
            endpoint.fragments, endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings);
    const auto directWarmPressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
            adjustedPrediction, pressureWarmStart,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric,
            endpoint.momentumCycle.pressureSettings);
    const auto warmWallPressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            exchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric, pressureWarmStart,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric, {},
            endpoint.momentumCycle.pressureSettings);
    validateSceneFluidRegionalOpeningMomentumWallPressureEpoch(
        warmWallPressureEpoch, exchange,
        endpoint.momentumCycle.transport, endpoint.openingMetric,
        pressureWarmStart, endpoint.pressureOperator,
        endpoint.basePressureOperator, endpoint.geometry, endpoint.sweep,
        endpoint.fragments, endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric);
    check(warmWallPressureEpoch.sourcePressureWarmStartFingerprint
              == pressureWarmStart.fingerprint
              && warmWallPressureEpoch.pressureEpoch
                  == directWarmPressureEpoch
              && warmWallPressureEpoch.pressureEpoch.diagnostics.accepted
              && warmWallPressureEpoch.pressureEpoch.diagnostics
                     .usedWarmPressureStart,
          "regional opening wall pressure epoch: pressure-only warm lineage composes without carrying stale flow");

    const auto baselineTargetFlow =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            baselinePressureEpoch.acceptedState,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric);
    const auto adjustedTargetFlow =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            adjustedPressureEpoch.acceptedState,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric);
    const auto zeroTargetFlow =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            zeroPressureEpoch.acceptedState, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric);
    const auto baselineNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            endpoint.momentumCycle.transport, endpoint.openingMetric,
            baselineTargetFlow, endpoint.openingMetric,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    const auto adjustedNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            adjustment, endpoint.openingMetric, adjustedTargetFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology,
            endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    const auto zeroNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            zeroAdjustment, endpoint.openingMetric, zeroTargetFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology,
            endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    validatePlanarPressureRegionFragmentOpeningMomentumTransport(
        adjustedNextTransport, adjustment, endpoint.openingMetric,
        adjustedTargetFlow, endpoint.openingMetric, endpoint.geometry,
        endpoint.sweep, endpoint.fragments, endpoint.topology,
        endpoint.volumeRates);
    check(baselineNextTransport.diagnostics.accepted
              && adjustedNextTransport.diagnostics.accepted
              && zeroNextTransport.diagnostics.accepted
              && baselineNextTransport.sourceTransportFingerprint
                  == endpoint.momentumCycle.transport.fingerprint
              && baselineNextTransport.sourceAdjustmentStateFingerprint
                  == 0
              && adjustedNextTransport.sourceTransportFingerprint == 0
              && adjustedNextTransport.sourceAdjustmentStateFingerprint
                  == adjustment.fingerprint
              && adjustedNextTransport.controls
                  != baselineNextTransport.controls
              && zeroNextTransport.controls
                  == baselineNextTransport.controls
              && zeroNextTransport.diagnostics
                  == baselineNextTransport.diagnostics
              && zeroNextTransport.fingerprint
                  != baselineNextTransport.fingerprint,
          "regional opening wall pressure epoch: adjusted collocated momentum propagates into the next ALE transport with zero-viscosity numeric identity");

    auto ambiguousNextTransport = adjustedNextTransport;
    ambiguousNextTransport.sourceTransportFingerprint =
        endpoint.momentumCycle.transport.fingerprint;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
                ambiguousNextTransport);
        },
        "regional opening wall pressure epoch: ambiguous adjusted/raw transport lineage rejects");
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningMomentumTransport(
                adjustedNextTransport, zeroAdjustment,
                endpoint.openingMetric, zeroTargetFlow,
                endpoint.openingMetric, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates);
        },
        "regional opening wall pressure epoch: foreign adjustment rejects next-transport validation");

    auto nextTransportLimits =
        PlanarPressureRegionFragmentOpeningMomentumTransportLimits{};
    nextTransportLimits.maximumFragments =
        adjustment.controls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                advancePlanarPressureRegionFragmentOpeningMomentum(
                    adjustment, endpoint.openingMetric,
                    adjustedTargetFlow, endpoint.openingMetric,
                    endpoint.geometry, endpoint.sweep,
                    endpoint.fragments, endpoint.topology,
                    endpoint.volumeRates,
                    endpoint.momentumCycle.transportSettings,
                    nextTransportLimits));
        },
        "regional opening wall pressure epoch: adjusted-source next-transport limit rejects");

    const auto wallCycleState =
        captureSceneFluidRegionalOpeningMomentumWallCycleState(
            wallPressureEpoch, exchange, endpoint.openingMetric,
            endpoint.openingMetric, fixture.quadrature);
    validateSceneFluidRegionalOpeningMomentumWallCycleState(
        wallCycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, fixture.quadrature);
    check(wallCycleState.sourceWallPressureEpochFingerprint
              == wallPressureEpoch.fingerprint
              && wallCycleState.sourceWallExchangeFingerprint
                  == exchange.fingerprint
              && wallCycleState.adjustedMomentum == adjustment
              && wallCycleState.acceptedPressure
                  == adjustedPressureEpoch.acceptedState
              && wallCycleState.wallTractions.wallExchangeFingerprint
                  == exchange.fingerprint
              && wallCycleState.wallTractions.tractions.size()
                  == exchange.samples.size(),
          "regional opening wall-cycle state: compact endpoint retains adjusted momentum, accepted pressure, and matching traction");

    using WallCyclePersistenceErrorCode =
        SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceErrorCode;
    using WallCyclePersistenceError =
        SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceError;
    using WallCyclePersistenceLimits =
        SceneFluidRegionalOpeningMomentumWallCycleStatePersistenceLimits;
    const auto serializeWallCycleState =
        [&](const SceneFluidRegionalOpeningMomentumWallCycleState& source,
            std::vector<std::uint8_t>& destination,
            WallCyclePersistenceError* error,
            const WallCyclePersistenceLimits& limits) {
            return serializeSceneFluidRegionalOpeningMomentumWallCycleState(
                source, endpoint.openingMetric,
                endpoint.pressureOperator,
                endpoint.basePressureOperator, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates, endpoint.openingDefinitions,
                endpoint.openings, endpoint.resistanceDefinitions,
                endpoint.baseMetric, endpoint.openingMetric,
                fixture.quadrature, destination, error, limits);
        };
    const auto deserializeWallCycleState =
        [&](const std::span<const std::uint8_t> source,
            SceneFluidRegionalOpeningMomentumWallCycleState& destination,
            WallCyclePersistenceError* error,
            const WallCyclePersistenceLimits& limits) {
            return deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
                source, endpoint.openingMetric,
                endpoint.pressureOperator,
                endpoint.basePressureOperator, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates, endpoint.openingDefinitions,
                endpoint.openings, endpoint.resistanceDefinitions,
                endpoint.baseMetric, endpoint.openingMetric,
                fixture.quadrature, destination, error, limits);
        };

    WallCyclePersistenceError wallCyclePersistenceError;
    const WallCyclePersistenceLimits defaultWallCyclePersistenceLimits;
    std::vector<std::uint8_t> wallCycleBytes;
    std::vector<std::uint8_t> repeatedWallCycleBytes;
    check(serializeWallCycleState(
              wallCycleState, wallCycleBytes,
              &wallCyclePersistenceError,
              defaultWallCyclePersistenceLimits)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::None
              && serializeWallCycleState(
                  wallCycleState, repeatedWallCycleBytes,
                  &wallCyclePersistenceError,
                  defaultWallCyclePersistenceLimits)
              && wallCycleBytes == repeatedWallCycleBytes
              && wallCycleBytes.size() > 24
              && wallCycleBytes[0] == 'S'
              && wallCycleBytes[1] == 'W'
              && wallCycleBytes[2] == 'R'
              && wallCycleBytes[3] == 'W',
          "regional opening wall-cycle persistence: SWRW encoding is deterministic");
    SceneFluidRegionalOpeningMomentumWallCycleState decodedWallCycleState;
    std::vector<std::uint8_t> decodedWallCycleBytes;
    check(deserializeWallCycleState(
              wallCycleBytes, decodedWallCycleState,
              &wallCyclePersistenceError,
              defaultWallCyclePersistenceLimits)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::None
              && decodedWallCycleState == wallCycleState
              && serializeWallCycleState(
                  decodedWallCycleState, decodedWallCycleBytes,
                  &wallCyclePersistenceError,
                  defaultWallCyclePersistenceLimits)
              && decodedWallCycleBytes == wallCycleBytes,
          "regional opening wall-cycle persistence: adjusted momentum, pressure, and traction round trip bit-exactly");

    const auto expectWallCyclePersistenceRejected =
        [&](const std::vector<std::uint8_t>& candidateBytes,
            const WallCyclePersistenceErrorCode expectedCode,
            const WallCyclePersistenceLimits& limits) {
            auto retainedState = wallCycleState;
            const auto before = retainedState;
            WallCyclePersistenceError rejectedError;
            const bool decoded = deserializeWallCycleState(
                candidateBytes, retainedState, &rejectedError, limits);
            return !decoded && rejectedError.code == expectedCode
                && retainedState == before;
        };
    auto corruptWallCycleBytes = wallCycleBytes;
    corruptWallCycleBytes[0] ^= 0x01U;
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::InvalidMagic,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: foreign magic rejects transactionally");
    corruptWallCycleBytes = wallCycleBytes;
    corruptWallCycleBytes[4] = 2;
    corruptWallCycleBytes[5] = 0;
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::UnsupportedVersion,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: unsupported protocol rejects transactionally");
    corruptWallCycleBytes = wallCycleBytes;
    corruptWallCycleBytes.back() ^= 0x01U;
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::ChecksumMismatch,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: payload corruption is detected");
    corruptWallCycleBytes = wallCycleBytes;
    constexpr std::size_t encodedWallCycleStateFingerprintOffset = 40;
    corruptWallCycleBytes[encodedWallCycleStateFingerprintOffset] ^= 0x01U;
    refreshOpeningMomentumWallCycleStatePersistenceChecksum(
        corruptWallCycleBytes);
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::InvalidData,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: recomputed-checksum state corruption rejects");
    corruptWallCycleBytes = wallCycleBytes;
    corruptWallCycleBytes.pop_back();
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::Truncated,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: truncation rejects transactionally");
    corruptWallCycleBytes = wallCycleBytes;
    corruptWallCycleBytes.push_back(0);
    check(expectWallCyclePersistenceRejected(
              corruptWallCycleBytes,
              WallCyclePersistenceErrorCode::TrailingData,
              defaultWallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: trailing bytes reject transactionally");

    auto wallCyclePersistenceLimits = defaultWallCyclePersistenceLimits;
    wallCyclePersistenceLimits.maximumAdjustmentControls =
        wallCycleState.adjustedMomentum.controls.size() - 1;
    check(expectWallCyclePersistenceRejected(
              wallCycleBytes,
              WallCyclePersistenceErrorCode::LimitExceeded,
              wallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: adjusted-control limit rejects");
    wallCyclePersistenceLimits = defaultWallCyclePersistenceLimits;
    wallCyclePersistenceLimits.maximumWallTractions =
        wallCycleState.wallTractions.tractions.size() - 1;
    check(expectWallCyclePersistenceRejected(
              wallCycleBytes,
              WallCyclePersistenceErrorCode::LimitExceeded,
              wallCyclePersistenceLimits),
          "regional opening wall-cycle persistence: traction-record limit rejects");

    SceneFixture foreignWallCycleFixture(false);
    auto foreignDecodedWallCycleState = wallCycleState;
    const auto beforeForeignDecode = foreignDecodedWallCycleState;
    check(!deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
              wallCycleBytes, endpoint.openingMetric,
              endpoint.pressureOperator,
              endpoint.basePressureOperator, endpoint.geometry,
              endpoint.sweep, endpoint.fragments, endpoint.topology,
              endpoint.volumeRates, endpoint.openingDefinitions,
              endpoint.openings, endpoint.resistanceDefinitions,
              endpoint.baseMetric, endpoint.openingMetric,
              foreignWallCycleFixture.quadrature,
              foreignDecodedWallCycleState, &wallCyclePersistenceError)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::SourceMismatch
              && foreignDecodedWallCycleState == beforeForeignDecode,
          "regional opening wall-cycle persistence: foreign quadrature rejects transactionally");
    const OpeningRegionalEndpoint foreignWallPressureEndpoint(true);
    foreignDecodedWallCycleState = wallCycleState;
    const auto beforeForeignPressureDecode = foreignDecodedWallCycleState;
    check(!deserializeSceneFluidRegionalOpeningMomentumWallCycleState(
              wallCycleBytes, endpoint.openingMetric,
              foreignWallPressureEndpoint.pressureOperator,
              foreignWallPressureEndpoint.basePressureOperator,
              foreignWallPressureEndpoint.geometry,
              foreignWallPressureEndpoint.sweep,
              foreignWallPressureEndpoint.fragments,
              foreignWallPressureEndpoint.topology,
              foreignWallPressureEndpoint.volumeRates,
              foreignWallPressureEndpoint.openingDefinitions,
              foreignWallPressureEndpoint.openings,
              foreignWallPressureEndpoint.resistanceDefinitions,
              endpoint.baseMetric, endpoint.openingMetric,
              fixture.quadrature, foreignDecodedWallCycleState,
              &wallCyclePersistenceError)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::SourceMismatch
              && foreignDecodedWallCycleState
                  == beforeForeignPressureDecode,
          "regional opening wall-cycle persistence: foreign accepted-pressure sources reject transactionally");

    wallCyclePersistenceLimits = defaultWallCyclePersistenceLimits;
    wallCyclePersistenceLimits.maximumEncodedBytes =
        wallCycleBytes.size() - 1;
    std::vector<std::uint8_t> retainedWallCycleEncoding{1, 2, 3};
    const auto retainedWallCycleEncodingBefore =
        retainedWallCycleEncoding;
    check(!serializeWallCycleState(
              wallCycleState, retainedWallCycleEncoding,
              &wallCyclePersistenceError, wallCyclePersistenceLimits)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::LimitExceeded
              && retainedWallCycleEncoding
                  == retainedWallCycleEncodingBefore,
          "regional opening wall-cycle persistence: encode byte limit retains destination");
    auto corruptEncodedWallCycleState = wallCycleState;
    corruptEncodedWallCycleState.wallTractions.tractions.front()
        .tractionPascals.x += 1.0;
    retainedWallCycleEncoding = retainedWallCycleEncodingBefore;
    check(!serializeWallCycleState(
              corruptEncodedWallCycleState, retainedWallCycleEncoding,
              &wallCyclePersistenceError,
              defaultWallCyclePersistenceLimits)
              && wallCyclePersistenceError.code
                  == WallCyclePersistenceErrorCode::InvalidData
              && retainedWallCycleEncoding
                  == retainedWallCycleEncodingBefore,
          "regional opening wall-cycle persistence: corrupt in-memory traction does not replace bytes");

    SceneFluidRegionalOpeningMomentumWallCycleOwner wallCycleOwner;
    check(wallCycleOwner.tryCommit(
              wallPressureEpoch, exchange, endpoint.openingMetric,
              endpoint.openingMetric, fixture.quadrature)
              && wallCycleOwner.hasState()
              && wallCycleOwner.checkpoint() == wallCycleState,
          "regional opening wall-cycle owner: accepted receipt commits atomically");
    SceneFluidRegionalOpeningMomentumWallCycleOwner restoredWallCycleOwner;
    restoredWallCycleOwner.restore(
        decodedWallCycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, fixture.quadrature);
    const auto restoredNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            restoredWallCycleOwner.state().adjustedMomentum,
            endpoint.openingMetric, adjustedTargetFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology,
            endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    check(restoredWallCycleOwner.state() == wallCycleState
              && restoredNextTransport == adjustedNextTransport,
          "regional opening wall-cycle owner: restored endpoint replays the exact next transport");

    SceneFixture wallLoadFixture(partialOpeningPressureScene(), false);
    const auto beforeWallLoad = wallLoadFixture.structure.checkpoint();
    const auto wallLoadApplication = applyOpeningMomentumWallLoads(
        decodedWallCycleState, endpoint, wallLoadFixture);
    validateOpeningMomentumWallLoads(
        wallLoadApplication, decodedWallCycleState, endpoint,
        wallLoadFixture);
    const auto afterWallLoad = wallLoadFixture.structure.checkpoint();
    const auto& actionReactionResidual =
        wallLoadApplication
            .actionReactionResidualKilogramMetersPerSecond;
    check(wallLoadApplication.applied
              && wallLoadApplication.sourceCycleStateFingerprint
                  == decodedWallCycleState.fingerprint
              && wallLoadApplication.sourceWallTractionFingerprint
                  == decodedWallCycleState.wallTractions.fingerprint
              && wallLoadApplication.nodeLoads.size()
                  == wallLoadFixture.transfer.nodes().size()
              && (std::abs(
                      wallLoadApplication.appliedWallForceNewtons.x)
                      > 0.0
                  || std::abs(
                      wallLoadApplication.appliedWallForceNewtons.y)
                      > 0.0
                  || std::abs(
                      wallLoadApplication.appliedWallForceNewtons.z)
                      > 0.0)
              && std::abs(actionReactionResidual.x) < 1.0e-12
              && std::abs(actionReactionResidual.y) < 1.0e-12
              && std::abs(actionReactionResidual.z) < 1.0e-12
              && afterWallLoad.pendingExternalForcesNewtons
                  != beforeWallLoad.pendingExternalForcesNewtons
              && afterWallLoad.nodes == beforeWallLoad.nodes
              && afterWallLoad.lastAppliedExternalForceNewtons
                  == beforeWallLoad.lastAppliedExternalForceNewtons,
          "regional opening wall-load application: restored SWRW traction closes action/reaction and changes only pending loads");

    SceneFixture repeatedWallLoadFixture(
        partialOpeningPressureScene(), false);
    const auto repeatedWallLoadApplication =
        applyOpeningMomentumWallLoads(
            decodedWallCycleState, endpoint, repeatedWallLoadFixture);
    check(repeatedWallLoadApplication == wallLoadApplication,
          "regional opening wall-load application: rebuilt target replay is deterministic");

    auto corruptWallLoadApplication = wallLoadApplication;
    corruptWallLoadApplication.nodeLoads.front()
        .appliedWallForceNewtons.x += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallLoadApplicationIntegrity(
                corruptWallLoadApplication);
        },
        "regional opening wall-load application: corrupt nodal force rejects");
    auto foreignWallLoadSettings =
        SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings{};
    foreignWallLoadSettings.transfer.momentReferenceMeters.x = 0.25;
    expectRejected(
        [&] {
            validateOpeningMomentumWallLoads(
                wallLoadApplication, decodedWallCycleState, endpoint,
                wallLoadFixture, foreignWallLoadSettings);
        },
        "regional opening wall-load application: foreign settings reject replay");

    SceneFixture limitedWallLoadFixture(
        partialOpeningPressureScene(), false);
    const auto beforeLimitedWallLoad =
        limitedWallLoadFixture.structure.checkpoint();
    auto wallLoadLimits =
        SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits{};
    wallLoadLimits.maximumNodeLoads =
        wallLoadApplication.nodeLoads.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningMomentumWallLoads(
                decodedWallCycleState, endpoint, limitedWallLoadFixture,
                {}, wallLoadLimits));
        },
        "regional opening wall-load application: node-load limit rejects before mutation");
    check(samePublicCheckpoint(
              beforeLimitedWallLoad,
              limitedWallLoadFixture.structure.checkpoint()),
          "regional opening wall-load application: limit rejection retains Structure checkpoint");

    const auto beforeForeignWallLoad =
        foreignWallCycleFixture.structure.checkpoint();
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningMomentumWallLoads(
                decodedWallCycleState, endpoint,
                foreignWallCycleFixture));
        },
        "regional opening wall-load application: foreign scene quadrature rejects before mutation");
    check(samePublicCheckpoint(
              beforeForeignWallLoad,
              foreignWallCycleFixture.structure.checkpoint()),
          "regional opening wall-load application: foreign-source rejection retains Structure checkpoint");

    SceneFixture combinedWallLoadFixture(
        partialOpeningPressureScene(), false);
    const auto beforeCombinedWallLoad =
        combinedWallLoadFixture.structure.checkpoint();
    const auto combinedWallLoadEpoch =
        applyOpeningMomentumWallLoadEpoch(
            decodedWallCycleState, endpoint, combinedWallLoadFixture);
    validateOpeningMomentumWallLoadEpoch(
        combinedWallLoadEpoch, decodedWallCycleState, endpoint,
        combinedWallLoadFixture);
    const auto afterCombinedWallLoad =
        combinedWallLoadFixture.structure.checkpoint();
    bool exactNodalHandoff =
        combinedWallLoadEpoch.pressureLoad.application.nodeLoads.size()
        == combinedWallLoadEpoch.wallLoad.nodeLoads.size();
    if (exactNodalHandoff) {
        for (std::size_t index = 0;
             index < combinedWallLoadEpoch.wallLoad.nodeLoads.size();
             ++index) {
            const auto& pressureNode =
                combinedWallLoadEpoch.pressureLoad.application
                    .nodeLoads[index];
            const auto& wallNode =
                combinedWallLoadEpoch.wallLoad.nodeLoads[index];
            exactNodalHandoff = exactNodalHandoff
                && pressureNode.stableId == wallNode.stableId
                && pressureNode.structureNode == wallNode.structureNode
                && pressureNode.resultingPendingForceNewtons
                    == wallNode.priorPendingForceNewtons;
        }
    }
    check(combinedWallLoadEpoch.applied
              && combinedWallLoadEpoch.sourceCycleStateFingerprint
                  == decodedWallCycleState.fingerprint
              && combinedWallLoadEpoch.pressureLoad.applied
              && combinedWallLoadEpoch.wallLoad.applied
              && exactNodalHandoff
              && combinedWallLoadEpoch
                     .pressureLoad.application
                     .resultingPendingForceNewtons
                  == combinedWallLoadEpoch
                         .wallLoad.priorPendingForceNewtons
              && combinedWallLoadEpoch.combinedAppliedForceNewtons
                  == StructureVector3{
                      combinedWallLoadEpoch.appliedPressureForceNewtons.x
                          + combinedWallLoadEpoch.appliedWallForceNewtons.x,
                      combinedWallLoadEpoch.appliedPressureForceNewtons.y
                          + combinedWallLoadEpoch.appliedWallForceNewtons.y,
                      combinedWallLoadEpoch.appliedPressureForceNewtons.z
                          + combinedWallLoadEpoch.appliedWallForceNewtons.z}
              && afterCombinedWallLoad.pendingExternalForcesNewtons
                  != beforeCombinedWallLoad.pendingExternalForcesNewtons
              && afterCombinedWallLoad.nodes
                  == beforeCombinedWallLoad.nodes
              && afterCombinedWallLoad.lastAppliedExternalForceNewtons
                  == beforeCombinedWallLoad
                         .lastAppliedExternalForceNewtons,
          "regional opening momentum wall-load epoch: pressure and wall loads hand off per node under one pending-load transaction");

    SceneFixture repeatedCombinedWallLoadFixture(
        partialOpeningPressureScene(), false);
    const auto repeatedCombinedWallLoadEpoch =
        applyOpeningMomentumWallLoadEpoch(
            decodedWallCycleState, endpoint,
            repeatedCombinedWallLoadFixture);
    check(repeatedCombinedWallLoadEpoch == combinedWallLoadEpoch,
          "regional opening momentum wall-load epoch: restored-state rebuilt-target replay is deterministic");

    auto corruptCombinedWallLoadEpoch = combinedWallLoadEpoch;
    corruptCombinedWallLoadEpoch.wallLoad.nodeLoads.front()
        .priorPendingForceNewtons.x += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallLoadEpochIntegrity(
                corruptCombinedWallLoadEpoch);
        },
        "regional opening momentum wall-load epoch: corrupt nodal handoff rejects");
    auto foreignCombinedSettings =
        SceneFluidRegionalOpeningMomentumWallLoadEpochSettings{};
    foreignCombinedSettings.transfer.momentReferenceMeters.y = 0.5;
    expectRejected(
        [&] {
            validateOpeningMomentumWallLoadEpoch(
                combinedWallLoadEpoch, decodedWallCycleState, endpoint,
                combinedWallLoadFixture, foreignCombinedSettings);
        },
        "regional opening momentum wall-load epoch: foreign settings reject replay");

    SceneFixture lateLimitedCombinedFixture(
        partialOpeningPressureScene(), false);
    const auto beforeLateLimitedCombined =
        lateLimitedCombinedFixture.structure.checkpoint();
    auto combinedWallLoadLimits =
        SceneFluidRegionalOpeningMomentumWallLoadEpochLimits{};
    combinedWallLoadLimits.maximumOwnedBytes =
        combinedWallLoadEpoch.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(applyOpeningMomentumWallLoadEpoch(
                decodedWallCycleState, endpoint,
                lateLimitedCombinedFixture, {},
                combinedWallLoadLimits));
        },
        "regional opening momentum wall-load epoch: late aggregate limit rejects after nested applications");
    check(samePublicCheckpoint(
              beforeLateLimitedCombined,
              lateLimitedCombinedFixture.structure.checkpoint()),
          "regional opening momentum wall-load epoch: late outer rejection restores both pressure and wall loads");

    SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings
        structureStepSettings;
    structureStepSettings.structure.timeStepSeconds =
        decodedWallCycleState.adjustedMomentum.timeStepSeconds;
    structureStepSettings.structure.substeps = 1;
    structureStepSettings.structure.constraintIterations = 8;
    structureStepSettings.structure.gravityMetersPerSecondSquared = {};
    structureStepSettings.structure.velocityDampingPerSecond = 0.0;
    structureStepSettings.structure.workerThreads = 0;

    SceneFixture structureStepFixture(
        partialOpeningPressureScene(), false);
    const auto beforeStructureStep =
        structureStepFixture.structure.checkpoint();
    const auto structureStepEpoch =
        advanceOpeningMomentumWallStructureStepEpoch(
            decodedWallCycleState, endpoint, structureStepFixture,
            structureStepSettings);
    const auto afterStructureStep =
        structureStepFixture.structure.checkpoint();
    validateOpeningMomentumWallStructureStepEpoch(
        structureStepEpoch, decodedWallCycleState, endpoint,
        structureStepFixture, structureStepSettings);
    const bool structureMoved =
        !std::ranges::equal(
            beforeStructureStep.nodes, afterStructureStep.nodes);
    const bool pendingLoadsConsumed = std::ranges::all_of(
        afterStructureStep.pendingExternalForcesNewtons,
        [](const StructureVector3& value) {
            return value == StructureVector3{};
        });
    check(structureStepEpoch.stepped
              && structureStepEpoch.version
                  == sceneFluidRegionalOpeningMomentumWallStructureStepEpochVersion
              && structureStepEpoch.fingerprint != 0
              && structureStepEpoch.loadEpoch.applied
              && structureStepEpoch.diagnostics.finite
              && structureStepEpoch.beforeStructureCheckpoint.size() > 24
              && structureStepEpoch.afterStructureCheckpoint.size() > 24
              && structureStepEpoch.beforeStructureCheckpoint
                  != structureStepEpoch.afterStructureCheckpoint
              && structureMoved && pendingLoadsConsumed
              && afterStructureStep.acceptedStepCount
                  == beforeStructureStep.acceptedStepCount + 1
              && afterStructureStep.simulationTimeSeconds
                  == beforeStructureStep.simulationTimeSeconds
                      + structureStepSettings.structure.timeStepSeconds
              && afterStructureStep.lastAppliedExternalForceNewtons
                  == structureStepEpoch.loadEpoch
                         .resultingPendingForceNewtons
              && structureStepEpoch.diagnostics
                     .lastAppliedExternalForceNewtons
                  == afterStructureStep.lastAppliedExternalForceNewtons,
          "regional opening wall structural-step epoch: combined loads deform XPBD once and retain complete endpoints");

    SceneFixture repeatedStructureStepFixture(
        partialOpeningPressureScene(), false);
    const auto repeatedStructureStepEpoch =
        advanceOpeningMomentumWallStructureStepEpoch(
            decodedWallCycleState, endpoint,
            repeatedStructureStepFixture, structureStepSettings);
    check(repeatedStructureStepEpoch.fingerprint
                  == structureStepEpoch.fingerprint
              && repeatedStructureStepEpoch.loadEpoch
                  == structureStepEpoch.loadEpoch
              && repeatedStructureStepEpoch.diagnostics
                  == structureStepEpoch.diagnostics
              && repeatedStructureStepEpoch.beforeStructureCheckpoint
                  == structureStepEpoch.beforeStructureCheckpoint
              && repeatedStructureStepEpoch.afterStructureCheckpoint
                  == structureStepEpoch.afterStructureCheckpoint,
          "regional opening wall structural-step epoch: rebuilt Structure replay is byte-identical");

    auto corruptStructureStepEpoch = structureStepEpoch;
    corruptStructureStepEpoch.afterStructureCheckpoint.back() ^= 0x01U;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
                corruptStructureStepEpoch);
        },
        "regional opening wall structural-step epoch: corrupt retained checkpoint rejects");
    corruptStructureStepEpoch = structureStepEpoch;
    corruptStructureStepEpoch.diagnostics.kineticEnergyJoules += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
                corruptStructureStepEpoch);
        },
        "regional opening wall structural-step epoch: corrupt diagnostics reject");

    auto foreignStructureStepSettings = structureStepSettings;
    foreignStructureStepSettings.structure.constraintIterations += 1;
    expectRejected(
        [&] {
            validateOpeningMomentumWallStructureStepEpoch(
                structureStepEpoch, decodedWallCycleState, endpoint,
                structureStepFixture, foreignStructureStepSettings);
        },
        "regional opening wall structural-step epoch: foreign structural settings reject replay");

    SceneFixture mismatchedStepFixture(
        partialOpeningPressureScene(), false);
    const auto beforeMismatchedStep =
        mismatchedStepFixture.structure.checkpoint();
    auto mismatchedStepSettings = structureStepSettings;
    mismatchedStepSettings.structure.timeStepSeconds *= 0.5;
    expectRejected(
        [&] {
            static_cast<void>(
                advanceOpeningMomentumWallStructureStepEpoch(
                    decodedWallCycleState, endpoint,
                    mismatchedStepFixture, mismatchedStepSettings));
        },
        "regional opening wall structural-step epoch: mismatched fluid and structural time steps reject");
    check(samePublicCheckpoint(
              beforeMismatchedStep,
              mismatchedStepFixture.structure.checkpoint()),
          "regional opening wall structural-step epoch: time-step rejection preserves target state");

    SceneFixture lateLimitedStructureStepFixture(
        partialOpeningPressureScene(), false);
    const auto beforeLateLimitedStructureStep =
        lateLimitedStructureStepFixture.structure.checkpoint();
    auto structureStepLimits =
        SceneFluidRegionalOpeningMomentumWallStructureStepEpochLimits{};
    structureStepLimits.maximumOwnedBytes =
        structureStepEpoch.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                advanceOpeningMomentumWallStructureStepEpoch(
                    decodedWallCycleState, endpoint,
                    lateLimitedStructureStepFixture,
                    structureStepSettings, structureStepLimits));
        },
        "regional opening wall structural-step epoch: late aggregate limit rejects after XPBD consumption");
    check(samePublicCheckpoint(
              beforeLateLimitedStructureStep,
              lateLimitedStructureStepFixture.structure.checkpoint()),
          "regional opening wall structural-step epoch: late outer rejection restores pre-load Structure state");

    SceneFixture coupledOwnerSources(
        partialOpeningPressureScene(), false);
    const auto advanceCoupledOwner = [&](
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& owner,
        const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits&
            limits) {
        owner.advance(
            decodedWallCycleState, endpoint.openingMetric,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric, coupledOwnerSources.surface.definition,
            coupledOwnerSources.state, coupledOwnerSources.transfer,
            coupledOwnerSources.quadrature, structureStepSettings, limits);
    };
    const auto restoreCoupledOwner = [&](
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& owner,
        const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
        const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits&
            limits) {
        owner.restore(
            state, endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.state, coupledOwnerSources.transfer,
            coupledOwnerSources.quadrature, structureStepSettings, limits);
    };
    const auto encodeOwnerStructure = [](
        const SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& owner) {
        std::vector<std::uint8_t> bytes;
        StructureCheckpointPersistenceError error;
        if (!serializeStructureCheckpoint(
                owner.structure(), owner.structure().checkpoint(), bytes,
                &error)) {
            throw std::runtime_error(
                "cannot encode coupled-owner Structure: " + error.message);
        }
        return bytes;
    };

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner coupledOwner(
        Structure(coupledOwnerSources.structureAssembly.definition));
    advanceCoupledOwner(coupledOwner, {});
    const auto coupledCheckpoint = coupledOwner.checkpoint();
    const auto coupledStructureBytes = encodeOwnerStructure(coupledOwner);
    check(coupledOwner.hasState()
              && coupledCheckpoint.version
                  == sceneFluidRegionalOpeningMomentumWallCoupledStateVersion
              && coupledCheckpoint.fingerprint != 0
              && coupledCheckpoint.cycleState
                  == decodedWallCycleState
              && coupledCheckpoint.structureStep.fingerprint
                  == structureStepEpoch.fingerprint
              && coupledCheckpoint.structureStep.afterStructureCheckpoint
                  == coupledStructureBytes,
          "regional opening wall coupled owner: fluid and post-XPBD endpoints publish atomically");

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner
        restoredCoupledOwner(
            Structure(coupledOwnerSources.structureAssembly.definition));
    restoreCoupledOwner(restoredCoupledOwner, coupledCheckpoint, {});
    const auto restoredCoupledStructureBytes =
        encodeOwnerStructure(restoredCoupledOwner);
    const auto restoredCoupledNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            restoredCoupledOwner.state().cycleState.adjustedMomentum,
            endpoint.openingMetric, adjustedTargetFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology,
            endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    check(restoredCoupledOwner.state().fingerprint
                  == coupledCheckpoint.fingerprint
              && restoredCoupledStructureBytes
                  == coupledStructureBytes
              && restoredCoupledNextTransport == adjustedNextTransport,
          "regional opening wall coupled owner: restored pair replays exact Structure and next fluid transport");

    const auto beforeRejectedCoupledRestore =
        encodeOwnerStructure(restoredCoupledOwner);
    const auto priorCoupledFingerprint =
        restoredCoupledOwner.state().fingerprint;
    auto corruptCoupledCheckpoint = coupledCheckpoint;
    corruptCoupledCheckpoint.structureStep.diagnostics
        .kineticEnergyJoules += 1.0;
    expectRejected(
        [&] {
            restoreCoupledOwner(
                restoredCoupledOwner, corruptCoupledCheckpoint, {});
        },
        "regional opening wall coupled owner: corrupt structural endpoint rejects restore");
    check(restoredCoupledOwner.state().fingerprint
                  == priorCoupledFingerprint
              && encodeOwnerStructure(restoredCoupledOwner)
                  == beforeRejectedCoupledRestore,
          "regional opening wall coupled owner: rejected restore preserves both live owners");

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner
        lateLimitedCoupledOwner(
            Structure(coupledOwnerSources.structureAssembly.definition));
    const auto beforeLateLimitedCoupledOwner =
        lateLimitedCoupledOwner.structure().checkpoint();
    auto coupledOwnerLimits =
        SceneFluidRegionalOpeningMomentumWallCoupledStateLimits{};
    coupledOwnerLimits.maximumOwnedBytes =
        coupledCheckpoint.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            advanceCoupledOwner(
                lateLimitedCoupledOwner, coupledOwnerLimits);
        },
        "regional opening wall coupled owner: late owner limit rejects after structural acceptance");
    check(!lateLimitedCoupledOwner.hasState()
              && samePublicCheckpoint(
                  beforeLateLimitedCoupledOwner,
                  lateLimitedCoupledOwner.structure().checkpoint()),
          "regional opening wall coupled owner: late rejection retains prior fluid absence and exact Structure");

    const auto serializeCoupledState = [&](
        const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
        const Structure& definitionOwner,
        std::vector<std::uint8_t>& bytes,
        SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError*
            error,
        const SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits&
            limits) {
        return serializeSceneFluidRegionalOpeningMomentumWallCoupledState(
            state, endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.state, coupledOwnerSources.transfer,
            coupledOwnerSources.quadrature, definitionOwner,
            structureStepSettings, bytes, error, limits);
    };
    const auto deserializeCoupledState = [&](
        const std::span<const std::uint8_t> bytes,
        const Structure& definitionOwner,
        const SceneFluidQuadratureDefinition& sourceQuadrature,
        SceneFluidRegionalOpeningMomentumWallCoupledState& state,
        SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError*
            error,
        const SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits&
            limits) {
        return deserializeSceneFluidRegionalOpeningMomentumWallCoupledState(
            bytes, endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.state, coupledOwnerSources.transfer,
            sourceQuadrature, definitionOwner, structureStepSettings,
            state, error, limits);
    };

    std::vector<std::uint8_t> coupledEncoding;
    SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceError
        coupledPersistenceError;
    check(serializeCoupledState(
              coupledCheckpoint, coupledOwner.structure(), coupledEncoding,
              &coupledPersistenceError, {})
              && coupledEncoding.size() > 64
              && coupledEncoding[0] == 'S'
              && coupledEncoding[1] == 'W'
              && coupledEncoding[2] == 'R'
              && coupledEncoding[3] == 'C',
          "regional opening wall coupled persistence: accepted pair encodes as distinct SWRC");
    SceneFluidRegionalOpeningMomentumWallCoupledState decodedCoupledState;
    check(deserializeCoupledState(
              coupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature, decodedCoupledState,
              &coupledPersistenceError, {})
              && decodedCoupledState.fingerprint
                  == coupledCheckpoint.fingerprint
              && decodedCoupledState.cycleState
                  == coupledCheckpoint.cycleState
              && decodedCoupledState.structureStep.fingerprint
                  == coupledCheckpoint.structureStep.fingerprint
              && decodedCoupledState.structureStep
                     .beforeStructureCheckpoint
                  == coupledCheckpoint.structureStep
                         .beforeStructureCheckpoint
              && decodedCoupledState.structureStep
                     .afterStructureCheckpoint
                  == coupledCheckpoint.structureStep
                         .afterStructureCheckpoint,
          "regional opening wall coupled persistence: decode replays the complete in-memory receipt");
    std::vector<std::uint8_t> repeatedCoupledEncoding;
    check(serializeCoupledState(
              decodedCoupledState, coupledOwnerSources.structure,
              repeatedCoupledEncoding, &coupledPersistenceError, {})
              && repeatedCoupledEncoding == coupledEncoding,
          "regional opening wall coupled persistence: round trip is bit-exact");

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner
        decodedCoupledOwner(
            Structure(coupledOwnerSources.structureAssembly.definition));
    restoreCoupledOwner(decodedCoupledOwner, decodedCoupledState, {});
    const auto decodedNextTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            decodedCoupledOwner.state().cycleState.adjustedMomentum,
            endpoint.openingMetric, adjustedTargetFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology,
            endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    check(encodeOwnerStructure(decodedCoupledOwner)
                  == coupledStructureBytes
              && decodedNextTransport == adjustedNextTransport,
          "regional opening wall coupled persistence: restored SWRC pair preserves both continuations");

    auto retainedDecodedCoupledState = decodedCoupledState;
    const auto decodedDestinationRetained = [&] {
        return retainedDecodedCoupledState.fingerprint
                == decodedCoupledState.fingerprint
            && retainedDecodedCoupledState.structureStep
                   .afterStructureCheckpoint
                == decodedCoupledState.structureStep
                       .afterStructureCheckpoint;
    };
    auto corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding[0] ^= 0x01U;
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::InvalidMagic
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: foreign magic retains destination");
    corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding.back() ^= 0x01U;
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::ChecksumMismatch
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: checksum corruption retains destination");
    corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding[32] ^= 0x01U;
    refreshOpeningMomentumWallCoupledStatePersistenceChecksum(
        corruptCoupledEncoding);
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::ReplayMismatch
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: recomputed-checksum receipt identity corruption rejects");
    corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding[64] ^= 0x01U;
    refreshOpeningMomentumWallCoupledStatePersistenceChecksum(
        corruptCoupledEncoding);
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: recomputed-checksum nested SWRW corruption rejects");
    corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding.pop_back();
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::Truncated
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: truncation retains destination");
    corruptCoupledEncoding = coupledEncoding;
    corruptCoupledEncoding.push_back(0);
    check(!deserializeCoupledState(
              corruptCoupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::TrailingData
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: trailing data retains destination");

    check(!deserializeCoupledState(
              coupledEncoding, coupledOwnerSources.structure,
              foreignWallCycleFixture.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::SourceMismatch
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: foreign quadrature retains destination");
    auto foreignStructureDefinition =
        coupledOwnerSources.structureAssembly.definition;
    foreignStructureDefinition.nodes.front().massKg *= 2.0;
    Structure foreignStructure(std::move(foreignStructureDefinition));
    check(!deserializeCoupledState(
              coupledEncoding, foreignStructure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError, {})
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::TopologyMismatch
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: foreign Structure definition retains destination");

    auto coupledPersistenceLimits =
        SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceLimits{};
    coupledPersistenceLimits.maximumEncodedBytes =
        coupledEncoding.size() - 1;
    check(!deserializeCoupledState(
              coupledEncoding, coupledOwnerSources.structure,
              coupledOwnerSources.quadrature,
              retainedDecodedCoupledState, &coupledPersistenceError,
              coupledPersistenceLimits)
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::LimitExceeded
              && decodedDestinationRetained(),
          "regional opening wall coupled persistence: decode byte limit retains destination");
    std::vector<std::uint8_t> retainedCoupledEncoding{1, 2, 3};
    check(!serializeCoupledState(
              coupledCheckpoint, coupledOwnerSources.structure,
              retainedCoupledEncoding, &coupledPersistenceError,
              coupledPersistenceLimits)
              && coupledPersistenceError.code
                  == SceneFluidRegionalOpeningMomentumWallCoupledStatePersistenceErrorCode::LimitExceeded
              && retainedCoupledEncoding
                  == std::vector<std::uint8_t>({1, 2, 3}),
          "regional opening wall coupled persistence: encode byte limit retains destination");

    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    coupledCheckpoint,
                    coupledOwnerSources.surface.definition,
                    coupledOwnerSources.structureAssembly.mappings, grid(),
                    coupledOwnerSources.transfer,
                    coupledOwner.structure()));
        },
        "regional opening wall post-step geometry: out-of-domain structural endpoint rejects without changing grids");

    auto fixedStructureDefinition =
        coupledOwnerSources.structureAssembly.definition;
    for (auto& node : fixedStructureDefinition.nodes) {
        node.fixed = true;
    }
    Structure fixedSourceStructure(fixedStructureDefinition);
    SceneFluidSurfaceTransfer fixedTransfer(
        coupledOwnerSources.surface.definition,
        coupledOwnerSources.structureAssembly.mappings,
        fixedSourceStructure);
    const auto fixedSourceState = captureSceneFluidSurfaceState(
        coupledOwnerSources.surface.definition,
        coupledOwnerSources.structureAssembly.mappings,
        fixedSourceStructure);
    const auto fixedSourceGridEpoch = buildSceneFluidGridEpoch(
        coupledOwnerSources.surface.definition, fixedSourceState, grid(),
        fixedTransfer);
    const auto fixedWallInput =
        captureSceneFluidRegionalOpeningMomentumWallInput(
            endpoint.momentumCycle.transport, endpoint.acceptedFlow,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric,
            coupledOwnerSources.surface.definition, fixedSourceState,
            fixedSourceGridEpoch.quadrature);
    const auto fixedWallExchange =
        exchangeSceneFluidRegionalOpeningMomentumWall(
            fixedWallInput, settings);
    const auto fixedWallPressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            fixedWallExchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric, {},
            endpoint.momentumCycle.pressureSettings);
    const auto fixedWallCycleState =
        captureSceneFluidRegionalOpeningMomentumWallCycleState(
            fixedWallPressureEpoch, fixedWallExchange,
            endpoint.openingMetric, endpoint.openingMetric,
            fixedSourceGridEpoch.quadrature);
    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner fixedCoupledOwner{
        Structure(fixedStructureDefinition)};
    fixedCoupledOwner.advance(
        fixedWallCycleState, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, coupledOwnerSources.surface.definition,
        fixedSourceState, fixedTransfer, fixedSourceGridEpoch.quadrature,
        structureStepSettings);
    const auto fixedCoupledCheckpoint = fixedCoupledOwner.checkpoint();
    const auto postStepGeometry =
        buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
            fixedCoupledCheckpoint,
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.structureAssembly.mappings, grid(),
            fixedTransfer, fixedCoupledOwner.structure());
    validateSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
        postStepGeometry, fixedCoupledCheckpoint,
        coupledOwnerSources.surface.definition,
        coupledOwnerSources.structureAssembly.mappings, grid(),
        fixedTransfer, fixedCoupledOwner.structure());
    const PeriodicCartesianGrid foreignPostStepGrid(
        {8, 2, 2}, {-2.0, -1.0, -1.0}, {2.0, 1.0, 1.0});
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                postStepGeometry, fixedCoupledCheckpoint,
                coupledOwnerSources.surface.definition,
                coupledOwnerSources.structureAssembly.mappings,
                foreignPostStepGrid, fixedTransfer,
                fixedCoupledOwner.structure());
        },
        "regional opening wall post-step geometry: foreign CFD grid rejects receipt replay");
    check(postStepGeometry.version
                  == sceneFluidRegionalOpeningMomentumWallPostStepGeometryVersion
              && postStepGeometry.fingerprint != 0
              && postStepGeometry.sourceCoupledStateFingerprint
                  == fixedCoupledCheckpoint.fingerprint
              && postStepGeometry.sourceSurfaceStateFingerprint
                  == fixedSourceState.fingerprint
              && postStepGeometry.previousAcceptedStepCount == 0
              && postStepGeometry.currentAcceptedStepCount == 1
              && postStepGeometry.surfaceState.acceptedStepCount == 1
              && postStepGeometry.surfaceState.simulationTimeSeconds
                  == structureStepSettings.structure.timeStepSeconds
              && postStepGeometry.surfaceState
                  != fixedSourceState
              && postStepGeometry.gridEpoch.fingerprint != 0
              && postStepGeometry.gridEpoch.quadrature.fingerprint != 0
              && postStepGeometry.ownedStorageBytes
                  > postStepGeometry.gridEpoch.ownedStorageBytes,
          "regional opening wall post-step geometry: accepted Structure publishes the next authoritative surface/grid epoch");

    SceneFluidRegionalOpeningMomentumWallCoupledStateOwner
        restoredFixedCoupledOwner{Structure(fixedStructureDefinition)};
    restoredFixedCoupledOwner.restore(
        fixedCoupledCheckpoint, endpoint.openingMetric,
        endpoint.pressureOperator, endpoint.basePressureOperator,
        endpoint.geometry, endpoint.sweep, endpoint.fragments,
        endpoint.topology, endpoint.volumeRates,
        endpoint.openingDefinitions, endpoint.openings,
        endpoint.resistanceDefinitions, endpoint.baseMetric,
        endpoint.openingMetric, coupledOwnerSources.surface.definition,
        fixedSourceState, fixedTransfer, fixedSourceGridEpoch.quadrature,
        structureStepSettings);
    const auto restoredPostStepGeometry =
        buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
            restoredFixedCoupledOwner.state(),
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.structureAssembly.mappings, grid(),
            fixedTransfer, restoredFixedCoupledOwner.structure());
    check(restoredPostStepGeometry == postStepGeometry,
          "regional opening wall post-step geometry: restored coupled owner rebuild is deterministic");

    auto corruptPostStepGeometry = postStepGeometry;
    corruptPostStepGeometry.surfaceState.vertices.front()
        .positionMeters.x += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallPostStepGeometryIntegrity(
                corruptPostStepGeometry);
        },
        "regional opening wall post-step geometry: corrupt surface state rejects");
    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    fixedCoupledCheckpoint,
                    coupledOwnerSources.surface.definition,
                    coupledOwnerSources.structureAssembly.mappings, grid(),
                    fixedTransfer, fixedSourceStructure));
        },
        "regional opening wall post-step geometry: stale pre-step Structure rejects");
    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    fixedCoupledCheckpoint,
                    foreignWallCycleFixture.surface.definition,
                    foreignWallCycleFixture.structureAssembly.mappings,
                    grid(), foreignWallCycleFixture.transfer,
                    fixedCoupledOwner.structure()));
        },
        "regional opening wall post-step geometry: foreign source surface rejects");

    auto postStepGeometryLimits =
        SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits{};
    postStepGeometryLimits.maximumSurfaceVertices =
        postStepGeometry.surfaceState.vertices.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    fixedCoupledCheckpoint,
                    coupledOwnerSources.surface.definition,
                    coupledOwnerSources.structureAssembly.mappings, grid(),
                    fixedTransfer, fixedCoupledOwner.structure(),
                    postStepGeometryLimits));
        },
        "regional opening wall post-step geometry: surface vertex limit rejects");
    postStepGeometryLimits = {};
    postStepGeometryLimits.structureCheckpoint.maximumEncodedBytes =
        fixedCoupledCheckpoint.structureStep
            .afterStructureCheckpoint.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    fixedCoupledCheckpoint,
                    coupledOwnerSources.surface.definition,
                    coupledOwnerSources.structureAssembly.mappings, grid(),
                    fixedTransfer, fixedCoupledOwner.structure(),
                    postStepGeometryLimits));
        },
        "regional opening wall post-step geometry: retained Structure byte limit rejects");
    postStepGeometryLimits = {};
    postStepGeometryLimits.maximumOwnedBytes =
        postStepGeometry.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
                    fixedCoupledCheckpoint,
                    coupledOwnerSources.surface.definition,
                    coupledOwnerSources.structureAssembly.mappings, grid(),
                    fixedTransfer, fixedCoupledOwner.structure(),
                    postStepGeometryLimits));
        },
        "regional opening wall post-step geometry: late aggregate limit rejects complete candidate");

    const auto fixedAcceptedFlow =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            fixedWallCycleState.acceptedPressure,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric);
    const auto fixedConsecutiveTransport =
        advancePlanarPressureRegionFragmentOpeningMomentum(
            fixedWallCycleState.adjustedMomentum,
            endpoint.openingMetric, fixedAcceptedFlow,
            endpoint.openingMetric, endpoint.geometry, endpoint.sweep,
            endpoint.fragments, endpoint.topology, endpoint.volumeRates,
            endpoint.momentumCycle.transportSettings);
    const auto fixedConsecutiveWallInput =
        captureSceneFluidRegionalOpeningMomentumWallInput(
            fixedConsecutiveTransport,
            fixedWallCycleState.acceptedPressure,
            endpoint.pressureOperator, endpoint.basePressureOperator,
            endpoint.geometry, endpoint.sweep, endpoint.fragments,
            endpoint.topology, endpoint.volumeRates,
            endpoint.openingDefinitions, endpoint.openings,
            endpoint.resistanceDefinitions, endpoint.baseMetric,
            endpoint.openingMetric,
            coupledOwnerSources.surface.definition,
            postStepGeometry.surfaceState,
            postStepGeometry.gridEpoch.quadrature);
    const auto fixedConsecutiveWallExchange =
        exchangeSceneFluidRegionalOpeningMomentumWall(
            fixedConsecutiveWallInput, settings);
    const auto fixedConsecutivePressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            fixedConsecutiveWallExchange, fixedConsecutiveTransport,
            endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric, {},
            endpoint.momentumCycle.pressureSettings);
    const auto fixedConsecutiveCycleState =
        captureSceneFluidRegionalOpeningMomentumWallCycleState(
            fixedConsecutivePressureEpoch,
            fixedConsecutiveWallExchange, endpoint.openingMetric,
            endpoint.openingMetric,
            postStepGeometry.gridEpoch.quadrature);
    const auto beforeRejectedConsecutiveState =
        fixedCoupledOwner.checkpoint();
    const auto beforeRejectedConsecutiveStructure =
        fixedCoupledOwner.structure().checkpoint();
    expectRejected(
        [&] {
            fixedCoupledOwner.advanceFixedMetricConsecutive(
                fixedConsecutiveCycleState, adjustedNextTransport,
                endpoint.openingMetric, endpoint.pressureOperator,
                endpoint.basePressureOperator, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates, endpoint.openingDefinitions,
                endpoint.openings, endpoint.resistanceDefinitions,
                endpoint.baseMetric, endpoint.openingMetric,
                coupledOwnerSources.surface.definition,
                postStepGeometry.surfaceState, fixedTransfer,
                postStepGeometry.gridEpoch.quadrature,
                structureStepSettings);
        },
        "regional opening wall consecutive owner: foreign transport lineage rejects");
    check(fixedCoupledOwner.state().fingerprint
                  == beforeRejectedConsecutiveState.fingerprint
              && samePublicCheckpoint(
                  beforeRejectedConsecutiveStructure,
                  fixedCoupledOwner.structure().checkpoint()),
          "regional opening wall consecutive owner: rejected lineage preserves both accepted owners");

    const auto advanceFixedConsecutive = [&](
        SceneFluidRegionalOpeningMomentumWallCoupledStateOwner& owner) {
        owner.advanceFixedMetricConsecutive(
            fixedConsecutiveCycleState, fixedConsecutiveTransport,
            endpoint.openingMetric, endpoint.pressureOperator,
            endpoint.basePressureOperator, endpoint.geometry,
            endpoint.sweep, endpoint.fragments, endpoint.topology,
            endpoint.volumeRates, endpoint.openingDefinitions,
            endpoint.openings, endpoint.resistanceDefinitions,
            endpoint.baseMetric, endpoint.openingMetric,
            coupledOwnerSources.surface.definition,
            postStepGeometry.surfaceState, fixedTransfer,
            postStepGeometry.gridEpoch.quadrature,
            structureStepSettings);
    };
    advanceFixedConsecutive(fixedCoupledOwner);
    advanceFixedConsecutive(restoredFixedCoupledOwner);
    const auto fixedConsecutiveCheckpoint =
        fixedCoupledOwner.checkpoint();
    const auto restoredFixedConsecutiveCheckpoint =
        restoredFixedCoupledOwner.checkpoint();
    const auto secondPostStepGeometry =
        buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
            fixedConsecutiveCheckpoint,
            coupledOwnerSources.surface.definition,
            coupledOwnerSources.structureAssembly.mappings, grid(),
            fixedTransfer, fixedCoupledOwner.structure());
    check(fixedConsecutiveCheckpoint.fingerprint
                  == restoredFixedConsecutiveCheckpoint.fingerprint
              && fixedConsecutiveCheckpoint.cycleState
                  == fixedConsecutiveCycleState
              && fixedConsecutiveCheckpoint.structureStep
                     .beforeAcceptedStepCount == 1
              && fixedConsecutiveCheckpoint.structureStep
                     .afterAcceptedStepCount == 2
              && secondPostStepGeometry.previousAcceptedStepCount == 1
              && secondPostStepGeometry.currentAcceptedStepCount == 2
              && secondPostStepGeometry.surfaceState.acceptedStepCount == 2,
          "regional opening wall consecutive owner: fixed-metric fluid lineage and XPBD advance repeat deterministically");

    const auto rejectingLayers = translatePlanarPressureJumpLayers(
        endpoint.geometry, endpoint.layers, 0.1).layers;
    const auto rejectingSweep = makePlanarPressureRegionSweepLedger(
        endpoint.geometry, endpoint.layers, rejectingLayers, 1.0);
    const auto rejectingFragments = buildPlanarPressureRegionFragments(
        endpoint.geometry, rejectingSweep);
    const auto rejectingTopology =
        buildPlanarPressureRegionFragmentTopology(
            endpoint.geometry, rejectingSweep, rejectingFragments);
    const auto rejectingVolumeRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            endpoint.geometry, rejectingSweep, rejectingFragments,
            rejectingTopology);
    auto rejectingOpeningDefinitions = endpoint.openingDefinitions;
    for (auto& definition : rejectingOpeningDefinitions) {
        if (definition.authoredWrappedCentroidMeters) {
            definition.authoredWrappedCentroidMeters->x += 0.1;
        }
    }
    const auto rejectingOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            endpoint.geometry, rejectingSweep, rejectingFragments,
            rejectingTopology, rejectingOpeningDefinitions);
    const auto rejectingBasePressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            endpoint.geometry, rejectingSweep, rejectingFragments,
            rejectingTopology);
    const auto rejectingPressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            rejectingBasePressureOperator, endpoint.geometry,
            rejectingSweep, rejectingFragments, rejectingTopology,
            rejectingOpeningDefinitions, rejectingOpenings);
    const auto rejectingBaseMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            endpoint.geometry, rejectingSweep, rejectingFragments,
            rejectingTopology);
    const auto rejectingMetric =
        buildPlanarPressureRegionFragmentOpeningVelocityMetric(
            endpoint.geometry, rejectingSweep, rejectingFragments,
            rejectingTopology, rejectingBaseMetric,
            rejectingOpeningDefinitions, rejectingOpenings);
    auto rejectingPressureSettings =
        endpoint.momentumCycle.pressureSettings;
    rejectingPressureSettings.projection.timeStepSeconds = 1.0;
    rejectingPressureSettings.projection.pressureSolve.maximumIterations = 1;
    rejectingPressureSettings.projection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-16;
    rejectingPressureSettings.projection.pressureSolve
        .relativeResidualTolerance = 0.0;
    const auto rejectedWallPressureEpoch =
        acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
            exchange, endpoint.momentumCycle.transport,
            endpoint.openingMetric, rejectingPressureOperator,
            rejectingBasePressureOperator, endpoint.geometry,
            rejectingSweep, rejectingFragments, rejectingTopology,
            rejectingVolumeRates, rejectingOpeningDefinitions,
            rejectingOpenings, endpoint.resistanceDefinitions,
            rejectingBaseMetric, rejectingMetric, {},
            rejectingPressureSettings);
    const auto beforeRejectedCommit = wallCycleOwner.checkpoint();
    const bool rejectedCommit = wallCycleOwner.tryCommit(
        rejectedWallPressureEpoch, exchange, endpoint.openingMetric,
        rejectingMetric, fixture.quadrature);
    check(!rejectedWallPressureEpoch.pressureEpoch.diagnostics.accepted,
          "regional opening wall-cycle owner: negative pressure oracle rejects");
    check(!rejectedCommit,
          "regional opening wall-cycle owner: rejected pressure does not commit");
    check(wallCycleOwner.checkpoint() == beforeRejectedCommit,
          "regional opening wall-cycle owner: rejected pressure retains the prior endpoint exactly");

    auto corruptWallCycleState = wallCycleState;
    corruptWallCycleState.wallTractions.tractions.front()
        .tractionPascals.y += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(
                corruptWallCycleState);
        },
        "regional opening wall-cycle state: corrupt retained traction rejects");
    auto wallCycleLimits =
        SceneFluidRegionalOpeningMomentumWallCycleStateLimits{};
    wallCycleLimits.maximumWallTractions =
        wallCycleState.wallTractions.tractions.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                captureSceneFluidRegionalOpeningMomentumWallCycleState(
                    wallPressureEpoch, exchange, endpoint.openingMetric,
                    endpoint.openingMetric, fixture.quadrature,
                    wallCycleLimits));
        },
        "regional opening wall-cycle state: traction limit rejects before publication");
    const auto beforeForeignRestore = restoredWallCycleOwner.checkpoint();
    expectRejected(
        [&] {
            restoredWallCycleOwner.restore(
                wallCycleState, endpoint.openingMetric,
                endpoint.pressureOperator,
                endpoint.basePressureOperator, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates, endpoint.openingDefinitions,
                endpoint.openings, endpoint.resistanceDefinitions,
                endpoint.baseMetric, endpoint.openingMetric,
                foreignWallCycleFixture.quadrature);
        },
        "regional opening wall-cycle owner: foreign quadrature restore rejects");
    check(restoredWallCycleOwner.checkpoint() == beforeForeignRestore,
          "regional opening wall-cycle owner: failed restore retains the prior state");

    auto corruptWallPressureEpoch = wallPressureEpoch;
    ++corruptWallPressureEpoch.pressureEpoch.sourcePredictionFingerprint;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallPressureEpochIntegrity(
                corruptWallPressureEpoch);
        },
        "regional opening wall pressure epoch: corrupt nested provenance rejects");
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallPressureEpoch(
                wallPressureEpoch, zeroExchange,
                endpoint.momentumCycle.transport, endpoint.openingMetric,
                endpoint.pressureOperator, endpoint.basePressureOperator,
                endpoint.geometry, endpoint.sweep, endpoint.fragments,
                endpoint.topology, endpoint.volumeRates,
                endpoint.openingDefinitions, endpoint.openings,
                endpoint.resistanceDefinitions, endpoint.baseMetric,
                endpoint.openingMetric);
        },
        "regional opening wall pressure epoch: foreign exchange rejects full validation");
    auto wallPressureLimits =
        SceneFluidRegionalOpeningMomentumWallPressureEpochLimits{};
    wallPressureLimits.maximumOwnedBytes =
        wallPressureEpoch.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
                    exchange, endpoint.momentumCycle.transport,
                    endpoint.openingMetric, endpoint.pressureOperator,
                    endpoint.basePressureOperator, endpoint.geometry,
                    endpoint.sweep, endpoint.fragments, endpoint.topology,
                    endpoint.volumeRates, endpoint.openingDefinitions,
                    endpoint.openings, endpoint.resistanceDefinitions,
                    endpoint.baseMetric, endpoint.openingMetric, {},
                    endpoint.momentumCycle.pressureSettings,
                    wallPressureLimits));
        },
        "regional opening wall pressure epoch: late aggregate limit rejects publication");

    auto corruptAdjustment = adjustment;
    corruptAdjustment.controls.front().velocityMetersPerSecond.y += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
                corruptAdjustment);
        },
        "regional opening wall exchange: corrupt adjusted momentum rejects");

    auto adjustmentLimits =
        PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits{};
    adjustmentLimits.maximumFragments = adjustment.controls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                captureSceneFluidRegionalOpeningMomentumAdjustmentState(
                    exchange, endpoint.momentumCycle.transport,
                    endpoint.openingMetric, {}, adjustmentLimits));
        },
        "regional opening wall exchange: adjusted-state fragment limit rejects");
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumAdjustmentState(
                adjustment, zeroExchange,
                endpoint.momentumCycle.transport,
                endpoint.openingMetric);
        },
        "regional opening wall exchange: foreign wall exchange rejects adjusted state");
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
                adjustedPrediction, zeroAdjustment,
                endpoint.openingMetric, endpoint.geometry,
                endpoint.sweep, endpoint.fragments, endpoint.topology,
                endpoint.volumeRates, endpoint.openingDefinitions,
                endpoint.openings, endpoint.baseMetric,
                endpoint.openingMetric);
        },
        "regional opening wall exchange: foreign adjustment rejects opt-in prediction");

    auto limitedSettings = settings;
    limitedSettings.kinematicViscositySquareMetersPerSecond = 100.0;
    limitedSettings.maximumSubsteps = 1;
    const auto limited = exchangeSceneFluidRegionalOpeningMomentumWall(
        input, limitedSettings);
    check(!limited.diagnostics.accepted
              && limited.diagnostics.failureStage
                  == SceneFluidRegionWallFailureStage::SubstepLimit
              && limited.controlVolumes.empty() && limited.samples.empty()
              && limited.sourceInput == input,
          "regional opening wall exchange: unsafe explicit step retains source but publishes no adjusted state");
    expectRejected(
        [&] {
            static_cast<void>(
                captureSceneFluidRegionalOpeningMomentumAdjustmentState(
                    limited, endpoint.momentumCycle.transport,
                    endpoint.openingMetric));
        },
        "regional opening wall exchange: rejected exchange cannot publish adjusted momentum");
    expectRejected(
        [&] {
            static_cast<void>(
                acceptSceneFluidRegionalOpeningMomentumWallPressureEpoch(
                    limited, endpoint.momentumCycle.transport,
                    endpoint.openingMetric, endpoint.pressureOperator,
                    endpoint.basePressureOperator, endpoint.geometry,
                    endpoint.sweep, endpoint.fragments, endpoint.topology,
                    endpoint.volumeRates, endpoint.openingDefinitions,
                    endpoint.openings, endpoint.resistanceDefinitions,
                    endpoint.baseMetric, endpoint.openingMetric, {},
                    endpoint.momentumCycle.pressureSettings));
        },
        "regional opening wall pressure epoch: rejected wall exchange cannot enter pressure acceptance");

    auto corrupt = exchange;
    corrupt.samples.front().structureTraction.tractionPascals.y += 1.0;
    expectRejected(
        [&] {
            validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(
                corrupt);
        },
        "regional opening wall exchange: corrupt traction rejects by independent kernel replay");

    auto limits = SceneFluidRegionalOpeningMomentumWallExchangeLimits{};
    limits.maximumOwnedBytes = exchange.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                exchangeSceneFluidRegionalOpeningMomentumWall(
                    input, settings, limits));
        },
        "regional opening wall exchange: aggregate owned limit rejects after private kernel execution");

    auto foreignSettings = settings;
    foreignSettings.timeStepSeconds *= 2.0;
    expectRejected(
        [&] {
            static_cast<void>(
                exchangeSceneFluidRegionalOpeningMomentumWall(
                    input, foreignSettings));
        },
        "regional opening wall exchange: foreign time step rejects");
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
        testMimeticGeometryEpochTransition();
        testStaticSamplingAndTransfer();
        testOpeningAwareSamplingAndApplication();
        testPartialOpeningSamplingAndApplication();
        testMovingSamplingAndPower();
        testTransactionalLoadApplication();
        testMovingLoadApplication();
        testApplicationRollbackAndLimits();
        testAtomicOpeningLoadEpoch();
        testOpeningMomentumWallInput();
        testOpeningMomentumWallExchange();
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
