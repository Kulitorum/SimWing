#include "scene_fluid_mimetic_trace_system.h"
#include "scene_fluid_mimetic_condensed_trace_system.h"
#include "scene_fluid_mimetic_trace_solve.h"
#include "scene_structure.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
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
                     "FAIL: %s (actual %.17g expected %.17g tolerance %.3g)\n",
                     message, actual, expected, tolerance);
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
    const std::array<Vec2, 3> chart{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
    }};
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
        "sha256:scene-fluid-mimetic-trace-system";
    scene.metadata.exporterVersion =
        "scene-fluid-mimetic-trace-system-test/1";
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

Scene tiltedOpenScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-mimetic-trace-tilted-opening";
    scene.metadata.exporterVersion =
        "scene-fluid-mimetic-trace-system-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 1.2, 2.0, 0.3, 0.3, 900, true);
    scene.vertices[1].positionMeters.x = 2.6;
    scene.vertices[2].positionMeters.x = 2.8;
    scene.vertices[3].positionMeters.x = 2.7;
    return scene;
}

Scene disconnectedSameRegionScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-mimetic-trace-disconnected";
    scene.metadata.exporterVersion =
        "scene-fluid-mimetic-trace-system-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "disconnected-cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    addTetra(scene, 10, 500, 2, 1, 0.2, 0.8, 0.2, 0.2, 900);
    addTetra(scene, 20, 600, 2, 1, 3.2, 3.8, 0.2, 0.2, 901);
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
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

    explicit Fixture(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
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
              surface.definition, volumes, connectivity)) {}

    SceneFluidMimeticControlCellSet shells(
        const SceneFluidPressureFaceLinkSettings& linkSettings = {}) const {
        const auto links = buildSceneFluidPressureFaceLinks(
            surface.definition, state, grid(), transfer, epoch, caps,
            openingQuadrature, openingPatches, openingFaceCrossings,
            cappedFacePartitions, volumes, connectivity, pressureVolumes,
            linkSettings);
        return buildSceneFluidMimeticControlCells(
            surface.definition, state, grid(), epoch, caps,
            openingQuadrature, openingPatches, pressureVolumes, links);
    }
};

double dot(const std::vector<double>& first,
           const std::vector<double>& second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result += first[index] * second[index];
    }
    return result;
}

SceneFluidMimeticTraceSolveSettings strictSolveSettings() {
    SceneFluidMimeticTraceSolveSettings settings;
    settings.absoluteResidualTolerancePascalsMeters = 1.0e-12;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-11;
    settings.maximumIterations = 4000;
    return settings;
}

template<typename System>
void normalizeTraceGauges(const System& system,
                          std::vector<double>& values) {
    std::vector<double> gauges(system.componentCount, 0.0);
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        gauges[component] = values[
            system.componentGaugeTraceIndices[component]];
    }
    for (const auto& trace : system.traces) {
        values[trace.traceIndex] -= gauges[trace.componentIndex];
    }
    for (const std::size_t gauge : system.componentGaugeTraceIndices) {
        values[gauge] = 0.0;
    }
}

double maximumError(const std::vector<double>& first,
                    const std::vector<double>& second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result = std::max(result, std::abs(first[index] - second[index]));
    }
    return result;
}

void testNestedAssemblyAndMatrixFreeAction() {
    Fixture fixture(nestedScene());
    const auto shells = fixture.shells();
    const auto first = buildSceneFluidMimeticTraceSystem(shells);
    const auto repeated = buildSceneFluidMimeticTraceSystem(shells);
    check(first == repeated
              && first.version == sceneFluidMimeticTraceSystemVersion
              && first.fingerprint != 0
              && first.mimeticControlCellFingerprint == shells.fingerprint
              && first.localOperators.size() == shells.controlCells.size()
              && first.incidences.size() == shells.halfFaces.size()
              && first.halfFaceTraceIndices.size()
                  == shells.halfFaces.size()
              && first.componentCount == fixture.connectivity.components.size()
              && first.componentGaugeTraceIndices.size()
                  == first.componentCount
              && first.materialWallTraceCount
                  == shells.materialWallHalfFaceCount
              && 2 * first.sharedTraceCount
                  + first.materialWallTraceCount == first.incidences.size()
              && first.minimumPositiveOperatorDiagonal > 0.0
              && first.maximumOperatorDiagonal
                  >= first.minimumPositiveOperatorDiagonal,
          "nested shells assemble one deterministic bounded global trace topology");

    std::vector<double> componentConstants(first.componentCount, 0.0);
    for (std::size_t component = 0;
         component < componentConstants.size(); ++component) {
        componentConstants[component] = 1.25
            + 0.75 * static_cast<double>(component);
    }
    std::vector<double> constantTraces(first.traces.size(), 0.0);
    for (const auto& trace : first.traces) {
        constantTraces[trace.traceIndex] =
            componentConstants[trace.componentIndex];
    }
    const auto nullAction = applySceneFluidMimeticTraceOperator(
        first, constantTraces);
    double maximumNullAction = 0.0;
    for (const double value : nullAction) {
        maximumNullAction = std::max(maximumNullAction, std::abs(value));
    }
    check(maximumNullAction < 2.0e-11,
          "one constant trace mode per pressure component remains null");

    std::vector<double> firstVector(first.traces.size(), 0.0);
    std::vector<double> secondVector(first.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < first.traces.size(); ++trace) {
        firstVector[trace] = std::sin(0.37 * static_cast<double>(trace + 1));
        secondVector[trace] = std::cos(0.23 * static_cast<double>(trace + 2));
    }
    const auto firstAction = applySceneFluidMimeticTraceOperator(
        first, firstVector);
    const auto secondAction = applySceneFluidMimeticTraceOperator(
        first, secondVector);
    const double forward = dot(firstVector, secondAction);
    const double reverse = dot(secondVector, firstAction);
    const double energy = dot(firstVector, firstAction);
    const double symmetryTolerance = 2.0e-10
        * std::max({1.0, std::abs(forward), std::abs(reverse)});
    checkNear(forward, reverse, symmetryTolerance,
              "matrix-free global trace action is symmetric");
    check(std::isfinite(energy) && energy >= -symmetryTolerance,
          "matrix-free global trace action is positive semidefinite");

    for (std::size_t trace = 0;
         trace < std::min<std::size_t>(first.traces.size(), 8); ++trace) {
        std::vector<double> basis(first.traces.size(), 0.0);
        basis[trace] = 1.0;
        const auto action = applySceneFluidMimeticTraceOperator(first, basis);
        checkNear(action[trace], first.traces[trace].operatorDiagonal,
                  2.0e-11 * std::max(
                      1.0, first.traces[trace].operatorDiagonal),
                  "stored trace diagonal matches matrix-free basis action");
    }

    std::vector<double> sources(first.localOperators.size(), 0.0);
    std::size_t donor = 0;
    std::size_t receiver = 1;
    while (receiver < shells.controlCells.size()
           && shells.controlCells[receiver].componentIndex
               != shells.controlCells[donor].componentIndex) {
        ++receiver;
    }
    check(receiver < shells.controlCells.size(),
          "nested trace fixture has two controls in one pressure component");
    if (receiver < shells.controlCells.size()) {
        sources[donor] = 0.2;
        sources[receiver] = -0.2;
    }
    const auto rightHandSide =
        buildSceneFluidMimeticTraceRightHandSide(first, sources);
    std::vector<double> componentSums(first.componentCount, 0.0);
    for (const auto& trace : first.traces) {
        componentSums[trace.componentIndex] +=
            rightHandSide[trace.traceIndex];
    }
    double maximumComponentSum = 0.0;
    for (const double value : componentSums) {
        maximumComponentSum = std::max(
            maximumComponentSum, std::abs(value));
    }
    check(maximumComponentSum < 2.0e-13
              && std::ranges::any_of(
                  rightHandSide,
                  [](const double value) { return value != 0.0; }),
          "compatible cell sources produce a nonzero gauge-compatible trace RHS");
    validateSceneFluidMimeticTraceSystem(first, shells);
}

void testRejectedTwoPointOpeningStillBuildsHybridTrace() {
    Fixture fixture(tiltedOpenScene());
    SceneFluidPressureFaceLinkSettings rejectedSettings;
    rejectedSettings.minimumCenterDistanceMeters = 10.0;
    const auto shells = fixture.shells(rejectedSettings);
    const auto system = buildSceneFluidMimeticTraceSystem(shells);
    const auto openingTraceCount = std::ranges::count(
        system.traces,
        SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace,
        &SceneFluidMimeticTrace::kind);
    check(shells.openingHalfFaceCount == 2
              && openingTraceCount == 1
              && system.sharedTraceCount > 0
              && system.localOperators.size() == shells.controlCells.size(),
          "rejected two-point aperture becomes one paired hybrid trace");
    const auto opening = std::ranges::find(
        system.traces,
        SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace,
        &SceneFluidMimeticTrace::kind);
    check(opening != system.traces.end()
              && opening->incidenceCount == 2
              && opening->operatorDiagonal > 0.0,
          "hybrid aperture trace has two owners and positive local scaling");
}

void testGaugeFixedJacobiPcgRecovery() {
    Fixture fixture(nestedScene());
    const auto shells = fixture.shells();
    const auto system = buildSceneFluidMimeticTraceSystem(shells);
    std::vector<double> expected(system.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < expected.size(); ++trace) {
        const double sample = static_cast<double>(trace + 1);
        expected[trace] = 0.7 * std::sin(0.17 * sample)
            + 0.2 * std::cos(0.31 * sample) + 0.003 * sample;
    }
    normalizeTraceGauges(system, expected);
    const auto rightHandSide = applySceneFluidMimeticTraceOperator(
        system, expected);
    std::vector<double> firstTrace(expected.size(), 0.0);
    std::vector<double> secondTrace(expected.size(), 0.0);
    const auto settings = strictSolveSettings();
    const auto first = solveSceneFluidMimeticTraceSystem(
        system, rightHandSide, firstTrace, settings);
    const auto second = solveSceneFluidMimeticTraceSystem(
        system, rightHandSide, secondTrace, settings);
    check(first == second
              && firstTrace == secondTrace
              && first.compatible && first.converged && first.finite
              && first.iterationCount > 0
              && first.componentCount == system.componentCount,
          "gauge-fixed Jacobi-PCG converges deterministically in every component");
    check(maximumError(firstTrace, expected) < 3.0e-9,
          "mimetic trace solve recovers the manufactured gauge-fixed field");
    check(first.finalResidualL2PascalsMeters < 2.0e-11
              && first.finalResidualMaximumPascalsMeters < 2.0e-10,
          "mimetic trace solve closes its explicitly recomputed full residual");
    for (const auto& component : first.components) {
        checkNear(component.rightHandSideSumPascalsMeters,
                  0.0, 2.0e-11,
                  "manufactured trace RHS satisfies component compatibility");
        check(component.traceGaugeAfterPascals == 0.0
                  && firstTrace[component.gaugeTraceIndex] == 0.0,
              "mimetic trace solve commits exact component gauges");
    }
}

void testGlobalMaterialWallCondensation() {
    Fixture fixture(nestedScene());
    const auto shells = fixture.shells();
    const auto full = buildSceneFluidMimeticTraceSystem(shells);
    const auto first = buildSceneFluidMimeticCondensedTraceSystem(full);
    const auto repeated = buildSceneFluidMimeticCondensedTraceSystem(full);
    check(first == repeated
              && first.version
                  == sceneFluidMimeticCondensedTraceSystemVersion
              && first.fingerprint != 0
              && first.fullTraceSystemFingerprint == full.fingerprint
              && first.traces.size() == full.sharedTraceCount
              && first.eliminatedMaterialWallTraceCount
                  == full.materialWallTraceCount
              && first.fullTraceReducedIndices.size()
                  == full.traces.size()
              && first.localCondensations.size()
                  == full.localOperators.size()
              && first.componentGaugeTraceIndices.size()
                  == full.componentCount
              && first.minimumPositiveOperatorDiagonal > 0.0
              && first.maximumOperatorDiagonal
                  >= first.minimumPositiveOperatorDiagonal,
          "global wall condensation deterministically retains only shared traces");

    std::vector<double> componentConstants(first.componentCount, 0.0);
    for (std::size_t component = 0;
         component < first.componentCount; ++component) {
        componentConstants[component] = 0.75
            + 0.5 * static_cast<double>(component);
    }
    std::vector<double> constant(first.traces.size(), 0.0);
    for (const auto& trace : first.traces) {
        constant[trace.traceIndex] =
            componentConstants[trace.componentIndex];
    }
    const auto nullAction =
        applySceneFluidMimeticCondensedTraceOperator(
            first, full, constant);
    check(std::ranges::all_of(
              nullAction,
              [](const double value) { return std::abs(value) < 3.0e-11; }),
          "global wall condensation preserves every component constant null mode");

    std::vector<double> firstVector(first.traces.size(), 0.0);
    std::vector<double> secondVector(first.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < first.traces.size(); ++trace) {
        firstVector[trace] = std::sin(
            0.19 * static_cast<double>(trace + 1));
        secondVector[trace] = std::cos(
            0.27 * static_cast<double>(trace + 2));
    }
    const auto firstAction =
        applySceneFluidMimeticCondensedTraceOperator(
            first, full, firstVector);
    const auto secondAction =
        applySceneFluidMimeticCondensedTraceOperator(
            first, full, secondVector);
    const double forward = dot(firstVector, secondAction);
    const double reverse = dot(secondVector, firstAction);
    const double energy = dot(firstVector, firstAction);
    const double tolerance = 5.0e-10
        * std::max({1.0, std::abs(forward), std::abs(reverse)});
    checkNear(forward, reverse, tolerance,
              "global wall-condensed action is symmetric");
    check(std::isfinite(energy) && energy >= -tolerance,
          "global wall-condensed action is positive semidefinite");
    for (std::size_t trace = 0;
         trace < std::min<std::size_t>(first.traces.size(), 8); ++trace) {
        std::vector<double> basis(first.traces.size(), 0.0);
        basis[trace] = 1.0;
        const auto action =
            applySceneFluidMimeticCondensedTraceOperator(
                first, full, basis);
        checkNear(action[trace], first.traces[trace].operatorDiagonal,
                  3.0e-10 * std::max(
                      1.0, first.traces[trace].operatorDiagonal),
                  "global condensed diagonal matches its basis action");
    }

    std::vector<double> expectedFull(full.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < expectedFull.size(); ++trace) {
        expectedFull[trace] = 0.4 * std::sin(
            0.13 * static_cast<double>(trace + 1))
            + 0.03 * static_cast<double>(trace + 1);
    }
    for (std::size_t component = 0;
         component < first.componentCount; ++component) {
        const std::size_t reducedGauge =
            first.componentGaugeTraceIndices[component];
        const std::size_t fullGauge =
            first.traces[reducedGauge].fullTraceIndex;
        const double gaugeValue = expectedFull[fullGauge];
        for (const auto& trace : full.traces) {
            if (trace.componentIndex == component) {
                expectedFull[trace.traceIndex] -= gaugeValue;
            }
        }
        expectedFull[fullGauge] = 0.0;
    }
    std::vector<double> expectedReduced(first.traces.size(), 0.0);
    for (const auto& trace : first.traces) {
        expectedReduced[trace.traceIndex] =
            expectedFull[trace.fullTraceIndex];
    }
    const auto fullRightHandSide = applySceneFluidMimeticTraceOperator(
        full, expectedFull);
    const auto condensedRightHandSide =
        condenseSceneFluidMimeticTraceRightHandSide(
            first, full, fullRightHandSide);
    const auto expectedReducedAction =
        applySceneFluidMimeticCondensedTraceOperator(
            first, full, expectedReduced);
    check(maximumError(
              condensedRightHandSide, expectedReducedAction) < 3.0e-10,
          "global RHS condensation matches the reduced manufactured action");
    std::vector<double> firstSolvedReduced(first.traces.size(), 0.0);
    std::vector<double> secondSolvedReduced(first.traces.size(), 0.0);
    const auto firstSolve = solveSceneFluidMimeticCondensedTraceSystem(
        first, full, condensedRightHandSide, firstSolvedReduced,
        strictSolveSettings());
    const auto secondSolve = solveSceneFluidMimeticCondensedTraceSystem(
        first, full, condensedRightHandSide, secondSolvedReduced,
        strictSolveSettings());
    check(firstSolve == secondSolve
              && firstSolvedReduced == secondSolvedReduced
              && firstSolve.compatible && firstSolve.converged
              && firstSolve.finite && firstSolve.iterationCount > 0
              && firstSolve.traceSystemFingerprint == first.fingerprint,
          "wall-condensed Jacobi-PCG converges transactionally and deterministically");
    check(maximumError(firstSolvedReduced, expectedReduced) < 4.0e-9,
          "wall-condensed solve recovers the manufactured shared traces");
    const auto reconstructed = reconstructSceneFluidMimeticFullTraces(
        first, full, fullRightHandSide, firstSolvedReduced);
    const auto exactReconstructed = reconstructSceneFluidMimeticFullTraces(
        first, full, fullRightHandSide, expectedReduced);
    check(maximumError(exactReconstructed, expectedFull) < 5.0e-9,
          "global wall reconstruction recovers every manufactured full trace");
    const auto reconstructedAction = applySceneFluidMimeticTraceOperator(
        full, reconstructed);
    check(maximumError(reconstructedAction, fullRightHandSide) < 4.0e-9,
          "solved and reconstructed traces close every original full equation");

    auto truncatedSettings = strictSolveSettings();
    truncatedSettings.absoluteResidualTolerancePascalsMeters = 1.0e-30;
    truncatedSettings.relativeResidualTolerance = 0.0;
    truncatedSettings.maximumIterations = 1;
    std::vector<double> truncatedWarmStart(first.traces.size(), 3.25);
    const auto originalTruncatedWarmStart = truncatedWarmStart;
    const auto truncated = solveSceneFluidMimeticCondensedTraceSystem(
        first, full, condensedRightHandSide, truncatedWarmStart,
        truncatedSettings);
    check(truncated.compatible && !truncated.converged
              && truncated.finite && truncated.iterationCount == 1
              && truncatedWarmStart == originalTruncatedWarmStart,
          "truncated wall-condensed solve rolls its warm start back exactly");
    auto incompatible = condensedRightHandSide;
    incompatible.front() += 1.0e-4;
    const auto incompatibleDiagnostics =
        solveSceneFluidMimeticCondensedTraceSystem(
            first, full, incompatible, truncatedWarmStart,
            strictSolveSettings());
    check(!incompatibleDiagnostics.compatible
              && !incompatibleDiagnostics.converged
              && incompatibleDiagnostics.finite
              && truncatedWarmStart == originalTruncatedWarmStart,
          "incompatible wall-condensed RHS is reported without state mutation");
    validateSceneFluidMimeticCondensedTraceSystem(first, full);
}

void testGlobalWallCondensationLimitsAndCorruption() {
    Fixture fixture(nestedScene());
    const auto full = buildSceneFluidMimeticTraceSystem(
        fixture.shells());
    const auto accepted =
        buildSceneFluidMimeticCondensedTraceSystem(full);
    auto corrupt = accepted;
    corrupt.traces.front().operatorDiagonal += 0.01;
    expectInvalid(
        [&] { validateSceneFluidMimeticCondensedTraceSystem(
            corrupt, full); },
        "global wall condensation rejects fingerprinted diagonal corruption");
    expectInvalid(
        [&] {
            std::vector<double> traces(corrupt.traces.size(), 0.0);
            static_cast<void>(solveSceneFluidMimeticCondensedTraceSystem(
                corrupt, full,
                std::vector<double>(corrupt.traces.size(), 0.0),
                traces, strictSolveSettings()));
        },
        "condensed PCG validates its immutable product before trusted iteration");

    SceneFluidMimeticCondensedTraceSystemLimits limits;
    limits.maximumReducedTraces = accepted.traces.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticCondensedTraceSystem(full, limits)); },
        "global wall condensation bounds reduced trace count");
    limits = {};
    limits.maximumLocalCondensations =
        accepted.localCondensations.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticCondensedTraceSystem(full, limits)); },
        "global wall condensation bounds local condensation count");
    limits = {};
    limits.maximumLocalCondensationBytes =
        accepted.localCondensationStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticCondensedTraceSystem(full, limits)); },
        "global wall condensation bounds nested linear storage");
    limits = {};
    limits.maximumOwnedBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticCondensedTraceSystem(full, limits)); },
        "global wall condensation bounds aggregate storage");

    expectInvalid(
        [&] { static_cast<void>(
            applySceneFluidMimeticCondensedTraceOperator(
                accepted, full,
                std::vector<double>(accepted.traces.size() - 1, 0.0))); },
        "global wall-condensed action rejects a short reduced field");
    expectInvalid(
        [&] { static_cast<void>(
            condenseSceneFluidMimeticTraceRightHandSide(
                accepted, full,
                std::vector<double>(
                    full.traces.size(),
                    std::numeric_limits<double>::infinity()))); },
        "global wall condensation rejects a non-finite full RHS");
}

void testSourceDrivenTraceBalance() {
    Fixture fixture(nestedScene());
    const auto shells = fixture.shells();
    const auto system = buildSceneFluidMimeticTraceSystem(shells);
    std::vector<double> sources(system.localOperators.size(), 0.0);
    std::size_t receiver = 1;
    while (receiver < shells.controlCells.size()
           && shells.controlCells[receiver].componentIndex
               != shells.controlCells[0].componentIndex) {
        ++receiver;
    }
    check(receiver < sources.size(),
          "source-driven solve fixture has two controls in one component");
    if (receiver < sources.size()) {
        sources[0] = 0.2;
        sources[receiver] = -0.2;
    }
    const auto rightHandSide =
        buildSceneFluidMimeticTraceRightHandSide(system, sources);
    std::vector<double> traces(system.traces.size(), 0.0);
    const auto diagnostics = solveSceneFluidMimeticTraceSystem(
        system, rightHandSide, traces, strictSolveSettings());
    const auto balance = evaluateSceneFluidMimeticTraceSystem(
        system, traces, sources);
    check(diagnostics.compatible && diagnostics.converged
              && diagnostics.finite
              && balance.maximumTraceFluxImbalance < 2.0e-10
              && balance.maximumCellConservationResidual < 2.0e-13,
          "source-driven trace solve closes local conservation and every hybrid trace");

    const auto condensed =
        buildSceneFluidMimeticCondensedTraceSystem(system);
    const auto reducedRightHandSide =
        condenseSceneFluidMimeticTraceRightHandSide(
            condensed, system, rightHandSide);
    std::vector<double> reducedTraces(condensed.traces.size(), 0.0);
    const auto reducedDiagnostics =
        solveSceneFluidMimeticCondensedTraceSystem(
            condensed, system, reducedRightHandSide, reducedTraces,
            strictSolveSettings());
    const auto reconstructed = reconstructSceneFluidMimeticFullTraces(
        condensed, system, rightHandSide, reducedTraces);
    const auto reducedBalance = evaluateSceneFluidMimeticTraceSystem(
        system, reconstructed, sources);
    check(reducedDiagnostics.compatible && reducedDiagnostics.converged
              && reducedDiagnostics.finite
              && reducedBalance.maximumTraceFluxImbalance < 3.0e-10
              && reducedBalance.maximumCellConservationResidual < 3.0e-13,
          "source-driven condensed solve reconstructs a conservative full trace field");
}

void testTraceSolveRollbackAndValidation() {
    Fixture fixture(nestedScene());
    const auto system = buildSceneFluidMimeticTraceSystem(
        fixture.shells());
    std::vector<double> expected(system.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < expected.size(); ++trace) {
        expected[trace] = std::sin(
            0.23 * static_cast<double>(trace + 1));
    }
    const auto compatible = applySceneFluidMimeticTraceOperator(
        system, expected);

    auto incompatible = compatible;
    incompatible.front() += 1.0e-4;
    std::vector<double> traces(system.traces.size(), 3.5);
    const auto originalTraces = traces;
    const auto incompatibleDiagnostics =
        solveSceneFluidMimeticTraceSystem(
            system, incompatible, traces, strictSolveSettings());
    check(!incompatibleDiagnostics.compatible
              && !incompatibleDiagnostics.converged
              && incompatibleDiagnostics.finite
              && incompatibleDiagnostics
                    .maximumAbsoluteComponentCompatibilityPascalsMeters
                    > 9.0e-5
              && traces == originalTraces,
          "incompatible trace RHS is reported without warm-state mutation");

    auto truncatedSettings = strictSolveSettings();
    truncatedSettings.absoluteResidualTolerancePascalsMeters = 1.0e-30;
    truncatedSettings.relativeResidualTolerance = 0.0;
    truncatedSettings.maximumIterations = 1;
    const auto truncated = solveSceneFluidMimeticTraceSystem(
        system, compatible, traces, truncatedSettings);
    check(truncated.compatible && !truncated.converged
              && truncated.finite && truncated.iterationCount == 1
              && traces == originalTraces,
          "truncated trace solve rolls its warm start back exactly");

    auto corrected = compatible;
    corrected.front() += 1.0e-13;
    auto correctedTraces = originalTraces;
    const auto correctedDiagnostics = solveSceneFluidMimeticTraceSystem(
        system, corrected, correctedTraces, strictSolveSettings());
    check(correctedDiagnostics.compatible
              && correctedDiagnostics.converged
              && correctedDiagnostics
                    .maximumAbsoluteComponentCompatibilityPascalsMeters
                    > 0.0,
          "trace solve removes only an admitted component roundoff defect");

    auto invalidSettings = strictSolveSettings();
    invalidSettings.relativeResidualTolerance = -1.0;
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidMimeticTraceSystem(
            system, compatible, traces, invalidSettings)); },
        "trace solve rejects invalid residual settings");
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidMimeticTraceSystem(
            system,
            std::span<const double>(compatible).first(
                compatible.size() - 1),
            traces, strictSolveSettings())); },
        "trace solve rejects a short right-hand side");
    auto nonFiniteTraces = originalTraces;
    nonFiniteTraces.front() =
        std::numeric_limits<double>::quiet_NaN();
    expectInvalid(
        [&] { static_cast<void>(solveSceneFluidMimeticTraceSystem(
            system, compatible, nonFiniteTraces,
            strictSolveSettings())); },
        "trace solve rejects a non-finite warm start");

    std::vector<double> overflowingWarmStart(
        system.traces.size(), std::numeric_limits<double>::max());
    for (const std::size_t gauge : system.componentGaugeTraceIndices) {
        overflowingWarmStart[gauge] =
            -std::numeric_limits<double>::max();
    }
    const auto originalOverflowingWarmStart = overflowingWarmStart;
    const auto overflowDiagnostics = solveSceneFluidMimeticTraceSystem(
        system, compatible, overflowingWarmStart, strictSolveSettings());
    check(overflowDiagnostics.compatible
              && !overflowDiagnostics.converged
              && !overflowDiagnostics.finite
              && overflowingWarmStart == originalOverflowingWarmStart,
          "non-finite trace-solve arithmetic preserves the warm start exactly");
}

void testLimitsCorruptionAndInputValidation() {
    Fixture fixture(nestedScene());
    const auto shells = fixture.shells();
    const auto accepted = buildSceneFluidMimeticTraceSystem(shells);
    auto corrupt = accepted;
    corrupt.traces.front().operatorDiagonal += 0.01;
    expectInvalid(
        [&] { validateSceneFluidMimeticTraceSystemIntegrity(corrupt); },
        "mimetic trace-system integrity rejects nested corruption");

    SceneFluidMimeticTraceSystemLimits limits;
    limits.maximumTraces = accepted.traces.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, {}, limits)); },
        "mimetic trace assembly bounds trace count");
    limits = {};
    limits.maximumIncidences = accepted.incidences.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, {}, limits)); },
        "mimetic trace assembly bounds incidence count");
    limits = {};
    limits.maximumLocalOperators = accepted.localOperators.size() - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, {}, limits)); },
        "mimetic trace assembly bounds local operator count");
    limits = {};
    limits.maximumLocalOperatorBytes =
        accepted.localOperatorStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, {}, limits)); },
        "mimetic trace assembly bounds compact local storage");
    limits = {};
    limits.maximumOwnedBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, {}, limits)); },
        "mimetic trace assembly bounds aggregate storage");

    SceneFluidMimeticTraceSystemSettings badSettings;
    badSettings.localCell.maximumHalfFaces = 3;
    expectInvalid(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceSystem(shells, badSettings)); },
        "mimetic trace assembly rejects invalid local settings");
    expectInvalid(
        [&] { static_cast<void>(applySceneFluidMimeticTraceOperator(
            accepted,
            std::vector<double>(accepted.traces.size() - 1, 0.0))); },
        "mimetic trace action rejects a short trace field");
    expectInvalid(
        [&] { static_cast<void>(
            buildSceneFluidMimeticTraceRightHandSide(
                accepted,
                std::vector<double>(
                    accepted.localOperators.size(),
                    std::numeric_limits<double>::infinity()))); },
        "mimetic trace RHS rejects a non-finite source field");

    Fixture disconnectedFixture(disconnectedSameRegionScene());
    const auto disconnectedShells = disconnectedFixture.shells();
    check(disconnectedShells.readyControlCellCount
              == disconnectedShells.controlCells.size(),
          "disconnected same-region fixture still has complete local shells");
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidMimeticTraceSystem(
            disconnectedShells)); },
        "one component gauge rejects disconnected trace topology");
}

} // namespace

int main() {
    try {
        testNestedAssemblyAndMatrixFreeAction();
        testRejectedTwoPointOpeningStillBuildsHybridTrace();
        testGaugeFixedJacobiPcgRecovery();
        testGlobalMaterialWallCondensation();
        testGlobalWallCondensationLimitsAndCorruption();
        testSourceDrivenTraceBalance();
        testTraceSolveRollbackAndValidation();
        testLimitsCorruptionAndInputValidation();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene fluid mimetic trace-system check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene fluid mimetic trace-system checks passed");
    return 0;
}
