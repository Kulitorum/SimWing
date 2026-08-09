#include "scene_fluid_mimetic_trace_system.h"
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
