#include "coupling.h"
#include "face_resolved_bridge.h"
#include "fluid_structure_bridge.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using simwing::fsi::ConservativeMacroStepCoupling;
using simwing::fsi::CouplingNodeKinematics;
using simwing::fsi::CouplingSurfaceNodeDefinition;
using simwing::fsi::CouplingSurfaceTriangleDefinition;
using simwing::fsi::Structure;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureStepSettings;
using simwing::fsi::StructureVector3;
using simwing::fsi::PlanarFaceResolvedFluidStructureBridge;
using simwing::fsi::PlanarFaceCorrespondenceMode;
using simwing::fsi::PlanarFaceResolvedBridgeSettings;
using simwing::fsi::UniformFluidStructureBridge;
using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionDiagnostics;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::PeriodicCartesianGrid;
using simwing::fsi::fluid::PlanarMovingControlVolume;
using simwing::fsi::fluid::evaluatePlanarCutSurfacePressure;
using simwing::fsi::fluid::projectVelocityWithMovingInterfaces;

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
                     "FAIL: %s (actual %.17g, expected %.17g, tolerance %.3g)\n",
                     message, actual, expected, tolerance);
        ++failures;
    }
}

void checkVectorNear(const StructureVector3& actual,
                     const StructureVector3& expected,
                     const double tolerance,
                     const char* message) {
    if (!std::isfinite(actual.x) || !std::isfinite(actual.y)
        || !std::isfinite(actual.z)
        || std::abs(actual.x - expected.x) > tolerance
        || std::abs(actual.y - expected.y) > tolerance
        || std::abs(actual.z - expected.z) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual [%.17g %.17g %.17g], expected [%.17g %.17g %.17g])\n",
                     message, actual.x, actual.y, actual.z,
                     expected.x, expected.y, expected.z);
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

StructureDefinition pistonDefinition(const double heightMeters = 3.0) {
    StructureDefinition definition;
    definition.nodes = {
        {{3.0, 0.0, 0.0}, 2.0, false},
        {{3.0, 2.0, 0.0}, 1.0, false},
        {{3.0, 2.0, heightMeters}, 2.0, false},
        {{3.0, 0.0, heightMeters}, 1.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

StructureDefinition movingPlaneDefinition(
    const GridFaceAxis axis,
    const double planeCoordinateMeters) {
    StructureDefinition definition;
    switch (axis) {
    case GridFaceAxis::X:
        definition.nodes = {
            {{planeCoordinateMeters, 0.0, 0.0}, 2.0, false},
            {{planeCoordinateMeters, 2.0, 0.0}, 1.0, false},
            {{planeCoordinateMeters, 2.0, 3.0}, 2.0, false},
            {{planeCoordinateMeters, 0.0, 3.0}, 1.0, false},
        };
        break;
    case GridFaceAxis::Y:
        definition.nodes = {
            {{0.0, planeCoordinateMeters, 0.0}, 2.0, false},
            {{0.0, planeCoordinateMeters, 3.0}, 1.0, false},
            {{4.0, planeCoordinateMeters, 3.0}, 2.0, false},
            {{4.0, planeCoordinateMeters, 0.0}, 1.0, false},
        };
        break;
    case GridFaceAxis::Z:
        definition.nodes = {
            {{0.0, 0.0, planeCoordinateMeters}, 2.0, false},
            {{4.0, 0.0, planeCoordinateMeters}, 1.0, false},
            {{4.0, 2.0, planeCoordinateMeters}, 2.0, false},
            {{0.0, 2.0, planeCoordinateMeters}, 1.0, false},
        };
        break;
    }
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

std::vector<CouplingSurfaceNodeDefinition> pistonNodes() {
    return {{40, 3}, {10, 0}, {30, 2}, {20, 1}};
}

std::vector<CouplingSurfaceTriangleDefinition> pistonTriangles() {
    return {
        {2000, {10, 30, 40}},
        {1000, {10, 20, 30}},
    };
}

PeriodicCartesianGrid pistonGrid() {
    return PeriodicCartesianGrid({8, 2, 3}, {}, {4.0, 2.0, 3.0});
}

std::vector<GridFaceMovingInterface> slabFaces(
    const PeriodicCartesianGrid& grid) {
    std::vector<GridFaceMovingInterface> result;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            result.push_back(
                {100, 1, 2, GridFaceAxis::X, 2, j, k, 0.25});
            result.push_back(
                {200, 2, 1, GridFaceAxis::X, 6, j, k, 0.25});
        }
    }
    return result;
}

MovingInterfaceProjectionDiagnostics fluidPistonDiagnostics(
    const bool disturbed = false) {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(grid, slabFaces(grid));
    MacVelocityField velocity(grid);
    if (disturbed) {
        for (std::size_t index = 0;
             index < velocity.xFaces().size(); ++index) {
            const double sample = static_cast<double>(index + 1);
            velocity.xFaces()[index] = 0.2 * std::sin(0.31 * sample);
            velocity.yFaces()[index] = 0.15 * std::cos(0.23 * sample);
            velocity.zFaces()[index] =
                0.1 * std::sin(0.17 * sample + 0.2);
        }
    } else {
        std::ranges::fill(velocity.xFaces(), 0.25);
    }
    CellScalarField pressure(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                pressure.values()[grid.cellIndex(i, j, k)] =
                    i >= 2 && i < 6 ? 110.0 : 0.0;
            }
        }
    }
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.4;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    return projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, settings);
}

simwing::fsi::fluid::PorousSurfaceTractionDiagnostics
porousSurfaceTraction(const double crossingFraction = 0.4) {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(grid, slabFaces(grid));
    const auto counts = grid.cellCounts();
    std::vector<simwing::fsi::fluid::PorousGridFaceCrossing> porous;
    std::vector<simwing::fsi::fluid::GridFacePressureJump> prescribed;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            porous.push_back({
                300, 10, 11, GridFaceAxis::X, 3, j, k,
                crossingFraction, 0.1, {10.0, 0.0}});
            prescribed.push_back({
                400, 11, 10, GridFaceAxis::X, 5, j, k,
                1.5, 0.6});
        }
    }
    simwing::fsi::fluid::MovingPorousProjectionSettings settings;
    settings.movingProjection.projection.densityKgPerCubicMeter = 1.2;
    settings.movingProjection.projection.timeStepSeconds = 0.4;
    settings.movingProjection.projection.absoluteResidualTolerance =
        1.0e-11;
    settings.movingProjection.projection.relativeResidualTolerance =
        1.0e-13;
    settings.movingProjection.projection.maximumIterations = 1000;
    settings.movingProjection
        .absoluteRegionVolumeRateToleranceCubicMetersPerSecond = 1.0e-12;
    settings.iteration.absoluteNormalVelocityToleranceMetersPerSecond =
        1.0e-12;
    settings.iteration.relativeNormalVelocityTolerance = 1.0e-12;
    settings.iteration.absolutePressureJumpTolerancePascals = 1.0e-11;
    settings.iteration.relativePressureJumpTolerance = 1.0e-12;
    settings.iteration.relaxation = 0.5;
    settings.iteration.maximumNonlinearIterations = 100;
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    const simwing::fsi::fluid::SharpPressureJumpField pressureSource(
        grid, std::move(prescribed));
    const auto diagnostics =
        simwing::fsi::fluid::projectVelocityWithMovingAndPorousInterfaces(
            grid, velocity, pressure, interfaces,
            porous, pressureSource, settings);
    return simwing::fsi::fluid::evaluatePorousSurfaceTraction(
        grid, diagnostics,
        settings.movingProjection.projection.timeStepSeconds);
}

FaceAlignedMovingInterface openPlaneInterfaces(
    const GridFaceAxis axis,
    const std::size_t movingPlaneCoordinate,
    const double speedMetersPerSecond = 0.25) {
    const auto grid = pistonGrid();
    const auto counts = grid.cellCounts();
    std::vector<GridFaceMovingInterface> faces;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t coordinate = axis == GridFaceAxis::X
                    ? i : (axis == GridFaceAxis::Y ? j : k);
                if (coordinate == movingPlaneCoordinate) {
                    faces.push_back({
                        300, 9, 9, axis, i, j, k,
                        speedMetersPerSecond,
                    });
                }
            }
        }
    }
    return FaceAlignedMovingInterface(grid, std::move(faces));
}

MovingInterfaceProjectionDiagnostics openPlaneDiagnostics(
    const GridFaceAxis axis,
    const std::size_t movingPlaneCoordinate,
    const double speedMetersPerSecond = 0.25) {
    const auto grid = pistonGrid();
    const auto interfaces = openPlaneInterfaces(
        axis, movingPlaneCoordinate, speedMetersPerSecond);
    MacVelocityField velocity(grid);
    CellScalarField pressure(grid);
    MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.4;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, settings);
}

MovingInterfaceProjectionDiagnostics openPistonDiagnostics(
    const std::size_t movingPlaneCoordinate = 6,
    const double speedMetersPerSecond = 0.25) {
    return openPlaneDiagnostics(
        GridFaceAxis::X, movingPlaneCoordinate,
        speedMetersPerSecond);
}

std::vector<CouplingNodeKinematics> translatingKinematics(
    const simwing::fsi::ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const double speedMetersPerSecond = 0.25) {
    auto result = transfer.captureKinematics(structure);
    for (auto& node : result) {
        node.velocityMetersPerSecond = {speedMetersPerSecond, 0.0, 0.0};
    }
    return result;
}

std::vector<CouplingNodeKinematics> movingPlaneKinematics(
    const simwing::fsi::ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const double physicalPlaneCoordinateMeters,
    const double speedMetersPerSecond = 0.25,
    const GridFaceAxis axis = GridFaceAxis::X) {
    auto result = transfer.captureKinematics(structure);
    for (auto& node : result) {
        switch (axis) {
        case GridFaceAxis::X:
            node.positionMeters.x = physicalPlaneCoordinateMeters;
            node.velocityMetersPerSecond = {
                speedMetersPerSecond, 0.0, 0.0};
            break;
        case GridFaceAxis::Y:
            node.positionMeters.y = physicalPlaneCoordinateMeters;
            node.velocityMetersPerSecond = {
                0.0, speedMetersPerSecond, 0.0};
            break;
        case GridFaceAxis::Z:
            node.positionMeters.z = physicalPlaneCoordinateMeters;
            node.velocityMetersPerSecond = {
                0.0, 0.0, speedMetersPerSecond};
            break;
        }
    }
    return result;
}

std::vector<CouplingNodeKinematics> translatingKinematics(
    const UniformFluidStructureBridge& bridge,
    const Structure& structure,
    const double speedMetersPerSecond = 0.25) {
    auto result = bridge.transfer().captureKinematics(structure);
    for (auto& node : result) {
        node.velocityMetersPerSecond = {speedMetersPerSecond, 0.0, 0.0};
    }
    return result;
}

void testAnalyticStableIdBridgeAndStructuralAcceptance() {
    Structure structure(pistonDefinition());
    UniformFluidStructureBridge bridge(
        structure, 200, pistonNodes(), pistonTriangles());
    const auto fluid = fluidPistonDiagnostics();
    const auto kinematics = translatingKinematics(bridge, structure);
    const auto first = bridge.evaluate(fluid, kinematics);
    const auto second = bridge.evaluate(fluid, kinematics);
    const auto& diagnostics = first.diagnostics();

    check(first == second,
          "bridge: identical stable-ID exchanges replay bit-for-bit");
    check(bridge.fluidSurfaceStableId() == 200
              && diagnostics.fluidSurfaceStableId == 200
              && diagnostics.fluidFaceCount == 6
              && diagnostics.structureTriangleCount == 2,
          "bridge: one canonical fluid surface maps to both structural triangles");
    checkNear(diagnostics.fluidAreaSquareMeters, 6.0, 0.0,
              "bridge: fluid surface has the analytic piston area");
    checkNear(diagnostics.structureAreaSquareMeters, 6.0, 1.0e-15,
              "bridge: structural surface has the analytic piston area");
    checkNear(diagnostics.areaResidualSquareMeters, 0.0, 1.0e-15,
              "bridge: independent fluid and structure area ledgers close");
    checkVectorNear(diagnostics.uniformPressureTractionPascals,
                    {110.0, 0.0, 0.0}, 0.0,
                    "bridge: stable surface force reconstructs the uniform traction");
    checkVectorNear(diagnostics.fluidPressureForceNewtons,
                    {660.0, 0.0, 0.0}, 0.0,
                    "bridge: source pressure force remains explicit");
    checkVectorNear(diagnostics.structureSurfaceForceNewtons,
                    {660.0, 0.0, 0.0}, 1.0e-13,
                    "bridge: structural transfer receives the source force");
    checkNear(diagnostics.fluidPressurePowerWatts, 165.0, 0.0,
              "bridge: source fluid pressure power is analytic");
    checkNear(diagnostics.structureSurfacePowerWatts, 165.0, 3.0e-14,
              "bridge: structural surface power matches the fluid ledger");
    check(diagnostics.forceResidualNormNewtons < 2.0e-13
              && std::abs(diagnostics.powerResidualWatts) < 4.0e-14
              && diagnostics.finite,
          "bridge: force and power residuals meet their budgets");

    ConservativeMacroStepCoupling coupling(bridge.transfer());
    const std::array<double, 2> offsets{0.0, 0.4};
    const std::array samples{
        first.transferResult(), second.transferResult()};
    const auto integrated = coupling.integrate(offsets, samples);
    checkVectorNear(
        integrated.diagnostics().integratedSurfaceImpulseNewtonSeconds,
        {264.0, 0.0, 0.0}, 2.0e-13,
        "bridge: fluid pressure reaches the temporal impulse exchange");
    checkNear(integrated.diagnostics().integratedSurfaceWorkJoules,
              66.0, 5.0e-14,
              "bridge: fluid pressure power reaches temporal work");

    StructureStepSettings step;
    step.timeStepSeconds = 0.4;
    step.substeps = 1;
    step.constraintIterations = 0;
    step.gravityMetersPerSecondSquared = {};
    step.velocityDampingPerSecond = 0.0;
    const auto accepted = coupling.advanceStructure(
        structure, integrated, step);
    checkVectorNear(accepted.linearMomentumKgMetersPerSecond,
                    {264.0, 0.0, 0.0}, 3.0e-13,
                    "bridge: the accepted XPBD state carries the fluid impulse");
    const auto states = structure.nodeStates();
    checkNear(states[0].velocityMetersPerSecond.x, 44.0, 5.0e-14,
              "bridge: tributary mass converts shared-node impulse correctly");
    checkNear(states[1].velocityMetersPerSecond.x, 44.0, 5.0e-14,
              "bridge: tributary mass converts boundary-node impulse correctly");
}

void testStrictRejectionContracts() {
    Structure structure(pistonDefinition());
    UniformFluidStructureBridge bridge(
        structure, 200, pistonNodes(), pistonTriangles());
    const auto fluid = fluidPistonDiagnostics();
    const auto kinematics = translatingKinematics(bridge, structure);

    auto nonuniform = fluid;
    nonuniform.surfaces[1].maximumPressureTractionDeviationPascals = 1.0;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(nonuniform, kinematics)); },
        "validation: unresolved nonuniform fluid traction is rejected");

    auto failed = fluid;
    failed.projection.converged = false;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(failed, kinematics)); },
        "validation: unaccepted fluid projection is rejected");

    auto wrongPowerKinematics = translatingKinematics(
        bridge, structure, 0.5);
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            fluid, wrongPowerKinematics)); },
        "validation: inconsistent fluid and structural interface power is rejected");

    Structure wrongAreaStructure(pistonDefinition(4.0));
    UniformFluidStructureBridge wrongAreaBridge(
        wrongAreaStructure, 200, pistonNodes(), pistonTriangles());
    const auto wrongAreaKinematics = translatingKinematics(
        wrongAreaBridge, wrongAreaStructure);
    expectRejected(
        [&] { static_cast<void>(wrongAreaBridge.evaluate(
            fluid, wrongAreaKinematics)); },
        "validation: fluid and structural area mismatch is rejected");

    expectRejected(
        [&] { UniformFluidStructureBridge missing(
            structure, 0, pistonNodes(), pistonTriangles()); },
        "validation: zero bound surface stable ID is rejected");

    UniformFluidStructureBridge absent(
        structure, 999, pistonNodes(), pistonTriangles());
    expectRejected(
        [&] { static_cast<void>(absent.evaluate(fluid, kinematics)); },
        "validation: absent fluid surface stable ID is rejected");
}

void testFaceResolvedNonuniformPressureTransfer() {
    Structure structure(pistonDefinition());
    const auto reference = fluidPistonDiagnostics();
    PlanarFaceResolvedFluidStructureBridge bridge(
        structure, 200, pistonNodes(), pistonTriangles(), reference.faces);
    const auto kinematics = translatingKinematics(
        bridge.transfer(), structure);
    const auto uniform = bridge.evaluate(reference, kinematics);
    checkNear(bridge.referenceAreaSquareMeters(), 6.0, 1.0e-15,
              "face bridge: clipped reference correspondence covers 6 square metres");
    check(bridge.overlapPatchCount() >= 6
              && uniform.diagnostics().overlapPatchCount
                  == bridge.overlapPatchCount(),
          "face bridge: fluid tiles are split into stable triangle overlap patches");
    check(uniform.transferResult().diagnostics().quadraturePointCount
              == bridge.overlapPatchCount(),
          "face bridge: every overlap becomes one conservative quadrature patch");

    const auto disturbed = fluidPistonDiagnostics(true);
    check(disturbed.surfaces[1]
              .maximumPressureTractionDeviationPascals > 0.0,
          "face bridge: source regression contains genuinely nonuniform pressure");
    const auto first = bridge.evaluate(disturbed, kinematics);
    const auto second = bridge.evaluate(disturbed, kinematics);
    const auto& diagnostics = first.diagnostics();
    check(first == second,
          "face bridge: nonuniform face correspondence replays bit-for-bit");
    check(diagnostics.forceResidualNormNewtons < 2.0e-10
              && diagnostics.momentResidualNormNewtonMeters < 2.0e-10
              && std::abs(diagnostics.powerResidualWatts) < 2.0e-10
              && diagnostics.maximumFacePowerResidualWatts < 2.0e-10
              && diagnostics.finite,
          "face bridge: nonuniform force, moment, global power, and per-face power close");
    checkVectorNear(diagnostics.structureSurfaceForceNewtons,
                    diagnostics.fluidPressureForceNewtons, 2.0e-10,
                    "face bridge: structural quadrature preserves fluid face force");
    checkVectorNear(diagnostics.structureSurfaceMomentNewtonMeters,
                    diagnostics.fluidPressureMomentNewtonMeters, 2.0e-10,
                    "face bridge: overlap centroids preserve fluid face moment");
    checkNear(diagnostics.structureSurfacePowerWatts,
              diagnostics.fluidPressurePowerWatts, 2.0e-10,
              "face bridge: barycentric velocity preserves fluid face power");

    auto changedGeometry = disturbed;
    changedGeometry.faces[1].lowerCornerMeters.x += 0.01;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            changedGeometry, kinematics)); },
        "face bridge validation: changed grid-face geometry is rejected");

    auto incompleteKinematics = kinematics;
    incompleteKinematics.pop_back();
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            disturbed, incompleteKinematics)); },
        "face bridge validation: incomplete kinematics are rejected before mapping");

    auto invalidFaceVelocity = disturbed;
    invalidFaceVelocity.faces[1].normalVelocityMetersPerSecond =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            invalidFaceVelocity, kinematics)); },
        "face bridge validation: non-finite face velocity is rejected");

    auto aggregateMismatch = disturbed;
    aggregateMismatch.surfaces[1].pressureForceNewtons.x += 1.0;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            aggregateMismatch, kinematics)); },
        "face bridge validation: aggregate force must equal its face ledger");

    auto cancellingPowerMismatch = disturbed;
    auto firstRight = std::ranges::find_if(
        cancellingPowerMismatch.faces,
        [](const auto& face) { return face.surfaceStableId == 200; });
    auto secondRight = std::ranges::find_if(
        std::next(firstRight), cancellingPowerMismatch.faces.end(),
        [](const auto& face) { return face.surfaceStableId == 200; });
    const double firstForce = firstRight->pressureForceNewtons.x;
    const double secondForce = secondRight->pressureForceNewtons.x;
    const double firstVelocityChange = 0.1;
    const double secondVelocityChange =
        -firstVelocityChange * firstForce / secondForce;
    firstRight->normalVelocityMetersPerSecond += firstVelocityChange;
    firstRight->pressurePowerWatts += firstForce * firstVelocityChange;
    secondRight->normalVelocityMetersPerSecond += secondVelocityChange;
    secondRight->pressurePowerWatts += secondForce * secondVelocityChange;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            cancellingPowerMismatch, kinematics)); },
        "face bridge validation: cancelling local power mismatch is rejected");

    auto incompleteReferenceFaces = reference.faces;
    const auto omitted = std::ranges::find_if(
        incompleteReferenceFaces,
        [](const auto& face) { return face.surfaceStableId == 200; });
    incompleteReferenceFaces.erase(omitted);
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge incomplete(
            structure, 200, pistonNodes(), pistonTriangles(),
            incompleteReferenceFaces); },
        "face bridge validation: incomplete reference coverage is rejected");

    auto reversedReferenceFaces = reference.faces;
    std::reverse(reversedReferenceFaces.begin(), reversedReferenceFaces.end());
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge reordered(
            structure, 200, pistonNodes(), pistonTriangles(),
            reversedReferenceFaces); },
        "face bridge validation: noncanonical reference face order is rejected");

    simwing::fsi::PlanarFaceResolvedBridgeSettings invalidSettings;
    invalidSettings.transfer.minimumQuadratureAreaSquareMeters = 1.0;
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge invalidConfiguration(
            structure, 200, pistonNodes(), pistonTriangles(),
            reference.faces, invalidSettings); },
        "face bridge validation: incompatible overlap and quadrature thresholds are rejected");

    simwing::fsi::PlanarFaceResolvedBridgeSettings invalidMode;
    invalidMode.correspondenceMode =
        static_cast<PlanarFaceCorrespondenceMode>(255);
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge invalidCorrespondenceMode(
            structure, 200, pistonNodes(), pistonTriangles(),
            reference.faces, invalidMode); },
        "face bridge validation: unknown correspondence modes are rejected");

    auto overlappingReferenceFaces = reference.faces;
    auto firstReferenceRight = std::ranges::find_if(
        overlappingReferenceFaces,
        [](const auto& face) { return face.surfaceStableId == 200; });
    auto secondReferenceRight = std::ranges::find_if(
        std::next(firstReferenceRight), overlappingReferenceFaces.end(),
        [](const auto& face) { return face.surfaceStableId == 200; });
    secondReferenceRight->lowerCornerMeters =
        firstReferenceRight->lowerCornerMeters;
    secondReferenceRight->upperCornerMeters =
        firstReferenceRight->upperCornerMeters;
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge overlappingFaces(
            structure, 200, pistonNodes(), pistonTriangles(),
            overlappingReferenceFaces); },
        "face bridge validation: overlapping fluid tiles are rejected");

    StructureDefinition reversedDefinition = pistonDefinition();
    reversedDefinition.triangles = {
        {{0, 2, 1}}, {{0, 3, 2}},
    };
    Structure reversedStructure(reversedDefinition);
    const std::vector<CouplingSurfaceTriangleDefinition> reversedTriangles{
        {1000, {10, 30, 20}},
        {2000, {10, 40, 30}},
    };
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge reversedOrientation(
            reversedStructure, 200, pistonNodes(), reversedTriangles,
            reference.faces); },
        "face bridge validation: reversed structural orientation is rejected");

    StructureDefinition overlappingDefinition;
    overlappingDefinition.nodes = {
        {{3.0, 0.0, 0.0}, 1.0, false},
        {{3.0, 2.0, 0.0}, 1.0, false},
        {{3.0, 2.0, 3.0}, 1.0, false},
        {{3.0, 0.0, 0.0}, 1.0, false},
        {{3.0, 2.0, 0.0}, 1.0, false},
        {{3.0, 2.0, 3.0}, 1.0, false},
    };
    overlappingDefinition.triangles = {
        {{0, 1, 2}}, {{3, 4, 5}},
    };
    Structure overlappingStructure(overlappingDefinition);
    const std::vector<CouplingSurfaceNodeDefinition> overlappingNodes{
        {10, 0}, {20, 1}, {30, 2},
        {40, 3}, {50, 4}, {60, 5},
    };
    const std::vector<CouplingSurfaceTriangleDefinition>
        overlappingTriangles{
            {1000, {10, 20, 30}},
            {2000, {40, 50, 60}},
        };
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge overlappingStructureFaces(
            overlappingStructure, 200, overlappingNodes,
            overlappingTriangles, reference.faces); },
        "face bridge validation: overlapping structural triangles are rejected");
}

void testMovingPlanarFaceCorrespondence() {
    Structure structure(pistonDefinition());
    const auto reference = openPistonDiagnostics();
    PlanarFaceResolvedBridgeSettings settings;
    settings.correspondenceMode =
        PlanarFaceCorrespondenceMode::RigidNormalTranslation;
    PlanarFaceResolvedFluidStructureBridge bridge(
        structure, 300, pistonNodes(), pistonTriangles(),
        reference.faces, settings);

    const auto translatedKinematics = movingPlaneKinematics(
        bridge.transfer(), structure, 3.2);
    const auto first = bridge.evaluateMovingPlane(
        reference, translatedKinematics, 3.2);
    const auto second = bridge.evaluateMovingPlane(
        reference, translatedKinematics, 3.2);
    const auto& diagnostics = first.diagnostics();
    check(first == second,
          "moving bridge: translated correspondence replays bit-for-bit");
    check(diagnostics.correspondenceMode
              == PlanarFaceCorrespondenceMode::RigidNormalTranslation
              && diagnostics.gridPlaneCoordinateMeters == 3.0
              && diagnostics.physicalPlaneCoordinateMeters == 3.2,
          "moving bridge: Eulerian grid plane and physical plane remain explicit");
    checkNear(diagnostics.normalTranslationFromReferenceMeters,
              0.2, 2.0e-16,
              "moving bridge: physical translation is measured from reference");
    checkNear(diagnostics.maximumRigidPositionResidualMeters,
              0.0, 0.0,
              "moving bridge: rigid translated nodes match the physical plane");
    checkNear(diagnostics.maximumRigidVelocityResidualMetersPerSecond,
              0.0, 0.0,
              "moving bridge: structural and fluid normal velocities match");
    check(diagnostics.forceResidualNormNewtons < 2.0e-10
              && diagnostics.momentResidualNormNewtonMeters < 2.0e-10
              && std::abs(diagnostics.powerResidualWatts) < 2.0e-10
              && diagnostics.maximumFacePowerResidualWatts < 2.0e-10,
          "moving bridge: translated force, moment, and power ledgers close");

    const auto grid = pistonGrid();
    const auto referenceInterfaces = openPlaneInterfaces(
        GridFaceAxis::X, 6);
    const PlanarMovingControlVolume referenceControlVolume(
        grid, referenceInterfaces, 300, 2);
    const auto cutSurface = evaluatePlanarCutSurfacePressure(
        grid, referenceControlVolume, reference, 0.2, 3.2);
    const auto cutTransfer = bridge.evaluateCutSurface(
        cutSurface, translatedKinematics);
    check(cutSurface.accepted && cutTransfer != first
              && cutTransfer.diagnostics().fluidPressureForceNewtons.x
                  == cutSurface.pressureForceNewtons.x,
          "moving bridge: cut geometry transfers the complete constraint reaction");
    auto corruptedCutGeometry = cutSurface;
    corruptedCutGeometry.faces.front().physicalLowerCornerMeters.y += 0.01;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateCutSurface(
            corruptedCutGeometry, translatedKinematics)); },
        "moving bridge validation: cut geometry is rebound before transfer");
    auto unacceptedCutSurface = cutSurface;
    unacceptedCutSurface.accepted = false;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateCutSurface(
            unacceptedCutSurface, translatedKinematics)); },
        "moving bridge validation: unaccepted cut geometry is rejected");
    auto nonfiniteCutLedger = cutSurface;
    nonfiniteCutLedger.pressureMomentNewtonMeters.x =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateCutSurface(
            nonfiniteCutLedger, translatedKinematics)); },
        "moving bridge validation: cut ledgers are rebound for finiteness");

    const auto rebasedFluid = openPistonDiagnostics(7);
    const auto rebasedKinematics = movingPlaneKinematics(
        bridge.transfer(), structure, 3.6);
    const auto rebased = bridge.evaluateMovingPlane(
        rebasedFluid, rebasedKinematics, 3.6);
    check(rebased.diagnostics().gridPlaneCoordinateMeters == 3.5
              && rebased.diagnostics().physicalPlaneCoordinateMeters == 3.6
              && rebased.diagnostics().overlapPatchCount
                  == bridge.overlapPatchCount(),
          "moving bridge: one material patch map survives a grid-plane rebase");
    checkNear(rebased.diagnostics().normalTranslationFromReferenceMeters,
              0.6, 5.0e-16,
              "moving bridge: translation remains unwrapped after rebase");

    expectRejected(
        [&] { static_cast<void>(bridge.evaluate(
            reference, translatedKinematics)); },
        "moving bridge validation: correspondence mode must be explicit");
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            reference, translatedKinematics, 3.19)); },
        "moving bridge validation: declared physical plane must match structure");
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            reference, translatedKinematics,
            std::numeric_limits<double>::quiet_NaN())); },
        "moving bridge validation: the declared physical plane must be finite");
    auto nonrigidKinematics = translatedKinematics;
    nonrigidKinematics.front().positionMeters.y += 0.01;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            reference, nonrigidKinematics, 3.2)); },
        "moving bridge validation: transverse structural motion is rejected");
    const auto wrongVelocityKinematics = movingPlaneKinematics(
        bridge.transfer(), structure, 3.2, 0.2);
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            reference, wrongVelocityKinematics, 3.2)); },
        "moving bridge validation: fluid and structural speeds must match");
    auto inconsistentZeroPressure = openPistonDiagnostics(6, 0.0);
    auto stationaryKinematics = movingPlaneKinematics(
        bridge.transfer(), structure, 3.2, 0.0);
    static_cast<void>(bridge.evaluateMovingPlane(
        inconsistentZeroPressure, stationaryKinematics, 3.2));
    inconsistentZeroPressure.faces.front().normalVelocityMetersPerSecond =
        0.01;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            inconsistentZeroPressure, stationaryKinematics, 3.2)); },
        "moving bridge validation: velocity correspondence is checked at zero pressure");
    auto changedTransverseTile = reference;
    changedTransverseTile.faces.front().lowerCornerMeters.y += 0.01;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluateMovingPlane(
            changedTransverseTile, translatedKinematics, 3.2)); },
        "moving bridge validation: transverse fluid tiling is immutable");

    PlanarFaceResolvedBridgeSettings fixedSettings;
    expectRejected(
        [&] { PlanarFaceResolvedFluidStructureBridge fixedNonseparating(
            structure, 300, pistonNodes(), pistonTriangles(),
            reference.faces, fixedSettings); },
        "moving bridge validation: nonseparating surfaces require moving mode");
}

void testMovingPlanarCorrespondenceAllAxes() {
    struct AxisCase {
        GridFaceAxis axis;
        std::size_t referencePlane;
        std::size_t rebasedPlane;
        double referenceCoordinateMeters;
        double rebasedGridCoordinateMeters;
        double rebasedPhysicalCoordinateMeters;
    };
    const std::array cases{
        AxisCase{GridFaceAxis::X, 6, 7, 3.0, 3.5, 3.7},
        AxisCase{GridFaceAxis::Y, 1, 0, 1.0, 0.0, 2.2},
        AxisCase{GridFaceAxis::Z, 2, 0, 2.0, 0.0, 3.2},
    };
    for (const auto& test : cases) {
        try {
            Structure structure(movingPlaneDefinition(
                test.axis, test.referenceCoordinateMeters));
            const auto reference = openPlaneDiagnostics(
                test.axis, test.referencePlane);
            PlanarFaceResolvedBridgeSettings settings;
            settings.correspondenceMode =
                PlanarFaceCorrespondenceMode::RigidNormalTranslation;
            PlanarFaceResolvedFluidStructureBridge bridge(
                structure, 300, pistonNodes(), pistonTriangles(),
                reference.faces, settings);
            const double firstPhysical =
                test.referenceCoordinateMeters + 0.2;
            const auto firstKinematics = movingPlaneKinematics(
                bridge.transfer(), structure, firstPhysical,
                0.25, test.axis);
            const auto first = bridge.evaluateMovingPlane(
                reference, firstKinematics, firstPhysical);
            check(first.diagnostics().correspondenceMode
                      == PlanarFaceCorrespondenceMode::RigidNormalTranslation
                      && first.diagnostics()
                             .maximumRigidPositionResidualMeters == 0.0
                      && first.diagnostics()
                             .maximumRigidVelocityResidualMetersPerSecond
                          == 0.0,
                  "moving axes: rigid correspondence accepts every orientation");

            const auto rebasedFluid = openPlaneDiagnostics(
                test.axis, test.rebasedPlane);
            const auto rebasedKinematics = movingPlaneKinematics(
                bridge.transfer(), structure,
                test.rebasedPhysicalCoordinateMeters,
                0.25, test.axis);
            const auto rebased = bridge.evaluateMovingPlane(
                rebasedFluid, rebasedKinematics,
                test.rebasedPhysicalCoordinateMeters);
            checkNear(rebased.diagnostics().gridPlaneCoordinateMeters,
                      test.rebasedGridCoordinateMeters, 0.0,
                      "moving axes: rebased Eulerian plane is orientation independent");
            checkNear(rebased.diagnostics().physicalPlaneCoordinateMeters,
                      test.rebasedPhysicalCoordinateMeters, 0.0,
                      "moving axes: unwrapped physical plane survives periodic rebase");
            check(rebased.diagnostics().forceResidualNormNewtons < 2.0e-10
                      && rebased.diagnostics()
                             .momentResidualNormNewtonMeters < 2.0e-10
                      && std::abs(
                             rebased.diagnostics().powerResidualWatts)
                          < 2.0e-10,
                  "moving axes: rebased force, moment, and power ledgers close");
        } catch (const std::exception& exception) {
            std::fprintf(
                stderr,
                "FAIL: moving axes: unexpected rejection for axis %d: %s\n",
                static_cast<int>(test.axis), exception.what());
            ++failures;
        }
    }
}

void testPorousFaceResolvedTransfer() {
    const auto traction = porousSurfaceTraction();
    Structure structure(movingPlaneDefinition(GridFaceAxis::X, 1.45));
    PlanarFaceResolvedFluidStructureBridge bridge(
        structure, 300, pistonNodes(), pistonTriangles(),
        traction.faces);
    const auto kinematics = translatingKinematics(
        bridge.transfer(), structure, 0.1);
    const auto first = bridge.evaluatePorousSurface(
        traction, kinematics);
    const auto second = bridge.evaluatePorousSurface(
        traction, kinematics);
    const auto& diagnostics = first.diagnostics();

    check(first == second && diagnostics.accepted && diagnostics.finite
              && diagnostics.version
                  == simwing::fsi::porousFaceResolvedBridgeVersion
              && diagnostics.fluidSurfaceStableId == 300
              && diagnostics.mapping.fluidFaceCount == 6
              && diagnostics.mapping.overlapPatchCount > 0,
          "porous bridge: stable-ID face transfer replays deterministically");
    checkVectorNear(diagnostics.pressureForceOnFluidNewtons,
                    {-9.0, 0.0, 0.0}, 2.0e-10,
                    "porous bridge: fluid-side porous force remains explicit");
    checkVectorNear(diagnostics.pressureForceOnSurfaceNewtons,
                    {9.0, 0.0, 0.0}, 2.0e-10,
                    "porous bridge: equal-and-opposite sheet load is selected for transfer");
    checkVectorNear(
        diagnostics.mapping.structureSurfaceForceNewtons,
        {9.0, 0.0, 0.0}, 2.0e-10,
        "porous bridge: structural quadrature receives the sheet reaction");
    checkVectorNear(
        diagnostics.transferredSurfaceImpulseNewtonSeconds,
        {3.6, 0.0, 0.0}, 2.0e-10,
        "porous bridge: mapped sheet force retains its macro-step impulse");
    check(diagnostics.impulseResidualNormNewtonSeconds < 3.0e-15,
          "porous bridge: source and transferred impulses close");
    checkNear(diagnostics.pressurePowerToFluidWatts,
              -2.25, 2.0e-10,
              "porous bridge: fluid pressure power is not reassigned to structure");
    checkNear(diagnostics.pressurePowerToSurfaceWatts,
              0.9, 2.0e-10,
              "porous bridge: authored sheet power remains explicit");
    checkNear(diagnostics.mapping.structureSurfacePowerWatts,
              0.9, 2.0e-10,
              "porous bridge: structural kinematics reproduce sheet power");
    checkNear(diagnostics.transferredSurfaceWorkJoules,
              0.36, 2.0e-10,
              "porous bridge: structural pressure work integrates over the macro step");
    checkNear(diagnostics.workResidualJoules,
              0.0, 3.0e-16,
              "porous bridge: source and transferred work close");
    checkNear(diagnostics.porousDissipatedEnergyJoules,
              0.54, 2.0e-10,
              "porous bridge: material dissipation remains separate from structure work");
    checkNear(diagnostics.sourceEnergyResidualJoules,
              0.0, 3.0e-16,
              "porous bridge: fluid, sheet, and dissipative energy identity survives transfer");

    ConservativeMacroStepCoupling coupling(bridge.transfer());
    const std::array<double, 2> offsets{0.0, 0.4};
    const std::array samples{
        first.transferResult(), second.transferResult()};
    const auto integrated = coupling.integrate(offsets, samples);
    checkVectorNear(
        integrated.diagnostics().integratedSurfaceImpulseNewtonSeconds,
        diagnostics.pressureImpulseOnSurfaceNewtonSeconds,
        2.0e-10,
        "porous bridge: sheet reaction reaches temporal coupling as the same impulse");
    checkNear(integrated.diagnostics().integratedSurfaceWorkJoules,
              diagnostics.pressureWorkToSurfaceJoules, 2.0e-10,
              "porous bridge: sheet power reaches temporal coupling as the same work");

    PlanarFaceResolvedBridgeSettings movingSettings;
    movingSettings.correspondenceMode =
        PlanarFaceCorrespondenceMode::RigidNormalTranslation;
    PlanarFaceResolvedFluidStructureBridge movingBridge(
        structure, 300, pistonNodes(), pistonTriangles(),
        traction.faces, movingSettings);
    const auto movedTraction = porousSurfaceTraction(0.6);
    const auto movedKinematics = movingPlaneKinematics(
        movingBridge.transfer(), structure, 1.55, 0.1);
    const auto moved = movingBridge.evaluateMovingPorousSurface(
        movedTraction, movedKinematics, 1.55);
    check(moved.diagnostics().accepted
              && moved.diagnostics().mapping.correspondenceMode
                  == PlanarFaceCorrespondenceMode::RigidNormalTranslation,
          "porous bridge: rigid sheet motion uses the explicit moving correspondence");
    checkNear(
        moved.diagnostics().mapping.normalTranslationFromReferenceMeters,
        0.1, 2.0e-16,
        "porous bridge: subcell sheet translation remains geometric, not a load shortcut");
    checkVectorNear(
        moved.diagnostics().mapping.structureSurfaceForceNewtons,
        {9.0, 0.0, 0.0}, 2.0e-10,
        "porous bridge: translated sheet retains conservative force transfer");

    auto rejected = traction;
    rejected.accepted = false;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluatePorousSurface(
            rejected, kinematics)); },
        "porous bridge validation: unaccepted traction is rejected");
    auto corruptedImpulse = traction;
    corruptedImpulse.surfaces.front()
        .pressureImpulseOnSurfaceNewtonSeconds.x += 1.0;
    expectRejected(
        [&] { static_cast<void>(bridge.evaluatePorousSurface(
            corruptedImpulse, kinematics)); },
        "porous bridge validation: source impulse must match mapped force");
    auto wrongVelocity = kinematics;
    for (auto& node : wrongVelocity) {
        node.velocityMetersPerSecond.x = 0.2;
    }
    expectRejected(
        [&] { static_cast<void>(bridge.evaluatePorousSurface(
            traction, wrongVelocity)); },
        "porous bridge validation: structural velocity must match sheet power");
    expectRejected(
        [&] { static_cast<void>(movingBridge.evaluatePorousSurface(
            traction, kinematics)); },
        "porous bridge validation: moving correspondence requires an explicit physical plane");
}

} // namespace

int main() {
    testAnalyticStableIdBridgeAndStructuralAcceptance();
    testStrictRejectionContracts();
    testFaceResolvedNonuniformPressureTransfer();
    testMovingPlanarFaceCorrespondence();
    testMovingPlanarCorrespondenceAllAxes();
    testPorousFaceResolvedTransfer();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid-structure bridge check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid-structure bridge checks passed");
    return 0;
}
