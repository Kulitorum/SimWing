#include "coupling.h"
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
using simwing::fsi::UniformFluidStructureBridge;
using simwing::fsi::fluid::CellScalarField;
using simwing::fsi::fluid::FaceAlignedMovingInterface;
using simwing::fsi::fluid::GridFaceAxis;
using simwing::fsi::fluid::GridFaceMovingInterface;
using simwing::fsi::fluid::MacVelocityField;
using simwing::fsi::fluid::MovingInterfaceProjectionDiagnostics;
using simwing::fsi::fluid::MovingInterfaceProjectionSettings;
using simwing::fsi::fluid::PeriodicCartesianGrid;
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

MovingInterfaceProjectionDiagnostics fluidPistonDiagnostics() {
    const auto grid = pistonGrid();
    const FaceAlignedMovingInterface interfaces(grid, slabFaces(grid));
    MacVelocityField velocity(grid);
    std::ranges::fill(velocity.xFaces(), 0.25);
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

} // namespace

int main() {
    testAnalyticStableIdBridgeAndStructuralAcceptance();
    testStrictRejectionContracts();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing fluid-structure bridge check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing fluid-structure bridge checks passed");
    return 0;
}
