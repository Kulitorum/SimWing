#include "transfer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using simwing::fsi::ConservativeSurfaceTransfer;
using simwing::fsi::ConservativeTransferSettings;
using simwing::fsi::CouplingNodeKinematics;
using simwing::fsi::CouplingSurfaceNodeDefinition;
using simwing::fsi::CouplingSurfaceTriangleDefinition;
using simwing::fsi::CouplingTriangleTraction;
using simwing::fsi::CouplingTriangleTractionQuadrature;
using simwing::fsi::Structure;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureStepSettings;
using simwing::fsi::StructureVector3;

int failures = 0;

template<typename Callback>
void expectRejected(Callback&& callback, const char* message);

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
        std::fprintf(
            stderr,
            "FAIL: %s (actual [%.17g %.17g %.17g], expected [%.17g %.17g %.17g])\n",
            message, actual.x, actual.y, actual.z,
            expected.x, expected.y, expected.z);
        ++failures;
    }
}

StructureDefinition rectangleDefinition(const double zOffset = 0.0) {
    StructureDefinition definition;
    definition.nodes = {
        {{1.0, -1.0, zOffset}, 1.0, false},
        {{3.0, -1.0, zOffset}, 1.0, false},
        {{3.0, 2.0, zOffset}, 1.0, false},
        {{1.0, 2.0, zOffset}, 1.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

std::vector<CouplingSurfaceNodeDefinition> scrambledNodes(
    const std::uint64_t offset = 0) {
    return {
        {40 + offset, 3},
        {10 + offset, 0},
        {30 + offset, 2},
        {20 + offset, 1},
    };
}

std::vector<CouplingSurfaceTriangleDefinition> scrambledTriangles(
    const std::uint64_t offset = 0) {
    return {
        {200 + offset, {10 + offset, 30 + offset, 40 + offset}},
        {100 + offset, {10 + offset, 20 + offset, 30 + offset}},
    };
}

std::vector<CouplingTriangleTraction> uniformTractions(
    const ConservativeSurfaceTransfer& transfer,
    const StructureVector3 traction = {2.0, -3.0, 5.0}) {
    std::vector<CouplingTriangleTraction> result;
    for (const auto& triangle : transfer.triangles()) {
        result.push_back({triangle.stableId, traction});
    }
    return result;
}

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

StructureVector3 cross(const StructureVector3& first,
                       const StructureVector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double dot(const StructureVector3& first, const StructureVector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

std::vector<CouplingNodeKinematics> rigidKinematics(
    const ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const StructureVector3 reference,
    const StructureVector3 translationVelocity,
    const StructureVector3 angularVelocity) {
    std::vector<CouplingNodeKinematics> result;
    for (const auto& node : transfer.nodes()) {
        const auto position =
            structure.definition().nodes[node.structureNode].positionMeters;
        result.push_back({
            node.stableId,
            position,
            add(translationVelocity,
                cross(angularVelocity, subtract(position, reference))),
        });
    }
    return result;
}

void testCanonicalTopologyAndCapture() {
    Structure structure(rectangleDefinition());
    auto nodes = scrambledNodes();
    auto triangles = scrambledTriangles();
    ConservativeSurfaceTransfer transfer(structure, nodes, triangles);
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(triangles.begin(), triangles.end());
    for (auto& triangle : triangles) {
        triangle.nodeStableIds = {
            triangle.nodeStableIds[1],
            triangle.nodeStableIds[2],
            triangle.nodeStableIds[0],
        };
    }
    ConservativeSurfaceTransfer reordered(structure, nodes, triangles);

    check(transfer.fingerprint() == reordered.fingerprint(),
          "topology: authored order and cyclic corners do not change the fingerprint");
    check(std::ranges::equal(transfer.nodes(), reordered.nodes())
              && std::ranges::equal(
                  transfer.triangles(), reordered.triangles()),
          "topology: node and triangle definitions are canonicalized by stable ID");
    check(transfer.nodes()[0].stableId == 10
              && transfer.nodes()[3].stableId == 40
              && transfer.triangles()[0].stableId == 100
              && transfer.triangles()[1].stableId == 200,
          "topology: public canonical order is explicit");
    const auto captured = transfer.captureKinematics(structure);
    check(captured.size() == 4 && captured[0].stableId == 10,
          "capture: stable node order matches the coupling surface");
    checkVectorNear(captured[0].positionMeters, {1.0, -1.0, 0.0}, 0.0,
                    "capture: structural positions cross the adapter exactly");
    checkVectorNear(captured[0].velocityMetersPerSecond, {}, 0.0,
                    "capture: structural velocities cross the adapter exactly");
}

void testAnalyticForceMomentAndTranslationPower() {
    Structure structure(rectangleDefinition());
    ConservativeSurfaceTransfer transfer(
        structure, scrambledNodes(), scrambledTriangles());
    const StructureVector3 reference{0.5, -0.5, 1.0};
    const StructureVector3 translation{0.4, -0.2, 0.1};
    const auto kinematics = rigidKinematics(
        transfer, structure, reference, translation, {});
    const auto tractions = uniformTractions(transfer);
    ConservativeTransferSettings settings;
    settings.momentReferenceMeters = reference;
    const auto result = transfer.evaluate(kinematics, tractions, settings);
    const auto& diagnostics = result.diagnostics();

    checkNear(diagnostics.surfaceAreaSquareMeters, 6.0, 1.0e-15,
              "translation: current triangle areas integrate to the rectangle area");
    checkVectorNear(diagnostics.integratedSurfaceForceNewtons,
                    {12.0, -18.0, 30.0}, 1.0e-14,
                    "translation: traction times area gives the analytic force");
    checkVectorNear(diagnostics.integratedSurfaceMomentNewtonMeters,
                    {12.0, -57.0, -39.0}, 2.0e-14,
                    "translation: centroid traction gives the analytic moment");
    checkNear(diagnostics.integratedSurfacePowerWatts, 11.4, 5.0e-15,
              "translation: surface power is resultant force dot rigid velocity");
    checkVectorNear(diagnostics.transferredNodalForceNewtons,
                    diagnostics.integratedSurfaceForceNewtons, 2.0e-14,
                    "translation: nodal and surface force ledgers agree");
    checkVectorNear(diagnostics.transferredNodalMomentNewtonMeters,
                    diagnostics.integratedSurfaceMomentNewtonMeters, 2.0e-14,
                    "translation: nodal and surface moment ledgers agree");
    checkNear(diagnostics.transferredNodalPowerWatts,
              diagnostics.integratedSurfacePowerWatts, 1.0e-14,
              "translation: nodal and surface power ledgers agree");
    check(diagnostics.forceResidualNormNewtons < 2.0e-14
              && diagnostics.momentResidualNormNewtonMeters < 3.0e-14
              && std::abs(diagnostics.powerResidualWatts) < 1.0e-14
              && diagnostics.finite,
          "translation: all conservative residuals meet their budgets");

    const auto loads = result.nodeLoads();
    check(loads.size() == 4
              && loads[0].stableId == 10
              && loads[0].structureNode == 0,
          "translation: nodal loads preserve canonical stable and target IDs");
    checkVectorNear(loads[0].forceNewtons, {4.0, -6.0, 10.0}, 1.0e-15,
                    "translation: shared diagonal vertex receives two barycentric shares");
    checkVectorNear(loads[1].forceNewtons, {2.0, -3.0, 5.0}, 1.0e-15,
                    "translation: boundary vertex receives one barycentric share");
    checkVectorNear(loads[2].forceNewtons, {4.0, -6.0, 10.0}, 1.0e-15,
                    "translation: second shared vertex receives two shares");
    checkVectorNear(loads[3].forceNewtons, {2.0, -3.0, 5.0}, 1.0e-15,
                    "translation: final boundary vertex receives one share");
}

void testRigidRotationPowerIdentityAndDeterminism() {
    Structure structure(rectangleDefinition());
    auto nodes = scrambledNodes();
    auto triangles = scrambledTriangles();
    ConservativeSurfaceTransfer first(structure, nodes, triangles);
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(triangles.begin(), triangles.end());
    ConservativeSurfaceTransfer second(structure, nodes, triangles);
    const StructureVector3 reference{0.5, -0.5, 1.0};
    const StructureVector3 translation{0.4, -0.2, 0.1};
    const StructureVector3 angular{0.3, -0.4, 0.2};
    const auto kinematics = rigidKinematics(
        first, structure, reference, translation, angular);
    const std::vector<CouplingTriangleTraction> tractions = {
        {100, {2.0, -3.0, 5.0}},
        {200, {-1.0, 4.0, 2.0}},
    };
    ConservativeTransferSettings settings;
    settings.momentReferenceMeters = reference;
    const auto firstResult = first.evaluate(kinematics, tractions, settings);
    const auto secondResult = second.evaluate(kinematics, tractions, settings);
    check(firstResult == secondResult,
          "rotation: canonicalized surfaces replay transfer bit-for-bit");

    const auto& diagnostics = firstResult.diagnostics();
    checkVectorNear(diagnostics.integratedSurfaceForceNewtons,
                    {3.0, 3.0, 21.0}, 1.0e-14,
                    "rotation: distinct triangle tractions give the analytic force");
    checkVectorNear(diagnostics.integratedSurfaceMomentNewtonMeters,
                    {19.5, -37.5, -1.0}, 2.0e-14,
                    "rotation: distinct triangle tractions give the analytic moment");
    const double analyticRigidPower =
        dot(diagnostics.integratedSurfaceForceNewtons, translation)
        + dot(diagnostics.integratedSurfaceMomentNewtonMeters, angular);
    checkNear(analyticRigidPower, 23.35, 1.0e-14,
              "rotation: analytic wrench identity has the expected value");
    checkNear(diagnostics.integratedSurfacePowerWatts,
              analyticRigidPower, 1.0e-14,
              "rotation: integrated power equals F dot V plus M dot omega");
    checkNear(diagnostics.transferredNodalPowerWatts,
              analyticRigidPower, 2.0e-14,
              "rotation: barycentric nodal power preserves rigid rotational work");
    check(std::abs(diagnostics.powerResidualWatts) < 2.0e-14,
          "rotation: independent surface and nodal power ledgers close");
}

void testLoadApplicationAndBinding() {
    Structure structure(rectangleDefinition());
    ConservativeSurfaceTransfer transfer(
        structure, scrambledNodes(), scrambledTriangles());
    const auto kinematics = transfer.captureKinematics(structure);
    const auto result = transfer.evaluate(
        kinematics, uniformTractions(transfer));
    structure.addExternalForce(0, {1.0, 0.0, 0.0});
    transfer.addLoadsTo(structure, result);
    const auto checkpoint = structure.checkpoint();
    checkVectorNear(checkpoint.pendingExternalForcesNewtons[0],
                    {5.0, -6.0, 10.0}, 1.0e-15,
                    "application: transfer adds to an existing nodal load");
    checkVectorNear(checkpoint.pendingExternalForcesNewtons[1],
                    {2.0, -3.0, 5.0}, 1.0e-15,
                    "application: canonical target index receives its load");
    checkVectorNear(structure.diagnostics().pendingExternalForceNewtons,
                    {13.0, -18.0, 30.0}, 2.0e-14,
                    "application: structural pending-force ledger includes the transfer");
    StructureStepSettings stepSettings;
    stepSettings.timeStepSeconds = 0.5;
    stepSettings.substeps = 1;
    stepSettings.constraintIterations = 0;
    stepSettings.gravityMetersPerSecondSquared = {};
    stepSettings.velocityDampingPerSecond = 0.0;
    const auto stepDiagnostics = structure.step(stepSettings);
    const auto steppedNodes = structure.nodeStates();
    checkVectorNear(stepDiagnostics.lastAppliedExternalForceNewtons,
                    {13.0, -18.0, 30.0}, 2.0e-14,
                    "application: accepted XPBD step reports the transferred resultant");
    checkVectorNear(steppedNodes[0].velocityMetersPerSecond,
                    {2.5, -3.0, 5.0}, 2.0e-14,
                    "application: transferred node load reaches structural integration");
    checkVectorNear(steppedNodes[1].velocityMetersPerSecond,
                    {1.0, -1.5, 2.5}, 2.0e-14,
                    "application: barycentric load follows the actual nodal path");

    Structure foreign(rectangleDefinition(0.25));
    const auto beforeForeign = foreign.checkpoint();
    bool rejected = false;
    try {
        transfer.addLoadsTo(foreign, result);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected
              && foreign.checkpoint().pendingExternalForcesNewtons
                  == beforeForeign.pendingExternalForcesNewtons,
          "binding: foreign structural definition is rejected before mutation");

    ConservativeSurfaceTransfer otherSurface(
        structure, scrambledNodes(1000), scrambledTriangles(1000));
    const auto otherResult = otherSurface.evaluate(
        rigidKinematics(otherSurface, structure, {}, {}, {}),
        uniformTractions(otherSurface));
    const auto beforeWrongSurface = structure.checkpoint();
    rejected = false;
    try {
        transfer.addLoadsTo(structure, otherResult);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected
              && structure.checkpoint().pendingExternalForcesNewtons
                  == beforeWrongSurface.pendingExternalForcesNewtons,
          "binding: result from another surface is rejected before mutation");
}

void testNonuniformBarycentricQuadrature() {
    Structure structure(rectangleDefinition());
    ConservativeSurfaceTransfer transfer(
        structure, scrambledNodes(), scrambledTriangles());
    auto kinematics = transfer.captureKinematics(structure);
    for (auto& node : kinematics) {
        node.velocityMetersPerSecond = {0.25, -0.5, 0.0};
    }
    const std::vector<CouplingTriangleTractionQuadrature> quadrature = {
        {1, 100, {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
         1.0, {2.0, 0.0, 0.0}},
        {2, 100, {1.0 / 6.0, 5.0 / 12.0, 5.0 / 12.0},
         2.0, {4.0, 0.0, 0.0}},
        {3, 200, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
         3.0, {0.0, 2.0, 0.0}},
    };
    const auto first = transfer.evaluateQuadrature(kinematics, quadrature);
    const auto second = transfer.evaluateQuadrature(kinematics, quadrature);
    const auto& diagnostics = first.diagnostics();
    check(first == second,
          "quadrature: identical nonuniform patches replay bit-for-bit");
    check(diagnostics.triangleCount == 2
              && diagnostics.quadraturePointCount == 3,
          "quadrature: triangle topology and patch count remain distinct");
    checkNear(diagnostics.surfaceAreaSquareMeters, 6.0, 1.0e-15,
              "quadrature: patch areas recover the analytic surface area");
    checkVectorNear(diagnostics.integratedSurfaceForceNewtons,
                    {10.0, 6.0, 0.0}, 1.0e-14,
                    "quadrature: nonuniform patches integrate analytic force");
    checkVectorNear(diagnostics.integratedSurfaceMomentNewtonMeters,
                    {0.0, 0.0, 9.0}, 2.0e-14,
                    "quadrature: patch centroids integrate analytic moment");
    checkNear(diagnostics.integratedSurfacePowerWatts, -0.5, 1.0e-14,
              "quadrature: patch work uses interpolated linear velocity");
    check(diagnostics.forceResidualNormNewtons < 2.0e-14
              && diagnostics.momentResidualNormNewtonMeters < 3.0e-14
              && std::abs(diagnostics.powerResidualWatts) < 2.0e-14,
          "quadrature: independent nodal force, moment, and power ledgers close");
    const auto loads = first.nodeLoads();
    checkVectorNear(loads[0].forceNewtons,
                    {8.0 / 3.0, 2.0, 0.0}, 2.0e-15,
                    "quadrature: first shared node receives geometric patch weights");
    checkVectorNear(loads[1].forceNewtons,
                    {11.0 / 3.0, 0.0, 0.0}, 2.0e-15,
                    "quadrature: first boundary node receives nonuniform load");
    checkVectorNear(loads[2].forceNewtons,
                    {11.0 / 3.0, 2.0, 0.0}, 2.0e-15,
                    "quadrature: second shared node receives both triangles");
    checkVectorNear(loads[3].forceNewtons,
                    {0.0, 2.0, 0.0}, 2.0e-15,
                    "quadrature: final boundary node receives its overlap load");

    auto invalid = quadrature;
    std::swap(invalid[0], invalid[1]);
    expectRejected(
        [&] { static_cast<void>(transfer.evaluateQuadrature(
            kinematics, invalid)); },
        "quadrature validation: noncanonical patch order is rejected");
    invalid = quadrature;
    invalid[0].barycentricCoordinates = {-0.1, 0.5, 0.6};
    expectRejected(
        [&] { static_cast<void>(transfer.evaluateQuadrature(
            kinematics, invalid)); },
        "quadrature validation: patch outside its triangle is rejected");
    invalid = quadrature;
    invalid[0].areaSquareMeters = 0.0;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluateQuadrature(
            kinematics, invalid)); },
        "quadrature validation: zero-area patch is rejected");
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

void testTopologyAndEvaluationValidation() {
    Structure structure(rectangleDefinition());
    auto nodes = scrambledNodes();
    auto triangles = scrambledTriangles();
    auto invalidNodes = nodes;
    invalidNodes[0].stableId = 0;
    expectRejected(
        [&] { ConservativeSurfaceTransfer invalid(
            structure, invalidNodes, triangles); },
        "validation: zero coupling node stable ID is rejected");
    invalidNodes = nodes;
    invalidNodes[0].structureNode = invalidNodes[1].structureNode;
    expectRejected(
        [&] { ConservativeSurfaceTransfer invalid(
            structure, invalidNodes, triangles); },
        "validation: duplicate structural node mapping is rejected");
    auto invalidTriangles = triangles;
    invalidTriangles[0].nodeStableIds = {10, 40, 30};
    expectRejected(
        [&] { ConservativeSurfaceTransfer invalid(
            structure, nodes, invalidTriangles); },
        "validation: reversed structural triangle orientation is rejected");
    invalidTriangles = triangles;
    invalidTriangles[0].nodeStableIds[0] = 999;
    expectRejected(
        [&] { ConservativeSurfaceTransfer invalid(
            structure, nodes, invalidTriangles); },
        "validation: unknown triangle node stable ID is rejected");

    ConservativeSurfaceTransfer transfer(structure, nodes, triangles);
    auto kinematics = transfer.captureKinematics(structure);
    auto tractions = uniformTractions(transfer);
    std::reverse(kinematics.begin(), kinematics.end());
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(kinematics, tractions)); },
        "validation: reordered kinematics are rejected explicitly");
    kinematics = transfer.captureKinematics(structure);
    std::reverse(tractions.begin(), tractions.end());
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(kinematics, tractions)); },
        "validation: reordered tractions are rejected explicitly");
    tractions = uniformTractions(transfer);
    kinematics[2].positionMeters = kinematics[1].positionMeters;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(kinematics, tractions)); },
        "validation: degenerate current triangle is rejected");
    kinematics = transfer.captureKinematics(structure);
    tractions[0].tractionPascals.x =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(kinematics, tractions)); },
        "validation: non-finite traction is rejected");
    tractions = uniformTractions(transfer);
    ConservativeTransferSettings invalidSettings;
    invalidSettings.minimumTriangleAreaSquareMeters = 0.0;
    expectRejected(
        [&] { static_cast<void>(transfer.evaluate(
            kinematics, tractions, invalidSettings)); },
        "validation: non-positive minimum triangle area is rejected");

    Structure foreign(rectangleDefinition(0.5));
    expectRejected(
        [&] { static_cast<void>(transfer.captureKinematics(foreign)); },
        "validation: kinematics capture rejects a foreign structure definition");
}

} // namespace

int main() {
    testCanonicalTopologyAndCapture();
    testAnalyticForceMomentAndTranslationPower();
    testRigidRotationPowerIdentityAndDeterminism();
    testLoadApplicationAndBinding();
    testNonuniformBarycentricQuadrature();
    testTopologyAndEvaluationValidation();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing transfer check(s) failed\n", failures);
        return 1;
    }
    std::puts("all SimWing transfer checks passed");
    return 0;
}
