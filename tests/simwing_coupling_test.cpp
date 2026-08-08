#include "coupling.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using simwing::fsi::ConservativeMacroStepCoupling;
using simwing::fsi::ConservativeSurfaceTransfer;
using simwing::fsi::ConservativeTransferResult;
using simwing::fsi::ConservativeTransferSettings;
using simwing::fsi::AitkenInterfaceRelaxation;
using simwing::fsi::AitkenRelaxationSettings;
using simwing::fsi::CouplingConvergenceSettings;
using simwing::fsi::CouplingResidualNorms;
using simwing::fsi::CouplingNodeKinematics;
using simwing::fsi::CouplingSurfaceNodeDefinition;
using simwing::fsi::CouplingSurfaceTriangleDefinition;
using simwing::fsi::CouplingTriangleTraction;
using simwing::fsi::Structure;
using simwing::fsi::StructureDefinition;
using simwing::fsi::StructureStepSettings;
using simwing::fsi::StructureVector3;

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
        std::fprintf(
            stderr,
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

StructureDefinition pistonDefinition(const double zOffset = 0.0) {
    StructureDefinition definition;
    // Mass follows the two-triangle barycentric tributary area. A uniform
    // pressure impulse therefore gives every unconstrained node the same
    // velocity increment in the structural acceptance check.
    definition.nodes = {
        {{1.0, -1.0, zOffset}, 2.0, false},
        {{3.0, -1.0, zOffset}, 1.0, false},
        {{3.0, 2.0, zOffset}, 2.0, false},
        {{1.0, 2.0, zOffset}, 1.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

std::vector<CouplingSurfaceNodeDefinition> pistonNodes(
    const std::uint64_t offset = 0) {
    return {
        {40 + offset, 3},
        {10 + offset, 0},
        {30 + offset, 2},
        {20 + offset, 1},
    };
}

std::vector<CouplingSurfaceTriangleDefinition> pistonTriangles(
    const std::uint64_t offset = 0) {
    return {
        {200 + offset, {10 + offset, 30 + offset, 40 + offset}},
        {100 + offset, {10 + offset, 20 + offset, 30 + offset}},
    };
}

struct PistonFixture {
    Structure structure{pistonDefinition()};
    ConservativeSurfaceTransfer transfer{
        structure, pistonNodes(), pistonTriangles()};
    ConservativeMacroStepCoupling coupling{transfer};
};

std::vector<CouplingNodeKinematics> pistonKinematics(
    const ConservativeSurfaceTransfer& transfer,
    const Structure& structure,
    const double timeSeconds) {
    constexpr double speedMetersPerSecond = 0.25;
    std::vector<CouplingNodeKinematics> result;
    for (const auto& node : transfer.nodes()) {
        auto position =
            structure.definition().nodes[node.structureNode].positionMeters;
        position.z += speedMetersPerSecond * timeSeconds;
        result.push_back(
            {node.stableId, position, {0.0, 0.0, speedMetersPerSecond}});
    }
    return result;
}

std::vector<CouplingTriangleTraction> pistonTractions(
    const ConservativeSurfaceTransfer& transfer,
    const double timeSeconds) {
    const double pressurePascals = 100.0 + 50.0 * timeSeconds;
    std::vector<CouplingTriangleTraction> result;
    for (const auto& triangle : transfer.triangles()) {
        result.push_back(
            {triangle.stableId, {0.0, 0.0, pressurePascals}});
    }
    return result;
}

std::vector<ConservativeTransferResult> pistonSamples(
    const PistonFixture& fixture,
    const std::span<const double> offsets,
    const StructureVector3 reference = {2.0, 0.5, 0.0}) {
    ConservativeTransferSettings settings;
    settings.momentReferenceMeters = reference;
    std::vector<ConservativeTransferResult> samples;
    samples.reserve(offsets.size());
    for (const double offset : offsets) {
        samples.push_back(fixture.transfer.evaluate(
            pistonKinematics(fixture.transfer, fixture.structure, offset),
            pistonTractions(fixture.transfer, offset), settings));
    }
    return samples;
}

void testMovingPistonImpulseVolumeAndWork() {
    PistonFixture fixture;
    constexpr std::array<double, 3> offsets{0.0, 0.1, 0.4};
    const auto samples = pistonSamples(fixture, offsets);
    const auto result = fixture.coupling.integrate(offsets, samples);
    const auto replay = fixture.coupling.integrate(offsets, samples);
    const auto& diagnostics = result.diagnostics();

    constexpr double areaSquareMeters = 6.0;
    constexpr double durationSeconds = 0.4;
    constexpr double meanPressurePascals = 110.0;
    const auto startKinematics = pistonKinematics(
        fixture.transfer, fixture.structure, offsets.front());
    const auto endKinematics = pistonKinematics(
        fixture.transfer, fixture.structure, offsets.back());
    const double sweptVolumeCubicMeters = areaSquareMeters
        * (endKinematics.front().positionMeters.z
           - startKinematics.front().positionMeters.z);
    constexpr double analyticImpulseNewtonSeconds =
        areaSquareMeters * meanPressurePascals * durationSeconds;
    const double analyticPressureWorkJoules =
        meanPressurePascals * sweptVolumeCubicMeters;

    check(result == replay,
          "piston: identical temporal samples integrate bit-for-bit");
    check(result.version()
              == simwing::fsi::interfaceImpulseExchangeVersion,
          "piston: immutable result carries the exchange contract version");
    check(result.surfaceFingerprint() == fixture.transfer.fingerprint()
              && result.targetDefinitionFingerprint()
                  == fixture.structure.definitionFingerprint(),
          "piston: immutable result remains bound to surface and Structure");
    check(diagnostics.sampleCount == offsets.size()
              && diagnostics.intervalCount == offsets.size() - 1
              && diagnostics.durationSeconds == durationSeconds,
          "piston: nonuniform sample intervals retain explicit macro-step time");
    checkVectorNear(diagnostics.integratedSurfaceImpulseNewtonSeconds,
                    {0.0, 0.0, analyticImpulseNewtonSeconds}, 2.0e-13,
                    "piston: linear pressure gives the analytic surface impulse");
    checkVectorNear(diagnostics.transferredNodalImpulseNewtonSeconds,
                    {0.0, 0.0, analyticImpulseNewtonSeconds}, 2.0e-13,
                    "piston: barycentric node impulses preserve total impulse");
    checkVectorNear(
        diagnostics.integratedSurfaceAngularImpulseNewtonMeterSeconds,
        {}, 2.0e-13,
        "piston: pressure through the moving centroid has zero angular impulse");
    checkVectorNear(
        diagnostics.transferredNodalAngularImpulseNewtonMeterSeconds,
        {}, 2.0e-13,
        "piston: nodal angular impulse preserves the zero analytic value");
    checkNear(sweptVolumeCubicMeters, 0.6, 2.0e-16,
              "piston: prescribed translation sweeps the analytic volume");
    checkNear(diagnostics.integratedSurfaceWorkJoules,
              analyticPressureWorkJoules, 5.0e-14,
              "piston: integrated surface power equals pressure-volume work");
    checkNear(diagnostics.transferredNodalWorkJoules,
              analyticPressureWorkJoules, 5.0e-14,
              "piston: nodal work preserves pressure-volume work");
    check(diagnostics.impulseResidualNormNewtonSeconds < 3.0e-13
              && diagnostics.angularImpulseResidualNormNewtonMeterSeconds
                  < 3.0e-13
              && std::abs(diagnostics.workResidualJoules) < 1.0e-13
              && diagnostics.finite,
          "piston: temporal impulse, angular impulse, and work ledgers close");

    const auto impulses = result.nodeImpulses();
    check(impulses.size() == 4
              && impulses[0].stableId == 10
              && impulses[3].stableId == 40,
          "piston: node impulses preserve canonical stable-ID order");
    checkVectorNear(impulses[0].impulseNewtonSeconds,
                    {0.0, 0.0, 88.0}, 5.0e-14,
                    "piston: shared diagonal node receives two temporal shares");
    checkVectorNear(impulses[1].impulseNewtonSeconds,
                    {0.0, 0.0, 44.0}, 3.0e-14,
                    "piston: boundary node receives one temporal share");
    checkVectorNear(impulses[2].impulseNewtonSeconds,
                    {0.0, 0.0, 88.0}, 5.0e-14,
                    "piston: second shared node receives two temporal shares");
    checkVectorNear(impulses[3].impulseNewtonSeconds,
                    {0.0, 0.0, 44.0}, 3.0e-14,
                    "piston: final boundary node receives one temporal share");
}

void testAcceptedImpulseReachesStructure() {
    PistonFixture fixture;
    constexpr std::array<double, 3> offsets{0.0, 0.1, 0.4};
    const auto samples = pistonSamples(fixture, offsets);
    const auto transfer = fixture.coupling.integrate(offsets, samples);
    StructureStepSettings settings;
    settings.timeStepSeconds = offsets.back();
    settings.substeps = 4;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.coupling.advanceStructure(
        fixture.structure, transfer, settings);

    check(fixture.structure.acceptedStepCount() == 1
              && fixture.structure.simulationTimeSeconds() == offsets.back(),
          "acceptance: one integrated exchange advances one macro-step");
    checkVectorNear(diagnostics.lastAppliedExternalForceNewtons,
                    {0.0, 0.0, 660.0}, 2.0e-13,
                    "acceptance: average load times duration equals the impulse");
    checkVectorNear(diagnostics.linearMomentumKgMetersPerSecond,
                    {0.0, 0.0, 264.0}, 3.0e-13,
                    "acceptance: actual XPBD momentum gain equals interface impulse");
    checkVectorNear(diagnostics.pendingExternalForceNewtons, {}, 0.0,
                    "acceptance: accepted Structure step consumes the average load");
    const auto states = fixture.structure.nodeStates();
    for (const auto& state : states) {
        checkVectorNear(state.velocityMetersPerSecond,
                        {0.0, 0.0, 44.0}, 8.0e-14,
                        "acceptance: tributary mass gives uniform piston velocity increment");
    }
}

void testValidationAndTransactionalFailure() {
    PistonFixture fixture;
    constexpr std::array<double, 3> offsets{0.0, 0.1, 0.4};
    auto samples = pistonSamples(fixture, offsets);
    const auto integrated = fixture.coupling.integrate(offsets, samples);

    const std::array<double, 1> oneOffset{0.0};
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(
            oneOffset, std::span<const ConservativeTransferResult>(samples).first(1))); },
        "validation: at least two temporal samples are required");
    const std::array<double, 3> nonzeroStart{0.01, 0.1, 0.4};
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(nonzeroStart, samples)); },
        "validation: macro-step-local offsets must start at zero");
    const std::array<double, 3> repeatedOffset{0.0, 0.1, 0.1};
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(repeatedOffset, samples)); },
        "validation: temporal samples must be strictly increasing");
    const std::array<double, 3> nonfiniteOffset{
        0.0, 0.1, std::numeric_limits<double>::quiet_NaN()};
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(nonfiniteOffset, samples)); },
        "validation: temporal sample offsets must be finite");

    ConservativeTransferSettings differentReference;
    differentReference.momentReferenceMeters = {2.0, 0.5, 1.0};
    samples.back() = fixture.transfer.evaluate(
        pistonKinematics(fixture.transfer, fixture.structure, offsets.back()),
        pistonTractions(fixture.transfer, offsets.back()),
        differentReference);
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(offsets, samples)); },
        "validation: one macro-step cannot mix moment references");

    ConservativeSurfaceTransfer otherSurface(
        fixture.structure, pistonNodes(1000), pistonTriangles(1000));
    ConservativeMacroStepCoupling otherCoupling(otherSurface);
    std::vector<ConservativeTransferResult> otherSamples;
    ConservativeTransferSettings reference;
    reference.momentReferenceMeters = {2.0, 0.5, 0.0};
    for (const double offset : offsets) {
        otherSamples.push_back(otherSurface.evaluate(
            pistonKinematics(otherSurface, fixture.structure, offset),
            pistonTractions(otherSurface, offset), reference));
    }
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.integrate(
            offsets, otherSamples)); },
        "validation: instantaneous samples from another surface are rejected");
    const auto otherIntegrated = otherCoupling.integrate(offsets, otherSamples);
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.advanceStructure(
            fixture.structure, otherIntegrated, StructureStepSettings{})); },
        "validation: an exchange from another surface is rejected before mutation");

    StructureStepSettings settings;
    settings.timeStepSeconds = offsets.back();
    settings.substeps = 0;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    fixture.structure.addExternalForce(0, {1.0, 2.0, 3.0});
    const auto before = fixture.structure.checkpoint();
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.advanceStructure(
            fixture.structure, integrated, settings)); },
        "rollback: rejected XPBD settings propagate to the coupling caller");
    const auto after = fixture.structure.checkpoint();
    check(after.acceptedStepCount == before.acceptedStepCount
              && after.simulationTimeSeconds == before.simulationTimeSeconds
              && after.nodes == before.nodes
              && after.pendingExternalForcesNewtons
                  == before.pendingExternalForcesNewtons
              && after.lastAppliedExternalForceNewtons
                  == before.lastAppliedExternalForceNewtons,
          "rollback: failed acceptance restores state from before impulse application");

    settings.substeps = 1;
    settings.timeStepSeconds = 0.2;
    const auto beforeDurationMismatch = fixture.structure.checkpoint();
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.advanceStructure(
            fixture.structure, integrated, settings)); },
        "validation: exchange duration must equal the structural macro-step");
    check(fixture.structure.checkpoint().pendingExternalForcesNewtons
              == beforeDurationMismatch.pendingExternalForcesNewtons,
          "validation: duration mismatch is rejected before pending-load mutation");

    Structure foreign(pistonDefinition(0.25));
    settings.timeStepSeconds = offsets.back();
    const auto foreignBefore = foreign.checkpoint();
    expectRejected(
        [&] { static_cast<void>(fixture.coupling.advanceStructure(
            foreign, integrated, settings)); },
        "validation: foreign Structure definition is rejected");
    check(foreign.checkpoint().nodes == foreignBefore.nodes,
          "validation: foreign Structure rejection is transactional");
}

void testAitkenLinearFixedPointAndBounds() {
    constexpr std::uint64_t definitionFingerprint = 0xa17e'0001ULL;
    AitkenRelaxationSettings settings;
    settings.initialRelaxation = 0.5;
    settings.minimumRelaxation = 0.1;
    settings.maximumRelaxation = 2.0;
    AitkenInterfaceRelaxation relaxation(
        definitionFingerprint, 1, settings);
    std::vector<double> iterate{0.0};
    const std::array<double, 1> firstCandidate{1.0};
    auto diagnostics = relaxation.relax(iterate, firstCandidate, iterate);
    checkNear(iterate.front(), 0.5, 0.0,
              "Aitken: first update uses the configured fixed factor");
    check(diagnostics.completedIterationCount == 1
              && diagnostics.relaxation == 0.5
              && diagnostics.residualL2 == 1.0
              && !diagnostics.usedDynamicRelaxation
              && !diagnostics.relaxationWasClipped
              && diagnostics.finite,
          "Aitken: first update reports its complete accepted state");

    const std::vector<double> secondCandidate{1.0 + 0.5 * iterate.front()};
    diagnostics = relaxation.relax(iterate, secondCandidate, iterate);
    checkNear(iterate.front(), 2.0, 0.0,
              "Aitken: scalar affine fixed point converges on the second update");
    check(diagnostics.completedIterationCount == 2
              && diagnostics.relaxation == 2.0
              && diagnostics.residualL2 == 0.75
              && diagnostics.residualChangeL2 == 0.25
              && diagnostics.denominator == 0.0625
              && diagnostics.usedDynamicRelaxation
              && !diagnostics.relaxationWasClipped,
          "Aitken: vector delta-squared diagnostics match the analytic update");

    AitkenRelaxationSettings boundedSettings;
    boundedSettings.initialRelaxation = 0.5;
    boundedSettings.minimumRelaxation = 0.1;
    boundedSettings.maximumRelaxation = 1.0;
    AitkenInterfaceRelaxation bounded(
        definitionFingerprint, 1, boundedSettings);
    std::vector<double> boundedIterate;
    const std::array<double, 1> zero{0.0};
    const std::array<double, 1> one{1.0};
    const std::array<double, 1> almostOneCandidate{1.49};
    static_cast<void>(bounded.relax(zero, one, boundedIterate));
    diagnostics = bounded.relax(
        boundedIterate, almostOneCandidate, boundedIterate);
    check(diagnostics.relaxation == 1.0
              && diagnostics.usedDynamicRelaxation
              && diagnostics.relaxationWasClipped
              && boundedIterate == std::vector<double>{1.49},
          "Aitken: a near-stationary residual clips the dynamic factor");

    AitkenInterfaceRelaxation degenerate(
        definitionFingerprint, 1, boundedSettings);
    std::vector<double> degenerateIterate;
    static_cast<void>(degenerate.relax(zero, one, degenerateIterate));
    const std::array<double, 1> repeatedResidualCandidate{1.5};
    diagnostics = degenerate.relax(
        degenerateIterate,
        repeatedResidualCandidate,
        degenerateIterate);
    check(diagnostics.relaxation == boundedSettings.initialRelaxation
              && diagnostics.residualChangeL2 == 0.0
              && diagnostics.denominator == 0.0
              && !diagnostics.usedDynamicRelaxation
              && !diagnostics.relaxationWasClipped
              && degenerateIterate == std::vector<double>{1.0},
          "Aitken: an unchanged residual deterministically retains its factor");
}

void testAitkenVectorCheckpointAndTransactionalFailure() {
    constexpr std::uint64_t definitionFingerprint = 0xa17e'0003ULL;
    AitkenRelaxationSettings settings;
    settings.initialRelaxation = 0.4;
    settings.minimumRelaxation = 0.05;
    settings.maximumRelaxation = 1.5;
    AitkenInterfaceRelaxation owner(
        definitionFingerprint, 3, settings);
    std::vector<double> first;
    const std::array<double, 3> initialCurrent{0.0, 1.0, -2.0};
    const std::array<double, 3> initialCandidate{1.0, 0.0, 2.0};
    const auto firstDiagnostics = owner.relax(
        initialCurrent, initialCandidate, first);
    const auto saved = owner.checkpoint();
    const std::array<double, 3> secondCandidate{1.2, 0.2, 1.5};
    std::vector<double> expected;
    const auto expectedDiagnostics = owner.relax(
        first, secondCandidate, expected);

    owner.restore(saved);
    std::vector<double> replay;
    const auto replayDiagnostics = owner.relax(
        first, secondCandidate, replay);
    check(saved.version
              == simwing::fsi::aitkenRelaxationCheckpointVersion
              && saved.interfaceDefinitionFingerprint
                  == definitionFingerprint
              && saved.valueCount == 3
              && saved.settings == settings
              && saved.completedIterationCount == 1
              && saved.relaxation == firstDiagnostics.relaxation
              && saved.previousResidual
                  == std::vector<double>({1.0, -1.0, 4.0})
              && replay == expected
              && replayDiagnostics == expectedDiagnostics,
          "Aitken: checkpoint restore reproduces the next vector update exactly");

    const auto beforeFailure = owner.checkpoint();
    std::vector<double> preserved{7.0, 8.0};
    const std::array<double, 3> nonfiniteCandidate{
        1.2, std::numeric_limits<double>::quiet_NaN(), 1.5};
    expectRejected(
        [&] {
            static_cast<void>(owner.relax(
                first, nonfiniteCandidate, preserved));
        },
        "Aitken: non-finite candidates are rejected");
    check(owner.checkpoint() == beforeFailure
              && preserved == std::vector<double>({7.0, 8.0}),
          "Aitken: rejected updates preserve state and output");

    auto foreign = saved;
    ++foreign.interfaceDefinitionFingerprint;
    const auto beforeRestore = owner.checkpoint();
    expectRejected(
        [&] { owner.restore(foreign); },
        "Aitken: foreign interface checkpoints are rejected");
    foreign = saved;
    ++foreign.valueCount;
    expectRejected(
        [&] { owner.restore(foreign); },
        "Aitken: foreign checkpoint dimensions are rejected");
    foreign = saved;
    foreign.previousResidual.pop_back();
    expectRejected(
        [&] { owner.restore(foreign); },
        "Aitken: truncated checkpoint residuals are rejected");
    check(owner.checkpoint() == beforeRestore,
          "Aitken: rejected checkpoint restores are transactional");

    owner.reset();
    const auto reset = owner.checkpoint();
    check(reset.completedIterationCount == 0
              && reset.relaxation == settings.initialRelaxation
              && reset.previousResidual.empty(),
          "Aitken: reset returns to the exact initial iteration state");

    expectRejected(
        [=] { AitkenInterfaceRelaxation invalid(definitionFingerprint, 0); },
        "Aitken: empty interface vectors are rejected");
    auto invalidSettings = settings;
    invalidSettings.minimumRelaxation = 0.0;
    expectRejected(
        [&] {
            AitkenInterfaceRelaxation invalid(
                definitionFingerprint, 3, invalidSettings);
        },
        "Aitken: non-positive relaxation bounds are rejected");
    expectRejected(
        [] { AitkenInterfaceRelaxation invalid(0, 3); },
        "Aitken: zero interface identities are rejected");
}

void testCouplingConvergenceDecision() {
    CouplingConvergenceSettings settings;
    settings.minimumIterations = 2;
    settings.maximumIterations = 4;
    settings.absoluteDisplacementToleranceMetres = 1.0e-4;
    settings.relativeDisplacementTolerance = 1.0e-2;
    settings.displacementReferenceFloorMetres = 1.0e-2;
    settings.absoluteVelocityToleranceMetersPerSecond = 2.0e-4;
    settings.relativeVelocityTolerance = 2.0e-2;
    settings.velocityReferenceFloorMetersPerSecond = 1.0e-2;
    settings.absoluteTractionToleranceNewtons = 0.5;
    settings.relativeTractionTolerance = 1.0e-2;
    settings.tractionReferenceFloorNewtons = 10.0;
    const CouplingResidualNorms convergedResiduals{
        5.0e-5, 2.0e-2,
        1.0e-4, 1.0e-2,
        0.25, 50.0,
    };

    const auto tooEarly = simwing::fsi::evaluateCouplingConvergence(
        1, convergedResiduals, settings);
    check(tooEarly.displacementConverged
              && tooEarly.velocityConverged
              && tooEarly.tractionConverged
              && !tooEarly.minimumIterationsSatisfied
              && !tooEarly.converged
              && !tooEarly.iterationLimitReached,
          "convergence: all residuals cannot bypass the minimum iteration count");

    const auto accepted = simwing::fsi::evaluateCouplingConvergence(
        2, convergedResiduals, settings);
    checkNear(accepted.relativeDisplacement, 2.5e-3, 1.0e-18,
              "convergence: displacement relative residual uses its reference");
    checkNear(accepted.relativeVelocity, 1.0e-2, 2.0e-18,
              "convergence: velocity relative residual uses its floor");
    checkNear(accepted.relativeTraction, 5.0e-3, 1.0e-18,
              "convergence: traction relative residual uses its reference");
    check(accepted.iteration == 2
              && accepted.residuals == convergedResiduals
              && accepted.displacementConverged
              && accepted.velocityConverged
              && accepted.tractionConverged
              && accepted.minimumIterationsSatisfied
              && accepted.converged
              && !accepted.iterationLimitReached
              && accepted.finite
              && accepted
                  == simwing::fsi::evaluateCouplingConvergence(
                      2, convergedResiduals, settings),
          "convergence: all absolute and relative channels accept deterministically");

    auto oneAbsoluteFailure = convergedResiduals;
    oneAbsoluteFailure.displacementMetres = 1.5e-4;
    oneAbsoluteFailure.displacementReferenceMetres = 1.0;
    const auto absoluteFailure =
        simwing::fsi::evaluateCouplingConvergence(
            2, oneAbsoluteFailure, settings);
    check(absoluteFailure.relativeDisplacement
              < settings.relativeDisplacementTolerance
              && !absoluteFailure.displacementConverged
              && !absoluteFailure.converged,
          "convergence: a relative pass cannot hide an absolute displacement failure");

    auto oneRelativeFailure = convergedResiduals;
    oneRelativeFailure.tractionNewtons = 0.2;
    oneRelativeFailure.tractionReferenceNewtons = 10.0;
    const auto relativeFailure =
        simwing::fsi::evaluateCouplingConvergence(
            settings.maximumIterations,
            oneRelativeFailure,
            settings);
    check(oneRelativeFailure.tractionNewtons
              < settings.absoluteTractionToleranceNewtons
              && relativeFailure.relativeTraction
                  > settings.relativeTractionTolerance
              && !relativeFailure.tractionConverged
              && !relativeFailure.converged
              && relativeFailure.iterationLimitReached,
          "convergence: a lone relative traction failure exhausts the iteration budget");

    auto invalidResiduals = convergedResiduals;
    invalidResiduals.velocityMetersPerSecond =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] {
            static_cast<void>(simwing::fsi::evaluateCouplingConvergence(
                2, invalidResiduals, settings));
        },
        "convergence: non-finite residuals are rejected");
    expectRejected(
        [&] {
            static_cast<void>(simwing::fsi::evaluateCouplingConvergence(
                0, convergedResiduals, settings));
        },
        "convergence: iteration zero is rejected");
    auto invalidSettings = settings;
    invalidSettings.minimumIterations = 5;
    expectRejected(
        [&] {
            static_cast<void>(simwing::fsi::evaluateCouplingConvergence(
                2, convergedResiduals, invalidSettings));
        },
        "convergence: inconsistent iteration limits are rejected");
}

void testTopologyBoundCouplingResidualReduction() {
    PistonFixture fixture;
    auto baseline = pistonKinematics(
        fixture.transfer, fixture.structure, 0.0);
    auto previous = baseline;
    auto current = baseline;
    for (std::size_t index = 0; index < baseline.size(); ++index) {
        previous[index].positionMeters.x += 0.01;
        current[index].positionMeters.x += 0.012;
        previous[index].velocityMetersPerSecond.y += 0.1;
        current[index].velocityMetersPerSecond.y += 0.12;
    }
    ConservativeTransferSettings transferSettings;
    const auto previousTraction = fixture.transfer.evaluate(
        previous, pistonTractions(fixture.transfer, 0.0), transferSettings);
    const auto currentTraction = fixture.transfer.evaluate(
        current, pistonTractions(fixture.transfer, 0.2), transferSettings);
    const auto residuals = fixture.coupling.measureResiduals(
        baseline, previous, current, previousTraction, currentTraction);
    checkNear(residuals.displacementMetres, 0.002, 3.0e-16,
              "residual reduction: maximum iterate displacement is physical");
    checkNear(residuals.displacementReferenceMetres, 0.012, 3.0e-16,
              "residual reduction: displacement reference uses macro-step motion");
    checkNear(residuals.velocityMetersPerSecond, 0.02, 3.0e-16,
              "residual reduction: maximum velocity update is physical");
    checkNear(residuals.velocityReferenceMetersPerSecond, 0.12, 3.0e-16,
              "residual reduction: velocity reference uses macro-step change");
    checkNear(residuals.tractionNewtons, 20.0, 1.0e-13,
              "residual reduction: shared node owns the maximum force update");
    checkNear(residuals.tractionReferenceNewtons, 220.0, 2.0e-13,
              "residual reduction: force reference uses the physical nodal load");

    for (std::size_t index = 0; index < baseline.size(); ++index) {
        baseline[index].positionMeters.z += 8.0;
        previous[index].positionMeters.z += 8.0;
        current[index].positionMeters.z += 8.0;
        baseline[index].velocityMetersPerSecond.z -= 4.0;
        previous[index].velocityMetersPerSecond.z -= 4.0;
        current[index].velocityMetersPerSecond.z -= 4.0;
    }
    const auto shifted = fixture.coupling.measureResiduals(
        baseline, previous, current, previousTraction, currentTraction);
    check(shifted == residuals,
          "residual reduction: rigid origin and bulk velocity shifts cancel exactly");

    auto invalidKinematics = current;
    ++invalidKinematics.front().stableId;
    expectRejected(
        [&] {
            static_cast<void>(fixture.coupling.measureResiduals(
                baseline, previous, invalidKinematics,
                previousTraction, currentTraction));
        },
        "residual reduction: edited stable-ID order is rejected");

    Structure otherStructure(pistonDefinition(0.25));
    ConservativeSurfaceTransfer otherSurface(
        otherStructure, pistonNodes(), pistonTriangles());
    const auto foreignTraction = otherSurface.evaluate(
        pistonKinematics(otherSurface, otherStructure, 0.2),
        pistonTractions(otherSurface, 0.2), transferSettings);
    expectRejected(
        [&] {
            static_cast<void>(fixture.coupling.measureResiduals(
                baseline, previous, current,
                previousTraction, foreignTraction));
        },
        "residual reduction: foreign structural definitions are rejected");
}

} // namespace

int main() {
    testMovingPistonImpulseVolumeAndWork();
    testAcceptedImpulseReachesStructure();
    testValidationAndTransactionalFailure();
    testAitkenLinearFixedPointAndBounds();
    testAitkenVectorCheckpointAndTransactionalFailure();
    testCouplingConvergenceDecision();
    testTopologyBoundCouplingResidualReduction();
    if (failures != 0) {
        std::fprintf(stderr, "%d SimWing coupling check(s) failed\n", failures);
        return 1;
    }
    std::puts("all SimWing coupling checks passed");
    return 0;
}
