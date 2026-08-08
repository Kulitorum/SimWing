#include "coupled_state.h"
#include "fluid/moving_interface.h"

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

StructureDefinition structureDefinition() {
    StructureDefinition definition;
    definition.nodes = {
        {{0.0, 0.0, 0.0}, 1.0, false},
        {{0.0, 1.0, 0.0}, 1.0, false},
        {{0.0, 1.0, 1.0}, 1.0, false},
        {{0.0, 0.0, 1.0}, 1.0, false},
    };
    definition.triangles = {{{0, 1, 2}}, {{0, 2, 3}}};
    return definition;
}

fluid::PeriodicCartesianGrid fluidGrid() {
    return fluid::PeriodicCartesianGrid(
        {4, 4, 4}, {}, {4.0, 4.0, 4.0});
}

fluid::FaceAlignedMovingInterface fluidInterfaces(
    const fluid::PeriodicCartesianGrid& grid) {
    std::vector<fluid::GridFaceMovingInterface> faces;
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            faces.push_back({
                300, 9, 9, fluid::GridFaceAxis::X,
                3, j, k, 0.125,
            });
        }
    }
    return fluid::FaceAlignedMovingInterface(grid, std::move(faces));
}

fluid::MovingInterfaceProjectionSettings fluidSettings() {
    fluid::MovingInterfaceProjectionSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.1;
    settings.projection.absoluteResidualTolerance = 1.0e-11;
    settings.projection.relativeResidualTolerance = 1.0e-13;
    settings.projection.maximumIterations = 1000;
    return settings;
}

fluid::MovingInterfaceFluidState acceptedFluidState(
    const fluid::PeriodicCartesianGrid& grid) {
    fluid::MacVelocityField velocity(grid);
    fluid::CellScalarField pressure(grid);
    auto interfaces = fluidInterfaces(grid);
    const auto diagnostics = fluid::projectVelocityWithMovingInterfaces(
        grid, velocity, pressure, interfaces, fluidSettings());
    return {
        std::move(velocity),
        std::move(pressure),
        std::move(interfaces),
        diagnostics,
    };
}

StrongCouplingIteration strongIteration(
    const std::uint64_t fingerprint) {
    const std::array<double, 2> initial{0.0, 0.0};
    CouplingConvergenceSettings settings;
    settings.minimumIterations = 2;
    settings.maximumIterations = 4;
    return StrongCouplingIteration(fingerprint, initial, {}, settings);
}

StrongCouplingIteration exhaustingIteration(
    const std::uint64_t fingerprint) {
    const std::array<double, 2> initial{0.0, 0.0};
    CouplingConvergenceSettings settings;
    settings.minimumIterations = 1;
    settings.maximumIterations = 1;
    return StrongCouplingIteration(fingerprint, initial, {}, settings);
}

bool sameStructureState(
    const StructureCheckpoint& first,
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

bool sameFluidCheckpoint(
    const fluid::MovingInterfaceFluidCheckpoint& first,
    const fluid::MovingInterfaceFluidCheckpoint& second) {
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> secondBytes;
    return fluid::serializeMovingInterfaceFluidCheckpoint(
               first, firstBytes)
        && fluid::serializeMovingInterfaceFluidCheckpoint(
               second, secondBytes)
        && firstBytes == secondBytes;
}

void mutate(StrongCouplingRollbackState& state) {
    state.structure().addExternalForce(0, {2.0, 0.0, 0.0});
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.1;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    static_cast<void>(state.structure().step(settings));
    state.fluidState().velocityMetersPerSecond.xFaces().front() += 3.0;
    state.fluidState().pressurePascals.values().front() -= 7.0;
    const std::array<double, 2> candidate{1.0, -1.0};
    const CouplingResidualNorms residuals{
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static_cast<void>(state.iteration().advance(candidate, residuals));
}

void mutatePhysicalSolvers(
    Structure& structure,
    fluid::MovingInterfaceFluidState& fluidState,
    const double timeStepSeconds) {
    structure.addExternalForce(0, {2.0, 0.0, 0.0});
    StructureStepSettings settings;
    settings.timeStepSeconds = timeStepSeconds;
    settings.substeps = 1;
    settings.constraintIterations = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    static_cast<void>(structure.step(settings));
    fluidState.velocityMetersPerSecond.xFaces().front()
        += timeStepSeconds;
    fluidState.pressurePascals.values().front()
        -= timeStepSeconds;
}

void testCompositeRollbackAndReplay() {
    constexpr std::uint64_t fingerprint = 0x5c02'0001ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        strongIteration(fingerprint));
    const auto baseline = state.checkpoint();
    const auto baselineStructure = state.structure().checkpoint();
    const auto baselineVelocity =
        state.fluidState().velocityMetersPerSecond;
    const auto baselinePressure = state.fluidState().pressurePascals;
    const auto baselineInterfaces = state.fluidState().interfaces;
    const auto baselineDiagnostics = state.fluidState().diagnostics;
    const auto baselineIteration = state.iteration().checkpoint();

    mutate(state);
    const auto firstMutatedStructure = state.structure().checkpoint();
    const auto firstMutatedVelocity =
        state.fluidState().velocityMetersPerSecond;
    const auto firstMutatedPressure = state.fluidState().pressurePascals;
    const auto firstMutatedIteration = state.iteration().checkpoint();
    state.restore(baseline);
    check(state.interfaceDefinitionFingerprint() == fingerprint
              && sameStructureState(
                  state.structure().checkpoint(), baselineStructure)
              && state.fluidState().velocityMetersPerSecond
                  == baselineVelocity
              && state.fluidState().pressurePascals == baselinePressure
              && state.fluidState().interfaces == baselineInterfaces
              && state.fluidState().diagnostics == baselineDiagnostics
              && state.iteration().checkpoint() == baselineIteration,
          "rollback: Structure, fluid, and iteration restore one exact baseline");

    mutate(state);
    check(sameStructureState(
              state.structure().checkpoint(), firstMutatedStructure)
              && state.fluidState().velocityMetersPerSecond
                  == firstMutatedVelocity
              && state.fluidState().pressurePascals
                  == firstMutatedPressure
              && state.iteration().checkpoint()
                  == firstMutatedIteration,
          "rollback: restored owners replay the same next mutation exactly");

    state.restore(baseline);
    state.fluidState().diagnostics.finite = false;
    expectRejected(
        [&] { static_cast<void>(state.checkpoint()); },
        "rollback: an unaccepted current fluid state is not checkpointable");
    state.restore(baseline);
}

void testCompositeRestoreIsTransactional() {
    constexpr std::uint64_t fingerprint = 0x5c02'0002ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        strongIteration(fingerprint));
    const auto baseline = state.checkpoint();
    const auto expectedStructure = state.structure().checkpoint();
    const auto expectedVelocity =
        state.fluidState().velocityMetersPerSecond;
    const auto expectedPressure = state.fluidState().pressurePascals;
    const auto expectedIteration = state.iteration().checkpoint();
    const auto unchanged = [&] {
        return sameStructureState(
                   state.structure().checkpoint(), expectedStructure)
            && state.fluidState().velocityMetersPerSecond
                == expectedVelocity
            && state.fluidState().pressurePascals == expectedPressure
            && state.iteration().checkpoint() == expectedIteration;
    };

    auto corrupt = baseline;
    ++corrupt.interfaceDefinitionFingerprint;
    expectRejected(
        [&] { state.restore(corrupt); },
        "rollback: a foreign composite identity is rejected");
    check(unchanged(),
          "rollback: identity rejection preserves every owner");

    corrupt = baseline;
    ++corrupt.structure.definitionFingerprint;
    expectRejected(
        [&] { state.restore(corrupt); },
        "rollback: a foreign Structure checkpoint is rejected");
    check(unchanged(),
          "rollback: Structure rejection preserves every owner");

    corrupt = baseline;
    ++corrupt.fluid.topologyFingerprint;
    expectRejected(
        [&] { state.restore(corrupt); },
        "rollback: edited fluid topology is rejected");
    check(unchanged(),
          "rollback: fluid rejection preserves every owner");

    corrupt = baseline;
    corrupt.iteration.currentInterface.front() =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { state.restore(corrupt); },
        "rollback: non-finite iteration state is rejected");
    check(unchanged(),
          "rollback: iteration rejection preserves every owner");
}

void testSolverOnlyRollbackPreservesIteration() {
    constexpr std::uint64_t fingerprint = 0x5c02'0003ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        strongIteration(fingerprint));
    const auto baseline = state.solverCheckpoint();
    mutate(state);
    const auto advancedIteration = state.iteration().checkpoint();
    const auto mutatedSolvers = state.solverCheckpoint();

    state.restoreSolvers(baseline);
    const auto restored = state.solverCheckpoint();
    check(sameStructureState(restored.structure, baseline.structure)
              && sameFluidCheckpoint(restored.fluid, baseline.fluid)
              && state.iteration().checkpoint() == advancedIteration,
          "solver rollback: physical owners rewind while Aitken state remains");

    auto corrupt = baseline;
    ++corrupt.interfaceDefinitionFingerprint;
    expectRejected(
        [&] { state.restoreSolvers(corrupt); },
        "solver rollback: foreign identity is rejected");
    corrupt = baseline;
    ++corrupt.structure.definitionFingerprint;
    expectRejected(
        [&] { state.restoreSolvers(corrupt); },
        "solver rollback: foreign Structure state is rejected");
    corrupt = baseline;
    ++corrupt.fluid.topologyFingerprint;
    expectRejected(
        [&] { state.restoreSolvers(corrupt); },
        "solver rollback: foreign fluid topology is rejected");
    const auto afterRejections = state.solverCheckpoint();
    check(sameStructureState(
              afterRejections.structure, baseline.structure)
              && sameFluidCheckpoint(
                  afterRejections.fluid, baseline.fluid)
              && state.iteration().checkpoint() == advancedIteration,
          "solver rollback: rejected replacements preserve every live owner");

    state.restoreSolvers(mutatedSolvers);
    check(state.iteration().checkpoint() == advancedIteration,
          "solver rollback: alternate valid physical epochs never alter iteration");
}

void testMacroStepRetryPolicyAndReplay() {
    constexpr std::uint64_t fingerprint = 0x5c03'0001ULL;
    CouplingMacroStepRetrySettings settings;
    settings.maximumRetries = 4;
    settings.reductionFactor = 0.5;
    settings.minimumTimeStepSeconds = 0.01;
    CouplingMacroStepRetry retry(fingerprint, 0.08, settings);

    const auto pending = retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted);
    check(pending.status == CouplingMacroStepRetryStatus::RetryPending
              && pending.attemptNumber == 1
              && pending.retryCount == 0
              && pending.timeStepSeconds == 0.08
              && pending.pendingTimeStepSeconds == 0.04,
          "retry: exhaustion proposes one bounded smaller macro-step");
    const auto savedPending = retry.checkpoint();
    const auto secondAttempt = retry.beginRetry();
    check(secondAttempt.status == CouplingMacroStepRetryStatus::Attempting
              && secondAttempt.attemptNumber == 2
              && secondAttempt.retryCount == 1
              && secondAttempt.timeStepSeconds == 0.04
              && secondAttempt.pendingTimeStepSeconds == 0.0,
          "retry: explicit acknowledgement activates the pending step");

    static_cast<void>(retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted));
    static_cast<void>(retry.beginRetry());
    static_cast<void>(retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted));
    static_cast<void>(retry.beginRetry());
    check(retry.timeStepSeconds() == settings.minimumTimeStepSeconds
              && retry.retryCount() == 3,
          "retry: reduction clamps exactly at the minimum step");
    const auto failed = retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted);
    check(failed.status == CouplingMacroStepRetryStatus::Failed
              && failed.timeStepSeconds
                  == settings.minimumTimeStepSeconds
              && failed.pendingTimeStepSeconds == 0.0,
          "retry: exhaustion at the minimum is terminal failure");

    retry.restore(savedPending);
    check(retry.checkpoint() == savedPending,
          "retry: pending restore reproduces the exact handshake state");
    const auto replaySecond = retry.beginRetry();
    check(replaySecond == secondAttempt,
          "retry: restored pending state replays the next attempt exactly");
    const auto accepted = retry.reportIteration(
        StrongCouplingIterationStatus::Converged);
    check(accepted.status == CouplingMacroStepRetryStatus::Accepted
              && accepted.timeStepSeconds == 0.04,
          "retry: converged iteration accepts its active reduced step");

    const auto terminal = retry.checkpoint();
    bool terminalRejected = false;
    try {
        static_cast<void>(retry.reportIteration(
            StrongCouplingIterationStatus::Converged));
    } catch (const std::logic_error&) {
        terminalRejected = true;
    }
    check(terminalRejected && retry.checkpoint() == terminal,
          "retry: accepted state is terminal and immutable");
}

void testMacroStepRetryValidationIsTransactional() {
    constexpr std::uint64_t fingerprint = 0x5c03'0002ULL;
    CouplingMacroStepRetrySettings settings;
    settings.maximumRetries = 1;
    settings.reductionFactor = 0.5;
    settings.minimumTimeStepSeconds = 0.01;
    CouplingMacroStepRetry retry(fingerprint, 0.08, settings);
    const auto initial = retry.checkpoint();

    bool activeRejected = false;
    try {
        static_cast<void>(retry.beginRetry());
    } catch (const std::logic_error&) {
        activeRejected = true;
    }
    expectRejected(
        [&] {
            static_cast<void>(retry.reportIteration(
                StrongCouplingIterationStatus::Iterating));
        },
        "retry: a nonterminal strong iteration cannot finish an attempt");
    check(activeRejected && retry.checkpoint() == initial,
          "retry: invalid sequencing preserves policy state");

    static_cast<void>(retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted));
    static_cast<void>(retry.beginRetry());
    const auto failed = retry.reportIteration(
        StrongCouplingIterationStatus::Exhausted);
    check(failed.status == CouplingMacroStepRetryStatus::Failed
              && failed.retryCount == settings.maximumRetries,
          "retry: the configured retry count is a hard terminal budget");
    const auto terminal = retry.checkpoint();

    auto corrupt = terminal;
    ++corrupt.macroStepDefinitionFingerprint;
    expectRejected(
        [&] { retry.restore(corrupt); },
        "retry: foreign policy identity is rejected");
    corrupt = terminal;
    corrupt.timeStepSeconds *= 0.75;
    expectRejected(
        [&] { retry.restore(corrupt); },
        "retry: a step inconsistent with retry history is rejected");
    corrupt = terminal;
    corrupt.status = CouplingMacroStepRetryStatus::RetryPending;
    corrupt.pendingTimeStepSeconds = settings.minimumTimeStepSeconds;
    expectRejected(
        [&] { retry.restore(corrupt); },
        "retry: a pending step beyond the retry budget is rejected");
    check(retry.checkpoint() == terminal,
          "retry: rejected restores preserve terminal policy state");

    auto invalidSettings = settings;
    invalidSettings.maximumRetries =
        maximumCouplingMacroStepRetries + 1;
    expectRejected(
        [&] {
            CouplingMacroStepRetry invalid(
                fingerprint, 0.08, invalidSettings);
        },
        "retry: excessive retry bounds are rejected");
}

void testMacroStepStateRestoresBeforeRetry() {
    constexpr std::uint64_t fingerprint = 0x5c04'0001ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        exhaustingIteration(fingerprint));
    CouplingMacroStepRetrySettings retrySettings;
    retrySettings.maximumRetries = 2;
    retrySettings.reductionFactor = 0.5;
    retrySettings.minimumTimeStepSeconds = 0.01;
    StrongCouplingMacroStepState macroStep(
        std::move(state), 0.08, retrySettings);
    const auto baseline = macroStep.rollbackState().checkpoint();

    const auto beforeEarlyReport = macroStep.decision();
    expectRejected(
        [&] {
            static_cast<void>(macroStep.reportTerminalIteration());
        },
        "macro-step: an active iteration cannot finish an attempt");
    check(macroStep.decision() == beforeEarlyReport,
          "macro-step: early reporting preserves retry state");

    mutate(macroStep.rollbackState());
    check(macroStep.rollbackState().iteration().status()
              == StrongCouplingIterationStatus::Exhausted,
          "macro-step: the synthetic first attempt exhausts its budget");
    const auto firstMutation = macroStep.rollbackState().checkpoint();
    const auto pending = macroStep.reportTerminalIteration();
    check(pending.status == CouplingMacroStepRetryStatus::RetryPending
              && pending.timeStepSeconds == 0.08
              && pending.pendingTimeStepSeconds == 0.04,
          "macro-step: exhaustion enters the restore-required pending state");

    const auto restarted = macroStep.restoreAndBeginRetry();
    const auto restored = macroStep.rollbackState().checkpoint();
    check(restarted.status == CouplingMacroStepRetryStatus::Attempting
              && restarted.attemptNumber == 2
              && restarted.timeStepSeconds == 0.04
              && sameStructureState(restored.structure, baseline.structure)
              && sameFluidCheckpoint(restored.fluid, baseline.fluid)
              && restored.iteration == baseline.iteration,
          "macro-step: retry activation first restores all three owners");

    mutate(macroStep.rollbackState());
    const auto replay = macroStep.rollbackState().checkpoint();
    check(sameStructureState(replay.structure, firstMutation.structure)
              && sameFluidCheckpoint(replay.fluid, firstMutation.fluid)
              && replay.iteration == firstMutation.iteration,
          "macro-step: restored retry replays the same attempted mutation");

    static_cast<void>(macroStep.reportTerminalIteration());
    const auto pendingState = macroStep.rollbackState().checkpoint();
    const auto pendingDecision = macroStep.decision();
    bool duplicateReportRejected = false;
    try {
        static_cast<void>(macroStep.reportTerminalIteration());
    } catch (const std::logic_error&) {
        duplicateReportRejected = true;
    }
    check(duplicateReportRejected
              && macroStep.decision() == pendingDecision
              && sameStructureState(
                  macroStep.rollbackState().checkpoint().structure,
                  pendingState.structure),
          "macro-step: duplicate terminal reports preserve pending state");
}

void testMacroStepStateAcceptanceAndFreshBaseline() {
    constexpr std::uint64_t fingerprint = 0x5c04'0002ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        strongIteration(fingerprint));
    StrongCouplingMacroStepState macroStep(std::move(state), 0.08);
    const std::array<double, 2> candidate{0.0, 0.0};
    static_cast<void>(macroStep.rollbackState().iteration().advance(
        candidate, {}));
    const auto converged = macroStep.rollbackState().iteration().advance(
        candidate, {});
    check(converged.status == StrongCouplingIterationStatus::Converged,
          "macro-step: synthetic zero residual reaches terminal convergence");
    const auto accepted = macroStep.reportTerminalIteration();
    check(accepted.status == CouplingMacroStepRetryStatus::Accepted
              && accepted.attemptNumber == 1
              && accepted.timeStepSeconds == 0.08,
          "macro-step: converged state accepts the active step directly");
    const auto acceptedState = macroStep.rollbackState().checkpoint();
    bool retryRejected = false;
    try {
        static_cast<void>(macroStep.restoreAndBeginRetry());
    } catch (const std::logic_error&) {
        retryRejected = true;
    }
    check(retryRejected
              && macroStep.decision() == accepted
              && macroStep.rollbackState().iteration().checkpoint()
                  == acceptedState.iteration,
          "macro-step: accepted state cannot roll back into a retry");

    expectRejected(
        [&] {
            StrongCouplingRollbackState advancedState(
                fingerprint + 1,
                Structure(structureDefinition()),
                grid,
                acceptedFluidState(grid),
                strongIteration(fingerprint + 1));
            static_cast<void>(advancedState.iteration().advance(
                candidate, {}));
            StrongCouplingMacroStepState invalid(
                std::move(advancedState), 0.08);
        },
        "macro-step: a partially advanced iteration is not a baseline");
}

void testMacroStepStateRewindsSolversBetweenIterations() {
    constexpr std::uint64_t fingerprint = 0x5c04'0003ULL;
    const auto grid = fluidGrid();
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        strongIteration(fingerprint));
    StrongCouplingMacroStepState macroStep(std::move(state), 0.08);
    const auto baseline = macroStep.rollbackState().solverCheckpoint();

    bool earlyRejected = false;
    try {
        macroStep.restoreSolversForNextIteration();
    } catch (const std::logic_error&) {
        earlyRejected = true;
    }
    check(earlyRejected,
          "macro-step iteration: an unadvanced baseline cannot be rewound");

    mutate(macroStep.rollbackState());
    const auto advancedIteration =
        macroStep.rollbackState().iteration().checkpoint();
    macroStep.restoreSolversForNextIteration();
    const auto restored = macroStep.rollbackState().solverCheckpoint();
    check(sameStructureState(restored.structure, baseline.structure)
              && sameFluidCheckpoint(restored.fluid, baseline.fluid)
              && macroStep.rollbackState().iteration().checkpoint()
                  == advancedIteration
              && macroStep.decision().status
                  == CouplingMacroStepRetryStatus::Attempting,
          "macro-step iteration: next solve starts at the physical baseline with advanced Aitken state");
}

StrongCouplingMacroStepState runnerMacroStep(
    const std::uint64_t fingerprint,
    const CouplingConvergenceSettings& convergenceSettings,
    const CouplingMacroStepRetrySettings& retrySettings) {
    const auto grid = fluidGrid();
    const std::array<double, 2> initial{0.0, 0.0};
    StrongCouplingRollbackState state(
        fingerprint,
        Structure(structureDefinition()),
        grid,
        acceptedFluidState(grid),
        StrongCouplingIteration(
            fingerprint, initial, {}, convergenceSettings));
    return StrongCouplingMacroStepState(
        std::move(state), 0.08, retrySettings);
}

void testStrongCouplingRunnerConvergesAfterRetry() {
    constexpr std::uint64_t fingerprint = 0x5c05'0001ULL;
    CouplingConvergenceSettings convergence;
    convergence.minimumIterations = 2;
    convergence.maximumIterations = 2;
    CouplingMacroStepRetrySettings retries;
    retries.maximumRetries = 2;
    retries.reductionFactor = 0.5;
    retries.minimumTimeStepSeconds = 0.01;
    auto macroStep = runnerMacroStep(fingerprint, convergence, retries);
    const auto baseline = macroStep.rollbackState().solverCheckpoint();
    const double baselineFluidVelocity = macroStep.rollbackState()
        .fluidState().velocityMetersPerSecond.xFaces().front();
    std::uint64_t largeStepRuns = 0;
    std::uint64_t reducedStepRuns = 0;
    bool everyRunStartedAtBaseline = true;

    const StrongCouplingSolverCallback solve =
        [&](Structure& structure,
            const fluid::PeriodicCartesianGrid&,
            fluid::MovingInterfaceFluidState& fluidState,
            const std::span<const double> current,
            const double timeStepSeconds) {
            const auto currentStructure = structure.checkpoint();
            everyRunStartedAtBaseline = everyRunStartedAtBaseline
                && sameStructureState(
                    currentStructure, baseline.structure)
                && fluidState.velocityMetersPerSecond.xFaces().front()
                    == baselineFluidVelocity;
            mutatePhysicalSolvers(
                structure, fluidState, timeStepSeconds);

            const bool largeStep = timeStepSeconds == 0.08;
            std::uint64_t& runs = largeStep
                ? largeStepRuns : reducedStepRuns;
            ++runs;
            StrongCouplingSolverResult result;
            result.unrelaxedInterface.assign(
                current.begin(), current.end());
            if (largeStep || runs == 1) {
                result.residuals = {
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            }
            return result;
        };

    const auto result = runStrongCouplingMacroStep(macroStep, solve);
    const auto& state = macroStep.rollbackState();
    check(result.decision.status
              == CouplingMacroStepRetryStatus::Accepted
              && result.decision.retryCount == 1
              && result.decision.timeStepSeconds == 0.04
              && result.lastIteration.status
                  == StrongCouplingIterationStatus::Converged
              && result.solverRunCount == 4
              && result.attempts.size() == 2
              && result.attempts[0].decision.status
                  == CouplingMacroStepRetryStatus::RetryPending
              && result.attempts[0].decision.timeStepSeconds == 0.08
              && result.attempts[0].terminalIteration.status
                  == StrongCouplingIterationStatus::Exhausted
              && result.attempts[0].solverRunCount == 2
              && result.attempts[1].decision.status
                  == CouplingMacroStepRetryStatus::Accepted
              && result.attempts[1].decision.timeStepSeconds == 0.04
              && result.attempts[1].terminalIteration.status
                  == StrongCouplingIterationStatus::Converged
              && result.attempts[1].solverRunCount == 2
              && largeStepRuns == 2
              && reducedStepRuns == 2,
          "runner: exhausted large step retries and converges at the reduced step");
    check(everyRunStartedAtBaseline
              && state.structure().acceptedStepCount() == 1
              && state.structure().simulationTimeSeconds() == 0.04
              && state.fluidState().velocityMetersPerSecond
                     .xFaces().front() == baselineFluidVelocity + 0.04
              && state.iteration().completedIterationCount() == 2,
          "runner: each solve rewinds physically while accepted final state is retained");
}

void testStrongCouplingRunnerFailureAndExceptionRollback() {
    constexpr std::uint64_t fingerprint = 0x5c05'0002ULL;
    CouplingConvergenceSettings convergence;
    convergence.minimumIterations = 1;
    convergence.maximumIterations = 1;
    CouplingMacroStepRetrySettings retries;
    retries.maximumRetries = 1;
    retries.reductionFactor = 0.5;
    retries.minimumTimeStepSeconds = 0.01;
    auto failedMacroStep = runnerMacroStep(
        fingerprint, convergence, retries);
    const auto failedBaseline =
        failedMacroStep.rollbackState().checkpoint();
    const StrongCouplingSolverCallback neverConverges =
        [](Structure& structure,
           const fluid::PeriodicCartesianGrid&,
           fluid::MovingInterfaceFluidState& fluidState,
           const std::span<const double> current,
           const double timeStepSeconds) {
            mutatePhysicalSolvers(
                structure, fluidState, timeStepSeconds);
            StrongCouplingSolverResult result;
            result.unrelaxedInterface.assign(
                current.begin(), current.end());
            result.residuals = {
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            return result;
        };
    const auto failed = runStrongCouplingMacroStep(
        failedMacroStep, neverConverges);
    const auto failedState = failedMacroStep.rollbackState().checkpoint();
    check(failed.decision.status
              == CouplingMacroStepRetryStatus::Failed
              && failed.solverRunCount == 2
              && failed.attempts.size() == 2
              && failed.attempts[0].decision.status
                  == CouplingMacroStepRetryStatus::RetryPending
              && failed.attempts[1].decision.status
                  == CouplingMacroStepRetryStatus::Failed
              && failed.attempts[0].solverRunCount == 1
              && failed.attempts[1].solverRunCount == 1
              && sameStructureState(
                  failedState.structure, failedBaseline.structure)
              && sameFluidCheckpoint(
                  failedState.fluid, failedBaseline.fluid)
              && failedState.iteration == failedBaseline.iteration,
          "runner: terminal retry failure leaves the accepted baseline intact");

    auto throwingMacroStep = runnerMacroStep(
        fingerprint + 1, convergence, retries);
    const auto throwingBaseline =
        throwingMacroStep.rollbackState().checkpoint();
    const StrongCouplingSolverCallback returnsUnacceptedFluid =
        [](Structure& structure,
           const fluid::PeriodicCartesianGrid&,
           fluid::MovingInterfaceFluidState& fluidState,
           const std::span<const double> current,
           const double timeStepSeconds) {
            mutatePhysicalSolvers(
                structure, fluidState, timeStepSeconds);
            fluidState.diagnostics.finite = false;
            StrongCouplingSolverResult result;
            result.unrelaxedInterface.assign(
                current.begin(), current.end());
            return result;
        };
    expectRejected(
        [&] {
            static_cast<void>(runStrongCouplingMacroStep(
                throwingMacroStep, returnsUnacceptedFluid));
        },
        "runner: an unaccepted callback fluid epoch is rejected");
    const auto afterInvalidFluid =
        throwingMacroStep.rollbackState().checkpoint();
    check(sameStructureState(
              afterInvalidFluid.structure, throwingBaseline.structure)
              && sameFluidCheckpoint(
                  afterInvalidFluid.fluid, throwingBaseline.fluid)
              && afterInvalidFluid.iteration
                  == throwingBaseline.iteration,
          "runner: rejected physical results restore the complete baseline");

    const StrongCouplingSolverCallback throwsAfterMutation =
        [](Structure& structure,
           const fluid::PeriodicCartesianGrid&,
           fluid::MovingInterfaceFluidState& fluidState,
           std::span<const double>,
           const double timeStepSeconds) -> StrongCouplingSolverResult {
            mutatePhysicalSolvers(
                structure, fluidState, timeStepSeconds);
            throw std::runtime_error("synthetic coupled solver failure");
        };
    bool exceptionPropagated = false;
    try {
        static_cast<void>(runStrongCouplingMacroStep(
            throwingMacroStep, throwsAfterMutation));
    } catch (const std::runtime_error&) {
        exceptionPropagated = true;
    }
    const auto throwingState =
        throwingMacroStep.rollbackState().checkpoint();
    check(exceptionPropagated
              && throwingMacroStep.decision().status
                  == CouplingMacroStepRetryStatus::Attempting
              && sameStructureState(
                  throwingState.structure, throwingBaseline.structure)
              && sameFluidCheckpoint(
                  throwingState.fluid, throwingBaseline.fluid)
              && throwingState.iteration == throwingBaseline.iteration,
          "runner: callback exceptions propagate only after full baseline rollback");

    const auto beforeEmpty = throwingMacroStep.rollbackState().checkpoint();
    expectRejected(
        [&] {
            static_cast<void>(runStrongCouplingMacroStep(
                throwingMacroStep, {}));
        },
        "runner: an empty solver callback is rejected");
    check(throwingMacroStep.rollbackState().iteration().checkpoint()
              == beforeEmpty.iteration,
          "runner: invalid startup preserves the fresh macro-step");
}

} // namespace

int main() {
    testCompositeRollbackAndReplay();
    testCompositeRestoreIsTransactional();
    testSolverOnlyRollbackPreservesIteration();
    testMacroStepRetryPolicyAndReplay();
    testMacroStepRetryValidationIsTransactional();
    testMacroStepStateRestoresBeforeRetry();
    testMacroStepStateAcceptanceAndFreshBaseline();
    testMacroStepStateRewindsSolversBetweenIterations();
    testStrongCouplingRunnerConvergesAfterRetry();
    testStrongCouplingRunnerFailureAndExceptionRollback();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d strong-coupling rollback check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all strong-coupling rollback checks passed");
    return 0;
}
