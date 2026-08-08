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

} // namespace

int main() {
    testCompositeRollbackAndReplay();
    testCompositeRestoreIsTransactional();
    testMacroStepRetryPolicyAndReplay();
    testMacroStepRetryValidationIsTransactional();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d strong-coupling rollback check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all strong-coupling rollback checks passed");
    return 0;
}
