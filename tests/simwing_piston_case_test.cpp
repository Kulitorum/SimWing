#include "piston_case.h"
#include "strong_piston_checkpoint_persistence.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using namespace simwing;

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

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "piston: deterministic frame serializes");
    return bytes;
}

std::uint64_t checkpointChecksum(
    const std::span<const std::uint8_t> bytes) {
    std::uint64_t result = 14695981039346656037ULL;
    for (const std::uint8_t value : bytes) {
        result ^= value;
        result *= 1099511628211ULL;
    }
    return result;
}

void writeU64(
    std::vector<std::uint8_t>& bytes,
    const std::size_t offset,
    std::uint64_t value) {
    for (std::size_t index = 0; index < sizeof(value); ++index) {
        bytes[offset + index] = static_cast<std::uint8_t>(value & 0xffU);
        value >>= 8U;
    }
}

void refreshStrongCheckpointChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 28;
    writeU64(
        bytes, 20,
        checkpointChecksum(std::span<const std::uint8_t>(bytes)
                               .subspan(envelopeBytes)));
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const std::string_view name) {
    const auto found = std::ranges::find(
        frame.scalarFields, name, &viewer::ScalarField::name);
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

void testDeterministicEndToEndPistonFrames() {
    fsi::CoupledPistonCase first;
    fsi::CoupledPistonCase observed;
    viewer::DiagnosticFrame finalFrame;
    constexpr std::uint64_t steps = 24;
    for (std::uint64_t step = 0; step < steps; ++step) {
        const auto firstFrame = first.advance();
        finalFrame = observed.advance();
        check(serialized(firstFrame) == serialized(finalFrame),
              "piston: frame observation cannot alter coupled arithmetic");
    }

    const auto firstCheckpoint = first.structure().checkpoint();
    const auto observedCheckpoint = observed.structure().checkpoint();
    check(firstCheckpoint.nodes == observedCheckpoint.nodes
              && firstCheckpoint.acceptedStepCount == steps
              && firstCheckpoint.version == observedCheckpoint.version
              && firstCheckpoint.definitionFingerprint
                  == observedCheckpoint.definitionFingerprint
              && firstCheckpoint.acceptedStepCount
                  == observedCheckpoint.acceptedStepCount
              && firstCheckpoint.simulationTimeSeconds
                  == observedCheckpoint.simulationTimeSeconds
              && firstCheckpoint.pendingExternalForcesNewtons
                  == observedCheckpoint.pendingExternalForcesNewtons
              && firstCheckpoint.lastAppliedExternalForceNewtons
                  == observedCheckpoint.lastAppliedExternalForceNewtons,
          "piston: repeated coupled runs commit bit-identical XPBD state");
    check(finalFrame.sceneChecksum == fsi::coupledPistonCaseChecksum
              && finalFrame.solverCommit == fsi::coupledPistonCaseSolverId
              && finalFrame.step == steps
              && finalFrame.couplingIteration == 1,
          "piston: frame provenance identifies the coupled verification path");
    check(finalFrame.vertices.size() == 4
              && finalFrame.triangles.size() == 2
              && finalFrame.lines.empty(),
          "piston: the visible two-triangle pressure surface is complete");
    check(finalFrame.triangles[0].negativeRegionId == 2
              && finalFrame.triangles[0].positiveRegionId == 1,
          "piston: viewer triangles retain the fluid side-region IDs");

    const double timeStep = observed.stepSettings().timeStepSeconds;
    const double expectedSpeed =
        static_cast<double>(steps) * 660.0 * timeStep / 6000.0;
    const double expectedMomentum =
        static_cast<double>(steps) * 660.0 * timeStep;
    for (const auto& node : observedCheckpoint.nodes) {
        checkNear(node.velocityMetersPerSecond.x, expectedSpeed, 5.0e-13,
                  "piston: every tributary-mass node translates at one speed");
        checkNear(node.velocityMetersPerSecond.y, 0.0, 0.0,
                  "piston: pressure creates no transverse velocity");
        checkNear(node.velocityMetersPerSecond.z, 0.0, 0.0,
                  "piston: pressure creates no vertical velocity");
    }
    checkNear(observed.structure().diagnostics()
                  .linearMomentumKgMetersPerSecond.x,
              expectedMomentum, 3.0e-9,
              "piston: accepted momentum equals accumulated fluid impulse");
    check(observedCheckpoint.nodes[0].positionMeters.x > 3.0,
          "piston: accepted pressure motion is visible in the published geometry");

    const auto* pressure = scalarField(
        finalFrame, "interface.pressure_traction");
    const auto* work = scalarField(finalFrame, "interface.step_work");
    check(pressure != nullptr
              && pressure->association == viewer::FieldAssociation::Triangle
              && pressure->values == std::vector<double>({110.0, 110.0}),
          "piston: replay exposes pressure traction on both triangles");
    check(work != nullptr && work->values.size() == 1
              && work->values.front() > 0.0,
          "piston: replay exposes positive pressure work during acceleration");
    check(finalFrame.couplingResiduals.tractionNewtons < 2.0e-13
              && finalFrame.couplingResiduals.fluid == 0.0
              && finalFrame.couplingResiduals.interfacePowerWatts < 1.0e-12,
          "piston: visible frame carries closed fluid-interface residuals");
    checkNear(finalFrame.conservation.fluidMassKilograms, 28.8, 0.0,
              "piston: visible frame carries the analytic fluid mass");
}

void testCompletedPistonTrace() {
    fsi::CoupledPistonCase simulation;
    std::stringstream trace(std::ios::in | std::ios::out | std::ios::binary);
    viewer::TraceWriter writer(trace);
    check(writer.writeHeader(simulation.traceHeader()),
          "piston trace: header writes");
    constexpr std::uint64_t steps = 5;
    for (std::uint64_t step = 0; step < steps; ++step) {
        check(writer.writeFrame(simulation.advance()),
              "piston trace: accepted frame writes");
    }
    check(writer.finish(), "piston trace: completion record writes");

    trace.seekg(0);
    viewer::TraceReader reader(trace);
    viewer::TraceHeader header;
    check(reader.readHeader(header)
              && header.sceneChecksum == fsi::coupledPistonCaseChecksum
              && header.solverCommit == fsi::coupledPistonCaseSolverId,
          "piston trace: provenance round-trips");
    std::uint64_t frameCount = 0;
    viewer::DiagnosticFrame decoded;
    for (;;) {
        const auto status = reader.readNext(decoded);
        if (status == viewer::TraceReadStatus::Frame) {
            ++frameCount;
            check(decoded.step == frameCount,
                  "piston trace: accepted steps remain consecutive");
        } else {
            check(status == viewer::TraceReadStatus::End,
                  "piston trace: completed replay terminates cleanly");
            break;
        }
    }
    check(frameCount == steps,
          "piston trace: replay contains every accepted frame");
}

void testStrongCoupledLightPiston() {
    fsi::StrongCoupledPistonCase first;
    fsi::StrongCoupledPistonCase observed;
    fsi::StrongCoupledPistonStepDiagnostics firstStep;
    fsi::StrongCoupledPistonStepDiagnostics finalStep;
    constexpr std::uint64_t steps = 8;
    for (std::uint64_t step = 0; step < steps; ++step) {
        const auto firstResult = first.advance();
        const auto observedResult = observed.advance();
        if (step == 0) {
            firstStep = firstResult;
        }
        finalStep = observedResult;
        check(firstResult == observedResult,
              "strong piston: repeated solver runs are deterministic");
        check(observedResult.finite
                  && observedResult.coupling.decision.status
                      == fsi::CouplingMacroStepRetryStatus::Accepted
                  && observedResult.coupling.solverRunCount >= 3
                  && observedResult.coupling.solverRunCount <= 30
                  && observedResult.coupling.attempts.size() == 1
                  && observedResult.coupling.attempts.front().decision.status
                      == fsi::CouplingMacroStepRetryStatus::Accepted
                  && observedResult.coupling.attempts.front().solverRunCount
                      == observedResult.coupling.solverRunCount
                  && observedResult.coupling.lastIteration.status
                      == fsi::StrongCouplingIterationStatus::Converged,
              "strong piston: the real callback converges inside bounded strong coupling");
    }

    const auto firstStructure = first.structure().checkpoint();
    const auto observedStructure = observed.structure().checkpoint();
    check(firstStructure.nodes == observedStructure.nodes
              && firstStructure.acceptedStepCount == steps
              && firstStructure.acceptedStepCount
                  == observedStructure.acceptedStepCount
              && firstStructure.simulationTimeSeconds
                  == observedStructure.simulationTimeSeconds
              && first.fluidState().velocityMetersPerSecond
                  == observed.fluidState().velocityMetersPerSecond
              && first.fluidState().pressurePascals
                  == observed.fluidState().pressurePascals
              && first.fluidState().interfaces
                  == observed.fluidState().interfaces
              && first.fluidState().diagnostics
                  == observed.fluidState().diagnostics,
          "strong piston: accepted Structure and persistent fluid state replay exactly");
    const double explicitWeakSpeed = 660.0
        * observed.stepSettings().timeStepSeconds
        / fsi::strongCoupledPistonStructuralMassKilograms;
    check(firstStep.acceptedSpeedMetersPerSecond > 0.0
              && firstStep.acceptedSpeedMetersPerSecond
                  < explicitWeakSpeed
              && firstStep.coupling.solverRunCount > 1,
          "strong piston: light added mass changes the weak one-pass response and requires iteration");
    checkNear(
        firstStep.measuredDiscreteAddedMassKilograms,
        fsi::strongCoupledPistonDiscreteAddedMassKilograms,
        2.0e-12,
        "strong piston: projection recovers the analytic discrete added mass");
    checkNear(
        firstStep.acceptedSpeedMetersPerSecond,
        firstStep.analyticAcceptedSpeedMetersPerSecond,
        2.0e-14,
        "strong piston: accepted speed matches the trapezoidal added-mass fixed point");
    check(finalStep.acceptedSpeedMetersPerSecond
              > firstStep.acceptedSpeedMetersPerSecond
              && finalStep.velocityClosureMetersPerSecond <= 1.0e-9,
          "strong piston: accepted interface and structure velocity close after repeated steps");
}

void testStrongCoupledPistonFrames() {
    fsi::StrongCoupledPistonWorkerCase first;
    fsi::StrongCoupledPistonWorkerCase observed;
    viewer::DiagnosticFrame finalFrame;
    constexpr std::uint64_t steps = 6;
    for (std::uint64_t step = 0; step < steps; ++step) {
        const auto firstFrame = first.advance();
        finalFrame = observed.advance();
        check(serialized(firstFrame) == serialized(finalFrame),
              "strong piston frame: observation cannot alter coupled arithmetic");
    }
    const auto& diagnostics = observed.diagnostics();
    check(finalFrame.sceneChecksum
              == fsi::strongCoupledPistonCaseChecksum
              && finalFrame.solverCommit
                  == fsi::strongCoupledPistonCaseSolverId
              && finalFrame.step == steps
              && finalFrame.couplingIteration
                  == diagnostics.coupling.lastIteration
                         .convergence.iteration
              && finalFrame.couplingIteration >= 3,
          "strong piston frame: provenance and true coupling iteration are published");
    check(finalFrame.vertices.size() == 4
              && finalFrame.triangles.size() == 2
              && scalarField(finalFrame, "coupling.solver_runs") != nullptr
              && scalarField(finalFrame, "coupling.retry_count") != nullptr
              && scalarField(finalFrame, "coupling.attempt_count") != nullptr
              && scalarField(finalFrame, "interface.speed") != nullptr
              && scalarField(
                  finalFrame, "interface.velocity_closure") != nullptr
              && scalarField(
                  finalFrame, "interface.mean_pressure_traction") != nullptr,
          "strong piston frame: accepted geometry and strong-coupling fields are complete");
    check(scalarField(
              finalFrame, "interface.discrete_added_mass") != nullptr
              && scalarField(
                  finalFrame, "coupling.analytic_speed_residual") != nullptr,
          "strong piston frame: discrete added-mass oracle is published");
    check(finalFrame.couplingResiduals.displacementMetres <= 1.0e-10,
          "strong piston frame: displacement residual meets its acceptance budget");
    check(finalFrame.couplingResiduals.tractionNewtons <= 1.0e-7,
          "strong piston frame: traction residual meets its acceptance budget");
    check(finalFrame.couplingResiduals.fluid <= 1.0e-12,
          "strong piston frame: accepted fluid projection is divergence free");
    check(finalFrame.conservation.fluidMassKilograms == 28.8,
          "strong piston frame: projected fluid mass is explicit");
}

void testStrongPistonAcceptedCheckpointReplay() {
    fsi::StrongCoupledPistonCase simulation;
    for (std::uint64_t step = 0; step < 4; ++step) {
        static_cast<void>(simulation.advance());
    }
    const auto saved = simulation.checkpoint();
    const auto expectedStep = simulation.advance();
    const auto expectedStructure = simulation.structure().checkpoint();
    const auto expectedVelocity =
        simulation.fluidState().velocityMetersPerSecond;
    const auto expectedPressure = simulation.fluidState().pressurePascals;
    const auto expectedInterfaces = simulation.fluidState().interfaces;
    const auto expectedFluidDiagnostics =
        simulation.fluidState().diagnostics;

    simulation.restore(saved);
    const auto replayedStep = simulation.advance();
    const auto replayedStructure = simulation.structure().checkpoint();
    check(replayedStep == expectedStep
              && replayedStructure.nodes == expectedStructure.nodes
              && replayedStructure.acceptedStepCount
                  == expectedStructure.acceptedStepCount
              && replayedStructure.simulationTimeSeconds
                  == expectedStructure.simulationTimeSeconds
              && simulation.fluidState().velocityMetersPerSecond
                  == expectedVelocity
              && simulation.fluidState().pressurePascals
                  == expectedPressure
              && simulation.fluidState().interfaces == expectedInterfaces
              && simulation.fluidState().diagnostics
                  == expectedFluidDiagnostics,
          "strong piston checkpoint: restore replays the exact next coupled step");

    const auto stableStructure = simulation.structure().checkpoint();
    const auto stableVelocity =
        simulation.fluidState().velocityMetersPerSecond;
    const auto stablePressure = simulation.fluidState().pressurePascals;
    const auto stableInterfaces = simulation.fluidState().interfaces;
    auto corrupt = saved;
    ++corrupt.version;
    expectRejected(
        [&] { simulation.restore(corrupt); },
        "strong piston checkpoint: foreign version is rejected");
    corrupt = saved;
    ++corrupt.interfaceDefinitionFingerprint;
    expectRejected(
        [&] { simulation.restore(corrupt); },
        "strong piston checkpoint: foreign interface identity is rejected");
    corrupt = saved;
    ++corrupt.structure.definitionFingerprint;
    expectRejected(
        [&] { simulation.restore(corrupt); },
        "strong piston checkpoint: foreign Structure is rejected");
    corrupt = saved;
    ++corrupt.fluid.topologyFingerprint;
    expectRejected(
        [&] { simulation.restore(corrupt); },
        "strong piston checkpoint: foreign fluid topology is rejected");
    check(simulation.structure().checkpoint().nodes
              == stableStructure.nodes
              && simulation.fluidState().velocityMetersPerSecond
                  == stableVelocity
              && simulation.fluidState().pressurePascals == stablePressure
              && simulation.fluidState().interfaces == stableInterfaces,
          "strong piston checkpoint: rejected restores preserve accepted state");

    fsi::StrongCoupledPistonWorkerCase worker;
    for (std::uint64_t step = 0; step < 3; ++step) {
        static_cast<void>(worker.advance());
    }
    const auto workerSaved = worker.checkpoint();
    const auto expectedFrame = serialized(worker.advance());
    worker.restore(workerSaved);
    const auto replayedFrame = serialized(worker.advance());
    check(replayedFrame == expectedFrame,
          "strong piston checkpoint: immutable next frame replays exactly");
}

void testStrongPistonPersistentCheckpoint() {
    fsi::StrongCoupledPistonCase simulation;
    for (std::uint64_t step = 0; step < 5; ++step) {
        static_cast<void>(simulation.advance());
    }
    const auto checkpoint = simulation.checkpoint();
    std::vector<std::uint8_t> firstBytes;
    std::vector<std::uint8_t> secondBytes;
    fsi::StrongPistonCheckpointPersistenceError error;
    check(fsi::serializeStrongPistonCheckpoint(
              checkpoint, firstBytes, &error)
              && fsi::serializeStrongPistonCheckpoint(
                  checkpoint, secondBytes, &error)
              && !firstBytes.empty()
              && firstBytes == secondBytes,
          "strong piston persistence: accepted state encodes deterministically");

    fsi::StrongCoupledPistonCase unchangedOwner;
    const auto unchanged = unchangedOwner.checkpoint();
    fsi::StrongCoupledPistonCheckpoint decoded = unchanged;
    check(fsi::deserializeStrongPistonCheckpoint(
              firstBytes, decoded, &error),
          "strong piston persistence: encoded state decodes");
    std::vector<std::uint8_t> reencoded;
    check(fsi::serializeStrongPistonCheckpoint(
              decoded, reencoded, &error)
              && reencoded == firstBytes,
          "strong piston persistence: decode and re-encode are byte exact");
    const auto expectedNext = simulation.advance();
    fsi::StrongCoupledPistonCase resumed;
    resumed.restore(decoded);
    const auto resumedNext = resumed.advance();
    check(resumedNext == expectedNext
              && resumed.structure().checkpoint().nodes
                  == simulation.structure().checkpoint().nodes
              && resumed.fluidState().velocityMetersPerSecond
                  == simulation.fluidState().velocityMetersPerSecond
              && resumed.fluidState().pressurePascals
                  == simulation.fluidState().pressurePascals,
          "strong piston persistence: rebuilt owner continues exactly");

    const auto outputUnchanged = [&] {
        return decoded.version == checkpoint.version
            && decoded.interfaceDefinitionFingerprint
                == checkpoint.interfaceDefinitionFingerprint
            && decoded.structure.acceptedStepCount
                == checkpoint.structure.acceptedStepCount
            && decoded.fluid.topologyFingerprint
                == checkpoint.fluid.topologyFingerprint;
    };
    const auto expectDecodeFailure = [&] (
        const std::vector<std::uint8_t>& bytes,
        const fsi::StrongPistonCheckpointPersistenceErrorCode code,
        const char* message,
        const fsi::StrongPistonCheckpointPersistenceLimits& limits = {}) {
        error = {};
        const bool accepted = fsi::deserializeStrongPistonCheckpoint(
            bytes, decoded, &error, limits);
        check(!accepted && error.code == code && outputUnchanged(), message);
    };

    auto corrupt = firstBytes;
    corrupt.front() ^= 0xffU;
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::InvalidMagic,
        "strong piston persistence: invalid magic is transactional");
    corrupt = firstBytes;
    corrupt[8] ^= 0x01U;
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::UnsupportedVersion,
        "strong piston persistence: protocol version is bounded");
    corrupt = firstBytes;
    corrupt[10] = 1U;
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::UnsupportedVersion,
        "strong piston persistence: reserved envelope bits are rejected");
    corrupt = firstBytes;
    corrupt.back() ^= 0x01U;
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::ChecksumMismatch,
        "strong piston persistence: checksum corruption is transactional");
    corrupt = firstBytes;
    corrupt.pop_back();
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::Truncated,
        "strong piston persistence: truncation is rejected");
    corrupt = firstBytes;
    corrupt.push_back(0U);
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::TrailingData,
        "strong piston persistence: trailing data is rejected");

    corrupt = firstBytes;
    corrupt[32] ^= 0x01U;
    refreshStrongCheckpointChecksum(corrupt);
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::InvalidData,
        "strong piston persistence: recomputed-checksum identity edits are rejected");
    corrupt = firstBytes;
    writeU64(corrupt, 40, std::numeric_limits<std::uint64_t>::max());
    refreshStrongCheckpointChecksum(corrupt);
    expectDecodeFailure(
        corrupt,
        fsi::StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
        "strong piston persistence: nested length is bounded before allocation");

    auto smallLimits = fsi::StrongPistonCheckpointPersistenceLimits{};
    smallLimits.maximumEncodedBytes = firstBytes.size() - 1;
    expectDecodeFailure(
        firstBytes,
        fsi::StrongPistonCheckpointPersistenceErrorCode::LimitExceeded,
        "strong piston persistence: outer byte limit is enforced",
        smallLimits);
}

} // namespace

int main() {
    testDeterministicEndToEndPistonFrames();
    testCompletedPistonTrace();
    testStrongCoupledLightPiston();
    testStrongCoupledPistonFrames();
    testStrongPistonAcceptedCheckpointReplay();
    testStrongPistonPersistentCheckpoint();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing coupled piston check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing coupled piston checks passed");
    return 0;
}
