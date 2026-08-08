#include "piston_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
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

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "piston: deterministic frame serializes");
    return bytes;
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

} // namespace

int main() {
    testDeterministicEndToEndPistonFrames();
    testCompletedPistonTrace();
    testStrongCoupledLightPiston();
    testStrongCoupledPistonFrames();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing coupled piston check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing coupled piston checks passed");
    return 0;
}
