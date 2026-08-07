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

} // namespace

int main() {
    testDeterministicEndToEndPistonFrames();
    testCompletedPistonTrace();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing coupled piston check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing coupled piston checks passed");
    return 0;
}
