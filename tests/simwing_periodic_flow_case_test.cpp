#include "periodic_flow_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>
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
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "periodic fluid frame serializes");
    return bytes;
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

void testDeterministicAcceptedFrames() {
    fsi::PeriodicFlowCase first;
    fsi::PeriodicFlowCase second;
    viewer::DiagnosticFrame firstFrame;
    viewer::DiagnosticFrame secondFrame;
    constexpr std::uint64_t steps = 12;
    for (std::uint64_t step = 0; step < steps; ++step) {
        firstFrame = first.advance();
        secondFrame = second.advance();
        check(serialized(firstFrame) == serialized(secondFrame),
              "periodic flow case and frame replay bit-for-bit");
        check(first.velocity() == second.velocity()
                  && first.pressure() == second.pressure()
                  && first.diagnostics() == second.diagnostics(),
              "frame observation cannot change committed fluid state");
    }

    check(firstFrame.step == steps
              && firstFrame.sceneChecksum == fsi::periodicFlowCaseChecksum
              && firstFrame.solverCommit == fsi::periodicFlowCaseSolverId,
          "periodic fluid frame provenance and step are stable");
    checkNear(firstFrame.simulationTimeSeconds,
              steps * first.stepSettings().flow.timeStepSeconds,
              1.0e-15,
              "periodic fluid frame time follows accepted outer intervals");
    check(firstFrame.vertices.size() == first.grid().cellCount()
              && firstFrame.triangles.empty()
              && firstFrame.lines.empty(),
          "periodic fluid frame publishes one unconnected point per cell");
    check(firstFrame.scalarFields.size() == 7
              && firstFrame.scalarFields[0].name == "pressure"
              && firstFrame.scalarFields[0].association
                  == viewer::FieldAssociation::Vertex
              && firstFrame.scalarFields[0].values.size()
                  == first.grid().cellCount()
              && firstFrame.vectorFields.size() == 2
              && firstFrame.vectorFields[0].name == "velocity"
              && firstFrame.vectorFields[0].values.size()
                  == first.grid().cellCount(),
          "periodic fluid frame retains pressure, speed, velocity, and ledgers");
    check(std::ranges::any_of(
              firstFrame.scalarFields[0].values,
              [](const double pressure) { return pressure != 0.0; })
              && first.diagnostics().accepted
              && firstFrame.couplingResiduals.fluid
                  == first.diagnostics().finalDivergenceL2PerSecond,
          "periodic fluid frame comes from an accepted nontrivial pressure solve");
    checkNear(firstFrame.conservation.fluidMassKilograms,
              first.stepSettings().flow.densityKgPerCubicMeter
                  * first.grid().cellVolumeCubicMeters()
                  * static_cast<double>(first.grid().cellCount()),
              0.0,
              "periodic fluid frame mass is the exact periodic domain mass");

    const auto ownedFrame = firstFrame;
    const auto ownedBytes = serialized(ownedFrame);
    static_cast<void>(first.advance());
    check(serialized(ownedFrame) == ownedBytes,
          "periodic fluid frame owns data independent of later solver steps");
}

void testFrameRejectsUnacceptedState() {
    fsi::PeriodicFlowCase simulation;
    viewer::PeriodicFluidFrameContext context;
    context.sceneChecksum = fsi::periodicFlowCaseChecksum;
    context.solverCommit = fsi::periodicFlowCaseSolverId;
    context.step = 1;
    context.simulationTimeSeconds =
        simulation.stepSettings().flow.timeStepSeconds;
    context.densityKgPerCubicMeter =
        simulation.stepSettings().flow.densityKgPerCubicMeter;
    const fsi::fluid::PeriodicFlowStrangSubcyclingDiagnostics unaccepted;
    expectRejected(
        [&] { static_cast<void>(viewer::buildPeriodicFluidFrame(
            simulation.grid(), simulation.velocity(), simulation.pressure(),
            unaccepted, context)); },
        "periodic fluid frame rejects an unaccepted solver state");
}

void testCompletedTrace() {
    fsi::PeriodicFlowCase simulation;
    std::stringstream trace(std::ios::in | std::ios::out | std::ios::binary);
    viewer::TraceWriter writer(trace);
    check(writer.writeHeader(simulation.traceHeader()),
          "periodic fluid trace header writes");
    constexpr std::uint64_t steps = 5;
    for (std::uint64_t step = 0; step < steps; ++step) {
        check(writer.writeFrame(simulation.advance()),
              "accepted periodic fluid frame writes");
    }
    check(writer.finish(), "periodic fluid trace receives an end marker");
    trace.seekg(0);
    viewer::TraceReader reader(trace);
    viewer::TraceHeader header;
    check(reader.readHeader(header)
              && header.sceneChecksum == fsi::periodicFlowCaseChecksum
              && header.solverCommit == fsi::periodicFlowCaseSolverId,
          "periodic fluid trace header replays");
    std::uint64_t frames = 0;
    viewer::DiagnosticFrame decoded;
    for (;;) {
        const viewer::TraceReadStatus status = reader.readNext(decoded);
        if (status == viewer::TraceReadStatus::Frame) {
            ++frames;
            check(decoded.step == frames
                      && decoded.vertices.size()
                          == simulation.grid().cellCount(),
                  "periodic fluid trace retains consecutive point snapshots");
        } else {
            check(status == viewer::TraceReadStatus::End,
                  "periodic fluid trace terminates cleanly");
            break;
        }
    }
    check(frames == steps,
          "periodic fluid trace contains exactly the accepted frames");
}

} // namespace

int main() {
    testDeterministicAcceptedFrames();
    testFrameRejectsUnacceptedState();
    testCompletedTrace();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d periodic fluid case check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all periodic fluid case checks passed");
    return 0;
}
