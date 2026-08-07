#include "canonical_case.h"
#include "viewer_protocol.h"

#include <cmath>
#include <cstdio>
#include <sstream>
#include <vector>

namespace {

using namespace simwing;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(double actual,
               double expected,
               double tolerance,
               const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message,
                     actual,
                     expected);
        ++failures;
    }
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "deterministic frame serializes");
    return bytes;
}

void testDeterministicAcceptedStateAndFrames() {
    fsi::CanonicalStructuralCase withoutObserver;
    fsi::CanonicalStructuralCase withObserver;
    viewer::DiagnosticFrame firstFinal;
    viewer::DiagnosticFrame secondFinal;
    std::vector<std::vector<std::uint8_t>> observedFrames;
    constexpr std::uint64_t steps = 24;
    for (std::uint64_t step = 0; step < steps; ++step) {
        firstFinal = withoutObserver.advance();
        secondFinal = withObserver.advance();
        observedFrames.push_back(serialized(secondFinal));
        check(serialized(firstFinal) == observedFrames.back(),
              "frame observation cannot change numerical output");
    }

    const fsi::StructureCheckpoint first =
        withoutObserver.structure().checkpoint();
    const fsi::StructureCheckpoint second =
        withObserver.structure().checkpoint();
    check(first.acceptedStepCount == steps
              && first.acceptedStepCount == second.acceptedStepCount
              && first.simulationTimeSeconds == second.simulationTimeSeconds
              && first.nodes == second.nodes
              && first.pendingExternalForcesNewtons
                     == second.pendingExternalForcesNewtons
              && first.lastAppliedExternalForceNewtons
                     == second.lastAppliedExternalForceNewtons,
          "viewer-side frame consumption leaves committed Structure state identical");
    check(firstFinal.step == steps && secondFinal.step == steps,
          "frames are published only after accepted steps");
    check(firstFinal.sceneChecksum == fsi::canonicalCaseChecksum
              && firstFinal.solverCommit == fsi::canonicalCaseSolverId
              && firstFinal.couplingIteration == 0,
          "canonical frame provenance and iteration metadata are deterministic");
    checkNear(firstFinal.simulationTimeSeconds,
              steps * withoutObserver.stepSettings().timeStepSeconds,
              1.0e-15,
              "canonical frame time matches accepted fixed steps");
    check(firstFinal.vertices.size() == 3
              && firstFinal.triangles.size() == 1
              && firstFinal.lines.size() == 3,
          "canonical diagnostic topology is complete");
}

void testCompletedHeadlessTrace() {
    fsi::CanonicalStructuralCase simulation;
    std::stringstream trace(std::ios::in | std::ios::out | std::ios::binary);
    viewer::TraceWriter writer(trace);
    check(writer.writeHeader(simulation.traceHeader()),
          "headless trace header writes");
    constexpr std::uint64_t steps = 7;
    for (std::uint64_t step = 0; step < steps; ++step) {
        const viewer::DiagnosticFrame frame = simulation.advance();
        check(writer.writeFrame(frame),
              "accepted headless frame writes");
        trace.flush();
        check(static_cast<bool>(trace),
              "completed headless frame record flushes");
    }
    check(writer.finish(), "headless trace receives an explicit completion record");
    trace.flush();
    check(static_cast<bool>(trace), "completed headless trace flushes");

    trace.seekg(0);
    viewer::TraceReader reader(trace);
    viewer::TraceHeader header;
    check(reader.readHeader(header), "completed headless trace header reads");
    check(header.sceneChecksum == fsi::canonicalCaseChecksum
              && header.solverCommit == fsi::canonicalCaseSolverId,
          "headless trace retains canonical provenance");
    std::uint64_t frames = 0;
    viewer::DiagnosticFrame decoded;
    for (;;) {
        const viewer::TraceReadStatus status = reader.readNext(decoded);
        if (status == viewer::TraceReadStatus::Frame) {
            ++frames;
            check(decoded.step == frames,
                  "headless trace contains consecutive accepted-step frames");
        } else {
            check(status == viewer::TraceReadStatus::End,
                  "completed headless trace terminates cleanly");
            break;
        }
    }
    check(frames == steps,
          "completed headless trace contains exactly the requested frames");
}

} // namespace

int main() {
    testDeterministicAcceptedStateAndFrames();
    testCompletedHeadlessTrace();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d SimWing FSI worker check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all SimWing FSI worker checks passed");
    return 0;
}
