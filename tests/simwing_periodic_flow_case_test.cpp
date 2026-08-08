#include "periodic_flow_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
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

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [&](const viewer::ScalarField& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

const viewer::VectorField* vectorField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.vectorFields,
        [&](const viewer::VectorField& field) { return field.name == name; });
    return found == frame.vectorFields.end() ? nullptr : &*found;
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

void testAnalyticCellDiagnostics() {
    fsi::PeriodicFlowCase simulation;
    const viewer::PeriodicFluidCellFields fields =
        viewer::buildPeriodicFluidCellFields(
            simulation.grid(), simulation.velocity());
    check(fields.velocityMetersPerSecond.size()
              == simulation.grid().cellCount()
              && fields.speedMetersPerSecond.size()
                  == simulation.grid().cellCount()
              && fields.divergencePerSecond.size()
                  == simulation.grid().cellCount()
              && fields.vorticityPerSecond.size()
                  == simulation.grid().cellCount()
              && fields.vorticityMagnitudePerSecond.size()
                  == simulation.grid().cellCount(),
          "periodic cell diagnostics publish one owning value per cell");

    constexpr std::size_t i = 3;
    constexpr std::size_t j = 5;
    constexpr std::size_t k = 0;
    const auto spacing = simulation.grid().cellSpacingMeters();
    const auto center = simulation.grid().cellCenterMeters(i, j, k);
    const std::size_t index = simulation.grid().cellIndex(i, j, k);
    const double expectedX = 0.35
        + std::cos(0.5 * spacing.x)
            * std::sin(center.x) * std::cos(center.y);
    const double expectedY = -0.2
        - std::cos(0.5 * spacing.y)
            * std::cos(center.x) * std::sin(center.y);
    const double expectedVorticityZ = std::sin(center.x)
        * std::sin(center.y)
        * (std::cos(0.5 * spacing.y) * std::sin(spacing.x) / spacing.x
           + std::cos(0.5 * spacing.x) * std::sin(spacing.y) / spacing.y);
    checkNear(fields.velocityMetersPerSecond[index].x,
              expectedX, 1.0e-15,
              "periodic cell diagnostics average X MAC faces analytically");
    checkNear(fields.velocityMetersPerSecond[index].y,
              expectedY, 1.0e-15,
              "periodic cell diagnostics average Y MAC faces analytically");
    checkNear(fields.velocityMetersPerSecond[index].z,
              0.0, 0.0,
              "periodic cell diagnostics retain zero Z velocity");
    checkNear(fields.vorticityPerSecond[index].x,
              0.0, 0.0,
              "periodic Taylor-Green vorticity has zero X component");
    checkNear(fields.vorticityPerSecond[index].y,
              0.0, 0.0,
              "periodic Taylor-Green vorticity has zero Y component");
    checkNear(fields.vorticityPerSecond[index].z,
              expectedVorticityZ, 3.0e-15,
              "periodic cell curl matches the discrete Taylor-Green oracle");
    checkNear(fields.vorticityMagnitudePerSecond[index],
              std::abs(expectedVorticityZ), 3.0e-15,
              "periodic vorticity magnitude matches the vector field");
    check(fields.maximumAbsoluteDivergencePerSecond < 2.0e-14
              && fields.maximumVorticityPerSecond > 1.0,
          "periodic diagnostic extrema retain solenoidal vortical structure");
}

void testCellDiagnosticsRejectInvalidInputs() {
    fsi::PeriodicFlowCase simulation;
    const fsi::fluid::PeriodicCartesianGrid otherGrid(
        {4, 4, 2}, {}, {1.0, 1.0, 1.0});
    const fsi::fluid::MacVelocityField wrongShape(otherGrid);
    expectRejected(
        [&] { static_cast<void>(viewer::buildPeriodicFluidCellFields(
            simulation.grid(), wrongShape)); },
        "periodic cell diagnostics reject a mismatched MAC field");

    fsi::fluid::MacVelocityField nonFinite = simulation.velocity();
    nonFinite.xFaces()[0] = std::numeric_limits<double>::infinity();
    expectRejected(
        [&] { static_cast<void>(viewer::buildPeriodicFluidCellFields(
            simulation.grid(), nonFinite)); },
        "periodic cell diagnostics reject a non-finite MAC field");

    const fsi::fluid::MacVelocityField overflowing(
        simulation.grid(), std::numeric_limits<double>::max());
    expectRejected(
        [&] { static_cast<void>(viewer::buildPeriodicFluidCellFields(
            simulation.grid(), overflowing)); },
        "periodic cell diagnostics reject non-finite derived values");
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
    check(firstFrame.scalarFields.size() == 11
              && firstFrame.scalarFields[0].name == "pressure"
              && firstFrame.scalarFields[0].association
                  == viewer::FieldAssociation::Vertex
              && firstFrame.scalarFields[0].values.size()
                  == first.grid().cellCount()
              && firstFrame.scalarFields[2].name == "divergence"
              && firstFrame.scalarFields[3].name == "vorticity magnitude"
              && firstFrame.vectorFields.size() == 3
              && firstFrame.vectorFields[0].name == "velocity"
              && firstFrame.vectorFields[0].values.size()
                  == first.grid().cellCount()
              && firstFrame.vectorFields[1].name == "vorticity"
              && firstFrame.vectorFields[1].values.size()
                  == first.grid().cellCount(),
          "periodic fluid frame retains scalar/vector flow diagnostics");
    check(std::ranges::any_of(
              firstFrame.scalarFields[0].values,
              [](const double pressure) { return pressure != 0.0; })
              && first.diagnostics().accepted
              && firstFrame.couplingResiduals.fluid
                  == first.diagnostics().finalDivergenceL2PerSecond,
          "periodic fluid frame comes from an accepted nontrivial pressure solve");

    const viewer::ScalarField* divergence =
        scalarField(firstFrame, "divergence");
    const viewer::ScalarField* vorticityMagnitude =
        scalarField(firstFrame, "vorticity magnitude");
    const viewer::ScalarField* maximumDivergence =
        scalarField(firstFrame, "divergence maximum absolute");
    const viewer::ScalarField* maximumVorticity =
        scalarField(firstFrame, "vorticity maximum");
    const viewer::VectorField* vorticity =
        vectorField(firstFrame, "vorticity");
    fsi::fluid::CellScalarField expectedDivergence(first.grid());
    fsi::fluid::computeDivergence(
        first.grid(), first.velocity(), expectedDivergence);
    check(divergence != nullptr
              && std::ranges::equal(
                  divergence->values, expectedDivergence.values())
              && fsi::fluid::l2Norm(expectedDivergence)
                  == first.diagnostics().finalDivergenceL2PerSecond,
          "periodic frame divergence is the exact accepted MAC operator");
    bool magnitudeMatches = vorticity != nullptr
        && vorticityMagnitude != nullptr
        && vorticity->values.size() == vorticityMagnitude->values.size();
    if (magnitudeMatches) {
        for (std::size_t index = 0;
             index < vorticity->values.size(); ++index) {
            const viewer::Vec3d value = vorticity->values[index];
            magnitudeMatches = magnitudeMatches
                && vorticityMagnitude->values[index]
                    == std::hypot(value.x, value.y, value.z);
        }
    }
    check(magnitudeMatches,
          "periodic frame vorticity magnitude derives exactly from curl vectors");
    const double expectedMaximumVorticity = vorticityMagnitude == nullptr
        ? 0.0
        : *std::ranges::max_element(vorticityMagnitude->values);
    check(maximumDivergence != nullptr
              && maximumDivergence->values[0]
                  == fsi::fluid::maximumAbsoluteValue(expectedDivergence)
              && maximumVorticity != nullptr
              && maximumVorticity->values[0] == expectedMaximumVorticity,
          "periodic frame publishes exact divergence/vorticity extrema");
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
    testAnalyticCellDiagnostics();
    testCellDiagnosticsRejectInvalidInputs();
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
