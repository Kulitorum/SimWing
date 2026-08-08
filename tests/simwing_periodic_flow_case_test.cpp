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

std::vector<std::uint8_t> serializedCheckpoint(
    const fsi::PeriodicFlowCaseCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    fsi::PeriodicFlowCaseCheckpointError error;
    check(fsi::serializePeriodicFlowCaseCheckpoint(
              checkpoint, bytes, &error),
          "periodic worker checkpoint serializes");
    check(!error,
          "successful periodic checkpoint serialization clears its error");
    return bytes;
}

void refreshCheckpointChecksum(std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 28;
    constexpr std::size_t checksumOffset = 20;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t checksum = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        checksum ^= bytes[index];
        checksum *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(checksum); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            checksum >> (8U * byte));
    }
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

void testCheckpointValidationAndReplay() {
    fsi::PeriodicFlowCase initialWorker;
    const auto initial = initialWorker.checkpoint();
    const auto expectedFirst = initialWorker.advance();
    initialWorker.restore(initial);
    const auto restoredFirst = initialWorker.advance();
    check(initial.version == fsi::periodicFlowCaseCheckpointVersion
              && initial.caseDefinitionFingerprint
                  == fsi::periodicFlowCaseDefinitionFingerprint
              && initial.acceptedStepCount == 0
              && initial.simulationTimeSeconds == 0.0
              && serialized(expectedFirst) == serialized(restoredFirst),
          "periodic checkpoint restores and replays the initial state");

    fsi::PeriodicFlowCase simulation;
    constexpr std::uint64_t savedSteps = 9;
    for (std::uint64_t step = 0; step < savedSteps; ++step) {
        static_cast<void>(simulation.advance());
    }
    const auto saved = simulation.checkpoint();
    const auto savedVelocity = simulation.velocity();
    const auto savedPressure = simulation.pressure();
    const auto savedDiagnostics = simulation.diagnostics();
    const auto expectedNext = simulation.advance();
    for (std::uint64_t step = 0; step < 4; ++step) {
        static_cast<void>(simulation.advance());
    }
    simulation.restore(saved);
    check(saved.cellCounts == simulation.grid().cellCounts()
              && saved.lowerMeters == simulation.grid().lowerMeters()
              && saved.upperMeters == simulation.grid().upperMeters()
              && saved.scalarSampleCount == simulation.grid().cellCount()
              && saved.acceptedStepCount == savedSteps
              && saved.simulationTimeSeconds
                  == simulation.simulationTimeSeconds()
              && simulation.velocity() == savedVelocity
              && simulation.pressure() == savedPressure
              && simulation.diagnostics() == savedDiagnostics,
          "periodic checkpoint restores fields, diagnostics, step, and time");
    const auto restoredNext = simulation.advance();
    check(serialized(expectedNext) == serialized(restoredNext),
          "periodic checkpoint continuation replays bit-for-bit");

    fsi::PeriodicFlowCase rebuiltWorker;
    rebuiltWorker.restore(saved);
    const auto rebuiltNext = rebuiltWorker.advance();
    check(serialized(expectedNext) == serialized(rebuiltNext),
          "equivalent rebuilt periodic worker resumes bit-for-bit");

    const auto beforeInvalidVelocity = simulation.velocity();
    const auto beforeInvalidPressure = simulation.pressure();
    const auto beforeInvalidDiagnostics = simulation.diagnostics();
    const std::uint64_t beforeInvalidStep = simulation.acceptedStepCount();
    const double beforeInvalidTime = simulation.simulationTimeSeconds();
    auto wrongVersion = saved;
    ++wrongVersion.version;
    expectRejected(
        [&] { simulation.restore(wrongVersion); },
        "periodic checkpoint rejects an unsupported version");
    auto wrongDefinition = saved;
    ++wrongDefinition.caseDefinitionFingerprint;
    expectRejected(
        [&] { simulation.restore(wrongDefinition); },
        "periodic checkpoint rejects a foreign case definition");
    auto wrongStep = saved;
    ++wrongStep.acceptedStepCount;
    expectRejected(
        [&] { simulation.restore(wrongStep); },
        "periodic checkpoint rejects corrupted step metadata");
    auto wrongTime = saved;
    wrongTime.simulationTimeSeconds += 1.0;
    expectRejected(
        [&] { simulation.restore(wrongTime); },
        "periodic checkpoint rejects corrupted time metadata");
    auto wrongGrid = saved;
    ++wrongGrid.cellCounts.x;
    expectRejected(
        [&] { simulation.restore(wrongGrid); },
        "periodic checkpoint rejects corrupted grid metadata");
    auto wrongSampleCount = saved;
    ++wrongSampleCount.scalarSampleCount;
    expectRejected(
        [&] { simulation.restore(wrongSampleCount); },
        "periodic checkpoint rejects corrupted sample metadata");
    const fsi::PeriodicFlowCaseCheckpoint empty;
    expectRejected(
        [&] { simulation.restore(empty); },
        "periodic checkpoint rejects a missing immutable payload");
    check(simulation.velocity() == beforeInvalidVelocity
              && simulation.pressure() == beforeInvalidPressure
              && simulation.diagnostics() == beforeInvalidDiagnostics
              && simulation.acceptedStepCount() == beforeInvalidStep
              && simulation.simulationTimeSeconds() == beforeInvalidTime,
          "rejected periodic checkpoint restores leave the worker untouched");
}

void testCheckpointBinaryProtocol() {
    fsi::PeriodicFlowCase simulation;
    constexpr std::uint64_t savedSteps = 9;
    for (std::uint64_t step = 0; step < savedSteps; ++step) {
        static_cast<void>(simulation.advance());
    }
    const auto saved = simulation.checkpoint();
    const auto firstEncoding = serializedCheckpoint(saved);
    const auto secondEncoding = serializedCheckpoint(saved);
    check(firstEncoding == secondEncoding && !firstEncoding.empty(),
          "periodic checkpoint encoding is byte deterministic");

    fsi::PeriodicFlowCaseCheckpointError error;
    fsi::PeriodicFlowCaseCheckpoint decoded;
    check(fsi::deserializePeriodicFlowCaseCheckpoint(
              firstEncoding, decoded, &error)
              && !error
              && serializedCheckpoint(decoded) == firstEncoding,
          "periodic checkpoint decode and re-encode preserves every byte");
    fsi::PeriodicFlowCase resumed;
    resumed.restore(decoded);
    check(resumed.acceptedStepCount() == savedSteps
              && resumed.simulationTimeSeconds()
                  == simulation.simulationTimeSeconds()
              && resumed.velocity() == simulation.velocity()
              && resumed.pressure() == simulation.pressure()
              && resumed.diagnostics() == simulation.diagnostics(),
          "decoded periodic checkpoint restores every committed field");
    const auto expectedNext = simulation.advance();
    const auto resumedNext = resumed.advance();
    check(serialized(expectedNext) == serialized(resumedNext),
          "decoded periodic checkpoint continues bit-for-bit");

    fsi::PeriodicFlowCase initialWorker;
    const auto initialBytes = serializedCheckpoint(initialWorker.checkpoint());
    fsi::PeriodicFlowCaseCheckpoint decodedInitial;
    check(fsi::deserializePeriodicFlowCaseCheckpoint(
              initialBytes, decodedInitial, &error),
          "initial periodic checkpoint decodes");
    const auto expectedInitialFirst = initialWorker.advance();
    fsi::PeriodicFlowCase resumedInitial;
    resumedInitial.restore(decodedInitial);
    check(serialized(expectedInitialFirst)
              == serialized(resumedInitial.advance()),
          "decoded initial periodic checkpoint replays its first step");

    fsi::PeriodicFlowCase preservedWorker;
    static_cast<void>(preservedWorker.advance());
    fsi::PeriodicFlowCaseCheckpoint output = preservedWorker.checkpoint();
    const auto preservedEncoding = serializedCheckpoint(output);
    const auto expectDecodeRejected = [&](
        const std::vector<std::uint8_t>& candidate,
        const fsi::PeriodicFlowCaseCheckpointErrorCode expectedCode,
        const fsi::PeriodicFlowCaseCheckpointLimits& limits,
        const char* message) {
        error = {};
        check(!fsi::deserializePeriodicFlowCaseCheckpoint(
                  candidate, output, &error, limits)
                  && error.code == expectedCode
                  && serializedCheckpoint(output) == preservedEncoding,
              message);
    };

    auto corrupt = firstEncoding;
    corrupt[0] ^= 0xffU;
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::InvalidMagic,
        {},
        "periodic checkpoint rejects bad magic without changing output");
    corrupt = firstEncoding;
    ++corrupt[8];
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::UnsupportedVersion,
        {},
        "periodic checkpoint rejects unsupported wire versions");
    corrupt = firstEncoding;
    ++corrupt[28];
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::UnsupportedVersion,
        {},
        "periodic checkpoint rejects unsupported state versions");
    corrupt = firstEncoding;
    corrupt[32] ^= 0x01U;
    refreshCheckpointChecksum(corrupt);
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::InvalidData,
        {},
        "periodic checkpoint rejects foreign case fingerprints");
    corrupt = firstEncoding;
    corrupt[10] = 1;
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::InvalidData,
        {},
        "periodic checkpoint rejects nonzero reserved bits");
    corrupt = firstEncoding;
    corrupt.pop_back();
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::Truncated,
        {},
        "periodic checkpoint rejects truncation transactionally");
    corrupt = firstEncoding;
    corrupt.push_back(0);
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::TrailingData,
        {},
        "periodic checkpoint rejects trailing data transactionally");
    corrupt = firstEncoding;
    corrupt.back() ^= 0x80U;
    expectDecodeRejected(
        corrupt,
        fsi::PeriodicFlowCaseCheckpointErrorCode::ChecksumMismatch,
        {},
        "periodic checkpoint detects payload corruption");

    fsi::PeriodicFlowCaseCheckpointLimits byteLimit;
    byteLimit.maximumBytes = firstEncoding.size() - 1;
    expectDecodeRejected(
        firstEncoding,
        fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
        byteLimit,
        "periodic checkpoint enforces its configured byte limit");
    fsi::PeriodicFlowCaseCheckpointLimits sampleLimit;
    sampleLimit.maximumScalarSamples = simulation.grid().cellCount() - 1;
    expectDecodeRejected(
        firstEncoding,
        fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
        sampleLimit,
        "periodic checkpoint enforces its configured sample limit");
    fsi::PeriodicFlowCaseCheckpointLimits substepLimit;
    substepLimit.maximumSubsteps = 0;
    expectDecodeRejected(
        firstEncoding,
        fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
        substepLimit,
        "periodic checkpoint enforces its configured substep limit");

    std::vector<std::uint8_t> rejectedBytes{1, 2, 3};
    check(!fsi::serializePeriodicFlowCaseCheckpoint(
              fsi::PeriodicFlowCaseCheckpoint{}, rejectedBytes, &error)
              && error.code
                  == fsi::PeriodicFlowCaseCheckpointErrorCode::InvalidData
              && rejectedBytes.empty(),
          "periodic checkpoint rejects missing payloads transactionally");
    rejectedBytes = {1, 2, 3};
    check(!fsi::serializePeriodicFlowCaseCheckpoint(
              saved, rejectedBytes, &error, byteLimit)
              && error.code
                  == fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded
              && rejectedBytes.empty(),
          "periodic checkpoint encoding honors its byte limit");
    rejectedBytes = {1, 2, 3};
    check(!fsi::serializePeriodicFlowCaseCheckpoint(
              saved, rejectedBytes, &error, sampleLimit)
              && error.code
                  == fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded
              && rejectedBytes.empty(),
          "periodic checkpoint encoding honors its sample limit");
    rejectedBytes = {1, 2, 3};
    check(!fsi::serializePeriodicFlowCaseCheckpoint(
              saved, rejectedBytes, &error, substepLimit)
              && error.code
                  == fsi::PeriodicFlowCaseCheckpointErrorCode::LimitExceeded
              && rejectedBytes.empty(),
          "periodic checkpoint encoding honors its substep limit");
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
    testCheckpointValidationAndReplay();
    testCheckpointBinaryProtocol();
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
