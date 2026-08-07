#include "viewer_protocol.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace simwing::viewer;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

DiagnosticFrame sampleFrame() {
    DiagnosticFrame frame;
    frame.sceneChecksum = "sha256:scene-fixture";
    frame.solverCommit = "0123456789abcdef-dirty";
    frame.step = 42;
    frame.simulationTimeSeconds = 0.7;
    frame.timeStepSeconds = 1.0 / 60.0;
    frame.couplingIteration = 3;
    frame.couplingResiduals = {1.25e-6, 2.5e-4, 3.0e-7, 4.0e-8,
                               -5.0e-5};
    frame.conservation.fluidMassKilograms = 12.5;
    frame.conservation.totalMomentumNewtonSeconds = {1.0, 2.0, 3.0};
    frame.conservation.totalEnergyJoules = 145.0;
    frame.conservation.interfaceForceResidualNewtons = {0.01, -0.02, 0.03};
    frame.conservation.interfaceMomentResidualNewtonMetres = {0.04, 0.05,
                                                               -0.06};
    frame.conservation.interfacePowerResidualWatts = -0.007;

    frame.vertices = {{100, {0.0, 0.0, 0.0}},
                      {101, {1.0, 0.0, 0.0}},
                      {102, {0.0, 1.0, 0.1}}};
    frame.triangles = {{200, 0, 1, 2, 1, 17}};
    frame.lines = {{300, 0, 2, 4}};
    frame.contacts = {{400,
                       {EntityKind::Vertex, 0},
                       {EntityKind::Triangle, 0},
                       {0.1, 0.2, 0.3},
                       {0.0, 0.0, 1.0},
                       -0.001}};
    frame.sealing = {{500,
                      1,
                      17,
                      {0.2, 0.3, 0.4},
                      {0.0, 1.0, 0.0},
                      0.0002,
                      true}};
    frame.scalarFields = {{"strain", "1", FieldAssociation::Triangle,
                           {0.012}},
                          {"line tension", "N", FieldAssociation::Line,
                           {123.0}},
                          {"CFL", "1", FieldAssociation::Global, {0.45}}};
    frame.vectorFields = {
        {"velocity",
         "m/s",
         FieldAssociation::Vertex,
         {{1.0, 0.0, 0.0}, {1.1, 0.0, 0.0}, {0.9, 0.1, 0.0}}},
        {"traction",
         "Pa",
         FieldAssociation::Triangle,
         {{0.0, 0.0, -18.0}}}};
    return frame;
}

void stampLittleDouble(
    std::vector<std::uint8_t>& bytes,
    std::size_t offset,
    double value) {
    std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    for (std::size_t i = 0; i < sizeof(bits); ++i) {
        bytes[offset + i] = static_cast<std::uint8_t>(bits & 0xffU);
        bits >>= 8U;
    }
}

void testDeterministicRoundTrip() {
    const DiagnosticFrame original = sampleFrame();
    ProtocolError error;
    std::vector<std::uint8_t> first;
    std::vector<std::uint8_t> second;
    check(serializeFrame(original, first, &error),
          "round trip: initial serialization succeeds");
    check(!error, "round trip: successful serialization clears error");
    check(serializeFrame(original, second, &error),
          "round trip: repeated serialization succeeds");
    check(first == second,
          "round trip: repeated serialization is byte deterministic");

    DiagnosticFrame decoded;
    check(deserializeFrame(first, decoded, &error),
          "round trip: deserialization succeeds");
    std::vector<std::uint8_t> reencoded;
    check(serializeFrame(decoded, reencoded, &error),
          "round trip: decoded frame can be serialized");
    check(first == reencoded,
          "round trip: decode and re-encode preserves every byte");
    check(decoded.sceneChecksum == original.sceneChecksum
              && decoded.solverCommit == original.solverCommit
              && decoded.step == original.step
              && decoded.vertices.size() == 3
              && decoded.triangles.size() == 1
              && decoded.lines.size() == 1
              && decoded.contacts.size() == 1
              && decoded.sealing.size() == 1
              && decoded.scalarFields.size() == 3
              && decoded.vectorFields.size() == 2,
          "round trip: frame content and provenance survive");
}

void testMultiFrameTrace() {
    DiagnosticFrame first = sampleFrame();
    DiagnosticFrame second = sampleFrame();
    second.step = 43;
    second.simulationTimeSeconds += second.timeStepSeconds;
    second.vertices[2].positionMetres.z = 0.15;
    second.scalarFields[0].values[0] = 0.02;

    std::ostringstream output(std::ios::binary);
    TraceWriter writer(output);
    const TraceHeader header{first.sceneChecksum, first.solverCommit};
    check(writer.writeHeader(header), "trace: header write succeeds");
    check(writer.writeFrame(first), "trace: first frame write succeeds");
    check(writer.writeFrame(second), "trace: second frame write succeeds");

    const std::string trace = output.str();
    std::istringstream input(trace, std::ios::binary);
    TraceReader reader(input);
    TraceHeader decodedHeader;
    check(reader.readHeader(decodedHeader), "trace: header read succeeds");
    check(decodedHeader.sceneChecksum == header.sceneChecksum
              && decodedHeader.solverCommit == header.solverCommit,
          "trace: header provenance survives");
    DiagnosticFrame frame;
    check(reader.readNext(frame) == TraceReadStatus::Frame && frame.step == 42,
          "trace: first frame replays in order");
    check(reader.readNext(frame) == TraceReadStatus::Frame && frame.step == 43
              && frame.vertices[2].positionMetres.z == 0.15,
          "trace: second frame replays in order");
    check(reader.readNext(frame) == TraceReadStatus::End,
          "trace: clean record boundary reports end of trace");

    std::ostringstream wrongOutput(std::ios::binary);
    TraceWriter wrongWriter(wrongOutput);
    check(wrongWriter.writeHeader(header),
          "trace provenance: header write succeeds");
    second.solverCommit = "another-commit";
    check(!wrongWriter.writeFrame(second)
              && wrongWriter.error().code == ProtocolErrorCode::InvalidData,
          "trace provenance: mismatched frame is rejected");
}

void testFrameCorruptionAndTruncation() {
    ProtocolError error;
    std::vector<std::uint8_t> bytes;
    check(serializeFrame(sampleFrame(), bytes, &error),
          "corruption: fixture serialization succeeds");

    std::vector<std::uint8_t> corrupt = bytes;
    corrupt[0] ^= 0xffU;
    DiagnosticFrame frame;
    check(!deserializeFrame(corrupt, frame, &error)
              && error.code == ProtocolErrorCode::InvalidMagic,
          "corruption: bad frame magic is rejected");

    corrupt = bytes;
    corrupt[4] = static_cast<std::uint8_t>(kFrameProtocolVersion + 1);
    check(!deserializeFrame(corrupt, frame, &error)
              && error.code == ProtocolErrorCode::UnsupportedVersion,
          "corruption: unsupported frame version is rejected");

    corrupt = bytes;
    corrupt.pop_back();
    frame = sampleFrame();
    check(!deserializeFrame(corrupt, frame, &error)
              && error.code == ProtocolErrorCode::Truncated,
          "corruption: truncated frame is rejected");
    check(frame.step == 42 && frame.vertices.size() == 3,
          "corruption: failed decode preserves the previous frame");

    corrupt = bytes;
    corrupt.push_back(0);
    check(!deserializeFrame(corrupt, frame, &error)
              && error.code == ProtocolErrorCode::TrailingData,
          "corruption: trailing frame bytes are rejected");

    corrupt = bytes;
    const std::size_t timeOffset = 16
                                   + 4 + sampleFrame().sceneChecksum.size()
                                   + 4 + sampleFrame().solverCommit.size() + 8;
    stampLittleDouble(corrupt, timeOffset,
                      std::numeric_limits<double>::quiet_NaN());
    check(!deserializeFrame(corrupt, frame, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "corruption: encoded non-finite values are rejected");
}

void testValidationAndLimits() {
    ProtocolError error;
    std::vector<std::uint8_t> bytes;
    DiagnosticFrame frame = sampleFrame();
    frame.triangles[0].vertex2 = 99;
    check(!serializeFrame(frame, bytes, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "validation: invalid triangle index is rejected");

    frame = sampleFrame();
    frame.vectorFields[0].values[1].x =
        std::numeric_limits<double>::infinity();
    check(!serializeFrame(frame, bytes, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "validation: non-finite field value is rejected");

    frame = sampleFrame();
    frame.scalarFields[0].values.push_back(0.0);
    check(!serializeFrame(frame, bytes, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "validation: field/association count mismatch is rejected");

    frame = sampleFrame();
    frame.vertices[0].stableId = 0;
    check(!serializeFrame(frame, bytes, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "validation: zero stable IDs are rejected");

    frame = sampleFrame();
    frame.triangles[0].positiveRegionId =
        frame.triangles[0].negativeRegionId;
    check(!serializeFrame(frame, bytes, &error)
              && error.code == ProtocolErrorCode::InvalidData,
          "validation: identical triangle side regions are rejected");

    frame = sampleFrame();
    ProtocolLimits tight;
    tight.maxVertices = 2;
    check(!serializeFrame(frame, bytes, &error, tight)
              && error.code == ProtocolErrorCode::LimitExceeded,
          "limits: writer rejects an unreasonable vertex count");

    check(serializeFrame(frame, bytes, &error),
          "limits: fixture serializes under default limits");
    DiagnosticFrame decoded;
    check(!deserializeFrame(bytes, decoded, &error, tight)
              && error.code == ProtocolErrorCode::LimitExceeded,
          "limits: reader rejects an unreasonable encoded vertex count");

    tight = {};
    tight.maxFrameBytes = 128;
    check(!serializeFrame(frame, bytes, &error, tight)
              && error.code == ProtocolErrorCode::LimitExceeded,
          "limits: writer enforces encoded byte cap while encoding");

    check(serializeFrame(frame, bytes, &error),
          "limits: fixture is restored for reader byte-cap test");
    check(!deserializeFrame(bytes, decoded, &error, tight)
              && error.code == ProtocolErrorCode::LimitExceeded,
          "limits: reader enforces encoded byte cap before decoding");
}

void testTraceCorruptionAndTruncation() {
    const DiagnosticFrame frame = sampleFrame();
    std::ostringstream output(std::ios::binary);
    TraceWriter writer(output);
    check(writer.writeHeader({frame.sceneChecksum, frame.solverCommit})
              && writer.writeFrame(frame),
          "trace corruption: fixture trace writes");
    const std::string trace = output.str();

    std::string corrupt = trace;
    corrupt[0] = 'X';
    std::istringstream badMagic(corrupt, std::ios::binary);
    TraceReader magicReader(badMagic);
    TraceHeader header;
    check(!magicReader.readHeader(header)
              && magicReader.error().code == ProtocolErrorCode::InvalidMagic,
          "trace corruption: bad trace magic is rejected");

    corrupt = trace;
    corrupt[4] = static_cast<char>(kTraceProtocolVersion + 1);
    std::istringstream badVersion(corrupt, std::ios::binary);
    TraceReader versionReader(badVersion);
    check(!versionReader.readHeader(header)
              && versionReader.error().code
                     == ProtocolErrorCode::UnsupportedVersion,
          "trace corruption: unsupported trace version is rejected");

    corrupt = trace.substr(0, trace.size() - 5);
    std::istringstream truncated(corrupt, std::ios::binary);
    TraceReader truncatedReader(truncated);
    DiagnosticFrame decoded;
    check(truncatedReader.readHeader(header),
          "trace truncation: intact header reads");
    check(truncatedReader.readNext(decoded) == TraceReadStatus::Error
              && truncatedReader.error().code == ProtocolErrorCode::Truncated,
          "trace truncation: incomplete frame is rejected");
}

} // namespace

int main() {
    testDeterministicRoundTrip();
    testMultiFrameTrace();
    testFrameCorruptionAndTruncation();
    testValidationAndLimits();
    testTraceCorruptionAndTruncation();
    if (failures == 0) {
        std::printf("simwing viewer protocol tests passed\n");
    }
    return failures == 0 ? 0 : 1;
}
