#include "scene_pressure_cell_case.h"
#include "scene_pressure_cell_checkpoint_persistence.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <ranges>
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

double norm(const fsi::StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

std::vector<std::uint8_t> serialized(
    const viewer::DiagnosticFrame& frame) {
    std::vector<std::uint8_t> bytes;
    viewer::ProtocolError error;
    check(viewer::serializeFrame(frame, bytes, &error),
          "scene pressure cell frame serializes");
    return bytes;
}

std::vector<std::uint8_t> serializedCheckpoint(
    const fsi::ScenePressureCellCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    fsi::ScenePressureCellCheckpointPersistenceError error;
    check(fsi::serializeScenePressureCellCheckpoint(
              checkpoint, bytes, &error),
          "scene pressure cell checkpoint serializes");
    return bytes;
}

const viewer::ScalarField* scalarField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.scalarFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.scalarFields.end() ? nullptr : &*found;
}

const viewer::VectorField* vectorField(
    const viewer::DiagnosticFrame& frame,
    const char* name) {
    const auto found = std::ranges::find_if(
        frame.vectorFields,
        [name](const auto& field) { return field.name == name; });
    return found == frame.vectorFields.end() ? nullptr : &*found;
}

void testVisibleStrongPressureCellAndReplay() {
    fsi::ScenePressureCellCase first;
    fsi::ScenePressureCellCase second;
    viewer::DiagnosticFrame frame;
    double peakDisplacement = 0.0;
    double peakPressure = 0.0;
    double peakPressureForce = 0.0;
    std::uint64_t peakIterations = 0;
    for (std::size_t step = 0; step < 120; ++step) {
        frame = first.advance();
        const auto repeated = second.advance();
        check(serialized(frame) == serialized(repeated),
              "scene pressure cell accepted frames are deterministic");
        const auto& diagnostics = first.diagnostics();
        peakDisplacement = std::max(
            peakDisplacement, diagnostics.maximumDisplacementMeters);
        peakPressure = std::max(
            peakPressure, diagnostics.maximumAbsolutePressurePascals);
        peakPressureForce = std::max(
            peakPressureForce, norm(diagnostics.pressureForceNewtons));
        peakIterations = std::max(
            peakIterations, diagnostics.coupling.solverRunCount);
        check(diagnostics.finite
                  && diagnostics.coupling.accepted
                  && diagnostics.coupling.iteration.status
                      == fsi::StrongCouplingIterationStatus::Converged
                  && diagnostics.coupling.pressureProjection.accepted
                  && diagnostics.coupling.pressureProjection
                         .correctedContinuityResidualMaximumCubicMetersPerSecond
                      < 2.0e-11
                  && diagnostics.coupling.pressureTransfer
                         .forceResidualNormNewtons < 1.0e-12
                  && diagnostics.coupling.pressureTransfer
                         .momentResidualNormNewtonMeters < 1.0e-12,
              "scene pressure cell closes every strong pressure-feedback step");
    }

    check(peakDisplacement > 1.0e-5
              && peakPressure > 1.0e-4
              && peakPressureForce > 1.0e-6
              && peakIterations >= 2,
          "visible cell develops motion, sparse pressure, and conservative feedback load");
    check(frame.step == 120
              && frame.vertices.size() == 4
              && frame.triangles.size() == 3
              && frame.sceneChecksum == fsi::scenePressureCellCaseChecksum
              && frame.solverCommit == fsi::scenePressureCellCaseSolverId,
          "scene pressure cell publishes its accepted scene-v2 shell");

    const auto* displacement = scalarField(
        frame, "pressure_cell.displacement");
    const auto* pressureJump = scalarField(
        frame, "pressure_cell.pressure_jump");
    const auto* maximumPressure = scalarField(
        frame, "pressure_cell.maximum_pressure");
    const auto* iterations = scalarField(
        frame, "pressure_cell.coupling_iterations");
    const auto* nodalForce = vectorField(
        frame, "pressure_cell.nodal_pressure_force");
    const auto* totalForce = vectorField(
        frame, "pressure_cell.total_pressure_force");
    const auto* actuator = vectorField(
        frame, "pressure_cell.actuator_force");
    check(displacement != nullptr
              && displacement->association
                  == viewer::FieldAssociation::Vertex
              && displacement->values.size() == frame.vertices.size()
              && pressureJump != nullptr
              && pressureJump->association
                  == viewer::FieldAssociation::Triangle
              && pressureJump->values.size() == frame.triangles.size()
              && maximumPressure != nullptr
              && maximumPressure->values.size() == 1
              && iterations != nullptr && iterations->values.size() == 1
              && nodalForce != nullptr
              && nodalForce->association
                  == viewer::FieldAssociation::Vertex
              && nodalForce->values.size() == frame.vertices.size()
              && totalForce != nullptr && totalForce->values.size() == 1
              && actuator != nullptr && actuator->values.size() == 1,
          "scene pressure cell frame exposes deformation, pressure, loads, and iteration count");

    const auto checkpoint = first.checkpoint();
    const auto expected = first.advance();
    const auto expectedDiagnostics = first.diagnostics();
    first.restore(checkpoint);
    const auto replay = first.advance();
    check(serialized(replay) == serialized(expected)
              && first.diagnostics() == expectedDiagnostics,
          "scene pressure cell checkpoint reproduces the exact next frame");
}

void testPersistentCheckpointAndRejection() {
    fsi::ScenePressureCellCase initial;
    const auto initialBytes = serializedCheckpoint(initial.checkpoint());
    fsi::ScenePressureCellCheckpoint initialDecoded;
    fsi::ScenePressureCellCheckpointPersistenceError error;
    check(fsi::deserializeScenePressureCellCheckpoint(
              initialBytes, initialDecoded, &error),
          "initial scene pressure cell checkpoint decodes without pressure state");
    fsi::ScenePressureCellCase initialReplay;
    initialReplay.restore(initialDecoded);
    check(serialized(initialReplay.advance()) == serialized(initial.advance()),
          "persisted initial scene pressure cell reproduces the first frame");

    fsi::ScenePressureCellCase source;
    for (std::size_t step = 0; step < 35; ++step) {
        static_cast<void>(source.advance());
    }
    const auto saved = source.checkpoint();
    const auto bytes = serializedCheckpoint(saved);
    check(!bytes.empty() && serializedCheckpoint(saved) == bytes,
          "scene pressure cell checkpoint encoding is deterministic");

    fsi::ScenePressureCellCheckpoint decoded;
    check(fsi::deserializeScenePressureCellCheckpoint(
              bytes, decoded, &error)
              && serializedCheckpoint(decoded) == bytes,
          "scene pressure cell checkpoint has a canonical bounded round trip");
    fsi::ScenePressureCellCase restored;
    restored.restore(decoded);
    const auto expected = source.advance();
    const auto replay = restored.advance();
    check(serialized(replay) == serialized(expected)
              && restored.diagnostics() == source.diagnostics(),
          "persisted scene pressure cell reproduces the exact next frame");

    const auto preserved = serializedCheckpoint(decoded);
    const auto rejects = [&](std::vector<std::uint8_t> damaged,
                             const fsi::ScenePressureCellCheckpointPersistenceErrorCode
                                 expectedCode,
                             const char* message,
                             const fsi::ScenePressureCellCheckpointPersistenceLimits&
                                 limits = {}) {
        fsi::ScenePressureCellCheckpointPersistenceError rejection;
        check(!fsi::deserializeScenePressureCellCheckpoint(
                  damaged, decoded, &rejection, limits)
                  && rejection.code == expectedCode
                  && serializedCheckpoint(decoded) == preserved,
              message);
    };

    auto damaged = bytes;
    damaged.front() ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidMagic,
            "scene pressure cell rejects foreign magic transactionally");
    damaged = bytes;
    damaged[8] ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::UnsupportedVersion,
            "scene pressure cell rejects an unsupported wire version transactionally");
    damaged = bytes;
    damaged.back() ^= 0xffU;
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::ChecksumMismatch,
            "scene pressure cell rejects payload corruption transactionally");
    damaged = bytes;
    damaged.pop_back();
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::Truncated,
            "scene pressure cell rejects truncation transactionally");
    damaged = bytes;
    damaged.push_back(0);
    rejects(std::move(damaged),
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::TrailingData,
            "scene pressure cell rejects trailing data transactionally");

    auto smallLimits =
        fsi::ScenePressureCellCheckpointPersistenceLimits{};
    smallLimits.maximumEncodedBytes = bytes.size() - 1;
    rejects(bytes,
            fsi::ScenePressureCellCheckpointPersistenceErrorCode::LimitExceeded,
            "scene pressure cell rejects an oversized checkpoint transactionally",
            smallLimits);

    auto foreign = saved;
    foreign.version += 1;
    std::vector<std::uint8_t> unchanged{1, 2, 3};
    const auto original = unchanged;
    check(!fsi::serializeScenePressureCellCheckpoint(
              foreign, unchanged, &error)
              && error.code
                  == fsi::ScenePressureCellCheckpointPersistenceErrorCode::InvalidData
              && unchanged == original,
          "scene pressure cell rejects foreign state before publishing bytes");
}

} // namespace

int main() {
    try {
        testVisibleStrongPressureCellAndReplay();
        testPersistentCheckpointAndRejection();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d scene pressure cell check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all scene pressure cell checks passed");
    return 0;
}
