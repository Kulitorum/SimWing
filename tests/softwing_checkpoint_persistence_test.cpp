#include <softwing/checkpoint_persistence.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <span>
#include <vector>

namespace {

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

softwing::OrthotropicMembraneMaterial material(
    const double stiffnessScale = 1.0) {
    softwing::OrthotropicMembraneMaterial result;
    result.warpStiffness = 800.0 * stiffnessScale;
    result.weftStiffness = 500.0;
    result.couplingStiffness = 100.0;
    result.shearStiffness = 180.0;
    result.warpPreTension = 0.02;
    result.weftPreTension = 0.01;
    result.dampingTime = 0.0;
    result.compressionStiffnessRatio = 0.2;
    return result;
}

softwing::SoftBody makeBody(const double stiffnessScale = 1.0) {
    softwing::SoftBody body;
    body.addFixedNode({0.0, 0.0, 2.0});
    body.addNode({1.08, 0.0, 2.0}, 1.0);
    body.addNode({0.0, 1.0, 2.0}, 1.0);
    body.addNode({1.0, -1.0, 2.25}, 1.0);
    body.addNode({-1.0, -1.0, 0.0}, 1.0);
    body.addNode({1.0, -1.0, 0.0}, 1.0);
    body.addNode({0.0, 1.0, 0.0}, 1.0);
    body.addNode({-1.0, -1.0, 0.01}, 1.0);
    body.addNode({0.0, 1.0, 0.01}, 1.0);
    body.addNode({1.0, -1.0, 0.01}, 1.0);

    const std::size_t membraneTriangle = body.addTriangle(0, 1, 2);
    body.addTriangle(1, 0, 3);
    body.addTriangle(4, 5, 6);
    body.addTriangle(7, 8, 9);
    body.addDistanceConstraint(0, 1, 1.0, 2.0e-5);
    body.addDihedralBendingConstraint(0, 1, 2, 3, 0.0, 1.0e-5);
    const softwing::MembraneElementDefinition membrane{
        membraneTriangle,
        {{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}},
        material(stiffnessScale),
        softwing::MaterialRole::Bulk,
    };
    static_cast<void>(body.addMembraneElements(
        std::span<const softwing::MembraneElementDefinition>{&membrane, 1}));

    const softwing::ContactSurfaceHandle contact = body.addContactSurface(
        body.surfaceGroup(2, 2), 0.01);
    static_cast<void>(body.addContactPair(
        contact.collider(), contact.collider(), {0.0, 0.4, 0.25}));
    body.setFacePressureDifference(0, 12.5);
    body.setFacePressureDifference(1, -3.0);
    body.addForce(1, {0.5, -0.25, 0.75});
    return body;
}

softwing::StepSettings settings() {
    softwing::StepSettings result;
    result.timeStep = 1.0 / 120.0;
    result.substeps = 2;
    result.constraintIterations = 4;
    result.gravity = {0.0, 0.0, -0.2};
    result.velocityDampingPerSecond = 0.05;
    return result;
}

std::vector<std::uint8_t> serialize(
    const softwing::SoftBodyCheckpoint& checkpoint) {
    std::vector<std::uint8_t> bytes;
    softwing::SoftBodyCheckpointPersistenceError error;
    check(softwing::serializeSoftBodyCheckpoint(
              checkpoint, bytes, &error),
          "valid body checkpoint serializes");
    return bytes;
}

std::uint64_t readU64(const std::vector<std::uint8_t>& bytes,
                      const std::size_t offset) {
    std::uint64_t result = 0;
    for (std::size_t index = 0; index < 8; ++index) {
        result |= static_cast<std::uint64_t>(bytes[offset + index])
            << (8U * index);
    }
    return result;
}

void writeU64(std::vector<std::uint8_t>& bytes,
              const std::size_t offset,
              std::uint64_t value) {
    for (std::size_t index = 0; index < 8; ++index) {
        bytes[offset + index] = static_cast<std::uint8_t>(value & 0xffU);
        value >>= 8U;
    }
}

std::uint64_t checksum(
    const std::span<const std::uint8_t> bytes) noexcept {
    std::uint64_t result = 14695981039346656037ULL;
    for (const std::uint8_t value : bytes) {
        result ^= value;
        result *= 1099511628211ULL;
    }
    return result;
}

void refreshChecksum(std::vector<std::uint8_t>& bytes) {
    const std::size_t payloadSize =
        static_cast<std::size_t>(readU64(bytes, 8));
    writeU64(bytes, 16, checksum(std::span<const std::uint8_t>{bytes}
                                    .subspan(24, payloadSize)));
}

void testRoundTripAndContinuation() {
    softwing::SoftBody source = makeBody();
    source.step(settings());
    check(!source.contactRecords().empty(),
          "fixture commits active contact records");
    const softwing::SoftBodyCheckpoint saved = source.checkpoint();
    const auto first = serialize(saved);
    const auto second = serialize(saved);
    check(first == second,
          "checkpoint encoding is byte deterministic");

    softwing::SoftBody rebuilt = makeBody();
    const auto topologyTemplate = rebuilt.checkpoint();
    softwing::SoftBodyCheckpoint decoded;
    softwing::SoftBodyCheckpointPersistenceError error;
    const bool decodedOk = softwing::deserializeSoftBodyCheckpoint(
        first, topologyTemplate, decoded, &error);
    if (!decodedOk) {
        std::fprintf(stderr, "decode error %d: %s\n",
                     static_cast<int>(error.code), error.message.c_str());
    }
    check(decodedOk,
          "checkpoint decodes against equivalent rebuilt topology");
    if (!decodedOk) return;
    check(serialize(decoded) == first,
          "decoded checkpoint re-encodes byte identically");
    rebuilt.restore(decoded);

    source.step(settings());
    rebuilt.step(settings());
    check(serialize(source.checkpoint())
              == serialize(rebuilt.checkpoint()),
          "decoded body state continues bit-for-bit including warm starts");

    softwing::SoftBody emptyBody;
    emptyBody.addNode({}, 1.0);
    const auto emptyBytes = serialize(emptyBody.checkpoint());
    softwing::SoftBodyCheckpoint emptyDecoded;
    check(softwing::deserializeSoftBodyCheckpoint(
              emptyBytes, emptyBody.checkpoint(), emptyDecoded, &error),
          "default positive-infinity contact gap round-trips");
}

void testTopologyBinding() {
    softwing::SoftBody source = makeBody();
    source.step(settings());
    const auto bytes = serialize(source.checkpoint());

    softwing::SoftBody foreign = makeBody(1.01);
    softwing::SoftBodyCheckpoint output = foreign.checkpoint();
    const auto preserved = serialize(output);
    softwing::SoftBodyCheckpointPersistenceError error;
    check(!softwing::deserializeSoftBodyCheckpoint(
              bytes, foreign.checkpoint(), output, &error)
              && error.code
                  == softwing::SoftBodyCheckpointPersistenceErrorCode::
                      TopologyMismatch
              && serialize(output) == preserved,
          "foreign materials are rejected transactionally by topology");
}

void testCorruptionAndLimits() {
    softwing::SoftBody body = makeBody();
    body.step(settings());
    const auto valid = serialize(body.checkpoint());
    const auto topologyTemplate = makeBody().checkpoint();
    softwing::SoftBodyCheckpoint preserved = topologyTemplate;
    const auto preservedBytes = serialize(preserved);

    const auto reject = [&](std::vector<std::uint8_t> bytes,
                            const softwing::SoftBodyCheckpointPersistenceErrorCode
                                expected,
                            const char* message,
                            const softwing::SoftBodyCheckpointPersistenceLimits&
                                limits = {}) {
        softwing::SoftBodyCheckpointPersistenceError error;
        check(!softwing::deserializeSoftBodyCheckpoint(
                  bytes, topologyTemplate, preserved, &error, limits)
                  && error.code == expected
                  && serialize(preserved) == preservedBytes,
              message);
    };

    auto corrupt = valid;
    corrupt[0] ^= 0xffU;
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::InvalidMagic,
           "bad magic is rejected transactionally");

    corrupt = valid;
    corrupt[4] = 2;
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::
               UnsupportedVersion,
           "unsupported envelope version is rejected");

    corrupt = valid;
    corrupt[6] = 1;
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::InvalidData,
           "nonzero reserved envelope bits are rejected");

    corrupt = valid;
    corrupt.pop_back();
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::Truncated,
           "truncated payload is rejected transactionally");

    corrupt = valid;
    corrupt.push_back(0);
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::TrailingData,
           "trailing envelope data is rejected transactionally");

    corrupt = valid;
    corrupt.back() ^= 1U;
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::ChecksumMismatch,
           "payload corruption is detected by checksum");

    corrupt = valid;
    writeU64(corrupt, 88,
             std::bit_cast<std::uint64_t>(
                 std::numeric_limits<double>::quiet_NaN()));
    refreshChecksum(corrupt);
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::InvalidData,
           "non-finite node state is rejected after checksum verification");

    softwing::SoftBodyCheckpointPersistenceLimits limits;
    limits.maximumEncodedBytes = valid.size() - 1;
    reject(valid,
           softwing::SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
           "configured byte limit is enforced", limits);

    limits = {};
    limits.maximumNodes = 9;
    reject(valid,
           softwing::SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
           "configured node limit is enforced", limits);

    // Payload header is 64 bytes, followed by fixed mutable records. The
    // multiplier count is therefore reached without trusting any wire count.
    constexpr std::size_t multiplierCountOffset =
        24 + 64 + 10 * 96 + 4 * 8 + 1 * 8 + 1 * 56 + 1 * 8;
    corrupt = valid;
    writeU64(corrupt, multiplierCountOffset, 5'000'001);
    refreshChecksum(corrupt);
    reject(corrupt,
           softwing::SoftBodyCheckpointPersistenceErrorCode::LimitExceeded,
           "impossible dynamic count is rejected before allocation");

    std::vector<std::uint8_t> output{1, 2, 3};
    softwing::SoftBodyCheckpointPersistenceError error;
    limits = {};
    limits.maximumNodes = 9;
    check(!softwing::serializeSoftBodyCheckpoint(
              body.checkpoint(), output, &error, limits)
              && error.code
                  == softwing::SoftBodyCheckpointPersistenceErrorCode::
                      LimitExceeded
              && output == std::vector<std::uint8_t>({1, 2, 3}),
          "failed serialization preserves caller output");

    const softwing::SoftBodyCheckpoint empty;
    check(!softwing::serializeSoftBodyCheckpoint(empty, output, &error)
              && error.code
                  == softwing::SoftBodyCheckpointPersistenceErrorCode::
                      InvalidData,
          "empty checkpoint is rejected");
}

} // namespace

int main() {
    try {
        testRoundTripAndContinuation();
        testTopologyBinding();
        testCorruptionAndLimits();
        if (failures != 0) {
            std::fprintf(
                stderr,
                "%d SoftBody checkpoint persistence check(s) failed\n",
                failures);
            return 1;
        }
        std::puts("SoftBody checkpoint persistence tests passed");
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "unexpected exception: %s\n", error.what());
        return 2;
    }
}
