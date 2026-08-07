#include "softwing/canopy.h"
#include "softwing/suspension.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

bool sameDouble(double first, double second) {
    return std::bit_cast<std::uint64_t>(first)
           == std::bit_cast<std::uint64_t>(second);
}

bool sameVec(const softwing::Vec3& first, const softwing::Vec3& second) {
    return sameDouble(first.x, second.x) && sameDouble(first.y, second.y)
           && sameDouble(first.z, second.z);
}

bool sameNode(const softwing::Node& first, const softwing::Node& second) {
    return sameVec(first.position, second.position)
           && sameVec(first.previousPosition, second.previousPosition)
           && sameVec(first.velocity, second.velocity)
           && sameVec(first.force, second.force)
           && sameDouble(first.inverseMass, second.inverseMass);
}

softwing::SuspensionDefinition makeDefinition(
    const softwing::CanopyMesh& canopy,
    double lineLengthOffset = 0.0) {
    const auto face = std::find_if(
        canopy.faces.begin(), canopy.faces.end(),
        [](const softwing::CanopyFaceRecord& value) {
            return !value.virtualClosure;
        });
    if (face == canopy.faces.end()) {
        throw std::runtime_error("checkpoint fixture has no material face");
    }
    const softwing::Triangle& triangle =
        canopy.body.triangles()[face->triangle];
    const softwing::Vec3 attachmentPosition =
        canopy.body.nodes()[triangle.a].position;

    softwing::SuspensionDefinition definition;
    definition.identifier = "checkpoint-fixture";
    definition.description = "minimal suspension checkpoint fixture";
    definition.unitsFrameTag = std::string(softwing::suspensionStage5FrameTag);
    definition.provenance = {{"fixture", "synthetic checkpoint test"}};
    definition.attachments = {
        {"attachment", face->panelId, face->chart[0],
         softwing::SuspensionSide::Centre, "", "fixture"}};

    definition.payload.mass = 80.0;
    definition.payload.centreOfMassLocal = {};
    definition.payload.inertiaBody = {{8.0, 0.0, 0.0,
                                       0.0, 10.0, 0.0,
                                       0.0, 0.0, 6.0}};
    definition.payload.initialState.centreOfMassWorld =
        attachmentPosition + softwing::Vec3{0.0, 0.0, -1.0};
    definition.payload.hangPoints = {
        {"hang", {0.0, 0.0, 0.0}, softwing::SuspensionSide::Centre,
         "fixture"}};
    definition.payload.supportPoints = {
        {"support", {0.0, 0.0, -0.25}, softwing::SuspensionSide::Centre,
         "fixture"}};
    definition.payload.provenanceId = "fixture";

    definition.segments = {
        {"main-line",
         {softwing::SuspensionEndpointKind::Attachment, "attachment"},
         {softwing::SuspensionEndpointKind::HangPoint, "hang"},
         0.82 + lineLengthOffset, 2500.0, 1.5,
         softwing::SuspensionSide::Centre, {"main"}, "fixture"}};

    softwing::SuspensionControlDefinition brake;
    brake.id = "brake";
    brake.kind = softwing::SuspensionControlKind::Brake;
    brake.side = softwing::SuspensionSide::Centre;
    brake.minimumCommand = 0.0;
    brake.maximumCommand = 1.0;
    brake.neutralCommand = 0.0;
    brake.maximumRate = 4.0;
    brake.targets = {{softwing::SuspensionControlTargetKind::SegmentRestLength,
                      "main-line", -0.08, {}}};
    brake.provenanceId = "fixture";

    softwing::SuspensionControlDefinition weightShift;
    weightShift.id = "weight-shift";
    weightShift.kind = softwing::SuspensionControlKind::WeightShift;
    weightShift.side = softwing::SuspensionSide::Centre;
    weightShift.minimumCommand = -1.0;
    weightShift.maximumCommand = 1.0;
    weightShift.neutralCommand = 0.0;
    weightShift.maximumRate = 3.0;
    weightShift.targets = {
        {softwing::SuspensionControlTargetKind::HangPointTravel,
         "hang", 0.06, {1.0, 0.0, 0.0}}};
    weightShift.provenanceId = "fixture";
    definition.controls = {brake, weightShift};

    definition.solver.lineIterations = 16;
    definition.solver.maximumLineResidual = 10.0;
    definition.solver.maximumControlWork = 1.0e9;
    definition.ground.mode = softwing::PayloadGroundMode::SupportPlane;
    definition.ground.planeNormal = {0.0, 0.0, 1.0};
    definition.ground.planeOffset = attachmentPosition.z - 0.10;
    definition.ground.compliance = 0.0;
    definition.ground.penetrationFraction = 1.0;
    return definition;
}

struct Rig {
    softwing::CanopyMesh canopy;
    softwing::SuspensionDefinition definition;
    softwing::SuspensionSystem suspension;

    explicit Rig(double lineLengthOffset = 0.0)
        : canopy(softwing::buildCanopy(
              softwing::makeStraightThreeCellDefinition())),
          definition(makeDefinition(canopy, lineLengthOffset)),
          suspension(softwing::SuspensionSystem::build(canopy, definition)) {}

    Rig(const Rig&) = delete;
    Rig& operator=(const Rig&) = delete;
};

softwing::StepSettings stepSettings(const Rig& rig) {
    softwing::StepSettings settings;
    settings.timeStep = 1.0 / 240.0;
    settings.substeps = 1;
    settings.constraintIterations = rig.definition.solver.lineIterations;
    settings.gravity = {0.0, 0.0, -9.80665};
    settings.velocityDampingPerSecond = 0.0;
    return settings;
}

void advance(Rig& rig, const softwing::StepSettings& settings, int count) {
    for (int step = 0; step < count; ++step) rig.suspension.step(settings);
}

void expectRejectedTransactionally(
    Rig& rig,
    const softwing::SuspensionCheckpoint& rejected,
    const char* message) {
    const std::uint64_t before = rig.suspension.checkpoint().stateFingerprint;
    bool threw = false;
    try {
        rig.suspension.restore(rejected);
    } catch (const softwing::SuspensionError& error) {
        threw = error.phase() == softwing::SuspensionPhase::Validation;
    }
    check(threw, message);
    check(rig.suspension.checkpoint().stateFingerprint == before,
          "rejected restore preserves every suspension state field");
}

void testReplayAndCompleteState() {
    Rig rig;
    const softwing::StepSettings settings = stepSettings(rig);

    softwing::RigidPayloadState state = rig.suspension.payloadState();
    state.centreOfMassWorld += {0.025, -0.015, 0.01};
    state.orientation = softwing::normalizedCanonical(
        softwing::rotationIncrement({0.08, -0.04, 0.03}));
    state.linearVelocity = {0.3, -0.12, 0.08};
    state.angularVelocity = {0.07, -0.05, 0.04};
    rig.suspension.setPayloadState(state);
    rig.suspension.setAppliedPayloadWrench(
        {{4.0, -1.0, 2.0}, {0.25, -0.15, 0.1}});
    rig.suspension.setControlTarget("brake", 0.8);
    rig.suspension.setControlTarget("weight-shift", -0.55);
    advance(rig, settings, 3);
    state = rig.suspension.payloadState();
    state.centreOfMassWorld.z = rig.definition.ground.planeOffset - 1.0;
    state.linearVelocity = {0.3, -0.12, -0.08};
    state.angularVelocity = {0.07, -0.05, 0.04};
    rig.suspension.setPayloadState(state);
    advance(rig, settings, 1);

    const softwing::SuspensionCheckpoint saved = rig.suspension.checkpoint();
    const std::vector<softwing::Node> savedNodes = rig.canopy.body.nodes();
    check(saved.stateFingerprint != 0 && saved.topologyFingerprint != 0,
          "checkpoint carries structure and full-state fingerprints");
    check(saved.payloadState.orientation.x != 0.0
              && saved.payloadState.linearVelocity.z != 0.0
              && saved.payloadState.angularVelocity.x != 0.0,
          "checkpoint covers rigid pose and velocities");
    check(saved.controls[0].actualCommand != 0.0
              || saved.controls[1].actualCommand != 0.0,
          "checkpoint covers rate-limited controls");
    check(!saved.groundMultipliers.empty()
              && saved.groundMultipliers.front() > 0.0,
          "checkpoint covers support-plane multiplier state");
    check(saved.appliedWrench.force.x == 4.0
              && saved.currentGravity.z == settings.gravity.z,
          "checkpoint covers persistent wrench and gravity");
    check(saved.diagnostics.registered
              && saved.committedDiagnostics.registered
              && saved.segmentDiagnostics.size() == saved.segments.size(),
          "checkpoint covers current and committed diagnostics");

    advance(rig, settings, 3);
    const softwing::SuspensionCheckpoint expected =
        rig.suspension.checkpoint();
    const std::vector<softwing::Node> expectedNodes = rig.canopy.body.nodes();

    rig.suspension.setControlTarget("brake", 0.0);
    rig.suspension.setControlTarget("weight-shift", 0.9);
    rig.suspension.setAppliedPayloadWrench(
        {{-8.0, 3.0, -4.0}, {-0.5, 0.2, 0.3}});
    advance(rig, settings, 5);

    rig.canopy.body.nodes() = savedNodes;
    rig.suspension.restore(saved);
    const softwing::SuspensionCheckpoint restored =
        rig.suspension.checkpoint();
    check(restored.stateFingerprint == saved.stateFingerprint,
          "restore reproduces all captured suspension state and diagnostics");

    advance(rig, settings, 3);
    const softwing::SuspensionCheckpoint replayed =
        rig.suspension.checkpoint();
    check(replayed.stateFingerprint == expected.stateFingerprint,
          "checkpoint restore replay is bit-identical");
    check(rig.canopy.body.nodes().size() == expectedNodes.size(),
          "replay retains structural node count");
    for (std::size_t index = 0;
         index < rig.canopy.body.nodes().size() &&
         index < expectedNodes.size(); ++index) {
        check(sameNode(rig.canopy.body.nodes()[index], expectedNodes[index]),
              "checkpoint replay body coupling is bit-identical");
    }
}

void testRejectionAndIntegrity() {
    Rig rig;
    const softwing::StepSettings settings = stepSettings(rig);
    rig.suspension.setControlTarget("brake", 0.6);
    advance(rig, settings, 2);
    const softwing::SuspensionCheckpoint saved = rig.suspension.checkpoint();

    softwing::SuspensionCheckpoint corrupt = saved;
    corrupt.stateFingerprint ^= 1U;
    expectRejectedTransactionally(rig, corrupt,
                                  "corrupt state fingerprint is rejected");

    corrupt = saved;
    corrupt.payloadState.linearVelocity.x =
        std::numeric_limits<double>::quiet_NaN();
    expectRejectedTransactionally(rig, corrupt,
                                  "non-finite rigid checkpoint is rejected");

    corrupt = saved;
    corrupt.segments.pop_back();
    expectRejectedTransactionally(rig, corrupt,
                                  "size-mismatched checkpoint is rejected");

    corrupt = saved;
    corrupt.controls.front().id = "foreign-control";
    expectRejectedTransactionally(rig, corrupt,
                                  "corrupt control identity is rejected");

    Rig equivalent;
    equivalent.suspension.restore(saved);
    check(equivalent.suspension.checkpoint().stateFingerprint
              == saved.stateFingerprint,
          "equivalent rebuilt topology accepts checkpoint without restoring "
          "owner pointers or lifetime tokens");

    Rig foreign(0.025);
    expectRejectedTransactionally(
        foreign, saved,
        "foreign normalized suspension definition is rejected");

    const double originalInverseMass =
        rig.canopy.body.nodes().front().inverseMass;
    rig.canopy.body.nodes().front().inverseMass *= 0.5;
    bool topologyRejected = false;
    try {
        rig.suspension.restore(saved);
    } catch (const softwing::SuspensionError& error) {
        topologyRejected = error.phase() == softwing::SuspensionPhase::Validation;
    }
    check(topologyRejected, "changed live structural topology is rejected");
    rig.canopy.body.nodes().front().inverseMass = originalInverseMass;
    check(rig.suspension.checkpoint().stateFingerprint
              == saved.stateFingerprint,
          "topology rejection preserves suspension state");

    rig.suspension.restore(saved);
    advance(rig, settings, 1);
    check(rig.suspension.diagnostics().registered,
          "restore leaves a no-active-substep solver-safe point");
}

} // namespace

int main() {
    testReplayAndCompleteState();
    testRejectionAndIntegrity();
    if (failures == 0) {
        std::printf("softwing suspension checkpoint tests passed\n");
    }
    return failures == 0 ? 0 : 1;
}
