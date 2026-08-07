#pragma once

#include <array>
#include <cstddef>
#include <compare>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace simwing::fsi {

// Structure owns SI-valued XPBD state without exposing softwing types. Scene-v2
// conversion deliberately remains a separate boundary so this adapter and its
// analytic tests do not depend on a serialization format.
struct StructureVector2 {
    double x = 0.0;
    double y = 0.0;
};

struct StructureVector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    auto operator<=>(const StructureVector3&) const = default;
};

struct StructureNodeDefinition {
    StructureVector3 positionMeters;
    double massKg = 0.0;
    bool fixed = false;
};

struct StructureTriangleDefinition {
    std::array<std::size_t, 3> nodes{};
};

enum class StructureConstraintKind : std::uint8_t {
    Distance = 1,
    Cable = 2,
    SuspensionTie = 3,
};

struct StructureConstraintDefinition {
    StructureConstraintKind kind = StructureConstraintKind::Distance;
    std::size_t firstNode = 0;
    std::size_t secondNode = 0;
    double restLengthMeters = 0.0;
    double complianceMetersPerNewton = 0.0;
};

enum class StructureMaterialRole : std::uint8_t {
    Bulk = 1,
    Seam = 2,
    Reinforcement = 3,
};

struct StructureMembraneMaterial {
    double warpStiffnessNewtonsPerMeter = 0.0;
    double weftStiffnessNewtonsPerMeter = 0.0;
    double couplingStiffnessNewtonsPerMeter = 0.0;
    double shearStiffnessNewtonsPerMeter = 0.0;
    double warpPreTensionNewtonsPerMeter = 0.0;
    double weftPreTensionNewtonsPerMeter = 0.0;
    double dampingSeconds = 0.0;
    double compressionStiffnessRatio = 1.0;
};

struct StructureMembraneDefinition {
    std::size_t triangle = 0;
    std::array<StructureVector2, 3> materialCoordinates{};
    StructureMembraneMaterial material;
    StructureMaterialRole role = StructureMaterialRole::Bulk;
};

struct StructureDihedralDefinition {
    std::array<std::size_t, 4> nodes{};
    double restAngleRadians = 0.0;
    double complianceRadiansPerNewtonMeter = 0.0;
};

// Optional exhaustive self-contact for the contiguous fabric triangle range.
// There is deliberately no implicit material default: scene assembly enables
// it only when the caller supplies explicit thickness/compliance/friction.
struct StructureFabricContactDefinition {
    double halfThicknessMeters = 0.0;
    double normalComplianceMetersPerNewton = 0.0;
    double staticFriction = 0.0;
    double dynamicFriction = 0.0;

    auto operator<=>(const StructureFabricContactDefinition&) const = default;
};

struct StructureQuaternion {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    auto operator<=>(const StructureQuaternion&) const = default;
};

enum class StructureSuspensionEndpointKind : std::uint8_t {
    SurfaceAttachment = 1,
    Junction = 2,
    PilotHarness = 3,
};

struct StructureSuspensionEndpointDefinition {
    StructureSuspensionEndpointKind kind =
        StructureSuspensionEndpointKind::SurfaceAttachment;
    std::uint64_t stableId = 0;
};

struct StructureSuspensionAttachmentDefinition {
    std::uint64_t stableId = 0;
    std::size_t node = 0;
};

struct StructureSuspensionJunctionDefinition {
    std::uint64_t stableId = 0;
    std::size_t node = 0;
};

struct StructurePilotHarnessDefinition {
    std::uint64_t stableId = 0;
    StructureVector3 localPositionMeters;
};

struct StructureSuspensionSegmentDefinition {
    std::uint64_t stableId = 0;
    StructureSuspensionEndpointDefinition from;
    StructureSuspensionEndpointDefinition to;
    double restLengthMeters = 0.0;
    double axialStiffnessNewtons = 0.0;
    double axialDampingNewtonSecondsPerMeter = 0.0;
    std::uint32_t role = 0;
};

struct StructureSuspensionDefinition {
    std::uint64_t pilotStableId = 0;
    double pilotMassKg = 0.0;
    StructureVector3 pilotCenterOfMassLocalMeters;
    StructureVector3 pilotInitialCenterOfMassWorldMeters;
    StructureVector3 pilotInitialLinearVelocityMetersPerSecond;
    StructureVector3 pilotInitialAngularVelocityRadiansPerSecond;
    StructureQuaternion pilotInitialBodyToWorld;
    StructureVector3 pilotPrincipalInertiaKgSquareMeters;
    std::vector<StructureSuspensionAttachmentDefinition> attachments;
    std::vector<StructureSuspensionJunctionDefinition> junctions;
    std::vector<StructurePilotHarnessDefinition> harnessPoints;
    std::vector<StructureSuspensionSegmentDefinition> segments;
    int solverIterations = 12;
    double attachmentTolerance = 1.0e-10;
    double minimumLineLengthMeters = 1.0e-10;
    double maximumLineResidualMeters = 2.0e-4;
    double maximumControlWorkJoules = 1.0e6;
};

struct StructureDefinition {
    std::vector<StructureNodeDefinition> nodes;
    std::vector<StructureTriangleDefinition> triangles;
    std::vector<StructureConstraintDefinition> constraints;
    std::vector<StructureMembraneDefinition> membranes;
    std::vector<StructureDihedralDefinition> dihedrals;
    std::optional<StructureFabricContactDefinition> fabricSelfContact;
    std::optional<StructureSuspensionDefinition> suspension;
};

struct StructureNodeState {
    StructureVector3 positionMeters;
    StructureVector3 previousPositionMeters;
    StructureVector3 velocityMetersPerSecond;

    auto operator<=>(const StructureNodeState&) const = default;
};

struct StructureRigidPayloadState {
    StructureVector3 centerOfMassWorldMeters;
    StructureQuaternion bodyToWorld;
    StructureVector3 linearVelocityMetersPerSecond;
    StructureVector3 angularVelocityRadiansPerSecond;

    auto operator<=>(const StructureRigidPayloadState&) const = default;
};

struct StructureSuspensionSegmentState {
    std::uint64_t stableId = 0;
    double currentLengthMeters = 0.0;
    double commandedRestLengthMeters = 0.0;
    double tensionNewtons = 0.0;

    auto operator<=>(const StructureSuspensionSegmentState&) const = default;
};

struct StructureSuspensionState {
    StructureRigidPayloadState payload;
    std::vector<StructureVector3> harnessPositionsMeters;
    std::vector<StructureVector3> harnessVelocitiesMetersPerSecond;
    std::vector<StructureSuspensionSegmentState> segments;

    auto operator<=>(const StructureSuspensionState&) const = default;
};

struct StructureStepSettings {
    double timeStepSeconds = 1.0 / 120.0;
    int substeps = 2;
    int constraintIterations = 12;
    int cableConstraintSweepPairs = 0;
    StructureVector3 gravityMetersPerSecondSquared{0.0, 0.0, -9.80665};
    double velocityDampingPerSecond = 0.25;
    StructureVector3 dampingReferenceVelocityMetersPerSecond;
    // Zero is the deterministic serial baseline. A nonzero value explicitly
    // selects softwing's reproducible parallel sweep.
    unsigned workerThreads = 0;
};

struct StructureDiagnostics {
    std::size_t nodeCount = 0;
    std::size_t dynamicNodeCount = 0;
    std::size_t triangleCount = 0;
    std::size_t constraintCount = 0;
    std::size_t membraneCount = 0;
    std::size_t dihedralCount = 0;
    std::size_t contactPairCount = 0;
    std::size_t activeContactCount = 0;
    std::size_t suspensionSegmentCount = 0;
    double totalDynamicMassKg = 0.0;
    StructureVector3 centerOfMassMeters;
    StructureVector3 linearMomentumKgMetersPerSecond;
    double kineticEnergyJoules = 0.0;
    double maximumDistanceErrorMeters = 0.0;
    double maximumCableExtensionMeters = 0.0;
    double maximumAbsoluteMembraneStrain = 0.0;
    double maximumMembraneResidual = 0.0;
    double maximumContactPenetrationMeters = 0.0;
    double maximumSuspensionResidualMeters = 0.0;
    StructureVector3 pendingExternalForceNewtons;
    StructureVector3 lastAppliedExternalForceNewtons;
    bool finite = true;

    auto operator<=>(const StructureDiagnostics&) const = default;
};

inline constexpr std::uint32_t structureCheckpointVersion = 2;

// A checkpoint is a committed macro-step state for the complete adapter. Its
// immutable private payload owns the soft-body/contact warm starts and the
// optional suspension/rigid-payload checkpoint. Pending nodal loads and public
// accounting remain visible here for orchestration and persistence adapters.
struct StructureCheckpoint {
    std::uint32_t version = structureCheckpointVersion;
    std::uint64_t definitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<StructureNodeState> nodes;
    std::vector<StructureVector3> pendingExternalForcesNewtons;
    StructureVector3 lastAppliedExternalForceNewtons;

private:
    friend class Structure;
    struct Detail;
    std::shared_ptr<const Detail> detail;
};

class Structure {
public:
    explicit Structure(StructureDefinition definition);
    ~Structure();

    Structure(Structure&&) noexcept;
    Structure& operator=(Structure&&) noexcept;
    Structure(const Structure&) = delete;
    Structure& operator=(const Structure&) = delete;

    [[nodiscard]] const StructureDefinition& definition() const noexcept;
    [[nodiscard]] std::uint64_t definitionFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] std::vector<StructureNodeState> nodeStates() const;
    [[nodiscard]] std::optional<StructureSuspensionState>
    suspensionState() const;

    void clearExternalForces() noexcept;
    void addExternalForce(std::size_t node,
                          const StructureVector3& forceNewtons);
    void setExternalForces(
        std::span<const StructureVector3> forcesNewtons);

    // A failed step restores the exact committed adapter checkpoint before
    // rethrowing. Loads are consumed only by an accepted step.
    [[nodiscard]] StructureDiagnostics step(
        const StructureStepSettings& settings);
    [[nodiscard]] StructureDiagnostics diagnostics() const;

    [[nodiscard]] StructureCheckpoint checkpoint() const;
    void restore(const StructureCheckpoint& checkpoint);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace simwing::fsi
