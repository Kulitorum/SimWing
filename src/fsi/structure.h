#pragma once

#include <array>
#include <cstddef>
#include <compare>
#include <cstdint>
#include <memory>
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

struct StructureDefinition {
    std::vector<StructureNodeDefinition> nodes;
    std::vector<StructureTriangleDefinition> triangles;
    std::vector<StructureConstraintDefinition> constraints;
    std::vector<StructureMembraneDefinition> membranes;
    std::vector<StructureDihedralDefinition> dihedrals;
};

struct StructureNodeState {
    StructureVector3 positionMeters;
    StructureVector3 previousPositionMeters;
    StructureVector3 velocityMetersPerSecond;

    auto operator<=>(const StructureNodeState&) const = default;
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
    double totalDynamicMassKg = 0.0;
    StructureVector3 centerOfMassMeters;
    StructureVector3 linearMomentumKgMetersPerSecond;
    double kineticEnergyJoules = 0.0;
    double maximumDistanceErrorMeters = 0.0;
    double maximumCableExtensionMeters = 0.0;
    double maximumAbsoluteMembraneStrain = 0.0;
    double maximumMembraneResidual = 0.0;
    StructureVector3 pendingExternalForceNewtons;
    StructureVector3 lastAppliedExternalForceNewtons;
    bool finite = true;

    auto operator<=>(const StructureDiagnostics&) const = default;
};

inline constexpr std::uint32_t structureCheckpointVersion = 1;

// A checkpoint is a committed macro-step state for the primitives this
// adapter currently wraps. Rebuilding on restore deliberately clears XPBD
// iteration multipliers, which softwing also clears at the beginning of every
// substep. Pending nodal loads are retained.
//
// Registered contact and the rigid-payload SuspensionSystem are intentionally
// absent: their persistent warm-start/private state has no production
// checkpoint API in softwing_core. They must not be added here until that API
// exists; using test friends would make rollback incomplete and unsafe.
struct StructureCheckpoint {
    std::uint32_t version = structureCheckpointVersion;
    std::uint64_t definitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<StructureNodeState> nodes;
    std::vector<StructureVector3> pendingExternalForcesNewtons;
    StructureVector3 lastAppliedExternalForceNewtons;
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
