#pragma once

#include "scene.h"
#include "structure.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace simwing::fsi {

// Assembly is deliberately bounded independently of Scene's in-memory
// representation. Lower limits are useful to applications accepting
// untrusted scene payloads and to tests of the transactional failure path.
struct SceneStructureLimits {
    std::size_t maximumNodes = 10'000'000;
    std::size_t maximumTriangles = 10'000'000;
    std::size_t maximumConstraints = 10'000'000;
    std::size_t maximumMappingBytes = 256ULL * 1024ULL * 1024ULL;
};

// Solver/contact policy is supplied by the worker rather than guessed from
// scene geometry. Scene-v2 carries authoritative fabric/line/pilot data, but
// it does not yet carry a verified contact material model.
struct SceneStructureSettings {
    // Paired seam vertices describe sewn coincidence. This tolerance accepts
    // only exporter roundoff; it is not a stitch-gap or compliance model.
    double seamCoincidenceToleranceMeters = 1.0e-9;
    int suspensionSolverIterations = 12;
    double suspensionAttachmentTolerance = 1.0e-10;
    double suspensionMinimumLineLengthMeters = 1.0e-10;
    double suspensionMaximumLineResidualMeters = 2.0e-4;
    double suspensionMaximumControlWorkJoules = 1.0e6;
    std::optional<StructureFabricContactDefinition> fabricSelfContact;

    auto operator<=>(const SceneStructureSettings&) const = default;
};

enum class SceneStructureDiagnosticSeverity : std::uint8_t {
    Error = 1,
};

enum class SceneStructureDiagnosticCode : std::uint16_t {
    InvalidScene = 1,
    MappingOverflow = 2,
    MaterialIncompatible = 3,
    DegenerateMaterialCoordinates = 4,
    ZeroMassDynamicNode = 5,
    UnsupportedPilot = 6,
    UnsupportedPilotHarness = 7,
    UnsupportedSuspensionTopology = 8,
    UnsupportedSeam = 9,
    UnsupportedBendingTopology = 10,
};

struct SceneStructureDiagnostic {
    SceneStructureDiagnosticSeverity severity =
        SceneStructureDiagnosticSeverity::Error;
    SceneStructureDiagnosticCode code =
        SceneStructureDiagnosticCode::InvalidScene;
    EntityKind entityKind = EntityKind::Scene;
    StableId entityId = invalidStableId;
    // Present when InvalidScene wraps a scene-v2 validation diagnostic.
    std::optional<ValidationCode> sceneValidationCode;
    std::string message;

    auto operator<=>(const SceneStructureDiagnostic&) const = default;
};

// The index-to-ID arrays follow StructureDefinition's deterministic ordering.
// Vertex nodes come first and junction nodes follow; each partition and each
// scalar-ID entity mapping is strictly increasing, so reverse lookups are
// deterministic binary searches rather than hash-table state. Dihedral pairs
// follow canonical undirected-edge order.
struct SceneStructureMappings {
    std::vector<StableId> nodeVertexIds;
    // Junction nodes follow all vertex nodes in StructureDefinition::nodes.
    std::vector<StableId> nodeSuspensionJunctionIds;
    std::vector<StableId> triangleIds;
    std::vector<StableId> membraneTriangleIds;
    std::vector<StableId> constraintSuspensionLineIds;
    struct SeamConstraintRange {
        StableId seamId = invalidStableId;
        std::size_t firstConstraint = 0;
        std::size_t constraintCount = 0;

        auto operator<=>(const SeamConstraintRange&) const = default;
    };
    // Sorted by seamId. Each range covers its pair stitches followed by the
    // two half-rigidity axial rails in StructureDefinition::constraints.
    std::vector<SeamConstraintRange> constraintSeamRanges;
    // Populated instead of constraintSuspensionLineIds when a rigid pilot
    // owns the directed suspension tree.
    std::vector<StableId> suspensionSegmentLineIds;
    // Aligned with StructureSuspensionDefinition::harnessPoints.
    std::vector<StableId> pilotHarnessAttachmentIds;
    // One stable triangle-ID pair per StructureDefinition::dihedrals entry.
    std::vector<std::array<StableId, 2>> dihedralTriangleIds;

    [[nodiscard]] std::optional<std::size_t> nodeIndex(
        StableId vertexId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> triangleIndex(
        StableId triangleId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> junctionNodeIndex(
        StableId junctionId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> membraneIndex(
        StableId triangleId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> constraintIndex(
        StableId suspensionLineId) const noexcept;
    [[nodiscard]] std::optional<SeamConstraintRange> seamConstraintRange(
        StableId seamId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> suspensionSegmentIndex(
        StableId suspensionLineId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> pilotHarnessIndex(
        StableId attachmentId) const noexcept;
};

struct SceneStructureAssembly {
    // assembled is false and all production data below is empty after any
    // error. Callers cannot accidentally run a partially converted scene.
    bool assembled = false;
    StructureDefinition definition;
    SceneStructureMappings mappings;
    SceneStructureSettings settings;
    double totalFabricMassKg = 0.0;
    double totalSeamMassKg = 0.0;
    std::vector<SceneStructureDiagnostic> diagnostics;

    [[nodiscard]] bool ok() const noexcept;
};

// Converts validated scene-v2 geometry into the currently checkpoint-safe
// Structure boundary. Vertices used by fabric triangles or surface
// attachments are dynamic; geometry used only to describe an opening loop is
// retained by scene-v2 but has no independent structural degree of freedom.
// Fabric mass is lumped from each triangle's undeformed material-chart area
// equally to its three vertices. In-plane warp, weft, shear, and damping properties are mapped to
// one Bulk membrane per triangle. A consistently oriented two-triangle
// manifold edge within the same authored sheet ID becomes one dihedral;
// merely welded edges between different sheets do not. If h0/h1 and L0/L1 are
// the altitudes and
// shared-edge lengths in each triangle's authoritative material chart, and
// D0/D1 are their bending rigidities in N*m, its XPBD angular compliance is
// `h0/(2*D0*L0) + h1/(2*D1*L1)` in radians/(N*m). This is the series
// compliance of two half-width fabric strips. Spatial geometry supplies only
// the signed rest angle. A zero rigidity produces no hinge.
//
// SurfaceVertex and SuspensionJunction attachments become explicit cable graph
// nodes; a dynamic junction must carry its scene-declared positive mass. A
// scene with one rigid pilot is oriented deterministically from surface leaves
// through junctions to harness roots and enters Structure's composite
// suspension solver. A scene without a pilot retains direct nodal cable
// constraints. Cable compliance is restLength/EA. Line mass and drag are
// intentionally not included in nodal mass or force ledgers by this bridge.
// A seam receives one rigid zero-rest stitch per authored vertex pair. Its
// assembled axial rigidity is split equally between matching distance-
// constraint rails on the two chains, and each centreline segment's mass is
// split equally across its four endpoint nodes. The parallel rails recover the
// authored EA while pair stitches preserve coincidence without welding.
[[nodiscard]] SceneStructureAssembly assembleSceneStructure(
    const Scene& scene,
    const SceneStructureLimits& limits = {},
    const SceneStructureSettings& settings = {});

} // namespace simwing::fsi
