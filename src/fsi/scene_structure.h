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

// The index-to-ID arrays are in the exact order of StructureDefinition's
// arrays. They are also strictly increasing, so the reverse lookups are
// deterministic binary searches rather than hash-table state.
struct SceneStructureMappings {
    std::vector<StableId> nodeVertexIds;
    std::vector<StableId> triangleIds;
    std::vector<StableId> membraneTriangleIds;
    std::vector<StableId> constraintSuspensionLineIds;

    [[nodiscard]] std::optional<std::size_t> nodeIndex(
        StableId vertexId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> triangleIndex(
        StableId triangleId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> membraneIndex(
        StableId triangleId) const noexcept;
    [[nodiscard]] std::optional<std::size_t> constraintIndex(
        StableId suspensionLineId) const noexcept;
};

struct SceneStructureAssembly {
    // assembled is false and all production data below is empty after any
    // error. Callers cannot accidentally run a partially converted scene.
    bool assembled = false;
    StructureDefinition definition;
    SceneStructureMappings mappings;
    double totalFabricMassKg = 0.0;
    std::vector<SceneStructureDiagnostic> diagnostics;

    [[nodiscard]] bool ok() const noexcept;
};

// Converts validated scene-v2 geometry into the currently checkpoint-safe
// Structure boundary. All scene vertices are dynamic. Fabric mass is lumped
// from each triangle's undeformed material-chart area equally to its three
// vertices. In-plane warp, weft, shear, and damping properties are mapped to
// one Bulk membrane per triangle; scene bending is not assembled in this
// slice.
//
// A SuspensionLine is currently representable only when both attachments are
// distinct SurfaceVertex targets. It becomes a unilateral cable with
// compliance restLength/EA. Line mass and drag are intentionally not included
// in nodal mass or force ledgers by this bridge. Pilot/harness topology is a
// hard unsupported diagnostic until the rigid-payload checkpoint API exists.
[[nodiscard]] SceneStructureAssembly assembleSceneStructure(
    const Scene& scene,
    const SceneStructureLimits& limits = {});

} // namespace simwing::fsi
