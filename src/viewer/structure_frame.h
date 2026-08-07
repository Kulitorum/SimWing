#pragma once

#include "structure.h"
#include "viewer_protocol.h"

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace simwing::fsi {
struct Scene;
struct SceneStructureAssembly;
}

namespace simwing::viewer {

struct StructureFrameTriangleMapping {
    std::uint64_t stableId = 0;
    std::uint64_t negativeRegionId = 0;
    std::uint64_t positiveRegionId = 0;
};

struct StructureFrameLineMapping {
    std::uint64_t stableId = 0;
    std::uint32_t role = 0;
    // Direct Structure constraints derive these from topology. Suspension
    // segments use explicit indices because one endpoint may be a rigid
    // harness vertex appended after all SoftBody nodes.
    std::uint32_t vertex0 = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t vertex1 = std::numeric_limits<std::uint32_t>::max();
};

struct StructureFrameMappingDefinition {
    std::vector<std::uint64_t> vertexStableIds;
    std::vector<StructureFrameTriangleMapping> triangles;
    std::vector<StructureFrameLineMapping> lines;
};

// This is an immutable, one-to-one mapping from a Structure definition to the
// stable IDs and side-region IDs used by diagnostic frames. Triangle and line
// Triangle entries follow StructureDefinition::triangles. Line entries contain
// all direct constraints followed by rigid-payload suspension segments.
// Construction rejects incomplete, ambiguous, or foreign topology before a
// run can publish a misleading frame.
class StructureFrameMapping final {
public:
    StructureFrameMapping(
        const fsi::Structure& structure,
        StructureFrameMappingDefinition definition);

    StructureFrameMapping(const StructureFrameMapping&) = default;
    StructureFrameMapping(StructureFrameMapping&&) noexcept = default;
    StructureFrameMapping& operator=(const StructureFrameMapping&) = delete;
    StructureFrameMapping& operator=(StructureFrameMapping&&) = delete;

    [[nodiscard]] std::uint64_t structureFingerprint() const noexcept;
    [[nodiscard]] const std::vector<std::uint64_t>&
    vertexStableIds() const noexcept;
    [[nodiscard]] const std::vector<StructureFrameTriangleMapping>&
    triangles() const noexcept;
    [[nodiscard]] const std::vector<StructureFrameLineMapping>&
    lines() const noexcept;

private:
    const std::uint64_t structureFingerprint_;
    const std::vector<std::uint64_t> vertexStableIds_;
    const std::vector<StructureFrameTriangleMapping> triangles_;
    const std::vector<StructureFrameLineMapping> lines_;
};

// Derives viewer identity exclusively from a successful scene-v2 structural
// assembly. The supplied assembly is checked against a canonical assembly of
// the validated scene, and both are checked against the live Structure before
// any mapping is returned. Scene collection order therefore cannot change the
// result, while stale, failed, or hand-edited assemblies are rejected.
[[nodiscard]] StructureFrameMapping makeStructureFrameMapping(
    const fsi::Scene& scene,
    const fsi::SceneStructureAssembly& assembly,
    const fsi::Structure& structure);

struct StructureFrameContext {
    std::string sceneChecksum;
    std::string solverCommit;
    double timeStepSeconds = 0.0;
    std::uint32_t couplingIteration = 0;
    CouplingResiduals couplingResiduals;
    ConservationValues conservation;
};

// Builds an owning frame from a committed Structure state. The initial field
// contract is deliberately limited to quantities exposed by Structure:
//
// - per-vertex velocity and pending external force;
// - per-line geometric length and unilateral/bilateral violation;
// - aggregate applied/pending force, mass, momentum, energy, and solver-error
//   diagnostics.
//
// Structure currently exposes membrane strain/residual only as global maxima,
// so no per-triangle strain field is emitted. It exposes the last applied load
// only as an aggregate, so no per-vertex applied-force field is invented.
// Direct-constraint multipliers remain private to softwing_core. Suspension
// tension is available, but a mixed field would invent values for direct
// constraints, so this adapter does not publish a line-tension field yet.
[[nodiscard]] DiagnosticFrame buildStructureFrame(
    const fsi::Structure& structure,
    const StructureFrameMapping& mapping,
    const StructureFrameContext& context);

} // namespace simwing::viewer
