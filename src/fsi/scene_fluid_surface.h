#pragma once

#include "scene_structure.h"

#include <array>
#include <cstddef>
#include <compare>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidSurfaceDefinitionVersion = 1;
inline constexpr std::uint32_t sceneFluidSurfaceStateVersion = 1;

struct SceneFluidSurfaceLimits {
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumMaterials = 1'000'000;
    std::size_t maximumVertices = 10'000'000;
    std::size_t maximumTriangles = 10'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumOpeningVertices = 10'000'000;
    std::size_t maximumMappingBytes = 256ULL * 1024ULL * 1024ULL;
};

enum class SceneFluidSurfaceDiagnosticCode : std::uint16_t {
    InvalidScene = 1,
    LimitExceeded = 2,
    OpeningVertexWithoutFabricMotion = 3,
    MappingOverflow = 4,
};

struct SceneFluidSurfaceDiagnostic {
    SceneFluidSurfaceDiagnosticCode code =
        SceneFluidSurfaceDiagnosticCode::InvalidScene;
    EntityKind entityKind = EntityKind::Scene;
    StableId entityId = invalidStableId;
    std::optional<ValidationCode> sceneValidationCode;
    std::string message;

    auto operator<=>(const SceneFluidSurfaceDiagnostic&) const = default;
};

struct SceneFluidSurfaceRegion {
    StableId id = invalidStableId;
    RegionKind kind = RegionKind::Cell;

    bool operator==(const SceneFluidSurfaceRegion&) const = default;
};

struct SceneFluidSurfaceMaterial {
    StableId id = invalidStableId;
    double porosityFraction = 0.0;
    double permeabilitySquareMeters = 0.0;

    bool operator==(const SceneFluidSurfaceMaterial&) const = default;
};

struct SceneFluidSurfaceVertex {
    StableId id = invalidStableId;
    Vec3 referencePositionMeters;

    bool operator==(const SceneFluidSurfaceVertex& other) const {
        return id == other.id
               && referencePositionMeters.x == other.referencePositionMeters.x
               && referencePositionMeters.y == other.referencePositionMeters.y
               && referencePositionMeters.z == other.referencePositionMeters.z;
    }
};

struct SceneFluidSurfaceTriangle {
    StableId id = invalidStableId;
    std::array<std::size_t, 3> vertexIndices{};
    std::size_t negativeSideRegionIndex = 0;
    std::size_t positiveSideRegionIndex = 0;
    std::size_t materialIndex = 0;
    StableId sheetId = invalidStableId;
    SurfaceRole role = SurfaceRole::Skin;

    bool operator==(const SceneFluidSurfaceTriangle&) const = default;
};

struct SceneFluidSurfaceOpening {
    StableId id = invalidStableId;
    std::vector<std::size_t> orderedVertexIndices;
    std::size_t negativeSideRegionIndex = 0;
    std::size_t positiveSideRegionIndex = 0;
    OpeningRole role = OpeningRole::Intake;

    bool operator==(const SceneFluidSurfaceOpening&) const = default;
};

struct SceneFluidSurfaceMappings {
    std::vector<StableId> regionIds;
    std::vector<StableId> materialIds;
    std::vector<StableId> vertexIds;
    std::vector<StableId> triangleIds;
    std::vector<StableId> openingIds;

    [[nodiscard]] std::optional<std::size_t> regionIndex(
        StableId id) const noexcept;
    [[nodiscard]] std::optional<std::size_t> materialIndex(
        StableId id) const noexcept;
    [[nodiscard]] std::optional<std::size_t> vertexIndex(
        StableId id) const noexcept;
    [[nodiscard]] std::optional<std::size_t> triangleIndex(
        StableId id) const noexcept;
    [[nodiscard]] std::optional<std::size_t> openingIndex(
        StableId id) const noexcept;

    bool operator==(const SceneFluidSurfaceMappings&) const = default;
};

// Canonical, compact surface ownership extracted from scene-v2. Triangle
// winding and the authored negative/positive region ordering are preserved;
// no grid crossing, cut-cell classification, or pressure model is invented.
struct SceneFluidSurfaceDefinition {
    std::uint32_t version = sceneFluidSurfaceDefinitionVersion;
    std::uint64_t fingerprint = 0;
    std::vector<SceneFluidSurfaceRegion> regions;
    std::vector<SceneFluidSurfaceMaterial> materials;
    std::vector<SceneFluidSurfaceVertex> vertices;
    std::vector<SceneFluidSurfaceTriangle> triangles;
    std::vector<SceneFluidSurfaceOpening> openings;
    SceneFluidSurfaceMappings mappings;

    bool operator==(const SceneFluidSurfaceDefinition&) const = default;
};

struct SceneFluidSurfaceAssembly {
    bool assembled = false;
    SceneFluidSurfaceDefinition definition;
    std::vector<SceneFluidSurfaceDiagnostic> diagnostics;

    [[nodiscard]] bool ok() const noexcept;
};

// Opening loops must use fabric-triangle vertices. Scene-v2 deliberately
// permits opening-only construction vertices, but their motion has no trusted
// Structure degree of freedom or interpolation rule yet, so this adapter
// rejects them instead of silently freezing fluid topology.
[[nodiscard]] SceneFluidSurfaceAssembly assembleSceneFluidSurface(
    const Scene& scene,
    const SceneFluidSurfaceLimits& limits = {});

struct SceneFluidSurfaceVertexState {
    StableId id = invalidStableId;
    Vec3 positionMeters;
    Vec3 velocityMetersPerSecond;

    bool operator==(const SceneFluidSurfaceVertexState& other) const {
        return id == other.id
               && positionMeters.x == other.positionMeters.x
               && positionMeters.y == other.positionMeters.y
               && positionMeters.z == other.positionMeters.z
               && velocityMetersPerSecond.x == other.velocityMetersPerSecond.x
               && velocityMetersPerSecond.y == other.velocityMetersPerSecond.y
               && velocityMetersPerSecond.z == other.velocityMetersPerSecond.z;
    }
};

struct SceneFluidSurfaceState {
    std::uint32_t version = sceneFluidSurfaceStateVersion;
    std::uint64_t definitionFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<SceneFluidSurfaceVertexState> vertices;

    bool operator==(const SceneFluidSurfaceState&) const = default;
};

// Captures one immutable accepted Structure epoch in the definition's stable
// vertex order. The scene-to-Structure mapping is checked in full before any
// state is returned.
[[nodiscard]] SceneFluidSurfaceState captureSceneFluidSurfaceState(
    const SceneFluidSurfaceDefinition& definition,
    const SceneStructureMappings& structureMappings,
    const Structure& structure);

} // namespace simwing::fsi
