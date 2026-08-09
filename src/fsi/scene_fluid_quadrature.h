#pragma once

#include "fluid/scene_surface_ownership.h"
#include "scene_fluid_surface_transfer.h"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidQuadratureVersion = 2;

enum class SceneFluidQuadratureOwnerKind : std::uint8_t {
    Cell = 1,
    Face = 2,
};

struct SceneFluidQuadraturePoint {
    std::uint64_t stableId = 0;
    StableId triangleId = invalidStableId;
    std::array<double, 3> barycentricCoordinates{};
    double areaSquareMeters = 0.0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    // Cell-owned patches sample both sides in their unique owner cell.
    // Face-owned patches retain the exact adjacent cell selected by authored
    // winding for each side, so downstream pressure lookup is unambiguous.
    std::size_t negativeSideCellIndex = 0;
    std::size_t positiveSideCellIndex = 0;
    StableId materialId = invalidStableId;
    StableId sheetId = invalidStableId;
    SurfaceRole role = SurfaceRole::Skin;
    SceneFluidQuadratureOwnerKind ownerKind =
        SceneFluidQuadratureOwnerKind::Cell;

    bool operator==(const SceneFluidQuadraturePoint&) const = default;
};

struct SceneFluidQuadratureDefinition {
    std::uint32_t version = sceneFluidQuadratureVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t patchOwnershipFingerprint = 0;
    std::uint64_t couplingSurfaceFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::vector<SceneFluidQuadraturePoint> points;

    bool operator==(const SceneFluidQuadratureDefinition&) const = default;
};

struct SceneFluidQuadratureTraction {
    std::uint64_t stableId = 0;
    StructureVector3 tractionPascals;

    bool operator==(const SceneFluidQuadratureTraction&) const = default;
};

struct SceneFluidQuadratureKinematics {
    std::uint64_t stableId = 0;
    StructureVector3 positionMeters;
    StructureVector3 velocityMetersPerSecond;

    bool operator==(const SceneFluidQuadratureKinematics&) const = default;
};

[[nodiscard]] SceneFluidQuadratureDefinition buildSceneFluidQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::SceneFluidGridCandidateSet& candidates,
    const fluid::SceneFluidGridIntersectionSet& intersections,
    const fluid::SceneFluidGridPatchSet& patches,
    const fluid::SceneFluidPatchOwnership& ownership,
    const SceneFluidSurfaceTransfer& transfer);

void validateSceneFluidQuadratureDefinition(
    const SceneFluidQuadratureDefinition& definition);

[[nodiscard]] std::vector<SceneFluidQuadratureKinematics>
sampleSceneFluidQuadratureKinematics(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& definition);

[[nodiscard]] ConservativeTransferResult evaluateSceneFluidQuadrature(
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& definition,
    std::span<const SceneFluidQuadratureTraction> tractions,
    const ConservativeTransferSettings& settings = {});

} // namespace simwing::fsi
