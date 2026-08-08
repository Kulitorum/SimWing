#pragma once

#include "fluid/interface_jump.h"
#include "fluid/scene_surface_clipping.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidPatchOwnershipVersion = 1;

struct SceneFluidOwnedCellPatch {
    std::size_t sourcePatchIndex = 0;
    std::size_t cellIndex = 0;
    GridCellCoordinate cell;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    double areaSquareMeters = 0.0;

    bool operator==(const SceneFluidOwnedCellPatch&) const = default;
};

struct SceneFluidOwnedFacePatch {
    std::size_t lowerCellSourcePatchIndex = 0;
    std::size_t upperCellSourcePatchIndex = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    // Existing MAC convention: (i,j,k) is the cell on the positive side.
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    // +1 means triangle winding points from the lower coordinate cell to the
    // upper coordinate cell; -1 means it points the other way.
    std::int8_t triangleNormalAxisSign = 0;
    double areaSquareMeters = 0.0;

    bool operator==(const SceneFluidOwnedFacePatch&) const = default;
};

struct SceneFluidPatchOwnershipLimits {
    std::size_t maximumOwnedAreaPatches = 5'000'000;
    std::size_t maximumOwnershipBytes = 256ULL * 1024ULL * 1024ULL;
};

struct SceneFluidPatchOwnership {
    std::uint32_t version = sceneFluidPatchOwnershipVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t patchSetFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t pointContactPatchCount = 0;
    std::size_t segmentContactPatchCount = 0;
    double ownedAreaSquareMeters = 0.0;
    std::vector<SceneFluidOwnedCellPatch> cellPatches;
    std::vector<SceneFluidOwnedFacePatch> facePatches;

    bool operator==(const SceneFluidPatchOwnership&) const = default;
};

// Resolves positive-area patches into unique owners. Ordinary clipped area is
// cell-owned. An internal grid-plane area must appear as one exact lower/upper
// pair and becomes one face-owned patch. Positive area on a periodic domain
// boundary is rejected until explicit image ownership exists.
[[nodiscard]] SceneFluidPatchOwnership ownSceneFluidSurfacePatches(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnershipLimits& limits = {});

void validateSceneFluidPatchOwnership(
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches);

} // namespace simwing::fsi::fluid
