#pragma once

#include "fluid/scene_surface_ownership.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFaceCrossingVersion = 1;

struct SceneFluidFaceCrossing {
    std::uint64_t stableId = 0;
    std::size_t lowerCellSourcePatchIndex = 0;
    std::size_t upperCellSourcePatchIndex = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    // Existing MAC convention: (i,j,k) is the cell on the positive side.
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t triangleIndex = 0;
    StableId triangleId = invalidStableId;
    SceneFluidClippedVertex first;
    SceneFluidClippedVertex second;
    Vec3 midpointMeters;
    // Unit vector in the face plane, following the triangle normal from its
    // authored negative-side region toward its positive-side region.
    Vec3 negativeToPositiveDirectionInFace;
    double lengthMeters = 0.0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    StableId materialId = invalidStableId;
    StableId sheetId = invalidStableId;
    SurfaceRole role = SurfaceRole::Skin;

    bool operator==(const SceneFluidFaceCrossing& other) const {
        return stableId == other.stableId
            && lowerCellSourcePatchIndex == other.lowerCellSourcePatchIndex
            && upperCellSourcePatchIndex == other.upperCellSourcePatchIndex
            && axis == other.axis
            && i == other.i && j == other.j && k == other.k
            && triangleIndex == other.triangleIndex
            && triangleId == other.triangleId
            && first == other.first && second == other.second
            && midpointMeters.x == other.midpointMeters.x
            && midpointMeters.y == other.midpointMeters.y
            && midpointMeters.z == other.midpointMeters.z
            && negativeToPositiveDirectionInFace.x
                == other.negativeToPositiveDirectionInFace.x
            && negativeToPositiveDirectionInFace.y
                == other.negativeToPositiveDirectionInFace.y
            && negativeToPositiveDirectionInFace.z
                == other.negativeToPositiveDirectionInFace.z
            && lengthMeters == other.lengthMeters
            && negativeSideRegionId == other.negativeSideRegionId
            && positiveSideRegionId == other.positiveSideRegionId
            && materialId == other.materialId
            && sheetId == other.sheetId
            && role == other.role;
    }
};

struct SceneFluidFaceCrossingLimits {
    std::size_t maximumCandidateSegments = 15'000'000;
    std::size_t maximumCrossings = 5'000'000;
    std::size_t maximumCrossingBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFaceCrossingSet {
    std::uint32_t version = sceneFluidFaceCrossingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t patchOwnershipFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t candidateSegmentCount = 0;
    std::size_t unpairedContactSegmentCount = 0;
    std::size_t coplanarAreaPatchCount = 0;
    double crossingLengthMeters = 0.0;
    std::vector<SceneFluidFaceCrossing> crossings;

    bool operator==(const SceneFluidFaceCrossingSet&) const = default;
};

// Extracts positive-length intersections between ordinary owned triangle
// patches and internal Cartesian faces. A true transverse crossing appears as
// one exact segment from each adjacent cell and becomes one canonical crossing.
// An unpaired segment is retained only as a contact count; triangle area that
// lies in a face remains owned by SceneFluidOwnedFacePatch and is not converted
// into a line crossing. Grid-edge-aligned paired segments are rejected until a
// dedicated edge owner exists.
[[nodiscard]] SceneFluidFaceCrossingSet buildSceneFluidFaceCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingLimits& limits = {});

void validateSceneFluidFaceCrossings(
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership);

} // namespace simwing::fsi::fluid
