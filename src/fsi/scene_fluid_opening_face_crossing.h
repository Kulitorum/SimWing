#pragma once

#include "fluid/interface_jump.h"
#include "scene_fluid_opening_patch.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidOpeningFaceCrossingVersion = 1;

struct SceneFluidOpeningFaceCrossingLimits {
    std::size_t maximumCandidateSegments = 30'000'000;
    std::size_t maximumCrossings = 10'000'000;
    std::size_t maximumCrossingBytes = 768ULL * 1024ULL * 1024ULL;
};

struct SceneFluidOpeningFaceCrossing {
    std::uint64_t stableId = 0;
    std::size_t lowerCellPatchIndex = 0;
    std::size_t upperCellPatchIndex = 0;
    fluid::GridFaceAxis axis = fluid::GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    std::size_t triangleOrdinal = 0;
    std::uint64_t sourcePointStableId = 0;
    fluid::SceneFluidClippedVertex first;
    fluid::SceneFluidClippedVertex second;
    Vec3 midpointMeters;
    // Projection of the authored opening-cap normal into this face. Together
    // with the positive Cartesian face normal it fixes the directed segment
    // whose left side is the authored positive region.
    Vec3 negativeToPositiveDirectionInFace;
    double lengthMeters = 0.0;
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;

    bool operator==(const SceneFluidOpeningFaceCrossing& other) const;
};

// Exact positive-length intersections of virtual opening-cap triangles with
// internal Cartesian faces. Each transverse segment must appear independently
// in the clipped polygons of both adjacent cells before one canonical crossing
// is published. Face-owned cap area remains an aperture patch and is counted,
// never converted into a crossing. Unpaired point/edge contact remains an
// explicit diagnostic; grid-edge-aligned paired segments reject until there is
// a dedicated edge owner.
struct SceneFluidOpeningFaceCrossingSet {
    std::uint32_t version = sceneFluidOpeningFaceCrossingVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t openingCapFingerprint = 0;
    std::uint64_t openingQuadratureFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t candidateSegmentCount = 0;
    std::size_t unpairedContactSegmentCount = 0;
    std::size_t faceOwnedPatchCount = 0;
    std::size_t ownedStorageBytes = 0;
    double crossingLengthMeters = 0.0;
    std::vector<SceneFluidOpeningFaceCrossing> crossings;

    bool operator==(const SceneFluidOpeningFaceCrossingSet&) const = default;
};

[[nodiscard]] SceneFluidOpeningFaceCrossingSet
buildSceneFluidOpeningFaceCrossings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidOpeningFaceCrossingLimits& limits = {});

void validateSceneFluidOpeningFaceCrossings(
    const SceneFluidOpeningFaceCrossingSet& crossings,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& quadrature,
    const SceneFluidOpeningGridPatchSet& patches,
    const fluid::PeriodicCartesianGrid& grid);

} // namespace simwing::fsi
