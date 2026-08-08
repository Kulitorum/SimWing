#pragma once

#include "fluid/scene_surface_crossings.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t sceneFluidFaceTopologyVersion = 1;

struct SceneFluidActiveFace {
    std::uint64_t stableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    // Existing MAC convention: (i,j,k) is the cell on the positive side.
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t firstCrossingReference = 0;
    std::size_t crossingReferenceCount = 0;
    std::size_t firstCoplanarPatchReference = 0;
    std::size_t coplanarPatchReferenceCount = 0;
    double summedCrossingLengthMeters = 0.0;
    double summedCoplanarAreaSquareMeters = 0.0;

    bool operator==(const SceneFluidActiveFace&) const = default;
};

struct SceneFluidFaceTopologyLimits {
    std::size_t maximumActiveFaces = 5'000'000;
    std::size_t maximumReferences = 15'000'000;
    std::size_t maximumTopologyBytes = 512ULL * 1024ULL * 1024ULL;
};

struct SceneFluidFaceTopology {
    std::uint32_t version = sceneFluidFaceTopologyVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t surfaceStateFingerprint = 0;
    std::uint64_t patchOwnershipFingerprint = 0;
    std::uint64_t faceCrossingFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    GridCellCounts cellCounts;
    std::vector<SceneFluidActiveFace> activeFaces;
    // Indices into SceneFluidFaceCrossingSet::crossings and
    // SceneFluidPatchOwnership::facePatches respectively.
    std::vector<std::size_t> crossingReferences;
    std::vector<std::size_t> coplanarPatchReferences;

    [[nodiscard]] const SceneFluidActiveFace* activeFace(
        GridFaceAxis axis,
        std::size_t i,
        std::size_t j,
        std::size_t k) const noexcept;
    [[nodiscard]] std::span<const std::size_t> crossingsForFace(
        const SceneFluidActiveFace& face) const;
    [[nodiscard]] std::span<const std::size_t> coplanarPatchesForFace(
        const SceneFluidActiveFace& face) const;

    bool operator==(const SceneFluidFaceTopology&) const = default;
};

// Builds a canonical sparse index over every internal MAC face touched by a
// transverse crossing or a coplanar owned area patch. Multiple sheets remain
// separate references; summed length/area are diagnostics, not union measures
// or fluid-region labels.
[[nodiscard]] SceneFluidFaceTopology buildSceneFluidFaceTopology(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopologyLimits& limits = {});

void validateSceneFluidFaceTopology(
    const SceneFluidFaceTopology& topology,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings);

} // namespace simwing::fsi::fluid
