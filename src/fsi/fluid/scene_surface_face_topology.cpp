#include "fluid/scene_surface_face_topology.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    template<typename Enum>
    void enumeration(const Enum value) {
        using Underlying = std::underlying_type_t<Enum>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

bool checkedMultiply(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

using FaceKey = std::tuple<std::uint8_t,
                           std::size_t,
                           std::size_t,
                           std::size_t>;

FaceKey faceKey(const GridFaceAxis axis,
                const std::size_t i,
                const std::size_t j,
                const std::size_t k) {
    return {static_cast<std::uint8_t>(axis), i, j, k};
}

bool internalFace(const GridCellCounts counts,
                  const GridFaceAxis axis,
                  const std::size_t i,
                  const std::size_t j,
                  const std::size_t k) {
    if (axis == GridFaceAxis::X) {
        return i > 0 && i < counts.x && j < counts.y && k < counts.z;
    }
    if (axis == GridFaceAxis::Y) {
        return j > 0 && j < counts.y && i < counts.x && k < counts.z;
    }
    return k > 0 && k < counts.z && i < counts.x && j < counts.y;
}

std::uint64_t activeFaceStableId(const PeriodicCartesianGrid& grid,
                                 const GridFaceAxis axis,
                                 const std::size_t i,
                                 const std::size_t j,
                                 const std::size_t k) {
    const GridCellCounts counts = grid.cellCounts();
    const Vector3 lower = grid.lowerMeters();
    const Vector3 upper = grid.upperMeters();
    Fingerprint fingerprint;
    fingerprint.integer(sceneFluidFaceTopologyVersion);
    fingerprint.integer(static_cast<std::uint64_t>(counts.x));
    fingerprint.integer(static_cast<std::uint64_t>(counts.y));
    fingerprint.integer(static_cast<std::uint64_t>(counts.z));
    fingerprint.real(lower.x);
    fingerprint.real(lower.y);
    fingerprint.real(lower.z);
    fingerprint.real(upper.x);
    fingerprint.real(upper.y);
    fingerprint.real(upper.z);
    fingerprint.enumeration(axis);
    fingerprint.integer(static_cast<std::uint64_t>(i));
    fingerprint.integer(static_cast<std::uint64_t>(j));
    fingerprint.integer(static_cast<std::uint64_t>(k));
    return fingerprint.value();
}

struct FaceReferences {
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::vector<std::size_t> crossings;
    std::vector<std::size_t> coplanarPatches;
};

std::uint64_t topologyFingerprint(const SceneFluidFaceTopology& topology) {
    Fingerprint fingerprint;
    fingerprint.integer(topology.version);
    fingerprint.integer(topology.surfaceDefinitionFingerprint);
    fingerprint.integer(topology.surfaceStateFingerprint);
    fingerprint.integer(topology.patchOwnershipFingerprint);
    fingerprint.integer(topology.faceCrossingFingerprint);
    fingerprint.integer(topology.structureDefinitionFingerprint);
    fingerprint.integer(topology.acceptedStepCount);
    fingerprint.real(topology.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(topology.cellCounts.z));
    fingerprint.integer(static_cast<std::uint64_t>(topology.activeFaces.size()));
    for (const auto& face : topology.activeFaces) {
        fingerprint.integer(face.stableId);
        fingerprint.enumeration(face.axis);
        fingerprint.integer(static_cast<std::uint64_t>(face.i));
        fingerprint.integer(static_cast<std::uint64_t>(face.j));
        fingerprint.integer(static_cast<std::uint64_t>(face.k));
        fingerprint.integer(static_cast<std::uint64_t>(
            face.firstCrossingReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            face.crossingReferenceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            face.firstCoplanarPatchReference));
        fingerprint.integer(static_cast<std::uint64_t>(
            face.coplanarPatchReferenceCount));
        fingerprint.real(face.summedCrossingLengthMeters);
        fingerprint.real(face.summedCoplanarAreaSquareMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        topology.crossingReferences.size()));
    for (const std::size_t reference : topology.crossingReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        topology.coplanarPatchReferences.size()));
    for (const std::size_t reference : topology.coplanarPatchReferences) {
        fingerprint.integer(static_cast<std::uint64_t>(reference));
    }
    return fingerprint.value();
}

SceneFluidFaceTopology buildTopology(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopologyLimits& limits) {
    SceneFluidFaceTopology result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.patchOwnershipFingerprint = ownership.fingerprint;
    result.faceCrossingFingerprint = crossings.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();

    std::size_t totalReferenceCount = 0;
    if (!checkedAdd(crossings.crossings.size(), ownership.facePatches.size(),
                    totalReferenceCount)
        || totalReferenceCount > limits.maximumReferences) {
        throw std::length_error(
            "scene fluid face topology exceeds its reference limit");
    }

    std::map<FaceKey, FaceReferences> referencesByFace;
    for (std::size_t index = 0; index < crossings.crossings.size(); ++index) {
        const auto& crossing = crossings.crossings[index];
        if (!internalFace(result.cellCounts, crossing.axis,
                          crossing.i, crossing.j, crossing.k)) {
            throw std::invalid_argument(
                "scene fluid crossing is not on an internal grid face");
        }
        auto& references = referencesByFace[faceKey(
            crossing.axis, crossing.i, crossing.j, crossing.k)];
        references.axis = crossing.axis;
        references.i = crossing.i;
        references.j = crossing.j;
        references.k = crossing.k;
        references.crossings.push_back(index);
    }
    for (std::size_t index = 0; index < ownership.facePatches.size(); ++index) {
        const auto& patch = ownership.facePatches[index];
        if (!internalFace(result.cellCounts, patch.axis,
                          patch.i, patch.j, patch.k)) {
            throw std::invalid_argument(
                "scene fluid coplanar patch is not on an internal grid face");
        }
        auto& references = referencesByFace[faceKey(
            patch.axis, patch.i, patch.j, patch.k)];
        references.axis = patch.axis;
        references.i = patch.i;
        references.j = patch.j;
        references.k = patch.k;
        references.coplanarPatches.push_back(index);
    }
    if (referencesByFace.size() > limits.maximumActiveFaces) {
        throw std::length_error(
            "scene fluid face topology exceeds its active-face limit");
    }

    result.activeFaces.reserve(referencesByFace.size());
    result.crossingReferences.reserve(crossings.crossings.size());
    result.coplanarPatchReferences.reserve(ownership.facePatches.size());
    std::set<std::uint64_t> stableIds;
    for (auto& [key, references] : referencesByFace) {
        static_cast<void>(key);
        std::sort(references.crossings.begin(), references.crossings.end(),
                  [&](const std::size_t first, const std::size_t second) {
                      return crossings.crossings[first].stableId
                          < crossings.crossings[second].stableId;
                  });
        std::sort(references.coplanarPatches.begin(),
                  references.coplanarPatches.end(),
                  [&](const std::size_t first, const std::size_t second) {
                      const auto& firstPatch = ownership.facePatches[first];
                      const auto& secondPatch = ownership.facePatches[second];
                      return std::pair{firstPatch.triangleId, first}
                          < std::pair{secondPatch.triangleId, second};
                  });

        SceneFluidActiveFace face;
        face.stableId = activeFaceStableId(
            grid, references.axis,
            references.i, references.j, references.k);
        if (!stableIds.insert(face.stableId).second) {
            throw std::invalid_argument(
                "scene fluid active-face stable-ID collision");
        }
        face.axis = references.axis;
        face.i = references.i;
        face.j = references.j;
        face.k = references.k;
        face.firstCrossingReference = result.crossingReferences.size();
        face.crossingReferenceCount = references.crossings.size();
        face.firstCoplanarPatchReference =
            result.coplanarPatchReferences.size();
        face.coplanarPatchReferenceCount =
            references.coplanarPatches.size();
        for (const std::size_t reference : references.crossings) {
            face.summedCrossingLengthMeters +=
                crossings.crossings[reference].lengthMeters;
            result.crossingReferences.push_back(reference);
        }
        for (const std::size_t reference : references.coplanarPatches) {
            face.summedCoplanarAreaSquareMeters +=
                ownership.facePatches[reference].areaSquareMeters;
            result.coplanarPatchReferences.push_back(reference);
        }
        if (!std::isfinite(face.summedCrossingLengthMeters)
            || !std::isfinite(face.summedCoplanarAreaSquareMeters)) {
            throw std::length_error(
                "scene fluid active-face metric is not finite");
        }
        result.activeFaces.push_back(face);
    }

    std::size_t faceBytes = 0;
    std::size_t referenceBytes = 0;
    std::size_t totalBytes = 0;
    if (!checkedMultiply(result.activeFaces.size(),
                         sizeof(SceneFluidActiveFace), faceBytes)
        || !checkedMultiply(totalReferenceCount,
                            sizeof(std::size_t), referenceBytes)
        || !checkedAdd(faceBytes, referenceBytes, totalBytes)
        || totalBytes > limits.maximumTopologyBytes) {
        throw std::length_error(
            "scene fluid face topology exceeds its storage limit");
    }
    result.fingerprint = topologyFingerprint(result);
    return result;
}

} // namespace

const SceneFluidActiveFace* SceneFluidFaceTopology::activeFace(
    const GridFaceAxis axis,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const noexcept {
    const FaceKey key = faceKey(axis, i, j, k);
    const auto found = std::lower_bound(
        activeFaces.begin(), activeFaces.end(), key,
        [](const SceneFluidActiveFace& face, const FaceKey& sought) {
            return faceKey(face.axis, face.i, face.j, face.k) < sought;
        });
    if (found == activeFaces.end()
        || faceKey(found->axis, found->i, found->j, found->k) != key) {
        return nullptr;
    }
    return &*found;
}

std::span<const std::size_t> SceneFluidFaceTopology::crossingsForFace(
    const SceneFluidActiveFace& face) const {
    if (face.firstCrossingReference > crossingReferences.size()
        || face.crossingReferenceCount
            > crossingReferences.size() - face.firstCrossingReference) {
        throw std::out_of_range(
            "scene fluid active-face crossing reference range is invalid");
    }
    return std::span<const std::size_t>(crossingReferences).subspan(
        face.firstCrossingReference, face.crossingReferenceCount);
}

std::span<const std::size_t> SceneFluidFaceTopology::coplanarPatchesForFace(
    const SceneFluidActiveFace& face) const {
    if (face.firstCoplanarPatchReference > coplanarPatchReferences.size()
        || face.coplanarPatchReferenceCount
            > coplanarPatchReferences.size()
                - face.firstCoplanarPatchReference) {
        throw std::out_of_range(
            "scene fluid active-face coplanar reference range is invalid");
    }
    return std::span<const std::size_t>(coplanarPatchReferences).subspan(
        face.firstCoplanarPatchReference,
        face.coplanarPatchReferenceCount);
}

SceneFluidFaceTopology buildSceneFluidFaceTopology(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings,
    const SceneFluidFaceTopologyLimits& limits) {
    validateSceneFluidFaceCrossings(
        crossings, surface, state, grid, candidates, intersections, patches,
        ownership);
    SceneFluidFaceTopology result = buildTopology(
        surface, state, grid, ownership, crossings, limits);
    validateSceneFluidFaceTopology(
        result, surface, state, grid, candidates, intersections, patches,
        ownership, crossings);
    return result;
}

void validateSceneFluidFaceTopology(
    const SceneFluidFaceTopology& topology,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const PeriodicCartesianGrid& grid,
    const SceneFluidGridCandidateSet& candidates,
    const SceneFluidGridIntersectionSet& intersections,
    const SceneFluidGridPatchSet& patches,
    const SceneFluidPatchOwnership& ownership,
    const SceneFluidFaceCrossingSet& crossings) {
    validateSceneFluidFaceCrossings(
        crossings, surface, state, grid, candidates, intersections, patches,
        ownership);
    if (topology.version != sceneFluidFaceTopologyVersion
        || topology.fingerprint == 0
        || topology.surfaceDefinitionFingerprint != surface.fingerprint
        || topology.surfaceStateFingerprint != state.fingerprint
        || topology.patchOwnershipFingerprint != ownership.fingerprint
        || topology.faceCrossingFingerprint != crossings.fingerprint
        || topology.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || topology.acceptedStepCount != state.acceptedStepCount
        || topology.simulationTimeSeconds != state.simulationTimeSeconds
        || topology.cellCounts != grid.cellCounts()) {
        throw std::invalid_argument(
            "scene fluid face-topology identity is invalid");
    }
    const SceneFluidFaceTopology expected = buildTopology(
        surface,
        state,
        grid,
        ownership,
        crossings,
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()});
    if (topology != expected) {
        throw std::invalid_argument(
            "scene fluid face topology does not match its source geometry");
    }
}

} // namespace simwing::fsi::fluid
