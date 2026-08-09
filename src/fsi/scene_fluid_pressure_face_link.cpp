#include "scene_fluid_pressure_face_link.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

namespace simwing::fsi {
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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
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

void validateSettings(const SceneFluidPressureFaceLinkSettings& settings) {
    if (!std::isfinite(settings.areaToleranceSquareMeters)
        || !(settings.areaToleranceSquareMeters > 0.0)) {
        throw std::invalid_argument(
            "scene fluid pressure-face-link settings are invalid");
    }
}

std::size_t faceOrdinal(const fluid::PeriodicCartesianGrid& grid,
                        const fluid::GridFaceAxis axis,
                        const std::size_t i,
                        const std::size_t j,
                        const std::size_t k) {
    return static_cast<std::size_t>(axis) * grid.cellCount()
        + grid.cellIndex(i, j, k);
}

fluid::GridFaceAxis gridFaceAxis(
    const SceneFluidOpeningPatchFaceAxis axis) {
    switch (axis) {
    case SceneFluidOpeningPatchFaceAxis::X:
        return fluid::GridFaceAxis::X;
    case SceneFluidOpeningPatchFaceAxis::Y:
        return fluid::GridFaceAxis::Y;
    case SceneFluidOpeningPatchFaceAxis::Z:
        return fluid::GridFaceAxis::Z;
    }
    throw std::invalid_argument(
        "scene fluid opening patch has invalid pressure-face axis");
}

std::uint64_t faceStableId(const fluid::GridFaceAxis axis,
                           const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x7072657373666163ULL});
    fingerprint.enumeration(axis);
    fingerprint.integer(static_cast<std::uint64_t>(i));
    fingerprint.integer(static_cast<std::uint64_t>(j));
    fingerprint.integer(static_cast<std::uint64_t>(k));
    return fingerprint.value();
}

std::uint64_t linkStableId(const std::uint64_t faceId,
                           const StableId regionId) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x70726573736c696eULL});
    fingerprint.integer(faceId);
    fingerprint.integer(regionId);
    return fingerprint.value();
}

double faceArea(const fluid::Vector3& spacing,
                const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return spacing.y * spacing.z;
    case fluid::GridFaceAxis::Y: return spacing.x * spacing.z;
    case fluid::GridFaceAxis::Z: return spacing.x * spacing.y;
    }
    throw std::invalid_argument("scene fluid pressure face has invalid axis");
}

double centerDistance(const fluid::Vector3& spacing,
                      const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return spacing.x;
    case fluid::GridFaceAxis::Y: return spacing.y;
    case fluid::GridFaceAxis::Z: return spacing.z;
    }
    throw std::invalid_argument("scene fluid pressure face has invalid axis");
}

std::size_t minusCellIndex(const fluid::PeriodicCartesianGrid& grid,
                           const fluid::GridFaceAxis axis,
                           const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) {
    const auto counts = grid.cellCounts();
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return grid.cellIndex(i == 0 ? counts.x - 1 : i - 1, j, k);
    case fluid::GridFaceAxis::Y:
        return grid.cellIndex(i, j == 0 ? counts.y - 1 : j - 1, k);
    case fluid::GridFaceAxis::Z:
        return grid.cellIndex(i, j, k == 0 ? counts.z - 1 : k - 1);
    }
    throw std::invalid_argument("scene fluid pressure face has invalid axis");
}

const SceneFluidPressureControlVolume* controlFor(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::size_t cellIndex,
    const StableId regionId) {
    const auto& cell = pressureVolumes.cells[cellIndex];
    for (std::size_t offset = 0; offset < cell.controlVolumeCount; ++offset) {
        const auto& control = pressureVolumes.controlVolumes[
            cell.firstControlVolume + offset];
        if (control.regionId == regionId) return &control;
    }
    return nullptr;
}

std::vector<StableId> commonRegions(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::size_t minusCell,
    const std::size_t plusCell) {
    std::vector<StableId> result;
    const auto& minus = pressureVolumes.cells[minusCell];
    const auto& plus = pressureVolumes.cells[plusCell];
    for (std::size_t minusOffset = 0;
         minusOffset < minus.controlVolumeCount; ++minusOffset) {
        const auto& minusControl = pressureVolumes.controlVolumes[
            minus.firstControlVolume + minusOffset];
        if (controlFor(pressureVolumes, plusCell, minusControl.regionId)) {
            result.push_back(minusControl.regionId);
        }
    }
    return result;
}

std::size_t storageBytesForCounts(const std::size_t faceCount,
                                  const std::size_t linkCount) {
    std::size_t faceBytes = 0;
    std::size_t linkBytes = 0;
    std::size_t total = 0;
    if (!checkedMultiply(faceCount, sizeof(SceneFluidPressureFace), faceBytes)
        || !checkedMultiply(linkCount, sizeof(SceneFluidPressureFaceLink),
                            linkBytes)
        || !checkedAdd(faceBytes, linkBytes, total)) {
        throw std::length_error(
            "scene fluid pressure-face-link storage size overflows");
    }
    return total;
}

std::size_t storageBytes(const SceneFluidPressureFaceLinkSet& faceLinks) {
    return storageBytesForCounts(faceLinks.faces.size(), faceLinks.links.size());
}

std::uint64_t faceLinkFingerprint(
    const SceneFluidPressureFaceLinkSet& faceLinks) {
    Fingerprint fingerprint;
    fingerprint.integer(faceLinks.version);
    for (const std::uint64_t value : {
             faceLinks.surfaceDefinitionFingerprint,
             faceLinks.surfaceStateFingerprint,
             faceLinks.gridEpochFingerprint,
             faceLinks.openingPatchFingerprint,
             faceLinks.pressureControlVolumeFingerprint,
             faceLinks.structureDefinitionFingerprint,
             faceLinks.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(faceLinks.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(faceLinks.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(faceLinks.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(faceLinks.cellCounts.z));
    for (const double value : {
             faceLinks.lowerMeters.x, faceLinks.lowerMeters.y,
             faceLinks.lowerMeters.z, faceLinks.upperMeters.x,
             faceLinks.upperMeters.y, faceLinks.upperMeters.z,
             faceLinks.settings.areaToleranceSquareMeters,
             faceLinks.totalFaceAreaSquareMeters,
             faceLinks.totalLinkedAreaSquareMeters,
             faceLinks.maximumResolvedAreaResidualSquareMeters}) {
        fingerprint.real(value);
    }
    for (const std::size_t value : {
             faceLinks.ownedStorageBytes,
             faceLinks.resolvedFullFaceCount,
             faceLinks.resolvedPartitionFaceCount,
             faceLinks.unresolvedActiveFaceCount,
             faceLinks.unresolvedAmbiguousFaceCount,
             faceLinks.unresolvedOpeningFaceCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(static_cast<std::uint64_t>(faceLinks.faces.size()));
    for (const auto& face : faceLinks.faces) {
        fingerprint.integer(static_cast<std::uint64_t>(face.faceIndex));
        fingerprint.integer(face.stableId);
        fingerprint.enumeration(face.axis);
        fingerprint.integer(static_cast<std::uint64_t>(face.i));
        fingerprint.integer(static_cast<std::uint64_t>(face.j));
        fingerprint.integer(static_cast<std::uint64_t>(face.k));
        fingerprint.integer(static_cast<std::uint64_t>(face.minusCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(face.plusCellIndex));
        fingerprint.enumeration(face.status);
        fingerprint.integer(static_cast<std::uint64_t>(face.firstLink));
        fingerprint.integer(static_cast<std::uint64_t>(face.linkCount));
        fingerprint.real(face.faceAreaSquareMeters);
        fingerprint.real(face.linkedAreaSquareMeters);
        fingerprint.real(face.areaResidualSquareMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(faceLinks.links.size()));
    for (const auto& link : faceLinks.links) {
        fingerprint.integer(static_cast<std::uint64_t>(link.linkIndex));
        fingerprint.integer(link.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(link.faceIndex));
        fingerprint.integer(link.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(link.regionIndex));
        fingerprint.integer(static_cast<std::uint64_t>(link.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            link.minusControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            link.plusControlVolumeIndex));
        fingerprint.real(link.areaSquareMeters);
        fingerprint.real(link.centerDistanceMeters);
        fingerprint.real(link.geometryWeightMeters);
    }
    return fingerprint.value();
}

SceneFluidPressureFaceLinkSet buildFaceLinks(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSettings& settings,
    const SceneFluidPressureFaceLinkLimits& limits) {
    validateSettings(settings);
    std::size_t faceCount = 0;
    if (!checkedMultiply(grid.cellCount(), std::size_t{3}, faceCount)) {
        throw std::length_error(
            "scene fluid pressure-face-link face count overflows");
    }
    if (faceCount > limits.maximumFaces) {
        throw std::length_error(
            "scene fluid pressure-face-link exceeds its face limit");
    }
    const std::size_t faceStorageBytes =
        storageBytesForCounts(faceCount, 0);
    if (faceStorageBytes > limits.maximumLinkBytes) {
        throw std::length_error(
            "scene fluid pressure-face-link exceeds its byte limit");
    }
    const std::size_t maximumLinksByBytes =
        (limits.maximumLinkBytes - faceStorageBytes)
        / sizeof(SceneFluidPressureFaceLink);
    if (pressureVolumes.surfaceDefinitionFingerprint != surface.fingerprint
        || pressureVolumes.surfaceStateFingerprint != state.fingerprint
        || pressureVolumes.cellCounts != grid.cellCounts()
        || pressureVolumes.lowerMeters != grid.lowerMeters()
        || pressureVolumes.upperMeters != grid.upperMeters()
        || pressureVolumes.cells.size() != grid.cellCount()
        || epoch.facePartitions.surfaceStateFingerprint != state.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-face-link source identity is invalid");
    }

    const auto spacing = grid.cellSpacingMeters();
    std::vector<std::size_t> activeFaceIndices(
        faceCount, std::numeric_limits<std::size_t>::max());
    for (std::size_t activeIndex = 0;
         activeIndex < epoch.faceTopology.activeFaces.size(); ++activeIndex) {
        const auto& active = epoch.faceTopology.activeFaces[activeIndex];
        const std::size_t ordinal = faceOrdinal(
            grid, active.axis, active.i, active.j, active.k);
        if (activeFaceIndices[ordinal]
            != std::numeric_limits<std::size_t>::max()) {
            throw std::invalid_argument(
                "scene fluid pressure-face-link active face is duplicated");
        }
        activeFaceIndices[ordinal] = activeIndex;
    }
    std::vector<std::size_t> partitionsByActiveFace(
        epoch.faceTopology.activeFaces.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::size_t partitionIndex = 0;
         partitionIndex < epoch.facePartitions.partitions.size();
         ++partitionIndex) {
        const auto& partition =
            epoch.facePartitions.partitions[partitionIndex];
        if (partition.activeFaceIndex >= partitionsByActiveFace.size()
            || partitionsByActiveFace[partition.activeFaceIndex]
                != std::numeric_limits<std::size_t>::max()) {
            throw std::invalid_argument(
                "scene fluid pressure-face-link partition index is invalid");
        }
        partitionsByActiveFace[partition.activeFaceIndex] = partitionIndex;
    }
    std::vector<std::uint8_t> openingFaces(faceCount, 0);
    for (const auto& patch : openingPatches.patches) {
        if (patch.ownerKind != SceneFluidOpeningPatchOwnerKind::Face) {
            continue;
        }
        const auto axis = gridFaceAxis(patch.faceAxis);
        const std::size_t ordinal = faceOrdinal(
            grid, axis, patch.faceI, patch.faceJ, patch.faceK);
        openingFaces[ordinal] = 1;
    }

    SceneFluidPressureFaceLinkSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.gridEpochFingerprint = epoch.fingerprint;
    result.openingPatchFingerprint = openingPatches.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    result.faces.reserve(faceCount);
    std::size_t reserveLinkCount = 0;
    if (!checkedMultiply(pressureVolumes.controlVolumes.size(),
                         std::size_t{3}, reserveLinkCount)) {
        reserveLinkCount = limits.maximumLinks;
    }
    result.links.reserve(std::min({
        limits.maximumLinks, maximumLinksByBytes, reserveLinkCount}));
    std::unordered_set<std::uint64_t> faceStableIds;
    std::unordered_set<std::uint64_t> linkStableIds;
    faceStableIds.reserve(faceCount);

    const auto appendLink = [&](SceneFluidPressureFace& face,
                                const StableId regionId,
                                const double areaSquareMeters) {
        if (!(areaSquareMeters > 0.0) || !std::isfinite(areaSquareMeters)) {
            throw std::invalid_argument(
                "scene fluid pressure face link has invalid area");
        }
        const auto* minus = controlFor(
            pressureVolumes, face.minusCellIndex, regionId);
        const auto* plus = controlFor(
            pressureVolumes, face.plusCellIndex, regionId);
        if (!minus || !plus || minus->regionIndex != plus->regionIndex
            || minus->componentIndex != plus->componentIndex) {
            throw std::invalid_argument(
                "scene fluid pressure face area has no matching control volumes");
        }
        if (result.links.size() == limits.maximumLinks) {
            throw std::length_error(
                "scene fluid pressure-face-link exceeds its link limit");
        }
        if (result.links.size() == maximumLinksByBytes) {
            throw std::length_error(
                "scene fluid pressure-face-link exceeds its byte limit");
        }
        SceneFluidPressureFaceLink link;
        link.linkIndex = result.links.size();
        link.stableId = linkStableId(face.stableId, regionId);
        link.faceIndex = face.faceIndex;
        link.regionId = regionId;
        link.regionIndex = minus->regionIndex;
        link.componentIndex = minus->componentIndex;
        link.minusControlVolumeIndex = minus->controlVolumeIndex;
        link.plusControlVolumeIndex = plus->controlVolumeIndex;
        link.areaSquareMeters = areaSquareMeters;
        link.centerDistanceMeters = centerDistance(spacing, face.axis);
        link.geometryWeightMeters =
            link.areaSquareMeters / link.centerDistanceMeters;
        if (!linkStableIds.insert(link.stableId).second) {
            throw std::invalid_argument(
                "scene fluid pressure face-link stable ID collides");
        }
        face.linkedAreaSquareMeters += areaSquareMeters;
        result.links.push_back(link);
    };

    const auto counts = grid.cellCounts();
    for (const fluid::GridFaceAxis axis : {
             fluid::GridFaceAxis::X,
             fluid::GridFaceAxis::Y,
             fluid::GridFaceAxis::Z}) {
        for (std::size_t k = 0; k < counts.z; ++k) {
            for (std::size_t j = 0; j < counts.y; ++j) {
                for (std::size_t i = 0; i < counts.x; ++i) {
                    SceneFluidPressureFace face;
                    face.faceIndex = result.faces.size();
                    face.stableId = faceStableId(axis, i, j, k);
                    face.axis = axis;
                    face.i = i;
                    face.j = j;
                    face.k = k;
                    face.minusCellIndex = minusCellIndex(
                        grid, axis, i, j, k);
                    face.plusCellIndex = grid.cellIndex(i, j, k);
                    face.firstLink = result.links.size();
                    face.faceAreaSquareMeters = faceArea(spacing, axis);
                    if (!faceStableIds.insert(face.stableId).second) {
                        throw std::invalid_argument(
                            "scene fluid pressure-face stable ID collides");
                    }
                    const std::size_t ordinal = faceOrdinal(
                        grid, axis, i, j, k);
                    const std::size_t activeIndex =
                        activeFaceIndices[ordinal];
                    if (openingFaces[ordinal] != 0) {
                        face.status = SceneFluidPressureFaceStatus::
                            UnresolvedOpening;
                        ++result.unresolvedOpeningFaceCount;
                    } else if (activeIndex
                        == std::numeric_limits<std::size_t>::max()) {
                        const auto common = commonRegions(
                            pressureVolumes,
                            face.minusCellIndex, face.plusCellIndex);
                        if (common.size() == 1) {
                            face.status =
                                SceneFluidPressureFaceStatus::ResolvedFull;
                            appendLink(face, common.front(),
                                       face.faceAreaSquareMeters);
                            ++result.resolvedFullFaceCount;
                        } else {
                            face.status = SceneFluidPressureFaceStatus::
                                UnresolvedAmbiguous;
                            ++result.unresolvedAmbiguousFaceCount;
                        }
                    } else {
                        const std::size_t partitionIndex =
                            partitionsByActiveFace[activeIndex];
                        if (partitionIndex
                            == std::numeric_limits<std::size_t>::max()) {
                            face.status =
                                SceneFluidPressureFaceStatus::UnresolvedActive;
                            ++result.unresolvedActiveFaceCount;
                        } else {
                            face.status = SceneFluidPressureFaceStatus::
                                ResolvedPartition;
                            const auto& partition =
                                epoch.facePartitions.partitions[partitionIndex];
                            for (std::size_t offset = 0;
                                 offset < partition.regionAreaCount; ++offset) {
                                const auto& regionArea =
                                    epoch.facePartitions.regionAreas[
                                        partition.firstRegionArea + offset];
                                if (regionArea.areaSquareMeters == 0.0) continue;
                                appendLink(face, regionArea.regionId,
                                           regionArea.areaSquareMeters);
                            }
                            ++result.resolvedPartitionFaceCount;
                        }
                    }
                    face.linkCount = result.links.size() - face.firstLink;
                    face.areaResidualSquareMeters =
                        face.linkedAreaSquareMeters
                        - face.faceAreaSquareMeters;
                    if (face.status
                            == SceneFluidPressureFaceStatus::ResolvedFull
                        || face.status
                            == SceneFluidPressureFaceStatus::ResolvedPartition) {
                        result.maximumResolvedAreaResidualSquareMeters =
                            std::max(
                                result.maximumResolvedAreaResidualSquareMeters,
                                std::abs(face.areaResidualSquareMeters));
                        if (face.linkCount == 0
                            || std::abs(face.areaResidualSquareMeters)
                                > settings.areaToleranceSquareMeters) {
                            throw std::invalid_argument(
                                "scene fluid pressure face links do not close area");
                        }
                    } else if (face.linkCount != 0) {
                        throw std::logic_error(
                            "unresolved scene fluid pressure face owns links");
                    }
                    result.totalFaceAreaSquareMeters +=
                        face.faceAreaSquareMeters;
                    result.totalLinkedAreaSquareMeters +=
                        face.linkedAreaSquareMeters;
                    result.faces.push_back(face);
                }
            }
        }
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumLinkBytes) {
        throw std::length_error(
            "scene fluid pressure-face-link exceeds its byte limit");
    }
    result.fingerprint = faceLinkFingerprint(result);
    return result;
}

} // namespace

SceneFluidPressureFaceLinkSet buildSceneFluidPressureFaceLinks(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSettings& settings,
    const SceneFluidPressureFaceLinkLimits& limits) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidOpeningGridPatches(
        openingPatches, surface, state, caps, openingQuadrature, grid);
    validateSceneFluidCellVolumes(
        volumes, surface, state, grid, transfer, epoch);
    validateSceneFluidPressureControlVolumes(
        pressureVolumes, surface, volumes, connectivity);
    auto result = buildFaceLinks(
        surface, state, grid, epoch, openingPatches,
        pressureVolumes, settings, limits);
    validateSceneFluidPressureFaceLinks(
        result, surface, state, grid, transfer, epoch, caps,
        openingQuadrature, openingPatches, volumes, connectivity,
        pressureVolumes);
    return result;
}

void validateSceneFluidPressureFaceLinks(
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidOpeningGridPatches(
        openingPatches, surface, state, caps, openingQuadrature, grid);
    validateSceneFluidCellVolumes(
        volumes, surface, state, grid, transfer, epoch);
    validateSceneFluidPressureControlVolumes(
        pressureVolumes, surface, volumes, connectivity);
    if (faceLinks.version != sceneFluidPressureFaceLinkVersion
        || faceLinks.fingerprint == 0
        || faceLinks.surfaceDefinitionFingerprint != surface.fingerprint
        || faceLinks.surfaceStateFingerprint != state.fingerprint
        || faceLinks.gridEpochFingerprint != epoch.fingerprint
        || faceLinks.openingPatchFingerprint != openingPatches.fingerprint
        || faceLinks.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-face-link identity is invalid");
    }
    const SceneFluidPressureFaceLinkLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildFaceLinks(
        surface, state, grid, epoch, openingPatches, pressureVolumes,
        faceLinks.settings, unlimited);
    if (faceLinks != expected
        || faceLinks.ownedStorageBytes != storageBytes(faceLinks)
        || faceLinks.fingerprint != faceLinkFingerprint(faceLinks)) {
        throw std::invalid_argument(
            "scene fluid pressure-face-link payload is invalid");
    }
}

} // namespace simwing::fsi
