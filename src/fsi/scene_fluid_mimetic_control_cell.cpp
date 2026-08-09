#include "scene_fluid_mimetic_control_cell.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    void integer(const std::uint64_t value) {
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value >> (8 * byte));
            value_ *= fnvPrime;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    template<typename Enum>
    void enumeration(const Enum value) {
        integer(static_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

fluid::Vector3 vector(const Vec3& value) {
    return {value.x, value.y, value.z};
}

fluid::Vector3 vector(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

fluid::Vector3 negate(const fluid::Vector3& value) {
    return {-value.x, -value.y, -value.z};
}

double dot(const fluid::Vector3& first,
           const fluid::Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

fluid::Vector3 cross(const fluid::Vector3& first,
                     const fluid::Vector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double component(const fluid::Vector3& value,
                 const std::size_t axis) {
    return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

void addComponent(fluid::Vector3& value,
                  const std::size_t axis,
                  const double increment) {
    if (axis == 0) {
        value.x += increment;
    } else if (axis == 1) {
        value.y += increment;
    } else {
        value.z += increment;
    }
}

bool sameGrid(const fluid::GridCellCounts& first,
              const fluid::GridCellCounts& second) {
    return first == second;
}

void validateSettings(
    const SceneFluidMimeticControlCellSettings& settings) {
    if (!std::isfinite(
            settings.absoluteAreaClosureToleranceSquareMeters)
        || settings.absoluteAreaClosureToleranceSquareMeters < 0.0
        || !std::isfinite(
            settings.absoluteDivergenceTheoremToleranceCubicMeters)
        || settings.absoluteDivergenceTheoremToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeGeometryTolerance)
        || settings.relativeGeometryTolerance < 0.0
        || !std::isfinite(settings.unitNormalTolerance)
        || settings.unitNormalTolerance < 0.0) {
        throw std::invalid_argument(
            "invalid scene fluid mimetic control-cell settings");
    }
}

std::size_t checkedStorageBytes(const std::size_t controlCount,
                                const std::size_t halfFaceCount) {
    if (controlCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticControlCell)
        || halfFaceCount > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidMimeticHalfFace)) {
        throw std::length_error(
            "scene fluid mimetic control-cell storage overflows");
    }
    const std::size_t controlBytes = controlCount
        * sizeof(SceneFluidMimeticControlCell);
    const std::size_t halfFaceBytes = halfFaceCount
        * sizeof(SceneFluidMimeticHalfFace);
    if (halfFaceBytes > std::numeric_limits<std::size_t>::max()
            - controlBytes) {
        throw std::length_error(
            "scene fluid mimetic control-cell storage overflows");
    }
    return controlBytes + halfFaceBytes;
}

std::uint64_t halfFaceStableId(
    const SceneFluidMimeticHalfFaceKind kind,
    const SceneFluidMimeticHalfFaceSide side,
    const std::uint64_t sourceStableId,
    const std::uint64_t controlStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(0x4d494d4554494346ULL);
    fingerprint.enumeration(kind);
    fingerprint.enumeration(side);
    fingerprint.integer(sourceStableId);
    fingerprint.integer(controlStableId);
    const std::uint64_t result = fingerprint.value();
    return result == 0 ? 1 : result;
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticControlCellSet& set) {
    Fingerprint fingerprint;
    fingerprint.integer(set.version);
    fingerprint.integer(set.surfaceDefinitionFingerprint);
    fingerprint.integer(set.surfaceStateFingerprint);
    fingerprint.integer(set.gridEpochFingerprint);
    fingerprint.integer(set.openingPatchFingerprint);
    fingerprint.integer(set.pressureControlVolumeFingerprint);
    fingerprint.integer(set.pressureFaceLinkFingerprint);
    fingerprint.integer(set.structureDefinitionFingerprint);
    fingerprint.integer(set.acceptedStepCount);
    fingerprint.real(set.simulationTimeSeconds);
    fingerprint.integer(set.cellCounts.x);
    fingerprint.integer(set.cellCounts.y);
    fingerprint.integer(set.cellCounts.z);
    fingerprint.real(set.lowerMeters.x);
    fingerprint.real(set.lowerMeters.y);
    fingerprint.real(set.lowerMeters.z);
    fingerprint.real(set.upperMeters.x);
    fingerprint.real(set.upperMeters.y);
    fingerprint.real(set.upperMeters.z);
    fingerprint.real(
        set.settings.absoluteAreaClosureToleranceSquareMeters);
    fingerprint.real(
        set.settings.absoluteDivergenceTheoremToleranceCubicMeters);
    fingerprint.real(set.settings.relativeGeometryTolerance);
    fingerprint.real(set.settings.unitNormalTolerance);
    fingerprint.integer(set.ownedStorageBytes);
    fingerprint.integer(set.readyControlCellCount);
    fingerprint.integer(set.incompleteTopologyControlCellCount);
    fingerprint.integer(set.nonclosingControlCellCount);
    fingerprint.integer(set.unresolvedCartesianFaceCount);
    fingerprint.integer(set.omittedZeroVolumeMaterialSideCount);
    fingerprint.integer(set.missingOpeningControlSideCount);
    fingerprint.integer(set.cartesianHalfFaceCount);
    fingerprint.integer(set.materialWallHalfFaceCount);
    fingerprint.integer(set.openingHalfFaceCount);
    fingerprint.integer(set.maximumHalfFaceCountPerControl);
    fingerprint.real(set.maximumAreaClosureErrorSquareMeters);
    fingerprint.real(set.maximumDivergenceTheoremErrorCubicMeters);
    fingerprint.integer(set.controlCells.size());
    for (const auto& cell : set.controlCells) {
        fingerprint.integer(cell.controlCellIndex);
        fingerprint.integer(cell.controlVolumeIndex);
        fingerprint.integer(cell.stableId);
        fingerprint.integer(cell.cellIndex);
        fingerprint.integer(cell.regionId);
        fingerprint.integer(cell.componentIndex);
        fingerprint.integer(cell.firstHalfFace);
        fingerprint.integer(cell.halfFaceCount);
        fingerprint.integer(cell.cartesianHalfFaceCount);
        fingerprint.integer(cell.materialWallHalfFaceCount);
        fingerprint.integer(cell.openingHalfFaceCount);
        fingerprint.integer(cell.unresolvedCartesianIncidentCount);
        fingerprint.integer(cell.missingOtherControlCount);
        fingerprint.real(cell.volumeCubicMeters);
        fingerprint.real(cell.centroidMeters.x);
        fingerprint.real(cell.centroidMeters.y);
        fingerprint.real(cell.centroidMeters.z);
        fingerprint.real(cell.areaClosureVectorSquareMeters.x);
        fingerprint.real(cell.areaClosureVectorSquareMeters.y);
        fingerprint.real(cell.areaClosureVectorSquareMeters.z);
        for (const double value : cell.normalMomentMatrixCubicMeters) {
            fingerprint.real(value);
        }
        fingerprint.real(cell.maximumAreaClosureErrorSquareMeters);
        fingerprint.real(cell.maximumDivergenceTheoremErrorCubicMeters);
        fingerprint.integer(cell.topologyComplete ? 1 : 0);
        fingerprint.integer(cell.areaVectorClosed ? 1 : 0);
        fingerprint.integer(cell.divergenceTheoremClosed ? 1 : 0);
        fingerprint.integer(cell.readyForMimeticOperator ? 1 : 0);
    }
    fingerprint.integer(set.halfFaces.size());
    for (const auto& face : set.halfFaces) {
        fingerprint.integer(face.halfFaceIndex);
        fingerprint.integer(face.stableId);
        fingerprint.enumeration(face.kind);
        fingerprint.enumeration(face.side);
        fingerprint.integer(face.sourceIndex);
        fingerprint.integer(face.sourceStableId);
        fingerprint.integer(face.traceStableId);
        fingerprint.integer(face.controlVolumeIndex);
        fingerprint.integer(face.otherControlVolumeIndex);
        fingerprint.real(face.areaSquareMeters);
        fingerprint.real(face.centroidMeters.x);
        fingerprint.real(face.centroidMeters.y);
        fingerprint.real(face.centroidMeters.z);
        fingerprint.real(face.periodicImageOffsetMeters.x);
        fingerprint.real(face.periodicImageOffsetMeters.y);
        fingerprint.real(face.periodicImageOffsetMeters.z);
        fingerprint.real(face.outwardUnitNormal.x);
        fingerprint.real(face.outwardUnitNormal.y);
        fingerprint.real(face.outwardUnitNormal.z);
    }
    return fingerprint.value();
}

fluid::Vector3 unwrapNear(const fluid::Vector3& position,
                          const fluid::Vector3& reference,
                          const fluid::Vector3& lower,
                          const fluid::Vector3& upper,
                          fluid::Vector3& imageOffset) {
    fluid::Vector3 result = position;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double length = component(upper, axis)
            - component(lower, axis);
        double delta = component(result, axis)
            - component(reference, axis);
        double offset = 0.0;
        if (delta > 0.5 * length) {
            offset = -length;
        } else if (delta < -0.5 * length) {
            offset = length;
        }
        addComponent(result, axis, offset);
        addComponent(imageOffset, axis, offset);
    }
    return result;
}

fluid::Vector3 materialNormal(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const StableId triangleId) {
    const auto triangleIndex = surface.mappings.triangleIndex(triangleId);
    if (!triangleIndex) {
        throw std::invalid_argument(
            "mimetic material half-face references a foreign triangle");
    }
    const auto& triangle = surface.triangles[*triangleIndex];
    const fluid::Vector3 first = vector(
        state.vertices[triangle.vertexIndices[0]].positionMeters);
    const fluid::Vector3 second = vector(
        state.vertices[triangle.vertexIndices[1]].positionMeters);
    const fluid::Vector3 third = vector(
        state.vertices[triangle.vertexIndices[2]].positionMeters);
    fluid::Vector3 normal = cross(
        subtract(second, first), subtract(third, first));
    const double magnitude = std::sqrt(dot(normal, normal));
    if (!std::isfinite(magnitude) || magnitude <= 0.0) {
        throw std::invalid_argument(
            "mimetic material half-face has a degenerate triangle");
    }
    normal.x /= magnitude;
    normal.y /= magnitude;
    normal.z /= magnitude;
    return normal;
}

bool resolvedFaceStatus(const SceneFluidPressureFaceStatus status) {
    return status == SceneFluidPressureFaceStatus::ResolvedFull
        || status == SceneFluidPressureFaceStatus::ResolvedPartition
        || status == SceneFluidPressureFaceStatus::ResolvedOpening;
}

void validateSourceBindings(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks) {
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidQuadratureDefinition(epoch.quadrature);
    validateSceneFluidOpeningQuadrature(
        openingQuadrature, surface, state, caps);
    validateSceneFluidOpeningGridPatches(
        openingPatches, surface, state, caps, openingQuadrature, grid);
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    if (epoch.fingerprint == 0
        || epoch.quadrature.surfaceDefinitionFingerprint
            != surface.fingerprint
        || epoch.quadrature.surfaceStateFingerprint != state.fingerprint
        || epoch.quadrature.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || epoch.quadrature.acceptedStepCount != state.acceptedStepCount
        || epoch.quadrature.simulationTimeSeconds
            != state.simulationTimeSeconds
        || openingPatches.fingerprint != faceLinks.openingPatchFingerprint
        || pressureVolumes.fingerprint
            != faceLinks.pressureControlVolumeFingerprint
        || faceLinks.gridEpochFingerprint != epoch.fingerprint
        || pressureVolumes.surfaceDefinitionFingerprint
            != surface.fingerprint
        || pressureVolumes.surfaceStateFingerprint != state.fingerprint
        || !sameGrid(pressureVolumes.cellCounts, grid.cellCounts())
        || pressureVolumes.lowerMeters != grid.lowerMeters()
        || pressureVolumes.upperMeters != grid.upperMeters()
        || !sameGrid(faceLinks.cellCounts, grid.cellCounts())
        || faceLinks.lowerMeters != grid.lowerMeters()
        || faceLinks.upperMeters != grid.upperMeters()
        || pressureVolumes.acceptedStepCount != state.acceptedStepCount
        || pressureVolumes.simulationTimeSeconds
            != state.simulationTimeSeconds
        || faceLinks.acceptedStepCount != state.acceptedStepCount
        || faceLinks.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic source binding is invalid");
    }
}

struct WorkingControlCell {
    std::vector<SceneFluidMimeticHalfFace> halfFaces;
    std::size_t unresolvedCartesianIncidentCount = 0;
    std::size_t missingOtherControlCount = 0;
};

} // namespace

SceneFluidMimeticControlCellSet buildSceneFluidMimeticControlCells(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticControlCellSettings& settings,
    const SceneFluidMimeticControlCellLimits& limits) {
    validateSettings(settings);
    validateSourceBindings(
        surface, state, grid, epoch, caps, openingQuadrature,
        openingPatches, pressureVolumes, faceLinks);
    const std::size_t controlCount = pressureVolumes.controlVolumes.size();
    if (controlCount > limits.maximumControlCells) {
        throw std::length_error(
            "scene fluid mimetic control-cell limit exceeded");
    }
    if (limits.maximumHalfFaces == 0 || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "invalid scene fluid mimetic control-cell limits");
    }
    if (checkedStorageBytes(controlCount, 0)
        > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic control-cell byte limit exceeded");
    }

    SceneFluidMimeticControlCellSet result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.gridEpochFingerprint = epoch.fingerprint;
    result.openingPatchFingerprint = openingPatches.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.structureDefinitionFingerprint = state.structureDefinitionFingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;

    std::vector<WorkingControlCell> working(controlCount);
    std::size_t appendedHalfFaceCount = 0;
    const auto controlIndex = [&](const std::size_t cellIndex,
                                  const StableId regionId)
        -> std::optional<std::size_t> {
        if (cellIndex >= pressureVolumes.cells.size()
            || pressureVolumes.cells[cellIndex].cellIndex != cellIndex) {
            throw std::invalid_argument(
                "mimetic control lookup has invalid pressure-cell ordering");
        }
        const auto& cell = pressureVolumes.cells[cellIndex];
        for (std::size_t offset = 0; offset < cell.controlVolumeCount;
             ++offset) {
            const std::size_t index = cell.firstControlVolume + offset;
            if (pressureVolumes.controlVolumes[index].regionId == regionId) {
                return index;
            }
        }
        return std::nullopt;
    };
    const auto append = [&](const std::size_t ownerIndex,
                            const std::optional<std::size_t> otherIndex,
                            const SceneFluidMimeticHalfFaceKind kind,
                            const SceneFluidMimeticHalfFaceSide side,
                            const std::size_t sourceIndex,
                            const std::uint64_t sourceStableId,
                            const double area,
                            const fluid::Vector3 sourceCentroid,
                            const fluid::Vector3 outwardNormal) {
        if (ownerIndex >= controlCount || sourceStableId == 0
            || !std::isfinite(area) || area <= 0.0
            || !finite(sourceCentroid) || !finite(outwardNormal)) {
            throw std::invalid_argument(
                "invalid scene fluid mimetic half-face source");
        }
        const double normalLength = std::sqrt(dot(
            outwardNormal, outwardNormal));
        if (!std::isfinite(normalLength)
            || std::abs(normalLength - 1.0) > settings.unitNormalTolerance) {
            throw std::invalid_argument(
                "scene fluid mimetic half-face normal is not unit length");
        }
        if (appendedHalfFaceCount >= limits.maximumHalfFaces) {
            throw std::length_error(
                "scene fluid mimetic half-face limit exceeded");
        }
        if (checkedStorageBytes(
                controlCount, appendedHalfFaceCount + 1)
            > limits.maximumOwnedBytes) {
            throw std::length_error(
                "scene fluid mimetic control-cell byte limit exceeded");
        }
        const auto& control = pressureVolumes.controlVolumes[ownerIndex];
        fluid::Vector3 imageOffset;
        const fluid::Vector3 centroid = unwrapNear(
            sourceCentroid, vector(control.centroidMeters),
            result.lowerMeters, result.upperMeters, imageOffset);
        SceneFluidMimeticHalfFace halfFace;
        halfFace.stableId = halfFaceStableId(
            kind, side, sourceStableId, control.stableId);
        halfFace.kind = kind;
        halfFace.side = side;
        halfFace.sourceIndex = sourceIndex;
        halfFace.sourceStableId = sourceStableId;
        halfFace.traceStableId =
            kind == SceneFluidMimeticHalfFaceKind::MaterialWall
            ? 0 : sourceStableId;
        halfFace.controlVolumeIndex = ownerIndex;
        halfFace.otherControlVolumeIndex = otherIndex.value_or(
            invalidSceneFluidMimeticControlVolumeIndex);
        halfFace.areaSquareMeters = area;
        halfFace.centroidMeters = centroid;
        halfFace.periodicImageOffsetMeters = imageOffset;
        halfFace.outwardUnitNormal = outwardNormal;
        working[ownerIndex].halfFaces.push_back(halfFace);
        ++appendedHalfFaceCount;
        if (!otherIndex
            && kind != SceneFluidMimeticHalfFaceKind::MaterialWall) {
            ++working[ownerIndex].missingOtherControlCount;
        }
    };

    for (const auto& face : faceLinks.faces) {
        if (resolvedFaceStatus(face.status)) {
            continue;
        }
        ++result.unresolvedCartesianFaceCount;
        const auto markCell = [&](const std::size_t cellIndex) {
            const auto& cell = pressureVolumes.cells[cellIndex];
            for (std::size_t offset = 0;
                 offset < cell.controlVolumeCount; ++offset) {
                ++working[cell.firstControlVolume + offset]
                    .unresolvedCartesianIncidentCount;
            }
        };
        markCell(face.minusCellIndex);
        if (face.plusCellIndex != face.minusCellIndex) {
            markCell(face.plusCellIndex);
        }
    }

    for (const auto& link : faceLinks.links) {
        if (link.geometryKind
            != SceneFluidPressureLinkGeometryKind::CartesianFace) {
            continue;
        }
        append(link.minusControlVolumeIndex,
               link.plusControlVolumeIndex,
               SceneFluidMimeticHalfFaceKind::CartesianTrace,
               SceneFluidMimeticHalfFaceSide::MinusOrNegative,
               link.linkIndex, link.stableId, link.areaSquareMeters,
               link.faceCentroidMeters, link.unitNormalMinusToPlus);
        append(link.plusControlVolumeIndex,
               link.minusControlVolumeIndex,
               SceneFluidMimeticHalfFaceKind::CartesianTrace,
               SceneFluidMimeticHalfFaceSide::PlusOrPositive,
               link.linkIndex, link.stableId, link.areaSquareMeters,
               link.faceCentroidMeters,
               negate(link.unitNormalMinusToPlus));
    }

    const auto materialKinematics = sampleSceneFluidQuadratureKinematics(
        surface, state, epoch.quadrature);
    if (materialKinematics.size() != epoch.quadrature.points.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic material kinematics are incomplete");
    }
    for (std::size_t pointIndex = 0;
         pointIndex < epoch.quadrature.points.size(); ++pointIndex) {
        const auto& point = epoch.quadrature.points[pointIndex];
        const auto& kinematics = materialKinematics[pointIndex];
        if (kinematics.stableId != point.stableId) {
            throw std::invalid_argument(
                "scene fluid mimetic material kinematics are out of order");
        }
        const auto negative = controlIndex(
            point.negativeSideCellIndex, point.negativeSideRegionId);
        const auto positive = controlIndex(
            point.positiveSideCellIndex, point.positiveSideRegionId);
        const fluid::Vector3 normal = materialNormal(
            surface, state, point.triangleId);
        const fluid::Vector3 centroid = vector(kinematics.positionMeters);
        if (negative) {
            append(*negative, positive,
                   SceneFluidMimeticHalfFaceKind::MaterialWall,
                   SceneFluidMimeticHalfFaceSide::MinusOrNegative,
                   pointIndex, point.stableId, point.areaSquareMeters,
                   centroid, normal);
        } else {
            ++result.omittedZeroVolumeMaterialSideCount;
        }
        if (positive) {
            append(*positive, negative,
                   SceneFluidMimeticHalfFaceKind::MaterialWall,
                   SceneFluidMimeticHalfFaceSide::PlusOrPositive,
                   pointIndex, point.stableId, point.areaSquareMeters,
                   centroid, negate(normal));
        } else {
            ++result.omittedZeroVolumeMaterialSideCount;
        }
    }

    for (std::size_t patchIndex = 0;
         patchIndex < openingPatches.patches.size(); ++patchIndex) {
        const auto& patch = openingPatches.patches[patchIndex];
        if (patch.ownerKind != SceneFluidOpeningPatchOwnerKind::Cell) {
            continue;
        }
        const auto negative = controlIndex(
            patch.cellIndex, patch.negativeSideRegionId);
        const auto positive = controlIndex(
            patch.cellIndex, patch.positiveSideRegionId);
        const fluid::Vector3 normal = vector(
            patch.unitNormalNegativeToPositive);
        const fluid::Vector3 centroid = vector(patch.centroidMeters);
        if (negative) {
            append(*negative, positive,
                   SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace,
                   SceneFluidMimeticHalfFaceSide::MinusOrNegative,
                   patchIndex, patch.stableId, patch.areaSquareMeters,
                   centroid, normal);
        } else {
            ++result.missingOpeningControlSideCount;
        }
        if (positive) {
            append(*positive, negative,
                   SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace,
                   SceneFluidMimeticHalfFaceSide::PlusOrPositive,
                   patchIndex, patch.stableId, patch.areaSquareMeters,
                   centroid, negate(normal));
        } else {
            ++result.missingOpeningControlSideCount;
        }
    }

    result.controlCells.reserve(controlCount);
    result.halfFaces.reserve(appendedHalfFaceCount);
    for (std::size_t controlIndexValue = 0;
         controlIndexValue < controlCount; ++controlIndexValue) {
        auto& source = working[controlIndexValue];
        std::ranges::sort(
            source.halfFaces,
            [](const auto& first, const auto& second) {
                return std::tie(first.kind, first.sourceStableId,
                                first.side, first.stableId)
                    < std::tie(second.kind, second.sourceStableId,
                               second.side, second.stableId);
            });
        const auto& control =
            pressureVolumes.controlVolumes[controlIndexValue];
        SceneFluidMimeticControlCell cell;
        cell.controlCellIndex = result.controlCells.size();
        cell.controlVolumeIndex = controlIndexValue;
        cell.stableId = control.stableId;
        cell.cellIndex = control.cellIndex;
        cell.regionId = control.regionId;
        cell.componentIndex = control.componentIndex;
        cell.firstHalfFace = result.halfFaces.size();
        cell.halfFaceCount = source.halfFaces.size();
        result.maximumHalfFaceCountPerControl = std::max(
            result.maximumHalfFaceCountPerControl, cell.halfFaceCount);
        cell.unresolvedCartesianIncidentCount =
            source.unresolvedCartesianIncidentCount;
        cell.missingOtherControlCount = source.missingOtherControlCount;
        cell.volumeCubicMeters = control.volumeCubicMeters;
        cell.centroidMeters = vector(control.centroidMeters);
        double summedArea = 0.0;
        for (auto& face : source.halfFaces) {
            face.halfFaceIndex = result.halfFaces.size();
            result.halfFaces.push_back(face);
            summedArea += face.areaSquareMeters;
            if (face.kind
                == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
                ++cell.cartesianHalfFaceCount;
                ++result.cartesianHalfFaceCount;
            } else if (face.kind
                       == SceneFluidMimeticHalfFaceKind::MaterialWall) {
                ++cell.materialWallHalfFaceCount;
                ++result.materialWallHalfFaceCount;
            } else {
                ++cell.openingHalfFaceCount;
                ++result.openingHalfFaceCount;
            }
            cell.areaClosureVectorSquareMeters.x +=
                face.areaSquareMeters * face.outwardUnitNormal.x;
            cell.areaClosureVectorSquareMeters.y +=
                face.areaSquareMeters * face.outwardUnitNormal.y;
            cell.areaClosureVectorSquareMeters.z +=
                face.areaSquareMeters * face.outwardUnitNormal.z;
            const fluid::Vector3 offset = subtract(
                face.centroidMeters, cell.centroidMeters);
            for (std::size_t row = 0; row < 3; ++row) {
                for (std::size_t column = 0; column < 3; ++column) {
                    cell.normalMomentMatrixCubicMeters[row * 3 + column] +=
                        component(face.outwardUnitNormal, row)
                        * face.areaSquareMeters
                        * component(offset, column);
                }
            }
        }
        cell.maximumAreaClosureErrorSquareMeters = std::max({
            std::abs(cell.areaClosureVectorSquareMeters.x),
            std::abs(cell.areaClosureVectorSquareMeters.y),
            std::abs(cell.areaClosureVectorSquareMeters.z),
        });
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = 0; column < 3; ++column) {
                const double expected = row == column
                    ? cell.volumeCubicMeters : 0.0;
                cell.maximumDivergenceTheoremErrorCubicMeters = std::max(
                    cell.maximumDivergenceTheoremErrorCubicMeters,
                    std::abs(cell.normalMomentMatrixCubicMeters[
                                 row * 3 + column] - expected));
            }
        }
        const double areaTolerance =
            settings.absoluteAreaClosureToleranceSquareMeters
            + settings.relativeGeometryTolerance * summedArea;
        const double volumeTolerance =
            settings.absoluteDivergenceTheoremToleranceCubicMeters
            + settings.relativeGeometryTolerance * cell.volumeCubicMeters;
        cell.topologyComplete =
            cell.unresolvedCartesianIncidentCount == 0
            && cell.missingOtherControlCount == 0;
        cell.areaVectorClosed =
            cell.maximumAreaClosureErrorSquareMeters <= areaTolerance;
        cell.divergenceTheoremClosed =
            cell.maximumDivergenceTheoremErrorCubicMeters
            <= volumeTolerance;
        cell.readyForMimeticOperator = cell.topologyComplete
            && cell.halfFaceCount >= 4
            && cell.areaVectorClosed
            && cell.divergenceTheoremClosed;
        if (cell.readyForMimeticOperator) {
            ++result.readyControlCellCount;
        }
        if (!cell.topologyComplete) {
            ++result.incompleteTopologyControlCellCount;
        }
        if (cell.halfFaceCount < 4 || !cell.areaVectorClosed
            || !cell.divergenceTheoremClosed) {
            ++result.nonclosingControlCellCount;
        }
        result.maximumAreaClosureErrorSquareMeters = std::max(
            result.maximumAreaClosureErrorSquareMeters,
            cell.maximumAreaClosureErrorSquareMeters);
        result.maximumDivergenceTheoremErrorCubicMeters = std::max(
            result.maximumDivergenceTheoremErrorCubicMeters,
            cell.maximumDivergenceTheoremErrorCubicMeters);
        result.controlCells.push_back(cell);
    }

    result.ownedStorageBytes = checkedStorageBytes(
        result.controlCells.size(), result.halfFaces.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic control-cell byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    validateSceneFluidMimeticControlCellIntegrity(result);
    return result;
}

void validateSceneFluidMimeticControlCellIntegrity(
    const SceneFluidMimeticControlCellSet& set) {
    validateSettings(set.settings);
    if (set.version != sceneFluidMimeticControlCellVersion
        || set.fingerprint == 0
        || set.surfaceDefinitionFingerprint == 0
        || set.surfaceStateFingerprint == 0
        || set.gridEpochFingerprint == 0
        || set.openingPatchFingerprint == 0
        || set.pressureControlVolumeFingerprint == 0
        || set.pressureFaceLinkFingerprint == 0
        || set.structureDefinitionFingerprint == 0
        || !std::isfinite(set.simulationTimeSeconds)
        || !finite(set.lowerMeters) || !finite(set.upperMeters)
        || set.cellCounts.x == 0 || set.cellCounts.y == 0
        || set.cellCounts.z == 0
        || set.ownedStorageBytes != checkedStorageBytes(
            set.controlCells.size(), set.halfFaces.size())
        || set.readyControlCellCount > set.controlCells.size()
        || set.incompleteTopologyControlCellCount > set.controlCells.size()
        || set.nonclosingControlCellCount > set.controlCells.size()
        || set.cartesianHalfFaceCount > set.halfFaces.size()
        || set.materialWallHalfFaceCount > set.halfFaces.size()
        || set.openingHalfFaceCount > set.halfFaces.size()
        || set.maximumHalfFaceCountPerControl > set.halfFaces.size()
        || set.cartesianHalfFaceCount + set.materialWallHalfFaceCount
            + set.openingHalfFaceCount != set.halfFaces.size()
        || !std::isfinite(set.maximumAreaClosureErrorSquareMeters)
        || set.maximumAreaClosureErrorSquareMeters < 0.0
        || !std::isfinite(
            set.maximumDivergenceTheoremErrorCubicMeters)
        || set.maximumDivergenceTheoremErrorCubicMeters < 0.0) {
        throw std::invalid_argument(
            "invalid scene fluid mimetic control-cell set");
    }
    std::size_t nextHalfFace = 0;
    std::size_t readyCount = 0;
    std::size_t incompleteCount = 0;
    std::size_t nonclosingCount = 0;
    std::size_t cartesianCount = 0;
    std::size_t materialCount = 0;
    std::size_t openingCount = 0;
    std::size_t maximumHalfFaceCount = 0;
    std::vector<std::uint64_t> stableIds;
    stableIds.reserve(set.halfFaces.size());
    for (std::size_t cellIndex = 0;
         cellIndex < set.controlCells.size(); ++cellIndex) {
        const auto& cell = set.controlCells[cellIndex];
        if (cell.controlCellIndex != cellIndex
            || cell.controlVolumeIndex != cellIndex
            || cell.stableId == 0 || cell.regionId == invalidStableId
            || cell.firstHalfFace != nextHalfFace
            || cell.firstHalfFace > set.halfFaces.size()
            || cell.halfFaceCount
                > set.halfFaces.size() - cell.firstHalfFace
            || cell.cartesianHalfFaceCount
                + cell.materialWallHalfFaceCount
                + cell.openingHalfFaceCount != cell.halfFaceCount
            || !std::isfinite(cell.volumeCubicMeters)
            || cell.volumeCubicMeters <= 0.0
            || !finite(cell.centroidMeters)
            || !finite(cell.areaClosureVectorSquareMeters)
            || !std::ranges::all_of(
                cell.normalMomentMatrixCubicMeters,
                [](const double value) { return std::isfinite(value); })
            || !std::isfinite(
                cell.maximumAreaClosureErrorSquareMeters)
            || cell.maximumAreaClosureErrorSquareMeters < 0.0
            || !std::isfinite(
                cell.maximumDivergenceTheoremErrorCubicMeters)
            || cell.maximumDivergenceTheoremErrorCubicMeters < 0.0
            || cell.topologyComplete
                != (cell.unresolvedCartesianIncidentCount == 0
                    && cell.missingOtherControlCount == 0)
            || cell.readyForMimeticOperator
                != (cell.topologyComplete && cell.halfFaceCount >= 4
                    && cell.areaVectorClosed
                    && cell.divergenceTheoremClosed)) {
            throw std::invalid_argument(
                "invalid scene fluid mimetic control-cell record");
        }
        nextHalfFace += cell.halfFaceCount;
        readyCount += cell.readyForMimeticOperator ? 1 : 0;
        incompleteCount += cell.topologyComplete ? 0 : 1;
        nonclosingCount += cell.halfFaceCount < 4
                || !cell.areaVectorClosed
                || !cell.divergenceTheoremClosed
            ? 1 : 0;
        cartesianCount += cell.cartesianHalfFaceCount;
        materialCount += cell.materialWallHalfFaceCount;
        openingCount += cell.openingHalfFaceCount;
        maximumHalfFaceCount = std::max(
            maximumHalfFaceCount, cell.halfFaceCount);
    }
    if (nextHalfFace != set.halfFaces.size()
        || readyCount != set.readyControlCellCount
        || incompleteCount != set.incompleteTopologyControlCellCount
        || nonclosingCount != set.nonclosingControlCellCount
        || cartesianCount != set.cartesianHalfFaceCount
        || materialCount != set.materialWallHalfFaceCount
        || openingCount != set.openingHalfFaceCount
        || maximumHalfFaceCount
            != set.maximumHalfFaceCountPerControl) {
        throw std::invalid_argument(
            "scene fluid mimetic aggregate counts are inconsistent");
    }
    for (std::size_t faceIndex = 0;
         faceIndex < set.halfFaces.size(); ++faceIndex) {
        const auto& face = set.halfFaces[faceIndex];
        const double normalLength = std::sqrt(dot(
            face.outwardUnitNormal, face.outwardUnitNormal));
        if (face.halfFaceIndex != faceIndex || face.stableId == 0
            || face.sourceStableId == 0
            || face.controlVolumeIndex >= set.controlCells.size()
            || (face.otherControlVolumeIndex
                    != invalidSceneFluidMimeticControlVolumeIndex
                && face.otherControlVolumeIndex >= set.controlCells.size())
            || !std::isfinite(face.areaSquareMeters)
            || face.areaSquareMeters <= 0.0
            || !finite(face.centroidMeters)
            || !finite(face.periodicImageOffsetMeters)
            || !finite(face.outwardUnitNormal)
            || !std::isfinite(normalLength)
            || std::abs(normalLength - 1.0)
                > set.settings.unitNormalTolerance
            || (face.kind == SceneFluidMimeticHalfFaceKind::MaterialWall
                ? face.traceStableId != 0
                : face.traceStableId == 0)) {
            throw std::invalid_argument(
                "invalid scene fluid mimetic half-face record");
        }
        stableIds.push_back(face.stableId);
    }
    std::ranges::sort(stableIds);
    if (std::ranges::adjacent_find(stableIds) != stableIds.end()) {
        throw std::invalid_argument(
            "scene fluid mimetic half-face stable IDs collide");
    }
    if (productFingerprint(set) != set.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic control-cell fingerprint mismatch");
    }
}

void validateSceneFluidMimeticControlCells(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks) {
    validateSourceBindings(
        surface, state, grid, epoch, caps, openingQuadrature,
        openingPatches, pressureVolumes, faceLinks);
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    if (controlCells.surfaceDefinitionFingerprint != surface.fingerprint
        || controlCells.surfaceStateFingerprint != state.fingerprint
        || controlCells.gridEpochFingerprint != epoch.fingerprint
        || controlCells.openingPatchFingerprint != openingPatches.fingerprint
        || controlCells.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || controlCells.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || controlCells.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || controlCells.acceptedStepCount != state.acceptedStepCount
        || controlCells.simulationTimeSeconds != state.simulationTimeSeconds
        || controlCells.cellCounts != grid.cellCounts()
        || controlCells.lowerMeters != grid.lowerMeters()
        || controlCells.upperMeters != grid.upperMeters()) {
        throw std::invalid_argument(
            "scene fluid mimetic control-cell binding is invalid");
    }
}

} // namespace simwing::fsi
