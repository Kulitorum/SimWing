#include "scene_fluid_region_momentum.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>

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

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double component(const fluid::Vector3& value, const std::size_t axis) {
    switch (axis) {
    case 0: return value.x;
    case 1: return value.y;
    case 2: return value.z;
    }
    throw std::invalid_argument(
        "scene fluid region momentum component is invalid");
}

double& component(fluid::Vector3& value, const std::size_t axis) {
    switch (axis) {
    case 0: return value.x;
    case 1: return value.y;
    case 2: return value.z;
    }
    throw std::invalid_argument(
        "scene fluid region momentum component is invalid");
}

std::size_t axisIndex(const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return 0;
    case fluid::GridFaceAxis::Y: return 1;
    case fluid::GridFaceAxis::Z: return 2;
    }
    throw std::invalid_argument(
        "scene fluid region momentum face axis is invalid");
}

double axisComponent(const Vec3& value, const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return value.x;
    case fluid::GridFaceAxis::Y: return value.y;
    case fluid::GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "scene fluid region momentum opening axis is invalid");
}

std::size_t storageBytesForControlVolumes(const std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max()
                    / sizeof(SceneFluidRegionMomentumControlVolume)) {
        throw std::length_error(
            "scene fluid region momentum storage size overflows");
    }
    return count * sizeof(SceneFluidRegionMomentumControlVolume);
}

std::uint64_t momentumFingerprint(
    const SceneFluidRegionMomentumState& momentum) {
    Fingerprint fingerprint;
    fingerprint.integer(momentum.version);
    fingerprint.integer(momentum.pressureProjectionFingerprint);
    fingerprint.integer(momentum.pressureControlVolumeFingerprint);
    fingerprint.integer(momentum.pressureFaceLinkFingerprint);
    fingerprint.integer(momentum.openingPatchFingerprint);
    fingerprint.integer(momentum.fallbackVelocityFingerprint);
    fingerprint.integer(momentum.acceptedStepCount);
    fingerprint.real(momentum.simulationTimeSeconds);
    fingerprint.real(momentum.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(momentum.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(momentum.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(momentum.cellCounts.z));
    for (const double value : {
             momentum.lowerMeters.x,
             momentum.lowerMeters.y,
             momentum.lowerMeters.z,
             momentum.upperMeters.x,
             momentum.upperMeters.y,
             momentum.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        momentum.ownedStorageBytes));
    const auto& diagnostics = momentum.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.linkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.openingLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.sampledComponentCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.fallbackComponentCount));
    fingerprint.real(diagnostics.totalMomentumKilogramMetersPerSecond.x);
    fingerprint.real(diagnostics.totalMomentumKilogramMetersPerSecond.y);
    fingerprint.real(diagnostics.totalMomentumKilogramMetersPerSecond.z);
    fingerprint.real(diagnostics.kineticEnergyJoules);
    fingerprint.real(diagnostics.maximumAbsoluteVelocityMetersPerSecond);
    fingerprint.real(
        diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond);
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.finite ? 1 : 0));
    fingerprint.integer(static_cast<std::uint64_t>(
        momentum.controlVolumes.size()));
    for (const auto& control : momentum.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(control.cellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(control.regionIndex));
        fingerprint.integer(control.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.real(control.volumeCubicMeters);
        for (const double value : {
                 control.velocityMetersPerSecond.x,
                 control.velocityMetersPerSecond.y,
                 control.velocityMetersPerSecond.z,
                 control.momentumKilogramMetersPerSecond.x,
                 control.momentumKilogramMetersPerSecond.y,
                 control.momentumKilogramMetersPerSecond.z}) {
            fingerprint.real(value);
        }
        for (const double area : control.sampledFaceAreaSquareMeters) {
            fingerprint.real(area);
        }
        for (const std::size_t count : control.sampledLinkCounts) {
            fingerprint.integer(static_cast<std::uint64_t>(count));
        }
    }
    return fingerprint.value();
}

void validateGridIdentity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::GridCellCounts counts,
    const fluid::Vector3 lower,
    const fluid::Vector3 upper,
    const char* message) {
    if (counts != grid.cellCounts()
        || lower != grid.lowerMeters()
        || upper != grid.upperMeters()) {
        throw std::invalid_argument(message);
    }
}

fluid::Vector3 cellCenteredVelocity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    const fluid::GridCellCoordinate cell) {
    const auto counts = grid.cellCounts();
    const std::size_t center = grid.cellIndex(cell.i, cell.j, cell.k);
    const std::size_t xPlus = grid.cellIndex(
        (cell.i + 1) % counts.x, cell.j, cell.k);
    const std::size_t yPlus = grid.cellIndex(
        cell.i, (cell.j + 1) % counts.y, cell.k);
    const std::size_t zPlus = grid.cellIndex(
        cell.i, cell.j, (cell.k + 1) % counts.z);
    return {
        0.5 * (velocity.xFaces()[center] + velocity.xFaces()[xPlus]),
        0.5 * (velocity.yFaces()[center] + velocity.yFaces()[yPlus]),
        0.5 * (velocity.zFaces()[center] + velocity.zFaces()[zPlus]),
    };
}

double absoluteLinkFlow(
    const SceneFluidPressureFaceLink& source,
    const SceneFluidPressureProjectedLink& projected,
    const SceneFluidPressureFace& face,
    const std::map<std::uint64_t, const SceneFluidOpeningGridPatch*>&
        patchById) {
    double result = projected
        .correctedRelativeVolumeFlowRateCubicMetersPerSecond;
    if (source.kind != SceneFluidPressureFaceLinkKind::AuthoredOpening) {
        return result;
    }
    const auto found = patchById.find(source.openingPatchStableId);
    if (found == patchById.end()) {
        throw std::invalid_argument(
            "scene fluid region momentum opening patch is missing");
    }
    const auto& patch = *found->second;
    const bool forwardRegions =
        patch.negativeSideRegionId == source.minusRegionId
        && patch.positiveSideRegionId == source.plusRegionId;
    const bool reverseRegions =
        patch.negativeSideRegionId == source.plusRegionId
        && patch.positiveSideRegionId == source.minusRegionId;
    const double normal = axisComponent(
        patch.unitNormalNegativeToPositive, face.axis);
    if (patch.openingId != source.openingId
        || patch.areaSquareMeters != source.areaSquareMeters
        || (!forwardRegions && !reverseRegions)
        || std::abs(std::abs(normal) - 1.0) > 1.0e-10) {
        throw std::invalid_argument(
            "scene fluid region momentum opening patch is foreign");
    }
    const double orientation = forwardRegions ? 1.0 : -1.0;
    if (orientation * normal < 1.0 - 1.0e-10) {
        throw std::invalid_argument(
            "scene fluid region momentum opening orientation is inconsistent");
    }
    return result
        + orientation * patch.surfaceSweepRateCubicMetersPerSecond;
}

} // namespace

SceneFluidRegionMomentumState reconstructSceneFluidRegionMomentumState(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond,
    const SceneFluidRegionMomentumLimits& limits) {
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidPressureProjectionIntegrity(projection);
    if (!fallbackVelocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(fallbackVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid region momentum fallback velocity is invalid");
    }
    validateGridIdentity(
        grid, pressureVolumes.cellCounts,
        pressureVolumes.lowerMeters, pressureVolumes.upperMeters,
        "scene fluid region momentum pressure grid is foreign");
    validateGridIdentity(
        grid, faceLinks.cellCounts,
        faceLinks.lowerMeters, faceLinks.upperMeters,
        "scene fluid region momentum face grid is foreign");
    validateGridIdentity(
        grid, openingPatches.cellCounts,
        openingPatches.lowerMeters, openingPatches.upperMeters,
        "scene fluid region momentum opening grid is foreign");
    const std::uint64_t fallbackFingerprint =
        sceneFluidOpeningFluxVelocityFingerprint(
            grid, fallbackVelocityMetersPerSecond);
    if (faceLinks.version != sceneFluidPressureFaceLinkVersion
        || faceLinks.fingerprint == 0
        || openingPatches.version != sceneFluidOpeningGridPatchVersion
        || openingPatches.fingerprint == 0
        || !projection.diagnostics.accepted
        || projection.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || projection.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || projection.velocityFingerprint != fallbackFingerprint
        || faceLinks.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || faceLinks.openingPatchFingerprint != openingPatches.fingerprint
        || projection.acceptedStepCount != pressureVolumes.acceptedStepCount
        || projection.acceptedStepCount != faceLinks.acceptedStepCount
        || projection.acceptedStepCount != openingPatches.acceptedStepCount
        || projection.simulationTimeSeconds
            != pressureVolumes.simulationTimeSeconds
        || projection.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || projection.simulationTimeSeconds
            != openingPatches.simulationTimeSeconds
        || projection.controlVolumes.size()
            != pressureVolumes.controlVolumes.size()
        || projection.links.size() != faceLinks.links.size()) {
        throw std::invalid_argument(
            "scene fluid region momentum identity is invalid");
    }
    const std::size_t storageBytes = storageBytesForControlVolumes(
        pressureVolumes.controlVolumes.size());
    if (pressureVolumes.controlVolumes.size() > limits.maximumControlVolumes
        || storageBytes > limits.maximumMomentumBytes) {
        throw std::length_error(
            "scene fluid region momentum exceeds its limits");
    }

    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*> patchById;
    for (const auto& patch : openingPatches.patches) {
        if (patch.stableId == 0
            || !patchById.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene fluid region momentum opening identity is invalid");
        }
    }

    struct Accumulator {
        std::array<double, 3> area{};
        std::array<double, 3> areaVelocity{};
        std::array<std::size_t, 3> linkCount{};
    };
    std::vector<Accumulator> accumulators(
        pressureVolumes.controlVolumes.size());
    std::vector<double> absoluteVelocities(faceLinks.links.size(), 0.0);
    std::size_t openingLinkCount = 0;
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        const auto& projected = projection.links[index];
        if (source.linkIndex != index
            || source.faceIndex >= faceLinks.faces.size()
            || source.minusControlVolumeIndex >= accumulators.size()
            || source.plusControlVolumeIndex >= accumulators.size()
            || source.minusControlVolumeIndex
                == source.plusControlVolumeIndex
            || projected.linkIndex != index
            || projected.stableId != source.stableId
            || projected.faceIndex != source.faceIndex
            || projected.kind != source.kind
            || projected.minusControlVolumeIndex
                != source.minusControlVolumeIndex
            || projected.plusControlVolumeIndex
                != source.plusControlVolumeIndex
            || projected.openingPatchStableId
                != source.openingPatchStableId
            || !(source.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region momentum link binding is invalid");
        }
        const auto& face = faceLinks.faces[source.faceIndex];
        const std::size_t axis = axisIndex(face.axis);
        const double velocity = absoluteLinkFlow(
            source, projected, face, patchById)
            / source.areaSquareMeters;
        if (!std::isfinite(velocity)) {
            throw std::overflow_error(
                "scene fluid region momentum link velocity is non-finite");
        }
        absoluteVelocities[index] = velocity;
        for (const std::size_t controlVolume : {
                 source.minusControlVolumeIndex,
                 source.plusControlVolumeIndex}) {
            auto& accumulator = accumulators[controlVolume];
            accumulator.area[axis] += source.areaSquareMeters;
            accumulator.areaVelocity[axis] +=
                source.areaSquareMeters * velocity;
            ++accumulator.linkCount[axis];
        }
        if (source.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++openingLinkCount;
        }
    }

    SceneFluidRegionMomentumState result;
    result.pressureProjectionFingerprint = projection.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.openingPatchFingerprint = openingPatches.fingerprint;
    result.fallbackVelocityFingerprint = fallbackFingerprint;
    result.acceptedStepCount = projection.acceptedStepCount;
    result.simulationTimeSeconds = projection.simulationTimeSeconds;
    result.densityKgPerCubicMeter = projection.settings.densityKgPerCubicMeter;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.ownedStorageBytes = storageBytes;
    result.controlVolumes.reserve(pressureVolumes.controlVolumes.size());
    auto& diagnostics = result.diagnostics;
    diagnostics.controlVolumeCount = pressureVolumes.controlVolumes.size();
    diagnostics.linkCount = faceLinks.links.size();
    diagnostics.openingLinkCount = openingLinkCount;

    for (const auto& source : pressureVolumes.controlVolumes) {
        if (source.controlVolumeIndex >= accumulators.size()
            || source.cellIndex >= pressureVolumes.cells.size()
            || pressureVolumes.cells[source.cellIndex].cellIndex
                != source.cellIndex
            || !(source.volumeCubicMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region momentum control volume is invalid");
        }
        const auto fallback = cellCenteredVelocity(
            grid, fallbackVelocityMetersPerSecond,
            pressureVolumes.cells[source.cellIndex].cell);
        const auto& accumulator = accumulators[source.controlVolumeIndex];
        SceneFluidRegionMomentumControlVolume control;
        control.controlVolumeIndex = source.controlVolumeIndex;
        control.stableId = source.stableId;
        control.cellIndex = source.cellIndex;
        control.regionIndex = source.regionIndex;
        control.regionId = source.regionId;
        control.componentIndex = source.componentIndex;
        control.volumeCubicMeters = source.volumeCubicMeters;
        control.sampledFaceAreaSquareMeters = accumulator.area;
        control.sampledLinkCounts = accumulator.linkCount;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            double velocity = component(fallback, axis);
            if (accumulator.linkCount[axis] != 0) {
                if (!(accumulator.area[axis] > 0.0)) {
                    throw std::invalid_argument(
                        "scene fluid region momentum sampled area is invalid");
                }
                velocity = accumulator.areaVelocity[axis]
                    / accumulator.area[axis];
                ++diagnostics.sampledComponentCount;
            } else {
                if (accumulator.area[axis] != 0.0) {
                    throw std::invalid_argument(
                        "scene fluid region momentum fallback area is invalid");
                }
                ++diagnostics.fallbackComponentCount;
            }
            component(control.velocityMetersPerSecond, axis) = velocity;
            component(control.momentumKilogramMetersPerSecond, axis) =
                result.densityKgPerCubicMeter
                * source.volumeCubicMeters * velocity;
            diagnostics.maximumAbsoluteVelocityMetersPerSecond = std::max(
                diagnostics.maximumAbsoluteVelocityMetersPerSecond,
                std::abs(velocity));
        }
        diagnostics.totalMomentumKilogramMetersPerSecond.x +=
            control.momentumKilogramMetersPerSecond.x;
        diagnostics.totalMomentumKilogramMetersPerSecond.y +=
            control.momentumKilogramMetersPerSecond.y;
        diagnostics.totalMomentumKilogramMetersPerSecond.z +=
            control.momentumKilogramMetersPerSecond.z;
        diagnostics.kineticEnergyJoules +=
            0.5 * result.densityKgPerCubicMeter
            * source.volumeCubicMeters
            * (control.velocityMetersPerSecond.x
                   * control.velocityMetersPerSecond.x
               + control.velocityMetersPerSecond.y
                   * control.velocityMetersPerSecond.y
               + control.velocityMetersPerSecond.z
                   * control.velocityMetersPerSecond.z);
        result.controlVolumes.push_back(control);
    }

    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& link = faceLinks.links[index];
        const auto& face = faceLinks.faces[link.faceIndex];
        const std::size_t axis = axisIndex(face.axis);
        const double reconstructed = 0.5
            * (component(
                   result.controlVolumes[link.minusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   axis)
               + component(
                   result.controlVolumes[link.plusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   axis));
        diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond =
            std::max(
                diagnostics
                    .maximumLinkNormalVelocityResidualMetersPerSecond,
                std::abs(absoluteVelocities[index] - reconstructed));
    }
    diagnostics.finite = finite(
            diagnostics.totalMomentumKilogramMetersPerSecond)
        && std::isfinite(diagnostics.kineticEnergyJoules)
        && diagnostics.kineticEnergyJoules >= 0.0
        && std::isfinite(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene fluid region momentum diagnostics are non-finite");
    }
    result.fingerprint = momentumFingerprint(result);
    validateSceneFluidRegionMomentumState(
        result, grid, pressureVolumes, faceLinks, openingPatches,
        projection, fallbackVelocityMetersPerSecond);
    return result;
}

void validateSceneFluidRegionMomentumStateIntegrity(
    const SceneFluidRegionMomentumState& momentum) {
    const auto& diagnostics = momentum.diagnostics;
    fluid::Vector3 totalMomentum;
    double kineticEnergy = 0.0;
    double maximumVelocity = 0.0;
    std::size_t sampledComponentCount = 0;
    std::size_t fallbackComponentCount = 0;
    bool controlsValid = true;
    for (std::size_t index = 0;
         index < momentum.controlVolumes.size(); ++index) {
        const auto& control = momentum.controlVolumes[index];
        controlsValid = controlsValid
            && control.controlVolumeIndex == index
            && control.stableId != 0
            && control.regionId != invalidStableId
            && control.regionIndex < momentum.controlVolumes.size()
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && finite(control.velocityMetersPerSecond)
            && finite(control.momentumKilogramMetersPerSecond);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double area = control.sampledFaceAreaSquareMeters[axis];
            const std::size_t count = control.sampledLinkCounts[axis];
            controlsValid = controlsValid
                && std::isfinite(area) && area >= 0.0
                && ((count == 0 && area == 0.0)
                    || (count != 0 && area > 0.0));
            if (count == 0) {
                ++fallbackComponentCount;
            } else {
                ++sampledComponentCount;
            }
            const double expectedMomentum = momentum.densityKgPerCubicMeter
                * control.volumeCubicMeters
                * component(control.velocityMetersPerSecond, axis);
            controlsValid = controlsValid
                && component(
                       control.momentumKilogramMetersPerSecond, axis)
                    == expectedMomentum;
            maximumVelocity = std::max(
                maximumVelocity,
                std::abs(component(control.velocityMetersPerSecond, axis)));
        }
        totalMomentum.x += control.momentumKilogramMetersPerSecond.x;
        totalMomentum.y += control.momentumKilogramMetersPerSecond.y;
        totalMomentum.z += control.momentumKilogramMetersPerSecond.z;
        kineticEnergy += 0.5 * momentum.densityKgPerCubicMeter
            * control.volumeCubicMeters
            * (control.velocityMetersPerSecond.x
                   * control.velocityMetersPerSecond.x
               + control.velocityMetersPerSecond.y
                   * control.velocityMetersPerSecond.y
               + control.velocityMetersPerSecond.z
                   * control.velocityMetersPerSecond.z);
    }
    if (momentum.version != sceneFluidRegionMomentumVersion
        || momentum.fingerprint == 0
        || momentum.pressureProjectionFingerprint == 0
        || momentum.pressureControlVolumeFingerprint == 0
        || momentum.pressureFaceLinkFingerprint == 0
        || momentum.openingPatchFingerprint == 0
        || momentum.fallbackVelocityFingerprint == 0
        || !std::isfinite(momentum.simulationTimeSeconds)
        || !std::isfinite(momentum.densityKgPerCubicMeter)
        || !(momentum.densityKgPerCubicMeter > 0.0)
        || momentum.cellCounts.x == 0
        || momentum.cellCounts.y == 0
        || momentum.cellCounts.z == 0
        || !finite(momentum.lowerMeters)
        || !finite(momentum.upperMeters)
        || !(momentum.upperMeters.x > momentum.lowerMeters.x)
        || !(momentum.upperMeters.y > momentum.lowerMeters.y)
        || !(momentum.upperMeters.z > momentum.lowerMeters.z)
        || momentum.ownedStorageBytes
            != storageBytesForControlVolumes(momentum.controlVolumes.size())
        || diagnostics.controlVolumeCount
            != momentum.controlVolumes.size()
        || diagnostics.linkCount == 0
        || diagnostics.openingLinkCount > diagnostics.linkCount
        || diagnostics.sampledComponentCount != sampledComponentCount
        || diagnostics.fallbackComponentCount != fallbackComponentCount
        || sampledComponentCount + fallbackComponentCount
            != 3 * momentum.controlVolumes.size()
        || diagnostics.totalMomentumKilogramMetersPerSecond != totalMomentum
        || diagnostics.kineticEnergyJoules != kineticEnergy
        || diagnostics.maximumAbsoluteVelocityMetersPerSecond
            != maximumVelocity
        || !std::isfinite(
            diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond)
        || diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond < 0.0
        || !diagnostics.finite
        || !controlsValid
        || momentum.fingerprint != momentumFingerprint(momentum)) {
        throw std::invalid_argument(
            "scene fluid region momentum integrity is invalid");
    }
}

void validateSceneFluidRegionMomentumState(
    const SceneFluidRegionMomentumState& momentum,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond) {
    validateSceneFluidRegionMomentumStateBinding(
        momentum, grid, pressureVolumes, faceLinks, openingPatches,
        projection);
    if (!fallbackVelocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(fallbackVelocityMetersPerSecond)
        || momentum.fallbackVelocityFingerprint
            != sceneFluidOpeningFluxVelocityFingerprint(
                grid, fallbackVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid region momentum fallback velocity is foreign");
    }
}

void validateSceneFluidRegionMomentumStateBinding(
    const SceneFluidRegionMomentumState& momentum,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection) {
    validateSceneFluidRegionMomentumStateIntegrity(momentum);
    validateSceneFluidPressureProjectionIntegrity(projection);
    validateGridIdentity(
        grid, momentum.cellCounts, momentum.lowerMeters,
        momentum.upperMeters,
        "scene fluid region momentum grid is foreign");
    if (momentum.pressureProjectionFingerprint != projection.fingerprint
        || momentum.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || momentum.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || momentum.openingPatchFingerprint != openingPatches.fingerprint
        || momentum.acceptedStepCount != projection.acceptedStepCount
        || momentum.simulationTimeSeconds
            != projection.simulationTimeSeconds
        || momentum.densityKgPerCubicMeter
            != projection.settings.densityKgPerCubicMeter
        || momentum.controlVolumes.size()
            != pressureVolumes.controlVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid region momentum binding is invalid");
    }
    for (std::size_t index = 0;
         index < momentum.controlVolumes.size(); ++index) {
        const auto& control = momentum.controlVolumes[index];
        const auto& source = pressureVolumes.controlVolumes[index];
        if (control.controlVolumeIndex != source.controlVolumeIndex
            || control.stableId != source.stableId
            || control.cellIndex != source.cellIndex
            || control.regionIndex != source.regionIndex
            || control.regionId != source.regionId
            || control.componentIndex != source.componentIndex
            || control.volumeCubicMeters != source.volumeCubicMeters) {
            throw std::invalid_argument(
                "scene fluid region momentum control binding is invalid");
        }
    }
}

} // namespace simwing::fsi
