#include "scene_fluid_region_momentum.h"

#include "scene_fluid_mimetic_control_cell.h"
#include "scene_fluid_mimetic_pressure_flow.h"

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

double dot(const fluid::Vector3& first, const fluid::Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

double dot(const Vec3& first, const fluid::Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

using Matrix3 = std::array<std::array<double, 3>, 3>;

fluid::Vector3 closestNormalEquationVelocity(
    Matrix3 matrix,
    const fluid::Vector3 rightHandSide,
    const fluid::Vector3 fallback) {
    Matrix3 eigenvectors{{
        {{1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 1.0}},
    }};
    constexpr std::array<std::array<std::size_t, 2>, 3> pairs{{
        {{0, 1}}, {{0, 2}}, {{1, 2}},
    }};
    for (std::size_t sweep = 0; sweep < 16; ++sweep) {
        for (const auto pair : pairs) {
            const std::size_t p = pair[0];
            const std::size_t q = pair[1];
            const double offDiagonal = matrix[p][q];
            if (offDiagonal == 0.0) {
                continue;
            }
            const double tau = (matrix[q][q] - matrix[p][p])
                / (2.0 * offDiagonal);
            const double tangent = std::copysign(1.0, tau)
                / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
            const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
            const double sine = tangent * cosine;
            matrix[p][p] -= tangent * offDiagonal;
            matrix[q][q] += tangent * offDiagonal;
            matrix[p][q] = 0.0;
            matrix[q][p] = 0.0;
            for (std::size_t row = 0; row < 3; ++row) {
                if (row == p || row == q) {
                    continue;
                }
                const double rowP = matrix[row][p];
                const double rowQ = matrix[row][q];
                matrix[row][p] = cosine * rowP - sine * rowQ;
                matrix[p][row] = matrix[row][p];
                matrix[row][q] = sine * rowP + cosine * rowQ;
                matrix[q][row] = matrix[row][q];
            }
            for (std::size_t row = 0; row < 3; ++row) {
                const double rowP = eigenvectors[row][p];
                const double rowQ = eigenvectors[row][q];
                eigenvectors[row][p] = cosine * rowP - sine * rowQ;
                eigenvectors[row][q] = sine * rowP + cosine * rowQ;
            }
        }
    }

    const std::array<double, 3> fallbackArray{
        fallback.x, fallback.y, fallback.z};
    const std::array<double, 3> rightArray{
        rightHandSide.x, rightHandSide.y, rightHandSide.z};
    // The Jacobi rotation overwrote the original matrix. In the eigenbasis,
    // M*fallback is simply lambda times the projected fallback, so solve the
    // closest-to-fallback correction directly there.
    std::array<double, 3> result = fallbackArray;
    const double maximumEigenvalue = std::max({
        0.0, matrix[0][0], matrix[1][1], matrix[2][2]});
    const double threshold = maximumEigenvalue * 1.0e-12;
    for (std::size_t eigen = 0; eigen < 3; ++eigen) {
        const double eigenvalue = matrix[eigen][eigen];
        if (!(eigenvalue > threshold)) {
            continue;
        }
        double projectedRight = 0.0;
        double projectedFallback = 0.0;
        for (std::size_t row = 0; row < 3; ++row) {
            projectedRight += eigenvectors[row][eigen] * rightArray[row];
            projectedFallback +=
                eigenvectors[row][eigen] * fallbackArray[row];
        }
        const double correction =
            projectedRight / eigenvalue - projectedFallback;
        for (std::size_t row = 0; row < 3; ++row) {
            result[row] += eigenvectors[row][eigen] * correction;
        }
    }
    return {result[0], result[1], result[2]};
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
        diagnostics.embeddedOpeningLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.normalEquationControlCount));
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
    const double normalAlignment = dot(
        patch.unitNormalNegativeToPositive,
        source.unitNormalMinusToPlus);
    if (patch.openingId != source.openingId
        || patch.areaSquareMeters != source.areaSquareMeters
        || (!forwardRegions && !reverseRegions)
        || std::abs(std::abs(normalAlignment) - 1.0) > 1.0e-10) {
        throw std::invalid_argument(
            "scene fluid region momentum opening patch is foreign");
    }
    const double orientation = forwardRegions ? 1.0 : -1.0;
    if (orientation * normalAlignment < 1.0 - 1.0e-10) {
        throw std::invalid_argument(
            "scene fluid region momentum opening orientation is inconsistent");
    }
    return result
        + orientation * patch.surfaceSweepRateCubicMetersPerSecond;
}

struct AdditionalRegionMomentumNormalSample {
    std::size_t minusControlVolumeIndex = 0;
    std::size_t plusControlVolumeIndex = 0;
    double areaSquareMeters = 0.0;
    fluid::Vector3 unitNormalMinusToPlus;
    double absoluteNormalVelocityMetersPerSecond = 0.0;
};

SceneFluidRegionMomentumState reconstructFromAbsoluteLinkVelocities(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const std::vector<double>& absoluteVelocities,
    const std::vector<AdditionalRegionMomentumNormalSample>&
        additionalNormalSamples,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond,
    const std::uint64_t pressureProjectionFingerprint,
    const std::uint64_t acceptedStepCount,
    const double simulationTimeSeconds,
    const double densityKgPerCubicMeter,
    const SceneFluidRegionMomentumLimits& limits) {
    const std::size_t storageBytes = storageBytesForControlVolumes(
        pressureVolumes.controlVolumes.size());
    if (pressureVolumes.controlVolumes.size() > limits.maximumControlVolumes
        || storageBytes > limits.maximumMomentumBytes) {
        throw std::length_error(
            "scene fluid region momentum exceeds its limits");
    }
    if (absoluteVelocities.size() != faceLinks.links.size()) {
        throw std::invalid_argument(
            "scene fluid region momentum link velocity count is invalid");
    }

    struct Accumulator {
        std::array<double, 3> area{};
        std::array<double, 3> areaVelocity{};
        std::array<std::size_t, 3> linkCount{};
        Matrix3 normalMatrix{};
        fluid::Vector3 normalRightHandSide;
        std::size_t embeddedLinkCount = 0;
    };
    std::vector<Accumulator> accumulators(
        pressureVolumes.controlVolumes.size());
    std::size_t openingLinkCount = 0;
    std::size_t embeddedOpeningLinkCount = 0;
    const auto accumulateNormalSample = [&accumulators](
        const std::size_t minusControlVolumeIndex,
        const std::size_t plusControlVolumeIndex,
        const double areaSquareMeters,
        const fluid::Vector3& normal,
        const double velocity,
        const bool embedded) {
        if (minusControlVolumeIndex >= accumulators.size()
            || plusControlVolumeIndex >= accumulators.size()
            || minusControlVolumeIndex == plusControlVolumeIndex
            || !(areaSquareMeters > 0.0)
            || !finite(normal) || !std::isfinite(velocity)) {
            throw std::invalid_argument(
                "scene fluid region momentum normal sample is invalid");
        }
        const std::array<double, 3> components{
            normal.x, normal.y, normal.z};
        for (const std::size_t controlVolume : {
                 minusControlVolumeIndex,
                 plusControlVolumeIndex}) {
            auto& accumulator = accumulators[controlVolume];
            if (embedded) {
                ++accumulator.embeddedLinkCount;
            }
            for (std::size_t row = 0; row < 3; ++row) {
                component(accumulator.normalRightHandSide, row) +=
                    areaSquareMeters * velocity * components[row];
                for (std::size_t column = 0; column < 3; ++column) {
                    accumulator.normalMatrix[row][column] +=
                        areaSquareMeters
                        * components[row] * components[column];
                }
            }
        }
    };
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        const bool cartesian = source.geometryKind
            == SceneFluidPressureLinkGeometryKind::CartesianFace;
        const bool embedded = source.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening;
        const double velocity = absoluteVelocities[index];
        if (source.linkIndex != index
            || (cartesian && source.faceIndex >= faceLinks.faces.size())
            || (embedded
                && (source.faceIndex != invalidSceneFluidPressureFaceIndex
                    || source.kind
                        != SceneFluidPressureFaceLinkKind::AuthoredOpening))
            || (!cartesian && !embedded)
            || source.minusControlVolumeIndex >= accumulators.size()
            || source.plusControlVolumeIndex >= accumulators.size()
            || source.minusControlVolumeIndex
                == source.plusControlVolumeIndex
            || !(source.areaSquareMeters > 0.0)
            || !std::isfinite(velocity)) {
            throw std::invalid_argument(
                "scene fluid region momentum link binding is invalid");
        }
        for (const std::size_t controlVolume : {
                 source.minusControlVolumeIndex,
                 source.plusControlVolumeIndex}) {
            auto& accumulator = accumulators[controlVolume];
            if (cartesian) {
                const std::size_t axis = axisIndex(
                    faceLinks.faces[source.faceIndex].axis);
                accumulator.area[axis] += source.areaSquareMeters;
                accumulator.areaVelocity[axis] +=
                    source.areaSquareMeters * velocity;
                ++accumulator.linkCount[axis];
            }
        }
        accumulateNormalSample(
            source.minusControlVolumeIndex,
            source.plusControlVolumeIndex,
            source.areaSquareMeters, source.unitNormalMinusToPlus,
            velocity, embedded);
        if (source.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++openingLinkCount;
            if (embedded) {
                ++embeddedOpeningLinkCount;
            }
        }
    }
    for (const auto& sample : additionalNormalSamples) {
        accumulateNormalSample(
            sample.minusControlVolumeIndex,
            sample.plusControlVolumeIndex,
            sample.areaSquareMeters,
            sample.unitNormalMinusToPlus,
            sample.absoluteNormalVelocityMetersPerSecond, true);
        ++openingLinkCount;
        ++embeddedOpeningLinkCount;
    }

    SceneFluidRegionMomentumState result;
    result.pressureProjectionFingerprint = pressureProjectionFingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.openingPatchFingerprint = openingPatches.fingerprint;
    result.fallbackVelocityFingerprint =
        sceneFluidOpeningFluxVelocityFingerprint(
            grid, fallbackVelocityMetersPerSecond);
    result.acceptedStepCount = acceptedStepCount;
    result.simulationTimeSeconds = simulationTimeSeconds;
    result.densityKgPerCubicMeter = densityKgPerCubicMeter;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.ownedStorageBytes = storageBytes;
    result.controlVolumes.reserve(pressureVolumes.controlVolumes.size());
    auto& diagnostics = result.diagnostics;
    diagnostics.controlVolumeCount = pressureVolumes.controlVolumes.size();
    diagnostics.linkCount =
        faceLinks.links.size() + additionalNormalSamples.size();
    diagnostics.openingLinkCount = openingLinkCount;
    diagnostics.embeddedOpeningLinkCount = embeddedOpeningLinkCount;

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
        control.velocityMetersPerSecond = fallback;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (accumulator.linkCount[axis] != 0) {
                if (!(accumulator.area[axis] > 0.0)) {
                    throw std::invalid_argument(
                        "scene fluid region momentum sampled area is invalid");
                }
                if (accumulator.embeddedLinkCount == 0) {
                    component(control.velocityMetersPerSecond, axis) =
                        accumulator.areaVelocity[axis]
                        / accumulator.area[axis];
                }
                ++diagnostics.sampledComponentCount;
            } else {
                if (accumulator.area[axis] != 0.0) {
                    throw std::invalid_argument(
                        "scene fluid region momentum fallback area is invalid");
                }
                ++diagnostics.fallbackComponentCount;
            }
        }
        if (accumulator.embeddedLinkCount != 0) {
            control.velocityMetersPerSecond = closestNormalEquationVelocity(
                accumulator.normalMatrix,
                accumulator.normalRightHandSide,
                fallback);
            ++diagnostics.normalEquationControlCount;
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double velocity = component(
                control.velocityMetersPerSecond, axis);
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
        double reconstructed = 0.0;
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::CartesianFace) {
            const std::size_t axis = axisIndex(
                faceLinks.faces[link.faceIndex].axis);
            reconstructed = 0.5
                * (component(
                       result.controlVolumes[link.minusControlVolumeIndex]
                           .velocityMetersPerSecond,
                       axis)
                   + component(
                       result.controlVolumes[link.plusControlVolumeIndex]
                           .velocityMetersPerSecond,
                       axis));
        } else {
            reconstructed = 0.5
                * (dot(
                       result.controlVolumes[link.minusControlVolumeIndex]
                           .velocityMetersPerSecond,
                       link.unitNormalMinusToPlus)
                   + dot(
                       result.controlVolumes[link.plusControlVolumeIndex]
                           .velocityMetersPerSecond,
                       link.unitNormalMinusToPlus));
        }
        diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond =
            std::max(
                diagnostics
                    .maximumLinkNormalVelocityResidualMetersPerSecond,
                std::abs(absoluteVelocities[index] - reconstructed));
    }
    for (const auto& sample : additionalNormalSamples) {
        const double reconstructed = 0.5
            * (dot(
                   result.controlVolumes[sample.minusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   sample.unitNormalMinusToPlus)
               + dot(
                   result.controlVolumes[sample.plusControlVolumeIndex]
                       .velocityMetersPerSecond,
                   sample.unitNormalMinusToPlus));
        diagnostics.maximumLinkNormalVelocityResidualMetersPerSecond =
            std::max(
                diagnostics
                    .maximumLinkNormalVelocityResidualMetersPerSecond,
                std::abs(
                    sample.absoluteNormalVelocityMetersPerSecond
                    - reconstructed));
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
    return result;
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
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
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
    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*> patchById;
    for (const auto& patch : openingPatches.patches) {
        if (patch.stableId == 0
            || !patchById.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene fluid region momentum opening identity is invalid");
        }
    }
    std::vector<double> absoluteVelocities(faceLinks.links.size(), 0.0);
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        const auto& projected = projection.links[index];
        if (projected.linkIndex != index
            || projected.stableId != source.stableId
            || projected.faceIndex != source.faceIndex
            || projected.kind != source.kind
            || projected.minusControlVolumeIndex
                != source.minusControlVolumeIndex
            || projected.plusControlVolumeIndex
                != source.plusControlVolumeIndex
            || projected.openingPatchStableId
                != source.openingPatchStableId) {
            throw std::invalid_argument(
                "scene fluid region momentum projected link is foreign");
        }
        absoluteVelocities[index] = absoluteLinkFlow(
            source, projected, patchById) / source.areaSquareMeters;
    }
    auto result = reconstructFromAbsoluteLinkVelocities(
        grid, pressureVolumes, faceLinks, openingPatches,
        absoluteVelocities, {}, fallbackVelocityMetersPerSecond,
        projection.fingerprint, projection.acceptedStepCount,
        projection.simulationTimeSeconds,
        projection.settings.densityKgPerCubicMeter, limits);
    validateSceneFluidRegionMomentumState(
        result, grid, pressureVolumes, faceLinks, openingPatches,
        projection, fallbackVelocityMetersPerSecond);
    return result;
}

SceneFluidRegionMomentumState reconstructSceneFluidRegionMomentumState(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidMimeticControlCellSet& mimeticControlCells,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond,
    const SceneFluidRegionMomentumLimits& limits) {
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningGridPatchIntegrity(openingPatches);
    validateSceneFluidMimeticControlCellIntegrity(mimeticControlCells);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(correctedFlow);
    if (!fallbackVelocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(fallbackVelocityMetersPerSecond)) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum fallback velocity is invalid");
    }
    validateGridIdentity(
        grid, pressureVolumes.cellCounts,
        pressureVolumes.lowerMeters, pressureVolumes.upperMeters,
        "scene fluid mimetic region momentum pressure grid is foreign");
    validateGridIdentity(
        grid, faceLinks.cellCounts,
        faceLinks.lowerMeters, faceLinks.upperMeters,
        "scene fluid mimetic region momentum face grid is foreign");
    validateGridIdentity(
        grid, openingPatches.cellCounts,
        openingPatches.lowerMeters, openingPatches.upperMeters,
        "scene fluid mimetic region momentum opening grid is foreign");
    if (!correctedFlow.accepted) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum corrected flow was rejected");
    }
    if (mimeticControlCells.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || mimeticControlCells.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || mimeticControlCells.openingPatchFingerprint
            != openingPatches.fingerprint
        || correctedFlow.mimeticControlCellFingerprint
            != mimeticControlCells.fingerprint
        || correctedFlow.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || correctedFlow.openingPatchFingerprint
            != openingPatches.fingerprint
        || correctedFlow.structureDefinitionFingerprint
            != pressureVolumes.structureDefinitionFingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum topology identity is invalid");
    }
    if (correctedFlow.acceptedStepCount
            != pressureVolumes.acceptedStepCount
        || correctedFlow.acceptedStepCount != faceLinks.acceptedStepCount
        || correctedFlow.acceptedStepCount
            != openingPatches.acceptedStepCount
        || correctedFlow.simulationTimeSeconds
            != pressureVolumes.simulationTimeSeconds
        || correctedFlow.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || correctedFlow.simulationTimeSeconds
            != openingPatches.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum epoch identity is invalid");
    }
    if (mimeticControlCells.controlCells.size()
        != pressureVolumes.controlVolumes.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum control count is invalid");
    }
    for (std::size_t index = 0;
         index < mimeticControlCells.controlCells.size(); ++index) {
        const auto& mimetic = mimeticControlCells.controlCells[index];
        const auto& pressure = pressureVolumes.controlVolumes[index];
        if (mimetic.controlCellIndex != index
            || mimetic.controlVolumeIndex != index
            || mimetic.stableId != pressure.stableId
            || mimetic.cellIndex != pressure.cellIndex
            || mimetic.regionId != pressure.regionId
            || mimetic.componentIndex != pressure.componentIndex
            || mimetic.volumeCubicMeters != pressure.volumeCubicMeters) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum control binding is invalid");
        }
    }

    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*> patchById;
    std::map<std::uint64_t, std::size_t> embeddedLinkByPatchId;
    for (const auto& patch : openingPatches.patches) {
        if (patch.stableId == 0
            || !patchById.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum opening identity is invalid");
        }
    }
    for (const auto& link : faceLinks.links) {
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening
            && (link.openingPatchStableId == 0
                || !embeddedLinkByPatchId.emplace(
                        link.openingPatchStableId,
                        link.linkIndex).second)) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum embedded identity is invalid");
        }
    }
    std::map<std::uint64_t, const SceneFluidMimeticHalfFace*>
        negativeOpeningHalfFaceByTraceId;
    for (const auto& halfFace : mimeticControlCells.halfFaces) {
        if (halfFace.kind
                == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace
            && halfFace.side
                == SceneFluidMimeticHalfFaceSide::MinusOrNegative
            && (halfFace.traceStableId == 0
                || !negativeOpeningHalfFaceByTraceId.emplace(
                        halfFace.traceStableId, &halfFace).second)) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum opening half-face identity is invalid");
        }
    }
    std::vector<const SceneFluidMimeticCorrectedTrace*> traceByLink(
        faceLinks.links.size(), nullptr);
    std::vector<AdditionalRegionMomentumNormalSample>
        additionalNormalSamples;
    for (const auto& trace : correctedFlow.traces) {
        std::size_t linkIndex = 0;
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            linkIndex = trace.sourceIndex;
        } else if (trace.kind
                   == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            const auto found = embeddedLinkByPatchId.find(
                trace.sourceStableId);
            if (found == embeddedLinkByPatchId.end()) {
                const auto halfFaceFound =
                    negativeOpeningHalfFaceByTraceId.find(
                        trace.stableId);
                const auto patchFound = patchById.find(
                    trace.sourceStableId);
                if (halfFaceFound
                        == negativeOpeningHalfFaceByTraceId.end()
                    || patchFound == patchById.end()) {
                    throw std::invalid_argument(
                        "scene fluid mimetic region momentum unresolved opening trace is missing");
                }
                const auto& halfFace = *halfFaceFound->second;
                const auto& patch = *patchFound->second;
                if (halfFace.sourceStableId != trace.sourceStableId
                    || halfFace.controlVolumeIndex
                        != trace.minusControlCellIndex
                    || halfFace.otherControlVolumeIndex
                        != trace.plusControlCellIndex
                    || halfFace.areaSquareMeters
                        != patch.areaSquareMeters
                    || mimeticControlCells.controlCells[
                            trace.minusControlCellIndex].regionId
                        != patch.negativeSideRegionId
                    || mimeticControlCells.controlCells[
                            trace.plusControlCellIndex].regionId
                        != patch.positiveSideRegionId
                    || (halfFace.outwardUnitNormal.x
                            * patch.unitNormalNegativeToPositive.x
                        + halfFace.outwardUnitNormal.y
                            * patch.unitNormalNegativeToPositive.y
                        + halfFace.outwardUnitNormal.z
                            * patch.unitNormalNegativeToPositive.z)
                        < 1.0 - 1.0e-10) {
                    throw std::invalid_argument(
                        "scene fluid mimetic region momentum unresolved opening binding is invalid");
                }
                AdditionalRegionMomentumNormalSample sample;
                sample.minusControlVolumeIndex =
                    trace.minusControlCellIndex;
                sample.plusControlVolumeIndex =
                    trace.plusControlCellIndex;
                sample.areaSquareMeters = halfFace.areaSquareMeters;
                sample.unitNormalMinusToPlus =
                    halfFace.outwardUnitNormal;
                sample.absoluteNormalVelocityMetersPerSecond =
                    (trace
                         .correctedRelativeVolumeFlowRateCubicMetersPerSecond
                     + patch.surfaceSweepRateCubicMetersPerSecond)
                    / halfFace.areaSquareMeters;
                additionalNormalSamples.push_back(sample);
                continue;
            }
            linkIndex = found->second;
        } else {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum trace kind is invalid");
        }
        if (linkIndex >= traceByLink.size()
            || traceByLink[linkIndex] != nullptr) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum trace ownership is invalid");
        }
        const auto& link = faceLinks.links[linkIndex];
        if (link.linkIndex != linkIndex
            || (trace.kind
                    == SceneFluidMimeticHalfFaceKind::CartesianTrace
                && (link.geometryKind
                        != SceneFluidPressureLinkGeometryKind::CartesianFace
                    || trace.sourceStableId != link.stableId))
            || (trace.kind
                    == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace
                && (link.geometryKind
                        != SceneFluidPressureLinkGeometryKind::EmbeddedOpening
                    || trace.sourceStableId
                        != link.openingPatchStableId))) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum trace binding is invalid");
        }
        traceByLink[linkIndex] = &trace;
    }

    std::vector<double> absoluteVelocities(faceLinks.links.size(), 0.0);
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        if (traceByLink[index] == nullptr) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum missed a link");
        }
        SceneFluidPressureProjectedLink projected;
        projected.correctedRelativeVolumeFlowRateCubicMetersPerSecond =
            traceByLink[index]
                ->correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        absoluteVelocities[index] = absoluteLinkFlow(
            faceLinks.links[index], projected, patchById)
            / faceLinks.links[index].areaSquareMeters;
    }
    auto result = reconstructFromAbsoluteLinkVelocities(
        grid, pressureVolumes, faceLinks, openingPatches,
        absoluteVelocities, additionalNormalSamples,
        fallbackVelocityMetersPerSecond,
        correctedFlow.fingerprint, correctedFlow.acceptedStepCount,
        correctedFlow.simulationTimeSeconds,
        correctedFlow.densityKgPerCubicMeter, limits);
    validateSceneFluidRegionMomentumState(
        result, grid, pressureVolumes, faceLinks, openingPatches,
        mimeticControlCells, correctedFlow,
        fallbackVelocityMetersPerSecond);
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
        || diagnostics.embeddedOpeningLinkCount
            > diagnostics.openingLinkCount
        || diagnostics.normalEquationControlCount
            > momentum.controlVolumes.size()
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
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidPressureProjectionIntegrity(projection);
    validateGridIdentity(
        grid, momentum.cellCounts, momentum.lowerMeters,
        momentum.upperMeters,
        "scene fluid region momentum grid is foreign");
    std::size_t openingLinkCount = 0;
    std::size_t embeddedOpeningLinkCount = 0;
    std::vector<bool> normalEquationControls(
        pressureVolumes.controlVolumes.size(), false);
    for (const auto& link : faceLinks.links) {
        if (link.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++openingLinkCount;
        }
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening) {
            if (link.minusControlVolumeIndex
                    >= normalEquationControls.size()
                || link.plusControlVolumeIndex
                    >= normalEquationControls.size()) {
                throw std::invalid_argument(
                    "scene fluid region momentum embedded binding is invalid");
            }
            ++embeddedOpeningLinkCount;
            normalEquationControls[link.minusControlVolumeIndex] = true;
            normalEquationControls[link.plusControlVolumeIndex] = true;
        }
    }
    const std::size_t normalEquationControlCount =
        static_cast<std::size_t>(std::count(
            normalEquationControls.begin(),
            normalEquationControls.end(), true));
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
            != pressureVolumes.controlVolumes.size()
        || momentum.diagnostics.linkCount != faceLinks.links.size()
        || momentum.diagnostics.openingLinkCount != openingLinkCount
        || momentum.diagnostics.embeddedOpeningLinkCount
            != embeddedOpeningLinkCount
        || momentum.diagnostics.normalEquationControlCount
            != normalEquationControlCount) {
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

void validateSceneFluidRegionMomentumState(
    const SceneFluidRegionMomentumState& momentum,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidMimeticControlCellSet& mimeticControlCells,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond) {
    validateSceneFluidRegionMomentumStateIntegrity(momentum);
    validateSceneFluidPressureControlVolumeIntegrity(pressureVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidOpeningGridPatchIntegrity(openingPatches);
    validateSceneFluidMimeticControlCellIntegrity(mimeticControlCells);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(correctedFlow);
    validateGridIdentity(
        grid, momentum.cellCounts, momentum.lowerMeters,
        momentum.upperMeters,
        "scene fluid mimetic region momentum grid is foreign");
    if (!fallbackVelocityMetersPerSecond.matches(grid)
        || !fluid::isFinite(fallbackVelocityMetersPerSecond)
        || momentum.fallbackVelocityFingerprint
            != sceneFluidOpeningFluxVelocityFingerprint(
                grid, fallbackVelocityMetersPerSecond)
        || momentum.pressureProjectionFingerprint
            != correctedFlow.fingerprint
        || momentum.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || momentum.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || momentum.openingPatchFingerprint != openingPatches.fingerprint
        || correctedFlow.mimeticControlCellFingerprint
            != mimeticControlCells.fingerprint
        || momentum.acceptedStepCount != correctedFlow.acceptedStepCount
        || momentum.simulationTimeSeconds
            != correctedFlow.simulationTimeSeconds
        || momentum.densityKgPerCubicMeter
            != correctedFlow.densityKgPerCubicMeter
        || momentum.controlVolumes.size()
            != pressureVolumes.controlVolumes.size()
        || momentum.controlVolumes.size()
            != mimeticControlCells.controlCells.size()
        || momentum.diagnostics.linkCount != correctedFlow.traces.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum binding is invalid");
    }
    std::size_t openingLinkCount = 0;
    std::size_t embeddedOpeningLinkCount = 0;
    std::vector<bool> normalEquationControls(
        pressureVolumes.controlVolumes.size(), false);
    std::map<std::uint64_t, std::size_t> embeddedLinkByPatchId;
    for (const auto& link : faceLinks.links) {
        if (link.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++openingLinkCount;
        }
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening) {
            if (link.minusControlVolumeIndex
                    >= normalEquationControls.size()
                || link.plusControlVolumeIndex
                    >= normalEquationControls.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic region momentum embedded binding is invalid");
            }
            ++embeddedOpeningLinkCount;
            if (link.openingPatchStableId == 0
                || !embeddedLinkByPatchId.emplace(
                        link.openingPatchStableId,
                        link.linkIndex).second) {
                throw std::invalid_argument(
                    "scene fluid mimetic region momentum embedded identity changed");
            }
            normalEquationControls[link.minusControlVolumeIndex] = true;
            normalEquationControls[link.plusControlVolumeIndex] = true;
        }
    }
    for (const auto& trace : correctedFlow.traces) {
        if (trace.kind
                == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace
            && !embeddedLinkByPatchId.contains(trace.sourceStableId)) {
            if (trace.minusControlCellIndex
                    >= normalEquationControls.size()
                || trace.plusControlCellIndex
                    >= normalEquationControls.size()) {
                throw std::invalid_argument(
                    "scene fluid mimetic region momentum unresolved opening binding changed");
            }
            ++openingLinkCount;
            ++embeddedOpeningLinkCount;
            normalEquationControls[trace.minusControlCellIndex] = true;
            normalEquationControls[trace.plusControlCellIndex] = true;
        }
    }
    if (momentum.diagnostics.openingLinkCount != openingLinkCount
        || momentum.diagnostics.embeddedOpeningLinkCount
            != embeddedOpeningLinkCount
        || momentum.diagnostics.normalEquationControlCount
            != static_cast<std::size_t>(std::count(
                normalEquationControls.begin(),
                normalEquationControls.end(), true))) {
        throw std::invalid_argument(
            "scene fluid mimetic region momentum diagnostics changed");
    }
    for (std::size_t index = 0;
         index < momentum.controlVolumes.size(); ++index) {
        const auto& control = momentum.controlVolumes[index];
        const auto& pressure = pressureVolumes.controlVolumes[index];
        const auto& mimetic = mimeticControlCells.controlCells[index];
        if (control.controlVolumeIndex != pressure.controlVolumeIndex
            || control.stableId != pressure.stableId
            || control.cellIndex != pressure.cellIndex
            || control.regionIndex != pressure.regionIndex
            || control.regionId != pressure.regionId
            || control.componentIndex != pressure.componentIndex
            || control.volumeCubicMeters != pressure.volumeCubicMeters
            || mimetic.controlCellIndex != index
            || mimetic.controlVolumeIndex != pressure.controlVolumeIndex
            || mimetic.stableId != pressure.stableId) {
            throw std::invalid_argument(
                "scene fluid mimetic region momentum control changed");
        }
    }
}

} // namespace simwing::fsi
