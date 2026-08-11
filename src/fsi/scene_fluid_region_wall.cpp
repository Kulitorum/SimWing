#include "scene_fluid_region_wall.h"

#include "scene_fluid_regional_opening_momentum_wall_exchange.h"
#include "scene_fluid_wall_exchange_kernel.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <ranges>
#include <span>
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
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        integer(static_cast<std::make_unsigned_t<Underlying>>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

fluid::Vector3 add(const fluid::Vector3& first,
                   const fluid::Vector3& second) {
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

fluid::Vector3 scale(const fluid::Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(const fluid::Vector3& first, const fluid::Vector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

double norm(const fluid::Vector3& value) {
    return std::sqrt(dot(value, value));
}

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double tolerance(const double absolute,
                 const double relative,
                 const double reference) {
    return absolute + relative * std::max(1.0, std::abs(reference));
}

double kineticEnergy(
    const std::span<const SceneFluidRegionWallControlVolume> controls,
    const double density) {
    double result = 0.0;
    for (const auto& control : controls) {
        result += 0.5 * density * control.volumeCubicMeters
            * dot(control.velocityMetersPerSecond,
                  control.velocityMetersPerSecond);
    }
    return result;
}

fluid::Vector3 totalMomentum(
    const std::span<const SceneFluidRegionWallControlVolume> controls) {
    fluid::Vector3 result;
    for (const auto& control : controls) {
        result = add(result, control.momentumKilogramMetersPerSecond);
    }
    return result;
}

bool validSettings(const SceneFluidRegionWallSettings& settings) {
    return std::isfinite(settings.timeStepSeconds)
        && settings.timeStepSeconds > 0.0
        && std::isfinite(settings.kinematicViscositySquareMetersPerSecond)
        && settings.kinematicViscositySquareMetersPerSecond >= 0.0
        && std::isfinite(settings.minimumWallDistanceMeters)
        && settings.minimumWallDistanceMeters > 0.0
        && std::isfinite(settings.maximumViscousNumber)
        && settings.maximumViscousNumber > 0.0
        && settings.maximumViscousNumber <= 1.0
        && settings.maximumSubsteps != 0
        && std::isfinite(
            settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && settings.absoluteMomentumToleranceKilogramMetersPerSecond >= 0.0
        && std::isfinite(settings.relativeMomentumTolerance)
        && settings.relativeMomentumTolerance >= 0.0
        && std::isfinite(settings.absoluteEnergyToleranceJoules)
        && settings.absoluteEnergyToleranceJoules >= 0.0
        && std::isfinite(settings.relativeEnergyTolerance)
        && settings.relativeEnergyTolerance >= 0.0;
}

std::size_t checkedStorageBytes(const std::size_t controls,
                                const std::size_t samples) {
    if (controls > std::numeric_limits<std::size_t>::max()
                       / sizeof(SceneFluidRegionWallControlVolume)
        || samples > std::numeric_limits<std::size_t>::max()
                       / sizeof(SceneFluidRegionWallSample)) {
        throw std::length_error(
            "scene fluid region wall storage size overflows");
    }
    const std::size_t controlBytes = controls
        * sizeof(SceneFluidRegionWallControlVolume);
    const std::size_t sampleBytes = samples
        * sizeof(SceneFluidRegionWallSample);
    if (sampleBytes > std::numeric_limits<std::size_t>::max()
                          - controlBytes) {
        throw std::length_error(
            "scene fluid region wall storage size overflows");
    }
    return controlBytes + sampleBytes;
}

std::uint64_t exchangeFingerprint(
    const SceneFluidRegionWallExchange& exchange) {
    Fingerprint fingerprint;
    fingerprint.integer(exchange.version);
    fingerprint.integer(exchange.sourceTransportFingerprint);
    fingerprint.integer(exchange.sourceRegionRebaseFingerprint);
    fingerprint.integer(exchange.currentPressureControlVolumeFingerprint);
    fingerprint.integer(exchange.quadratureFingerprint);
    fingerprint.integer(exchange.surfaceDefinitionFingerprint);
    fingerprint.integer(exchange.surfaceStateFingerprint);
    fingerprint.integer(exchange.structureDefinitionFingerprint);
    fingerprint.integer(exchange.acceptedStepCount);
    fingerprint.real(exchange.simulationTimeSeconds);
    fingerprint.real(exchange.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(exchange.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(exchange.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(exchange.cellCounts.z));
    for (const double value : {
             exchange.lowerMeters.x, exchange.lowerMeters.y,
             exchange.lowerMeters.z, exchange.upperMeters.x,
             exchange.upperMeters.y, exchange.upperMeters.z}) {
        fingerprint.real(value);
    }
    const auto& settings = exchange.settings;
    fingerprint.real(settings.timeStepSeconds);
    fingerprint.real(settings.kinematicViscositySquareMetersPerSecond);
    fingerprint.real(settings.minimumWallDistanceMeters);
    fingerprint.real(settings.maximumViscousNumber);
    fingerprint.integer(static_cast<std::uint64_t>(settings.maximumSubsteps));
    fingerprint.real(
        settings.absoluteMomentumToleranceKilogramMetersPerSecond);
    fingerprint.real(settings.relativeMomentumTolerance);
    fingerprint.real(settings.absoluteEnergyToleranceJoules);
    fingerprint.real(settings.relativeEnergyTolerance);
    fingerprint.integer(static_cast<std::uint64_t>(exchange.ownedStorageBytes));
    const auto& diagnostics = exchange.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.quadraturePointCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.substepCount));
    for (const double value : {
             diagnostics.maximumFullStepViscousNumber,
             diagnostics.maximumAcceptedSubstepViscousNumber,
             diagnostics.maximumWallDistanceMeters,
             diagnostics.maximumRelativeTangentialSpeedMetersPerSecond,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.x,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.y,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.z,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.x,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.y,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.z,
             diagnostics.fluidImpulseKilogramMetersPerSecond.x,
             diagnostics.fluidImpulseKilogramMetersPerSecond.y,
             diagnostics.fluidImpulseKilogramMetersPerSecond.z,
             diagnostics.structureImpulseKilogramMetersPerSecond.x,
             diagnostics.structureImpulseKilogramMetersPerSecond.y,
             diagnostics.structureImpulseKilogramMetersPerSecond.z,
             diagnostics.momentumResidualKilogramMetersPerSecond.x,
             diagnostics.momentumResidualKilogramMetersPerSecond.y,
             diagnostics.momentumResidualKilogramMetersPerSecond.z,
             diagnostics.momentumResidualNormKilogramMetersPerSecond,
             diagnostics.kineticEnergyBeforeJoules,
             diagnostics.kineticEnergyAfterJoules,
             diagnostics.wallWorkOnFluidJoules,
             diagnostics.viscousDissipationJoules}) {
        fingerprint.real(value);
    }
    fingerprint.enumeration(diagnostics.failureStage);
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.finite));
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.accepted));
    fingerprint.integer(static_cast<std::uint64_t>(
        exchange.controlVolumes.size()));
    for (const auto& control : exchange.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        for (const double value : {
                 control.volumeCubicMeters,
                 control.incidentWallAreaSquareMeters,
                 control.velocityMetersPerSecond.x,
                 control.velocityMetersPerSecond.y,
                 control.velocityMetersPerSecond.z,
                 control.momentumKilogramMetersPerSecond.x,
                 control.momentumKilogramMetersPerSecond.y,
                 control.momentumKilogramMetersPerSecond.z}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(exchange.samples.size()));
    for (const auto& sample : exchange.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(sample.stableId);
        fingerprint.integer(sample.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.negativeSideControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.positiveSideControlVolumeIndex));
        fingerprint.real(sample.areaSquareMeters);
        for (const double value : {
                 sample.unitNormalNegativeToPositive.x,
                 sample.unitNormalNegativeToPositive.y,
                 sample.unitNormalNegativeToPositive.z,
                 sample.wallVelocityMetersPerSecond.x,
                 sample.wallVelocityMetersPerSecond.y,
                 sample.wallVelocityMetersPerSecond.z,
                 sample.negativeSideFluidImpulseKilogramMetersPerSecond.x,
                 sample.negativeSideFluidImpulseKilogramMetersPerSecond.y,
                 sample.negativeSideFluidImpulseKilogramMetersPerSecond.z,
                 sample.positiveSideFluidImpulseKilogramMetersPerSecond.x,
                 sample.positiveSideFluidImpulseKilogramMetersPerSecond.y,
                 sample.positiveSideFluidImpulseKilogramMetersPerSecond.z,
                 sample.structureTraction.tractionPascals.x,
                 sample.structureTraction.tractionPascals.y,
                 sample.structureTraction.tractionPascals.z}) {
            fingerprint.real(value);
        }
        fingerprint.integer(sample.structureTraction.stableId);
    }
    return fingerprint.value();
}

std::uint64_t acceptedTractionFingerprint(
    const SceneFluidAcceptedWallTractionSet& tractions) {
    Fingerprint fingerprint;
    fingerprint.integer(tractions.version);
    fingerprint.integer(tractions.wallExchangeFingerprint);
    fingerprint.integer(tractions.quadratureFingerprint);
    fingerprint.integer(tractions.surfaceDefinitionFingerprint);
    fingerprint.integer(tractions.surfaceStateFingerprint);
    fingerprint.integer(tractions.structureDefinitionFingerprint);
    fingerprint.integer(tractions.acceptedStepCount);
    fingerprint.real(tractions.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        tractions.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        tractions.tractions.size()));
    for (const auto& traction : tractions.tractions) {
        fingerprint.integer(traction.stableId);
        fingerprint.real(traction.tractionPascals.x);
        fingerprint.real(traction.tractionPascals.y);
        fingerprint.real(traction.tractionPascals.z);
    }
    return fingerprint.value();
}

fluid::Vector3 triangleNormal(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const StableId triangleId) {
    const auto triangleIndex = surface.mappings.triangleIndex(triangleId);
    if (!triangleIndex) {
        throw std::invalid_argument(
            "scene fluid region wall triangle is missing");
    }
    const auto& triangle = surface.triangles[*triangleIndex];
    const auto& first = state.vertices[triangle.vertexIndices[0]].positionMeters;
    const auto& second = state.vertices[triangle.vertexIndices[1]].positionMeters;
    const auto& third = state.vertices[triangle.vertexIndices[2]].positionMeters;
    const fluid::Vector3 a{
        second.x - first.x, second.y - first.y, second.z - first.z};
    const fluid::Vector3 b{
        third.x - first.x, third.y - first.y, third.z - first.z};
    const fluid::Vector3 cross{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
    const double length = norm(cross);
    if (!(length > 0.0) || !std::isfinite(length)) {
        throw std::invalid_argument(
            "scene fluid region wall triangle normal is invalid");
    }
    return scale(cross, 1.0 / length);
}

double periodicSquaredDistance(
    const fluid::Vector3& point,
    const Vec3& controlCentroid,
    const fluid::Vector3& lower,
    const fluid::Vector3& upper) {
    double squared = 0.0;
    const std::array<double, 3> pointValues{
        point.x, point.y, point.z};
    const std::array<double, 3> controlValues{
        controlCentroid.x, controlCentroid.y, controlCentroid.z};
    const std::array<double, 3> lowerValues{
        lower.x, lower.y, lower.z};
    const std::array<double, 3> upperValues{
        upper.x, upper.y, upper.z};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double length = upperValues[axis] - lowerValues[axis];
        double delta = std::abs(pointValues[axis] - controlValues[axis]);
        delta = std::fmod(delta, length);
        delta = std::min(delta, length - delta);
        squared += delta * delta;
    }
    return squared;
}

std::size_t resolveWallControl(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const std::map<std::pair<std::size_t, StableId>, std::size_t>&
        controlByOwner,
    const std::size_t cellIndex,
    const StableId regionId,
    const fluid::Vector3& sampleCentroidMeters,
    const SceneFluidQuadraturePoint& point,
    const char* const side) {
    if (const auto exact = controlByOwner.find({cellIndex, regionId});
        exact != controlByOwner.end()) {
        return exact->second;
    }
    if (cellIndex >= pressureVolumes.cells.size()) {
        throw std::invalid_argument(
            "scene fluid region wall sample is outside the pressure grid");
    }

    const auto counts = pressureVolumes.cellCounts;
    const auto& cell = pressureVolumes.cells[cellIndex].cell;
    const std::array<std::size_t, 6> neighborCells{
        (cell.i + counts.x - 1) % counts.x
            + counts.x * (cell.j + counts.y * cell.k),
        (cell.i + 1) % counts.x
            + counts.x * (cell.j + counts.y * cell.k),
        cell.i + counts.x
            * ((cell.j + counts.y - 1) % counts.y + counts.y * cell.k),
        cell.i + counts.x
            * ((cell.j + 1) % counts.y + counts.y * cell.k),
        cell.i + counts.x
            * (cell.j + counts.y * ((cell.k + counts.z - 1) % counts.z)),
        cell.i + counts.x
            * (cell.j + counts.y * ((cell.k + 1) % counts.z)),
    };
    const SceneFluidPressureControlVolume* nearest = nullptr;
    double nearestSquaredDistance =
        std::numeric_limits<double>::infinity();
    std::array<std::size_t, 6> visited{};
    std::size_t visitedCount = 0;
    for (const std::size_t neighborCell : neighborCells) {
        if (neighborCell == cellIndex
            || std::find(
                   visited.begin(), visited.begin() + visitedCount,
                   neighborCell)
                != visited.begin() + visitedCount) {
            continue;
        }
        visited[visitedCount++] = neighborCell;
        const auto candidateOwner = controlByOwner.find(
            {neighborCell, regionId});
        if (candidateOwner == controlByOwner.end()) {
            continue;
        }
        const auto& candidate = pressureVolumes.controlVolumes[
            candidateOwner->second];
        const double squaredDistance = periodicSquaredDistance(
            sampleCentroidMeters, candidate.centroidMeters,
            pressureVolumes.lowerMeters, pressureVolumes.upperMeters);
        if (squaredDistance < nearestSquaredDistance
            || (squaredDistance == nearestSquaredDistance
                && nearest != nullptr
                && candidate.stableId < nearest->stableId)) {
            nearest = &candidate;
            nearestSquaredDistance = squaredDistance;
        }
    }
    if (nearest == nullptr || !std::isfinite(nearestSquaredDistance)) {
        throw std::invalid_argument(
            "scene fluid region wall sample has no matching "
            + std::string(side) + " control in its cell or immediate "
              "neighbors: sample=" + std::to_string(point.stableId)
            + ", triangle=" + std::to_string(point.triangleId)
            + ", cell=" + std::to_string(cellIndex)
            + ", region=" + std::to_string(regionId));
    }
    return nearest->controlVolumeIndex;
}

} // namespace

static SceneFluidRegionWallExchange exchangeSceneFluidRegionWallMomentumImpl(
    const SceneFluidRegionTransport* const transport,
    const SceneFluidRegionRebase* const rebase,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionWallSettings& settings,
    const SceneFluidRegionWallLimits& limits) {
    const bool usesRebase = rebase != nullptr;
    if (usesRebase == (transport != nullptr)) {
        throw std::invalid_argument(
            "scene fluid region wall source selection is invalid");
    }
    if (usesRebase) {
        validateSceneFluidRegionRebaseIntegrity(*rebase);
    } else {
        validateSceneFluidRegionTransportIntegrity(*transport);
    }
    validateSceneFluidPressureControlVolumeIntegrity(currentPressureVolumes);
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidSurfaceState(surface, currentState);
    validateSceneFluidQuadratureDefinition(quadrature);
    const std::size_t sourceControlCount = usesRebase
        ? rebase->controlVolumes.size() : transport->controlVolumes.size();
    const std::uint64_t sourceAcceptedStepCount = usesRebase
        ? rebase->currentAcceptedStepCount : transport->acceptedStepCount + 1;
    const double sourceSimulationTimeSeconds = usesRebase
        ? rebase->currentSimulationTimeSeconds
        : transport->targetSimulationTimeSeconds;
    const double sourceDensityKgPerCubicMeter = usesRebase
        ? rebase->densityKgPerCubicMeter : transport->densityKgPerCubicMeter;
    if ((!usesRebase && !transport->diagnostics.accepted)
        || !validSettings(settings)
        || grid.cellCounts() != currentPressureVolumes.cellCounts
        || grid.lowerMeters() != currentPressureVolumes.lowerMeters
        || grid.upperMeters() != currentPressureVolumes.upperMeters
        || currentPressureVolumes.surfaceDefinitionFingerprint
            != surface.fingerprint
        || currentPressureVolumes.surfaceStateFingerprint
            != currentState.fingerprint
        || quadrature.surfaceDefinitionFingerprint != surface.fingerprint
        || quadrature.surfaceStateFingerprint != currentState.fingerprint
        || quadrature.structureDefinitionFingerprint
            != currentState.structureDefinitionFingerprint
        || currentPressureVolumes.acceptedStepCount
            != sourceAcceptedStepCount
        || currentState.acceptedStepCount
            != currentPressureVolumes.acceptedStepCount
        || quadrature.acceptedStepCount
            != currentPressureVolumes.acceptedStepCount
        || currentPressureVolumes.simulationTimeSeconds
            != sourceSimulationTimeSeconds
        || currentState.simulationTimeSeconds
            != currentPressureVolumes.simulationTimeSeconds
        || quadrature.simulationTimeSeconds
            != currentPressureVolumes.simulationTimeSeconds
        || sourceControlCount
            != currentPressureVolumes.controlVolumes.size()
        || (usesRebase
            && rebase->currentPressureControlVolumeFingerprint
                != currentPressureVolumes.fingerprint)) {
        throw std::invalid_argument(
            "scene fluid region wall identity is invalid");
    }
    const std::size_t storageBytes = checkedStorageBytes(
        currentPressureVolumes.controlVolumes.size(), quadrature.points.size());
    if (currentPressureVolumes.controlVolumes.size()
            > limits.maximumControlVolumes
        || quadrature.points.size() > limits.maximumQuadraturePoints
        || storageBytes > limits.maximumWallBytes) {
        throw std::length_error(
            "scene fluid region wall exceeds its limits");
    }

    SceneFluidRegionWallExchange result;
    result.sourceTransportFingerprint = usesRebase
        ? rebase->sourceTransportFingerprint : transport->fingerprint;
    result.sourceRegionRebaseFingerprint = usesRebase
        ? rebase->fingerprint : 0;
    result.currentPressureControlVolumeFingerprint =
        currentPressureVolumes.fingerprint;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = currentState.fingerprint;
    result.structureDefinitionFingerprint =
        currentState.structureDefinitionFingerprint;
    result.acceptedStepCount = currentState.acceptedStepCount;
    result.simulationTimeSeconds = currentState.simulationTimeSeconds;
    result.densityKgPerCubicMeter = sourceDensityKgPerCubicMeter;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    std::vector<SceneFluidRegionWallControlVolume> controls;
    controls.reserve(currentPressureVolumes.controlVolumes.size());
    std::map<std::pair<std::size_t, StableId>, std::size_t> controlByOwner;
    for (std::size_t index = 0;
         index < currentPressureVolumes.controlVolumes.size(); ++index) {
        const auto& current = currentPressureVolumes.controlVolumes[index];
        const auto sourceControlVolumeIndex = usesRebase
            ? rebase->controlVolumes[index].controlVolumeIndex
            : transport->controlVolumes[index].controlVolumeIndex;
        const auto sourceStableId = usesRebase
            ? rebase->controlVolumes[index].stableId
            : transport->controlVolumes[index].stableId;
        const auto sourceVelocity = usesRebase
            ? rebase->controlVolumes[index].velocityMetersPerSecond
            : transport->controlVolumes[index].velocityMetersPerSecond;
        if (current.controlVolumeIndex != index
            || sourceControlVolumeIndex != index
            || current.stableId != sourceStableId
            || !(current.volumeCubicMeters > 0.0)
            || !controlByOwner.emplace(
                    std::pair{current.cellIndex, current.regionId}, index)
                    .second) {
            throw std::invalid_argument(
                "scene fluid region wall control topology changed");
        }
        const double mass = result.densityKgPerCubicMeter
            * current.volumeCubicMeters;
        controls.push_back({
            index, current.stableId, current.volumeCubicMeters, 0.0,
            sourceVelocity,
            scale(sourceVelocity, mass),
        });
    }

    const auto kinematics = sampleSceneFluidQuadratureKinematics(
        surface, currentState, quadrature);
    std::vector<SceneFluidRegionWallSample> samples;
    samples.reserve(quadrature.points.size());
    for (std::size_t index = 0; index < quadrature.points.size(); ++index) {
        const auto& point = quadrature.points[index];
        const auto& motion = kinematics[index];
        if (motion.stableId != point.stableId
            || !(point.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region wall sample ownership is invalid");
        }
        const fluid::Vector3 sampleCentroidMeters{
            motion.positionMeters.x,
            motion.positionMeters.y,
            motion.positionMeters.z,
        };
        const std::size_t negative = resolveWallControl(
            currentPressureVolumes, controlByOwner,
            point.negativeSideCellIndex, point.negativeSideRegionId,
            sampleCentroidMeters, point, "negative-side");
        const std::size_t positive = resolveWallControl(
            currentPressureVolumes, controlByOwner,
            point.positiveSideCellIndex, point.positiveSideRegionId,
            sampleCentroidMeters, point, "positive-side");
        const auto normal = triangleNormal(
            surface, currentState, point.triangleId);
        SceneFluidRegionWallSample sample;
        sample.sampleIndex = index;
        sample.stableId = point.stableId;
        sample.triangleId = point.triangleId;
        sample.negativeSideControlVolumeIndex = negative;
        sample.positiveSideControlVolumeIndex = positive;
        sample.areaSquareMeters = point.areaSquareMeters;
        sample.unitNormalNegativeToPositive = normal;
        sample.wallVelocityMetersPerSecond = {
            motion.velocityMetersPerSecond.x,
            motion.velocityMetersPerSecond.y,
            motion.velocityMetersPerSecond.z,
        };
        sample.structureTraction.stableId = point.stableId;
        samples.push_back(sample);
        controls[negative].incidentWallAreaSquareMeters +=
            point.areaSquareMeters;
        controls[positive].incidentWallAreaSquareMeters +=
            point.areaSquareMeters;
    }

    auto kernel = detail::exchangeSceneFluidWallMomentumKernel(
        result.densityKgPerCubicMeter, std::move(controls),
        std::move(samples), settings, storageBytes);
    result.ownedStorageBytes = kernel.ownedStorageBytes;
    result.diagnostics = std::move(kernel.diagnostics);
    result.controlVolumes = std::move(kernel.controlVolumes);
    result.samples = std::move(kernel.samples);
    result.fingerprint = exchangeFingerprint(result);
    validateSceneFluidRegionWallExchangeIntegrity(result);
    return result;
}

SceneFluidRegionWallExchange exchangeSceneFluidRegionWallMomentum(
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionWallSettings& settings,
    const SceneFluidRegionWallLimits& limits) {
    return exchangeSceneFluidRegionWallMomentumImpl(
        &transport, nullptr, grid, currentPressureVolumes, surface,
        currentState, quadrature, settings, limits);
}

SceneFluidRegionWallExchange exchangeSceneFluidRegionWallMomentum(
    const SceneFluidRegionRebase& rebase,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionWallSettings& settings,
    const SceneFluidRegionWallLimits& limits) {
    return exchangeSceneFluidRegionWallMomentumImpl(
        nullptr, &rebase, grid, currentPressureVolumes, surface,
        currentState, quadrature, settings, limits);
}

void validateSceneFluidRegionWallExchangeIntegrity(
    const SceneFluidRegionWallExchange& exchange) {
    const auto& diagnostics = exchange.diagnostics;
    bool controlsValid = true;
    for (std::size_t index = 0; index < exchange.controlVolumes.size(); ++index) {
        const auto& control = exchange.controlVolumes[index];
        controlsValid = controlsValid
            && control.controlVolumeIndex == index
            && control.stableId != 0
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && std::isfinite(control.incidentWallAreaSquareMeters)
            && control.incidentWallAreaSquareMeters >= 0.0
            && finite(control.velocityMetersPerSecond)
            && finite(control.momentumKilogramMetersPerSecond);
    }
    bool samplesValid = true;
    for (std::size_t index = 0; index < exchange.samples.size(); ++index) {
        const auto& sample = exchange.samples[index];
        samplesValid = samplesValid
            && sample.sampleIndex == index
            && sample.stableId != 0
            && sample.triangleId != invalidStableId
            && sample.negativeSideControlVolumeIndex
                < exchange.controlVolumes.size()
            && sample.positiveSideControlVolumeIndex
                < exchange.controlVolumes.size()
            && std::isfinite(sample.areaSquareMeters)
            && sample.areaSquareMeters > 0.0
            && finite(sample.unitNormalNegativeToPositive)
            && std::abs(norm(sample.unitNormalNegativeToPositive) - 1.0)
                <= 1.0e-12
            && finite(sample.wallVelocityMetersPerSecond)
            && finite(sample.negativeSideFluidImpulseKilogramMetersPerSecond)
            && finite(sample.positiveSideFluidImpulseKilogramMetersPerSecond)
            && sample.structureTraction.stableId == sample.stableId
            && std::isfinite(sample.structureTraction.tractionPascals.x)
            && std::isfinite(sample.structureTraction.tractionPascals.y)
            && std::isfinite(sample.structureTraction.tractionPascals.z);
    }
    const bool acceptedShape = diagnostics.accepted
        && diagnostics.failureStage == SceneFluidRegionWallFailureStage::None
        && exchange.controlVolumes.size() == diagnostics.controlVolumeCount
        && exchange.samples.size() == diagnostics.quadraturePointCount
        && exchange.ownedStorageBytes == checkedStorageBytes(
            exchange.controlVolumes.size(), exchange.samples.size());
    const bool rejectedShape = !diagnostics.accepted
        && diagnostics.failureStage != SceneFluidRegionWallFailureStage::None
        && exchange.controlVolumes.empty() && exchange.samples.empty()
        && exchange.ownedStorageBytes == 0;
    bool ledgersValid = true;
    if (acceptedShape) {
        fluid::Vector3 sampleFluidImpulse;
        fluid::Vector3 sampleStructureImpulse;
        for (const auto& sample : exchange.samples) {
            const auto fluidImpulse = add(
                sample.negativeSideFluidImpulseKilogramMetersPerSecond,
                sample.positiveSideFluidImpulseKilogramMetersPerSecond);
            const auto structureImpulse = scale(fluidImpulse, -1.0);
            sampleFluidImpulse = add(sampleFluidImpulse, fluidImpulse);
            sampleStructureImpulse = add(
                sampleStructureImpulse, structureImpulse);
            const auto expectedTraction = scale(
                structureImpulse,
                1.0 / (sample.areaSquareMeters
                       * exchange.settings.timeStepSeconds));
            ledgersValid = ledgersValid
                && sample.structureTraction.tractionPascals.x
                    == expectedTraction.x
                && sample.structureTraction.tractionPascals.y
                    == expectedTraction.y
                && sample.structureTraction.tractionPascals.z
                    == expectedTraction.z;
        }
        const auto controlMomentum = totalMomentum(exchange.controlVolumes);
        const double controlEnergy = kineticEnergy(
            exchange.controlVolumes, exchange.densityKgPerCubicMeter);
        const auto expectedFluidImpulse = subtract(
            diagnostics.fluidMomentumAfterKilogramMetersPerSecond,
            diagnostics.fluidMomentumBeforeKilogramMetersPerSecond);
        const auto expectedResidual = add(
            diagnostics.fluidImpulseKilogramMetersPerSecond,
            diagnostics.structureImpulseKilogramMetersPerSecond);
        const double momentumTolerance = tolerance(
            exchange.settings
                .absoluteMomentumToleranceKilogramMetersPerSecond,
            exchange.settings.relativeMomentumTolerance,
            std::max(
                norm(diagnostics.fluidImpulseKilogramMetersPerSecond),
                norm(sampleFluidImpulse)));
        const double energyResidual =
            diagnostics.kineticEnergyBeforeJoules
            + diagnostics.wallWorkOnFluidJoules
            - diagnostics.kineticEnergyAfterJoules
            - diagnostics.viscousDissipationJoules;
        const double energyTolerance = tolerance(
            exchange.settings.absoluteEnergyToleranceJoules,
            exchange.settings.relativeEnergyTolerance,
            std::max({
                diagnostics.kineticEnergyBeforeJoules,
                diagnostics.kineticEnergyAfterJoules,
                std::abs(diagnostics.wallWorkOnFluidJoules),
                diagnostics.viscousDissipationJoules,
            }));
        ledgersValid = ledgersValid
            && controlMomentum
                == diagnostics.fluidMomentumAfterKilogramMetersPerSecond
            && controlEnergy == diagnostics.kineticEnergyAfterJoules
            && expectedFluidImpulse
                == diagnostics.fluidImpulseKilogramMetersPerSecond
            && sampleStructureImpulse
                == diagnostics.structureImpulseKilogramMetersPerSecond
            && expectedResidual
                == diagnostics.momentumResidualKilogramMetersPerSecond
            && norm(expectedResidual)
                == diagnostics
                    .momentumResidualNormKilogramMetersPerSecond
            && norm(subtract(
                   sampleFluidImpulse,
                   diagnostics.fluidImpulseKilogramMetersPerSecond))
                <= momentumTolerance
            && diagnostics.viscousDissipationJoules >= 0.0
            && std::abs(energyResidual) <= energyTolerance
            && diagnostics.substepCount <= exchange.settings.maximumSubsteps
            && diagnostics.maximumAcceptedSubstepViscousNumber
                <= exchange.settings.maximumViscousNumber;
    }
    if (exchange.version != sceneFluidRegionWallExchangeVersion
        || exchange.fingerprint == 0
        || exchange.sourceTransportFingerprint == 0
        || exchange.currentPressureControlVolumeFingerprint == 0
        || exchange.quadratureFingerprint == 0
        || exchange.surfaceDefinitionFingerprint == 0
        || exchange.surfaceStateFingerprint == 0
        || exchange.structureDefinitionFingerprint == 0
        || !std::isfinite(exchange.simulationTimeSeconds)
        || !(exchange.densityKgPerCubicMeter > 0.0)
        || exchange.cellCounts.x == 0 || exchange.cellCounts.y == 0
        || exchange.cellCounts.z == 0
        || !finite(exchange.lowerMeters) || !finite(exchange.upperMeters)
        || !validSettings(exchange.settings)
        || diagnostics.controlVolumeCount == 0
        || diagnostics.quadraturePointCount == 0
        || diagnostics.substepCount == 0
        || !std::isfinite(diagnostics.maximumFullStepViscousNumber)
        || !std::isfinite(diagnostics.maximumAcceptedSubstepViscousNumber)
        || !std::isfinite(diagnostics.maximumWallDistanceMeters)
        || !std::isfinite(
            diagnostics.maximumRelativeTangentialSpeedMetersPerSecond)
        || !finite(diagnostics.fluidMomentumBeforeKilogramMetersPerSecond)
        || !finite(diagnostics.fluidMomentumAfterKilogramMetersPerSecond)
        || !finite(diagnostics.fluidImpulseKilogramMetersPerSecond)
        || !finite(diagnostics.structureImpulseKilogramMetersPerSecond)
        || !finite(diagnostics.momentumResidualKilogramMetersPerSecond)
        || !std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        || !std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        || !std::isfinite(diagnostics.kineticEnergyAfterJoules)
        || !std::isfinite(diagnostics.wallWorkOnFluidJoules)
        || !std::isfinite(diagnostics.viscousDissipationJoules)
        || !diagnostics.finite || (!acceptedShape && !rejectedShape)
        || !controlsValid || !samplesValid || !ledgersValid
        || exchange.fingerprint != exchangeFingerprint(exchange)) {
        throw std::invalid_argument(
            "scene fluid region wall exchange integrity is invalid");
    }
}

void validateSceneFluidRegionWallExchange(
    const SceneFluidRegionWallExchange& exchange,
    const SceneFluidRegionTransport& transport,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature) {
    validateSceneFluidRegionWallExchangeIntegrity(exchange);
    if (exchange.sourceTransportFingerprint != transport.fingerprint
        || exchange.sourceRegionRebaseFingerprint != 0
        || exchange.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || exchange.quadratureFingerprint != quadrature.fingerprint
        || exchange.surfaceDefinitionFingerprint != surface.fingerprint
        || exchange.surfaceStateFingerprint != currentState.fingerprint
        || exchange.structureDefinitionFingerprint
            != currentState.structureDefinitionFingerprint
        || exchange.acceptedStepCount != currentState.acceptedStepCount
        || exchange.simulationTimeSeconds != currentState.simulationTimeSeconds
        || exchange.cellCounts != grid.cellCounts()
        || exchange.lowerMeters != grid.lowerMeters()
        || exchange.upperMeters != grid.upperMeters()
        || exchange.densityKgPerCubicMeter
            != transport.densityKgPerCubicMeter) {
        throw std::invalid_argument(
            "scene fluid region wall exchange binding is invalid");
    }
}

void validateSceneFluidRegionWallExchange(
    const SceneFluidRegionWallExchange& exchange,
    const SceneFluidRegionRebase& rebase,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& currentState,
    const SceneFluidQuadratureDefinition& quadrature) {
    validateSceneFluidRegionWallExchangeIntegrity(exchange);
    validateSceneFluidRegionRebaseIntegrity(rebase);
    if (exchange.sourceTransportFingerprint
            != rebase.sourceTransportFingerprint
        || exchange.sourceRegionRebaseFingerprint != rebase.fingerprint
        || exchange.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || rebase.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || exchange.quadratureFingerprint != quadrature.fingerprint
        || exchange.surfaceDefinitionFingerprint != surface.fingerprint
        || exchange.surfaceStateFingerprint != currentState.fingerprint
        || exchange.structureDefinitionFingerprint
            != currentState.structureDefinitionFingerprint
        || exchange.acceptedStepCount != currentState.acceptedStepCount
        || exchange.simulationTimeSeconds != currentState.simulationTimeSeconds
        || exchange.cellCounts != grid.cellCounts()
        || exchange.lowerMeters != grid.lowerMeters()
        || exchange.upperMeters != grid.upperMeters()
        || exchange.densityKgPerCubicMeter
            != rebase.densityKgPerCubicMeter) {
        throw std::invalid_argument(
            "scene fluid rebased wall exchange binding is invalid");
    }
}

SceneFluidAcceptedWallTractionSet captureSceneFluidAcceptedWallTractions(
    const SceneFluidRegionWallExchange& exchange) {
    validateSceneFluidRegionWallExchangeIntegrity(exchange);
    if (!exchange.diagnostics.accepted) {
        throw std::invalid_argument(
            "scene fluid accepted wall traction source was not accepted");
    }
    SceneFluidAcceptedWallTractionSet result;
    result.wallExchangeFingerprint = exchange.fingerprint;
    result.quadratureFingerprint = exchange.quadratureFingerprint;
    result.surfaceDefinitionFingerprint =
        exchange.surfaceDefinitionFingerprint;
    result.surfaceStateFingerprint = exchange.surfaceStateFingerprint;
    result.structureDefinitionFingerprint =
        exchange.structureDefinitionFingerprint;
    result.acceptedStepCount = exchange.acceptedStepCount;
    result.simulationTimeSeconds = exchange.simulationTimeSeconds;
    if (exchange.samples.size() > std::numeric_limits<std::size_t>::max()
                                    / sizeof(SceneFluidQuadratureTraction)) {
        throw std::length_error(
            "scene fluid accepted wall traction storage overflows");
    }
    result.ownedStorageBytes = exchange.samples.size()
        * sizeof(SceneFluidQuadratureTraction);
    result.tractions.reserve(exchange.samples.size());
    for (const auto& sample : exchange.samples) {
        result.tractions.push_back(sample.structureTraction);
    }
    result.fingerprint = acceptedTractionFingerprint(result);
    validateSceneFluidAcceptedWallTractionSetIntegrity(result);
    return result;
}

SceneFluidAcceptedWallTractionSet captureSceneFluidAcceptedWallTractions(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange) {
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(exchange);
    if (!exchange.diagnostics.accepted) {
        throw std::invalid_argument(
            "scene fluid accepted regional wall traction source was not accepted");
    }
    const auto& input = exchange.sourceInput;
    SceneFluidAcceptedWallTractionSet result;
    result.wallExchangeFingerprint = exchange.fingerprint;
    result.quadratureFingerprint = input.quadratureFingerprint;
    result.surfaceDefinitionFingerprint = input.surfaceDefinitionFingerprint;
    result.surfaceStateFingerprint = input.surfaceStateFingerprint;
    result.structureDefinitionFingerprint =
        input.structureDefinitionFingerprint;
    result.acceptedStepCount = input.acceptedStepCount;
    result.simulationTimeSeconds = input.simulationTimeSeconds;
    if (exchange.samples.size() > std::numeric_limits<std::size_t>::max()
                                    / sizeof(SceneFluidQuadratureTraction)) {
        throw std::length_error(
            "scene fluid accepted regional wall traction storage overflows");
    }
    result.ownedStorageBytes = exchange.samples.size()
        * sizeof(SceneFluidQuadratureTraction);
    result.tractions.reserve(exchange.samples.size());
    for (const auto& sample : exchange.samples) {
        result.tractions.push_back(sample.structureTraction);
    }
    result.fingerprint = acceptedTractionFingerprint(result);
    validateSceneFluidAcceptedWallTractionSetIntegrity(result);
    return result;
}

void validateSceneFluidAcceptedWallTractionSetIntegrity(
    const SceneFluidAcceptedWallTractionSet& tractions) {
    bool valuesValid = true;
    std::map<std::uint64_t, bool> stableIds;
    for (const auto& traction : tractions.tractions) {
        valuesValid = valuesValid
            && traction.stableId != 0
            && stableIds.emplace(traction.stableId, true).second
            && std::isfinite(traction.tractionPascals.x)
            && std::isfinite(traction.tractionPascals.y)
            && std::isfinite(traction.tractionPascals.z);
    }
    if (tractions.version != sceneFluidAcceptedWallTractionVersion
        || tractions.fingerprint == 0
        || tractions.wallExchangeFingerprint == 0
        || tractions.quadratureFingerprint == 0
        || tractions.surfaceDefinitionFingerprint == 0
        || tractions.surfaceStateFingerprint == 0
        || tractions.structureDefinitionFingerprint == 0
        || !std::isfinite(tractions.simulationTimeSeconds)
        || tractions.tractions.empty()
        || tractions.ownedStorageBytes
            != tractions.tractions.size()
                * sizeof(SceneFluidQuadratureTraction)
        || !valuesValid
        || tractions.fingerprint != acceptedTractionFingerprint(tractions)) {
        throw std::invalid_argument(
            "scene fluid accepted wall traction integrity is invalid");
    }
}

void validateSceneFluidAcceptedWallTractions(
    const SceneFluidAcceptedWallTractionSet& tractions,
    const SceneFluidQuadratureDefinition& quadrature,
    const std::uint64_t expectedWallExchangeFingerprint) {
    validateSceneFluidAcceptedWallTractionSetIntegrity(tractions);
    validateSceneFluidQuadratureDefinition(quadrature);
    if (tractions.wallExchangeFingerprint != expectedWallExchangeFingerprint
        || tractions.quadratureFingerprint != quadrature.fingerprint
        || tractions.surfaceDefinitionFingerprint
            != quadrature.surfaceDefinitionFingerprint
        || tractions.surfaceStateFingerprint
            != quadrature.surfaceStateFingerprint
        || tractions.structureDefinitionFingerprint
            != quadrature.structureDefinitionFingerprint
        || tractions.acceptedStepCount != quadrature.acceptedStepCount
        || tractions.simulationTimeSeconds
            != quadrature.simulationTimeSeconds
        || tractions.tractions.size() != quadrature.points.size()) {
        throw std::invalid_argument(
            "scene fluid accepted wall traction binding is invalid");
    }
    for (std::size_t index = 0; index < tractions.tractions.size(); ++index) {
        if (tractions.tractions[index].stableId
            != quadrature.points[index].stableId) {
            throw std::invalid_argument(
                "scene fluid accepted wall traction sample is foreign");
        }
    }
}

} // namespace simwing::fsi
