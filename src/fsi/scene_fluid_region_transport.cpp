#include "scene_fluid_region_transport.h"

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

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

fluid::Vector3 add(const fluid::Vector3& first,
                   const fluid::Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

fluid::Vector3 scale(const fluid::Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double norm(const fluid::Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

double maximumAbsoluteComponent(const fluid::Vector3& value) {
    return std::max({std::abs(value.x), std::abs(value.y),
                     std::abs(value.z)});
}

fluid::Vector3 cellCenteredVelocity(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity,
    const std::size_t cellIndex) {
    const auto counts = grid.cellCounts();
    const std::size_t i = cellIndex % counts.x;
    const std::size_t j = (cellIndex / counts.x) % counts.y;
    const std::size_t k = cellIndex / (counts.x * counts.y);
    return {
        0.5 * (velocity.xFaces()[grid.cellIndex(i, j, k)]
               + velocity.xFaces()[grid.cellIndex(
                   (i + 1) % counts.x, j, k)]),
        0.5 * (velocity.yFaces()[grid.cellIndex(i, j, k)]
               + velocity.yFaces()[grid.cellIndex(
                   i, (j + 1) % counts.y, k)]),
        0.5 * (velocity.zFaces()[grid.cellIndex(i, j, k)]
               + velocity.zFaces()[grid.cellIndex(
                   i, j, (k + 1) % counts.z)]),
    };
}

double tolerance(const double absolute,
                 const double relative,
                 const double reference) {
    return absolute + relative * std::max(1.0, std::abs(reference));
}

double kineticEnergy(
    const std::span<const SceneFluidRegionTransportControlVolume> controls,
    const double density) {
    double result = 0.0;
    for (const auto& control : controls) {
        result += 0.5 * density * control.volumeCubicMeters
            * (control.velocityMetersPerSecond.x
                   * control.velocityMetersPerSecond.x
               + control.velocityMetersPerSecond.y
                   * control.velocityMetersPerSecond.y
               + control.velocityMetersPerSecond.z
                   * control.velocityMetersPerSecond.z);
    }
    return result;
}

fluid::Vector3 totalMomentum(
    const std::span<const SceneFluidRegionTransportControlVolume> controls) {
    fluid::Vector3 result;
    for (const auto& control : controls) {
        result = add(result, control.momentumKilogramMetersPerSecond);
    }
    return result;
}

std::size_t storageBytesForControlVolumes(const std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max()
                    / sizeof(SceneFluidRegionTransportControlVolume)) {
        throw std::length_error(
            "scene fluid region transport storage size overflows");
    }
    return count * sizeof(SceneFluidRegionTransportControlVolume);
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "scene fluid region transport working storage size overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "scene fluid region transport working storage size overflows");
    }
    return first + second;
}

std::size_t workingStorageBytes(const std::size_t controlVolumeCount,
                                const std::size_t linkCount) {
    const std::size_t perControl =
        sizeof(SceneFluidRegionTransportControlVolume)
        + 4 * sizeof(double) + sizeof(fluid::Vector3);
    return checkedAdd(
        checkedMultiply(controlVolumeCount, perControl),
        checkedMultiply(linkCount, sizeof(double)));
}

std::uint64_t transportFingerprint(
    const SceneFluidRegionTransport& transport) {
    Fingerprint fingerprint;
    fingerprint.integer(transport.version);
    fingerprint.integer(transport.sourceMomentumFingerprint);
    fingerprint.integer(transport.pressureProjectionFingerprint);
    fingerprint.integer(transport.pressureFaceLinkFingerprint);
    fingerprint.integer(transport.pressureVolumeRateFingerprint);
    fingerprint.integer(transport.previousBulkVelocityFingerprint);
    fingerprint.integer(transport.currentBulkVelocityFingerprint);
    fingerprint.integer(transport.acceptedStepCount);
    fingerprint.real(transport.sourceSimulationTimeSeconds);
    fingerprint.real(transport.targetSimulationTimeSeconds);
    fingerprint.real(transport.densityKgPerCubicMeter);
    const auto& settings = transport.settings;
    fingerprint.real(settings.timeStepSeconds);
    fingerprint.real(settings.kinematicViscositySquareMetersPerSecond);
    fingerprint.real(settings.maximumOutgoingCourantNumber);
    fingerprint.real(settings.maximumViscousNumber);
    fingerprint.integer(static_cast<std::uint64_t>(
        settings.maximumSubsteps));
    fingerprint.real(
        settings.absoluteMomentumToleranceKilogramMetersPerSecond);
    fingerprint.real(settings.relativeMomentumTolerance);
    fingerprint.real(settings.absoluteEnergyToleranceJoules);
    fingerprint.real(settings.relativeEnergyTolerance);
    fingerprint.integer(static_cast<std::uint64_t>(
        transport.ownedStorageBytes));
    const auto& diagnostics = transport.diagnostics;
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.linkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.openingLinkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.substepCount));
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.usesMovingVolumeRates ? 1 : 0));
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.usesBulkVelocityIncrement ? 1 : 0));
    fingerprint.real(
        diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond);
    fingerprint.real(
        diagnostics.maximumAbsoluteGeometryVolumeChangeCubicMeters);
    fingerprint.real(
        diagnostics.maximumBulkVelocityIncrementMetersPerSecond);
    fingerprint.real(
        diagnostics.bulkVelocityIncrementImpulseKilogramMetersPerSecond.x);
    fingerprint.real(
        diagnostics.bulkVelocityIncrementImpulseKilogramMetersPerSecond.y);
    fingerprint.real(
        diagnostics.bulkVelocityIncrementImpulseKilogramMetersPerSecond.z);
    fingerprint.real(diagnostics.bulkVelocityIncrementWorkJoules);
    fingerprint.real(
        diagnostics.maximumCorrectedContinuityResidualCubicMetersPerSecond);
    fingerprint.real(
        diagnostics.maximumFullStepOutgoingCourantNumber);
    fingerprint.real(
        diagnostics.maximumAcceptedSubstepOutgoingCourantNumber);
    fingerprint.real(diagnostics.maximumFullStepViscousNumber);
    fingerprint.real(diagnostics.maximumAcceptedSubstepViscousNumber);
    for (const double value : {
             diagnostics.momentumBeforeKilogramMetersPerSecond.x,
             diagnostics.momentumBeforeKilogramMetersPerSecond.y,
             diagnostics.momentumBeforeKilogramMetersPerSecond.z,
             diagnostics.momentumAfterKilogramMetersPerSecond.x,
             diagnostics.momentumAfterKilogramMetersPerSecond.y,
             diagnostics.momentumAfterKilogramMetersPerSecond.z,
             diagnostics.momentumResidualKilogramMetersPerSecond.x,
             diagnostics.momentumResidualKilogramMetersPerSecond.y,
             diagnostics.momentumResidualKilogramMetersPerSecond.z,
             diagnostics.momentumResidualNormKilogramMetersPerSecond,
             diagnostics.kineticEnergyBeforeJoules,
             diagnostics.kineticEnergyAfterAdvectionJoules,
             diagnostics.kineticEnergyAfterJoules,
             diagnostics.advectiveEnergyLossJoules,
             diagnostics.viscousEnergyLossJoules,
             diagnostics.maximumVelocityChangeMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.enumeration(diagnostics.failureStage);
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.finite ? 1 : 0));
    fingerprint.integer(static_cast<std::uint8_t>(
        diagnostics.accepted ? 1 : 0));
    fingerprint.integer(static_cast<std::uint64_t>(
        transport.controlVolumes.size()));
    for (const auto& control : transport.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
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
    }
    return fingerprint.value();
}

bool validSettings(const SceneFluidRegionTransportSettings& settings) {
    return std::isfinite(settings.timeStepSeconds)
        && settings.timeStepSeconds > 0.0
        && std::isfinite(
            settings.kinematicViscositySquareMetersPerSecond)
        && settings.kinematicViscositySquareMetersPerSecond >= 0.0
        && std::isfinite(settings.maximumOutgoingCourantNumber)
        && settings.maximumOutgoingCourantNumber > 0.0
        && settings.maximumOutgoingCourantNumber <= 1.0
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

std::size_t requiredSubsteps(const double fullStepValue,
                             const double maximumValue) {
    if (!(fullStepValue > maximumValue)) {
        return 1;
    }
    const double required = std::ceil(fullStepValue / maximumValue);
    if (!std::isfinite(required)
        || required > static_cast<double>(
            std::numeric_limits<std::size_t>::max())) {
        return std::numeric_limits<std::size_t>::max();
    }
    return static_cast<std::size_t>(required);
}

void updateVelocities(
    std::vector<SceneFluidRegionTransportControlVolume>& controls,
    const double density) {
    for (auto& control : controls) {
        const double inverseMass =
            1.0 / (density * control.volumeCubicMeters);
        control.velocityMetersPerSecond = scale(
            control.momentumKilogramMetersPerSecond, inverseMass);
    }
}

struct RegionTransportControlSource {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t componentIndex = 0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;
};

struct RegionTransportLinkSource {
    std::size_t linkIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidPressureFaceLinkKind kind =
        SceneFluidPressureFaceLinkKind::SameRegion;
    std::size_t minusControlVolumeIndex = 0;
    std::size_t plusControlVolumeIndex = 0;
    double correctedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
};

struct RegionTransportFlowSource {
    std::uint64_t fingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t pressureVolumeRateFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    double absoluteContinuityToleranceCubicMetersPerSecond = 0.0;
    double relativeContinuityTolerance = 0.0;
    double predictedContinuityResidualMaximumCubicMetersPerSecond = 0.0;
    bool usesMovingVolumeRates = false;
    std::vector<RegionTransportControlSource> controls;
    std::vector<RegionTransportLinkSource> links;
};

RegionTransportFlowSource regionTransportFlowSource(
    const SceneFluidPressureProjection& projection) {
    validateSceneFluidPressureProjectionIntegrity(projection);
    RegionTransportFlowSource result;
    result.fingerprint = projection.fingerprint;
    result.pressureFaceLinkFingerprint =
        projection.pressureFaceLinkFingerprint;
    result.pressureVolumeRateFingerprint =
        projection.pressureVolumeRateFingerprint;
    result.acceptedStepCount = projection.acceptedStepCount;
    result.simulationTimeSeconds = projection.simulationTimeSeconds;
    result.densityKgPerCubicMeter =
        projection.settings.densityKgPerCubicMeter;
    result.absoluteContinuityToleranceCubicMetersPerSecond =
        projection.settings
            .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond;
    result.relativeContinuityTolerance =
        projection.settings.relativeCorrectedVolumeRateTolerance;
    result.predictedContinuityResidualMaximumCubicMetersPerSecond =
        projection.diagnostics
            .predictedContinuityResidualMaximumCubicMetersPerSecond;
    result.usesMovingVolumeRates =
        projection.diagnostics.usesMovingVolumeRates;
    result.controls.reserve(projection.controlVolumes.size());
    for (const auto& source : projection.controlVolumes) {
        result.controls.push_back({
            source.controlVolumeIndex,
            source.stableId,
            source.componentIndex,
            source.geometryVolumeChangeRateCubicMetersPerSecond,
        });
    }
    result.links.reserve(projection.links.size());
    for (const auto& source : projection.links) {
        result.links.push_back({
            source.linkIndex,
            source.stableId,
            source.kind,
            source.minusControlVolumeIndex,
            source.plusControlVolumeIndex,
            source.correctedRelativeVolumeFlowRateCubicMetersPerSecond,
        });
    }
    return result;
}

RegionTransportFlowSource regionTransportFlowSource(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow) {
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(correctedFlow);
    if (!correctedFlow.accepted
        || sourceMomentum.pressureProjectionFingerprint
            != correctedFlow.fingerprint
        || sourceMomentum.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || correctedFlow.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || correctedFlow.acceptedStepCount
            != sourceMomentum.acceptedStepCount
        || correctedFlow.simulationTimeSeconds
            != sourceMomentum.simulationTimeSeconds
        || correctedFlow.densityKgPerCubicMeter
            != sourceMomentum.densityKgPerCubicMeter
        || correctedFlow.traces.size() != faceLinks.links.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic region transport source is foreign");
    }

    std::map<std::uint64_t, std::size_t> embeddedLinkByPatchId;
    for (const auto& link : faceLinks.links) {
        if (link.geometryKind
            == SceneFluidPressureLinkGeometryKind::EmbeddedOpening
            && (link.openingPatchStableId == 0
                || !embeddedLinkByPatchId.emplace(
                        link.openingPatchStableId,
                        link.linkIndex).second)) {
            throw std::invalid_argument(
                "scene fluid mimetic region transport embedded identity is invalid");
        }
    }
    std::vector<const SceneFluidMimeticCorrectedTrace*> traceByLink(
        faceLinks.links.size(), nullptr);
    for (const auto& trace : correctedFlow.traces) {
        std::size_t linkIndex = trace.sourceIndex;
        if (trace.kind
            == SceneFluidMimeticHalfFaceKind::AuthoredOpeningTrace) {
            const auto found = embeddedLinkByPatchId.find(
                trace.sourceStableId);
            if (found == embeddedLinkByPatchId.end()) {
                throw std::invalid_argument(
                    "scene fluid mimetic region transport embedded trace is missing");
            }
            linkIndex = found->second;
        } else if (trace.kind
                   != SceneFluidMimeticHalfFaceKind::CartesianTrace) {
            throw std::invalid_argument(
                "scene fluid mimetic region transport trace kind is invalid");
        }
        if (linkIndex >= traceByLink.size()
            || traceByLink[linkIndex] != nullptr) {
            throw std::invalid_argument(
                "scene fluid mimetic region transport trace ownership is invalid");
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
                "scene fluid mimetic region transport trace binding is invalid");
        }
        traceByLink[linkIndex] = &trace;
    }

    RegionTransportFlowSource result;
    result.fingerprint = correctedFlow.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.acceptedStepCount = correctedFlow.acceptedStepCount;
    result.simulationTimeSeconds = correctedFlow.simulationTimeSeconds;
    result.densityKgPerCubicMeter = correctedFlow.densityKgPerCubicMeter;
    result.absoluteContinuityToleranceCubicMetersPerSecond =
        correctedFlow.correctedContinuityToleranceCubicMetersPerSecond;
    result.controls.reserve(sourceMomentum.controlVolumes.size());
    for (const auto& source : sourceMomentum.controlVolumes) {
        result.controls.push_back({
            source.controlVolumeIndex,
            source.stableId,
            source.componentIndex,
            0.0,
        });
    }
    result.links.reserve(faceLinks.links.size());
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        const auto* trace = traceByLink[index];
        if (trace == nullptr) {
            throw std::invalid_argument(
                "scene fluid mimetic region transport missed a link");
        }
        result.links.push_back({
            source.linkIndex,
            source.stableId,
            source.kind,
            source.minusControlVolumeIndex,
            source.plusControlVolumeIndex,
            trace->correctedRelativeVolumeFlowRateCubicMetersPerSecond,
        });
    }
    return result;
}

} // namespace

static SceneFluidRegionTransport advanceSceneFluidRegionMomentumImpl(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const RegionTransportFlowSource& flowSource,
    const fluid::PeriodicCartesianGrid* const bulkGrid,
    const fluid::MacVelocityField* const previousBulkVelocity,
    const fluid::MacVelocityField* const currentBulkVelocity,
    const SceneFluidRegionTransportSettings& settings,
    const SceneFluidRegionTransportLimits& limits) {
    validateSceneFluidRegionMomentumStateIntegrity(sourceMomentum);
    if (!validSettings(settings)) {
        throw std::invalid_argument(
            "scene fluid region transport settings are invalid");
    }
    const bool usesBulkVelocityIncrement = bulkGrid != nullptr;
    if (usesBulkVelocityIncrement
            != (previousBulkVelocity != nullptr
                && currentBulkVelocity != nullptr)
        || (usesBulkVelocityIncrement
            && (bulkGrid->cellCounts() != sourceMomentum.cellCounts
                || bulkGrid->lowerMeters() != sourceMomentum.lowerMeters
                || bulkGrid->upperMeters() != sourceMomentum.upperMeters
                || !previousBulkVelocity->matches(*bulkGrid)
                || !currentBulkVelocity->matches(*bulkGrid)
                || !fluid::isFinite(*previousBulkVelocity)
                || !fluid::isFinite(*currentBulkVelocity)))) {
        throw std::invalid_argument(
            "scene fluid region transport bulk velocity increment is foreign");
    }
    if (faceLinks.version != sceneFluidPressureFaceLinkVersion
        || faceLinks.fingerprint == 0
        || sourceMomentum.pressureProjectionFingerprint
            != flowSource.fingerprint
        || sourceMomentum.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || flowSource.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || sourceMomentum.acceptedStepCount
            != flowSource.acceptedStepCount
        || sourceMomentum.simulationTimeSeconds
            != flowSource.simulationTimeSeconds
        || sourceMomentum.densityKgPerCubicMeter
            != flowSource.densityKgPerCubicMeter
        || sourceMomentum.controlVolumes.size()
            != flowSource.controls.size()
        || flowSource.links.size() != faceLinks.links.size()) {
        throw std::invalid_argument(
            "scene fluid region transport identity is invalid");
    }
    const std::size_t storageBytes = storageBytesForControlVolumes(
        sourceMomentum.controlVolumes.size());
    const std::size_t workingBytes = workingStorageBytes(
        sourceMomentum.controlVolumes.size(), faceLinks.links.size());
    if (sourceMomentum.controlVolumes.size() > limits.maximumControlVolumes
        || faceLinks.links.size() > limits.maximumLinks
        || workingBytes > limits.maximumTransportBytes) {
        throw std::length_error(
            "scene fluid region transport exceeds its limits");
    }

    SceneFluidRegionTransport result;
    result.sourceMomentumFingerprint = sourceMomentum.fingerprint;
    result.pressureProjectionFingerprint = flowSource.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.pressureVolumeRateFingerprint =
        flowSource.pressureVolumeRateFingerprint;
    result.previousBulkVelocityFingerprint = usesBulkVelocityIncrement
        ? sceneFluidOpeningFluxVelocityFingerprint(
            *bulkGrid, *previousBulkVelocity) : 0;
    result.currentBulkVelocityFingerprint = usesBulkVelocityIncrement
        ? sceneFluidOpeningFluxVelocityFingerprint(
            *bulkGrid, *currentBulkVelocity) : 0;
    result.acceptedStepCount = sourceMomentum.acceptedStepCount;
    result.sourceSimulationTimeSeconds = sourceMomentum.simulationTimeSeconds;
    result.targetSimulationTimeSeconds = sourceMomentum.simulationTimeSeconds
        + settings.timeStepSeconds;
    result.densityKgPerCubicMeter = sourceMomentum.densityKgPerCubicMeter;
    result.settings = settings;
    auto& diagnostics = result.diagnostics;
    diagnostics.controlVolumeCount = sourceMomentum.controlVolumes.size();
    diagnostics.linkCount = faceLinks.links.size();
    diagnostics.usesMovingVolumeRates =
        flowSource.usesMovingVolumeRates;
    diagnostics.usesBulkVelocityIncrement = usesBulkVelocityIncrement;
    std::vector<SceneFluidRegionTransportControlVolume> forcedControls;
    forcedControls.reserve(sourceMomentum.controlVolumes.size());
    for (const auto& source : sourceMomentum.controlVolumes) {
        auto velocity = source.velocityMetersPerSecond;
        if (usesBulkVelocityIncrement) {
            const auto increment = subtract(
                cellCenteredVelocity(
                    *bulkGrid, *currentBulkVelocity, source.cellIndex),
                cellCenteredVelocity(
                    *bulkGrid, *previousBulkVelocity, source.cellIndex));
            velocity = add(velocity, increment);
            diagnostics.maximumBulkVelocityIncrementMetersPerSecond =
                std::max(
                    diagnostics.maximumBulkVelocityIncrementMetersPerSecond,
                    norm(increment));
        }
        const double mass = result.densityKgPerCubicMeter
            * source.volumeCubicMeters;
        forcedControls.push_back({
            source.controlVolumeIndex,
            source.stableId,
            source.volumeCubicMeters,
            velocity,
            scale(velocity, mass),
        });
    }
    diagnostics.momentumBeforeKilogramMetersPerSecond =
        totalMomentum(forcedControls);
    diagnostics.kineticEnergyBeforeJoules = kineticEnergy(
        forcedControls, result.densityKgPerCubicMeter);
    diagnostics.bulkVelocityIncrementImpulseKilogramMetersPerSecond =
        subtract(
            diagnostics.momentumBeforeKilogramMetersPerSecond,
            sourceMomentum.diagnostics.totalMomentumKilogramMetersPerSecond);
    diagnostics.bulkVelocityIncrementWorkJoules =
        diagnostics.kineticEnergyBeforeJoules
        - sourceMomentum.diagnostics.kineticEnergyJoules;

    std::vector<double> correctedFlows(faceLinks.links.size(), 0.0);
    std::vector<double> outwardRates(sourceMomentum.controlVolumes.size(), 0.0);
    std::vector<double> netOutwardRates(
        sourceMomentum.controlVolumes.size(), 0.0);
    std::vector<double> viscousGeometryWeights(
        sourceMomentum.controlVolumes.size(), 0.0);
    std::vector<double> geometryVolumeRates(
        sourceMomentum.controlVolumes.size(), 0.0);
    for (std::size_t index = 0;
         index < sourceMomentum.controlVolumes.size(); ++index) {
        const auto& source = sourceMomentum.controlVolumes[index];
        const auto& projected = flowSource.controls[index];
        if (projected.controlVolumeIndex != index
            || projected.stableId != source.stableId
            || projected.componentIndex != source.componentIndex) {
            throw std::invalid_argument(
                "scene fluid region transport control binding is invalid");
        }
        const double rate = projected
            .geometryVolumeChangeRateCubicMetersPerSecond;
        const double change = settings.timeStepSeconds * rate;
        const double targetVolume = source.volumeCubicMeters + change;
        if (!std::isfinite(rate) || !std::isfinite(targetVolume)) {
            throw std::overflow_error(
                "scene fluid region transport geometry volume is non-finite");
        }
        geometryVolumeRates[index] = rate;
        diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
                std::abs(rate));
        diagnostics.maximumAbsoluteGeometryVolumeChangeCubicMeters =
            std::max(
                diagnostics.maximumAbsoluteGeometryVolumeChangeCubicMeters,
                std::abs(change));
        if (!(targetVolume > 0.0)) {
            diagnostics.failureStage =
                SceneFluidRegionTransportFailureStage::GeometryVolume;
        }
    }
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        const auto& projected = flowSource.links[index];
        if (source.linkIndex != index
            || source.minusControlVolumeIndex >= outwardRates.size()
            || source.plusControlVolumeIndex >= outwardRates.size()
            || projected.linkIndex != index
            || projected.stableId != source.stableId
            || projected.minusControlVolumeIndex
                != source.minusControlVolumeIndex
            || projected.plusControlVolumeIndex
                != source.plusControlVolumeIndex
            || projected.kind != source.kind
            || !(source.areaSquareMeters > 0.0)
            || !(source.centerDistanceMeters > 0.0)) {
            throw std::invalid_argument(
                "scene fluid region transport link binding is invalid");
        }
        const double flow = projected
            .correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        if (!std::isfinite(flow)) {
            throw std::overflow_error(
                "scene fluid region transport flow is non-finite");
        }
        correctedFlows[index] = flow;
        netOutwardRates[source.minusControlVolumeIndex] += flow;
        netOutwardRates[source.plusControlVolumeIndex] -= flow;
        outwardRates[source.minusControlVolumeIndex] += std::max(0.0, flow);
        outwardRates[source.plusControlVolumeIndex] += std::max(0.0, -flow);
        viscousGeometryWeights[source.minusControlVolumeIndex] +=
            source.geometryWeightMeters;
        viscousGeometryWeights[source.plusControlVolumeIndex] +=
            source.geometryWeightMeters;
        if (source.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++diagnostics.openingLinkCount;
        }
    }
    for (std::size_t index = 0; index < outwardRates.size(); ++index) {
        const double volume =
            sourceMomentum.controlVolumes[index].volumeCubicMeters;
        const double targetVolume = volume
            + settings.timeStepSeconds * geometryVolumeRates[index];
        const double stabilityVolume = targetVolume > 0.0
            ? std::min(volume, targetVolume) : volume;
        diagnostics.maximumFullStepOutgoingCourantNumber = std::max(
            diagnostics.maximumFullStepOutgoingCourantNumber,
            settings.timeStepSeconds * outwardRates[index]
                / stabilityVolume);
        diagnostics.maximumCorrectedContinuityResidualCubicMetersPerSecond =
            std::max(
                diagnostics
                    .maximumCorrectedContinuityResidualCubicMetersPerSecond,
                std::abs(
                    geometryVolumeRates[index]
                    + netOutwardRates[index]));
        diagnostics.maximumFullStepViscousNumber = std::max(
            diagnostics.maximumFullStepViscousNumber,
            settings.timeStepSeconds
                * settings.kinematicViscositySquareMetersPerSecond
                * viscousGeometryWeights[index] / stabilityVolume);
    }
    const double continuityTolerance = std::max(
        flowSource.absoluteContinuityToleranceCubicMetersPerSecond,
        flowSource.relativeContinuityTolerance
            * flowSource
                .predictedContinuityResidualMaximumCubicMetersPerSecond);
    if (diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None
        && diagnostics.maximumCorrectedContinuityResidualCubicMetersPerSecond
            > continuityTolerance) {
        diagnostics.failureStage =
            SceneFluidRegionTransportFailureStage::FlowContinuity;
    }
    const std::size_t advectionSubsteps = requiredSubsteps(
        diagnostics.maximumFullStepOutgoingCourantNumber,
        settings.maximumOutgoingCourantNumber);
    const std::size_t viscositySubsteps = requiredSubsteps(
        diagnostics.maximumFullStepViscousNumber,
        settings.maximumViscousNumber);
    diagnostics.substepCount = std::max(advectionSubsteps, viscositySubsteps);
    if (diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None
        && diagnostics.substepCount > settings.maximumSubsteps) {
        diagnostics.failureStage =
            SceneFluidRegionTransportFailureStage::SubstepLimit;
    }
    if (diagnostics.failureStage
        != SceneFluidRegionTransportFailureStage::None) {
        diagnostics.finite = true;
        result.fingerprint = transportFingerprint(result);
        validateSceneFluidRegionTransportIntegrity(result);
        return result;
    }
    diagnostics.maximumAcceptedSubstepOutgoingCourantNumber =
        diagnostics.maximumFullStepOutgoingCourantNumber
        / static_cast<double>(diagnostics.substepCount);
    diagnostics.maximumAcceptedSubstepViscousNumber =
        diagnostics.maximumFullStepViscousNumber
        / static_cast<double>(diagnostics.substepCount);

    std::vector<SceneFluidRegionTransportControlVolume> candidate =
        forcedControls;
    const double substepSeconds = settings.timeStepSeconds
        / static_cast<double>(diagnostics.substepCount);
    double energyAfterAdvection = diagnostics.kineticEnergyBeforeJoules;
    double energyAfterViscosity = diagnostics.kineticEnergyBeforeJoules;
    double advectiveEnergyLoss = 0.0;
    double viscousEnergyLoss = 0.0;
    std::vector<fluid::Vector3> impulse(candidate.size());
    for (std::size_t substep = 0;
         substep < diagnostics.substepCount; ++substep) {
        std::ranges::fill(impulse, fluid::Vector3{});
        for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
            const auto& link = faceLinks.links[index];
            const double flow = correctedFlows[index];
            const auto& donor = flow >= 0.0
                ? candidate[link.minusControlVolumeIndex]
                      .velocityMetersPerSecond
                : candidate[link.plusControlVolumeIndex]
                      .velocityMetersPerSecond;
            const auto transported = scale(
                donor, result.densityKgPerCubicMeter
                    * flow * substepSeconds);
            impulse[link.minusControlVolumeIndex] = subtract(
                impulse[link.minusControlVolumeIndex], transported);
            impulse[link.plusControlVolumeIndex] = add(
                impulse[link.plusControlVolumeIndex], transported);
        }
        for (std::size_t index = 0; index < candidate.size(); ++index) {
            candidate[index].momentumKilogramMetersPerSecond = add(
                candidate[index].momentumKilogramMetersPerSecond,
                impulse[index]);
            const double elapsedSeconds =
                substep + 1 == diagnostics.substepCount
                ? settings.timeStepSeconds
                : substepSeconds * static_cast<double>(substep + 1);
            candidate[index].volumeCubicMeters =
                sourceMomentum.controlVolumes[index].volumeCubicMeters
                + geometryVolumeRates[index] * elapsedSeconds;
        }
        updateVelocities(candidate, result.densityKgPerCubicMeter);
        energyAfterAdvection = kineticEnergy(
            candidate, result.densityKgPerCubicMeter);
        const double advectionReference = std::max(
            diagnostics.kineticEnergyBeforeJoules, energyAfterViscosity);
        if (!std::isfinite(energyAfterAdvection)
            || energyAfterAdvection > energyAfterViscosity
                + tolerance(
                    settings.absoluteEnergyToleranceJoules,
                    settings.relativeEnergyTolerance,
                    advectionReference)) {
            diagnostics.failureStage =
                SceneFluidRegionTransportFailureStage::AdvectionEnergy;
            break;
        }
        advectiveEnergyLoss +=
            energyAfterViscosity - energyAfterAdvection;

        std::ranges::fill(impulse, fluid::Vector3{});
        for (const auto& link : faceLinks.links) {
            const auto difference = subtract(
                candidate[link.plusControlVolumeIndex]
                    .velocityMetersPerSecond,
                candidate[link.minusControlVolumeIndex]
                    .velocityMetersPerSecond);
            const double coefficient = result.densityKgPerCubicMeter
                * settings.kinematicViscositySquareMetersPerSecond
                * link.geometryWeightMeters * substepSeconds;
            const auto exchanged = scale(difference, coefficient);
            impulse[link.minusControlVolumeIndex] = add(
                impulse[link.minusControlVolumeIndex], exchanged);
            impulse[link.plusControlVolumeIndex] = subtract(
                impulse[link.plusControlVolumeIndex], exchanged);
        }
        for (std::size_t index = 0; index < candidate.size(); ++index) {
            candidate[index].momentumKilogramMetersPerSecond = add(
                candidate[index].momentumKilogramMetersPerSecond,
                impulse[index]);
        }
        updateVelocities(candidate, result.densityKgPerCubicMeter);
        energyAfterViscosity = kineticEnergy(
            candidate, result.densityKgPerCubicMeter);
        if (!std::isfinite(energyAfterViscosity)
            || energyAfterViscosity > energyAfterAdvection
                + tolerance(
                    settings.absoluteEnergyToleranceJoules,
                    settings.relativeEnergyTolerance,
                    energyAfterAdvection)) {
            diagnostics.failureStage =
                SceneFluidRegionTransportFailureStage::ViscosityEnergy;
            break;
        }
        viscousEnergyLoss +=
            energyAfterAdvection - energyAfterViscosity;
    }
    diagnostics.kineticEnergyAfterAdvectionJoules = energyAfterAdvection;
    diagnostics.kineticEnergyAfterJoules = energyAfterViscosity;
    diagnostics.advectiveEnergyLossJoules = advectiveEnergyLoss;
    diagnostics.viscousEnergyLossJoules = viscousEnergyLoss;
    diagnostics.momentumAfterKilogramMetersPerSecond =
        totalMomentum(candidate);
    diagnostics.momentumResidualKilogramMetersPerSecond = subtract(
        diagnostics.momentumAfterKilogramMetersPerSecond,
        diagnostics.momentumBeforeKilogramMetersPerSecond);
    diagnostics.momentumResidualNormKilogramMetersPerSecond = norm(
        diagnostics.momentumResidualKilogramMetersPerSecond);
    const double momentumReference = std::max(
        norm(diagnostics.momentumBeforeKilogramMetersPerSecond),
        norm(diagnostics.momentumAfterKilogramMetersPerSecond));
    if (diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None
        && diagnostics.momentumResidualNormKilogramMetersPerSecond
            > tolerance(
                settings.absoluteMomentumToleranceKilogramMetersPerSecond,
                settings.relativeMomentumTolerance, momentumReference)) {
        diagnostics.failureStage =
            SceneFluidRegionTransportFailureStage::Conservation;
    }
    for (std::size_t index = 0; index < candidate.size(); ++index) {
        diagnostics.maximumVelocityChangeMetersPerSecond = std::max(
            diagnostics.maximumVelocityChangeMetersPerSecond,
            norm(subtract(
                candidate[index].velocityMetersPerSecond,
                sourceMomentum.controlVolumes[index]
                    .velocityMetersPerSecond)));
    }
    diagnostics.finite = finite(
            diagnostics.momentumBeforeKilogramMetersPerSecond)
        && finite(diagnostics
            .bulkVelocityIncrementImpulseKilogramMetersPerSecond)
        && finite(diagnostics.momentumAfterKilogramMetersPerSecond)
        && finite(diagnostics.momentumResidualKilogramMetersPerSecond)
        && std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterAdvectionJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.advectiveEnergyLossJoules)
        && std::isfinite(diagnostics.viscousEnergyLossJoules)
        && std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumBulkVelocityIncrementMetersPerSecond)
        && std::isfinite(diagnostics.bulkVelocityIncrementWorkJoules)
        && std::ranges::all_of(candidate, [](const auto& control) {
            return finite(control.velocityMetersPerSecond)
                && finite(control.momentumKilogramMetersPerSecond);
        });
    if (!diagnostics.finite
        && diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None) {
        diagnostics.failureStage =
            SceneFluidRegionTransportFailureStage::NonFinite;
    }
    if (diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None
        && diagnostics.finite) {
        diagnostics.accepted = true;
        result.ownedStorageBytes = storageBytes;
        result.controlVolumes = std::move(candidate);
    }
    result.fingerprint = transportFingerprint(result);
    validateSceneFluidRegionTransportIntegrity(result);
    return result;
}

SceneFluidRegionTransport advanceSceneFluidRegionMomentum(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureProjection& correctedProjection,
    const SceneFluidRegionTransportSettings& settings,
    const SceneFluidRegionTransportLimits& limits) {
    const auto flowSource = regionTransportFlowSource(
        correctedProjection);
    return advanceSceneFluidRegionMomentumImpl(
        sourceMomentum, faceLinks, flowSource,
        nullptr, nullptr, nullptr, settings, limits);
}

SceneFluidRegionTransport advanceSceneFluidRegionMomentum(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const SceneFluidRegionTransportSettings& settings,
    const SceneFluidRegionTransportLimits& limits) {
    const auto flowSource = regionTransportFlowSource(
        sourceMomentum, faceLinks, correctedFlow);
    auto result = advanceSceneFluidRegionMomentumImpl(
        sourceMomentum, faceLinks, flowSource,
        nullptr, nullptr, nullptr, settings, limits);
    validateSceneFluidRegionTransport(
        result, sourceMomentum, faceLinks, correctedFlow);
    return result;
}

SceneFluidRegionTransport advanceSceneFluidRegionMomentum(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureProjection& correctedProjection,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& previousBulkVelocityMetersPerSecond,
    const fluid::MacVelocityField& currentBulkVelocityMetersPerSecond,
    const SceneFluidRegionTransportSettings& settings,
    const SceneFluidRegionTransportLimits& limits) {
    const auto flowSource = regionTransportFlowSource(
        correctedProjection);
    return advanceSceneFluidRegionMomentumImpl(
        sourceMomentum, faceLinks, flowSource,
        &grid, &previousBulkVelocityMetersPerSecond,
        &currentBulkVelocityMetersPerSecond, settings, limits);
}

SceneFluidRegionTransport advanceSceneFluidRegionMomentum(
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& previousBulkVelocityMetersPerSecond,
    const fluid::MacVelocityField& currentBulkVelocityMetersPerSecond,
    const SceneFluidRegionTransportSettings& settings,
    const SceneFluidRegionTransportLimits& limits) {
    const auto flowSource = regionTransportFlowSource(
        sourceMomentum, faceLinks, correctedFlow);
    auto result = advanceSceneFluidRegionMomentumImpl(
        sourceMomentum, faceLinks, flowSource,
        &grid, &previousBulkVelocityMetersPerSecond,
        &currentBulkVelocityMetersPerSecond, settings, limits);
    validateSceneFluidRegionTransport(
        result, sourceMomentum, faceLinks, correctedFlow);
    return result;
}

void validateSceneFluidRegionTransportIntegrity(
    const SceneFluidRegionTransport& transport) {
    const auto& diagnostics = transport.diagnostics;
    bool controlsValid = true;
    for (std::size_t index = 0;
         index < transport.controlVolumes.size(); ++index) {
        const auto& control = transport.controlVolumes[index];
        controlsValid = controlsValid
            && control.controlVolumeIndex == index
            && control.stableId != 0
            && std::isfinite(control.volumeCubicMeters)
            && control.volumeCubicMeters > 0.0
            && finite(control.velocityMetersPerSecond)
            && finite(control.momentumKilogramMetersPerSecond);
        const auto reconstructedMomentum = scale(
            control.velocityMetersPerSecond,
            transport.densityKgPerCubicMeter
                * control.volumeCubicMeters);
        controlsValid = controlsValid
            && norm(subtract(
                   control.momentumKilogramMetersPerSecond,
                   reconstructedMomentum))
                <= tolerance(
                    transport.settings
                        .absoluteMomentumToleranceKilogramMetersPerSecond,
                    transport.settings.relativeMomentumTolerance,
                    norm(control.momentumKilogramMetersPerSecond));
    }
    const bool acceptedShape = diagnostics.accepted
        && diagnostics.failureStage
            == SceneFluidRegionTransportFailureStage::None
        && transport.controlVolumes.size()
            == diagnostics.controlVolumeCount
        && transport.ownedStorageBytes
            == storageBytesForControlVolumes(
                transport.controlVolumes.size());
    const bool rejectedShape = !diagnostics.accepted
        && diagnostics.failureStage
            != SceneFluidRegionTransportFailureStage::None
        && transport.controlVolumes.empty()
        && transport.ownedStorageBytes == 0;
    if (transport.version != sceneFluidRegionTransportVersion
        || transport.fingerprint == 0
        || transport.sourceMomentumFingerprint == 0
        || transport.pressureProjectionFingerprint == 0
        || transport.pressureFaceLinkFingerprint == 0
        || diagnostics.usesMovingVolumeRates
            != (transport.pressureVolumeRateFingerprint != 0)
        || diagnostics.usesBulkVelocityIncrement
            != (transport.previousBulkVelocityFingerprint != 0)
        || diagnostics.usesBulkVelocityIncrement
            != (transport.currentBulkVelocityFingerprint != 0)
        || !std::isfinite(transport.sourceSimulationTimeSeconds)
        || !std::isfinite(transport.targetSimulationTimeSeconds)
        || transport.targetSimulationTimeSeconds
            != transport.sourceSimulationTimeSeconds
                + transport.settings.timeStepSeconds
        || !std::isfinite(transport.densityKgPerCubicMeter)
        || !(transport.densityKgPerCubicMeter > 0.0)
        || !validSettings(transport.settings)
        || diagnostics.controlVolumeCount == 0
        || diagnostics.linkCount == 0
        || diagnostics.openingLinkCount > diagnostics.linkCount
        || diagnostics.substepCount == 0
        || !std::isfinite(diagnostics
            .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond)
        || diagnostics
               .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond < 0.0
        || !std::isfinite(diagnostics
            .maximumAbsoluteGeometryVolumeChangeCubicMeters)
        || diagnostics.maximumAbsoluteGeometryVolumeChangeCubicMeters < 0.0
        || !std::isfinite(
            diagnostics.maximumBulkVelocityIncrementMetersPerSecond)
        || diagnostics.maximumBulkVelocityIncrementMetersPerSecond < 0.0
        || !finite(diagnostics
            .bulkVelocityIncrementImpulseKilogramMetersPerSecond)
        || !std::isfinite(diagnostics.bulkVelocityIncrementWorkJoules)
        || !std::isfinite(diagnostics
            .maximumCorrectedContinuityResidualCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumFullStepOutgoingCourantNumber)
        || !std::isfinite(
            diagnostics.maximumAcceptedSubstepOutgoingCourantNumber)
        || !std::isfinite(diagnostics.maximumFullStepViscousNumber)
        || !std::isfinite(
            diagnostics.maximumAcceptedSubstepViscousNumber)
        || !finite(diagnostics.momentumBeforeKilogramMetersPerSecond)
        || !finite(diagnostics.momentumAfterKilogramMetersPerSecond)
        || !finite(diagnostics.momentumResidualKilogramMetersPerSecond)
        || !std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        || !std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        || !std::isfinite(diagnostics.kineticEnergyAfterAdvectionJoules)
        || !std::isfinite(diagnostics.kineticEnergyAfterJoules)
        || !std::isfinite(diagnostics.advectiveEnergyLossJoules)
        || !std::isfinite(diagnostics.viscousEnergyLossJoules)
        || !std::isfinite(
            diagnostics.maximumVelocityChangeMetersPerSecond)
        || !diagnostics.finite
        || (!acceptedShape && !rejectedShape)
        || !controlsValid
        || transport.fingerprint != transportFingerprint(transport)) {
        throw std::invalid_argument(
            "scene fluid region transport integrity is invalid");
    }
}

void validateSceneFluidRegionTransport(
    const SceneFluidRegionTransport& transport,
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureProjection& correctedProjection) {
    validateSceneFluidRegionTransportIntegrity(transport);
    validateSceneFluidRegionMomentumStateIntegrity(sourceMomentum);
    validateSceneFluidPressureProjectionIntegrity(correctedProjection);
    if (transport.sourceMomentumFingerprint != sourceMomentum.fingerprint
        || transport.pressureProjectionFingerprint
            != correctedProjection.fingerprint
        || transport.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || transport.pressureVolumeRateFingerprint
            != correctedProjection.pressureVolumeRateFingerprint
        || transport.acceptedStepCount != sourceMomentum.acceptedStepCount
        || transport.sourceSimulationTimeSeconds
            != sourceMomentum.simulationTimeSeconds
        || transport.densityKgPerCubicMeter
            != sourceMomentum.densityKgPerCubicMeter) {
        throw std::invalid_argument(
            "scene fluid region transport binding is invalid");
    }
    if (transport.diagnostics.accepted) {
        for (std::size_t index = 0;
             index < transport.controlVolumes.size(); ++index) {
            const auto& control = transport.controlVolumes[index];
            const auto& source = sourceMomentum.controlVolumes[index];
            const double expectedVolume = source.volumeCubicMeters
                + transport.settings.timeStepSeconds
                    * correctedProjection.controlVolumes[index]
                        .geometryVolumeChangeRateCubicMetersPerSecond;
            if (control.controlVolumeIndex != source.controlVolumeIndex
                || control.stableId != source.stableId
                || control.volumeCubicMeters != expectedVolume) {
                throw std::invalid_argument(
                    "scene fluid region transport control binding is invalid");
            }
        }
    }
}

void validateSceneFluidRegionTransport(
    const SceneFluidRegionTransport& transport,
    const SceneFluidRegionMomentumState& sourceMomentum,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow) {
    validateSceneFluidRegionTransportIntegrity(transport);
    validateSceneFluidRegionMomentumStateIntegrity(sourceMomentum);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    validateSceneFluidMimeticCorrectedTraceFlowIntegrity(correctedFlow);
    if (!correctedFlow.accepted
        || transport.sourceMomentumFingerprint != sourceMomentum.fingerprint
        || transport.pressureProjectionFingerprint
            != correctedFlow.fingerprint
        || transport.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || transport.pressureVolumeRateFingerprint != 0
        || sourceMomentum.pressureProjectionFingerprint
            != correctedFlow.fingerprint
        || correctedFlow.pressureFaceLinkFingerprint
            != faceLinks.fingerprint
        || transport.acceptedStepCount != sourceMomentum.acceptedStepCount
        || transport.sourceSimulationTimeSeconds
            != sourceMomentum.simulationTimeSeconds
        || transport.densityKgPerCubicMeter
            != sourceMomentum.densityKgPerCubicMeter) {
        throw std::invalid_argument(
            "scene fluid mimetic region transport binding is invalid");
    }
    if (transport.diagnostics.accepted) {
        for (std::size_t index = 0;
             index < transport.controlVolumes.size(); ++index) {
            const auto& control = transport.controlVolumes[index];
            const auto& source = sourceMomentum.controlVolumes[index];
            if (control.controlVolumeIndex != source.controlVolumeIndex
                || control.stableId != source.stableId
                || control.volumeCubicMeters != source.volumeCubicMeters) {
                throw std::invalid_argument(
                    "scene fluid mimetic region transport control binding is invalid");
            }
        }
    }
}

} // namespace simwing::fsi
