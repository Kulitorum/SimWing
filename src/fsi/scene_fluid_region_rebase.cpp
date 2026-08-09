#include "scene_fluid_region_rebase.h"

#include <algorithm>
#include <bit>
#include <cmath>
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

fluid::Vector3 add(const fluid::Vector3& first,
                   const fluid::Vector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

fluid::Vector3 scale(const fluid::Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double norm(const fluid::Vector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

std::size_t storageBytesForControls(const std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidRegionRebaseControlVolume)) {
        throw std::length_error(
            "scene fluid region rebase storage size overflows");
    }
    return count * sizeof(SceneFluidRegionRebaseControlVolume);
}

std::uint64_t rebaseFingerprint(const SceneFluidRegionRebase& rebase) {
    Fingerprint fingerprint;
    fingerprint.integer(rebase.version);
    for (const std::uint64_t value : {
             rebase.sourceTransportFingerprint,
             rebase.sourceTopologyTransitionFingerprint,
             rebase.previousPressureControlVolumeFingerprint,
             rebase.currentPressureControlVolumeFingerprint,
             rebase.surfaceDefinitionFingerprint,
             rebase.structureDefinitionFingerprint,
             rebase.previousAcceptedStepCount,
             rebase.currentAcceptedStepCount}) {
        fingerprint.integer(value);
    }
    for (const double value : {
             rebase.previousSimulationTimeSeconds,
             rebase.currentSimulationTimeSeconds,
             rebase.densityKgPerCubicMeter,
             rebase.lowerMeters.x,
             rebase.lowerMeters.y,
             rebase.lowerMeters.z,
             rebase.upperMeters.x,
             rebase.upperMeters.y,
             rebase.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(rebase.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(rebase.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(rebase.cellCounts.z));
    fingerprint.integer(static_cast<std::uint64_t>(
        rebase.ownedStorageBytes));
    const auto& diagnostics = rebase.diagnostics;
    for (const std::size_t value : {
             diagnostics.previousControlVolumeCount,
             diagnostics.currentControlVolumeCount,
             diagnostics.retainedControlVolumeCount,
             diagnostics.appearedControlVolumeCount,
             diagnostics.disappearedControlVolumeCount,
             diagnostics.maximumDonorControlVolumeCount,
             diagnostics.maximumRetirementRecipientCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    for (const double value : {
             diagnostics.sourceVolumeCubicMeters,
             diagnostics.mappedSourceVolumeCubicMeters,
             diagnostics.sourceVolumeMappingResidualCubicMeters,
             diagnostics.sourceMomentumKilogramMetersPerSecond.x,
             diagnostics.sourceMomentumKilogramMetersPerSecond.y,
             diagnostics.sourceMomentumKilogramMetersPerSecond.z,
             diagnostics.mappedSourceMomentumKilogramMetersPerSecond.x,
             diagnostics.mappedSourceMomentumKilogramMetersPerSecond.y,
             diagnostics.mappedSourceMomentumKilogramMetersPerSecond.z,
             diagnostics.sourceMomentumMappingResidualKilogramMetersPerSecond.x,
             diagnostics.sourceMomentumMappingResidualKilogramMetersPerSecond.y,
             diagnostics.sourceMomentumMappingResidualKilogramMetersPerSecond.z,
             diagnostics.sourceMomentumMappingResidualNormKilogramMetersPerSecond,
             diagnostics.rebasedMomentumKilogramMetersPerSecond.x,
             diagnostics.rebasedMomentumKilogramMetersPerSecond.y,
             diagnostics.rebasedMomentumKilogramMetersPerSecond.z,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.x,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.y,
             diagnostics.geometricMomentumChangeKilogramMetersPerSecond.z,
             diagnostics.sourceKineticEnergyJoules,
             diagnostics.rebasedKineticEnergyJoules,
             diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
             diagnostics.maximumAppearedVolumeCubicMeters,
             diagnostics.maximumAbsoluteVelocityMetersPerSecond}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.finite));
    fingerprint.integer(static_cast<std::uint64_t>(
        rebase.controlVolumes.size()));
    for (const auto& control : rebase.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(control.cellIndex));
        fingerprint.integer(control.regionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.integer(static_cast<std::uint8_t>(
            control.appearedThisEpoch));
        fingerprint.integer(static_cast<std::uint64_t>(
            control.sourceControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            control.donorControlVolumeCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            control.retiredSourceControlVolumeCount));
        for (const double value : {
                 control.donorLinkAreaSquareMeters,
                 control.retiredSourceLinkAreaSquareMeters,
                 control.mappedSourceVolumeCubicMeters,
                 control.mappedSourceMomentumKilogramMetersPerSecond.x,
                 control.mappedSourceMomentumKilogramMetersPerSecond.y,
                 control.mappedSourceMomentumKilogramMetersPerSecond.z,
                 control.volumeCubicMeters,
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

void validateSources(
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureTopologyTransition& topologyTransition) {
    validateSceneFluidRegionTransportIntegrity(transport);
    validateSceneFluidPressureControlVolumeIntegrity(previousPressureVolumes);
    validateSceneFluidPressureControlVolumeIntegrity(currentPressureVolumes);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        topologyTransition);
    if (!transport.diagnostics.accepted
        || previousPressureVolumes.surfaceDefinitionFingerprint == 0
        || previousPressureVolumes.surfaceDefinitionFingerprint
            != currentPressureVolumes.surfaceDefinitionFingerprint
        || previousPressureVolumes.structureDefinitionFingerprint
            != currentPressureVolumes.structureDefinitionFingerprint
        || previousPressureVolumes.acceptedStepCount
            != transport.acceptedStepCount
        || previousPressureVolumes.simulationTimeSeconds
            != transport.sourceSimulationTimeSeconds
        || previousPressureVolumes.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentPressureVolumes.acceptedStepCount
            != previousPressureVolumes.acceptedStepCount + 1
        || currentPressureVolumes.simulationTimeSeconds
            != transport.targetSimulationTimeSeconds
        || previousPressureVolumes.cellCounts
            != currentPressureVolumes.cellCounts
        || previousPressureVolumes.lowerMeters
            != currentPressureVolumes.lowerMeters
        || previousPressureVolumes.upperMeters
            != currentPressureVolumes.upperMeters
        || !(currentPressureVolumes.simulationTimeSeconds
            > previousPressureVolumes.simulationTimeSeconds)
        || previousPressureVolumes.controlVolumes.size()
            != transport.controlVolumes.size()
        || topologyTransition.previousPressureControlVolumeFingerprint
            != previousPressureVolumes.fingerprint
        || topologyTransition.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || topologyTransition.previousAcceptedStepCount
            != previousPressureVolumes.acceptedStepCount
        || topologyTransition.currentAcceptedStepCount
            != currentPressureVolumes.acceptedStepCount
        || topologyTransition.previousSimulationTimeSeconds
            != previousPressureVolumes.simulationTimeSeconds
        || topologyTransition.currentSimulationTimeSeconds
            != currentPressureVolumes.simulationTimeSeconds
        || topologyTransition.cellCounts != currentPressureVolumes.cellCounts
        || topologyTransition.lowerMeters != currentPressureVolumes.lowerMeters
        || topologyTransition.upperMeters != currentPressureVolumes.upperMeters) {
        throw std::invalid_argument(
            "scene fluid region rebase source identity is invalid");
    }
    for (std::size_t index = 0;
         index < previousPressureVolumes.controlVolumes.size(); ++index) {
        const auto& previous = previousPressureVolumes.controlVolumes[index];
        const auto& source = transport.controlVolumes[index];
        if (previous.controlVolumeIndex != index
            || source.controlVolumeIndex != index
            || source.stableId != previous.stableId) {
            throw std::invalid_argument(
                "scene fluid region rebase source binding is invalid");
        }
    }
}

SceneFluidRegionRebase buildRebase(
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidRegionRebaseLimits& limits) {
    if (currentPressureVolumes.controlVolumes.size()
            > limits.maximumControlVolumes
        || topologyTransition.appearanceDonors.size()
            > limits.maximumMappings
        || topologyTransition.retirementRecipients.size()
            > limits.maximumMappings
                - topologyTransition.appearanceDonors.size()) {
        throw std::length_error(
            "scene fluid region rebase exceeds its count limit");
    }
    const std::size_t storageBytes = storageBytesForControls(
        currentPressureVolumes.controlVolumes.size());
    if (storageBytes > limits.maximumRebaseBytes) {
        throw std::length_error(
            "scene fluid region rebase exceeds its byte limit");
    }

    SceneFluidRegionRebase result;
    result.sourceTransportFingerprint = transport.fingerprint;
    result.sourceTopologyTransitionFingerprint =
        topologyTransition.fingerprint;
    result.previousPressureControlVolumeFingerprint =
        previousPressureVolumes.fingerprint;
    result.currentPressureControlVolumeFingerprint =
        currentPressureVolumes.fingerprint;
    result.surfaceDefinitionFingerprint =
        currentPressureVolumes.surfaceDefinitionFingerprint;
    result.structureDefinitionFingerprint =
        currentPressureVolumes.structureDefinitionFingerprint;
    result.previousAcceptedStepCount =
        previousPressureVolumes.acceptedStepCount;
    result.currentAcceptedStepCount = currentPressureVolumes.acceptedStepCount;
    result.previousSimulationTimeSeconds =
        previousPressureVolumes.simulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        currentPressureVolumes.simulationTimeSeconds;
    result.densityKgPerCubicMeter = transport.densityKgPerCubicMeter;
    result.cellCounts = currentPressureVolumes.cellCounts;
    result.lowerMeters = currentPressureVolumes.lowerMeters;
    result.upperMeters = currentPressureVolumes.upperMeters;
    result.ownedStorageBytes = storageBytes;
    auto& diagnostics = result.diagnostics;
    diagnostics.previousControlVolumeCount =
        previousPressureVolumes.controlVolumes.size();
    diagnostics.currentControlVolumeCount =
        currentPressureVolumes.controlVolumes.size();
    diagnostics.sourceMomentumKilogramMetersPerSecond =
        transport.diagnostics.momentumAfterKilogramMetersPerSecond;
    diagnostics.sourceKineticEnergyJoules =
        transport.diagnostics.kineticEnergyAfterJoules;

    std::vector<std::size_t> sourceIndexByCurrent(
        currentPressureVolumes.controlVolumes.size(),
        std::numeric_limits<std::size_t>::max());
    for (const auto& retained : topologyTransition.retainedControls) {
        sourceIndexByCurrent[retained.currentControlVolumeIndex] =
            retained.previousControlVolumeIndex;
    }
    struct WorkingControl {
        SceneFluidRegionRebaseControlVolume control;
        fluid::Vector3 geometricSeedVelocityMetersPerSecond;
    };
    std::vector<WorkingControl> working;
    working.reserve(currentPressureVolumes.controlVolumes.size());
    std::size_t appearanceDonorCursor = 0;
    for (const auto& current : currentPressureVolumes.controlVolumes) {
        WorkingControl entry;
        auto& control = entry.control;
        control.controlVolumeIndex = current.controlVolumeIndex;
        control.stableId = current.stableId;
        control.cellIndex = current.cellIndex;
        control.regionId = current.regionId;
        control.componentIndex = current.componentIndex;
        control.volumeCubicMeters = current.volumeCubicMeters;
        const std::size_t retained =
            sourceIndexByCurrent[current.controlVolumeIndex];
        control.appearedThisEpoch =
            retained == std::numeric_limits<std::size_t>::max();
        if (!control.appearedThisEpoch) {
            control.sourceControlVolumeIndex = retained;
            control.donorControlVolumeCount = 1;
            entry.geometricSeedVelocityMetersPerSecond =
                transport.controlVolumes[retained]
                    .velocityMetersPerSecond;
            control.mappedSourceVolumeCubicMeters =
                transport.controlVolumes[retained]
                    .volumeCubicMeters;
            control.mappedSourceMomentumKilogramMetersPerSecond =
                transport.controlVolumes[retained]
                    .momentumKilogramMetersPerSecond;
            ++diagnostics.retainedControlVolumeCount;
        } else {
            fluid::Vector3 weightedVelocity;
            while (appearanceDonorCursor
                    < topologyTransition.appearanceDonors.size()
                && topologyTransition.appearanceDonors[
                       appearanceDonorCursor]
                       .appearedCurrentControlVolumeIndex
                    == current.controlVolumeIndex) {
                const auto& donor = topologyTransition.appearanceDonors[
                    appearanceDonorCursor++];
                weightedVelocity = add(
                    weightedVelocity,
                    scale(
                        transport.controlVolumes[
                            donor.retainedPreviousControlVolumeIndex]
                            .velocityMetersPerSecond,
                        donor.linkAreaSquareMeters));
                control.donorLinkAreaSquareMeters +=
                    donor.linkAreaSquareMeters;
                ++control.donorControlVolumeCount;
            }
            if (control.donorControlVolumeCount == 0
                || !(control.donorLinkAreaSquareMeters > 0.0)) {
                throw std::invalid_argument(
                    "scene fluid region rebase appeared control lacks a "
                    "retained same-region donor");
            }
            entry.geometricSeedVelocityMetersPerSecond = scale(
                weightedVelocity,
                1.0 / control.donorLinkAreaSquareMeters);
            ++diagnostics.appearedControlVolumeCount;
            diagnostics.maximumAppearedVolumeCubicMeters = std::max(
                diagnostics.maximumAppearedVolumeCubicMeters,
                current.volumeCubicMeters);
            diagnostics.maximumDonorControlVolumeCount = std::max(
                diagnostics.maximumDonorControlVolumeCount,
                control.donorControlVolumeCount);
        }
        working.push_back(entry);
    }

    diagnostics.disappearedControlVolumeCount =
        topologyTransition.disappearedControlVolumeCount;
    diagnostics.maximumRetirementRecipientCount =
        topologyTransition.maximumRetirementRecipientCount;
    for (const auto& recipient :
         topologyTransition.retirementRecipients) {
        const auto& source =
            transport.controlVolumes[
                recipient.disappearedPreviousControlVolumeIndex];
        auto& control = working[
            recipient.retainedCurrentControlVolumeIndex].control;
        control.mappedSourceVolumeCubicMeters +=
            recipient.normalizedWeight * source.volumeCubicMeters;
        control.mappedSourceMomentumKilogramMetersPerSecond = add(
            control.mappedSourceMomentumKilogramMetersPerSecond,
            scale(
                source.momentumKilogramMetersPerSecond,
                recipient.normalizedWeight));
        ++control.retiredSourceControlVolumeCount;
        control.retiredSourceLinkAreaSquareMeters +=
            recipient.linkAreaSquareMeters;
    }

    result.controlVolumes.reserve(working.size());
    for (auto& entry : working) {
        auto& control = entry.control;
        const double geometricVolumeChange = control.volumeCubicMeters
            - control.mappedSourceVolumeCubicMeters;
        const double mass = result.densityKgPerCubicMeter
            * control.volumeCubicMeters;
        if (!control.appearedThisEpoch
            && control.retiredSourceControlVolumeCount == 0) {
            control.velocityMetersPerSecond =
                entry.geometricSeedVelocityMetersPerSecond;
            control.momentumKilogramMetersPerSecond = scale(
                control.velocityMetersPerSecond, mass);
        } else {
            control.momentumKilogramMetersPerSecond = add(
                control.mappedSourceMomentumKilogramMetersPerSecond,
                scale(
                    entry.geometricSeedVelocityMetersPerSecond,
                    result.densityKgPerCubicMeter
                        * geometricVolumeChange));
            control.velocityMetersPerSecond = scale(
                control.momentumKilogramMetersPerSecond, 1.0 / mass);
        }
        diagnostics.maximumAbsoluteVolumeChangeCubicMeters = std::max(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters,
            std::abs(geometricVolumeChange));
        diagnostics.mappedSourceVolumeCubicMeters +=
            control.mappedSourceVolumeCubicMeters;
        diagnostics.mappedSourceMomentumKilogramMetersPerSecond = add(
            diagnostics.mappedSourceMomentumKilogramMetersPerSecond,
            control.mappedSourceMomentumKilogramMetersPerSecond);
        diagnostics.rebasedMomentumKilogramMetersPerSecond = add(
            diagnostics.rebasedMomentumKilogramMetersPerSecond,
            control.momentumKilogramMetersPerSecond);
        diagnostics.rebasedKineticEnergyJoules += 0.5 * mass
            * (control.velocityMetersPerSecond.x
                    * control.velocityMetersPerSecond.x
               + control.velocityMetersPerSecond.y
                    * control.velocityMetersPerSecond.y
               + control.velocityMetersPerSecond.z
                    * control.velocityMetersPerSecond.z);
        diagnostics.maximumAbsoluteVelocityMetersPerSecond = std::max(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond,
            norm(control.velocityMetersPerSecond));
        result.controlVolumes.push_back(control);
    }
    diagnostics.sourceVolumeCubicMeters = 0.0;
    for (const auto& source : transport.controlVolumes) {
        diagnostics.sourceVolumeCubicMeters += source.volumeCubicMeters;
    }
    diagnostics.sourceVolumeMappingResidualCubicMeters =
        diagnostics.mappedSourceVolumeCubicMeters
        - diagnostics.sourceVolumeCubicMeters;
    diagnostics.sourceMomentumMappingResidualKilogramMetersPerSecond =
        subtract(
            diagnostics.mappedSourceMomentumKilogramMetersPerSecond,
            diagnostics.sourceMomentumKilogramMetersPerSecond);
    diagnostics.sourceMomentumMappingResidualNormKilogramMetersPerSecond =
        norm(diagnostics
            .sourceMomentumMappingResidualKilogramMetersPerSecond);
    diagnostics.geometricMomentumChangeKilogramMetersPerSecond = subtract(
        diagnostics.rebasedMomentumKilogramMetersPerSecond,
        diagnostics.sourceMomentumKilogramMetersPerSecond);
    const double volumeMappingTolerance = 1.0e-12
        + 64.0 * std::numeric_limits<double>::epsilon()
            * std::max(
                diagnostics.sourceVolumeCubicMeters,
                diagnostics.mappedSourceVolumeCubicMeters);
    const double momentumMappingTolerance =
        transport.settings.absoluteMomentumToleranceKilogramMetersPerSecond
        + transport.settings.relativeMomentumTolerance
            * std::max(
                norm(diagnostics.sourceMomentumKilogramMetersPerSecond),
                norm(diagnostics
                    .mappedSourceMomentumKilogramMetersPerSecond));
    diagnostics.finite = finite(
            diagnostics.sourceMomentumKilogramMetersPerSecond)
        && finite(diagnostics.mappedSourceMomentumKilogramMetersPerSecond)
        && finite(diagnostics
            .sourceMomentumMappingResidualKilogramMetersPerSecond)
        && finite(diagnostics.rebasedMomentumKilogramMetersPerSecond)
        && finite(diagnostics.geometricMomentumChangeKilogramMetersPerSecond)
        && std::isfinite(diagnostics.sourceVolumeCubicMeters)
        && std::isfinite(diagnostics.mappedSourceVolumeCubicMeters)
        && std::isfinite(
            diagnostics.sourceVolumeMappingResidualCubicMeters)
        && std::isfinite(diagnostics
            .sourceMomentumMappingResidualNormKilogramMetersPerSecond)
        && std::isfinite(diagnostics.sourceKineticEnergyJoules)
        && std::isfinite(diagnostics.rebasedKineticEnergyJoules)
        && std::isfinite(
            diagnostics.maximumAbsoluteVolumeChangeCubicMeters)
        && std::isfinite(diagnostics.maximumAppearedVolumeCubicMeters)
        && std::isfinite(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene fluid region rebase is non-finite");
    }
    if (std::abs(diagnostics.sourceVolumeMappingResidualCubicMeters)
            > volumeMappingTolerance
        || diagnostics.sourceMomentumMappingResidualNormKilogramMetersPerSecond
            > momentumMappingTolerance) {
        throw std::runtime_error(
            "scene fluid region rebase source mapping is not conservative");
    }
    result.fingerprint = rebaseFingerprint(result);
    return result;
}

} // namespace

SceneFluidRegionRebase rebaseSceneFluidRegionTransport(
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const SceneFluidRegionRebaseLimits& limits) {
    validateSources(
        transport, previousPressureVolumes, currentPressureVolumes,
        topologyTransition);
    auto result = buildRebase(
        transport, previousPressureVolumes, currentPressureVolumes,
        topologyTransition, limits);
    validateSceneFluidRegionRebase(
        result, transport, previousPressureVolumes,
        currentPressureVolumes, topologyTransition);
    return result;
}

void validateSceneFluidRegionRebaseIntegrity(
    const SceneFluidRegionRebase& rebase) {
    if (rebase.version != sceneFluidRegionRebaseVersion
        || rebase.fingerprint == 0
        || !(rebase.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(rebase.densityKgPerCubicMeter)
        || rebase.controlVolumes.empty()
        || rebase.ownedStorageBytes
            != storageBytesForControls(rebase.controlVolumes.size())
        || rebase.diagnostics.previousControlVolumeCount == 0
        || rebase.diagnostics.retainedControlVolumeCount
                + rebase.diagnostics.disappearedControlVolumeCount
            != rebase.diagnostics.previousControlVolumeCount
        || rebase.diagnostics.retainedControlVolumeCount
                + rebase.diagnostics.appearedControlVolumeCount
            != rebase.controlVolumes.size()
        || !rebase.diagnostics.finite
        || rebase.fingerprint != rebaseFingerprint(rebase)) {
        throw std::invalid_argument(
            "scene fluid region rebase integrity is invalid");
    }
}

void validateSceneFluidRegionRebase(
    const SceneFluidRegionRebase& rebase,
    const SceneFluidRegionTransport& transport,
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureTopologyTransition& topologyTransition) {
    validateSources(
        transport, previousPressureVolumes, currentPressureVolumes,
        topologyTransition);
    validateSceneFluidRegionRebaseIntegrity(rebase);
    const SceneFluidRegionRebaseLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildRebase(
        transport, previousPressureVolumes, currentPressureVolumes,
        topologyTransition, unlimited);
    if (rebase != expected) {
        throw std::invalid_argument(
            "scene fluid region rebase payload is invalid");
    }
}

std::vector<double> rebaseSceneFluidPressureWarmStart(
    const SceneFluidPressureControlVolumeSet& previousPressureVolumes,
    const SceneFluidPressureControlVolumeSet& currentPressureVolumes,
    const SceneFluidPressureTopologyTransition& topologyTransition,
    const std::span<const double> previousPressurePascals) {
    validateSceneFluidPressureControlVolumeIntegrity(previousPressureVolumes);
    validateSceneFluidPressureControlVolumeIntegrity(currentPressureVolumes);
    validateSceneFluidPressureTopologyTransitionIntegrity(
        topologyTransition);
    if (previousPressurePascals.size()
            != previousPressureVolumes.controlVolumes.size()
        || previousPressureVolumes.surfaceDefinitionFingerprint
            != currentPressureVolumes.surfaceDefinitionFingerprint
        || previousPressureVolumes.structureDefinitionFingerprint
            != currentPressureVolumes.structureDefinitionFingerprint
        || previousPressureVolumes.cellCounts
            != currentPressureVolumes.cellCounts
        || previousPressureVolumes.lowerMeters
            != currentPressureVolumes.lowerMeters
        || previousPressureVolumes.upperMeters
            != currentPressureVolumes.upperMeters
        || previousPressureVolumes.acceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || currentPressureVolumes.acceptedStepCount
            != previousPressureVolumes.acceptedStepCount + 1
        || !(currentPressureVolumes.simulationTimeSeconds
            > previousPressureVolumes.simulationTimeSeconds)
        || topologyTransition.previousPressureControlVolumeFingerprint
            != previousPressureVolumes.fingerprint
        || topologyTransition.currentPressureControlVolumeFingerprint
            != currentPressureVolumes.fingerprint
        || topologyTransition.previousAcceptedStepCount
            != previousPressureVolumes.acceptedStepCount
        || topologyTransition.currentAcceptedStepCount
            != currentPressureVolumes.acceptedStepCount
        || !std::ranges::all_of(
            previousPressurePascals,
            [](const double pressure) { return std::isfinite(pressure); })) {
        throw std::invalid_argument(
            "scene fluid pressure warm-start rebase identity is invalid");
    }

    std::vector<double> result(
        currentPressureVolumes.controlVolumes.size(), 0.0);
    std::vector<bool> initialized(result.size(), false);
    std::vector<double> appearanceDonorArea(result.size(), 0.0);
    for (const auto& retained : topologyTransition.retainedControls) {
        result[retained.currentControlVolumeIndex] =
            previousPressurePascals[
                retained.previousControlVolumeIndex];
        initialized[retained.currentControlVolumeIndex] = true;
    }
    for (const auto& donor : topologyTransition.appearanceDonors) {
        result[donor.appearedCurrentControlVolumeIndex] +=
            donor.linkAreaSquareMeters
            * previousPressurePascals[
                donor.retainedPreviousControlVolumeIndex];
        appearanceDonorArea[donor.appearedCurrentControlVolumeIndex] +=
            donor.linkAreaSquareMeters;
        initialized[donor.appearedCurrentControlVolumeIndex] = true;
    }
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (appearanceDonorArea[index] > 0.0) {
            result[index] /= appearanceDonorArea[index];
        }
    }
    if (!std::ranges::all_of(initialized, [](const bool value) {
            return value;
        })
        || !std::ranges::all_of(
            result,
            [](const double pressure) { return std::isfinite(pressure); })) {
        throw std::overflow_error(
            "scene fluid pressure warm-start rebase is non-finite");
    }
    return result;
}

} // namespace simwing::fsi
