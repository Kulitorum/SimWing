#include "scene_fluid_regional_opening_load_epoch.h"

#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening load-epoch storage overflows");
    }
    return first + second;
}

void validateLimits(
    const SceneFluidRegionalOpeningLoadEpochLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening load-epoch limits are invalid");
    }
}

std::uint64_t settingsFingerprint(
    const SceneFluidRegionalOpeningLoadEpochSettings& settings) {
    Fingerprint fingerprint;
    for (const double value : {
             settings.pressureState
                 .absolutePressureResidualTolerancePascals,
             settings.pressureState.relativePressureResidualTolerance,
             settings.pressureState.absoluteForceResidualToleranceNewtons,
             settings.pressureState.relativeForceResidualTolerance,
             settings.pressureState.absoluteWorkResidualToleranceJoules,
             settings.pressureState.relativeWorkResidualTolerance,
             settings.transfer.momentReferenceMeters.x,
             settings.transfer.momentReferenceMeters.y,
             settings.transfer.momentReferenceMeters.z,
             settings.transfer.minimumTriangleAreaSquareMeters,
             settings.transfer.minimumQuadratureAreaSquareMeters,
             settings.transfer.barycentricTolerance}) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningLoadEpoch& epoch) {
    return checkedAdd(
        checkedAdd(epoch.loadState.ownedStorageBytes,
                   epoch.samples.ownedStorageBytes),
        epoch.application.ownedStorageBytes);
}

std::size_t workingStorageBytes(
    const SceneFluidRegionalOpeningLoadEpoch& epoch) {
    std::size_t result = epoch.loadState.pressure.workingStorageBytes;
    result = checkedAdd(
        result, epoch.loadState.surfaceLoads.workingStorageBytes);
    result = checkedAdd(
        result, epoch.loadState.openingSurfaceLoads.workingStorageBytes);
    result = checkedAdd(result, epoch.samples.workingStorageBytes);
    return checkedAdd(result, epoch.application.workingStorageBytes);
}

std::uint64_t epochFingerprint(
    const SceneFluidRegionalOpeningLoadEpoch& epoch) {
    Fingerprint fingerprint;
    fingerprint.integer(epoch.version);
    for (const std::uint64_t value : {
             epoch.sourceAcceptedStateFingerprint,
             epoch.sourcePressureOperatorFingerprint,
             epoch.sourceBasePressureOperatorFingerprint,
             epoch.sourceOpeningFingerprint,
             epoch.sourceFragmentFingerprint,
             epoch.sourceTopologyFingerprint,
             epoch.sourceVolumeRateFingerprint,
             epoch.surfaceDefinitionFingerprint,
             epoch.surfaceStateFingerprint,
             epoch.quadratureFingerprint,
             epoch.couplingSurfaceFingerprint,
             epoch.targetDefinitionFingerprint,
             epoch.sourceSettingsFingerprint,
             epoch.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(epoch.simulationTimeSeconds);
    fingerprint.integer(epoch.loadState.fingerprint);
    fingerprint.integer(epoch.samples.fingerprint);
    fingerprint.integer(epoch.application.fingerprint);
    fingerprint.boolean(epoch.applied);
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(epoch.workingStorageBytes));
    return fingerprint.value();
}

fluid::PlanarPressureRegionFragmentOpeningLoadState buildLoadState(
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        acceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningLoadEpochLimits& limits) {
    const auto pressure =
        fluid::composePlanarPressureRegionFragmentOpeningPressureState(
            acceptedState, pressureOperator, basePressureOperator, grid,
            sweep, fragments, topology, volumeRates, openingDefinitions,
            openings, resistanceDefinitions, settings.pressureState,
            limits.loadState.pressureStateLimits);
    const auto surfaceLoads =
        fluid::capturePlanarPressureRegionFragmentSurfaceLoads(
            pressure,
            limits.loadState.surfaceLoadLimits.surfaceLoadLimits);
    const auto openingSurfaceLoads =
        fluid::capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressure, grid, sweep, fragments, topology,
            openingDefinitions, openings,
            limits.loadState.surfaceLoadLimits);
    return fluid::capturePlanarPressureRegionFragmentOpeningLoadState(
        acceptedState, pressure, surfaceLoads, openingSurfaceLoads,
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits.loadState);
}

void validateAggregateLimits(
    const SceneFluidRegionalOpeningLoadEpoch& epoch,
    const SceneFluidRegionalOpeningLoadEpochLimits& limits) {
    if (epoch.ownedStorageBytes > limits.maximumOwnedBytes
        || epoch.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening load epoch exceeds aggregate limits");
    }
}

} // namespace

SceneFluidRegionalOpeningLoadEpoch
applySceneFluidRegionalOpeningLoadEpoch(
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        acceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningLoadEpochLimits& limits) {
    validateLimits(limits);
    auto loadState = buildLoadState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, settings, limits);
    auto samples = sampleSceneFluidRegionalOpeningPressure(
        loadState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, surface, surfaceState, quadrature,
        limits.loadState, limits.sampling);

    const StructureCheckpoint before = target.checkpoint();
    try {
        auto application =
            applySceneFluidRegionalAcceptedPressureLoads(
                surface, surfaceState, transfer, quadrature, samples,
                target, settings.transfer, limits.application);

        SceneFluidRegionalOpeningLoadEpoch result;
        result.sourceAcceptedStateFingerprint = acceptedState.fingerprint;
        result.sourcePressureOperatorFingerprint =
            pressureOperator.fingerprint;
        result.sourceBasePressureOperatorFingerprint =
            basePressureOperator.fingerprint;
        result.sourceOpeningFingerprint = openings.fingerprint;
        result.sourceFragmentFingerprint = fragments.fingerprint;
        result.sourceTopologyFingerprint = topology.fingerprint;
        result.sourceVolumeRateFingerprint = volumeRates.fingerprint;
        result.surfaceDefinitionFingerprint = surface.fingerprint;
        result.surfaceStateFingerprint = surfaceState.fingerprint;
        result.quadratureFingerprint = quadrature.fingerprint;
        result.couplingSurfaceFingerprint =
            transfer.couplingSurfaceFingerprint();
        result.targetDefinitionFingerprint =
            transfer.targetDefinitionFingerprint();
        result.sourceSettingsFingerprint =
            settingsFingerprint(settings);
        result.acceptedStepCount = surfaceState.acceptedStepCount;
        result.simulationTimeSeconds =
            surfaceState.simulationTimeSeconds;
        result.loadState = std::move(loadState);
        result.samples = std::move(samples);
        result.application = std::move(application);
        result.applied = true;
        result.ownedStorageBytes = ownedStorageBytes(result);
        result.workingStorageBytes = workingStorageBytes(result);
        result.fingerprint = epochFingerprint(result);
        validateSceneFluidRegionalOpeningLoadEpoch(
            result, acceptedState, pressureOperator, basePressureOperator,
            grid, sweep, fragments, topology, volumeRates,
            openingDefinitions, openings, resistanceDefinitions, surface,
            surfaceState, transfer, quadrature, settings, limits);
        return result;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalOpeningLoadEpochIntegrity(
    const SceneFluidRegionalOpeningLoadEpoch& epoch) {
    fluid::validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
        epoch.loadState);
    validateSceneFluidRegionalPressureSampleIntegrity(epoch.samples);
    validateSceneFluidRegionalPressureLoadApplicationIntegrity(
        epoch.application);
    if (epoch.version != sceneFluidRegionalOpeningLoadEpochVersion
        || epoch.fingerprint == 0
        || epoch.sourceAcceptedStateFingerprint == 0
        || epoch.sourcePressureOperatorFingerprint == 0
        || epoch.sourceBasePressureOperatorFingerprint == 0
        || epoch.sourceOpeningFingerprint == 0
        || epoch.sourceFragmentFingerprint == 0
        || epoch.sourceTopologyFingerprint == 0
        || epoch.sourceVolumeRateFingerprint == 0
        || epoch.surfaceDefinitionFingerprint == 0
        || epoch.surfaceStateFingerprint == 0
        || epoch.quadratureFingerprint == 0
        || epoch.couplingSurfaceFingerprint == 0
        || epoch.targetDefinitionFingerprint == 0
        || epoch.sourceSettingsFingerprint == 0
        || !std::isfinite(epoch.simulationTimeSeconds)
        || epoch.sourceAcceptedStateFingerprint
            != epoch.loadState.sourceAcceptedStateFingerprint
        || epoch.sourcePressureOperatorFingerprint
            != epoch.loadState.pressure
                   .sourcePressureOperatorFingerprint
        || epoch.sourceBasePressureOperatorFingerprint
            != epoch.loadState.pressure
                   .sourceBasePressureOperatorFingerprint
        || epoch.sourceOpeningFingerprint
            != epoch.loadState.sourceOpeningFingerprint
        || epoch.sourceFragmentFingerprint
            != epoch.loadState.sourceFragmentFingerprint
        || epoch.sourceTopologyFingerprint
            != epoch.loadState.sourceTopologyFingerprint
        || epoch.sourceVolumeRateFingerprint
            != epoch.loadState.sourceVolumeRateFingerprint
        || epoch.surfaceDefinitionFingerprint
            != epoch.samples.surfaceDefinitionFingerprint
        || epoch.surfaceStateFingerprint
            != epoch.samples.surfaceStateFingerprint
        || epoch.quadratureFingerprint
            != epoch.samples.quadratureFingerprint
        || epoch.couplingSurfaceFingerprint
            != epoch.application.couplingSurfaceFingerprint
        || epoch.targetDefinitionFingerprint
            != epoch.application.targetDefinitionFingerprint
        || epoch.acceptedStepCount != epoch.samples.acceptedStepCount
        || epoch.acceptedStepCount
            != epoch.application.acceptedStepCount
        || epoch.simulationTimeSeconds
            != epoch.samples.simulationTimeSeconds
        || epoch.simulationTimeSeconds
            != epoch.application.simulationTimeSeconds
        || epoch.samples.regionalOpeningLoadStateFingerprint
            != epoch.loadState.fingerprint
        || epoch.application.sourceSamplingFingerprint
            != epoch.samples.fingerprint
        || !epoch.applied || !epoch.application.applied
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.workingStorageBytes != workingStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "regional opening load-epoch integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningLoadEpoch(
    const SceneFluidRegionalOpeningLoadEpoch& epoch,
    const fluid::PlanarPressureRegionFragmentOpeningAcceptedState&
        acceptedState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningLoadEpochSettings& settings,
    const SceneFluidRegionalOpeningLoadEpochLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningLoadEpochIntegrity(epoch);
    if (epoch.sourceAcceptedStateFingerprint != acceptedState.fingerprint
        || epoch.sourcePressureOperatorFingerprint
            != pressureOperator.fingerprint
        || epoch.sourceBasePressureOperatorFingerprint
            != basePressureOperator.fingerprint
        || epoch.sourceOpeningFingerprint != openings.fingerprint
        || epoch.sourceFragmentFingerprint != fragments.fingerprint
        || epoch.sourceTopologyFingerprint != topology.fingerprint
        || epoch.sourceVolumeRateFingerprint != volumeRates.fingerprint
        || epoch.surfaceDefinitionFingerprint != surface.fingerprint
        || epoch.surfaceStateFingerprint != surfaceState.fingerprint
        || epoch.quadratureFingerprint != quadrature.fingerprint
        || epoch.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || epoch.targetDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || epoch.sourceSettingsFingerprint
            != settingsFingerprint(settings)) {
        throw std::invalid_argument(
            "regional opening load epoch is foreign to its source");
    }
    validateAggregateLimits(epoch, limits);

    const auto expectedLoadState = buildLoadState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, settings, limits);
    if (epoch.loadState != expectedLoadState) {
        throw std::invalid_argument(
            "regional opening load epoch changed its load state");
    }
    const auto expectedSamples = sampleSceneFluidRegionalOpeningPressure(
        expectedLoadState, pressureOperator, basePressureOperator, grid,
        sweep, fragments, topology, volumeRates, openingDefinitions,
        openings, resistanceDefinitions, surface, surfaceState, quadrature,
        limits.loadState, limits.sampling);
    if (epoch.samples != expectedSamples) {
        throw std::invalid_argument(
            "regional opening load epoch changed its scene samples");
    }
    validateSceneFluidRegionalPressureLoadApplication(
        epoch.application, surface, surfaceState, transfer, quadrature,
        epoch.samples, settings.transfer, limits.application);
}

} // namespace simwing::fsi
