#pragma once

#include "fluid/planar_region_fragment_opening_continuation.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningPressureEpochVersion = 1;

enum class PlanarPressureRegionFragmentOpeningPressureEpochFailureStage
    : std::uint8_t {
    None = 0,
    Resistance = 1,
    PressureProjection = 2,
    AggregateEnergy = 3,
};

struct PlanarPressureRegionFragmentOpeningPressureEpochLimits {
    PlanarPressureRegionFragmentOpeningContinuationLimits continuation;
    PlanarPressureRegionFragmentOpeningPressureStepLimits pressureStep;
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits acceptedState;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
};

struct PlanarPressureRegionFragmentOpeningPressureEpochDiagnostics {
    bool accepted = false;
    bool usedConsecutiveContinuation = false;
    PlanarPressureRegionFragmentOpeningPressureEpochFailureStage
        failureStage =
            PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::None;
    PlanarPressureRegionFragmentOpeningPressureStepDiagnostics pressureStep;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureEpochDiagnostics&)
        const = default;
};

// Atomic continuation-to-acceptance transaction for one topology-stable
// active-aperture epoch. The stable-ID continuation is built privately, its
// four mutable fields are advanced through resistance and augmented pressure
// projection, and an immutable accepted endpoint is captured only after the
// complete step accepts. A rejected numerical step returns diagnostics and no
// endpoint; neither source epoch can be altered.
//
// Input/source errors still throw before numerical publication. This is an
// opt-in regional harness: direct velocity continuation is not conservative
// advection, and this result owns no restart codec, topology rebase, Structure
// load application, or production worker selection.
struct PlanarPressureRegionFragmentOpeningPressureEpochResult {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningPressureEpochVersion;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t continuationFingerprint = 0;
    std::uint64_t continuationOpeningFluxFingerprint = 0;
    std::uint64_t currentPressureOperatorFingerprint = 0;
    std::uint64_t currentBasePressureOperatorFingerprint = 0;
    std::uint64_t currentOpeningFingerprint = 0;
    std::uint64_t currentFragmentFingerprint = 0;
    std::uint64_t currentTopologyFingerprint = 0;
    std::uint64_t currentVolumeRateFingerprint = 0;
    std::uint64_t currentResistanceDefinitionFingerprint = 0;
    PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
    PlanarPressureRegionFragmentOpeningPressureEpochDiagnostics diagnostics;
    PlanarPressureRegionFragmentOpeningAcceptedState acceptedState;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureEpochResult&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningPressureEpochResult
acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningPressureEpochLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureEpochResult& result);

void validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
    const PlanarPressureRegionFragmentOpeningPressureEpochResult& result,
    const PlanarPressureRegionFragmentOpeningAcceptedState& previousState,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        previousPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        previousBasePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& previousSweep,
    const PlanarPressureRegionFragmentSet& previousFragments,
    const PlanarPressureRegionFragmentTopology& previousTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& previousVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        currentResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningPressureEpochLimits& limits = {});

} // namespace simwing::fsi::fluid
