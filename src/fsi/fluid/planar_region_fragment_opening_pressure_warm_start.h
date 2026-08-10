#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningPressureWarmStartVersion = 1;

struct PlanarPressureRegionFragmentOpeningPressureWarmStartLimits {
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits
        previousAcceptedStateLimits;
    PlanarPressureRegionFragmentOpeningPressureOperatorLimits
        currentPressureOperatorLimits;
    PlanarPressureRegionFragmentVolumeRateLimits currentVolumeRateLimits;
    std::size_t maximumPressureCorrections = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 1024ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 1024ULL * 1024ULL * 1024ULL;
};

// Pressure-only history for one consecutive topology-stable opening epoch.
// The previous accepted correction is mapped to current fragments by exact
// stable identity and shifted to the current opening-connected volume gauge.
// It deliberately carries no Cartesian velocity, aperture velocity, opening
// flux, momentum, or accepted endpoint.
struct PlanarPressureRegionFragmentOpeningPressureWarmStart {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningPressureWarmStartVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t previousPressureOperatorFingerprint = 0;
    std::uint64_t previousBasePressureOperatorFingerprint = 0;
    std::uint64_t previousOpeningFingerprint = 0;
    std::uint64_t previousFragmentFingerprint = 0;
    std::uint64_t previousTopologyFingerprint = 0;
    std::uint64_t previousVolumeRateFingerprint = 0;
    std::uint64_t currentPressureOperatorFingerprint = 0;
    std::uint64_t currentBasePressureOperatorFingerprint = 0;
    std::uint64_t currentOpeningFingerprint = 0;
    std::uint64_t currentFragmentFingerprint = 0;
    std::uint64_t currentTopologyFingerprint = 0;
    std::uint64_t currentVolumeRateFingerprint = 0;
    std::size_t pressureCorrectionCount = 0;
    std::size_t connectedComponentCount = 0;
    double maximumAbsoluteGaugeShiftPascals = 0.0;
    double maximumAbsolutePressureCorrectionPascals = 0.0;
    std::vector<double> componentGaugeShiftsPascals;
    std::vector<double> pressureCorrectionPascals;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureWarmStart&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningPressureWarmStart
buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
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
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart);

void validatePlanarPressureRegionFragmentOpeningPressureWarmStart(
    const PlanarPressureRegionFragmentOpeningPressureWarmStart& warmStart,
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
    const PlanarPressureRegionFragmentOpeningPressureWarmStartLimits& limits =
        {});

} // namespace simwing::fsi::fluid
