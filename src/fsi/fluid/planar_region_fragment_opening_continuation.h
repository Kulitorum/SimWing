#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningContinuationVersion = 1;

struct PlanarPressureRegionFragmentOpeningContinuationLimits {
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits
        previousAcceptedStateLimits;
    PlanarPressureRegionFragmentOpeningPressureOperatorLimits
        currentPressureOperatorLimits;
    PlanarPressureRegionFragmentVolumeRateLimits currentVolumeRateLimits;
    PlanarPressureRegionFragmentOpeningFluxLimits currentOpeningFluxLimits;
    std::size_t maximumTopologyLinkVelocities = 140'000'000;
    std::size_t maximumOpeningSamples = 20'000'000;
    std::size_t maximumPressureCorrections = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Immutable initialization for the next topology-stable moving aperture
// epoch. Every Cartesian-link velocity, material-relative aperture velocity,
// and correction-pressure sample is matched by stable identity. The pressure
// field is then shifted to the current opening-connected volume gauge and the
// current opening-flux ledger is rebuilt from the mapped samples.
//
// The endpoint profiles must be exactly consecutive and all identities must
// have one-to-one coverage. Appearance, retirement, changed aperture
// semantics, and topology rebases reject; there is deliberately no donor or
// repair rule. This is a bounded warm-state transfer, not momentum advection,
// an accepted next pressure step, persistence, or production selection.
struct PlanarPressureRegionFragmentOpeningContinuation {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningContinuationVersion;
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
    std::uint64_t currentOpeningFluxFingerprint = 0;
    std::size_t topologyLinkCount = 0;
    std::size_t openingPatchCount = 0;
    std::size_t pressureCorrectionCount = 0;
    std::size_t connectedComponentCount = 0;
    double maximumAbsoluteGaugeShiftPascals = 0.0;
    double maximumAbsolutePressureCorrectionPascals = 0.0;
    std::vector<double> componentGaugeShiftsPascals;
    std::vector<double>
        orientedTopologyLinkVelocityMetersPerSecond;
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        openingVelocitySamples;
    PlanarPressureRegionFragmentOpeningFluxState openingFlux;
    std::vector<double> pressureCorrectionPascals;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningContinuation&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningContinuation
buildPlanarPressureRegionFragmentOpeningContinuation(
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
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation);

void validatePlanarPressureRegionFragmentOpeningContinuation(
    const PlanarPressureRegionFragmentOpeningContinuation& continuation,
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
    const PlanarPressureRegionFragmentOpeningContinuationLimits& limits = {});

} // namespace simwing::fsi::fluid
