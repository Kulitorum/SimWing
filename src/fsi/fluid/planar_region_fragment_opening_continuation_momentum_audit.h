#pragma once

#include "fluid/planar_region_fragment_opening_continuation.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningContinuationMomentumAuditVersion = 1;

enum class PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind
    : std::uint8_t {
    SameRegionGrid = 1,
    OpeningPatch = 2,
};

struct PlanarPressureRegionFragmentOpeningContinuationMomentumSample {
    std::size_t sampleIndex = 0;
    PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind kind =
        PlanarPressureRegionFragmentOpeningContinuationMomentumDofKind::
            SameRegionGrid;
    std::uint64_t stableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t previousSourceIndex = 0;
    std::size_t currentSourceIndex = 0;
    double previousDualVolumeCubicMeters = 0.0;
    double currentDualVolumeCubicMeters = 0.0;
    double carriedVelocityMetersPerSecond = 0.0;
    double previousMassKilograms = 0.0;
    double currentMassKilograms = 0.0;
    double previousMomentumKilogramMetersPerSecond = 0.0;
    double carriedMomentumKilogramMetersPerSecond = 0.0;
    double momentumChangeKilogramMetersPerSecond = 0.0;
    double previousKineticEnergyJoules = 0.0;
    double carriedKineticEnergyJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    bool metricChanged = false;
    bool momentumChanged = false;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningContinuationMomentumSample&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits {
    PlanarPressureRegionFragmentOpeningContinuationLimits continuation;
    std::size_t maximumSamples = 160'000'000;
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Read-only conservation audit for the current stable-ID velocity carry.
// Every same-region Cartesian degree uses rho*area*centre-distance mass; every
// aperture degree uses the identical plug mass owned by resistance/projection.
// The audit evaluates the previous and current diagonal momentum/energy at the
// unchanged carried velocity and publishes every metric-induced change.
//
// This deliberately does not repair the change. Local mass rescaling would
// preserve diagonal momentum but can violate ALE free-stream/GCL behavior;
// only a swept-volume transport operator can decide where exchanged momentum
// belongs. The artifact therefore supplies blocker evidence, not advection,
// a replacement continuation, pressure acceptance, or worker selection.
struct PlanarPressureRegionFragmentOpeningContinuationMomentumAudit {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningContinuationMomentumAuditVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceContinuationFingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t previousTopologyFingerprint = 0;
    std::uint64_t previousOpeningFingerprint = 0;
    std::uint64_t currentTopologyFingerprint = 0;
    std::uint64_t currentOpeningFingerprint = 0;
    double densityKgPerCubicMeter = 0.0;
    std::vector<
        PlanarPressureRegionFragmentOpeningContinuationMomentumSample>
        samples;
    std::size_t sameRegionGridSampleCount = 0;
    std::size_t openingPatchSampleCount = 0;
    std::size_t metricChangedSampleCount = 0;
    std::size_t momentumChangedSampleCount = 0;
    Vector3 previousMassByAxisKilograms;
    Vector3 currentMassByAxisKilograms;
    Vector3 previousMomentumKilogramMetersPerSecond;
    Vector3 carriedMomentumKilogramMetersPerSecond;
    Vector3 momentumChangeKilogramMetersPerSecond;
    double previousKineticEnergyJoules = 0.0;
    double carriedKineticEnergyJoules = 0.0;
    double kineticEnergyChangeJoules = 0.0;
    double maximumAbsoluteDualVolumeChangeCubicMeters = 0.0;
    double maximumRelativeDualVolumeChange = 0.0;
    double maximumAbsoluteMomentumChangeKilogramMetersPerSecond = 0.0;
    double maximumAbsoluteKineticEnergyChangeJoules = 0.0;
    bool metricChanged = false;
    bool warmCarryChangesMomentum = false;
    bool audited = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit&)
        const = default;
};

[[nodiscard]]
PlanarPressureRegionFragmentOpeningContinuationMomentumAudit
auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
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
    std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits = {});

void validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit);

void validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAudit& audit,
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
    std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        previousOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& previousOpenings,
    std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        previousResistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureOperator&
        currentPressureOperator,
    const PlanarPressureRegionFragmentPressureOperator&
        currentBasePressureOperator,
    const PlanarPressureRegionSweepLedger& currentSweep,
    const PlanarPressureRegionFragmentSet& currentFragments,
    const PlanarPressureRegionFragmentTopology& currentTopology,
    const PlanarPressureRegionFragmentVolumeRateSet& currentVolumeRates,
    std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        currentOpeningDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& currentOpenings,
    const PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits&
        limits = {});

} // namespace simwing::fsi::fluid
