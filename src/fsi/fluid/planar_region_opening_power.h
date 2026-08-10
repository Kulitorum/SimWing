#pragma once

#include "fluid/planar_region_opening_flow.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t planarPressureRegionOpeningPowerVersion = 1;

struct PlanarPressureRegionOpeningPowerSettings {
    double absolutePowerToleranceWatts = 1.0e-12;
    double relativePowerTolerance = 1.0e-12;

    bool operator==(
        const PlanarPressureRegionOpeningPowerSettings&) const = default;
};

struct PlanarPressureRegionOpeningPowerLimits {
    std::size_t maximumOpenings = 16'384;
    std::size_t maximumRegions = 4096;
    std::size_t maximumOwnedBytes = 256ULL * 1024ULL * 1024ULL;
    PlanarPressureRegionOpeningFlowLimits sourceLimits;
};

struct PlanarPressureRegionOpeningPower {
    std::uint64_t openingStableId = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    double midpointNegativePressurePascals = 0.0;
    double midpointPositivePressurePascals = 0.0;
    double negativeToPositivePressureDropPascals = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double pressurePowerWatts = 0.0;
    double minimumExternalPowerWatts = 0.0;
    double passiveToleranceWatts = 0.0;
    bool passiveWithinTolerance = false;

    bool operator==(
        const PlanarPressureRegionOpeningPower&) const = default;
};

struct PlanarPressureRegionPressureVolumePower {
    std::uint64_t regionStableId = 0;
    double previousPressurePascals = 0.0;
    double currentPressurePascals = 0.0;
    double midpointPressurePascals = 0.0;
    double geometryVolumeRateCubicMetersPerSecond = 0.0;
    double pressureVolumePowerWatts = 0.0;

    bool operator==(
        const PlanarPressureRegionPressureVolumePower&) const = default;
};

// Pressure power is positive when pressure drops along the oriented opening
// flow and can supply fluid energy. Negative power is uphill flow and exposes
// the minimum external power needed by that opening. For compatible regional
// continuity, midpoint pressure work closes as
//
//   sum(opening pressure power) + sum(region p * dV/dt) = 0.
//
// The audit requires an already feasible opening allocation. It diagnoses
// energetic direction only and does not prescribe an opening law, pressure,
// momentum, or external energy source.
struct PlanarPressureRegionOpeningPowerAudit {
    std::uint32_t version = planarPressureRegionOpeningPowerVersion;
    std::uint64_t fingerprint = 0;
    std::uint32_t sourceSweepVersion = 0;
    std::uint32_t sourceOpeningFlowVersion = 0;
    std::uint64_t sourceOpeningFlowFingerprint = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    double durationSeconds = 0.0;
    PlanarPressureRegionOpeningPowerSettings settings;
    std::vector<PlanarPressureRegionOpeningPower> openings;
    std::vector<PlanarPressureRegionPressureVolumePower> regions;
    std::size_t failedPassiveOpeningCount = 0;
    double maximumOpeningExternalPowerWatts = 0.0;
    double summedOpeningExternalPowerDeficitWatts = 0.0;
    double totalOpeningPressurePowerWatts = 0.0;
    double minimumNetExternalPowerWatts = 0.0;
    double totalRegionPressureVolumePowerWatts = 0.0;
    double pressurePowerClosureResidualWatts = 0.0;
    double pressurePowerClosureToleranceWatts = 0.0;
    bool allOpeningsPassiveWithinTolerance = false;
    bool pressurePowerClosesWithinTolerance = false;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionOpeningPowerAudit&) const = default;
};

[[nodiscard]] PlanarPressureRegionOpeningPowerAudit
auditPlanarPressureRegionOpeningPower(
    const PlanarPressureRegionSweepLedger& sweep,
    std::span<const PlanarPressureRegionOpeningDefinition> openingDefinitions,
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionOpeningPowerSettings& settings = {},
    const PlanarPressureRegionOpeningPowerLimits& limits = {});

void validatePlanarPressureRegionOpeningPower(
    const PlanarPressureRegionOpeningPowerAudit& audit,
    const PlanarPressureRegionSweepLedger& sweep,
    std::span<const PlanarPressureRegionOpeningDefinition> openingDefinitions,
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const PlanarPressureRegionOpeningPowerLimits& limits = {});

} // namespace simwing::fsi::fluid
