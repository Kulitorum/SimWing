#pragma once

#include "fluid/planar_region_fragment_opening.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningFluxVersion = 1;

struct PlanarPressureRegionFragmentOpeningVelocitySample {
    std::uint64_t patchStableId = 0;
    double relativeNormalVelocityMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocitySample&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPatchFlux {
    std::size_t patchIndex = 0;
    std::uint64_t patchStableId = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::size_t minusFragmentIndex = 0;
    std::size_t plusFragmentIndex = 0;
    std::size_t minusBaseComponentIndex = 0;
    std::size_t plusBaseComponentIndex = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    double areaSquareMeters = 0.0;
    double relativeNormalVelocityMetersPerSecond = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPatchFlux&) const = default;
};

struct PlanarPressureRegionFragmentOpeningFluxSummary {
    std::size_t openingIndex = 0;
    std::uint64_t openingStableId = 0;
    std::uint64_t surfaceStableId = 0;
    std::uint64_t negativeSideRegionStableId = 0;
    std::uint64_t positiveSideRegionStableId = 0;
    std::size_t patchCount = 0;
    double areaSquareMeters = 0.0;
    double relativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double areaWeightedRelativeNormalVelocityMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningFluxSummary&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningFragmentFlux {
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t baseComponentIndex = 0;
    std::size_t connectedComponentIndex = 0;
    std::size_t incidentOpeningPatchCount = 0;
    double outwardRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningFragmentFlux&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningBaseComponentFlux {
    std::size_t baseComponentIndex = 0;
    std::uint64_t baseComponentStableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t connectedComponentIndex = 0;
    std::uint64_t connectedComponentStableId = 0;
    std::size_t incidentOpeningPatchCount = 0;
    double outwardRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningBaseComponentFlux&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningConnectedComponentFlux {
    std::size_t connectedComponentIndex = 0;
    std::uint64_t connectedComponentStableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t openingPatchCount = 0;
    double outwardRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningConnectedComponentFlux&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningFluxLimits {
    PlanarPressureRegionFragmentOpeningLimits openingLimits;
    std::size_t maximumPatches = 20'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumBaseComponents = 20'000'000;
    std::size_t maximumConnectedComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Immutable oriented relative flow through every patch of one exact opening
// overlay. Positive velocity and volume flow travel from the patch's authored
// negative-side fragment to its positive-side fragment. Fragment and sealed-
// base-component ledgers therefore store positive outward flow on the minus
// side and the exact opposite contribution on the plus side.
//
// This state accepts prescribed kinematics only. It owns no conductance,
// pressure-drop law, aperture momentum/inertia, pressure RHS, geometry dV/dt,
// fabric-load subtraction, or worker state.
struct PlanarPressureRegionFragmentOpeningFluxState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningFluxVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    std::vector<PlanarPressureRegionFragmentOpeningPatchFlux> patches;
    std::vector<PlanarPressureRegionFragmentOpeningFluxSummary> openings;
    std::vector<PlanarPressureRegionFragmentOpeningFragmentFlux> fragments;
    std::vector<PlanarPressureRegionFragmentOpeningBaseComponentFlux>
        baseComponents;
    std::vector<PlanarPressureRegionFragmentOpeningConnectedComponentFlux>
        connectedComponents;
    double maximumAbsoluteRelativeNormalVelocityMetersPerSecond = 0.0;
    double maximumAbsolutePatchVolumeFlowRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteBaseComponentOutwardFlowRateCubicMetersPerSecond =
        0.0;
    double totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond = 0.0;
    double globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond =
        0.0;
    double conservationToleranceCubicMetersPerSecond = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningFluxState&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningFluxState
buildPlanarPressureRegionFragmentOpeningFluxState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningVelocitySample> samples,
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningFluxStateIntegrity(
    const PlanarPressureRegionFragmentOpeningFluxState& state);

void validatePlanarPressureRegionFragmentOpeningFluxState(
    const PlanarPressureRegionFragmentOpeningFluxState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningVelocitySample> samples,
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits = {});

} // namespace simwing::fsi::fluid
