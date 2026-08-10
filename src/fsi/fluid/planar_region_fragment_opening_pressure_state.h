#pragma once

#include "fluid/planar_region_fragment_opening_accepted_state.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningPressureStateVersion = 2;

struct PlanarPressureRegionFragmentOpeningPressureStateSettings {
    double absolutePressureResidualTolerancePascals = 1.0e-12;
    double relativePressureResidualTolerance = 1.0e-10;
    double absoluteForceResidualToleranceNewtons = 1.0e-12;
    double relativeForceResidualTolerance = 1.0e-10;
    double absoluteWorkResidualToleranceJoules = 1.0e-12;
    double relativeWorkResidualTolerance = 1.0e-10;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStateSettings&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureStateControl {
    std::size_t fragmentIndex = 0;
    std::uint64_t fragmentStableId = 0;
    std::size_t baseComponentIndex = 0;
    std::size_t connectedComponentIndex = 0;
    std::uint64_t connectedComponentStableId = 0;
    std::uint64_t regionStableId = 0;
    double volumeCubicMeters = 0.0;
    double authoredPressurePascals = 0.0;
    double correctionPressurePascals = 0.0;
    double totalPressurePascals = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStateControl&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureStateWall {
    std::size_t wallIndex = 0;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::uint64_t surfaceStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::size_t minusFragmentIndex = 0;
    std::size_t plusFragmentIndex = 0;
    std::size_t minusBaseComponentIndex = 0;
    std::size_t plusBaseComponentIndex = 0;
    std::size_t minusConnectedComponentIndex = 0;
    std::size_t plusConnectedComponentIndex = 0;
    std::uint64_t minusRegionStableId = 0;
    std::uint64_t plusRegionStableId = 0;
    double areaSquareMeters = 0.0;
    Vector3 wrappedCentroidMeters;
    Vector3 unitNormalMinusToPlus;
    double minusTotalPressurePascals = 0.0;
    double plusTotalPressurePascals = 0.0;
    double authoredPressureJumpPascals = 0.0;
    double correctionPressureJumpPascals = 0.0;
    double totalPressureJumpPascals = 0.0;
    double pressureSplitResidualPascals = 0.0;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 pressureForceSplitResidualNewtons;
    Vector3 totalPressureImpulseOnSheetNewtonSeconds;
    double materialWallVelocityMetersPerSecond = 0.0;
    double authoredPressureWorkToFluidJoules = 0.0;
    double correctionPressureWorkToFluidJoules = 0.0;
    double totalPressureWorkToFluidJoules = 0.0;
    double pressureWorkSplitResidualJoules = 0.0;
    double totalPressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStateWall&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureStateComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t fragmentCount = 0;
    double volumeCubicMeters = 0.0;
    double authoredVolumeMeanPressurePascals = 0.0;
    double correctionVolumeMeanPressurePascals = 0.0;
    double totalVolumeMeanPressurePascals = 0.0;
    double pressureMeanSplitResidualPascals = 0.0;
    double authoredGeometryPressureWorkJoules = 0.0;
    double correctionGeometryPressureWorkJoules = 0.0;
    double totalGeometryPressureWorkJoules = 0.0;
    double totalWallPressureWorkToFluidJoules = 0.0;
    double wallGeometryWorkResidualJoules = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureStateComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningPressureStateLimits {
    PlanarPressureRegionFragmentOpeningAcceptedStateLimits acceptedStateLimits;
    std::size_t maximumControls = 20'000'000;
    std::size_t maximumWalls = 60'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 2048ULL * 1024ULL * 1024ULL;
};

// Immutable authored-plus-correction pressure composition over one accepted
// active-opening endpoint. Correction gauges follow the opening-connected
// pressure operator rather than the sealed base components. Every full source
// wall publishes pressure jump, sheet force/impulse, and moving material work;
// connected-component and global wall work close independently to
// -dt*sum(p*dV/dt).
//
// The full-wall loads are diagnostic source data. This state neither removes
// aperture area from fabric traction nor applies a load to Structure. It owns
// no pressure relaxation, opening resistance selection, topology rebase, or
// production worker state.
struct PlanarPressureRegionFragmentOpeningPressureState {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningPressureStateVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceAcceptedStateFingerprint = 0;
    std::uint64_t sourcePressureOperatorFingerprint = 0;
    std::uint64_t sourceBasePressureOperatorFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    std::uint64_t sourceVolumeRateFingerprint = 0;
    PlanarPressureRegionFragmentOpeningPressureStateSettings settings;
    bool staticGeometry = false;
    bool usesMovingVolumeRates = false;
    double timeStepSeconds = 0.0;
    std::vector<PlanarPressureRegionFragmentOpeningPressureStateControl>
        controls;
    std::vector<PlanarPressureRegionFragmentOpeningPressureStateWall> walls;
    std::vector<PlanarPressureRegionFragmentOpeningPressureStateComponent>
        components;
    double maximumAbsoluteCorrectionGaugePascals = 0.0;
    double maximumAbsolutePressureSplitResidualPascals = 0.0;
    double maximumAbsoluteForceSplitResidualNewtons = 0.0;
    double maximumAbsoluteWorkResidualJoules = 0.0;
    Vector3 authoredPressureForceOnSheetNewtons;
    Vector3 correctionPressureForceOnSheetNewtons;
    Vector3 totalPressureForceOnSheetNewtons;
    Vector3 pressureForceSplitResidualNewtons;
    Vector3 totalPressureImpulseOnSheetNewtonSeconds;
    double authoredPressureWorkToFluidJoules = 0.0;
    double correctionPressureWorkToFluidJoules = 0.0;
    double totalPressureWorkToFluidJoules = 0.0;
    double pressureWorkSplitResidualJoules = 0.0;
    double totalPressureWorkToSheetJoules = 0.0;
    double actionReactionWorkResidualJoules = 0.0;
    double authoredGeometryPressureWorkJoules = 0.0;
    double correctionGeometryPressureWorkJoules = 0.0;
    double totalGeometryPressureWorkJoules = 0.0;
    double wallGeometryWorkResidualJoules = 0.0;
    bool accepted = false;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningPressureState&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningPressureState
composePlanarPressureRegionFragmentOpeningPressureState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStateSettings& settings = {},
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits = {});

void validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(
    const PlanarPressureRegionFragmentOpeningPressureState& state);

void validatePlanarPressureRegionFragmentOpeningPressureState(
    const PlanarPressureRegionFragmentOpeningPressureState& state,
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    std::span<const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStateLimits& limits = {});

} // namespace simwing::fsi::fluid
