#include "scene_fluid_regional_opening_momentum_wall_coupled_state.h"

#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening wall coupled-state storage overflows");
    }
    return first + second;
}

std::uint64_t fingerprint(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state) {
    std::uint64_t result = fnvOffsetBasis;
    const auto append = [&](std::uint64_t value) {
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            result ^= static_cast<std::uint8_t>(value & 0xffU);
            result *= fnvPrime;
            value >>= 8U;
        }
    };
    append(state.version);
    append(state.cycleState.fingerprint);
    append(state.structureStep.fingerprint);
    append(static_cast<std::uint64_t>(state.ownedStorageBytes));
    return result == 0 ? 1 : result;
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state) {
    return checkedAdd(
        state.cycleState.ownedStorageBytes,
        state.structureStep.ownedStorageBytes);
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    if (limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "regional opening wall coupled-state limits are invalid");
    }
}

void validateOwnedLimit(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    if (state.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "regional opening wall coupled-state exceeds its owned-byte limit");
    }
}

StructureCheckpoint decodePostStepCheckpoint(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    const Structure& owner,
    const StructureCheckpointPersistenceLimits& limits) {
    StructureCheckpoint result;
    StructureCheckpointPersistenceError error;
    if (!deserializeStructureCheckpoint(
            state.structureStep.afterStructureCheckpoint, owner, result,
            &error, limits)) {
        const std::string message =
            "cannot restore regional opening wall coupled Structure: "
            + error.message;
        if (error.code
            == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
            throw std::length_error(message);
        }
        throw std::invalid_argument(message);
    }
    return result;
}

std::optional<SceneFluidRegionalOpeningMomentumWallCoupledState>
advanceCandidate(
    Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    validateLimits(limits);
    std::optional<SceneFluidRegionalOpeningMomentumWallCoupledState> candidate;
    candidate.emplace();
    candidate->cycleState = cycleState;

    const StructureCheckpoint before = structure.checkpoint();
    try {
        candidate->structureStep =
            advanceSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
                candidate->cycleState, transportMetric,
                acceptedPressureOperator, acceptedBasePressureOperator,
                grid, acceptedSweep, acceptedFragments, acceptedTopology,
                acceptedVolumeRates, acceptedOpeningDefinitions,
                acceptedOpenings, acceptedResistanceDefinitions,
                acceptedBaseMetric, acceptedMetric, surface, surfaceState,
                transfer, quadrature, structure, settings,
                limits.structureStep);
        candidate->ownedStorageBytes = ownedStorageBytes(*candidate);
        validateOwnedLimit(*candidate, limits);
        candidate->fingerprint = fingerprint(*candidate);
        validateSceneFluidRegionalOpeningMomentumWallCoupledState(
            *candidate, transportMetric, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, surface, surfaceState, transfer, quadrature,
            structure, settings, limits);
        return candidate;
    } catch (...) {
        structure.restore(before);
        throw;
    }
}

} // namespace

void validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state) {
    validateSceneFluidRegionalOpeningMomentumWallCycleStateIntegrity(
        state.cycleState);
    validateSceneFluidRegionalOpeningMomentumWallStructureStepEpochIntegrity(
        state.structureStep);
    if (state.version
            != sceneFluidRegionalOpeningMomentumWallCoupledStateVersion
        || state.fingerprint == 0
        || state.structureStep.sourceCycleStateFingerprint
            != state.cycleState.fingerprint
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != fingerprint(state)) {
        throw std::invalid_argument(
            "regional opening wall coupled-state integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallCoupledState(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const Structure& structureDefinitionOwner,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(state);
    validateOwnedLimit(state, limits);
    validateSceneFluidRegionalOpeningMomentumWallStructureStepEpoch(
        state.structureStep, state.cycleState, transportMetric,
        acceptedPressureOperator, acceptedBasePressureOperator, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, acceptedOpeningDefinitions,
        acceptedOpenings, acceptedResistanceDefinitions,
        acceptedBaseMetric, acceptedMetric, surface, surfaceState,
        transfer, quadrature, structureDefinitionOwner, settings,
        limits.structureStep);
}

SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::
SceneFluidRegionalOpeningMomentumWallCoupledStateOwner(Structure structure)
    : structure_(std::move(structure)) {}

const Structure& SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::
structure() const noexcept {
    return structure_;
}

bool SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::hasState()
    const noexcept {
    return state_.has_value();
}

const SceneFluidRegionalOpeningMomentumWallCoupledState&
SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::state() const {
    if (!state_) {
        throw std::logic_error(
            "regional opening wall coupled-state owner has no accepted state");
    }
    return *state_;
}

SceneFluidRegionalOpeningMomentumWallCoupledState
SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::checkpoint() const {
    const auto& current = state();
    validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
        current);
    return current;
}

void SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::advance(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    if (state_) {
        throw std::logic_error(
            "regional opening wall coupled-state bootstrap requires an empty owner");
    }
    auto candidate = advanceCandidate(
        structure_, cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric, acceptedMetric,
        surface, surfaceState, transfer, quadrature, settings, limits);
    static_assert(std::is_nothrow_swappable_v<
        SceneFluidRegionalOpeningMomentumWallCoupledState>);
    state_.swap(candidate);
}

void SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::
advanceFixedMetricConsecutive(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        consecutiveTransport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    validateLimits(limits);
    const auto& prior = state();
    validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
        prior);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        consecutiveTransport);
    if (prior.cycleState.transportMetricFingerprint
            != transportMetric.fingerprint
        || prior.cycleState.acceptedMetricFingerprint
            != transportMetric.fingerprint
        || acceptedMetric.fingerprint != transportMetric.fingerprint) {
        throw std::invalid_argument(
            "consecutive regional opening wall advance requires one fixed fluid metric");
    }

    const auto priorAcceptedFlow =
        fluid::capturePlanarPressureRegionFragmentOpeningVelocityState(
            prior.cycleState.acceptedPressure, acceptedPressureOperator,
            acceptedBasePressureOperator, grid, acceptedSweep,
            acceptedFragments, acceptedTopology, acceptedVolumeRates,
            acceptedOpeningDefinitions, acceptedOpenings,
            acceptedResistanceDefinitions, acceptedBaseMetric,
            acceptedMetric, limits.consecutiveFlow);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumTransport(
        consecutiveTransport, prior.cycleState.adjustedMomentum,
        transportMetric, priorAcceptedFlow, transportMetric, grid,
        acceptedSweep, acceptedFragments, acceptedTopology,
        acceptedVolumeRates, limits.consecutiveTransport);
    if (!consecutiveTransport.diagnostics.accepted
        || consecutiveTransport.sourceAdjustmentStateFingerprint
            != prior.cycleState.adjustedMomentum.fingerprint
        || consecutiveTransport.targetFlowStateFingerprint
            != priorAcceptedFlow.fingerprint
        || cycleState.adjustedMomentum.sourceTransportFingerprint
            != consecutiveTransport.fingerprint
        || cycleState.transportMetricFingerprint
            != consecutiveTransport.targetMetricFingerprint) {
        throw std::invalid_argument(
            "consecutive regional opening wall cycle lineage is invalid");
    }

    auto candidate = advanceCandidate(
        structure_, cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric, acceptedMetric,
        surface, surfaceState, transfer, quadrature, settings, limits);
    static_assert(std::is_nothrow_swappable_v<
        SceneFluidRegionalOpeningMomentumWallCoupledState>);
    state_.swap(candidate);
}

void SceneFluidRegionalOpeningMomentumWallCoupledStateOwner::restore(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& state,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        transportMetric,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        acceptedPressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        acceptedBasePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& acceptedSweep,
    const fluid::PlanarPressureRegionFragmentSet& acceptedFragments,
    const fluid::PlanarPressureRegionFragmentTopology& acceptedTopology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet&
        acceptedVolumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        acceptedOpeningDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& acceptedOpenings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        acceptedResistanceDefinitions,
    const fluid::PlanarPressureRegionFragmentVelocityMetric&
        acceptedBaseMetric,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric&
        acceptedMetric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallStructureStepEpochSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallCoupledStateLimits& limits) {
    std::optional<SceneFluidRegionalOpeningMomentumWallCoupledState> candidate{
        state};
    validateSceneFluidRegionalOpeningMomentumWallCoupledState(
        *candidate, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric, acceptedMetric,
        surface, surfaceState, transfer, quadrature, structure_, settings,
        limits);
    const StructureCheckpoint restoredCheckpoint = decodePostStepCheckpoint(
        *candidate, structure_,
        limits.structureStep.checkpointPersistence);
    Structure restoredStructure(structure_.definition());
    restoredStructure.restore(restoredCheckpoint);

    static_assert(std::is_nothrow_move_assignable_v<Structure>);
    static_assert(std::is_nothrow_swappable_v<
        SceneFluidRegionalOpeningMomentumWallCoupledState>);
    structure_ = std::move(restoredStructure);
    state_.swap(candidate);
}

} // namespace simwing::fsi
