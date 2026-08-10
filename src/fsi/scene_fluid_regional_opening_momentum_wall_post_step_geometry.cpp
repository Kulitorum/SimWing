#include "scene_fluid_regional_opening_momentum_wall_post_step_geometry.h"

#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
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
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening wall post-step geometry storage overflows");
    }
    return first + second;
}

std::size_t surfaceStateStorageBytes(
    const SceneFluidSurfaceState& state) {
    if (state.vertices.size() > std::numeric_limits<std::size_t>::max()
            / sizeof(SceneFluidSurfaceVertexState)) {
        throw std::length_error(
            "regional opening wall post-step surface storage overflows");
    }
    return state.vertices.size() * sizeof(SceneFluidSurfaceVertexState);
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry) {
    return checkedAdd(
        surfaceStateStorageBytes(geometry.surfaceState),
        geometry.gridEpoch.ownedStorageBytes);
}

std::uint64_t geometryFingerprint(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry) {
    Fingerprint fingerprint;
    fingerprint.integer(geometry.version);
    for (const std::uint64_t value : {
             geometry.sourceCoupledStateFingerprint,
             geometry.sourceStructureStepFingerprint,
             geometry.sourcePostStepCheckpointFingerprint,
             geometry.sourceSurfaceStateFingerprint,
             geometry.surfaceDefinitionFingerprint,
             geometry.structureDefinitionFingerprint,
             geometry.previousAcceptedStepCount,
             geometry.currentAcceptedStepCount,
             geometry.surfaceState.fingerprint,
             geometry.gridEpoch.fingerprint}) {
        fingerprint.integer(value);
    }
    fingerprint.real(geometry.previousSimulationTimeSeconds);
    fingerprint.real(geometry.currentSimulationTimeSeconds);
    fingerprint.integer(
        static_cast<std::uint64_t>(geometry.ownedStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits&
        limits) {
    if (limits.maximumSurfaceVertices == 0
        || limits.maximumOwnedBytes == 0
        || limits.structureCheckpoint.maximumEncodedBytes == 0
        || limits.structureCheckpoint.maximumNodes == 0) {
        throw std::invalid_argument(
            "regional opening wall post-step geometry limits are invalid");
    }
}

void validateOwnedLimits(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits&
        limits) {
    if (geometry.surfaceState.vertices.size()
            > limits.maximumSurfaceVertices
        || geometry.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "regional opening wall post-step geometry exceeds its aggregate limit");
    }
}

std::vector<std::uint8_t> encodeStructure(
    const Structure& structure,
    const StructureCheckpointPersistenceLimits& limits) {
    std::vector<std::uint8_t> bytes;
    StructureCheckpointPersistenceError error;
    if (!serializeStructureCheckpoint(
            structure, structure.checkpoint(), bytes, &error, limits)) {
        const std::string message =
            "cannot bind regional opening wall post-step Structure: "
            + error.message;
        if (error.code
            == StructureCheckpointPersistenceErrorCode::LimitExceeded) {
            throw std::length_error(message);
        }
        throw std::invalid_argument(message);
    }
    return bytes;
}

void validateLiveStructure(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const Structure& structure,
    const StructureCheckpointPersistenceLimits& limits) {
    const auto liveBytes = encodeStructure(structure, limits);
    if (liveBytes
        != coupledState.structureStep.afterStructureCheckpoint) {
        throw std::invalid_argument(
            "regional opening wall post-step Structure is not the accepted coupled endpoint");
    }
}

void validateSourceGeometryBinding(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure) {
    const auto& loadEpoch = coupledState.structureStep.loadEpoch;
    if (surface.fingerprint
            != loadEpoch.pressureLoad.surfaceDefinitionFingerprint
        || transfer.surfaceDefinitionFingerprint() != surface.fingerprint
        || transfer.couplingSurfaceFingerprint()
            != loadEpoch.couplingSurfaceFingerprint
        || structure.definitionFingerprint()
            != coupledState.structureStep.targetDefinitionFingerprint
        || transfer.targetDefinitionFingerprint()
            != structure.definitionFingerprint()) {
        throw std::invalid_argument(
            "regional opening wall post-step geometry sources are foreign");
    }
}

SceneFluidRegionalOpeningMomentumWallPostStepGeometry buildGeometry(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits&
        limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
        coupledState);
    validateLiveStructure(
        coupledState, structure, limits.structureCheckpoint);
    validateSourceGeometryBinding(
        coupledState, surface, transfer, structure);
    if (surface.vertices.size() > limits.maximumSurfaceVertices) {
        throw std::length_error(
            "regional opening wall post-step surface exceeds its vertex limit");
    }

    SceneFluidRegionalOpeningMomentumWallPostStepGeometry result;
    result.sourceCoupledStateFingerprint = coupledState.fingerprint;
    result.sourceStructureStepFingerprint =
        coupledState.structureStep.fingerprint;
    result.sourcePostStepCheckpointFingerprint =
        coupledState.structureStep.afterCheckpointFingerprint;
    result.sourceSurfaceStateFingerprint =
        coupledState.structureStep.loadEpoch.surfaceStateFingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.structureDefinitionFingerprint =
        structure.definitionFingerprint();
    result.previousAcceptedStepCount =
        coupledState.structureStep.beforeAcceptedStepCount;
    result.currentAcceptedStepCount =
        coupledState.structureStep.afterAcceptedStepCount;
    result.previousSimulationTimeSeconds =
        coupledState.structureStep.beforeSimulationTimeSeconds;
    result.currentSimulationTimeSeconds =
        coupledState.structureStep.afterSimulationTimeSeconds;
    result.surfaceState = captureSceneFluidSurfaceState(
        surface, structureMappings, structure);
    if (result.surfaceState.acceptedStepCount
            != result.currentAcceptedStepCount
        || result.surfaceState.simulationTimeSeconds
            != result.currentSimulationTimeSeconds
        || result.surfaceState.definitionFingerprint
            != result.surfaceDefinitionFingerprint
        || result.surfaceState.structureDefinitionFingerprint
            != result.structureDefinitionFingerprint) {
        throw std::invalid_argument(
            "regional opening wall post-step surface epoch is inconsistent");
    }
    result.gridEpoch = buildSceneFluidGridEpoch(
        surface, result.surfaceState, grid, transfer, {},
        limits.gridEpoch);
    result.ownedStorageBytes = ownedStorageBytes(result);
    validateOwnedLimits(result, limits);
    result.fingerprint = geometryFingerprint(result);
    return result;
}

} // namespace

SceneFluidRegionalOpeningMomentumWallPostStepGeometry
buildSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits& limits) {
    auto result = buildGeometry(
        coupledState, surface, structureMappings, grid, transfer,
        structure, limits);
    validateSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
        result, coupledState, surface, structureMappings, grid, transfer,
        structure, limits);
    return result;
}

void validateSceneFluidRegionalOpeningMomentumWallPostStepGeometryIntegrity(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry) {
    validateSceneFluidSurfaceState(geometry.surfaceState);
    if (geometry.version
            != sceneFluidRegionalOpeningMomentumWallPostStepGeometryVersion
        || geometry.fingerprint == 0
        || geometry.sourceCoupledStateFingerprint == 0
        || geometry.sourceStructureStepFingerprint == 0
        || geometry.sourcePostStepCheckpointFingerprint == 0
        || geometry.sourceSurfaceStateFingerprint == 0
        || geometry.surfaceDefinitionFingerprint == 0
        || geometry.structureDefinitionFingerprint == 0
        || geometry.previousAcceptedStepCount
            == std::numeric_limits<std::uint64_t>::max()
        || geometry.currentAcceptedStepCount
            != geometry.previousAcceptedStepCount + 1
        || !std::isfinite(geometry.previousSimulationTimeSeconds)
        || !std::isfinite(geometry.currentSimulationTimeSeconds)
        || geometry.previousSimulationTimeSeconds < 0.0
        || !(geometry.currentSimulationTimeSeconds
            > geometry.previousSimulationTimeSeconds)
        || geometry.surfaceState.fingerprint == 0
        || geometry.gridEpoch.fingerprint == 0
        || geometry.surfaceState.acceptedStepCount
            != geometry.currentAcceptedStepCount
        || geometry.surfaceState.simulationTimeSeconds
            != geometry.currentSimulationTimeSeconds
        || geometry.surfaceState.definitionFingerprint
            != geometry.surfaceDefinitionFingerprint
        || geometry.surfaceState.structureDefinitionFingerprint
            != geometry.structureDefinitionFingerprint
        || geometry.ownedStorageBytes != ownedStorageBytes(geometry)
        || geometry.fingerprint != geometryFingerprint(geometry)) {
        throw std::invalid_argument(
            "regional opening wall post-step geometry integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallPostStepGeometry(
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometry& geometry,
    const SceneFluidRegionalOpeningMomentumWallCoupledState& coupledState,
    const SceneFluidSurfaceDefinition& surface,
    const SceneStructureMappings& structureMappings,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const Structure& structure,
    const SceneFluidRegionalOpeningMomentumWallPostStepGeometryLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallPostStepGeometryIntegrity(
        geometry);
    validateSceneFluidRegionalOpeningMomentumWallCoupledStateIntegrity(
        coupledState);
    validateLiveStructure(
        coupledState, structure, limits.structureCheckpoint);
    validateSceneFluidSurfaceState(surface, geometry.surfaceState);
    validateSceneFluidGridEpoch(
        geometry.gridEpoch, surface, geometry.surfaceState, grid, transfer);
    validateOwnedLimits(geometry, limits);
    const auto expected = buildGeometry(
        coupledState, surface, structureMappings, grid, transfer,
        structure, limits);
    if (geometry.sourceCoupledStateFingerprint
            != coupledState.fingerprint
        || geometry.sourceStructureStepFingerprint
            != coupledState.structureStep.fingerprint
        || geometry.sourcePostStepCheckpointFingerprint
            != coupledState.structureStep.afterCheckpointFingerprint
        || geometry.sourceSurfaceStateFingerprint
            != coupledState.structureStep.loadEpoch.surfaceStateFingerprint
        || geometry.surfaceDefinitionFingerprint != surface.fingerprint
        || geometry.structureDefinitionFingerprint
            != structure.definitionFingerprint()
        || geometry.surfaceState != expected.surfaceState
        || geometry.gridEpoch != expected.gridEpoch
        || geometry != expected) {
        throw std::invalid_argument(
            "regional opening wall post-step geometry sources are foreign");
    }
}

} // namespace simwing::fsi
