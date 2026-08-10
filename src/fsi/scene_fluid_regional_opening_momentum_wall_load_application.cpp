#include "scene_fluid_regional_opening_momentum_wall_load_application.h"

#include <algorithm>
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
        for (std::size_t index = 0; index < sizeof(value); ++index) {
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

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

StructureVector3 converted(const fluid::Vector3& value) {
    return {value.x, value.y, value.z};
}

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

StructureVector3 scaled(const StructureVector3& value,
                        const double scale) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

double maximumAbsolute(const StructureVector3& value) {
    return std::max({std::abs(value.x),
                     std::abs(value.y),
                     std::abs(value.z)});
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "regional opening wall-load application storage overflows");
    }
    return first * second;
}

void fingerprintVector(Fingerprint& fingerprint,
                       const StructureVector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void validateSettings(
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings) {
    if (!finite(settings.transfer.momentReferenceMeters)
        || !std::isfinite(
            settings.transfer.minimumTriangleAreaSquareMeters)
        || !(settings.transfer.minimumTriangleAreaSquareMeters > 0.0)
        || !std::isfinite(
            settings.transfer.minimumQuadratureAreaSquareMeters)
        || !(settings.transfer.minimumQuadratureAreaSquareMeters > 0.0)
        || !std::isfinite(settings.transfer.barycentricTolerance)
        || settings.transfer.barycentricTolerance < 0.0
        || !std::isfinite(
            settings
                .absoluteActionReactionToleranceKilogramMetersPerSecond)
        || settings
               .absoluteActionReactionToleranceKilogramMetersPerSecond
            < 0.0
        || !std::isfinite(settings.relativeActionReactionTolerance)
        || settings.relativeActionReactionTolerance < 0.0) {
        throw std::invalid_argument(
            "regional opening wall-load application settings are invalid");
    }
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits) {
    if (limits.maximumNodeLoads == 0
        || limits.maximumStructureNodes == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening wall-load application limits are invalid");
    }
}

std::uint64_t settingsFingerprint(
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings) {
    Fingerprint fingerprint;
    for (const double value : {
             settings.transfer.momentReferenceMeters.x,
             settings.transfer.momentReferenceMeters.y,
             settings.transfer.momentReferenceMeters.z,
             settings.transfer.minimumTriangleAreaSquareMeters,
             settings.transfer.minimumQuadratureAreaSquareMeters,
             settings.transfer.barycentricTolerance,
             settings
                 .absoluteActionReactionToleranceKilogramMetersPerSecond,
             settings.relativeActionReactionTolerance}) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application) {
    return checkedMultiply(
        application.nodeLoads.size(),
        sizeof(SceneFluidRegionalOpeningMomentumWallAppliedNodeLoad));
}

std::size_t workingStorageBytes(const std::size_t structureNodeCount) {
    return checkedMultiply(
        structureNodeCount, sizeof(StructureVector3));
}

std::uint64_t applicationFingerprint(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application) {
    Fingerprint fingerprint;
    fingerprint.integer(application.version);
    for (const std::uint64_t value : {
             application.sourceCycleStateFingerprint,
             application.sourceAdjustmentStateFingerprint,
             application.sourceWallTractionFingerprint,
             application.sourceWallExchangeFingerprint,
             application.transportMetricFingerprint,
             application.acceptedMetricFingerprint,
             application.quadratureFingerprint,
             application.surfaceStateFingerprint,
             application.couplingSurfaceFingerprint,
             application.targetDefinitionFingerprint,
             application.sourceSettingsFingerprint,
             application.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(application.simulationTimeSeconds);
    fingerprint.real(application.timeStepSeconds);
    fingerprint.integer(
        static_cast<std::uint64_t>(application.structureNodeCount));
    fingerprint.integer(
        static_cast<std::uint64_t>(application.nodeLoads.size()));
    for (const auto& load : application.nodeLoads) {
        fingerprint.integer(static_cast<std::uint64_t>(load.loadIndex));
        fingerprint.integer(load.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(load.structureNode));
        fingerprintVector(fingerprint, load.priorPendingForceNewtons);
        fingerprintVector(fingerprint, load.appliedWallForceNewtons);
        fingerprintVector(fingerprint, load.resultingPendingForceNewtons);
        fingerprintVector(fingerprint, load.applicationResidualNewtons);
    }
    fingerprintVector(
        fingerprint,
        application.fluidAdjustmentImpulseKilogramMetersPerSecond);
    fingerprintVector(fingerprint, application.transferredWallForceNewtons);
    fingerprintVector(
        fingerprint,
        application.structureWallImpulseKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint,
        application.actionReactionResidualKilogramMetersPerSecond);
    fingerprintVector(fingerprint, application.priorPendingForceNewtons);
    fingerprintVector(fingerprint, application.appliedWallForceNewtons);
    fingerprintVector(fingerprint, application.resultingPendingForceNewtons);
    fingerprintVector(fingerprint, application.applicationResidualNewtons);
    fingerprint.boolean(application.applied);
    fingerprint.integer(
        static_cast<std::uint64_t>(application.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(application.workingStorageBytes));
    return fingerprint.value();
}

double actionReactionTolerance(
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings,
    const StructureVector3& fluidImpulse,
    const StructureVector3& structureImpulse) {
    const double scale = std::max({
        maximumAbsolute(fluidImpulse),
        maximumAbsolute(structureImpulse), 1.0});
    return std::max(
        settings.absoluteActionReactionToleranceKilogramMetersPerSecond,
        settings.relativeActionReactionTolerance * scale);
}

void validateActionReaction(
    const StructureVector3& fluidImpulse,
    const StructureVector3& structureImpulse,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings) {
    const auto residual = add(fluidImpulse, structureImpulse);
    if (!finite(residual)
        || maximumAbsolute(residual)
            > actionReactionTolerance(
                settings, fluidImpulse, structureImpulse)) {
        throw std::invalid_argument(
            "regional opening wall load does not close action and reaction");
    }
}

ConservativeTransferResult evaluateWallTransfer(
    const SceneFluidRegionalOpeningMomentumWallCycleState& cycleState,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const ConservativeTransferSettings& settings) {
    validateSceneFluidAcceptedWallTractions(
        cycleState.wallTractions, quadrature,
        cycleState.sourceWallExchangeFingerprint);
    return evaluateSceneFluidQuadrature(
        transfer, surfaceState, quadrature,
        cycleState.wallTractions.tractions, settings);
}

void validateSourceState(
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
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallCycleStateLimits& limits) {
    validateSceneFluidRegionalOpeningMomentumWallCycleState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits);
}

} // namespace

SceneFluidRegionalOpeningMomentumWallLoadApplication
applySceneFluidRegionalOpeningMomentumWallLoads(
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
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    Structure& target,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    validateSourceState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits.cycleState);
    validateSceneFluidSurfaceState(surfaceState);

    const std::size_t nodeLoadCount = transfer.nodes().size();
    const std::size_t structureNodeCount = target.definition().nodes.size();
    if (nodeLoadCount == 0 || structureNodeCount == 0
        || nodeLoadCount > limits.maximumNodeLoads
        || structureNodeCount > limits.maximumStructureNodes) {
        throw std::length_error(
            "regional opening wall-load application count limit exceeded");
    }
    const std::size_t expectedOwnedBytes = checkedMultiply(
        nodeLoadCount,
        sizeof(SceneFluidRegionalOpeningMomentumWallAppliedNodeLoad));
    const std::size_t expectedWorkingBytes =
        workingStorageBytes(structureNodeCount);
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || expectedWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening wall-load application byte limit exceeded");
    }
    if (target.definitionFingerprint()
            != transfer.targetDefinitionFingerprint()
        || target.acceptedStepCount() != surfaceState.acceptedStepCount
        || target.simulationTimeSeconds()
            != surfaceState.simulationTimeSeconds) {
        throw std::invalid_argument(
            "regional opening wall-load target epoch is stale");
    }
    const auto sampledKinematics = transfer.kinematics(surfaceState);
    const auto targetKinematics =
        transfer.conservativeTransfer().captureKinematics(target);
    if (sampledKinematics != targetKinematics) {
        throw std::invalid_argument(
            "regional opening wall-load target kinematics are stale");
    }

    const auto transferred = evaluateWallTransfer(
        cycleState, surfaceState, transfer, quadrature,
        settings.transfer);
    const auto transferredLoads = transferred.nodeLoads();
    if (transferredLoads.size() != nodeLoadCount) {
        throw std::logic_error(
            "regional opening wall transfer changed node count");
    }
    const StructureVector3 fluidImpulse = converted(
        cycleState.adjustedMomentum
            .adjustmentImpulseKilogramMetersPerSecond);
    const StructureVector3 wallForce =
        transferred.diagnostics().integratedSurfaceForceNewtons;
    const double timeStepSeconds =
        cycleState.adjustedMomentum.timeStepSeconds;
    const StructureVector3 structureImpulse =
        scaled(wallForce, timeStepSeconds);
    validateActionReaction(fluidImpulse, structureImpulse, settings);

    const StructureCheckpoint before = target.checkpoint();
    if (before.definitionFingerprint != target.definitionFingerprint()
        || before.acceptedStepCount != surfaceState.acceptedStepCount
        || before.simulationTimeSeconds != surfaceState.simulationTimeSeconds
        || before.nodes.size() != structureNodeCount
        || before.pendingExternalForcesNewtons.size()
            != structureNodeCount) {
        throw std::invalid_argument(
            "regional opening wall-load target checkpoint is incompatible");
    }
    std::vector<StructureVector3> expectedPending =
        before.pendingExternalForcesNewtons;

    SceneFluidRegionalOpeningMomentumWallLoadApplication application;
    application.sourceCycleStateFingerprint = cycleState.fingerprint;
    application.sourceAdjustmentStateFingerprint =
        cycleState.adjustedMomentum.fingerprint;
    application.sourceWallTractionFingerprint =
        cycleState.wallTractions.fingerprint;
    application.sourceWallExchangeFingerprint =
        cycleState.sourceWallExchangeFingerprint;
    application.transportMetricFingerprint =
        cycleState.transportMetricFingerprint;
    application.acceptedMetricFingerprint =
        cycleState.acceptedMetricFingerprint;
    application.quadratureFingerprint = quadrature.fingerprint;
    application.surfaceStateFingerprint = surfaceState.fingerprint;
    application.couplingSurfaceFingerprint =
        transfer.couplingSurfaceFingerprint();
    application.targetDefinitionFingerprint =
        transfer.targetDefinitionFingerprint();
    application.sourceSettingsFingerprint = settingsFingerprint(settings);
    application.acceptedStepCount = surfaceState.acceptedStepCount;
    application.simulationTimeSeconds = surfaceState.simulationTimeSeconds;
    application.timeStepSeconds = timeStepSeconds;
    application.structureNodeCount = structureNodeCount;
    application.fluidAdjustmentImpulseKilogramMetersPerSecond =
        fluidImpulse;
    application.transferredWallForceNewtons = wallForce;
    application.structureWallImpulseKilogramMetersPerSecond =
        structureImpulse;
    application.actionReactionResidualKilogramMetersPerSecond =
        add(fluidImpulse, structureImpulse);
    application.nodeLoads.reserve(nodeLoadCount);
    for (const auto& force : before.pendingExternalForcesNewtons) {
        application.priorPendingForceNewtons = add(
            application.priorPendingForceNewtons, force);
    }
    for (std::size_t index = 0; index < nodeLoadCount; ++index) {
        const auto& load = transferredLoads[index];
        if (load.structureNode >= structureNodeCount
            || !finite(load.forceNewtons)) {
            throw std::invalid_argument(
                "regional opening wall node load is invalid");
        }
        const StructureVector3 prior = expectedPending[load.structureNode];
        const StructureVector3 resulting = add(prior, load.forceNewtons);
        if (!finite(resulting)) {
            throw std::overflow_error(
                "regional opening wall pending load is not finite");
        }
        expectedPending[load.structureNode] = resulting;
        application.appliedWallForceNewtons = add(
            application.appliedWallForceNewtons, load.forceNewtons);
        application.nodeLoads.push_back({
            index,
            load.stableId,
            load.structureNode,
            prior,
            load.forceNewtons,
            resulting,
            subtract(resulting, add(prior, load.forceNewtons)),
        });
    }
    for (const auto& force : expectedPending) {
        application.resultingPendingForceNewtons = add(
            application.resultingPendingForceNewtons, force);
    }
    application.applicationResidualNewtons = subtract(
        application.resultingPendingForceNewtons,
        add(application.priorPendingForceNewtons,
            application.appliedWallForceNewtons));
    application.ownedStorageBytes = ownedStorageBytes(application);
    application.workingStorageBytes = expectedWorkingBytes;
    if (application.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "regional opening wall-load application storage changed");
    }

    try {
        transfer.addLoadsTo(target, transferred);
        const StructureCheckpoint after = target.checkpoint();
        if (after.definitionFingerprint != before.definitionFingerprint
            || after.acceptedStepCount != before.acceptedStepCount
            || after.simulationTimeSeconds != before.simulationTimeSeconds
            || after.nodes != before.nodes
            || after.pendingExternalForcesNewtons != expectedPending
            || after.lastAppliedExternalForceNewtons
                != before.lastAppliedExternalForceNewtons) {
            throw std::logic_error(
                "regional opening wall-load application changed non-load state");
        }
        application.applied = true;
        application.fingerprint = applicationFingerprint(application);
        validateSceneFluidRegionalOpeningMomentumWallLoadApplication(
            application, cycleState, transportMetric,
            acceptedPressureOperator, acceptedBasePressureOperator, grid,
            acceptedSweep, acceptedFragments, acceptedTopology,
            acceptedVolumeRates, acceptedOpeningDefinitions,
            acceptedOpenings, acceptedResistanceDefinitions,
            acceptedBaseMetric, acceptedMetric, surfaceState, transfer,
            quadrature, settings, limits);
        return application;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalOpeningMomentumWallLoadApplicationIntegrity(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application) {
    if (application.version
            != sceneFluidRegionalOpeningMomentumWallLoadApplicationVersion
        || application.fingerprint == 0
        || application.sourceCycleStateFingerprint == 0
        || application.sourceAdjustmentStateFingerprint == 0
        || application.sourceWallTractionFingerprint == 0
        || application.sourceWallExchangeFingerprint == 0
        || application.transportMetricFingerprint == 0
        || application.acceptedMetricFingerprint == 0
        || application.quadratureFingerprint == 0
        || application.surfaceStateFingerprint == 0
        || application.couplingSurfaceFingerprint == 0
        || application.targetDefinitionFingerprint == 0
        || application.sourceSettingsFingerprint == 0
        || !std::isfinite(application.simulationTimeSeconds)
        || !std::isfinite(application.timeStepSeconds)
        || !(application.timeStepSeconds > 0.0)
        || application.structureNodeCount == 0
        || application.nodeLoads.empty()
        || !finite(
            application
                .fluidAdjustmentImpulseKilogramMetersPerSecond)
        || !finite(application.transferredWallForceNewtons)
        || !finite(
            application
                .structureWallImpulseKilogramMetersPerSecond)
        || !finite(
            application
                .actionReactionResidualKilogramMetersPerSecond)
        || !finite(application.priorPendingForceNewtons)
        || !finite(application.appliedWallForceNewtons)
        || !finite(application.resultingPendingForceNewtons)
        || !finite(application.applicationResidualNewtons)
        || !application.applied
        || application.ownedStorageBytes != ownedStorageBytes(application)
        || application.workingStorageBytes
            != workingStorageBytes(application.structureNodeCount)
        || application.structureWallImpulseKilogramMetersPerSecond
            != scaled(
                application.transferredWallForceNewtons,
                application.timeStepSeconds)
        || application.actionReactionResidualKilogramMetersPerSecond
            != add(
                application
                    .fluidAdjustmentImpulseKilogramMetersPerSecond,
                application
                    .structureWallImpulseKilogramMetersPerSecond)
        || application.fingerprint != applicationFingerprint(application)) {
        throw std::invalid_argument(
            "regional opening wall-load application integrity is invalid");
    }

    StructureVector3 appliedForce;
    std::uint64_t previousStableId = 0;
    std::size_t previousStructureNode = 0;
    bool havePrevious = false;
    for (std::size_t index = 0;
         index < application.nodeLoads.size(); ++index) {
        const auto& load = application.nodeLoads[index];
        const StructureVector3 reconstructed = add(
            load.priorPendingForceNewtons,
            load.appliedWallForceNewtons);
        const StructureVector3 residual = subtract(
            load.resultingPendingForceNewtons, reconstructed);
        if (load.loadIndex != index || load.stableId == 0
            || load.structureNode >= application.structureNodeCount
            || !finite(load.priorPendingForceNewtons)
            || !finite(load.appliedWallForceNewtons)
            || !finite(load.resultingPendingForceNewtons)
            || !finite(load.applicationResidualNewtons)
            || load.resultingPendingForceNewtons != reconstructed
            || load.applicationResidualNewtons != residual
            || (havePrevious
                && (load.stableId <= previousStableId
                    || load.structureNode <= previousStructureNode))) {
            throw std::invalid_argument(
                "regional opening wall applied node load is invalid");
        }
        havePrevious = true;
        previousStableId = load.stableId;
        previousStructureNode = load.structureNode;
        appliedForce = add(appliedForce, load.appliedWallForceNewtons);
    }
    const StructureVector3 aggregateResidual = subtract(
        application.resultingPendingForceNewtons,
        add(application.priorPendingForceNewtons,
            application.appliedWallForceNewtons));
    const double forceScale = std::max({
        maximumAbsolute(application.priorPendingForceNewtons),
        maximumAbsolute(application.appliedWallForceNewtons),
        maximumAbsolute(application.resultingPendingForceNewtons), 1.0});
    const double forceTolerance =
        16384.0 * std::numeric_limits<double>::epsilon() * forceScale;
    if (appliedForce != application.appliedWallForceNewtons
        || aggregateResidual != application.applicationResidualNewtons
        || maximumAbsolute(application.applicationResidualNewtons)
            > forceTolerance) {
        throw std::invalid_argument(
            "regional opening wall-load application closure is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallLoadApplication(
    const SceneFluidRegionalOpeningMomentumWallLoadApplication& application,
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
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationSettings&
        settings,
    const SceneFluidRegionalOpeningMomentumWallLoadApplicationLimits& limits) {
    validateSettings(settings);
    validateLimits(limits);
    validateSourceState(
        cycleState, transportMetric, acceptedPressureOperator,
        acceptedBasePressureOperator, grid, acceptedSweep,
        acceptedFragments, acceptedTopology, acceptedVolumeRates,
        acceptedOpeningDefinitions, acceptedOpenings,
        acceptedResistanceDefinitions, acceptedBaseMetric,
        acceptedMetric, quadrature, limits.cycleState);
    validateSceneFluidSurfaceState(surfaceState);
    validateSceneFluidRegionalOpeningMomentumWallLoadApplicationIntegrity(
        application);
    if (application.nodeLoads.size() > limits.maximumNodeLoads
        || application.structureNodeCount > limits.maximumStructureNodes
        || application.ownedStorageBytes > limits.maximumOwnedBytes
        || application.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening wall-load application validation limit exceeded");
    }
    if (application.sourceCycleStateFingerprint != cycleState.fingerprint
        || application.sourceAdjustmentStateFingerprint
            != cycleState.adjustedMomentum.fingerprint
        || application.sourceWallTractionFingerprint
            != cycleState.wallTractions.fingerprint
        || application.sourceWallExchangeFingerprint
            != cycleState.sourceWallExchangeFingerprint
        || application.transportMetricFingerprint
            != cycleState.transportMetricFingerprint
        || application.acceptedMetricFingerprint
            != cycleState.acceptedMetricFingerprint
        || application.quadratureFingerprint != quadrature.fingerprint
        || application.surfaceStateFingerprint != surfaceState.fingerprint
        || application.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || application.targetDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || application.sourceSettingsFingerprint
            != settingsFingerprint(settings)
        || application.acceptedStepCount != surfaceState.acceptedStepCount
        || application.simulationTimeSeconds
            != surfaceState.simulationTimeSeconds
        || application.timeStepSeconds
            != cycleState.adjustedMomentum.timeStepSeconds
        || application.nodeLoads.size() != transfer.nodes().size()) {
        throw std::invalid_argument(
            "regional opening wall-load application is foreign to its source");
    }

    const auto transferred = evaluateWallTransfer(
        cycleState, surfaceState, transfer, quadrature,
        settings.transfer);
    const StructureVector3 expectedFluidImpulse = converted(
        cycleState.adjustedMomentum
            .adjustmentImpulseKilogramMetersPerSecond);
    const StructureVector3 expectedWallForce =
        transferred.diagnostics().integratedSurfaceForceNewtons;
    const StructureVector3 expectedStructureImpulse = scaled(
        expectedWallForce,
        cycleState.adjustedMomentum.timeStepSeconds);
    validateActionReaction(
        expectedFluidImpulse, expectedStructureImpulse, settings);
    if (application.fluidAdjustmentImpulseKilogramMetersPerSecond
            != expectedFluidImpulse
        || application.transferredWallForceNewtons != expectedWallForce
        || application.structureWallImpulseKilogramMetersPerSecond
            != expectedStructureImpulse
        || application.actionReactionResidualKilogramMetersPerSecond
            != add(expectedFluidImpulse, expectedStructureImpulse)) {
        throw std::invalid_argument(
            "regional opening wall-load action/reaction ledger is foreign");
    }

    const auto expectedLoads = transferred.nodeLoads();
    if (expectedLoads.size() != application.nodeLoads.size()) {
        throw std::invalid_argument(
            "regional opening wall-load application changed node count");
    }
    for (std::size_t index = 0; index < expectedLoads.size(); ++index) {
        const auto& expected = expectedLoads[index];
        const auto& actual = application.nodeLoads[index];
        if (actual.loadIndex != index
            || actual.stableId != expected.stableId
            || actual.structureNode != expected.structureNode
            || actual.appliedWallForceNewtons != expected.forceNewtons) {
            throw std::invalid_argument(
                "regional opening wall applied node load is foreign to its source");
        }
    }
}

} // namespace simwing::fsi
