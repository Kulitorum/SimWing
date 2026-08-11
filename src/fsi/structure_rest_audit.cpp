#include "structure_rest_audit.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace simwing::fsi {
namespace {

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

bool zero(const StructureVector3& value) {
    return value == StructureVector3{};
}

double length(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

StructureVector3 difference(const StructureVector3& first,
                            const StructureVector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

double quaternionDistance(const StructureQuaternion& first,
                          const StructureQuaternion& second) {
    const double firstNorm = std::sqrt(
        first.w * first.w + first.x * first.x
        + first.y * first.y + first.z * first.z);
    const double secondNorm = std::sqrt(
        second.w * second.w + second.x * second.x
        + second.y * second.y + second.z * second.z);
    if (!std::isfinite(firstNorm) || !(firstNorm > 0.0)
        || !std::isfinite(secondNorm) || !(secondNorm > 0.0)) {
        throw std::logic_error(
            "Structure rest audit found an invalid payload orientation");
    }
    const double dot = (first.w * second.w + first.x * second.x
                        + first.y * second.y + first.z * second.z)
        / (firstNorm * secondNorm);
    return 2.0 * std::acos(std::clamp(std::abs(dot), 0.0, 1.0));
}

void validateSettings(const StructureStepSettings& settings,
                      const StructureRestAuditLimits& limits) {
    if (!std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || settings.substeps <= 0
        || settings.constraintIterations < 0
        || settings.cableConstraintSweepPairs < 0
        || !finite(settings.gravityMetersPerSecondSquared)
        || !zero(settings.gravityMetersPerSecondSquared)
        || !std::isfinite(settings.velocityDampingPerSecond)
        || settings.velocityDampingPerSecond < 0.0
        || !finite(settings.dampingReferenceVelocityMetersPerSecond)
        || !zero(settings.dampingReferenceVelocityMetersPerSecond)
        || !std::isfinite(limits.maximumLinearDisplacementMeters)
        || limits.maximumLinearDisplacementMeters < 0.0
        || !std::isfinite(limits.maximumLinearSpeedMetersPerSecond)
        || limits.maximumLinearSpeedMetersPerSecond < 0.0
        || !std::isfinite(limits.maximumAngularDisplacementRadians)
        || limits.maximumAngularDisplacementRadians < 0.0
        || !std::isfinite(limits.maximumAngularSpeedRadiansPerSecond)
        || limits.maximumAngularSpeedRadiansPerSecond < 0.0) {
        throw std::invalid_argument(
            "Structure rest-audit settings are invalid");
    }
}

void validateInitialState(
    const StructureCheckpoint& checkpoint,
    const std::optional<StructureSuspensionState>& suspension) {
    for (const auto& force : checkpoint.pendingExternalForcesNewtons) {
        if (!zero(force)) {
            throw std::invalid_argument(
                "Structure rest audit requires zero pending loads");
        }
    }
    for (const auto& node : checkpoint.nodes) {
        if (node.positionMeters != node.previousPositionMeters
            || !zero(node.velocityMetersPerSecond)) {
            throw std::invalid_argument(
                "Structure rest audit requires stationary node state");
        }
    }
    if (!suspension) {
        return;
    }
    if (!zero(suspension->payload.linearVelocityMetersPerSecond)
        || !zero(suspension->payload.angularVelocityRadiansPerSecond)) {
        throw std::invalid_argument(
            "Structure rest audit requires a stationary payload");
    }
    for (const auto& velocity : suspension->harnessVelocitiesMetersPerSecond) {
        if (!zero(velocity)) {
            throw std::invalid_argument(
                "Structure rest audit requires stationary harness points");
        }
    }
}

} // namespace

StructureRestAuditDiagnostics auditStructureRestState(
    Structure& structure,
    const StructureStepSettings& stepSettings,
    const StructureRestAuditLimits& limits) {
    validateSettings(stepSettings, limits);
    const auto before = structure.checkpoint();
    const auto suspensionBefore = structure.suspensionState();
    validateInitialState(before, suspensionBefore);

    try {
        const auto stepDiagnostics = structure.step(stepSettings);
        if (!stepDiagnostics.finite) {
            throw std::runtime_error(
                "Structure rest audit produced non-finite diagnostics");
        }
        const auto after = structure.nodeStates();
        const auto suspensionAfter = structure.suspensionState();
        if (after.size() != before.nodes.size()
            || suspensionAfter.has_value() != suspensionBefore.has_value()) {
            throw std::logic_error(
                "Structure rest audit changed state topology");
        }

        StructureRestAuditDiagnostics result;
        result.nodeCount = after.size();
        long double squaredDisplacementSum = 0.0L;
        for (std::size_t node = 0; node < after.size(); ++node) {
            const double displacement = length(difference(
                after[node].positionMeters,
                before.nodes[node].positionMeters));
            squaredDisplacementSum += static_cast<long double>(displacement)
                * static_cast<long double>(displacement);
            if (displacement > result.maximumNodeDisplacementMeters) {
                result.maximumNodeDisplacementMeters = displacement;
                result.maximumDisplacementNode = node;
            }
            result.maximumNodeSpeedMetersPerSecond = std::max(
                result.maximumNodeSpeedMetersPerSecond,
                length(after[node].velocityMetersPerSecond));
        }
        if (!after.empty()) {
            result.rmsNodeDisplacementMeters = std::sqrt(
                static_cast<double>(squaredDisplacementSum
                    / static_cast<long double>(after.size())));
        }

        if (suspensionBefore) {
            if (suspensionBefore->harnessPositionsMeters.size()
                    != suspensionAfter->harnessPositionsMeters.size()
                || suspensionBefore->harnessVelocitiesMetersPerSecond.size()
                    != suspensionAfter->harnessVelocitiesMetersPerSecond.size()) {
                throw std::logic_error(
                    "Structure rest audit changed suspension topology");
            }
            result.payloadDisplacementMeters = length(difference(
                suspensionAfter->payload.centerOfMassWorldMeters,
                suspensionBefore->payload.centerOfMassWorldMeters));
            result.payloadRotationRadians = quaternionDistance(
                suspensionAfter->payload.bodyToWorld,
                suspensionBefore->payload.bodyToWorld);
            result.maximumPayloadOrHarnessSpeedMetersPerSecond = length(
                suspensionAfter->payload.linearVelocityMetersPerSecond);
            result.payloadAngularSpeedRadiansPerSecond = length(
                suspensionAfter->payload.angularVelocityRadiansPerSecond);
            for (std::size_t point = 0;
                 point < suspensionAfter->harnessPositionsMeters.size();
                 ++point) {
                result.maximumHarnessDisplacementMeters = std::max(
                    result.maximumHarnessDisplacementMeters,
                    length(difference(
                        suspensionAfter->harnessPositionsMeters[point],
                        suspensionBefore->harnessPositionsMeters[point])));
                result.maximumPayloadOrHarnessSpeedMetersPerSecond = std::max(
                    result.maximumPayloadOrHarnessSpeedMetersPerSecond,
                    length(suspensionAfter->harnessVelocitiesMetersPerSecond[
                        point]));
            }
        }
        result.maximumSuspensionResidualMeters =
            stepDiagnostics.maximumSuspensionResidualMeters;
        result.stationary =
            result.maximumNodeDisplacementMeters
                <= limits.maximumLinearDisplacementMeters
            && result.payloadDisplacementMeters
                <= limits.maximumLinearDisplacementMeters
            && result.maximumHarnessDisplacementMeters
                <= limits.maximumLinearDisplacementMeters
            && result.maximumNodeSpeedMetersPerSecond
                <= limits.maximumLinearSpeedMetersPerSecond
            && result.maximumPayloadOrHarnessSpeedMetersPerSecond
                <= limits.maximumLinearSpeedMetersPerSecond
            && result.payloadRotationRadians
                <= limits.maximumAngularDisplacementRadians
            && result.payloadAngularSpeedRadiansPerSecond
                <= limits.maximumAngularSpeedRadiansPerSecond;
        structure.restore(before);
        return result;
    } catch (...) {
        structure.restore(before);
        throw;
    }
}

} // namespace simwing::fsi
