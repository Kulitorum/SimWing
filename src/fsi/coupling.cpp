#include "coupling.h"

#include <cmath>
#include <stdexcept>

namespace simwing::fsi {
namespace {

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
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

StructureVector3 scale(const StructureVector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double norm(const StructureVector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

} // namespace

std::uint32_t TimeIntegratedTransferResult::version() const noexcept {
    return version_;
}

std::uint64_t
TimeIntegratedTransferResult::surfaceFingerprint() const noexcept {
    return surfaceFingerprint_;
}

std::uint64_t
TimeIntegratedTransferResult::targetDefinitionFingerprint() const noexcept {
    return targetDefinitionFingerprint_;
}

std::span<const CouplingNodeImpulse>
TimeIntegratedTransferResult::nodeImpulses() const noexcept {
    return nodeImpulses_;
}

const TimeIntegratedTransferDiagnostics&
TimeIntegratedTransferResult::diagnostics() const noexcept {
    return diagnostics_;
}

ConservativeMacroStepCoupling::ConservativeMacroStepCoupling(
    const ConservativeSurfaceTransfer& transfer)
    : surfaceFingerprint_(transfer.fingerprint()),
      targetDefinitionFingerprint_(transfer.targetDefinitionFingerprint()),
      nodes_(transfer.nodes().begin(), transfer.nodes().end()) {}

std::uint64_t
ConservativeMacroStepCoupling::surfaceFingerprint() const noexcept {
    return surfaceFingerprint_;
}

std::uint64_t
ConservativeMacroStepCoupling::targetDefinitionFingerprint() const noexcept {
    return targetDefinitionFingerprint_;
}

TimeIntegratedTransferResult ConservativeMacroStepCoupling::integrate(
    const std::span<const double> sampleOffsetsSeconds,
    const std::span<const ConservativeTransferResult> samples) const {
    if (sampleOffsetsSeconds.size() < 2
        || sampleOffsetsSeconds.size() != samples.size()) {
        throw std::invalid_argument(
            "macro-step transfer requires at least two equally counted samples");
    }
    if (sampleOffsetsSeconds.front() != 0.0) {
        throw std::invalid_argument(
            "macro-step transfer sample offsets must start exactly at zero");
    }
    for (std::size_t index = 0; index < sampleOffsetsSeconds.size(); ++index) {
        if (!std::isfinite(sampleOffsetsSeconds[index])
            || (index != 0
                && !(sampleOffsetsSeconds[index]
                     > sampleOffsetsSeconds[index - 1]))) {
            throw std::invalid_argument(
                "macro-step transfer sample offsets must be finite and increasing");
        }
    }

    const StructureVector3 momentReference =
        samples.front().diagnostics().momentReferenceMeters;
    for (const ConservativeTransferResult& sample : samples) {
        const auto loads = sample.nodeLoads();
        const auto& diagnostics = sample.diagnostics();
        if (sample.surfaceFingerprint() != surfaceFingerprint_
            || sample.targetDefinitionFingerprint()
                != targetDefinitionFingerprint_
            || loads.size() != nodes_.size()
            || !diagnostics.finite
            || diagnostics.momentReferenceMeters != momentReference) {
            throw std::invalid_argument(
                "instantaneous transfer sample does not belong to this macro-step");
        }
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            if (loads[node].stableId != nodes_[node].stableId
                || loads[node].structureNode != nodes_[node].structureNode
                || !finite(loads[node].forceNewtons)) {
                throw std::invalid_argument(
                    "instantaneous transfer sample has invalid nodal bindings");
            }
        }
    }

    TimeIntegratedTransferResult result;
    result.surfaceFingerprint_ = surfaceFingerprint_;
    result.targetDefinitionFingerprint_ = targetDefinitionFingerprint_;
    result.nodeImpulses_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        result.nodeImpulses_.push_back(
            {node.stableId, node.structureNode, {}});
    }
    auto& diagnostics = result.diagnostics_;
    diagnostics.sampleCount = samples.size();
    diagnostics.intervalCount = samples.size() - 1;
    diagnostics.durationSeconds = sampleOffsetsSeconds.back();
    diagnostics.momentReferenceMeters = momentReference;

    for (std::size_t interval = 0;
         interval + 1 < samples.size(); ++interval) {
        const double halfInterval = 0.5
            * (sampleOffsetsSeconds[interval + 1]
               - sampleOffsetsSeconds[interval]);
        const auto& firstDiagnostics = samples[interval].diagnostics();
        const auto& secondDiagnostics = samples[interval + 1].diagnostics();
        diagnostics.integratedSurfaceImpulseNewtonSeconds = add(
            diagnostics.integratedSurfaceImpulseNewtonSeconds,
            scale(add(firstDiagnostics.integratedSurfaceForceNewtons,
                      secondDiagnostics.integratedSurfaceForceNewtons),
                  halfInterval));
        diagnostics.integratedSurfaceAngularImpulseNewtonMeterSeconds = add(
            diagnostics.integratedSurfaceAngularImpulseNewtonMeterSeconds,
            scale(add(
                      firstDiagnostics.integratedSurfaceMomentNewtonMeters,
                      secondDiagnostics.integratedSurfaceMomentNewtonMeters),
                  halfInterval));
        diagnostics.transferredNodalAngularImpulseNewtonMeterSeconds = add(
            diagnostics.transferredNodalAngularImpulseNewtonMeterSeconds,
            scale(add(
                      firstDiagnostics.transferredNodalMomentNewtonMeters,
                      secondDiagnostics.transferredNodalMomentNewtonMeters),
                  halfInterval));
        diagnostics.integratedSurfaceWorkJoules += halfInterval
            * (firstDiagnostics.integratedSurfacePowerWatts
               + secondDiagnostics.integratedSurfacePowerWatts);
        diagnostics.transferredNodalWorkJoules += halfInterval
            * (firstDiagnostics.transferredNodalPowerWatts
               + secondDiagnostics.transferredNodalPowerWatts);

        const auto firstLoads = samples[interval].nodeLoads();
        const auto secondLoads = samples[interval + 1].nodeLoads();
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            result.nodeImpulses_[node].impulseNewtonSeconds = add(
                result.nodeImpulses_[node].impulseNewtonSeconds,
                scale(add(firstLoads[node].forceNewtons,
                          secondLoads[node].forceNewtons),
                      halfInterval));
        }
    }

    for (const auto& node : result.nodeImpulses_) {
        diagnostics.transferredNodalImpulseNewtonSeconds = add(
            diagnostics.transferredNodalImpulseNewtonSeconds,
            node.impulseNewtonSeconds);
    }
    diagnostics.impulseResidualNewtonSeconds = subtract(
        diagnostics.transferredNodalImpulseNewtonSeconds,
        diagnostics.integratedSurfaceImpulseNewtonSeconds);
    diagnostics.impulseResidualNormNewtonSeconds = norm(
        diagnostics.impulseResidualNewtonSeconds);
    diagnostics.angularImpulseResidualNewtonMeterSeconds = subtract(
        diagnostics.transferredNodalAngularImpulseNewtonMeterSeconds,
        diagnostics.integratedSurfaceAngularImpulseNewtonMeterSeconds);
    diagnostics.angularImpulseResidualNormNewtonMeterSeconds = norm(
        diagnostics.angularImpulseResidualNewtonMeterSeconds);
    diagnostics.workResidualJoules =
        diagnostics.transferredNodalWorkJoules
        - diagnostics.integratedSurfaceWorkJoules;
    diagnostics.finite = std::isfinite(diagnostics.durationSeconds)
        && diagnostics.durationSeconds > 0.0
        && finite(diagnostics.momentReferenceMeters)
        && finite(diagnostics.integratedSurfaceImpulseNewtonSeconds)
        && finite(diagnostics.transferredNodalImpulseNewtonSeconds)
        && finite(diagnostics.impulseResidualNewtonSeconds)
        && std::isfinite(diagnostics.impulseResidualNormNewtonSeconds)
        && finite(
            diagnostics.integratedSurfaceAngularImpulseNewtonMeterSeconds)
        && finite(
            diagnostics.transferredNodalAngularImpulseNewtonMeterSeconds)
        && finite(diagnostics.angularImpulseResidualNewtonMeterSeconds)
        && std::isfinite(
            diagnostics.angularImpulseResidualNormNewtonMeterSeconds)
        && std::isfinite(diagnostics.integratedSurfaceWorkJoules)
        && std::isfinite(diagnostics.transferredNodalWorkJoules)
        && std::isfinite(diagnostics.workResidualJoules);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "time-integrated interface transfer produced non-finite ledgers");
    }
    return result;
}

StructureDiagnostics ConservativeMacroStepCoupling::advanceStructure(
    Structure& target,
    const TimeIntegratedTransferResult& transfer,
    const StructureStepSettings& settings) const {
    const auto impulses = transfer.nodeImpulses();
    const auto& diagnostics = transfer.diagnostics();
    if (transfer.version() != interfaceImpulseExchangeVersion
        || target.definitionFingerprint() != targetDefinitionFingerprint_
        || transfer.surfaceFingerprint() != surfaceFingerprint_
        || transfer.targetDefinitionFingerprint()
            != targetDefinitionFingerprint_
        || impulses.size() != nodes_.size()
        || !diagnostics.finite
        || settings.timeStepSeconds != diagnostics.durationSeconds) {
        throw std::invalid_argument(
            "time-integrated transfer does not match this structure macro-step");
    }

    std::vector<StructureVector3> averageForces;
    averageForces.reserve(impulses.size());
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        if (impulses[node].stableId != nodes_[node].stableId
            || impulses[node].structureNode != nodes_[node].structureNode
            || !finite(impulses[node].impulseNewtonSeconds)) {
            throw std::invalid_argument(
                "time-integrated transfer has invalid nodal bindings");
        }
        const auto averageForce = scale(
            impulses[node].impulseNewtonSeconds,
            1.0 / diagnostics.durationSeconds);
        if (!finite(averageForce)) {
            throw std::overflow_error(
                "time-integrated transfer average load is non-finite");
        }
        averageForces.push_back(averageForce);
    }

    const StructureCheckpoint before = target.checkpoint();
    try {
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            target.addExternalForce(
                nodes_[node].structureNode, averageForces[node]);
        }
        return target.step(settings);
    } catch (...) {
        target.restore(before);
        throw;
    }
}

} // namespace simwing::fsi
