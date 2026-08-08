#include "coupling.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

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

bool valid(const AitkenRelaxationSettings& settings) {
    return std::isfinite(settings.initialRelaxation)
        && std::isfinite(settings.minimumRelaxation)
        && std::isfinite(settings.maximumRelaxation)
        && settings.minimumRelaxation > 0.0
        && settings.initialRelaxation >= settings.minimumRelaxation
        && settings.initialRelaxation <= settings.maximumRelaxation;
}

bool finite(const std::span<const double> values) {
    for (const double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

} // namespace

AitkenInterfaceRelaxation::AitkenInterfaceRelaxation(
    const std::uint64_t interfaceDefinitionFingerprint,
    const std::size_t valueCount,
    const AitkenRelaxationSettings& settings)
    : interfaceDefinitionFingerprint_(interfaceDefinitionFingerprint),
      valueCount_(valueCount),
      settings_(settings),
      relaxation_(settings.initialRelaxation) {
    if (interfaceDefinitionFingerprint_ == 0
        || valueCount_ == 0 || !valid(settings_)) {
        throw std::invalid_argument(
            "Aitken relaxation identity, size, and factor bounds must be valid");
    }
}

std::uint64_t AitkenInterfaceRelaxation::
interfaceDefinitionFingerprint() const noexcept {
    return interfaceDefinitionFingerprint_;
}

std::size_t AitkenInterfaceRelaxation::valueCount() const noexcept {
    return valueCount_;
}

const AitkenRelaxationSettings&
AitkenInterfaceRelaxation::settings() const noexcept {
    return settings_;
}

std::uint64_t
AitkenInterfaceRelaxation::completedIterationCount() const noexcept {
    return completedIterationCount_;
}

double AitkenInterfaceRelaxation::relaxation() const noexcept {
    return relaxation_;
}

AitkenRelaxationCheckpoint
AitkenInterfaceRelaxation::checkpoint() const {
    AitkenRelaxationCheckpoint result;
    result.interfaceDefinitionFingerprint = interfaceDefinitionFingerprint_;
    result.valueCount = valueCount_;
    result.settings = settings_;
    result.completedIterationCount = completedIterationCount_;
    result.relaxation = relaxation_;
    result.previousResidual = previousResidual_;
    return result;
}

void AitkenInterfaceRelaxation::restore(
    const AitkenRelaxationCheckpoint& checkpoint) {
    const bool initial = checkpoint.completedIterationCount == 0;
    if (checkpoint.version != aitkenRelaxationCheckpointVersion
        || checkpoint.interfaceDefinitionFingerprint
            != interfaceDefinitionFingerprint_
        || checkpoint.valueCount != valueCount_
        || checkpoint.settings != settings_
        || !std::isfinite(checkpoint.relaxation)
        || checkpoint.relaxation < settings_.minimumRelaxation
        || checkpoint.relaxation > settings_.maximumRelaxation
        || (initial
            && (checkpoint.relaxation != settings_.initialRelaxation
                || !checkpoint.previousResidual.empty()))
        || (!initial
            && (checkpoint.previousResidual.size() != valueCount_
                || !finite(checkpoint.previousResidual)))) {
        throw std::invalid_argument(
            "Aitken relaxation checkpoint is incompatible or invalid");
    }

    std::vector<double> restoredResidual = checkpoint.previousResidual;
    completedIterationCount_ = checkpoint.completedIterationCount;
    relaxation_ = checkpoint.relaxation;
    previousResidual_.swap(restoredResidual);
}

void AitkenInterfaceRelaxation::reset() noexcept {
    completedIterationCount_ = 0;
    relaxation_ = settings_.initialRelaxation;
    previousResidual_.clear();
}

AitkenRelaxationDiagnostics AitkenInterfaceRelaxation::relax(
    const std::span<const double> current,
    const std::span<const double> candidate,
    std::vector<double>& relaxed) {
    if (current.size() != valueCount_
        || candidate.size() != valueCount_
        || !finite(current) || !finite(candidate)) {
        throw std::invalid_argument(
            "Aitken relaxation vectors must match and contain finite values");
    }
    if (completedIterationCount_
        == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error(
            "Aitken relaxation iteration counter overflowed");
    }

    std::vector<double> residual(valueCount_);
    double residualSquared = 0.0;
    for (std::size_t index = 0; index < valueCount_; ++index) {
        residual[index] = candidate[index] - current[index];
        residualSquared += residual[index] * residual[index];
    }

    double nextRelaxation = relaxation_;
    double residualChangeSquared = 0.0;
    double denominator = 0.0;
    bool usedDynamic = false;
    bool clipped = false;
    if (completedIterationCount_ != 0) {
        double previousResidualSquared = 0.0;
        double numerator = 0.0;
        for (std::size_t index = 0; index < valueCount_; ++index) {
            const double change = residual[index] - previousResidual_[index];
            residualChangeSquared += change * change;
            numerator += previousResidual_[index] * change;
            previousResidualSquared +=
                previousResidual_[index] * previousResidual_[index];
        }
        denominator = residualChangeSquared;
        const double scaleSquared = std::max(
            residualSquared, previousResidualSquared);
        constexpr double reliableDifferenceFactor =
            64.0 * std::numeric_limits<double>::epsilon();
        if (denominator > reliableDifferenceFactor * scaleSquared) {
            const double unbounded = -relaxation_ * numerator / denominator;
            if (!std::isfinite(unbounded)) {
                throw std::overflow_error(
                    "Aitken relaxation factor became non-finite");
            }
            nextRelaxation = std::clamp(
                unbounded,
                settings_.minimumRelaxation,
                settings_.maximumRelaxation);
            usedDynamic = true;
            clipped = nextRelaxation != unbounded;
        }
    }

    std::vector<double> next(valueCount_);
    for (std::size_t index = 0; index < valueCount_; ++index) {
        next[index] = current[index] + nextRelaxation * residual[index];
    }
    const double residualL2 = std::sqrt(residualSquared);
    const double residualChangeL2 = std::sqrt(residualChangeSquared);
    if (!finite(next) || !std::isfinite(residualL2)
        || !std::isfinite(residualChangeL2)
        || !std::isfinite(denominator)) {
        throw std::overflow_error(
            "Aitken relaxation produced non-finite state");
    }

    relaxed.swap(next);
    previousResidual_ = std::move(residual);
    relaxation_ = nextRelaxation;
    ++completedIterationCount_;
    return {
        completedIterationCount_,
        relaxation_,
        residualL2,
        residualChangeL2,
        denominator,
        usedDynamic,
        clipped,
        true,
    };
}

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
