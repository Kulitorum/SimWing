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

bool finiteNonnegative(const double value) {
    return std::isfinite(value) && value >= 0.0;
}

bool valid(const CouplingConvergenceSettings& settings) {
    return settings.minimumIterations > 0
        && settings.minimumIterations <= settings.maximumIterations
        && finiteNonnegative(
            settings.absoluteDisplacementToleranceMetres)
        && finiteNonnegative(settings.relativeDisplacementTolerance)
        && std::isfinite(settings.displacementReferenceFloorMetres)
        && settings.displacementReferenceFloorMetres > 0.0
        && finiteNonnegative(
            settings.absoluteVelocityToleranceMetersPerSecond)
        && finiteNonnegative(settings.relativeVelocityTolerance)
        && std::isfinite(settings.velocityReferenceFloorMetersPerSecond)
        && settings.velocityReferenceFloorMetersPerSecond > 0.0
        && finiteNonnegative(settings.absoluteTractionToleranceNewtons)
        && finiteNonnegative(settings.relativeTractionTolerance)
        && std::isfinite(settings.tractionReferenceFloorNewtons)
        && settings.tractionReferenceFloorNewtons > 0.0;
}

} // namespace

CouplingConvergenceDecision evaluateCouplingConvergence(
    const std::uint64_t iteration,
    const CouplingResidualNorms& residuals,
    const CouplingConvergenceSettings& settings) {
    if (!valid(settings)) {
        throw std::invalid_argument(
            "coupling convergence iteration and tolerance settings are invalid");
    }
    if (iteration == 0 || iteration > settings.maximumIterations
        || !finiteNonnegative(residuals.displacementMetres)
        || !finiteNonnegative(residuals.displacementReferenceMetres)
        || !finiteNonnegative(residuals.velocityMetersPerSecond)
        || !finiteNonnegative(
            residuals.velocityReferenceMetersPerSecond)
        || !finiteNonnegative(residuals.tractionNewtons)
        || !finiteNonnegative(residuals.tractionReferenceNewtons)) {
        throw std::invalid_argument(
            "coupling convergence iteration and residual norms are invalid");
    }

    CouplingConvergenceDecision decision;
    decision.iteration = iteration;
    decision.residuals = residuals;
    decision.relativeDisplacement = residuals.displacementMetres
        / std::max(residuals.displacementReferenceMetres,
                   settings.displacementReferenceFloorMetres);
    decision.relativeVelocity = residuals.velocityMetersPerSecond
        / std::max(residuals.velocityReferenceMetersPerSecond,
                   settings.velocityReferenceFloorMetersPerSecond);
    decision.relativeTraction = residuals.tractionNewtons
        / std::max(residuals.tractionReferenceNewtons,
                   settings.tractionReferenceFloorNewtons);
    decision.finite = std::isfinite(decision.relativeDisplacement)
        && std::isfinite(decision.relativeVelocity)
        && std::isfinite(decision.relativeTraction);
    if (!decision.finite) {
        throw std::overflow_error(
            "coupling convergence relative residual overflowed");
    }
    decision.displacementConverged =
        residuals.displacementMetres
            <= settings.absoluteDisplacementToleranceMetres
        && decision.relativeDisplacement
            <= settings.relativeDisplacementTolerance;
    decision.velocityConverged =
        residuals.velocityMetersPerSecond
            <= settings.absoluteVelocityToleranceMetersPerSecond
        && decision.relativeVelocity
            <= settings.relativeVelocityTolerance;
    decision.tractionConverged =
        residuals.tractionNewtons
            <= settings.absoluteTractionToleranceNewtons
        && decision.relativeTraction
            <= settings.relativeTractionTolerance;
    decision.minimumIterationsSatisfied =
        iteration >= settings.minimumIterations;
    decision.converged = decision.minimumIterationsSatisfied
        && decision.displacementConverged
        && decision.velocityConverged
        && decision.tractionConverged;
    decision.iterationLimitReached =
        !decision.converged && iteration == settings.maximumIterations;
    return decision;
}

StrongCouplingIteration::StrongCouplingIteration(
    const std::uint64_t interfaceDefinitionFingerprint,
    const std::span<const double> initialInterface,
    const AitkenRelaxationSettings& relaxationSettings,
    const CouplingConvergenceSettings& convergenceSettings)
    : convergenceSettings_(convergenceSettings),
      relaxation_(interfaceDefinitionFingerprint,
                  initialInterface.size(),
                  relaxationSettings),
      currentInterface_(initialInterface.begin(), initialInterface.end()) {
    if (!valid(convergenceSettings_) || !finite(currentInterface_)) {
        throw std::invalid_argument(
            "strong-coupling initial interface or convergence settings are invalid");
    }
}

std::span<const double>
StrongCouplingIteration::currentInterface() const noexcept {
    return currentInterface_;
}

StrongCouplingIterationStatus
StrongCouplingIteration::status() const noexcept {
    return status_;
}

std::uint64_t
StrongCouplingIteration::completedIterationCount() const noexcept {
    return relaxation_.completedIterationCount();
}

const AitkenRelaxationDiagnostics&
StrongCouplingIteration::lastRelaxation() const noexcept {
    return lastRelaxation_;
}

const CouplingConvergenceDecision&
StrongCouplingIteration::lastConvergence() const noexcept {
    return lastConvergence_;
}

StrongCouplingIterationCheckpoint
StrongCouplingIteration::checkpoint() const {
    StrongCouplingIterationCheckpoint result;
    result.convergenceSettings = convergenceSettings_;
    result.relaxation = relaxation_.checkpoint();
    result.currentInterface = currentInterface_;
    result.status = status_;
    result.lastRelaxation = lastRelaxation_;
    result.lastConvergence = lastConvergence_;
    return result;
}

StrongCouplingIterationResult StrongCouplingIteration::advance(
    const std::span<const double> unrelaxedCandidate,
    const CouplingResidualNorms& residuals) {
    if (status_ != StrongCouplingIterationStatus::Iterating) {
        throw std::logic_error(
            "strong-coupling iteration is already terminal");
    }
    const std::uint64_t nextIteration =
        relaxation_.completedIterationCount() + 1;
    const CouplingConvergenceDecision convergence =
        evaluateCouplingConvergence(
            nextIteration, residuals, convergenceSettings_);
    std::vector<double> nextInterface;
    const AitkenRelaxationDiagnostics relaxation = relaxation_.relax(
        currentInterface_, unrelaxedCandidate, nextInterface);

    StrongCouplingIterationStatus nextStatus =
        StrongCouplingIterationStatus::Iterating;
    if (convergence.converged) {
        nextStatus = StrongCouplingIterationStatus::Converged;
    } else if (convergence.iterationLimitReached) {
        nextStatus = StrongCouplingIterationStatus::Exhausted;
    }
    currentInterface_.swap(nextInterface);
    status_ = nextStatus;
    lastRelaxation_ = relaxation;
    lastConvergence_ = convergence;
    return {status_, lastRelaxation_, lastConvergence_};
}

void StrongCouplingIteration::restore(
    const StrongCouplingIterationCheckpoint& checkpoint) {
    if (checkpoint.version != strongCouplingIterationCheckpointVersion
        || checkpoint.convergenceSettings != convergenceSettings_
        || checkpoint.currentInterface.size() != currentInterface_.size()
        || !finite(checkpoint.currentInterface)) {
        throw std::invalid_argument(
            "strong-coupling iteration checkpoint is incompatible or invalid");
    }

    const std::uint64_t iteration =
        checkpoint.relaxation.completedIterationCount;
    if (iteration == 0) {
        if (checkpoint.status != StrongCouplingIterationStatus::Iterating
            || checkpoint.lastRelaxation
                != AitkenRelaxationDiagnostics{}
            || checkpoint.lastConvergence
                != CouplingConvergenceDecision{}) {
            throw std::invalid_argument(
                "initial strong-coupling checkpoint has committed diagnostics");
        }
    } else {
        double storedResidualSquared = 0.0;
        for (const double value : checkpoint.relaxation.previousResidual) {
            storedResidualSquared += value * value;
        }
        const double storedResidualL2 = std::sqrt(storedResidualSquared);
        const double storedChangeL2 = std::sqrt(
            checkpoint.lastRelaxation.denominator);
        const auto& relaxationSettings = relaxation_.settings();
        if (!checkpoint.lastRelaxation.finite
            || checkpoint.lastRelaxation.completedIterationCount
                != iteration
            || checkpoint.lastRelaxation.relaxation
                != checkpoint.relaxation.relaxation
            || !finiteNonnegative(
                checkpoint.lastRelaxation.residualL2)
            || !finiteNonnegative(
                checkpoint.lastRelaxation.residualChangeL2)
            || !finiteNonnegative(
                checkpoint.lastRelaxation.denominator)
            || !std::isfinite(storedResidualL2)
            || checkpoint.lastRelaxation.residualL2 != storedResidualL2
            || checkpoint.lastRelaxation.residualChangeL2
                != storedChangeL2
            || (iteration == 1
                && (checkpoint.lastRelaxation.relaxation
                        != relaxationSettings.initialRelaxation
                    || checkpoint.lastRelaxation.denominator != 0.0
                    || checkpoint.lastRelaxation.usedDynamicRelaxation
                    || checkpoint.lastRelaxation.relaxationWasClipped))
            || (!checkpoint.lastRelaxation.usedDynamicRelaxation
                && checkpoint.lastRelaxation.relaxationWasClipped)
            || (checkpoint.lastRelaxation.usedDynamicRelaxation
                && !(checkpoint.lastRelaxation.denominator > 0.0))
            || (checkpoint.lastRelaxation.relaxationWasClipped
                && checkpoint.lastRelaxation.relaxation
                    != relaxationSettings.minimumRelaxation
                && checkpoint.lastRelaxation.relaxation
                    != relaxationSettings.maximumRelaxation)) {
            throw std::invalid_argument(
                "strong-coupling checkpoint relaxation diagnostics are invalid");
        }
        const CouplingConvergenceDecision expected =
            evaluateCouplingConvergence(
                iteration,
                checkpoint.lastConvergence.residuals,
                convergenceSettings_);
        const StrongCouplingIterationStatus expectedStatus =
            expected.converged
            ? StrongCouplingIterationStatus::Converged
            : (expected.iterationLimitReached
                ? StrongCouplingIterationStatus::Exhausted
                : StrongCouplingIterationStatus::Iterating);
        if (checkpoint.lastConvergence != expected
            || checkpoint.status != expectedStatus) {
            throw std::invalid_argument(
                "strong-coupling checkpoint convergence state is invalid");
        }
    }

    std::vector<double> restoredInterface = checkpoint.currentInterface;
    relaxation_.restore(checkpoint.relaxation);
    currentInterface_.swap(restoredInterface);
    status_ = checkpoint.status;
    lastRelaxation_ = checkpoint.lastRelaxation;
    lastConvergence_ = checkpoint.lastConvergence;
}

void StrongCouplingIteration::reset(
    const std::span<const double> initialInterface) {
    if (initialInterface.size() != currentInterface_.size()
        || !finite(initialInterface)) {
        throw std::invalid_argument(
            "strong-coupling reset interface is invalid");
    }
    std::vector<double> resetInterface(
        initialInterface.begin(), initialInterface.end());
    relaxation_.reset();
    currentInterface_.swap(resetInterface);
    status_ = StrongCouplingIterationStatus::Iterating;
    lastRelaxation_ = {};
    lastConvergence_ = {};
}

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

CouplingResidualNorms ConservativeMacroStepCoupling::measureResiduals(
    const std::span<const CouplingNodeKinematics> baselineKinematics,
    const std::span<const CouplingNodeKinematics> previousKinematics,
    const std::span<const CouplingNodeKinematics> currentKinematics,
    const ConservativeTransferResult& previousTraction,
    const ConservativeTransferResult& currentTraction) const {
    const auto previousLoads = previousTraction.nodeLoads();
    const auto currentLoads = currentTraction.nodeLoads();
    if (baselineKinematics.size() != nodes_.size()
        || previousKinematics.size() != nodes_.size()
        || currentKinematics.size() != nodes_.size()
        || previousLoads.size() != nodes_.size()
        || currentLoads.size() != nodes_.size()
        || previousTraction.surfaceFingerprint() != surfaceFingerprint_
        || currentTraction.surfaceFingerprint() != surfaceFingerprint_
        || previousTraction.targetDefinitionFingerprint()
            != targetDefinitionFingerprint_
        || currentTraction.targetDefinitionFingerprint()
            != targetDefinitionFingerprint_
        || !previousTraction.diagnostics().finite
        || !currentTraction.diagnostics().finite) {
        throw std::invalid_argument(
            "coupling residual inputs do not belong to this interface");
    }

    CouplingResidualNorms result;
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const auto& node = nodes_[index];
        const auto& baseline = baselineKinematics[index];
        const auto& previous = previousKinematics[index];
        const auto& current = currentKinematics[index];
        const auto& previousLoad = previousLoads[index];
        const auto& currentLoad = currentLoads[index];
        if (baseline.stableId != node.stableId
            || previous.stableId != node.stableId
            || current.stableId != node.stableId
            || previousLoad.stableId != node.stableId
            || currentLoad.stableId != node.stableId
            || previousLoad.structureNode != node.structureNode
            || currentLoad.structureNode != node.structureNode
            || !finite(baseline.positionMeters)
            || !finite(previous.positionMeters)
            || !finite(current.positionMeters)
            || !finite(baseline.velocityMetersPerSecond)
            || !finite(previous.velocityMetersPerSecond)
            || !finite(current.velocityMetersPerSecond)
            || !finite(previousLoad.forceNewtons)
            || !finite(currentLoad.forceNewtons)) {
            throw std::invalid_argument(
                "coupling residual inputs have invalid nodal bindings");
        }

        result.displacementMetres = std::max(
            result.displacementMetres,
            norm(subtract(
                current.positionMeters, previous.positionMeters)));
        result.displacementReferenceMetres = std::max({
            result.displacementReferenceMetres,
            norm(subtract(
                previous.positionMeters, baseline.positionMeters)),
            norm(subtract(
                current.positionMeters, baseline.positionMeters)),
        });
        result.velocityMetersPerSecond = std::max(
            result.velocityMetersPerSecond,
            norm(subtract(
                current.velocityMetersPerSecond,
                previous.velocityMetersPerSecond)));
        result.velocityReferenceMetersPerSecond = std::max({
            result.velocityReferenceMetersPerSecond,
            norm(subtract(
                previous.velocityMetersPerSecond,
                baseline.velocityMetersPerSecond)),
            norm(subtract(
                current.velocityMetersPerSecond,
                baseline.velocityMetersPerSecond)),
        });
        result.tractionNewtons = std::max(
            result.tractionNewtons,
            norm(subtract(
                currentLoad.forceNewtons, previousLoad.forceNewtons)));
        result.tractionReferenceNewtons = std::max({
            result.tractionReferenceNewtons,
            norm(previousLoad.forceNewtons),
            norm(currentLoad.forceNewtons),
        });
    }

    if (!finiteNonnegative(result.displacementMetres)
        || !finiteNonnegative(result.displacementReferenceMetres)
        || !finiteNonnegative(result.velocityMetersPerSecond)
        || !finiteNonnegative(
            result.velocityReferenceMetersPerSecond)
        || !finiteNonnegative(result.tractionNewtons)
        || !finiteNonnegative(result.tractionReferenceNewtons)) {
        throw std::overflow_error(
            "coupling residual reduction produced non-finite norms");
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
