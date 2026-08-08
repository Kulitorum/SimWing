#include "coupled_state.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

fluid::MovingInterfaceFluidCheckpoint checkpointFluid(
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MovingInterfaceFluidState& state) {
    return fluid::checkpointMovingInterfaceFluidState(
        grid,
        state.velocityMetersPerSecond,
        state.pressurePascals,
        state.interfaces,
        state.diagnostics);
}

bool valid(const CouplingMacroStepRetrySettings& settings) {
    return settings.maximumRetries <= maximumCouplingMacroStepRetries
        && std::isfinite(settings.reductionFactor)
        && settings.reductionFactor > 0.0
        && settings.reductionFactor < 1.0
        && std::isfinite(settings.minimumTimeStepSeconds)
        && settings.minimumTimeStepSeconds > 0.0;
}

double reducedTimeStep(
    const double current,
    const CouplingMacroStepRetrySettings& settings) {
    return std::max(
        settings.minimumTimeStepSeconds,
        current * settings.reductionFactor);
}

} // namespace

static_assert(std::is_nothrow_move_assignable_v<Structure>);
static_assert(
    std::is_nothrow_move_assignable_v<fluid::MovingInterfaceFluidState>);
static_assert(std::is_nothrow_move_assignable_v<StrongCouplingIteration>);

CouplingMacroStepRetry::CouplingMacroStepRetry(
    const std::uint64_t macroStepDefinitionFingerprint,
    const double initialTimeStepSeconds,
    const CouplingMacroStepRetrySettings& settings)
    : macroStepDefinitionFingerprint_(macroStepDefinitionFingerprint),
      settings_(settings),
      initialTimeStepSeconds_(initialTimeStepSeconds),
      timeStepSeconds_(initialTimeStepSeconds) {
    if (macroStepDefinitionFingerprint_ == 0
        || !valid(settings_)
        || !std::isfinite(initialTimeStepSeconds_)
        || initialTimeStepSeconds_ < settings_.minimumTimeStepSeconds) {
        throw std::invalid_argument(
            "coupling macro-step retry identity, step, or settings are invalid");
    }
}

CouplingMacroStepRetryStatus
CouplingMacroStepRetry::status() const noexcept {
    return status_;
}

std::uint32_t CouplingMacroStepRetry::retryCount() const noexcept {
    return retryCount_;
}

double CouplingMacroStepRetry::timeStepSeconds() const noexcept {
    return timeStepSeconds_;
}

double CouplingMacroStepRetry::pendingTimeStepSeconds() const noexcept {
    return pendingTimeStepSeconds_;
}

CouplingMacroStepRetryDecision
CouplingMacroStepRetry::decision() const noexcept {
    return {
        status_,
        retryCount_ + 1,
        retryCount_,
        timeStepSeconds_,
        pendingTimeStepSeconds_,
    };
}

CouplingMacroStepRetryCheckpoint
CouplingMacroStepRetry::checkpoint() const noexcept {
    return {
        couplingMacroStepRetryCheckpointVersion,
        macroStepDefinitionFingerprint_,
        settings_,
        initialTimeStepSeconds_,
        timeStepSeconds_,
        pendingTimeStepSeconds_,
        retryCount_,
        status_,
    };
}

CouplingMacroStepRetryDecision CouplingMacroStepRetry::reportIteration(
    const StrongCouplingIterationStatus iterationStatus) {
    if (status_ != CouplingMacroStepRetryStatus::Attempting) {
        throw std::logic_error(
            "coupling macro-step attempt is not active");
    }
    if (iterationStatus == StrongCouplingIterationStatus::Converged) {
        status_ = CouplingMacroStepRetryStatus::Accepted;
        return decision();
    }
    if (iterationStatus != StrongCouplingIterationStatus::Exhausted) {
        throw std::invalid_argument(
            "coupling macro-step requires a terminal iteration result");
    }

    if (retryCount_ >= settings_.maximumRetries
        || timeStepSeconds_ <= settings_.minimumTimeStepSeconds) {
        status_ = CouplingMacroStepRetryStatus::Failed;
        return decision();
    }
    pendingTimeStepSeconds_ = reducedTimeStep(
        timeStepSeconds_, settings_);
    status_ = CouplingMacroStepRetryStatus::RetryPending;
    return decision();
}

CouplingMacroStepRetryDecision CouplingMacroStepRetry::beginRetry() {
    if (status_ != CouplingMacroStepRetryStatus::RetryPending) {
        throw std::logic_error(
            "coupling macro-step has no pending retry");
    }
    timeStepSeconds_ = pendingTimeStepSeconds_;
    pendingTimeStepSeconds_ = 0.0;
    ++retryCount_;
    status_ = CouplingMacroStepRetryStatus::Attempting;
    return decision();
}

void CouplingMacroStepRetry::restore(
    const CouplingMacroStepRetryCheckpoint& checkpointValue) {
    if (checkpointValue.version
            != couplingMacroStepRetryCheckpointVersion
        || checkpointValue.macroStepDefinitionFingerprint
            != macroStepDefinitionFingerprint_
        || checkpointValue.settings != settings_
        || !std::isfinite(checkpointValue.initialTimeStepSeconds)
        || checkpointValue.initialTimeStepSeconds
            < settings_.minimumTimeStepSeconds
        || checkpointValue.initialTimeStepSeconds
            != initialTimeStepSeconds_
        || !std::isfinite(checkpointValue.timeStepSeconds)
        || checkpointValue.timeStepSeconds
            < settings_.minimumTimeStepSeconds
        || !std::isfinite(checkpointValue.pendingTimeStepSeconds)
        || checkpointValue.pendingTimeStepSeconds < 0.0
        || checkpointValue.retryCount > settings_.maximumRetries) {
        throw std::invalid_argument(
            "coupling macro-step retry checkpoint is incompatible or invalid");
    }

    double expectedTimeStep = initialTimeStepSeconds_;
    for (std::uint32_t retry = 0;
         retry < checkpointValue.retryCount; ++retry) {
        expectedTimeStep = reducedTimeStep(expectedTimeStep, settings_);
    }
    const bool canRetry =
        checkpointValue.retryCount < settings_.maximumRetries
        && expectedTimeStep > settings_.minimumTimeStepSeconds;
    const double expectedPending = canRetry
        ? reducedTimeStep(expectedTimeStep, settings_)
        : 0.0;
    const bool validStatus =
        (checkpointValue.status == CouplingMacroStepRetryStatus::Attempting
            && checkpointValue.pendingTimeStepSeconds == 0.0)
        || (checkpointValue.status
                == CouplingMacroStepRetryStatus::RetryPending
            && canRetry
            && checkpointValue.pendingTimeStepSeconds == expectedPending)
        || (checkpointValue.status == CouplingMacroStepRetryStatus::Accepted
            && checkpointValue.pendingTimeStepSeconds == 0.0)
        || (checkpointValue.status == CouplingMacroStepRetryStatus::Failed
            && !canRetry
            && checkpointValue.pendingTimeStepSeconds == 0.0);
    if (checkpointValue.timeStepSeconds != expectedTimeStep
        || !validStatus) {
        throw std::invalid_argument(
            "coupling macro-step retry checkpoint state is inconsistent");
    }

    timeStepSeconds_ = checkpointValue.timeStepSeconds;
    pendingTimeStepSeconds_ = checkpointValue.pendingTimeStepSeconds;
    retryCount_ = checkpointValue.retryCount;
    status_ = checkpointValue.status;
}

StrongCouplingRollbackState::StrongCouplingRollbackState(
    const std::uint64_t interfaceDefinitionFingerprint,
    Structure structure,
    fluid::PeriodicCartesianGrid grid,
    fluid::MovingInterfaceFluidState fluidState,
    StrongCouplingIteration iteration)
    : interfaceDefinitionFingerprint_(interfaceDefinitionFingerprint),
      structure_(std::move(structure)),
      grid_(std::move(grid)),
      fluidState_(std::move(fluidState)),
      iteration_(std::move(iteration)) {
    if (interfaceDefinitionFingerprint_ == 0
        || iteration_.interfaceDefinitionFingerprint()
            != interfaceDefinitionFingerprint_) {
        throw std::invalid_argument(
            "strong-coupling rollback identity is invalid");
    }
    static_cast<void>(checkpointFluid(grid_, fluidState_));
}

std::uint64_t StrongCouplingRollbackState::
interfaceDefinitionFingerprint() const noexcept {
    return interfaceDefinitionFingerprint_;
}

const fluid::PeriodicCartesianGrid&
StrongCouplingRollbackState::grid() const noexcept {
    return grid_;
}

Structure& StrongCouplingRollbackState::structure() noexcept {
    return structure_;
}

const Structure& StrongCouplingRollbackState::structure() const noexcept {
    return structure_;
}

fluid::MovingInterfaceFluidState&
StrongCouplingRollbackState::fluidState() noexcept {
    return fluidState_;
}

const fluid::MovingInterfaceFluidState&
StrongCouplingRollbackState::fluidState() const noexcept {
    return fluidState_;
}

StrongCouplingIteration& StrongCouplingRollbackState::iteration() noexcept {
    return iteration_;
}

const StrongCouplingIteration&
StrongCouplingRollbackState::iteration() const noexcept {
    return iteration_;
}

StrongCouplingRollbackCheckpoint
StrongCouplingRollbackState::checkpoint() const {
    StrongCouplingRollbackCheckpoint result;
    result.interfaceDefinitionFingerprint =
        interfaceDefinitionFingerprint_;
    result.structure = structure_.checkpoint();
    result.fluid = checkpointFluid(grid_, fluidState_);
    result.iteration = iteration_.checkpoint();
    return result;
}

void StrongCouplingRollbackState::restore(
    const StrongCouplingRollbackCheckpoint& checkpointValue) {
    if (checkpointValue.version
            != strongCouplingRollbackCheckpointVersion
        || checkpointValue.interfaceDefinitionFingerprint
            != interfaceDefinitionFingerprint_
        || checkpointValue.iteration.relaxation
               .interfaceDefinitionFingerprint
            != interfaceDefinitionFingerprint_) {
        throw std::invalid_argument(
            "strong-coupling rollback checkpoint identity is invalid");
    }

    Structure restoredStructure(structure_.definition());
    restoredStructure.restore(checkpointValue.structure);
    fluid::MovingInterfaceFluidState restoredFluid =
        fluid::restoreMovingInterfaceFluidState(
            grid_, checkpointValue.fluid);

    StrongCouplingIteration restoredIteration(
        interfaceDefinitionFingerprint_,
        iteration_.currentInterface(),
        iteration_.relaxationSettings(),
        iteration_.convergenceSettings());
    restoredIteration.restore(checkpointValue.iteration);

    structure_ = std::move(restoredStructure);
    fluidState_ = std::move(restoredFluid);
    iteration_ = std::move(restoredIteration);
}

StrongCouplingSolverCheckpoint
StrongCouplingRollbackState::solverCheckpoint() const {
    StrongCouplingSolverCheckpoint result;
    result.interfaceDefinitionFingerprint =
        interfaceDefinitionFingerprint_;
    result.structure = structure_.checkpoint();
    result.fluid = checkpointFluid(grid_, fluidState_);
    return result;
}

void StrongCouplingRollbackState::restoreSolvers(
    const StrongCouplingSolverCheckpoint& checkpointValue) {
    if (checkpointValue.version != strongCouplingSolverCheckpointVersion
        || checkpointValue.interfaceDefinitionFingerprint
            != interfaceDefinitionFingerprint_) {
        throw std::invalid_argument(
            "strong-coupling solver checkpoint identity is invalid");
    }

    Structure restoredStructure(structure_.definition());
    restoredStructure.restore(checkpointValue.structure);
    fluid::MovingInterfaceFluidState restoredFluid =
        fluid::restoreMovingInterfaceFluidState(
            grid_, checkpointValue.fluid);

    structure_ = std::move(restoredStructure);
    fluidState_ = std::move(restoredFluid);
}

StrongCouplingMacroStepState::StrongCouplingMacroStepState(
    StrongCouplingRollbackState state,
    const double initialTimeStepSeconds,
    const CouplingMacroStepRetrySettings& retrySettings)
    : state_(std::move(state)),
      baseline_(state_.checkpoint()),
      solverBaseline_(state_.solverCheckpoint()),
      retry_(
          state_.interfaceDefinitionFingerprint(),
          initialTimeStepSeconds,
          retrySettings) {
    if (state_.iteration().status()
            != StrongCouplingIterationStatus::Iterating
        || state_.iteration().completedIterationCount() != 0) {
        throw std::invalid_argument(
            "strong-coupling macro-step requires a fresh iteration baseline");
    }
}

StrongCouplingRollbackState&
StrongCouplingMacroStepState::rollbackState() noexcept {
    return state_;
}

const StrongCouplingRollbackState&
StrongCouplingMacroStepState::rollbackState() const noexcept {
    return state_;
}

CouplingMacroStepRetryDecision
StrongCouplingMacroStepState::decision() const noexcept {
    return retry_.decision();
}

CouplingMacroStepRetryDecision
StrongCouplingMacroStepState::reportTerminalIteration() {
    return retry_.reportIteration(state_.iteration().status());
}

void StrongCouplingMacroStepState::restoreSolversForNextIteration() {
    if (retry_.status() != CouplingMacroStepRetryStatus::Attempting
        || state_.iteration().status()
            != StrongCouplingIterationStatus::Iterating
        || state_.iteration().completedIterationCount() == 0) {
        throw std::logic_error(
            "strong-coupling macro-step has no advanced next iteration");
    }
    state_.restoreSolvers(solverBaseline_);
}

CouplingMacroStepRetryDecision
StrongCouplingMacroStepState::restoreAndBeginRetry() {
    if (retry_.status() != CouplingMacroStepRetryStatus::RetryPending) {
        throw std::logic_error(
            "strong-coupling macro-step has no pending retry to restore");
    }
    state_.restore(baseline_);
    return retry_.beginRetry();
}

} // namespace simwing::fsi
