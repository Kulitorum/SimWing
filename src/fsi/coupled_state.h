#pragma once

#include "coupling.h"
#include "fluid/checkpoint.h"

#include <cstdint>
#include <functional>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t strongCouplingRollbackCheckpointVersion = 1;
inline constexpr std::uint32_t strongCouplingSolverCheckpointVersion = 1;
inline constexpr std::uint32_t couplingMacroStepRetryCheckpointVersion = 1;
inline constexpr std::uint32_t maximumCouplingMacroStepRetries = 64;

struct CouplingMacroStepRetrySettings {
    std::uint32_t maximumRetries = 4;
    double reductionFactor = 0.5;
    double minimumTimeStepSeconds = 1.0e-6;

    bool operator==(const CouplingMacroStepRetrySettings&) const = default;
};

enum class CouplingMacroStepRetryStatus : std::uint8_t {
    Attempting = 1,
    RetryPending = 2,
    Accepted = 3,
    Failed = 4,
};

struct CouplingMacroStepRetryDecision {
    CouplingMacroStepRetryStatus status =
        CouplingMacroStepRetryStatus::Attempting;
    std::uint32_t attemptNumber = 1;
    std::uint32_t retryCount = 0;
    double timeStepSeconds = 0.0;
    double pendingTimeStepSeconds = 0.0;

    bool operator==(const CouplingMacroStepRetryDecision&) const = default;
};

struct CouplingMacroStepRetryCheckpoint {
    std::uint32_t version = couplingMacroStepRetryCheckpointVersion;
    std::uint64_t macroStepDefinitionFingerprint = 0;
    CouplingMacroStepRetrySettings settings;
    double initialTimeStepSeconds = 0.0;
    double timeStepSeconds = 0.0;
    double pendingTimeStepSeconds = 0.0;
    std::uint32_t retryCount = 0;
    CouplingMacroStepRetryStatus status =
        CouplingMacroStepRetryStatus::Attempting;

    bool operator==(
        const CouplingMacroStepRetryCheckpoint&) const = default;
};

// Bounded macro-step retry policy. Exhaustion first creates a RetryPending
// state; the caller restores its composite baseline before acknowledging that
// retry. Only then does the smaller time step become the active next attempt.
// Convergence and unrecoverable exhaustion are distinct terminal states.
class CouplingMacroStepRetry final {
public:
    CouplingMacroStepRetry(
        std::uint64_t macroStepDefinitionFingerprint,
        double initialTimeStepSeconds,
        const CouplingMacroStepRetrySettings& settings = {});

    [[nodiscard]] CouplingMacroStepRetryStatus status() const noexcept;
    [[nodiscard]] std::uint32_t retryCount() const noexcept;
    [[nodiscard]] double timeStepSeconds() const noexcept;
    [[nodiscard]] double pendingTimeStepSeconds() const noexcept;
    [[nodiscard]] CouplingMacroStepRetryDecision decision() const noexcept;
    [[nodiscard]] CouplingMacroStepRetryCheckpoint checkpoint() const noexcept;

    [[nodiscard]] CouplingMacroStepRetryDecision reportIteration(
        StrongCouplingIterationStatus iterationStatus);
    [[nodiscard]] CouplingMacroStepRetryDecision beginRetry();
    void restore(const CouplingMacroStepRetryCheckpoint& checkpoint);

private:
    std::uint64_t macroStepDefinitionFingerprint_ = 0;
    CouplingMacroStepRetrySettings settings_;
    double initialTimeStepSeconds_ = 0.0;
    double timeStepSeconds_ = 0.0;
    double pendingTimeStepSeconds_ = 0.0;
    std::uint32_t retryCount_ = 0;
    CouplingMacroStepRetryStatus status_ =
        CouplingMacroStepRetryStatus::Attempting;
};

struct StrongCouplingRollbackCheckpoint {
    std::uint32_t version = strongCouplingRollbackCheckpointVersion;
    std::uint64_t interfaceDefinitionFingerprint = 0;
    StructureCheckpoint structure;
    fluid::MovingInterfaceFluidCheckpoint fluid;
    StrongCouplingIterationCheckpoint iteration;
};

struct StrongCouplingSolverCheckpoint {
    std::uint32_t version = strongCouplingSolverCheckpointVersion;
    std::uint64_t interfaceDefinitionFingerprint = 0;
    StructureCheckpoint structure;
    fluid::MovingInterfaceFluidCheckpoint fluid;
};

// Composite owner for the three states that must return to one macro-step
// baseline together: Structure, one accepted moving-interface fluid epoch,
// and the strong-iteration algorithm. Restore validates and allocates complete
// replacements first, then commits all three through no-throw moves.
class StrongCouplingRollbackState final {
public:
    StrongCouplingRollbackState(
        std::uint64_t interfaceDefinitionFingerprint,
        Structure structure,
        fluid::PeriodicCartesianGrid grid,
        fluid::MovingInterfaceFluidState fluidState,
        StrongCouplingIteration iteration);

    StrongCouplingRollbackState(const StrongCouplingRollbackState&) = delete;
    StrongCouplingRollbackState& operator=(
        const StrongCouplingRollbackState&) = delete;
    StrongCouplingRollbackState(
        StrongCouplingRollbackState&&) noexcept = default;
    StrongCouplingRollbackState& operator=(
        StrongCouplingRollbackState&&) noexcept = default;

    [[nodiscard]] std::uint64_t
    interfaceDefinitionFingerprint() const noexcept;
    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] Structure& structure() noexcept;
    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] fluid::MovingInterfaceFluidState& fluidState() noexcept;
    [[nodiscard]] const fluid::MovingInterfaceFluidState&
    fluidState() const noexcept;
    [[nodiscard]] StrongCouplingIteration& iteration() noexcept;
    [[nodiscard]] const StrongCouplingIteration& iteration() const noexcept;

    [[nodiscard]] StrongCouplingRollbackCheckpoint checkpoint() const;
    void restore(const StrongCouplingRollbackCheckpoint& checkpoint);
    [[nodiscard]] StrongCouplingSolverCheckpoint solverCheckpoint() const;
    void restoreSolvers(const StrongCouplingSolverCheckpoint& checkpoint);

private:
    std::uint64_t interfaceDefinitionFingerprint_ = 0;
    Structure structure_;
    fluid::PeriodicCartesianGrid grid_;
    fluid::MovingInterfaceFluidState fluidState_;
    StrongCouplingIteration iteration_;
};

// Owns the complete mutable state and immutable rollback baseline for one
// strong-coupling macro-step. A terminal iteration is reported to the bounded
// retry policy. Pending retries can only be activated through an operation
// that first restores Structure, fluid, and iteration to that baseline.
class StrongCouplingMacroStepState final {
public:
    StrongCouplingMacroStepState(
        StrongCouplingRollbackState state,
        double initialTimeStepSeconds,
        const CouplingMacroStepRetrySettings& retrySettings = {});

    StrongCouplingMacroStepState(
        const StrongCouplingMacroStepState&) = delete;
    StrongCouplingMacroStepState& operator=(
        const StrongCouplingMacroStepState&) = delete;
    StrongCouplingMacroStepState(
        StrongCouplingMacroStepState&&) noexcept = default;
    StrongCouplingMacroStepState& operator=(
        StrongCouplingMacroStepState&&) noexcept = default;

    [[nodiscard]] StrongCouplingRollbackState& rollbackState() noexcept;
    [[nodiscard]] const StrongCouplingRollbackState&
    rollbackState() const noexcept;
    [[nodiscard]] CouplingMacroStepRetryDecision decision() const noexcept;

    [[nodiscard]] CouplingMacroStepRetryDecision reportTerminalIteration();
    void restoreSolversForNextIteration();
    [[nodiscard]] CouplingMacroStepRetryDecision restoreAndBeginRetry();
    void restoreAttemptBaseline();

private:
    StrongCouplingRollbackState state_;
    StrongCouplingRollbackCheckpoint baseline_;
    StrongCouplingSolverCheckpoint solverBaseline_;
    CouplingMacroStepRetry retry_;
};

struct StrongCouplingSolverResult {
    std::vector<double> unrelaxedInterface;
    CouplingResidualNorms residuals;

    bool operator==(const StrongCouplingSolverResult&) const = default;
};

using StrongCouplingSolverCallback = std::function<StrongCouplingSolverResult(
    Structure&,
    const fluid::PeriodicCartesianGrid&,
    fluid::MovingInterfaceFluidState&,
    std::span<const double>,
    double)>;

struct StrongCouplingMacroStepAttempt {
    CouplingMacroStepRetryDecision decision;
    StrongCouplingIterationResult terminalIteration;
    std::uint64_t solverRunCount = 0;

    bool operator==(const StrongCouplingMacroStepAttempt&) const = default;
};

struct StrongCouplingMacroStepRunResult {
    CouplingMacroStepRetryDecision decision;
    StrongCouplingIterationResult lastIteration;
    std::uint64_t solverRunCount = 0;
    // One terminal record per attempted time step, bounded by the retry policy
    // to maximumCouplingMacroStepRetries + 1 entries.
    std::vector<StrongCouplingMacroStepAttempt> attempts;

    bool operator==(const StrongCouplingMacroStepRunResult&) const = default;
};

// Runs one complete bounded macro-step. The callback mutates only the physical
// solver owners and returns one unrelaxed interface candidate plus already
// reduced residuals. Nonterminal fixed-point iterations rewind Structure and
// fluid while retaining Aitken history. Exhausted attempts rewind all state
// before a smaller retry; terminal failure also leaves the accepted baseline
// intact. The returned attempt history retains each exhausted/accepted/failed
// terminal decision. Callback and validation failures restore the active
// baseline before propagating the exception.
[[nodiscard]] StrongCouplingMacroStepRunResult runStrongCouplingMacroStep(
    StrongCouplingMacroStepState& macroStep,
    const StrongCouplingSolverCallback& solve);

} // namespace simwing::fsi
