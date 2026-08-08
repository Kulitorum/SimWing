#pragma once

#include "coupling.h"
#include "fluid/checkpoint.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t strongCouplingRollbackCheckpointVersion = 1;
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

private:
    std::uint64_t interfaceDefinitionFingerprint_ = 0;
    Structure structure_;
    fluid::PeriodicCartesianGrid grid_;
    fluid::MovingInterfaceFluidState fluidState_;
    StrongCouplingIteration iteration_;
};

} // namespace simwing::fsi
