#pragma once

#include "coupling.h"
#include "fluid/checkpoint.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr std::uint32_t strongCouplingRollbackCheckpointVersion = 1;

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
