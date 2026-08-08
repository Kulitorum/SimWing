#include "coupled_state.h"

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

} // namespace

static_assert(std::is_nothrow_move_assignable_v<Structure>);
static_assert(
    std::is_nothrow_move_assignable_v<fluid::MovingInterfaceFluidState>);
static_assert(std::is_nothrow_move_assignable_v<StrongCouplingIteration>);

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

} // namespace simwing::fsi
