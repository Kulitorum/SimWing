#pragma once

#include "scene_fluid_regional_opening_momentum_wall_input.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidRegionalOpeningMomentumWallExchangeVersion = 1;

struct SceneFluidRegionalOpeningMomentumWallExchangeLimits {
    std::size_t maximumOwnedBytes = 16384ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 16384ULL * 1024ULL * 1024ULL;
};

// Immutable opt-in execution of the shared two-sided tangential wall law on
// one exact same-epoch regional opening input. The complete zero-exchange
// input is retained beside the adjusted collocated controls, per-sample fluid
// impulses, equal-and-opposite Structure tractions, and energy/conservation
// diagnostics so integrity can rerun the numerical kernel independently.
//
// Accepted adjusted controls are not fed into the next pressure prediction and
// tractions are not applied to Structure here. No cycle owner, topology rebase,
// or production worker is mutated or selected.
struct SceneFluidRegionalOpeningMomentumWallExchange {
    std::uint32_t version =
        sceneFluidRegionalOpeningMomentumWallExchangeVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceWallInputFingerprint = 0;
    SceneFluidRegionWallSettings settings;
    SceneFluidRegionalOpeningMomentumWallInput sourceInput;
    SceneFluidRegionWallDiagnostics diagnostics;
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes;
    std::vector<SceneFluidRegionWallSample> samples;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const SceneFluidRegionalOpeningMomentumWallExchange&) const = default;
};

[[nodiscard]] SceneFluidRegionalOpeningMomentumWallExchange
exchangeSceneFluidRegionalOpeningMomentumWall(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const SceneFluidRegionWallSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallExchangeLimits& limits = {});

void validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange);

void validateSceneFluidRegionalOpeningMomentumWallExchange(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const SceneFluidRegionWallSettings& settings = {},
    const SceneFluidRegionalOpeningMomentumWallExchangeLimits& limits = {});

} // namespace simwing::fsi
