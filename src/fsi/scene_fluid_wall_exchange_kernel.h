#pragma once

#include "scene_fluid_region_wall.h"

#include <cstddef>
#include <vector>

namespace simwing::fsi::detail {

// Source-neutral numerical core for two-sided material-wall exchange. Callers
// remain responsible for proving how collocated controls and quadrature sides
// map to their authoritative topology; this kernel owns only the bounded
// tangential impulse, energy, and action/reaction arithmetic.
struct SceneFluidWallExchangeKernelResult {
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionWallDiagnostics diagnostics;
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes;
    std::vector<SceneFluidRegionWallSample> samples;

    bool operator==(const SceneFluidWallExchangeKernelResult&) const = default;
};

[[nodiscard]] SceneFluidWallExchangeKernelResult
exchangeSceneFluidWallMomentumKernel(
    double densityKgPerCubicMeter,
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes,
    std::vector<SceneFluidRegionWallSample> samples,
    const SceneFluidRegionWallSettings& settings,
    std::size_t acceptedStorageBytes);

} // namespace simwing::fsi::detail
