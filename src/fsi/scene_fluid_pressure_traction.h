#pragma once

#include "scene_fluid_quadrature.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

struct SceneFluidQuadraturePressure {
    std::uint64_t stableId = 0;
    double negativeSidePressurePascals = 0.0;
    double positiveSidePressurePascals = 0.0;

    bool operator==(const SceneFluidQuadraturePressure&) const = default;
};

// Converts explicit one-sided CFD pressure samples into net sheet traction:
// (p_negative - p_positive) times the current oriented triangle normal. It
// does not invent pressure, shear, a polar force, or a second load path.
[[nodiscard]] std::vector<SceneFluidQuadratureTraction>
buildSceneFluidPressureTractions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& quadrature,
    std::span<const SceneFluidQuadraturePressure> pressures);

[[nodiscard]] ConservativeTransferResult evaluateSceneFluidPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    std::span<const SceneFluidQuadraturePressure> pressures,
    const ConservativeTransferSettings& settings = {});

} // namespace simwing::fsi
