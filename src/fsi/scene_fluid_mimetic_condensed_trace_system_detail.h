#pragma once

#include "scene_fluid_mimetic_condensed_trace_system.h"

#include <span>
#include <vector>

namespace simwing::fsi::detail {

// Internal repeated-application path. The caller must have validated the
// complete condensed/full binding and finite reduced field once for the
// enclosing transaction.
[[nodiscard]] std::vector<double>
applySceneFluidMimeticCondensedTraceOperatorAssumingValidated(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    std::span<const double> reducedTraceScalars);

} // namespace simwing::fsi::detail
