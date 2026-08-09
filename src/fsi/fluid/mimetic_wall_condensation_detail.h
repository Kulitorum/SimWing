#pragma once

#include "mimetic_wall_condensation.h"

#include <span>
#include <vector>

namespace simwing::fsi::fluid::detail {

// Internal repeated-application path. The caller must have validated the
// condensation/operator binding and the finite, zero-on-wall field once for
// the enclosing transaction.
[[nodiscard]] std::vector<double>
applyMimeticWallCondensedTraceOperatorAssumingValidated(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    std::span<const double> activeTraceValues);

} // namespace simwing::fsi::fluid::detail
