#pragma once

#include "scene_fluid_mimetic_pressure_source.h"
#include "scene_fluid_mimetic_trace_solve.h"

#include <span>
#include <vector>

namespace simwing::fsi {

struct SceneFluidMimeticPressureSolveDiagnostics {
    bool accepted = false;
    bool reconstructedFullResidualConverged = false;
    double reconstructedFullResidualTolerancePascalsMeters = 0.0;
    double reconstructedFullResidualL2PascalsMeters = 0.0;
    double reconstructedFullResidualMaximumPascalsMeters = 0.0;
    double maximumCellConservationResidual = 0.0;
    SceneFluidMimeticTraceSolveDiagnostics reducedTraceSolve;

    bool operator==(
        const SceneFluidMimeticPressureSolveDiagnostics&) const = default;
};

// One atomic source-to-pressure audit result. State fields remain empty unless
// the reduced solve converges and the reconstructed full trace equations close
// inside the larger of the declared RMS and admitted component-compatibility
// tolerances. The caller-supplied reduced field is a read-only warm start; no
// partial candidate can escape on failure.
struct SceneFluidMimeticPressureSolveResult {
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t pressureSourceFingerprint = 0;
    SceneFluidMimeticPressureSolveDiagnostics diagnostics;
    std::vector<double> reducedTracePascals;
    std::vector<double> fullTracePascals;
    SceneFluidMimeticTraceEvaluation evaluation;
};

[[nodiscard]] SceneFluidMimeticPressureSolveResult
solveSceneFluidMimeticPressureSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    std::span<const double> integratedCellSources,
    std::span<const double> reducedTraceWarmStartPascals,
    const SceneFluidMimeticTraceSolveSettings& settings = {});

[[nodiscard]] SceneFluidMimeticPressureSolveResult
solveSceneFluidMimeticPressureSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    std::span<const double> reducedTraceWarmStartPascals,
    const SceneFluidMimeticTraceSolveSettings& settings = {});

} // namespace simwing::fsi
