#pragma once

#include "scene_fluid_mimetic_trace_system.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

struct SceneFluidMimeticTraceSolveSettings {
    double absoluteResidualTolerancePascalsMeters = 1.0e-12;
    double relativeResidualTolerance = 1.0e-10;
    double absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-12;
    std::size_t maximumIterations = 4000;

    bool operator==(
        const SceneFluidMimeticTraceSolveSettings&) const = default;
};

struct SceneFluidMimeticTraceSolveComponentDiagnostics {
    std::size_t componentIndex = 0;
    std::size_t traceCount = 0;
    std::size_t gaugeTraceIndex = 0;
    double rightHandSideSumPascalsMeters = 0.0;
    double compatibilityCorrectionPascalsMeters = 0.0;
    double traceGaugeBeforePascals = 0.0;
    double traceGaugeAfterPascals = 0.0;

    bool operator==(
        const SceneFluidMimeticTraceSolveComponentDiagnostics&) const =
        default;
};

struct SceneFluidMimeticTraceSolveDiagnostics {
    bool compatible = false;
    bool converged = false;
    bool finite = false;
    std::uint64_t traceSystemFingerprint = 0;
    std::size_t traceCount = 0;
    std::size_t componentCount = 0;
    std::size_t iterationCount = 0;
    double maximumAbsoluteComponentCompatibilityPascalsMeters = 0.0;
    double initialResidualL2PascalsMeters = 0.0;
    double finalResidualL2PascalsMeters = 0.0;
    double finalResidualMaximumPascalsMeters = 0.0;
    std::vector<SceneFluidMimeticTraceSolveComponentDiagnostics> components;

    bool operator==(
        const SceneFluidMimeticTraceSolveDiagnostics&) const = default;
};

// Solves H*lambda=b with the retained component gauge traces fixed exactly to
// zero. Only a component compatibility defect within the declared absolute
// tolerance is removed. Jacobi-preconditioned conjugate gradients uses the
// exact stored condensed diagonal and the matrix-free trace action. The caller
// field is committed only after a freshly recomputed full residual converges;
// incompatibility, non-finite arithmetic, or iteration exhaustion leaves its
// warm start bit-for-bit unchanged.
[[nodiscard]] SceneFluidMimeticTraceSolveDiagnostics
solveSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& tracePascals,
    const SceneFluidMimeticTraceSolveSettings& settings = {});

} // namespace simwing::fsi
