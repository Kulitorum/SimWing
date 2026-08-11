#include "scene_fluid_mimetic_pressure_solve.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double next = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) {
            correction_ += (sum_ - next) + value;
        } else {
            correction_ += (value - next) + sum_;
        }
        sum_ = next;
    }

    [[nodiscard]] double value() const noexcept {
        return sum_ + correction_;
    }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

bool finiteField(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

double rootMeanSquare(const std::span<const double> values) {
    CompensatedSum squared;
    for (const double value : values) squared.add(value * value);
    return std::sqrt(
        squared.value() / static_cast<double>(values.size()));
}

} // namespace

SceneFluidMimeticPressureSolveResult
solveSceneFluidMimeticPressureSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const std::span<const double> integratedCellSources,
    const std::span<const double> reducedTraceWarmStartPascals,
    const SceneFluidMimeticTraceSolveSettings& settings) {
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedSystem, fullSystem);
    if (integratedCellSources.size() != fullSystem.localOperators.size()
        || reducedTraceWarmStartPascals.size()
            != condensedSystem.traces.size()
        || !finiteField(integratedCellSources)
        || !finiteField(reducedTraceWarmStartPascals)) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure-solve fields are invalid");
    }

    SceneFluidMimeticPressureSolveResult result;
    result.fullTraceSystemFingerprint = fullSystem.fingerprint;
    result.condensedTraceSystemFingerprint = condensedSystem.fingerprint;
    const auto fullRightHandSide =
        buildSceneFluidMimeticTraceRightHandSide(
            fullSystem, integratedCellSources);
    const auto reducedRightHandSide =
        condenseSceneFluidMimeticTraceRightHandSide(
            condensedSystem, fullSystem, fullRightHandSide);
    std::vector<double> reducedCandidate(
        reducedTraceWarmStartPascals.begin(),
        reducedTraceWarmStartPascals.end());
    result.diagnostics.reducedTraceSolve =
        solveSceneFluidMimeticCondensedTraceSystem(
            condensedSystem, fullSystem, reducedRightHandSide,
            reducedCandidate, settings);
    if (!result.diagnostics.reducedTraceSolve.converged) return result;

    auto fullCandidate = reconstructSceneFluidMimeticFullTraces(
        condensedSystem, fullSystem, fullRightHandSide,
        reducedCandidate);
    auto evaluation = evaluateSceneFluidMimeticTraceSystem(
        fullSystem, fullCandidate, integratedCellSources);
    result.diagnostics.reconstructedFullResidualL2PascalsMeters =
        rootMeanSquare(evaluation.traceIntegratedOutwardFluxSums);
    result.diagnostics.reconstructedFullResidualMaximumPascalsMeters =
        evaluation.maximumTraceFluxImbalance;
    result.diagnostics.maximumCellConservationResidual =
        evaluation.maximumCellConservationResidual;
    result.diagnostics.reconstructedFullResidualTolerancePascalsMeters =
        std::max({
            settings.absoluteResidualTolerancePascalsMeters,
            settings.relativeResidualTolerance
                * result.diagnostics.reducedTraceSolve
                    .initialResidualL2PascalsMeters,
            settings
                .absoluteReconstructedResidualTolerancePascalsMeters,
            settings
                .absoluteComponentCompatibilityTolerancePascalsMeters});
    result.diagnostics.reconstructedFullResidualConverged =
        std::isfinite(result.diagnostics
            .reconstructedFullResidualL2PascalsMeters)
        && result.diagnostics.reconstructedFullResidualL2PascalsMeters
            <= result.diagnostics
                .reconstructedFullResidualTolerancePascalsMeters;
    if (!result.diagnostics.reconstructedFullResidualConverged) {
        return result;
    }

    result.diagnostics.accepted = true;
    result.reducedTracePascals = std::move(reducedCandidate);
    result.fullTracePascals = std::move(fullCandidate);
    result.evaluation = std::move(evaluation);
    return result;
}

SceneFluidMimeticPressureSolveResult
solveSceneFluidMimeticPressureSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticPressureSourceSet& sources,
    const std::span<const double> reducedTraceWarmStartPascals,
    const SceneFluidMimeticTraceSolveSettings& settings) {
    validateSceneFluidMimeticPressureSourceIntegrity(sources);
    if (sources.mimeticControlCellFingerprint
            != fullSystem.mimeticControlCellFingerprint
        || sources.structureDefinitionFingerprint
            != fullSystem.structureDefinitionFingerprint
        || sources.acceptedStepCount != fullSystem.acceptedStepCount
        || sources.simulationTimeSeconds != fullSystem.simulationTimeSeconds
        || sources.controls.size() != fullSystem.localOperators.size()) {
        throw std::invalid_argument(
            "scene fluid mimetic pressure source is foreign");
    }
    const auto integratedSources =
        sceneFluidMimeticIntegratedCellSources(sources);
    auto result = solveSceneFluidMimeticPressureSystem(
        condensedSystem, fullSystem, integratedSources,
        reducedTraceWarmStartPascals, settings);
    result.pressureSourceFingerprint = sources.fingerprint;
    return result;
}

} // namespace simwing::fsi
