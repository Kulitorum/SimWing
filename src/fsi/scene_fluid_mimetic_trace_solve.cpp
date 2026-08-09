#include "scene_fluid_mimetic_trace_solve.h"
#include "scene_fluid_mimetic_condensed_trace_system_detail.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

void validateSettings(
    const SceneFluidMimeticTraceSolveSettings& settings) {
    if (!std::isfinite(
            settings.absoluteResidualTolerancePascalsMeters)
        || settings.absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(settings.relativeResidualTolerance)
        || settings.relativeResidualTolerance < 0.0
        || (settings.absoluteResidualTolerancePascalsMeters == 0.0
            && settings.relativeResidualTolerance == 0.0)
        || !std::isfinite(
            settings.absoluteComponentCompatibilityTolerancePascalsMeters)
        || settings.absoluteComponentCompatibilityTolerancePascalsMeters
            < 0.0) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-solve settings are invalid");
    }
}

double vectorDot(const std::span<const double> first,
                 const std::span<const double> second) {
    CompensatedSum sum;
    for (std::size_t index = 0; index < first.size(); ++index) {
        sum.add(first[index] * second[index]);
    }
    return sum.value();
}

double vectorL2(const std::span<const double> values) {
    const double squared = vectorDot(values, values);
    return std::sqrt(squared / static_cast<double>(values.size()));
}

double vectorMaximumAbsolute(const std::span<const double> values) {
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

template<typename System>
void shiftComponentGauges(
    const System& system,
    std::vector<double>& traceValues) {
    std::vector<double> gaugeValues(system.componentCount, 0.0);
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        gaugeValues[component] = traceValues[
            system.componentGaugeTraceIndices[component]];
    }
    for (const auto& trace : system.traces) {
        traceValues[trace.traceIndex] -=
            gaugeValues[trace.componentIndex];
    }
    for (const std::size_t gauge : system.componentGaugeTraceIndices) {
        traceValues[gauge] = 0.0;
    }
}

bool allFinite(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

template<typename ApplyOperator>
bool applyOperatorFinite(
    ApplyOperator&& applyOperator,
    const std::span<const double> values,
    std::vector<double>& result) {
    try {
        result = applyOperator(values);
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::overflow_error&) {
        return false;
    }
    return allFinite(result);
}

template<typename System, typename ApplyOperator>
SceneFluidMimeticTraceSolveDiagnostics solveTraceSystem(
    const System& system,
    const std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& tracePascals,
    const SceneFluidMimeticTraceSolveSettings& settings,
    ApplyOperator&& applyOperator) {
    validateSettings(settings);
    if (integratedRightHandSidePascalsMeters.size()
            != system.traces.size()
        || tracePascals.size() != system.traces.size()
        || !allFinite(integratedRightHandSidePascalsMeters)
        || !allFinite(tracePascals)) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-solve fields are invalid");
    }

    SceneFluidMimeticTraceSolveDiagnostics diagnostics;
    diagnostics.finite = true;
    diagnostics.traceSystemFingerprint = system.fingerprint;
    diagnostics.traceCount = system.traces.size();
    diagnostics.componentCount = system.componentCount;
    diagnostics.components.resize(system.componentCount);
    std::vector<CompensatedSum> componentSums(system.componentCount);
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        auto& componentDiagnostics = diagnostics.components[component];
        componentDiagnostics.componentIndex = component;
        componentDiagnostics.gaugeTraceIndex =
            system.componentGaugeTraceIndices[component];
        componentDiagnostics.traceGaugeBeforePascals = tracePascals[
            componentDiagnostics.gaugeTraceIndex];
        componentDiagnostics.traceGaugeAfterPascals =
            componentDiagnostics.traceGaugeBeforePascals;
    }
    for (const auto& trace : system.traces) {
        auto& component = diagnostics.components[trace.componentIndex];
        ++component.traceCount;
        componentSums[trace.componentIndex].add(
            integratedRightHandSidePascalsMeters[trace.traceIndex]);
    }

    diagnostics.compatible = true;
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        auto& componentDiagnostics = diagnostics.components[component];
        componentDiagnostics.rightHandSideSumPascalsMeters =
            componentSums[component].value();
        diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters =
            std::max(
                diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters,
                std::abs(componentDiagnostics
                    .rightHandSideSumPascalsMeters));
        if (!std::isfinite(
                componentDiagnostics.rightHandSideSumPascalsMeters)) {
            diagnostics.finite = false;
            diagnostics.compatible = false;
        } else if (std::abs(
                       componentDiagnostics.rightHandSideSumPascalsMeters)
                   > settings
                       .absoluteComponentCompatibilityTolerancePascalsMeters) {
            diagnostics.compatible = false;
        }
        componentDiagnostics.compatibilityCorrectionPascalsMeters =
            componentDiagnostics.rightHandSideSumPascalsMeters
            / static_cast<double>(componentDiagnostics.traceCount);
    }
    if (!diagnostics.compatible) return diagnostics;

    std::vector<double> rightHandSide(
        integratedRightHandSidePascalsMeters.begin(),
        integratedRightHandSidePascalsMeters.end());
    for (const auto& trace : system.traces) {
        rightHandSide[trace.traceIndex] -=
            diagnostics.components[trace.componentIndex]
                .compatibilityCorrectionPascalsMeters;
    }
    std::vector<CompensatedSum> correctedFreeSums(system.componentCount);
    for (const auto& trace : system.traces) {
        if (!trace.isGauge) {
            correctedFreeSums[trace.componentIndex].add(
                rightHandSide[trace.traceIndex]);
        }
    }
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        rightHandSide[system.componentGaugeTraceIndices[component]] =
            -correctedFreeSums[component].value();
    }

    std::vector<double> candidateTrace = tracePascals;
    shiftComponentGauges(system, candidateTrace);
    std::vector<double> operatorTrace;
    if (!applyOperatorFinite(
            applyOperator, candidateTrace, operatorTrace)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    std::vector<double> residual(system.traces.size(), 0.0);
    for (std::size_t trace = 0; trace < residual.size(); ++trace) {
        residual[trace] = rightHandSide[trace] - operatorTrace[trace];
    }
    if (!allFinite(residual)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    diagnostics.initialResidualL2PascalsMeters = vectorL2(residual);
    diagnostics.finalResidualL2PascalsMeters =
        diagnostics.initialResidualL2PascalsMeters;
    diagnostics.finalResidualMaximumPascalsMeters =
        vectorMaximumAbsolute(residual);
    if (!std::isfinite(diagnostics.initialResidualL2PascalsMeters)
        || !std::isfinite(
            diagnostics.finalResidualMaximumPascalsMeters)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    const double convergenceThreshold = std::max(
        settings.absoluteResidualTolerancePascalsMeters,
        settings.relativeResidualTolerance
            * diagnostics.initialResidualL2PascalsMeters);
    diagnostics.converged =
        diagnostics.finalResidualL2PascalsMeters <= convergenceThreshold;

    std::vector<double> preconditionedResidual(system.traces.size(), 0.0);
    std::vector<double> direction(system.traces.size(), 0.0);
    std::vector<double> operatorDirection;
    double residualPreconditionedResidual = 0.0;
    if (!diagnostics.converged) {
        for (const auto& trace : system.traces) {
            if (!trace.isGauge) {
                preconditionedResidual[trace.traceIndex] =
                    residual[trace.traceIndex] / trace.operatorDiagonal;
                direction[trace.traceIndex] =
                    preconditionedResidual[trace.traceIndex];
            }
        }
        residualPreconditionedResidual = vectorDot(
            residual, preconditionedResidual);
        if (!std::isfinite(residualPreconditionedResidual)) {
            diagnostics.finite = false;
            return diagnostics;
        }
    }

    while (!diagnostics.converged
           && diagnostics.iterationCount < settings.maximumIterations) {
        if (!(residualPreconditionedResidual > 0.0)) break;
        if (!applyOperatorFinite(
                applyOperator, direction, operatorDirection)) {
            diagnostics.finite = false;
            break;
        }
        const double denominator = vectorDot(
            direction, operatorDirection);
        if (!std::isfinite(denominator)) {
            diagnostics.finite = false;
            break;
        }
        if (!(denominator > 0.0)) break;
        const double alpha = residualPreconditionedResidual / denominator;
        if (!std::isfinite(alpha)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t trace = 0; trace < residual.size(); ++trace) {
            candidateTrace[trace] += alpha * direction[trace];
            residual[trace] -= alpha * operatorDirection[trace];
        }
        ++diagnostics.iterationCount;
        for (const std::size_t gauge :
             system.componentGaugeTraceIndices) {
            candidateTrace[gauge] = 0.0;
        }
        if (!allFinite(candidateTrace) || !allFinite(residual)) {
            diagnostics.finite = false;
            break;
        }
        diagnostics.finalResidualL2PascalsMeters = vectorL2(residual);
        diagnostics.finalResidualMaximumPascalsMeters =
            vectorMaximumAbsolute(residual);
        diagnostics.converged = std::isfinite(
            diagnostics.finalResidualL2PascalsMeters)
            && diagnostics.finalResidualL2PascalsMeters
                <= convergenceThreshold;
        if (diagnostics.converged) break;

        for (const auto& trace : system.traces) {
            preconditionedResidual[trace.traceIndex] = trace.isGauge
                ? 0.0
                : residual[trace.traceIndex] / trace.operatorDiagonal;
        }
        const double nextResidualPreconditionedResidual = vectorDot(
            residual, preconditionedResidual);
        if (!std::isfinite(nextResidualPreconditionedResidual)) {
            diagnostics.finite = false;
            break;
        }
        if (!(nextResidualPreconditionedResidual > 0.0)) break;
        const double beta = nextResidualPreconditionedResidual
            / residualPreconditionedResidual;
        if (!std::isfinite(beta)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t trace = 0; trace < direction.size(); ++trace) {
            direction[trace] = preconditionedResidual[trace]
                + beta * direction[trace];
        }
        for (const std::size_t gauge :
             system.componentGaugeTraceIndices) {
            direction[gauge] = 0.0;
        }
        residualPreconditionedResidual =
            nextResidualPreconditionedResidual;
    }

    if (diagnostics.converged) {
        if (!applyOperatorFinite(
                applyOperator, candidateTrace, operatorTrace)) {
            diagnostics.finite = false;
            diagnostics.converged = false;
            return diagnostics;
        }
        for (std::size_t trace = 0; trace < residual.size(); ++trace) {
            residual[trace] = rightHandSide[trace] - operatorTrace[trace];
        }
        diagnostics.finalResidualL2PascalsMeters = vectorL2(residual);
        diagnostics.finalResidualMaximumPascalsMeters =
            vectorMaximumAbsolute(residual);
        diagnostics.converged = allFinite(residual)
            && std::isfinite(diagnostics.finalResidualL2PascalsMeters)
            && diagnostics.finalResidualL2PascalsMeters
                <= convergenceThreshold;
    }
    if (!diagnostics.converged) return diagnostics;

    shiftComponentGauges(system, candidateTrace);
    for (auto& component : diagnostics.components) {
        component.traceGaugeAfterPascals = candidateTrace[
            component.gaugeTraceIndex];
    }
    tracePascals = std::move(candidateTrace);
    return diagnostics;
}

} // namespace

SceneFluidMimeticTraceSolveDiagnostics
solveSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    const std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& tracePascals,
    const SceneFluidMimeticTraceSolveSettings& settings) {
    validateSceneFluidMimeticTraceSystemIntegrity(system);
    const auto applyOperator = [&](const std::span<const double> values) {
        return applySceneFluidMimeticTraceOperator(system, values);
    };
    return solveTraceSystem(
        system, integratedRightHandSidePascalsMeters, tracePascals,
        settings, applyOperator);
}

SceneFluidMimeticTraceSolveDiagnostics
solveSceneFluidMimeticCondensedTraceSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const std::span<const double>
        integratedReducedRightHandSidePascalsMeters,
    std::vector<double>& reducedTracePascals,
    const SceneFluidMimeticTraceSolveSettings& settings) {
    validateSceneFluidMimeticCondensedTraceSystem(
        condensedSystem, fullSystem);
    const auto applyOperator = [&](const std::span<const double> values) {
        return detail::
            applySceneFluidMimeticCondensedTraceOperatorAssumingValidated(
            condensedSystem, fullSystem, values);
    };
    return solveTraceSystem(
        condensedSystem,
        integratedReducedRightHandSidePascalsMeters,
        reducedTracePascals, settings, applyOperator);
}

} // namespace simwing::fsi
