#include "fluid/planar_region_fragment_pressure_solve.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

void validateSettings(
    const PlanarPressureRegionFragmentPressureSolveSettings& settings) {
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
            < 0.0
        || settings.maximumIterations == 0) {
        throw std::invalid_argument(
            "planar regional pressure-correction solve settings are invalid");
    }
}

double vectorDot(const std::span<const double> first,
                 const std::span<const double> second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result += first[index] * second[index];
    }
    return result;
}

double vectorL2(const std::span<const double> values) {
    return std::sqrt(
        vectorDot(values, values) / static_cast<double>(values.size()));
}

double vectorMaximumAbsolute(const std::span<const double> values) {
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

bool applyOperator(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const std::span<const double> pressure,
    std::vector<double>& result) {
    result.assign(pressureOperator.rows.size(), 0.0);
    for (const auto& row : pressureOperator.rows) {
        for (std::size_t offset = 0; offset < row.entryCount; ++offset) {
            const auto& entry = pressureOperator.entries[
                row.firstEntry + offset];
            result[row.rowIndex] += entry.geometryWeightMeters
                * (pressure[row.rowIndex]
                   - pressure[entry.columnFragmentIndex]);
        }
        if (!std::isfinite(result[row.rowIndex])) return false;
    }
    return true;
}

void subtractComponentArithmeticMeans(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    std::vector<double>& values) {
    for (const auto& component : pressureOperator.components) {
        double sum = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            sum += values[
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset]];
        }
        const double mean = sum
            / static_cast<double>(component.fragmentCount);
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            values[pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset]] -= mean;
        }
    }
}

double componentVolumeMean(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentPressureOperatorComponent& component,
    const std::span<const double> values) {
    double moment = 0.0;
    for (std::size_t offset = 0;
         offset < component.fragmentCount; ++offset) {
        const std::size_t index =
            pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset];
        moment += fragments.fragments[index].volumeCubicMeters
            * values[index];
    }
    return moment / component.totalVolumeCubicMeters;
}

void subtractComponentVolumeMeans(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    std::vector<double>& values) {
    for (const auto& component : pressureOperator.components) {
        const double mean = componentVolumeMean(
            pressureOperator, fragments, component, values);
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            values[pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset]] -= mean;
        }
    }
}

bool allFinite(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

} // namespace

PlanarPressureRegionFragmentPressureSolveDiagnostics
solvePlanarPressureRegionFragmentPressureCorrection(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& correctionPascals,
    const PlanarPressureRegionFragmentPressureSolveSettings& settings) {
    validateSettings(settings);
    validatePlanarPressureRegionFragmentPressureOperator(
        pressureOperator, grid, sweep, fragments, topology);
    if (integratedRightHandSidePascalsMeters.size()
            != pressureOperator.rows.size()
        || correctionPascals.size() != pressureOperator.rows.size()
        || !allFinite(integratedRightHandSidePascalsMeters)
        || !allFinite(correctionPascals)) {
        throw std::invalid_argument(
            "planar regional pressure-correction fields are invalid");
    }

    PlanarPressureRegionFragmentPressureSolveDiagnostics diagnostics;
    diagnostics.finite = true;
    diagnostics.pressureOperatorFingerprint = pressureOperator.fingerprint;
    diagnostics.fragmentFingerprint = fragments.fingerprint;
    diagnostics.rowCount = pressureOperator.rows.size();
    diagnostics.componentCount = pressureOperator.components.size();
    diagnostics.components.reserve(pressureOperator.components.size());
    std::vector<double> rightHandSide(
        integratedRightHandSidePascalsMeters.begin(),
        integratedRightHandSidePascalsMeters.end());
    diagnostics.compatible = true;
    for (const auto& component : pressureOperator.components) {
        PlanarPressureRegionFragmentPressureSolveComponentDiagnostics
            componentDiagnostics;
        componentDiagnostics.componentIndex = component.componentIndex;
        componentDiagnostics.fragmentCount = component.fragmentCount;
        componentDiagnostics.totalVolumeCubicMeters =
            component.totalVolumeCubicMeters;
        componentDiagnostics.correctionVolumeMeanBeforePascals =
            componentVolumeMean(
                pressureOperator, fragments, component,
                correctionPascals);
        componentDiagnostics.correctionVolumeMeanAfterPascals =
            componentDiagnostics.correctionVolumeMeanBeforePascals;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            componentDiagnostics.rightHandSideSumPascalsMeters +=
                rightHandSide[
                    pressureOperator.componentFragmentIndices[
                        component.firstFragmentMember + offset]];
        }
        diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters =
            std::max(
                diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters,
                std::abs(componentDiagnostics
                    .rightHandSideSumPascalsMeters));
        if (!std::isfinite(
                componentDiagnostics.rightHandSideSumPascalsMeters)
            || !std::isfinite(
                componentDiagnostics.correctionVolumeMeanBeforePascals)) {
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
            / static_cast<double>(component.fragmentCount);
        diagnostics.components.push_back(componentDiagnostics);
    }
    if (!diagnostics.compatible) return diagnostics;

    for (const auto& componentDiagnostics : diagnostics.components) {
        const auto& component = pressureOperator.components[
            componentDiagnostics.componentIndex];
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            rightHandSide[
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset]]
                -= componentDiagnostics
                    .compatibilityCorrectionPascalsMeters;
        }
    }

    std::vector<double> candidate = correctionPascals;
    subtractComponentVolumeMeans(
        pressureOperator, fragments, candidate);
    std::vector<double> operatorCandidate;
    if (!applyOperator(
            pressureOperator, candidate, operatorCandidate)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    std::vector<double> residual(pressureOperator.rows.size(), 0.0);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        residual[row] = rightHandSide[row] - operatorCandidate[row];
    }
    subtractComponentArithmeticMeans(pressureOperator, residual);
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

    std::vector<double> direction = residual;
    std::vector<double> operatorDirection;
    double residualSquared = vectorDot(residual, residual);
    diagnostics.converged =
        diagnostics.finalResidualL2PascalsMeters <= convergenceThreshold;
    while (!diagnostics.converged
           && diagnostics.iterationCount < settings.maximumIterations) {
        if (!applyOperator(
                pressureOperator, direction, operatorDirection)) {
            diagnostics.finite = false;
            break;
        }
        const double denominator = vectorDot(
            direction, operatorDirection);
        if (!std::isfinite(denominator) || !(denominator > 0.0)
            || !std::isfinite(residualSquared)
            || !(residualSquared > 0.0)) {
            diagnostics.finite = false;
            break;
        }
        const double alpha = residualSquared / denominator;
        if (!std::isfinite(alpha)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t row = 0; row < residual.size(); ++row) {
            candidate[row] += alpha * direction[row];
            residual[row] -= alpha * operatorDirection[row];
        }
        ++diagnostics.iterationCount;
        subtractComponentVolumeMeans(
            pressureOperator, fragments, candidate);
        subtractComponentArithmeticMeans(pressureOperator, residual);
        const double nextResidualSquared = vectorDot(residual, residual);
        diagnostics.finalResidualL2PascalsMeters = std::sqrt(
            nextResidualSquared / static_cast<double>(residual.size()));
        diagnostics.finalResidualMaximumPascalsMeters =
            vectorMaximumAbsolute(residual);
        diagnostics.converged = std::isfinite(nextResidualSquared)
            && std::isfinite(diagnostics.finalResidualL2PascalsMeters)
            && diagnostics.finalResidualL2PascalsMeters
                <= convergenceThreshold;
        if (diagnostics.converged) {
            residualSquared = nextResidualSquared;
            break;
        }
        if (!std::isfinite(nextResidualSquared)
            || !(residualSquared > 0.0)) {
            diagnostics.finite = false;
            break;
        }
        const double beta = nextResidualSquared / residualSquared;
        if (!std::isfinite(beta)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t row = 0; row < residual.size(); ++row) {
            direction[row] = residual[row] + beta * direction[row];
        }
        subtractComponentArithmeticMeans(pressureOperator, direction);
        residualSquared = nextResidualSquared;
    }

    if (diagnostics.converged) {
        if (!applyOperator(
                pressureOperator, candidate, operatorCandidate)) {
            diagnostics.finite = false;
            diagnostics.converged = false;
        } else {
            for (std::size_t row = 0; row < residual.size(); ++row) {
                residual[row] =
                    rightHandSide[row] - operatorCandidate[row];
            }
            subtractComponentArithmeticMeans(
                pressureOperator, residual);
            diagnostics.finalResidualL2PascalsMeters = vectorL2(residual);
            diagnostics.finalResidualMaximumPascalsMeters =
                vectorMaximumAbsolute(residual);
            diagnostics.converged = std::isfinite(
                diagnostics.finalResidualL2PascalsMeters)
                && diagnostics.finalResidualL2PascalsMeters
                    <= convergenceThreshold;
        }
    }
    if (!diagnostics.converged) return diagnostics;

    subtractComponentVolumeMeans(
        pressureOperator, fragments, candidate);
    for (auto& componentDiagnostics : diagnostics.components) {
        componentDiagnostics.correctionVolumeMeanAfterPascals =
            componentVolumeMean(
                pressureOperator, fragments,
                pressureOperator.components[
                    componentDiagnostics.componentIndex],
                candidate);
        diagnostics.maximumAbsoluteCorrectionVolumeMeanPascals =
            std::max(
                diagnostics.maximumAbsoluteCorrectionVolumeMeanPascals,
                std::abs(componentDiagnostics
                    .correctionVolumeMeanAfterPascals));
    }
    if (!std::isfinite(
            diagnostics.maximumAbsoluteCorrectionVolumeMeanPascals)) {
        diagnostics.finite = false;
        diagnostics.converged = false;
        return diagnostics;
    }
    correctionPascals = std::move(candidate);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
