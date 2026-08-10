#include "scene_fluid_wall_exchange_kernel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>

namespace simwing::fsi::detail {
namespace {

fluid::Vector3 add(const fluid::Vector3& first,
                   const fluid::Vector3& second) {
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

fluid::Vector3 subtract(const fluid::Vector3& first,
                        const fluid::Vector3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

fluid::Vector3 scale(const fluid::Vector3& value, const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(const fluid::Vector3& first, const fluid::Vector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

double norm(const fluid::Vector3& value) {
    return std::sqrt(dot(value, value));
}

bool finite(const fluid::Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double tolerance(const double absolute,
                 const double relative,
                 const double reference) {
    return absolute + relative * std::max(1.0, std::abs(reference));
}

double kineticEnergy(
    const std::span<const SceneFluidRegionWallControlVolume> controls,
    const double density) {
    double result = 0.0;
    for (const auto& control : controls) {
        result += 0.5 * density * control.volumeCubicMeters
            * dot(control.velocityMetersPerSecond,
                  control.velocityMetersPerSecond);
    }
    return result;
}

fluid::Vector3 totalMomentum(
    const std::span<const SceneFluidRegionWallControlVolume> controls) {
    fluid::Vector3 result;
    for (const auto& control : controls) {
        result = add(result, control.momentumKilogramMetersPerSecond);
    }
    return result;
}

bool validSettings(const SceneFluidRegionWallSettings& settings) {
    return std::isfinite(settings.timeStepSeconds)
        && settings.timeStepSeconds > 0.0
        && std::isfinite(settings.kinematicViscositySquareMetersPerSecond)
        && settings.kinematicViscositySquareMetersPerSecond >= 0.0
        && std::isfinite(settings.minimumWallDistanceMeters)
        && settings.minimumWallDistanceMeters > 0.0
        && std::isfinite(settings.maximumViscousNumber)
        && settings.maximumViscousNumber > 0.0
        && settings.maximumViscousNumber <= 1.0
        && settings.maximumSubsteps != 0
        && std::isfinite(
            settings.absoluteMomentumToleranceKilogramMetersPerSecond)
        && settings.absoluteMomentumToleranceKilogramMetersPerSecond >= 0.0
        && std::isfinite(settings.relativeMomentumTolerance)
        && settings.relativeMomentumTolerance >= 0.0
        && std::isfinite(settings.absoluteEnergyToleranceJoules)
        && settings.absoluteEnergyToleranceJoules >= 0.0
        && std::isfinite(settings.relativeEnergyTolerance)
        && settings.relativeEnergyTolerance >= 0.0;
}

void validateInput(
    const double density,
    const std::span<const SceneFluidRegionWallControlVolume> controls,
    const std::span<const SceneFluidRegionWallSample> samples,
    const SceneFluidRegionWallSettings& settings,
    const std::size_t acceptedStorageBytes) {
    if (!(density > 0.0) || !std::isfinite(density)
        || controls.empty() || samples.empty() || !validSettings(settings)
        || acceptedStorageBytes == 0) {
        throw std::invalid_argument(
            "scene fluid wall-exchange kernel input is invalid");
    }
    for (std::size_t index = 0; index < controls.size(); ++index) {
        const auto& control = controls[index];
        if (control.controlVolumeIndex != index || control.stableId == 0
            || !(control.volumeCubicMeters > 0.0)
            || !std::isfinite(control.volumeCubicMeters)
            || !(control.incidentWallAreaSquareMeters >= 0.0)
            || !std::isfinite(control.incidentWallAreaSquareMeters)
            || !finite(control.velocityMetersPerSecond)
            || !finite(control.momentumKilogramMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid wall-exchange kernel control is invalid");
        }
    }
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& sample = samples[index];
        if (sample.sampleIndex != index || sample.stableId == 0
            || sample.triangleId == invalidStableId
            || sample.negativeSideControlVolumeIndex >= controls.size()
            || sample.positiveSideControlVolumeIndex >= controls.size()
            || !(sample.areaSquareMeters > 0.0)
            || !std::isfinite(sample.areaSquareMeters)
            || !finite(sample.unitNormalNegativeToPositive)
            || std::abs(norm(sample.unitNormalNegativeToPositive) - 1.0)
                > 1.0e-12
            || !finite(sample.wallVelocityMetersPerSecond)
            || sample.structureTraction.stableId != sample.stableId) {
            throw std::invalid_argument(
                "scene fluid wall-exchange kernel sample is invalid");
        }
    }
}

} // namespace

SceneFluidWallExchangeKernelResult exchangeSceneFluidWallMomentumKernel(
    const double densityKgPerCubicMeter,
    std::vector<SceneFluidRegionWallControlVolume> controlVolumes,
    std::vector<SceneFluidRegionWallSample> samples,
    const SceneFluidRegionWallSettings& settings,
    const std::size_t acceptedStorageBytes) {
    validateInput(
        densityKgPerCubicMeter, controlVolumes, samples, settings,
        acceptedStorageBytes);

    SceneFluidWallExchangeKernelResult result;
    auto& diagnostics = result.diagnostics;
    diagnostics.controlVolumeCount = controlVolumes.size();
    diagnostics.quadraturePointCount = samples.size();

    std::vector<double> distances(controlVolumes.size(), 0.0);
    std::vector<double> viscousRows(controlVolumes.size(), 0.0);
    for (std::size_t index = 0; index < controlVolumes.size(); ++index) {
        if (!(controlVolumes[index].incidentWallAreaSquareMeters > 0.0)) {
            continue;
        }
        distances[index] = std::max(
            settings.minimumWallDistanceMeters,
            0.5 * controlVolumes[index].volumeCubicMeters
                / controlVolumes[index].incidentWallAreaSquareMeters);
        diagnostics.maximumWallDistanceMeters = std::max(
            diagnostics.maximumWallDistanceMeters, distances[index]);
    }
    for (const auto& sample : samples) {
        for (const std::size_t controlIndex : {
                 sample.negativeSideControlVolumeIndex,
                 sample.positiveSideControlVolumeIndex}) {
            viscousRows[controlIndex] +=
                settings.kinematicViscositySquareMetersPerSecond
                * sample.areaSquareMeters
                / (controlVolumes[controlIndex].volumeCubicMeters
                   * distances[controlIndex]);
        }
    }
    for (const double row : viscousRows) {
        diagnostics.maximumFullStepViscousNumber = std::max(
            diagnostics.maximumFullStepViscousNumber,
            settings.timeStepSeconds * row);
    }
    const double requiredSubsteps = std::ceil(
        diagnostics.maximumFullStepViscousNumber
        / settings.maximumViscousNumber);
    if (!std::isfinite(requiredSubsteps)) {
        throw std::overflow_error(
            "scene fluid region wall substep count is non-finite");
    }
    const double maximumRepresentableSubsteps = static_cast<double>(
        std::numeric_limits<std::size_t>::max());
    if (requiredSubsteps >= maximumRepresentableSubsteps
        || requiredSubsteps
            > static_cast<double>(settings.maximumSubsteps)) {
        diagnostics.substepCount = settings.maximumSubsteps;
        diagnostics.failureStage =
            SceneFluidRegionWallFailureStage::SubstepLimit;
    } else {
        diagnostics.substepCount = std::max<std::size_t>(
            1, static_cast<std::size_t>(requiredSubsteps));
        if (diagnostics.substepCount > settings.maximumSubsteps) {
            diagnostics.substepCount = settings.maximumSubsteps;
            diagnostics.failureStage =
                SceneFluidRegionWallFailureStage::SubstepLimit;
        }
    }
    diagnostics.maximumAcceptedSubstepViscousNumber =
        diagnostics.maximumFullStepViscousNumber
        / static_cast<double>(diagnostics.substepCount);
    diagnostics.fluidMomentumBeforeKilogramMetersPerSecond =
        totalMomentum(controlVolumes);
    diagnostics.kineticEnergyBeforeJoules = kineticEnergy(
        controlVolumes, densityKgPerCubicMeter);
    if (diagnostics.failureStage != SceneFluidRegionWallFailureStage::None) {
        diagnostics.fluidMomentumAfterKilogramMetersPerSecond =
            diagnostics.fluidMomentumBeforeKilogramMetersPerSecond;
        diagnostics.kineticEnergyAfterJoules =
            diagnostics.kineticEnergyBeforeJoules;
        diagnostics.finite = true;
        return result;
    }

    const double substepSeconds = settings.timeStepSeconds
        / static_cast<double>(diagnostics.substepCount);
    std::vector<fluid::Vector3> impulses(controlVolumes.size());
    for (std::size_t substep = 0;
         substep < diagnostics.substepCount; ++substep) {
        std::ranges::fill(impulses, fluid::Vector3{});
        double wallWork = 0.0;
        const double energyBefore = kineticEnergy(
            controlVolumes, densityKgPerCubicMeter);
        for (auto& sample : samples) {
            const auto exchangeSide = [&](
                const std::size_t controlIndex,
                fluid::Vector3& accumulatedSampleImpulse) {
                const auto relative = subtract(
                    controlVolumes[controlIndex].velocityMetersPerSecond,
                    sample.wallVelocityMetersPerSecond);
                const auto tangential = subtract(
                    relative,
                    scale(sample.unitNormalNegativeToPositive,
                          dot(relative,
                              sample.unitNormalNegativeToPositive)));
                diagnostics.maximumRelativeTangentialSpeedMetersPerSecond =
                    std::max(
                        diagnostics
                            .maximumRelativeTangentialSpeedMetersPerSecond,
                        norm(tangential));
                const double coefficient = densityKgPerCubicMeter
                    * settings.kinematicViscositySquareMetersPerSecond
                    * sample.areaSquareMeters / distances[controlIndex];
                const auto impulse = scale(
                    tangential, -coefficient * substepSeconds);
                impulses[controlIndex] = add(
                    impulses[controlIndex], impulse);
                accumulatedSampleImpulse = add(
                    accumulatedSampleImpulse, impulse);
                wallWork += dot(
                    impulse, sample.wallVelocityMetersPerSecond);
            };
            exchangeSide(
                sample.negativeSideControlVolumeIndex,
                sample.negativeSideFluidImpulseKilogramMetersPerSecond);
            exchangeSide(
                sample.positiveSideControlVolumeIndex,
                sample.positiveSideFluidImpulseKilogramMetersPerSecond);
        }
        if (settings.kinematicViscositySquareMetersPerSecond > 0.0) {
            for (std::size_t index = 0;
                 index < controlVolumes.size(); ++index) {
                controlVolumes[index].momentumKilogramMetersPerSecond = add(
                    controlVolumes[index].momentumKilogramMetersPerSecond,
                    impulses[index]);
                controlVolumes[index].velocityMetersPerSecond = scale(
                    controlVolumes[index].momentumKilogramMetersPerSecond,
                    1.0 / (densityKgPerCubicMeter
                           * controlVolumes[index].volumeCubicMeters));
            }
        }
        const double energyAfter = kineticEnergy(
            controlVolumes, densityKgPerCubicMeter);
        const double dissipation = energyBefore + wallWork - energyAfter;
        const double energyTolerance = tolerance(
            settings.absoluteEnergyToleranceJoules,
            settings.relativeEnergyTolerance,
            std::max({energyBefore, energyAfter, std::abs(wallWork)}));
        if (!std::isfinite(energyAfter) || !std::isfinite(wallWork)
            || !std::isfinite(dissipation)
            || dissipation < -energyTolerance) {
            diagnostics.failureStage =
                SceneFluidRegionWallFailureStage::Energy;
            break;
        }
        diagnostics.wallWorkOnFluidJoules += wallWork;
        diagnostics.viscousDissipationJoules += std::max(0.0, dissipation);
    }

    diagnostics.fluidMomentumAfterKilogramMetersPerSecond =
        totalMomentum(controlVolumes);
    diagnostics.fluidImpulseKilogramMetersPerSecond = subtract(
        diagnostics.fluidMomentumAfterKilogramMetersPerSecond,
        diagnostics.fluidMomentumBeforeKilogramMetersPerSecond);
    diagnostics.kineticEnergyAfterJoules = kineticEnergy(
        controlVolumes, densityKgPerCubicMeter);
    for (auto& sample : samples) {
        const auto fluidImpulse = add(
            sample.negativeSideFluidImpulseKilogramMetersPerSecond,
            sample.positiveSideFluidImpulseKilogramMetersPerSecond);
        const auto structureImpulse = scale(fluidImpulse, -1.0);
        diagnostics.structureImpulseKilogramMetersPerSecond = add(
            diagnostics.structureImpulseKilogramMetersPerSecond,
            structureImpulse);
        const auto traction = scale(
            structureImpulse,
            1.0 / (sample.areaSquareMeters * settings.timeStepSeconds));
        sample.structureTraction.tractionPascals = {
            traction.x, traction.y, traction.z};
    }
    diagnostics.momentumResidualKilogramMetersPerSecond = add(
        diagnostics.fluidImpulseKilogramMetersPerSecond,
        diagnostics.structureImpulseKilogramMetersPerSecond);
    diagnostics.momentumResidualNormKilogramMetersPerSecond = norm(
        diagnostics.momentumResidualKilogramMetersPerSecond);
    const double momentumReference = std::max(
        norm(diagnostics.fluidImpulseKilogramMetersPerSecond),
        norm(diagnostics.structureImpulseKilogramMetersPerSecond));
    if (diagnostics.failureStage == SceneFluidRegionWallFailureStage::None
        && diagnostics.momentumResidualNormKilogramMetersPerSecond
            > tolerance(
                settings.absoluteMomentumToleranceKilogramMetersPerSecond,
                settings.relativeMomentumTolerance, momentumReference)) {
        diagnostics.failureStage =
            SceneFluidRegionWallFailureStage::Conservation;
    }
    diagnostics.finite = finite(
            diagnostics.fluidMomentumBeforeKilogramMetersPerSecond)
        && finite(diagnostics.fluidMomentumAfterKilogramMetersPerSecond)
        && finite(diagnostics.fluidImpulseKilogramMetersPerSecond)
        && finite(diagnostics.structureImpulseKilogramMetersPerSecond)
        && finite(diagnostics.momentumResidualKilogramMetersPerSecond)
        && std::isfinite(
            diagnostics.momentumResidualNormKilogramMetersPerSecond)
        && std::isfinite(diagnostics.kineticEnergyBeforeJoules)
        && std::isfinite(diagnostics.kineticEnergyAfterJoules)
        && std::isfinite(diagnostics.wallWorkOnFluidJoules)
        && std::isfinite(diagnostics.viscousDissipationJoules)
        && std::ranges::all_of(controlVolumes, [](const auto& control) {
            return finite(control.velocityMetersPerSecond)
                && finite(control.momentumKilogramMetersPerSecond);
        });
    if (!diagnostics.finite
        && diagnostics.failureStage == SceneFluidRegionWallFailureStage::None) {
        diagnostics.failureStage = SceneFluidRegionWallFailureStage::NonFinite;
    }
    if (diagnostics.finite
        && diagnostics.failureStage == SceneFluidRegionWallFailureStage::None) {
        diagnostics.accepted = true;
        result.ownedStorageBytes = acceptedStorageBytes;
        result.controlVolumes = std::move(controlVolumes);
        result.samples = std::move(samples);
    }
    return result;
}

} // namespace simwing::fsi::detail
