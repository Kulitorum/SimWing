#include "scene_fluid_regional_opening_momentum_wall_exchange.h"

#include "scene_fluid_wall_exchange_kernel.h"

#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        integer(static_cast<std::make_unsigned_t<Underlying>>(value));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "regional opening momentum wall-exchange storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "regional opening momentum wall-exchange storage overflows");
    }
    return first * second;
}

void validateLimits(
    const SceneFluidRegionalOpeningMomentumWallExchangeLimits& limits) {
    if (limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "regional opening momentum wall-exchange limits are invalid");
    }
}

std::size_t kernelWorkingStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallInput& input) {
    std::size_t result = input.ownedStorageBytes;
    result = checkedAdd(
        result,
        checkedMultiply(
            input.controlVolumes.size(), 2 * sizeof(double)));
    result = checkedAdd(
        result,
        checkedMultiply(
            input.controlVolumes.size(), sizeof(fluid::Vector3)));
    return result;
}

std::size_t aggregateOwnedStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const detail::SceneFluidWallExchangeKernelResult& kernel) {
    return checkedAdd(input.ownedStorageBytes, kernel.ownedStorageBytes);
}

std::size_t aggregateWorkingStorageBytes(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const detail::SceneFluidWallExchangeKernelResult& kernel) {
    std::size_t result = checkedAdd(
        input.workingStorageBytes, kernelWorkingStorageBytes(input));
    result = checkedAdd(result, kernel.ownedStorageBytes);
    return result;
}

void fingerprintVector(Fingerprint& fingerprint,
                       const fluid::Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

void fingerprintSettings(Fingerprint& fingerprint,
                         const SceneFluidRegionWallSettings& settings) {
    fingerprint.real(settings.timeStepSeconds);
    fingerprint.real(settings.kinematicViscositySquareMetersPerSecond);
    fingerprint.real(settings.minimumWallDistanceMeters);
    fingerprint.real(settings.maximumViscousNumber);
    fingerprint.integer(static_cast<std::uint64_t>(settings.maximumSubsteps));
    fingerprint.real(
        settings.absoluteMomentumToleranceKilogramMetersPerSecond);
    fingerprint.real(settings.relativeMomentumTolerance);
    fingerprint.real(settings.absoluteEnergyToleranceJoules);
    fingerprint.real(settings.relativeEnergyTolerance);
}

void fingerprintDiagnostics(
    Fingerprint& fingerprint,
    const SceneFluidRegionWallDiagnostics& diagnostics) {
    for (const std::size_t value : {
             diagnostics.controlVolumeCount,
             diagnostics.quadraturePointCount,
             diagnostics.substepCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    for (const double value : {
             diagnostics.maximumFullStepViscousNumber,
             diagnostics.maximumAcceptedSubstepViscousNumber,
             diagnostics.maximumWallDistanceMeters,
             diagnostics.maximumRelativeTangentialSpeedMetersPerSecond,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.x,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.y,
             diagnostics.fluidMomentumBeforeKilogramMetersPerSecond.z,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.x,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.y,
             diagnostics.fluidMomentumAfterKilogramMetersPerSecond.z,
             diagnostics.fluidImpulseKilogramMetersPerSecond.x,
             diagnostics.fluidImpulseKilogramMetersPerSecond.y,
             diagnostics.fluidImpulseKilogramMetersPerSecond.z,
             diagnostics.structureImpulseKilogramMetersPerSecond.x,
             diagnostics.structureImpulseKilogramMetersPerSecond.y,
             diagnostics.structureImpulseKilogramMetersPerSecond.z,
             diagnostics.momentumResidualKilogramMetersPerSecond.x,
             diagnostics.momentumResidualKilogramMetersPerSecond.y,
             diagnostics.momentumResidualKilogramMetersPerSecond.z,
             diagnostics.momentumResidualNormKilogramMetersPerSecond,
             diagnostics.kineticEnergyBeforeJoules,
             diagnostics.kineticEnergyAfterJoules,
             diagnostics.wallWorkOnFluidJoules,
             diagnostics.viscousDissipationJoules}) {
        fingerprint.real(value);
    }
    fingerprint.enumeration(diagnostics.failureStage);
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.finite));
    fingerprint.integer(static_cast<std::uint8_t>(diagnostics.accepted));
}

std::uint64_t exchangeFingerprint(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange) {
    Fingerprint fingerprint;
    fingerprint.integer(exchange.version);
    fingerprint.integer(exchange.sourceWallInputFingerprint);
    fingerprintSettings(fingerprint, exchange.settings);
    fingerprint.integer(exchange.sourceInput.fingerprint);
    fingerprintDiagnostics(fingerprint, exchange.diagnostics);
    fingerprint.integer(
        static_cast<std::uint64_t>(exchange.controlVolumes.size()));
    for (const auto& control : exchange.controlVolumes) {
        fingerprint.integer(
            static_cast<std::uint64_t>(control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.real(control.volumeCubicMeters);
        fingerprint.real(control.incidentWallAreaSquareMeters);
        fingerprintVector(fingerprint, control.velocityMetersPerSecond);
        fingerprintVector(fingerprint, control.momentumKilogramMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(exchange.samples.size()));
    for (const auto& sample : exchange.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(sample.stableId);
        fingerprint.integer(sample.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.negativeSideControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.positiveSideControlVolumeIndex));
        fingerprint.real(sample.areaSquareMeters);
        fingerprintVector(
            fingerprint, sample.unitNormalNegativeToPositive);
        fingerprintVector(
            fingerprint, sample.wallVelocityMetersPerSecond);
        fingerprintVector(
            fingerprint,
            sample.negativeSideFluidImpulseKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint,
            sample.positiveSideFluidImpulseKilogramMetersPerSecond);
        fingerprint.integer(sample.structureTraction.stableId);
        fingerprint.real(sample.structureTraction.tractionPascals.x);
        fingerprint.real(sample.structureTraction.tractionPascals.y);
        fingerprint.real(sample.structureTraction.tractionPascals.z);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(exchange.ownedStorageBytes));
    fingerprint.integer(
        static_cast<std::uint64_t>(exchange.workingStorageBytes));
    return fingerprint.value();
}

detail::SceneFluidWallExchangeKernelResult rerunKernel(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const SceneFluidRegionWallSettings& settings) {
    return detail::exchangeSceneFluidWallMomentumKernel(
        input.densityKgPerCubicMeter, input.controlVolumes, input.samples,
        settings, input.ownedStorageBytes);
}

} // namespace

SceneFluidRegionalOpeningMomentumWallExchange
exchangeSceneFluidRegionalOpeningMomentumWall(
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const SceneFluidRegionWallSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallExchangeLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(input);
    if (settings.timeStepSeconds != input.timeStepSeconds) {
        throw std::invalid_argument(
            "regional opening momentum wall-exchange time step is foreign");
    }
    auto kernel = rerunKernel(input, settings);

    SceneFluidRegionalOpeningMomentumWallExchange result;
    result.sourceWallInputFingerprint = input.fingerprint;
    result.settings = settings;
    result.sourceInput = input;
    result.diagnostics = kernel.diagnostics;
    result.controlVolumes = kernel.controlVolumes;
    result.samples = kernel.samples;
    result.ownedStorageBytes = aggregateOwnedStorageBytes(input, kernel);
    result.workingStorageBytes = aggregateWorkingStorageBytes(input, kernel);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "regional opening momentum wall-exchange aggregate limit exceeded");
    }
    result.fingerprint = exchangeFingerprint(result);
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(result);
    return result;
}

void validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange) {
    validateSceneFluidRegionalOpeningMomentumWallInputIntegrity(
        exchange.sourceInput);
    if (exchange.settings.timeStepSeconds
        != exchange.sourceInput.timeStepSeconds) {
        throw std::invalid_argument(
            "regional opening momentum wall-exchange integrity is invalid");
    }
    const auto expected = rerunKernel(exchange.sourceInput, exchange.settings);
    const bool outputMatches = exchange.diagnostics == expected.diagnostics
        && exchange.controlVolumes == expected.controlVolumes
        && exchange.samples == expected.samples;
    if (exchange.version
            != sceneFluidRegionalOpeningMomentumWallExchangeVersion
        || exchange.fingerprint == 0
        || exchange.sourceWallInputFingerprint == 0
        || exchange.sourceWallInputFingerprint
            != exchange.sourceInput.fingerprint
        || !outputMatches
        || exchange.ownedStorageBytes
            != aggregateOwnedStorageBytes(exchange.sourceInput, expected)
        || exchange.workingStorageBytes
            != aggregateWorkingStorageBytes(exchange.sourceInput, expected)
        || exchange.fingerprint != exchangeFingerprint(exchange)) {
        throw std::invalid_argument(
            "regional opening momentum wall-exchange integrity is invalid");
    }
}

void validateSceneFluidRegionalOpeningMomentumWallExchange(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const SceneFluidRegionalOpeningMomentumWallInput& input,
    const SceneFluidRegionWallSettings& settings,
    const SceneFluidRegionalOpeningMomentumWallExchangeLimits& limits) {
    validateLimits(limits);
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(exchange);
    if (exchange.sourceInput != input || exchange.settings != settings
        || exchange.ownedStorageBytes > limits.maximumOwnedBytes
        || exchange.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::invalid_argument(
            "regional opening momentum wall-exchange sources are foreign");
    }
}

fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState
captureSceneFluidRegionalOpeningMomentumAdjustmentState(
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateSettings&
        settings,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    validateSceneFluidRegionalOpeningMomentumWallExchangeIntegrity(exchange);
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
        transport);
    fluid::validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        metric);
    const auto& input = exchange.sourceInput;
    if (!exchange.diagnostics.accepted
        || exchange.controlVolumes.empty()
        || exchange.controlVolumes.size() != transport.controls.size()
        || input.controlVolumes.size() != transport.controls.size()
        || input.sourceTransportFingerprint != transport.fingerprint
        || input.sourceTransportMetricFingerprint != metric.fingerprint
        || transport.targetMetricFingerprint != metric.fingerprint
        || input.densityKgPerCubicMeter
            != transport.densityKgPerCubicMeter
        || input.timeStepSeconds != transport.timeStepSeconds) {
        throw std::invalid_argument(
            "regional opening momentum adjustment source is foreign");
    }

    std::vector<
        fluid::PlanarPressureRegionFragmentOpeningMomentumTransportControl>
        adjustedControls;
    adjustedControls.reserve(transport.controls.size());
    for (std::size_t index = 0; index < transport.controls.size(); ++index) {
        const auto& source = transport.controls[index];
        const auto& inputControl = input.controlVolumes[index];
        const auto& adjusted = exchange.controlVolumes[index];
        if (source.fragmentIndex != index
            || inputControl.controlVolumeIndex != index
            || adjusted.controlVolumeIndex != index
            || inputControl.stableId != source.stableId
            || adjusted.stableId != source.stableId
            || inputControl.volumeCubicMeters != source.volumeCubicMeters
            || adjusted.volumeCubicMeters != source.volumeCubicMeters
            || inputControl.velocityMetersPerSecond
                != source.velocityMetersPerSecond
            || inputControl.momentumKilogramMetersPerSecond
                != source.momentumKilogramMetersPerSecond) {
            throw std::invalid_argument(
                "regional opening momentum adjustment control is foreign");
        }
        auto control = source;
        control.velocityMetersPerSecond = adjusted.velocityMetersPerSecond;
        control.momentumKilogramMetersPerSecond =
            adjusted.momentumKilogramMetersPerSecond;
        adjustedControls.push_back(control);
    }
    return fluid::capturePlanarPressureRegionFragmentOpeningMomentumAdjustmentState(
        transport, metric, exchange.fingerprint, adjustedControls, settings,
        limits);
}

void validateSceneFluidRegionalOpeningMomentumAdjustmentState(
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentState&
        state,
    const SceneFluidRegionalOpeningMomentumWallExchange& exchange,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumTransport&
        transport,
    const fluid::PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const fluid::PlanarPressureRegionFragmentOpeningMomentumAdjustmentStateLimits&
        limits) {
    fluid::validatePlanarPressureRegionFragmentOpeningMomentumAdjustmentStateIntegrity(
        state);
    if (state
        != captureSceneFluidRegionalOpeningMomentumAdjustmentState(
            exchange, transport, metric, state.settings, limits)) {
        throw std::invalid_argument(
            "regional opening momentum adjustment state is foreign");
    }
}

} // namespace simwing::fsi
