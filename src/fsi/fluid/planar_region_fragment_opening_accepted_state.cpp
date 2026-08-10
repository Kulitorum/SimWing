#include "fluid/planar_region_fragment_opening_accepted_state.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned> void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
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

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening accepted-state storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening accepted-state storage overflows");
    }
    return first + second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits) {
    if (limits.maximumTopologyLinkVelocities == 0
        || limits.maximumOpeningSamples == 0
        || limits.maximumPressureCorrections == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening accepted-state limits are invalid");
    }
}

void validateSettings(
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings) {
    const auto& projection = settings.projection;
    const auto& solve = projection.pressureSolve;
    if (!std::isfinite(projection.densityKgPerCubicMeter)
        || !(projection.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(projection.timeStepSeconds)
        || !(projection.timeStepSeconds > 0.0)
        || !std::isfinite(
            projection.absoluteContinuityToleranceCubicMetersPerSecond)
        || projection.absoluteContinuityToleranceCubicMetersPerSecond < 0.0
        || !std::isfinite(projection.relativeContinuityTolerance)
        || projection.relativeContinuityTolerance < 0.0
        || (projection.absoluteContinuityToleranceCubicMetersPerSecond == 0.0
            && projection.relativeContinuityTolerance == 0.0)
        || !std::isfinite(
            projection.absoluteMomentumResidualToleranceKilogramMetersPerSecond)
        || projection
               .absoluteMomentumResidualToleranceKilogramMetersPerSecond < 0.0
        || !std::isfinite(projection.relativeMomentumResidualTolerance)
        || projection.relativeMomentumResidualTolerance < 0.0
        || (projection
                    .absoluteMomentumResidualToleranceKilogramMetersPerSecond
                == 0.0
            && projection.relativeMomentumResidualTolerance == 0.0)
        || !std::isfinite(
            projection.absoluteEnergyResidualToleranceJoules)
        || projection.absoluteEnergyResidualToleranceJoules < 0.0
        || !std::isfinite(projection.relativeEnergyResidualTolerance)
        || projection.relativeEnergyResidualTolerance < 0.0
        || (projection.absoluteEnergyResidualToleranceJoules == 0.0
            && projection.relativeEnergyResidualTolerance == 0.0)
        || !std::isfinite(solve.absoluteResidualTolerancePascalsMeters)
        || solve.absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(solve.relativeResidualTolerance)
        || solve.relativeResidualTolerance < 0.0
        || (solve.absoluteResidualTolerancePascalsMeters == 0.0
            && solve.relativeResidualTolerance == 0.0)
        || !std::isfinite(
            solve.absoluteComponentCompatibilityTolerancePascalsMeters)
        || solve.absoluteComponentCompatibilityTolerancePascalsMeters < 0.0
        || solve.maximumIterations == 0) {
        throw std::invalid_argument(
            "opening accepted-state settings are invalid");
    }
}

bool allFinite(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
canonicalSamples(
    const std::span<
        const PlanarPressureRegionFragmentOpeningVelocitySample> samples) {
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample> result(
        samples.begin(), samples.end());
    std::ranges::sort(result, {},
        &PlanarPressureRegionFragmentOpeningVelocitySample::patchStableId);
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (result[index].patchStableId == 0
            || !std::isfinite(
                result[index].relativeNormalVelocityMetersPerSecond)
            || (index != 0
                && result[index - 1].patchStableId
                    == result[index].patchStableId)) {
            throw std::invalid_argument(
                "opening accepted-state samples are invalid");
        }
    }
    return result;
}

std::uint64_t resistanceDefinitionFingerprint(
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    if (definitions.size() != openings.patches.size()) {
        throw std::invalid_argument(
            "opening accepted-state resistance coverage is incomplete");
    }
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        canonical(definitions.begin(), definitions.end());
    std::ranges::sort(canonical, {},
        &PlanarPressureRegionFragmentOpeningResistanceDefinition::patchStableId);
    std::vector<std::uint64_t> patchIds;
    patchIds.reserve(openings.patches.size());
    for (const auto& patch : openings.patches)
        patchIds.push_back(patch.patchStableId);
    std::ranges::sort(patchIds);

    Fingerprint fingerprint;
    fingerprint.integer(static_cast<std::uint64_t>(canonical.size()));
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto& definition = canonical[index];
        const auto& resistance = definition.resistance;
        if (definition.patchStableId == 0
            || definition.patchStableId != patchIds[index]
            || (index != 0
                && canonical[index - 1].patchStableId
                    == definition.patchStableId)
            || !std::isfinite(resistance.linearPascalSecondsPerMeter)
            || resistance.linearPascalSecondsPerMeter < 0.0
            || !std::isfinite(
                resistance.quadraticPascalSecondsSquaredPerSquareMeter)
            || resistance.quadraticPascalSecondsSquaredPerSquareMeter < 0.0) {
            throw std::invalid_argument(
                "opening accepted-state resistance identity is invalid");
        }
        fingerprint.integer(definition.patchStableId);
        fingerprint.real(resistance.linearPascalSecondsPerMeter);
        fingerprint.real(
            resistance.quadraticPascalSecondsSquaredPerSquareMeter);
    }
    return fingerprint.value();
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state) {
    return checkedAdd(
        checkedAdd(
            checkedMultiply(
                state.orientedTopologyLinkVelocityMetersPerSecond.size(),
                sizeof(double)),
            checkedMultiply(
                state.openingVelocitySamples.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocitySample))),
        checkedAdd(
            state.openingFlux.ownedStorageBytes,
            checkedMultiply(
                state.pressureCorrectionPascals.size(), sizeof(double))));
}

std::size_t workingStorageBytes(const std::size_t fragmentCount,
                                const std::size_t sampleCount,
                                const std::size_t resistanceCount) {
    return checkedAdd(
        checkedMultiply(fragmentCount, sizeof(double)),
        checkedAdd(
            checkedMultiply(
                sampleCount,
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocitySample)),
            checkedMultiply(
                resistanceCount,
                sizeof(
                    PlanarPressureRegionFragmentOpeningResistanceDefinition))));
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings) {
    const auto& projection = settings.projection;
    const auto& solve = projection.pressureSolve;
    fingerprint.real(projection.densityKgPerCubicMeter);
    fingerprint.real(projection.timeStepSeconds);
    fingerprint.real(
        projection.absoluteContinuityToleranceCubicMetersPerSecond);
    fingerprint.real(projection.relativeContinuityTolerance);
    fingerprint.real(
        projection
            .absoluteMomentumResidualToleranceKilogramMetersPerSecond);
    fingerprint.real(projection.relativeMomentumResidualTolerance);
    fingerprint.real(projection.absoluteEnergyResidualToleranceJoules);
    fingerprint.real(projection.relativeEnergyResidualTolerance);
    fingerprint.real(solve.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(solve.relativeResidualTolerance);
    fingerprint.real(
        solve.absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(solve.maximumIterations));
    fingerprint.boolean(settings.useAuthoredPressureDrive);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    for (const std::uint64_t value : {
             state.sourcePressureOperatorFingerprint,
             state.sourceBasePressureOperatorFingerprint,
             state.sourceOpeningFingerprint,
             state.sourceFragmentFingerprint,
             state.sourceTopologyFingerprint,
             state.sourceVolumeRateFingerprint,
             state.sourceOpeningFluxFingerprint,
             state.resultOpeningFluxFingerprint,
             state.resistanceDefinitionFingerprint}) {
        fingerprint.integer(value);
    }
    fingerprintSettings(fingerprint, state.settings);
    fingerprint.integer(static_cast<std::uint64_t>(
        state.orientedTopologyLinkVelocityMetersPerSecond.size()));
    for (const double value
         : state.orientedTopologyLinkVelocityMetersPerSecond)
        fingerprint.real(value);
    fingerprint.integer(static_cast<std::uint64_t>(
        state.openingVelocitySamples.size()));
    for (const auto& sample : state.openingVelocitySamples) {
        fingerprint.integer(sample.patchStableId);
        fingerprint.real(sample.relativeNormalVelocityMetersPerSecond);
    }
    fingerprint.integer(state.openingFlux.fingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        state.pressureCorrectionPascals.size()));
    for (const double value : state.pressureCorrectionPascals)
        fingerprint.real(value);
    for (const double value : {
             state.correctedContinuityResidualL2CubicMetersPerSecond,
             state.correctedContinuityResidualMaximumCubicMetersPerSecond,
             state.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond,
             state.continuityToleranceCubicMetersPerSecond,
             state.maximumAbsoluteCorrectionVolumeMeanPascals,
             state.maximumAbsolutePressureCorrectionPascals,
             state.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
             state.momentumResidualToleranceKilogramMetersPerSecond,
             state.kineticEnergyBeforeJoules,
             state.kineticEnergyAfterJoules,
             state.kineticEnergyChangeJoules,
             state.authoredPressureWorkJoules,
             state.geometryPressureWorkJoules,
             state.correctionKineticEnergyJoules,
             state.dissipatedEnergyJoules,
             state.energyResidualJoules,
             state.energyToleranceJoules,
             state.pressureSolveFinalResidualL2PascalsMeters,
             state.pressureSolveFinalResidualMaximumPascalsMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        state.pressureSolveIterationCount));
    fingerprint.boolean(state.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(state.workingStorageBytes));
    return fingerprint.value();
}

struct DerivedEndpoint {
    double continuityL2 = 0.0;
    double continuityMaximum = 0.0;
    double componentContinuityMaximum = 0.0;
    double correctionGaugeMaximum = 0.0;
    double pressureMaximum = 0.0;
    double kineticEnergyAfter = 0.0;
};

DerivedEndpoint deriveEndpoint(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<const double> topologyVelocity,
    const PlanarPressureRegionFragmentOpeningFluxState& openingFlux,
    const std::span<const double> pressureCorrection,
    const double densityKgPerCubicMeter) {
    std::vector<double> continuity(fragments.fragments.size(), 0.0);
    DerivedEndpoint result;
    for (const auto& link : topology.links) {
        const double velocity = topologyVelocity[link.linkIndex];
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            if (velocity != 0.0) {
                throw std::invalid_argument(
                    "opening accepted-state pressure wall carries flow");
            }
            continue;
        }
        const double flow = link.areaSquareMeters * velocity;
        continuity[link.minusFragmentIndex] += flow;
        continuity[link.plusFragmentIndex] -= flow;
        const double mass = densityKgPerCubicMeter
            * link.areaSquareMeters * link.centerDistanceMeters;
        result.kineticEnergyAfter += 0.5 * mass * velocity * velocity;
    }
    for (const auto& fragment : openingFlux.fragments) {
        continuity[fragment.fragmentIndex] +=
            fragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond;
    }
    for (const auto& fragment : volumeRates.fragments) {
        continuity[fragment.fragmentIndex] +=
            fragment.geometryVolumeChangeRateCubicMetersPerSecond;
    }
    double squared = 0.0;
    for (const double residual : continuity) {
        squared += residual * residual;
        result.continuityMaximum = std::max(
            result.continuityMaximum, std::abs(residual));
    }
    result.continuityL2 = std::sqrt(
        squared / static_cast<double>(continuity.size()));
    for (const auto& component : pressureOperator.components) {
        double sum = 0.0;
        double pressureMoment = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            const std::size_t fragmentIndex =
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset];
            sum += continuity[fragmentIndex];
            pressureMoment +=
                fragments.fragments[fragmentIndex].volumeCubicMeters
                * pressureCorrection[fragmentIndex];
        }
        result.componentContinuityMaximum = std::max(
            result.componentContinuityMaximum, std::abs(sum));
        result.correctionGaugeMaximum = std::max(
            result.correctionGaugeMaximum,
            std::abs(pressureMoment / component.totalVolumeCubicMeters));
    }
    for (const double pressure : pressureCorrection)
        result.pressureMaximum = std::max(
            result.pressureMaximum, std::abs(pressure));
    for (const auto& patch : openings.patches) {
        if (patch.patchIndex >= openingFlux.patches.size()
            || openingFlux.patches[patch.patchIndex].patchStableId
                != patch.patchStableId) {
            throw std::logic_error(
                "opening accepted-state flux ordering is inconsistent");
        }
        const double velocity = openingFlux.patches[patch.patchIndex]
            .relativeNormalVelocityMetersPerSecond;
        const double mass = densityKgPerCubicMeter
            * patch.areaSquareMeters * patch.centerDistanceMeters;
        result.kineticEnergyAfter += 0.5 * mass * velocity * velocity;
    }
    if (!std::ranges::all_of(
            std::initializer_list<double>{
                result.continuityL2,
                result.continuityMaximum,
                result.componentContinuityMaximum,
                result.correctionGaugeMaximum,
                result.pressureMaximum,
                result.kineticEnergyAfter},
            [](const double value) { return std::isfinite(value); })) {
        throw std::overflow_error(
            "opening accepted-state derived endpoint is non-finite");
    }
    return result;
}

void validateSources(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits) {
    validateLimits(limits);
    validateSettings(state.settings);
    validatePlanarPressureRegionFragmentOpeningPressureOperator(
        pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, openingDefinitions, openings,
        limits.pressureOperatorLimits);
    validatePlanarPressureRegionFragmentVolumeRates(
        volumeRates, grid, sweep, fragments, topology,
        limits.volumeRateLimits);
    validatePlanarPressureRegionFragmentOpeningFluxState(
        state.openingFlux, grid, sweep, fragments, topology,
        openingDefinitions, openings, state.openingVelocitySamples,
        limits.openingFluxLimits);
    if (state.orientedTopologyLinkVelocityMetersPerSecond.size()
            != topology.links.size()
        || state.pressureCorrectionPascals.size()
            != fragments.fragments.size()
        || state.openingVelocitySamples.size() != openings.patches.size()
        || topology.links.size() > limits.maximumTopologyLinkVelocities
        || openings.patches.size() > limits.maximumOpeningSamples
        || fragments.fragments.size() > limits.maximumPressureCorrections
        || !allFinite(
            state.orientedTopologyLinkVelocityMetersPerSecond)
        || !allFinite(state.pressureCorrectionPascals)
        || state.sourcePressureOperatorFingerprint
            != pressureOperator.fingerprint
        || state.sourceBasePressureOperatorFingerprint
            != basePressureOperator.fingerprint
        || state.sourceOpeningFingerprint != openings.fingerprint
        || state.sourceFragmentFingerprint != fragments.fingerprint
        || state.sourceTopologyFingerprint != topology.fingerprint
        || state.sourceVolumeRateFingerprint != volumeRates.fingerprint
        || state.resultOpeningFluxFingerprint
            != state.openingFlux.fingerprint
        || state.resistanceDefinitionFingerprint
            != resistanceDefinitionFingerprint(
                resistanceDefinitions, openings)
        || state.settings.projection.timeStepSeconds
            != volumeRates.durationSeconds
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.workingStorageBytes != workingStorageBytes(
            fragments.fragments.size(), openings.patches.size(),
            resistanceDefinitions.size())
        || state.ownedStorageBytes > limits.maximumOwnedBytes
        || state.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::invalid_argument(
            "opening accepted-state sources are incompatible");
    }
    const auto derived = deriveEndpoint(
        pressureOperator, fragments, topology, volumeRates, openings,
        state.orientedTopologyLinkVelocityMetersPerSecond,
        state.openingFlux, state.pressureCorrectionPascals,
        state.settings.projection.densityKgPerCubicMeter);
    if (state.correctedContinuityResidualL2CubicMetersPerSecond
            != derived.continuityL2
        || state.correctedContinuityResidualMaximumCubicMetersPerSecond
            != derived.continuityMaximum
        || state.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
            != derived.componentContinuityMaximum
        || state.maximumAbsoluteCorrectionVolumeMeanPascals
            != derived.correctionGaugeMaximum
        || state.maximumAbsolutePressureCorrectionPascals
            != derived.pressureMaximum
        || state.kineticEnergyAfterJoules
            != derived.kineticEnergyAfter) {
        throw std::invalid_argument(
            "opening accepted-state derived ledgers disagree");
    }
}

} // namespace

PlanarPressureRegionFragmentOpeningAcceptedState
capturePlanarPressureRegionFragmentOpeningAcceptedState(
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningPressureStepDiagnostics&
        diagnostics,
    const std::span<const double>
        orientedTopologyLinkVelocityMetersPerSecond,
    const std::span<
        const PlanarPressureRegionFragmentOpeningVelocitySample>
        openingVelocitySamples,
    const PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    const std::span<const double> pressureCorrectionPascals,
    const PlanarPressureRegionFragmentOpeningPressureStepSettings& settings,
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits) {
    PlanarPressureRegionFragmentOpeningAcceptedState result;
    result.sourcePressureOperatorFingerprint = pressureOperator.fingerprint;
    result.sourceBasePressureOperatorFingerprint =
        basePressureOperator.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.sourceVolumeRateFingerprint = volumeRates.fingerprint;
    result.sourceOpeningFluxFingerprint =
        diagnostics.sourceOpeningFluxFingerprint;
    result.resultOpeningFluxFingerprint = openingFluxState.fingerprint;
    result.resistanceDefinitionFingerprint = resistanceDefinitionFingerprint(
        resistanceDefinitions, openings);
    result.settings = settings;
    result.orientedTopologyLinkVelocityMetersPerSecond.assign(
        orientedTopologyLinkVelocityMetersPerSecond.begin(),
        orientedTopologyLinkVelocityMetersPerSecond.end());
    result.openingVelocitySamples = canonicalSamples(openingVelocitySamples);
    result.openingFlux = openingFluxState;
    result.pressureCorrectionPascals.assign(
        pressureCorrectionPascals.begin(), pressureCorrectionPascals.end());

    const PlanarPressureRegionFragmentOpeningResistanceSettings
        expectedResistanceSettings{
            settings.projection.densityKgPerCubicMeter,
            settings.projection.timeStepSeconds,
            settings.useAuthoredPressureDrive};
    const auto& projection = diagnostics.projection;
    const auto& solve = projection.pressureSolve;
    if (!diagnostics.accepted || !diagnostics.finite
        || !diagnostics.energyAccepted
        || !diagnostics.resistance.accepted
        || !diagnostics.resistance.finite
        || !projection.accepted || !projection.finite
        || !projection.energyAccepted
        || !solve.compatible || !solve.converged || !solve.finite
        || !solve.usesOpeningPressureOperator
        || diagnostics.sourceOpeningFluxFingerprint == 0
        || diagnostics.resultOpeningFluxFingerprint
            != openingFluxState.fingerprint
        || diagnostics.resistance.sourceOpeningFingerprint
            != openings.fingerprint
        || diagnostics.resistance.sourceOpeningFluxFingerprint
            != diagnostics.sourceOpeningFluxFingerprint
        || diagnostics.resistance.resultOpeningFluxFingerprint
            != projection.predictedOpeningFluxFingerprint
        || diagnostics.resistance.resistanceDefinitionFingerprint
            != result.resistanceDefinitionFingerprint
        || diagnostics.resistance.settings
            != expectedResistanceSettings
        || diagnostics.resistance.usesAuthoredPressureDrive
            != settings.useAuthoredPressureDrive
        || projection.pressureOperatorFingerprint
            != pressureOperator.fingerprint
        || projection.basePressureOperatorFingerprint
            != basePressureOperator.fingerprint
        || projection.topologyFingerprint != topology.fingerprint
        || projection.fragmentFingerprint != fragments.fingerprint
        || projection.volumeRateFingerprint != volumeRates.fingerprint
        || projection.openingFingerprint != openings.fingerprint
        || projection.correctedOpeningFluxFingerprint
            != openingFluxState.fingerprint
        || solve.pressureOperatorFingerprint != pressureOperator.fingerprint
        || solve.basePressureOperatorFingerprint
            != basePressureOperator.fingerprint
        || solve.openingFingerprint != openings.fingerprint
        || solve.fragmentFingerprint != fragments.fingerprint
        || projection.fragmentCount != fragments.fragments.size()
        || projection.topologyLinkCount != topology.links.size()
        || projection.openingPatchCount != openings.patches.size()) {
        throw std::invalid_argument(
            "opening accepted-state step evidence is incompatible");
    }

    result.correctedContinuityResidualL2CubicMetersPerSecond =
        projection.correctedContinuityResidualL2CubicMetersPerSecond;
    result.correctedContinuityResidualMaximumCubicMetersPerSecond =
        projection.correctedContinuityResidualMaximumCubicMetersPerSecond;
    result.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond =
        projection
            .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond;
    result.continuityToleranceCubicMetersPerSecond =
        projection.continuityToleranceCubicMetersPerSecond;
    result.maximumAbsoluteCorrectionVolumeMeanPascals =
        solve.maximumAbsoluteCorrectionVolumeMeanPascals;
    for (const double value : result.pressureCorrectionPascals) {
        result.maximumAbsolutePressureCorrectionPascals = std::max(
            result.maximumAbsolutePressureCorrectionPascals,
            std::abs(value));
    }
    result.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond =
        projection
            .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond;
    result.momentumResidualToleranceKilogramMetersPerSecond =
        projection.momentumResidualToleranceKilogramMetersPerSecond;
    result.kineticEnergyBeforeJoules =
        diagnostics.kineticEnergyBeforeJoules;
    result.kineticEnergyAfterJoules =
        diagnostics.kineticEnergyAfterJoules;
    result.kineticEnergyChangeJoules =
        diagnostics.kineticEnergyChangeJoules;
    result.authoredPressureWorkJoules =
        diagnostics.authoredPressureWorkJoules;
    result.geometryPressureWorkJoules =
        diagnostics.geometryPressureWorkJoules;
    result.correctionKineticEnergyJoules =
        diagnostics.correctionKineticEnergyJoules;
    result.dissipatedEnergyJoules = diagnostics.dissipatedEnergyJoules;
    result.energyResidualJoules = diagnostics.energyResidualJoules;
    result.energyToleranceJoules = diagnostics.energyToleranceJoules;
    result.pressureSolveIterationCount = solve.iterationCount;
    result.pressureSolveFinalResidualL2PascalsMeters =
        solve.finalResidualL2PascalsMeters;
    result.pressureSolveFinalResidualMaximumPascalsMeters =
        solve.finalResidualMaximumPascalsMeters;
    result.ownedStorageBytes = ownedStorageBytes(result);
    result.workingStorageBytes = workingStorageBytes(
        fragments.fragments.size(), openings.patches.size(),
        resistanceDefinitions.size());
    result.accepted = true;

    validateSources(
        result, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits);
    if (result.kineticEnergyAfterJoules
            != projection.kineticEnergyAfterJoules
        || result.correctedContinuityResidualL2CubicMetersPerSecond
            != projection.correctedContinuityResidualL2CubicMetersPerSecond
        || result.maximumAbsoluteCorrectionVolumeMeanPascals
            != solve.maximumAbsoluteCorrectionVolumeMeanPascals) {
        throw std::invalid_argument(
            "opening accepted-state output evidence disagrees");
    }
    result.fingerprint = stateFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(result);
    return result;
}

void validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state) {
    validatePlanarPressureRegionFragmentOpeningFluxStateIntegrity(
        state.openingFlux);
    validateSettings(state.settings);
    const double reconstructedEnergyResidual =
        state.kineticEnergyChangeJoules
        - state.authoredPressureWorkJoules
        - state.geometryPressureWorkJoules
        + state.correctionKineticEnergyJoules
        + state.dissipatedEnergyJoules;
    bool sampleFluxIdentity =
        state.openingVelocitySamples.size()
        == state.openingFlux.patches.size();
    for (std::size_t index = 0;
         sampleFluxIdentity
             && index < state.openingVelocitySamples.size(); ++index) {
        const auto& sample = state.openingVelocitySamples[index];
        sampleFluxIdentity = sample.patchStableId != 0
            && std::isfinite(
                sample.relativeNormalVelocityMetersPerSecond)
            && (index == 0
                || state.openingVelocitySamples[index - 1].patchStableId
                    < sample.patchStableId)
            && sample.patchStableId
                == state.openingFlux.patches[index].patchStableId
            && sample.relativeNormalVelocityMetersPerSecond
                == state.openingFlux.patches[index]
                       .relativeNormalVelocityMetersPerSecond;
    }
    if (state.version
            != planarPressureRegionFragmentOpeningAcceptedStateVersion
        || state.fingerprint == 0 || !state.accepted
        || state.sourcePressureOperatorFingerprint == 0
        || state.sourceBasePressureOperatorFingerprint == 0
        || state.sourceOpeningFingerprint == 0
        || state.sourceFragmentFingerprint == 0
        || state.sourceTopologyFingerprint == 0
        || state.sourceVolumeRateFingerprint == 0
        || state.sourceOpeningFluxFingerprint == 0
        || state.resultOpeningFluxFingerprint == 0
        || state.resistanceDefinitionFingerprint == 0
        || state.resultOpeningFluxFingerprint
            != state.openingFlux.fingerprint
        || state.sourceOpeningFingerprint
            != state.openingFlux.sourceOpeningFingerprint
        || state.sourceFragmentFingerprint
            != state.openingFlux.sourceFragmentFingerprint
        || state.sourceTopologyFingerprint
            != state.openingFlux.sourceTopologyFingerprint
        || !sampleFluxIdentity
        || !allFinite(
            state.orientedTopologyLinkVelocityMetersPerSecond)
        || !allFinite(state.pressureCorrectionPascals)
        || !std::ranges::all_of(
            std::initializer_list<double>{
                state.correctedContinuityResidualL2CubicMetersPerSecond,
                state.correctedContinuityResidualMaximumCubicMetersPerSecond,
                state.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond,
                state.continuityToleranceCubicMetersPerSecond,
                state.maximumAbsoluteCorrectionVolumeMeanPascals,
                state.maximumAbsolutePressureCorrectionPascals,
                state.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond,
                state.momentumResidualToleranceKilogramMetersPerSecond,
                state.kineticEnergyBeforeJoules,
                state.kineticEnergyAfterJoules,
                state.kineticEnergyChangeJoules,
                state.authoredPressureWorkJoules,
                state.geometryPressureWorkJoules,
                state.correctionKineticEnergyJoules,
                state.dissipatedEnergyJoules,
                state.energyResidualJoules,
                state.energyToleranceJoules,
                state.pressureSolveFinalResidualL2PascalsMeters,
                state.pressureSolveFinalResidualMaximumPascalsMeters,
                reconstructedEnergyResidual},
            [](const double value) { return std::isfinite(value); })
        || state.correctedContinuityResidualMaximumCubicMetersPerSecond
            > state.continuityToleranceCubicMetersPerSecond
        || state.maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
            > state.continuityToleranceCubicMetersPerSecond
        || state.maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond
            > state.momentumResidualToleranceKilogramMetersPerSecond
        || state.dissipatedEnergyJoules < -state.energyToleranceJoules
        || state.kineticEnergyBeforeJoules < 0.0
        || state.kineticEnergyAfterJoules < 0.0
        || state.correctionKineticEnergyJoules < 0.0
        || state.pressureSolveIterationCount
            > state.settings.projection.pressureSolve.maximumIterations
        || state.kineticEnergyChangeJoules
            != state.kineticEnergyAfterJoules
                - state.kineticEnergyBeforeJoules
        || state.energyResidualJoules != reconstructedEnergyResidual
        || std::abs(state.energyResidualJoules)
            > state.energyToleranceJoules
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening accepted-state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningAcceptedState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& state,
    const PlanarPressureRegionFragmentOpeningPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentPressureOperator& basePressureOperator,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const PlanarPressureRegionFragmentOpeningAcceptedStateLimits& limits) {
    validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(state);
    validateSources(
        state, pressureOperator, basePressureOperator, grid, sweep, fragments,
        topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits);
}

} // namespace simwing::fsi::fluid
