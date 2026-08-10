#include "fluid/planar_region_fragment_opening_resistance.h"

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
        throw std::length_error("opening-resistance storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error("opening-resistance storage overflows");
    }
    return first + second;
}

bool zeroResistance(const DarcyForchheimerResistance& resistance) {
    return resistance.linearPascalSecondsPerMeter == 0.0
        && resistance.quadraticPascalSecondsSquaredPerSquareMeter == 0.0;
}

void validateResistance(const DarcyForchheimerResistance& resistance) {
    if (!std::isfinite(resistance.linearPascalSecondsPerMeter)
        || !std::isfinite(
            resistance.quadraticPascalSecondsSquaredPerSquareMeter)
        || resistance.linearPascalSecondsPerMeter < 0.0
        || resistance.quadraticPascalSecondsSquaredPerSquareMeter < 0.0) {
        throw std::invalid_argument(
            "opening resistance coefficients are invalid");
    }
}

std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
canonicalDefinitions(
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings) {
    if (definitions.size() != openings.patches.size()) {
        throw std::invalid_argument(
            "opening resistance requires exactly one definition per patch");
    }
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        canonical(definitions.begin(), definitions.end());
    std::ranges::sort(canonical, {},
        &PlanarPressureRegionFragmentOpeningResistanceDefinition::patchStableId);
    std::vector<std::uint64_t> openingPatchIds;
    openingPatchIds.reserve(openings.patches.size());
    for (const auto& patch : openings.patches)
        openingPatchIds.push_back(patch.patchStableId);
    std::ranges::sort(openingPatchIds);
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        validateResistance(canonical[index].resistance);
        if (canonical[index].patchStableId == 0
            || (index != 0
                && canonical[index - 1].patchStableId
                    == canonical[index].patchStableId)
            || canonical[index].patchStableId != openingPatchIds[index]) {
            throw std::invalid_argument(
                "opening resistance definition identity is invalid");
        }
    }
    return canonical;
}

std::uint64_t definitionFingerprint(
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        definitions) {
    Fingerprint fingerprint;
    fingerprint.integer(static_cast<std::uint64_t>(definitions.size()));
    for (const auto& definition : definitions) {
        fingerprint.integer(definition.patchStableId);
        fingerprint.real(
            definition.resistance.linearPascalSecondsPerMeter);
        fingerprint.real(
            definition.resistance
                .quadraticPascalSecondsSquaredPerSquareMeter);
    }
    return fingerprint.value();
}

PorousPlugFlowDiagnostics identityDiagnostics(
    const double velocity,
    const double density,
    const double length,
    const double area) {
    PorousPlugFlowDiagnostics result;
    result.accepted = true;
    result.velocityBeforeMetersPerSecond = velocity;
    result.midpointVelocityMetersPerSecond = velocity;
    result.velocityAfterMetersPerSecond = velocity;
    result.fluidMassKilograms = density * length * area;
    result.momentumBeforeNewtonSeconds =
        result.fluidMassKilograms * velocity;
    result.momentumAfterNewtonSeconds =
        result.momentumBeforeNewtonSeconds;
    result.kineticEnergyBeforeJoules = 0.5 * result.fluidMassKilograms
        * velocity * velocity;
    result.kineticEnergyAfterJoules = result.kineticEnergyBeforeJoules;
    if (!std::isfinite(result.fluidMassKilograms)
        || !std::isfinite(result.momentumBeforeNewtonSeconds)
        || !std::isfinite(result.kineticEnergyBeforeJoules)) {
        throw std::overflow_error(
            "zero opening resistance identity is non-finite");
    }
    return result;
}

double roundoffTolerance(const std::initializer_list<double> values) {
    double scale = 1.0;
    for (const double value : values)
        scale = std::max(scale, std::abs(value));
    return 1024.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

PlanarPressureRegionFragmentOpeningResistanceDiagnostics
advancePlanarPressureRegionFragmentOpeningResistance(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>&
        openingVelocitySamples,
    PlanarPressureRegionFragmentOpeningFluxState& openingFluxState,
    const PlanarPressureRegionFragmentOpeningResistanceSettings& settings,
    const PlanarPressureRegionFragmentOpeningResistanceLimits& limits) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || limits.maximumPatches == 0 || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening-resistance settings or limits are invalid");
    }
    validatePlanarPressureRegionFragmentOpeningFluxState(
        openingFluxState, grid, sweep, fragments, topology,
        openingDefinitions, openings, openingVelocitySamples,
        limits.openingFluxLimits);
    if (openings.patches.size() > limits.maximumPatches) {
        throw std::length_error(
            "opening-resistance patch limit exceeded");
    }
    auto canonical = canonicalDefinitions(resistanceDefinitions, openings);
    const std::size_t ownedBytes = checkedMultiply(
        openings.patches.size(),
        sizeof(
            PlanarPressureRegionFragmentOpeningResistancePatchDiagnostics));
    const std::size_t workingBytes = checkedAdd(
        checkedAdd(
            checkedMultiply(
                canonical.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningResistanceDefinition)),
            checkedMultiply(canonical.size(), sizeof(std::uint64_t))),
        checkedAdd(
            checkedMultiply(
                openingVelocitySamples.size(),
                sizeof(PlanarPressureRegionFragmentOpeningVelocitySample)),
            checkedMultiply(
                openingVelocitySamples.size(), sizeof(std::size_t))));
    if (ownedBytes > limits.maximumOwnedBytes
        || workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening-resistance storage limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningResistanceDiagnostics diagnostics;
    diagnostics.finite = true;
    diagnostics.sourceOpeningFingerprint = openings.fingerprint;
    diagnostics.sourceOpeningFluxFingerprint = openingFluxState.fingerprint;
    diagnostics.resistanceDefinitionFingerprint =
        definitionFingerprint(canonical);
    diagnostics.settings = settings;
    diagnostics.ownedStorageBytes = ownedBytes;
    diagnostics.workingStorageBytes = workingBytes;
    diagnostics.patches.reserve(openings.patches.size());

    auto candidateSamples = openingVelocitySamples;
    std::vector<std::size_t> sampleOrder(candidateSamples.size());
    std::iota(sampleOrder.begin(), sampleOrder.end(), 0);
    std::ranges::sort(sampleOrder, [&](const std::size_t first,
                                      const std::size_t second) {
        return candidateSamples[first].patchStableId
            < candidateSamples[second].patchStableId;
    });
    for (const auto& patch : openings.patches) {
        const auto definition = std::ranges::lower_bound(
            canonical, patch.patchStableId, {},
            &PlanarPressureRegionFragmentOpeningResistanceDefinition::
                patchStableId);
        const auto orderedSample = std::ranges::lower_bound(
            sampleOrder, patch.patchStableId, {},
            [&](const std::size_t index) {
                return candidateSamples[index].patchStableId;
            });
        if (definition == canonical.end()
            || definition->patchStableId != patch.patchStableId
            || orderedSample == sampleOrder.end()
            || candidateSamples[*orderedSample].patchStableId
                != patch.patchStableId) {
            throw std::logic_error(
                "opening-resistance source mapping is incomplete");
        }
        auto& sample = candidateSamples[*orderedSample];
        PorousPlugFlowDiagnostics plugFlow;
        const bool identity = zeroResistance(definition->resistance);
        if (identity) {
            ++diagnostics.zeroResistancePatchCount;
            plugFlow = identityDiagnostics(
                sample.relativeNormalVelocityMetersPerSecond,
                settings.densityKgPerCubicMeter,
                patch.centerDistanceMeters, patch.areaSquareMeters);
        } else {
            PorousPlugFlowSettings plugSettings;
            plugSettings.resistance = definition->resistance;
            plugSettings.densityKgPerCubicMeter =
                settings.densityKgPerCubicMeter;
            plugSettings.flowLengthMeters = patch.centerDistanceMeters;
            plugSettings.crossSectionAreaSquareMeters =
                patch.areaSquareMeters;
            plugSettings.drivingPressureRisePascals = 0.0;
            plugSettings.timeStepSeconds = settings.timeStepSeconds;
            plugFlow = advancePorousPlugFlow(
                sample.relativeNormalVelocityMetersPerSecond,
                plugSettings);
        }
        diagnostics.maximumAbsoluteMidpointPressureDropPascals = std::max(
            diagnostics.maximumAbsoluteMidpointPressureDropPascals,
            std::abs(plugFlow.midpointPressureDropPascals));
        diagnostics.maximumAbsoluteMomentumResidualKilogramMetersPerSecond =
            std::max(
                diagnostics
                    .maximumAbsoluteMomentumResidualKilogramMetersPerSecond,
                std::abs(plugFlow.momentumResidualNewtonSeconds));
        diagnostics.kineticEnergyBeforeJoules +=
            plugFlow.kineticEnergyBeforeJoules;
        diagnostics.kineticEnergyAfterJoules +=
            plugFlow.kineticEnergyAfterJoules;
        diagnostics.dissipatedEnergyJoules +=
            plugFlow.porousDissipationJoules;
        diagnostics.patches.push_back({
            patch.patchIndex,
            patch.patchStableId,
            patch.openingStableId,
            patch.surfaceStableId,
            patch.areaSquareMeters,
            patch.centerDistanceMeters,
            definition->resistance,
            identity,
            plugFlow,
        });
    }

    auto candidateFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        grid, sweep, fragments, topology, openingDefinitions, openings,
        candidateSamples, limits.openingFluxLimits);
    diagnostics.resultOpeningFluxFingerprint = candidateFlux.fingerprint;
    diagnostics.kineticEnergyChangeJoules =
        diagnostics.kineticEnergyAfterJoules
        - diagnostics.kineticEnergyBeforeJoules;
    diagnostics.energyResidualJoules =
        diagnostics.kineticEnergyChangeJoules
        + diagnostics.dissipatedEnergyJoules;
    diagnostics.energyToleranceJoules = roundoffTolerance({
        diagnostics.kineticEnergyBeforeJoules,
        diagnostics.kineticEnergyAfterJoules,
        diagnostics.dissipatedEnergyJoules,
    });
    diagnostics.finite = std::ranges::all_of(
        std::initializer_list<double>{
            diagnostics.maximumAbsoluteMidpointPressureDropPascals,
            diagnostics
                .maximumAbsoluteMomentumResidualKilogramMetersPerSecond,
            diagnostics.kineticEnergyBeforeJoules,
            diagnostics.kineticEnergyAfterJoules,
            diagnostics.kineticEnergyChangeJoules,
            diagnostics.dissipatedEnergyJoules,
            diagnostics.energyResidualJoules,
            diagnostics.energyToleranceJoules,
        }, [](const double value) { return std::isfinite(value); });
    diagnostics.nonIncreasingKineticEnergy = diagnostics.finite
        && diagnostics.kineticEnergyAfterJoules
            <= diagnostics.kineticEnergyBeforeJoules
                + diagnostics.energyToleranceJoules;
    diagnostics.accepted = diagnostics.finite
        && diagnostics.nonIncreasingKineticEnergy
        && diagnostics.dissipatedEnergyJoules
            >= -diagnostics.energyToleranceJoules
        && std::abs(diagnostics.energyResidualJoules)
            <= diagnostics.energyToleranceJoules;
    if (!diagnostics.accepted) return diagnostics;

    openingVelocitySamples = std::move(candidateSamples);
    openingFluxState = std::move(candidateFlux);
    return diagnostics;
}

} // namespace simwing::fsi::fluid
