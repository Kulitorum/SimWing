#include "fluid/planar_region_fragment_opening_flux.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        integer(static_cast<std::underlying_type_t<Enumeration>>(value));
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
            "planar regional fragment-opening flux storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional fragment-opening flux storage overflows");
    }
    return first + second;
}

template<typename... Counts>
std::size_t summedStorage(const Counts... counts) {
    std::size_t result = 0;
    ((result = checkedAdd(result, counts)), ...);
    return result;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningFluxState& state) {
    return summedStorage(
        checkedMultiply(
            state.patches.size(),
            sizeof(PlanarPressureRegionFragmentOpeningPatchFlux)),
        checkedMultiply(
            state.openings.size(),
            sizeof(PlanarPressureRegionFragmentOpeningFluxSummary)),
        checkedMultiply(
            state.fragments.size(),
            sizeof(PlanarPressureRegionFragmentOpeningFragmentFlux)),
        checkedMultiply(
            state.baseComponents.size(),
            sizeof(PlanarPressureRegionFragmentOpeningBaseComponentFlux)),
        checkedMultiply(
            state.connectedComponents.size(),
            sizeof(
                PlanarPressureRegionFragmentOpeningConnectedComponentFlux)));
}

std::uint64_t openingFluxFingerprint(
    const PlanarPressureRegionFragmentOpeningFluxState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceOpeningFingerprint);
    fingerprint.integer(state.sourceFragmentFingerprint);
    fingerprint.integer(state.sourceTopologyFingerprint);
    fingerprint.enumeration(state.profileAxis);
    fingerprint.integer(static_cast<std::uint64_t>(state.patches.size()));
    for (const auto& patch : state.patches) {
        for (const std::size_t value : {
                 patch.patchIndex, patch.sourceFaceLinkIndex,
                 patch.minusFragmentIndex, patch.plusFragmentIndex,
                 patch.minusBaseComponentIndex,
                 patch.plusBaseComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        for (const std::uint64_t value : {
                 patch.patchStableId, patch.openingStableId,
                 patch.surfaceStableId, patch.sourceFaceLinkStableId,
                 patch.negativeSideRegionStableId,
                 patch.positiveSideRegionStableId}) {
            fingerprint.integer(value);
        }
        fingerprint.real(patch.areaSquareMeters);
        fingerprint.real(patch.relativeNormalVelocityMetersPerSecond);
        fingerprint.real(
            patch.relativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.openings.size()));
    for (const auto& opening : state.openings) {
        fingerprint.integer(static_cast<std::uint64_t>(opening.openingIndex));
        fingerprint.integer(opening.openingStableId);
        fingerprint.integer(opening.surfaceStableId);
        fingerprint.integer(opening.negativeSideRegionStableId);
        fingerprint.integer(opening.positiveSideRegionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(opening.patchCount));
        fingerprint.real(opening.areaSquareMeters);
        fingerprint.real(
            opening.relativeVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            opening.areaWeightedRelativeNormalVelocityMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.fragments.size()));
    for (const auto& fragment : state.fragments) {
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.fragmentIndex));
        fingerprint.integer(fragment.fragmentStableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.baseComponentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.connectedComponentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.incidentOpeningPatchCount));
        fingerprint.real(
            fragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        state.baseComponents.size()));
    for (const auto& component : state.baseComponents) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.baseComponentIndex));
        fingerprint.integer(component.baseComponentStableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.connectedComponentIndex));
        fingerprint.integer(component.connectedComponentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.incidentOpeningPatchCount));
        fingerprint.real(
            component.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        state.connectedComponents.size()));
    for (const auto& component : state.connectedComponents) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.connectedComponentIndex));
        fingerprint.integer(component.connectedComponentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.baseComponentCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.openingPatchCount));
        fingerprint.real(
            component.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.real(
        state.maximumAbsoluteRelativeNormalVelocityMetersPerSecond);
    fingerprint.real(
        state.maximumAbsolutePatchVolumeFlowRateCubicMetersPerSecond);
    fingerprint.real(
        state.maximumAbsoluteBaseComponentOutwardFlowRateCubicMetersPerSecond);
    fingerprint.real(
        state.totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond);
    fingerprint.real(
        state.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond);
    fingerprint.real(
        state.maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond);
    fingerprint.real(state.conservationToleranceCubicMetersPerSecond);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(state.workingStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits) {
    if (limits.maximumPatches == 0 || limits.maximumOpenings == 0
        || limits.maximumFragments == 0
        || limits.maximumBaseComponents == 0
        || limits.maximumConnectedComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional fragment-opening flux limits are invalid");
    }
}

void addFinite(double& target, const double value) {
    target += value;
    if (!std::isfinite(target)) {
        throw std::overflow_error(
            "planar regional fragment-opening flux aggregate is non-finite");
    }
}

PlanarPressureRegionFragmentOpeningFluxState buildFluxState(
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningVelocitySample> samples,
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits) {
    validateLimits(limits);
    const std::size_t patchCount = openings.patches.size();
    const std::size_t openingCount = openings.openings.size();
    const std::size_t fragmentCount = topology.fragments.size();
    const std::size_t baseComponentCount = topology.components.size();
    const std::size_t connectedComponentCount =
        openings.connectedComponents.size();
    if (samples.size() != patchCount
        || patchCount > limits.maximumPatches
        || openingCount > limits.maximumOpenings
        || fragmentCount > limits.maximumFragments
        || baseComponentCount > limits.maximumBaseComponents
        || connectedComponentCount > limits.maximumConnectedComponents) {
        throw std::length_error(
            "planar regional fragment-opening flux count limit exceeded");
    }
    const std::size_t expectedOwnedBytes = summedStorage(
        checkedMultiply(
            patchCount,
            sizeof(PlanarPressureRegionFragmentOpeningPatchFlux)),
        checkedMultiply(
            openingCount,
            sizeof(PlanarPressureRegionFragmentOpeningFluxSummary)),
        checkedMultiply(
            fragmentCount,
            sizeof(PlanarPressureRegionFragmentOpeningFragmentFlux)),
        checkedMultiply(
            baseComponentCount,
            sizeof(PlanarPressureRegionFragmentOpeningBaseComponentFlux)),
        checkedMultiply(
            connectedComponentCount,
            sizeof(
                PlanarPressureRegionFragmentOpeningConnectedComponentFlux)));
    const std::size_t expectedWorkingBytes = checkedMultiply(
        patchCount,
        sizeof(PlanarPressureRegionFragmentOpeningVelocitySample));
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || expectedWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional fragment-opening flux byte limit exceeded");
    }

    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample> canonical(
        samples.begin(), samples.end());
    std::ranges::sort(
        canonical, {},
        &PlanarPressureRegionFragmentOpeningVelocitySample::patchStableId);
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        if (canonical[index].patchStableId == 0
            || !std::isfinite(
                canonical[index].relativeNormalVelocityMetersPerSecond)
            || (index != 0
                && canonical[index - 1].patchStableId
                    == canonical[index].patchStableId)) {
            throw std::invalid_argument(
                "planar regional fragment-opening velocity sample is invalid");
        }
    }

    PlanarPressureRegionFragmentOpeningFluxState result;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.profileAxis = openings.profileAxis;
    result.patches.reserve(patchCount);
    result.openings.resize(openingCount);
    result.fragments.resize(fragmentCount);
    result.baseComponents.resize(baseComponentCount);
    result.connectedComponents.resize(connectedComponentCount);

    for (std::size_t index = 0; index < openingCount; ++index) {
        const auto& source = openings.openings[index];
        result.openings[index] = {
            source.openingIndex,
            source.openingStableId,
            source.surfaceStableId,
            source.negativeSideRegionStableId,
            source.positiveSideRegionStableId,
            source.patchCount,
            source.areaSquareMeters,
            0.0,
            0.0,
        };
    }
    for (std::size_t index = 0; index < fragmentCount; ++index) {
        const auto& source = topology.fragments[index];
        result.fragments[index] = {
            source.fragmentIndex,
            source.fragmentStableId,
            source.regionStableId,
            source.componentIndex,
            openings.baseComponents[source.componentIndex]
                .connectedComponentIndex,
            0,
            0.0,
        };
    }
    for (std::size_t index = 0; index < baseComponentCount; ++index) {
        const auto& source = topology.components[index];
        const auto& mapping = openings.baseComponents[index];
        result.baseComponents[index] = {
            source.componentIndex,
            source.stableId,
            source.regionStableId,
            mapping.connectedComponentIndex,
            mapping.connectedComponentStableId,
            0,
            0.0,
        };
    }
    for (std::size_t index = 0; index < connectedComponentCount; ++index) {
        const auto& source = openings.connectedComponents[index];
        result.connectedComponents[index] = {
            source.componentIndex,
            source.stableId,
            source.baseComponentCount,
            source.openingPatchCount,
            0.0,
        };
    }

    for (const auto& source : openings.patches) {
        const auto sample = std::ranges::lower_bound(
            canonical, source.patchStableId, {},
            &PlanarPressureRegionFragmentOpeningVelocitySample::patchStableId);
        if (sample == canonical.end()
            || sample->patchStableId != source.patchStableId) {
            throw std::invalid_argument(
                "planar regional fragment-opening velocity sample is missing");
        }
        const double flow = source.areaSquareMeters
            * sample->relativeNormalVelocityMetersPerSecond;
        if (!std::isfinite(flow)) {
            throw std::overflow_error(
                "planar regional fragment-opening patch flow is non-finite");
        }
        result.patches.push_back({
            source.patchIndex,
            source.patchStableId,
            source.openingStableId,
            source.surfaceStableId,
            source.sourceFaceLinkIndex,
            source.sourceFaceLinkStableId,
            source.minusFragmentIndex,
            source.plusFragmentIndex,
            source.minusBaseComponentIndex,
            source.plusBaseComponentIndex,
            source.negativeSideRegionStableId,
            source.positiveSideRegionStableId,
            source.areaSquareMeters,
            sample->relativeNormalVelocityMetersPerSecond,
            flow,
        });
        result.maximumAbsoluteRelativeNormalVelocityMetersPerSecond =
            std::max(
                result.maximumAbsoluteRelativeNormalVelocityMetersPerSecond,
                std::abs(sample->relativeNormalVelocityMetersPerSecond));
        result.maximumAbsolutePatchVolumeFlowRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsolutePatchVolumeFlowRateCubicMetersPerSecond,
                std::abs(flow));
        addFinite(
            result.totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond,
            std::abs(flow));

        const auto opening = std::ranges::lower_bound(
            result.openings, source.openingStableId, {},
            &PlanarPressureRegionFragmentOpeningFluxSummary::openingStableId);
        if (opening == result.openings.end()
            || opening->openingStableId != source.openingStableId) {
            throw std::logic_error(
                "planar regional fragment-opening flux lost its opening");
        }
        addFinite(
            opening->relativeVolumeFlowRateCubicMetersPerSecond, flow);

        auto& minusFragment = result.fragments[source.minusFragmentIndex];
        auto& plusFragment = result.fragments[source.plusFragmentIndex];
        ++minusFragment.incidentOpeningPatchCount;
        ++plusFragment.incidentOpeningPatchCount;
        addFinite(
            minusFragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            flow);
        addFinite(
            plusFragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            -flow);

        auto& minusComponent =
            result.baseComponents[source.minusBaseComponentIndex];
        auto& plusComponent =
            result.baseComponents[source.plusBaseComponentIndex];
        ++minusComponent.incidentOpeningPatchCount;
        ++plusComponent.incidentOpeningPatchCount;
        addFinite(
            minusComponent.outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            flow);
        addFinite(
            plusComponent.outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            -flow);
    }

    for (auto& opening : result.openings) {
        opening.areaWeightedRelativeNormalVelocityMetersPerSecond =
            opening.relativeVolumeFlowRateCubicMetersPerSecond
            / opening.areaSquareMeters;
        if (!std::isfinite(
                opening.areaWeightedRelativeNormalVelocityMetersPerSecond)) {
            throw std::overflow_error(
                "planar regional fragment-opening summary is non-finite");
        }
    }
    for (const auto& component : result.baseComponents) {
        result.maximumAbsoluteBaseComponentOutwardFlowRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteBaseComponentOutwardFlowRateCubicMetersPerSecond,
                std::abs(
                    component.outwardRelativeVolumeFlowRateCubicMetersPerSecond));
        auto& connected = result.connectedComponents[
            component.connectedComponentIndex];
        addFinite(
            connected.outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            component.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    for (const auto& fragment : result.fragments) {
        addFinite(
            result.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond,
            fragment.outwardRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    for (const auto& component : result.connectedComponents) {
        result.maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond =
            std::max(
                result.maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond,
                std::abs(
                    component.outwardRelativeVolumeFlowRateCubicMetersPerSecond));
    }
    result.conservationToleranceCubicMetersPerSecond =
        4096.0 * std::numeric_limits<double>::epsilon()
        * result.totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond;
    if (!std::isfinite(result.conservationToleranceCubicMetersPerSecond)
        || std::abs(
               result.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond)
            > result.conservationToleranceCubicMetersPerSecond
        || result.maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond
            > result.conservationToleranceCubicMetersPerSecond) {
        throw std::overflow_error(
            "planar regional fragment-opening flux does not conserve volume");
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    result.workingStorageBytes = expectedWorkingBytes;
    if (result.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "planar regional fragment-opening flux storage is inconsistent");
    }
    result.fingerprint = openingFluxFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningFluxState
buildPlanarPressureRegionFragmentOpeningFluxState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningVelocitySample> samples,
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits) {
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, definitions,
        limits.openingLimits);
    auto result = buildFluxState(
        fragments, topology, openings, samples, limits);
    validatePlanarPressureRegionFragmentOpeningFluxState(
        result, grid, sweep, fragments, topology, definitions, openings,
        samples, limits);
    return result;
}

void validatePlanarPressureRegionFragmentOpeningFluxState(
    const PlanarPressureRegionFragmentOpeningFluxState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition> definitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const PlanarPressureRegionFragmentOpeningVelocitySample> samples,
    const PlanarPressureRegionFragmentOpeningFluxLimits& limits) {
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology, definitions,
        limits.openingLimits);
    validateLimits(limits);
    if (state.version != planarPressureRegionFragmentOpeningFluxVersion
        || state.fingerprint == 0
        || state.sourceOpeningFingerprint != openings.fingerprint
        || state.sourceFragmentFingerprint != fragments.fingerprint
        || state.sourceTopologyFingerprint != topology.fingerprint
        || state.profileAxis != topology.profileAxis
        || !std::isfinite(
            state.maximumAbsoluteRelativeNormalVelocityMetersPerSecond)
        || !std::isfinite(
            state.maximumAbsolutePatchVolumeFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            state.maximumAbsoluteBaseComponentOutwardFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            state.totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            state.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond)
        || !std::isfinite(
            state.maximumAbsoluteConnectedComponentOutwardFlowRateCubicMetersPerSecond)
        || !std::isfinite(state.conservationToleranceCubicMetersPerSecond)
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.fingerprint != openingFluxFingerprint(state)
        || state != buildFluxState(
                        fragments, topology, openings, samples, limits)) {
        throw std::invalid_argument(
            "planar regional fragment-opening flux state is invalid");
    }
}

} // namespace simwing::fsi::fluid
