#include "fluid/planar_region_fragment_velocity_state.h"

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
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

void validateLimits(
    const PlanarPressureRegionFragmentVelocityStateLimits& limits) {
    if (limits.maximumSamples == 0 || limits.maximumFragments == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional velocity-state limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional velocity-state storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional velocity-state storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t sampleCount,
                              const std::size_t fragmentCount,
                              const std::size_t componentCount) {
    return checkedSum(
        checkedSum(
            checkedProduct(
                sampleCount,
                sizeof(PlanarPressureRegionFragmentVelocitySample)),
            checkedProduct(
                fragmentCount,
                sizeof(
                    PlanarPressureRegionFragmentVelocityStateFragment))),
        checkedProduct(
            componentCount,
            sizeof(
                PlanarPressureRegionFragmentVelocityStateComponent)));
}

double& vectorCoordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "planar regional velocity-state axis is invalid");
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

template<typename Ledger>
bool finiteLedger(const Ledger& ledger) {
    return finiteVector(ledger.massByAxisKilograms)
        && finiteVector(ledger.momentumKilogramMetersPerSecond)
        && finiteVector(ledger.kineticEnergyByAxisJoules)
        && std::isfinite(ledger.kineticEnergyJoules);
}

double closureTolerance(const double scale) {
    return 4096.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale);
}

template<typename Ledger>
void addContribution(Ledger& ledger,
                     const GridFaceAxis axis,
                     const double massKilograms,
                     const double velocityMetersPerSecond) {
    const double momentum =
        massKilograms * velocityMetersPerSecond;
    const double kineticEnergy = 0.5 * massKilograms
        * velocityMetersPerSecond * velocityMetersPerSecond;
    vectorCoordinate(ledger.massByAxisKilograms, axis) += massKilograms;
    vectorCoordinate(ledger.momentumKilogramMetersPerSecond, axis) +=
        momentum;
    vectorCoordinate(ledger.kineticEnergyByAxisJoules, axis) +=
        kineticEnergy;
    ledger.kineticEnergyJoules += kineticEnergy;
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentVelocityState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceMetricFingerprint);
    fingerprint.integer(state.sourceFragmentFingerprint);
    fingerprint.integer(state.sourceTopologyFingerprint);
    fingerprint.enumeration(state.profileAxis);
    fingerprint.real(state.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(state.samples.size()));
    for (const auto& sample : state.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.dofIndex));
        fingerprint.integer(sample.stableId);
        fingerprint.enumeration(sample.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.sourceFaceLinkIndex));
        fingerprint.integer(sample.sourceFaceLinkStableId);
        fingerprint.enumeration(sample.axis);
        fingerprint.integer(sample.surfaceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.componentIndex));
        fingerprint.integer(sample.regionStableId);
        fingerprint.real(sample.dualVolumeCubicMeters);
        fingerprint.real(sample.normalVelocityMetersPerSecond);
        fingerprint.real(sample.normalMomentumKilogramMetersPerSecond);
        fingerprint.real(sample.kineticEnergyJoules);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.fragments.size()));
    for (const auto& fragment : state.fragments) {
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.fragmentIndex));
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.componentIndex));
        fingerprintVector(fingerprint, fragment.massByAxisKilograms);
        fingerprintVector(
            fingerprint, fragment.momentumKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint, fragment.kineticEnergyByAxisJoules);
        fingerprint.real(fragment.kineticEnergyJoules);
        fingerprintVector(
            fingerprint, fragment.massClosureResidualByAxisKilograms);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.components.size()));
    for (const auto& component : state.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprintVector(fingerprint, component.massByAxisKilograms);
        fingerprintVector(
            fingerprint, component.momentumKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint, component.kineticEnergyByAxisJoules);
        fingerprint.real(component.kineticEnergyJoules);
        fingerprintVector(
            fingerprint, component.massClosureResidualByAxisKilograms);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        state.sharedRegionGridSampleCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        state.pressureLayerTraceSampleCount));
    fingerprintVector(fingerprint, state.massByAxisKilograms);
    fingerprintVector(
        fingerprint, state.momentumKilogramMetersPerSecond);
    fingerprintVector(fingerprint, state.kineticEnergyByAxisJoules);
    fingerprint.real(state.kineticEnergyJoules);
    fingerprint.real(state.maximumAbsoluteVelocityMetersPerSecond);
    fingerprint.real(
        state.maximumAbsoluteFragmentMassClosureResidualKilograms);
    fingerprint.real(
        state.maximumAbsoluteComponentMassClosureResidualKilograms);
    fingerprintVector(
        fingerprint, state.domainMassClosureResidualByAxisKilograms);
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentVelocityState buildState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::vector<double>& normalVelocityMetersPerSecond,
    const double densityKgPerCubicMeter,
    const PlanarPressureRegionFragmentVelocityStateLimits& limits) {
    validateLimits(limits);
    if (!std::isfinite(densityKgPerCubicMeter)
        || !(densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument(
            "planar regional velocity-state density is invalid");
    }
    validatePlanarPressureRegionFragmentVelocityMetric(
        metric, grid, sweep, fragments, topology, limits.metricLimits);
    if (normalVelocityMetersPerSecond.size() != metric.dofs.size()) {
        throw std::invalid_argument(
            "planar regional velocity-state sample count is invalid");
    }
    if (metric.dofs.size() > limits.maximumSamples
        || metric.fragments.size() > limits.maximumFragments
        || metric.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional velocity-state count limit exceeded");
    }

    PlanarPressureRegionFragmentVelocityState result;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.sourceFragmentFingerprint = metric.sourceFragmentFingerprint;
    result.sourceTopologyFingerprint = metric.sourceTopologyFingerprint;
    result.profileAxis = metric.profileAxis;
    result.densityKgPerCubicMeter = densityKgPerCubicMeter;
    result.ownedStorageBytes = ownedStorageBytes(
        metric.dofs.size(), metric.fragments.size(),
        metric.components.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional velocity-state storage limit exceeded");
    }

    result.fragments.resize(metric.fragments.size());
    for (std::size_t index = 0; index < metric.fragments.size(); ++index) {
        const auto& source = metric.fragments[index];
        auto& fragment = result.fragments[index];
        fragment.fragmentIndex = source.fragmentIndex;
        fragment.stableId = source.stableId;
        fragment.regionStableId = source.regionStableId;
        fragment.componentIndex = source.componentIndex;
    }
    result.components.resize(metric.components.size());
    for (std::size_t index = 0; index < metric.components.size(); ++index) {
        const auto& source = metric.components[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.regionStableId = source.regionStableId;
    }

    result.samples.reserve(metric.dofs.size());
    for (std::size_t index = 0; index < metric.dofs.size(); ++index) {
        const auto& dof = metric.dofs[index];
        const double velocity = normalVelocityMetersPerSecond[index];
        if (!std::isfinite(velocity)) {
            throw std::invalid_argument(
                "planar regional velocity-state sample is not finite");
        }
        const double mass =
            densityKgPerCubicMeter * dof.dualVolumeCubicMeters;
        const double momentum = mass * velocity;
        const double kineticEnergy = 0.5 * mass * velocity * velocity;
        if (!std::isfinite(mass) || !(mass > 0.0)
            || !std::isfinite(momentum)
            || !std::isfinite(kineticEnergy)) {
            throw std::invalid_argument(
                "planar regional velocity-state inertia is invalid");
        }
        result.samples.push_back({
            dof.dofIndex,
            dof.stableId,
            dof.kind,
            dof.sourceFaceLinkIndex,
            dof.sourceFaceLinkStableId,
            dof.axis,
            dof.surfaceStableId,
            dof.componentIndex,
            dof.regionStableId,
            dof.dualVolumeCubicMeters,
            velocity,
            momentum,
            kineticEnergy,
        });
        result.maximumAbsoluteVelocityMetersPerSecond = std::max(
            result.maximumAbsoluteVelocityMetersPerSecond,
            std::abs(velocity));
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            ++result.sharedRegionGridSampleCount;
        } else {
            ++result.pressureLayerTraceSampleCount;
        }

        addContribution(
            result, dof.axis, mass, velocity);
        addContribution(
            result.components[dof.componentIndex], dof.axis, mass,
            velocity);
        addContribution(
            result.fragments[dof.ownerFragmentIndex], dof.axis,
            densityKgPerCubicMeter * dof.ownerDualVolumeCubicMeters,
            velocity);
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            addContribution(
                result.fragments[dof.oppositeFragmentIndex], dof.axis,
                densityKgPerCubicMeter
                    * dof.oppositeDualVolumeCubicMeters,
                velocity);
        }
    }

    for (std::size_t index = 0; index < result.fragments.size(); ++index) {
        auto& fragment = result.fragments[index];
        if (!finiteLedger(fragment)) {
            throw std::invalid_argument(
                "planar regional velocity-state fragment ledger is "
                "invalid");
        }
        const double expectedMass = densityKgPerCubicMeter
            * metric.fragments[index].sourceVolumeCubicMeters;
        fragment.massClosureResidualByAxisKilograms = {
            fragment.massByAxisKilograms.x - expectedMass,
            fragment.massByAxisKilograms.y - expectedMass,
            fragment.massByAxisKilograms.z - expectedMass,
        };
        const double residual = maximumAbsoluteComponent(
            fragment.massClosureResidualByAxisKilograms);
        result.maximumAbsoluteFragmentMassClosureResidualKilograms =
            std::max(
                result
                    .maximumAbsoluteFragmentMassClosureResidualKilograms,
                residual);
        if (residual > closureTolerance(expectedMass)) {
            throw std::invalid_argument(
                "planar regional velocity-state fragment mass closure "
                "failed");
        }
    }
    for (std::size_t index = 0; index < result.components.size(); ++index) {
        auto& component = result.components[index];
        if (!finiteLedger(component)) {
            throw std::invalid_argument(
                "planar regional velocity-state component ledger is "
                "invalid");
        }
        const double expectedMass = densityKgPerCubicMeter
            * metric.components[index].sourceVolumeCubicMeters;
        component.massClosureResidualByAxisKilograms = {
            component.massByAxisKilograms.x - expectedMass,
            component.massByAxisKilograms.y - expectedMass,
            component.massByAxisKilograms.z - expectedMass,
        };
        const double residual = maximumAbsoluteComponent(
            component.massClosureResidualByAxisKilograms);
        result.maximumAbsoluteComponentMassClosureResidualKilograms =
            std::max(
                result
                    .maximumAbsoluteComponentMassClosureResidualKilograms,
                residual);
        if (residual > closureTolerance(expectedMass)) {
            throw std::invalid_argument(
                "planar regional velocity-state component mass closure "
                "failed");
        }
    }
    const double domainMass = densityKgPerCubicMeter
        * grid.cellVolumeCubicMeters() * static_cast<double>(grid.cellCount());
    result.domainMassClosureResidualByAxisKilograms = {
        result.massByAxisKilograms.x - domainMass,
        result.massByAxisKilograms.y - domainMass,
        result.massByAxisKilograms.z - domainMass,
    };
    if (!finiteLedger(result)
        || !std::isfinite(
            result.maximumAbsoluteVelocityMetersPerSecond)
        || maximumAbsoluteComponent(
            result.domainMassClosureResidualByAxisKilograms)
        > closureTolerance(domainMass)) {
        throw std::invalid_argument(
            "planar regional velocity-state domain mass closure failed");
    }
    result.fingerprint = stateFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentVelocityState
buildPlanarPressureRegionFragmentVelocityState(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::vector<double>& normalVelocityMetersPerSecond,
    const double densityKgPerCubicMeter,
    const PlanarPressureRegionFragmentVelocityStateLimits& limits) {
    return buildState(
        grid, sweep, fragments, topology, metric,
        normalVelocityMetersPerSecond, densityKgPerCubicMeter, limits);
}

void validatePlanarPressureRegionFragmentVelocityState(
    const PlanarPressureRegionFragmentVelocityState& state,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PlanarPressureRegionFragmentVelocityStateLimits& limits) {
    validateLimits(limits);
    if (state.samples.size() > limits.maximumSamples
        || state.fragments.size() > limits.maximumFragments
        || state.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional velocity-state validation limit exceeded");
    }
    if (checkedProduct(state.samples.size(), sizeof(double))
        > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional velocity-state validation storage limit "
            "exceeded");
    }
    std::vector<double> velocity;
    velocity.reserve(state.samples.size());
    for (const auto& sample : state.samples) {
        velocity.push_back(sample.normalVelocityMetersPerSecond);
    }
    if (state != buildState(
                     grid, sweep, fragments, topology, metric, velocity,
                     state.densityKgPerCubicMeter, limits)) {
        throw std::invalid_argument(
            "planar regional velocity state is corrupted");
    }
}

} // namespace simwing::fsi::fluid
