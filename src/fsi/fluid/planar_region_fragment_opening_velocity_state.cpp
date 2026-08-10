#include "fluid/planar_region_fragment_opening_velocity_state.h"

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

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening velocity-state storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening velocity-state storage overflows");
    }
    return first * second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    if (limits.maximumSamples == 0 || limits.maximumFragments == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening velocity-state limits are invalid");
    }
}

double& coordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "opening velocity-state axis is invalid");
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 scaledVector(const Vector3& value, const double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

void addVector(Vector3& target, const Vector3& value) {
    target.x += value.x;
    target.y += value.y;
    target.z += value.z;
}

double norm(const Vector3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

double energyTolerance(const std::initializer_list<double> values) {
    double scale = 1.0;
    for (const double value : values)
        scale = std::max(scale, std::abs(value));
    return 16384.0 * std::numeric_limits<double>::epsilon() * scale;
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningVelocityState& state) {
    return checkedAdd(
        checkedMultiply(
            state.samples.size(),
            sizeof(PlanarPressureRegionFragmentOpeningVelocityStateSample)),
        checkedAdd(
            checkedMultiply(
                state.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateFragment)),
            checkedMultiply(
                state.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateComponent))));
}

std::uint64_t stateFingerprint(
    const PlanarPressureRegionFragmentOpeningVelocityState& state) {
    Fingerprint fingerprint;
    fingerprint.integer(state.version);
    fingerprint.integer(state.sourceMetricFingerprint);
    fingerprint.integer(state.sourceAcceptedStateFingerprint);
    fingerprint.integer(state.sourceVolumeRateFingerprint);
    fingerprint.boolean(state.mappedFromAcceptedEndpoint);
    fingerprint.real(state.densityKgPerCubicMeter);
    fingerprint.integer(static_cast<std::uint64_t>(state.samples.size()));
    for (const auto& sample : state.samples) {
        for (const std::size_t value : {
                 sample.dofIndex,
                 sample.sourceFaceLinkIndex,
                 sample.sourceOpeningPatchIndex,
                 sample.connectedComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(sample.stableId);
        fingerprint.enumeration(sample.kind);
        fingerprint.integer(sample.sourceFaceLinkStableId);
        fingerprint.integer(sample.sourceOpeningPatchStableId);
        fingerprint.enumeration(sample.axis);
        for (const double value : {
                 sample.dualVolumeCubicMeters,
                 sample.materialNormalVelocityMetersPerSecond,
                 sample.relativeNormalVelocityMetersPerSecond,
                 sample.normalVelocityMetersPerSecond,
                 sample.diagonalMassKilograms,
                 sample.normalMomentumKilogramMetersPerSecond,
                 sample.diagonalKineticEnergyJoules}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.fragments.size()));
    for (const auto& fragment : state.fragments) {
        for (const std::size_t value : {
                 fragment.fragmentIndex,
                 fragment.baseComponentIndex,
                 fragment.connectedComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.real(fragment.volumeCubicMeters);
        fingerprint.real(fragment.massKilograms);
        fingerprintVector(
            fingerprint, fragment.collocatedVelocityMetersPerSecond);
        fingerprintVector(
            fingerprint, fragment.momentumKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint, fragment.diagonalKineticEnergyByAxisJoules);
        fingerprint.real(fragment.diagonalKineticEnergyJoules);
        fingerprint.real(fragment.collocatedKineticEnergyJoules);
        fingerprint.real(fragment.staggeringKineticEnergyJoules);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.components.size()));
    for (const auto& component : state.components) {
        for (const std::size_t value : {
                 component.componentIndex,
                 component.baseComponentCount,
                 component.fragmentCount}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(component.stableId);
        fingerprint.real(component.volumeCubicMeters);
        fingerprint.real(component.massKilograms);
        fingerprintVector(
            fingerprint, component.momentumKilogramMetersPerSecond);
        fingerprintVector(
            fingerprint, component.diagonalKineticEnergyByAxisJoules);
        fingerprint.real(component.diagonalKineticEnergyJoules);
        fingerprint.real(component.collocatedKineticEnergyJoules);
        fingerprint.real(component.staggeringKineticEnergyJoules);
    }
    for (const std::size_t value : {
             state.sharedRegionGridSampleCount,
             state.solidWallTraceSampleCount,
             state.openingPatchSampleCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprintVector(fingerprint, state.diagonalMassByAxisKilograms);
    fingerprint.real(state.physicalMassKilograms);
    fingerprintVector(
        fingerprint, state.momentumKilogramMetersPerSecond);
    fingerprintVector(
        fingerprint, state.diagonalKineticEnergyByAxisJoules);
    for (const double value : {
             state.diagonalKineticEnergyJoules,
             state.collocatedKineticEnergyJoules,
             state.staggeringKineticEnergyJoules,
             state.maximumAbsoluteNormalVelocityMetersPerSecond,
             state.maximumCollocatedSpeedMetersPerSecond,
             state.maximumAbsoluteVelocityCompositionResidualMetersPerSecond,
             state.minimumFragmentStaggeringKineticEnergyJoules}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(state.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(state.workingStorageBytes));
    return fingerprint.value();
}

void addFragmentContribution(
    PlanarPressureRegionFragmentOpeningVelocityStateFragment& fragment,
    const GridFaceAxis axis,
    const double massKilograms,
    const double normalVelocityMetersPerSecond) {
    coordinate(fragment.momentumKilogramMetersPerSecond, axis) +=
        massKilograms * normalVelocityMetersPerSecond;
    coordinate(fragment.diagonalKineticEnergyByAxisJoules, axis) +=
        0.5 * massKilograms * normalVelocityMetersPerSecond
            * normalVelocityMetersPerSecond;
}

PlanarPressureRegionFragmentOpeningVelocityState buildState(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const std::span<const double> normalVelocityMetersPerSecond,
    const std::span<const double> materialNormalVelocityMetersPerSecond,
    const std::span<const double> relativeNormalVelocityMetersPerSecond,
    const double densityKgPerCubicMeter,
    const bool mappedFromAcceptedEndpoint,
    const std::uint64_t sourceAcceptedStateFingerprint,
    const std::uint64_t sourceVolumeRateFingerprint,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(metric);
    if (!std::isfinite(densityKgPerCubicMeter)
        || !(densityKgPerCubicMeter > 0.0)
        || normalVelocityMetersPerSecond.size() != metric.dofs.size()
        || materialNormalVelocityMetersPerSecond.size()
            != metric.dofs.size()
        || relativeNormalVelocityMetersPerSecond.size()
            != metric.dofs.size()
        || mappedFromAcceptedEndpoint
            != (sourceAcceptedStateFingerprint != 0)
        || mappedFromAcceptedEndpoint
            != (sourceVolumeRateFingerprint != 0)) {
        throw std::invalid_argument(
            "opening velocity-state source is invalid");
    }
    if (metric.dofs.size() > limits.maximumSamples
        || metric.fragments.size() > limits.maximumFragments
        || metric.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "opening velocity-state entity limit exceeded");
    }
    const std::size_t expectedOwnedBytes = checkedAdd(
        checkedMultiply(
            metric.dofs.size(),
            sizeof(PlanarPressureRegionFragmentOpeningVelocityStateSample)),
        checkedAdd(
            checkedMultiply(
                metric.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateFragment)),
            checkedMultiply(
                metric.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateComponent))));
    const std::size_t workingStorageBytes = checkedMultiply(
        checkedMultiply(metric.dofs.size(), 3), sizeof(double));
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening velocity-state storage limit exceeded");
    }

    PlanarPressureRegionFragmentOpeningVelocityState result;
    result.sourceMetricFingerprint = metric.fingerprint;
    result.sourceAcceptedStateFingerprint =
        sourceAcceptedStateFingerprint;
    result.sourceVolumeRateFingerprint = sourceVolumeRateFingerprint;
    result.mappedFromAcceptedEndpoint = mappedFromAcceptedEndpoint;
    result.densityKgPerCubicMeter = densityKgPerCubicMeter;
    result.workingStorageBytes = workingStorageBytes;
    result.samples.reserve(metric.dofs.size());
    result.fragments.resize(metric.fragments.size());
    for (std::size_t index = 0; index < metric.fragments.size(); ++index) {
        const auto& source = metric.fragments[index];
        auto& fragment = result.fragments[index];
        fragment.fragmentIndex = source.fragmentIndex;
        fragment.stableId = source.stableId;
        fragment.regionStableId = source.regionStableId;
        fragment.baseComponentIndex = source.baseComponentIndex;
        fragment.connectedComponentIndex = source.connectedComponentIndex;
        fragment.volumeCubicMeters = source.sourceVolumeCubicMeters;
        fragment.massKilograms = densityKgPerCubicMeter
            * fragment.volumeCubicMeters;
    }
    result.components.resize(metric.components.size());
    for (std::size_t index = 0; index < metric.components.size(); ++index) {
        const auto& source = metric.components[index];
        auto& component = result.components[index];
        component.componentIndex = source.componentIndex;
        component.stableId = source.stableId;
        component.baseComponentCount = source.baseComponentCount;
        component.fragmentCount = source.fragmentCount;
        component.volumeCubicMeters = source.sourceVolumeCubicMeters;
        component.massKilograms = densityKgPerCubicMeter
            * component.volumeCubicMeters;
    }

    for (std::size_t index = 0; index < metric.dofs.size(); ++index) {
        const auto& dof = metric.dofs[index];
        const double normal = normalVelocityMetersPerSecond[index];
        const double material =
            materialNormalVelocityMetersPerSecond[index];
        const double relative =
            relativeNormalVelocityMetersPerSecond[index];
        const double compositionResidual = normal - (material + relative);
        const bool gridDof = dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid;
        const bool openingDof = dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                OpeningPatch;
        if (!std::isfinite(normal) || !std::isfinite(material)
            || !std::isfinite(relative)
            || !std::isfinite(compositionResidual)
            || compositionResidual != 0.0
            || (gridDof && material != 0.0)
            || (!gridDof && !openingDof && relative != 0.0)) {
            throw std::invalid_argument(
                "opening velocity-state velocity composition is invalid");
        }
        const double mass = densityKgPerCubicMeter
            * dof.dualVolumeCubicMeters;
        const double momentum = mass * normal;
        const double energy = 0.5 * mass * normal * normal;
        if (!std::isfinite(mass) || !(mass > 0.0)
            || !std::isfinite(momentum) || !std::isfinite(energy)) {
            throw std::invalid_argument(
                "opening velocity-state inertia is invalid");
        }
        result.samples.push_back({
            index,
            dof.stableId,
            dof.kind,
            dof.sourceFaceLinkIndex,
            dof.sourceFaceLinkStableId,
            dof.sourceOpeningPatchIndex,
            dof.sourceOpeningPatchStableId,
            dof.axis,
            dof.connectedComponentIndex,
            dof.dualVolumeCubicMeters,
            material,
            relative,
            normal,
            mass,
            momentum,
            energy,
        });
        result.maximumAbsoluteNormalVelocityMetersPerSecond = std::max(
            result.maximumAbsoluteNormalVelocityMetersPerSecond,
            std::abs(normal));
        result.maximumAbsoluteVelocityCompositionResidualMetersPerSecond =
            std::max(
                result
                    .maximumAbsoluteVelocityCompositionResidualMetersPerSecond,
                std::abs(compositionResidual));
        coordinate(result.diagonalMassByAxisKilograms, dof.axis) += mass;
        coordinate(result.diagonalKineticEnergyByAxisJoules, dof.axis) +=
            energy;
        result.diagonalKineticEnergyJoules += energy;
        if (gridDof) {
            ++result.sharedRegionGridSampleCount;
        } else if (openingDof) {
            ++result.openingPatchSampleCount;
        } else {
            ++result.solidWallTraceSampleCount;
        }

        addFragmentContribution(
            result.fragments[dof.ownerFragmentIndex], dof.axis,
            densityKgPerCubicMeter * dof.ownerDualVolumeCubicMeters,
            normal);
        if (gridDof || openingDof) {
            addFragmentContribution(
                result.fragments[dof.oppositeFragmentIndex], dof.axis,
                densityKgPerCubicMeter
                    * dof.oppositeDualVolumeCubicMeters,
                normal);
        }
    }

    double minimumFragmentStaggering =
        std::numeric_limits<double>::infinity();
    for (auto& fragment : result.fragments) {
        fragment.collocatedVelocityMetersPerSecond = scaledVector(
            fragment.momentumKilogramMetersPerSecond,
            1.0 / fragment.massKilograms);
        fragment.diagonalKineticEnergyJoules =
            fragment.diagonalKineticEnergyByAxisJoules.x
            + fragment.diagonalKineticEnergyByAxisJoules.y
            + fragment.diagonalKineticEnergyByAxisJoules.z;
        const double speed = norm(
            fragment.collocatedVelocityMetersPerSecond);
        fragment.collocatedKineticEnergyJoules =
            0.5 * fragment.massKilograms * speed * speed;
        fragment.staggeringKineticEnergyJoules =
            fragment.diagonalKineticEnergyJoules
            - fragment.collocatedKineticEnergyJoules;
        if (!std::isfinite(fragment.massKilograms)
            || !(fragment.massKilograms > 0.0)
            || !finiteVector(fragment.momentumKilogramMetersPerSecond)
            || !finiteVector(fragment.collocatedVelocityMetersPerSecond)
            || !finiteVector(fragment.diagonalKineticEnergyByAxisJoules)
            || !std::isfinite(fragment.diagonalKineticEnergyJoules)
            || !std::isfinite(fragment.collocatedKineticEnergyJoules)
            || !std::isfinite(fragment.staggeringKineticEnergyJoules)
            || fragment.staggeringKineticEnergyJoules
                < -energyTolerance({
                    fragment.diagonalKineticEnergyJoules,
                    fragment.collocatedKineticEnergyJoules})) {
            throw std::invalid_argument(
                "opening velocity-state fragment energy is invalid");
        }
        minimumFragmentStaggering = std::min(
            minimumFragmentStaggering,
            fragment.staggeringKineticEnergyJoules);
        result.maximumCollocatedSpeedMetersPerSecond = std::max(
            result.maximumCollocatedSpeedMetersPerSecond, speed);
        result.physicalMassKilograms += fragment.massKilograms;
        addVector(
            result.momentumKilogramMetersPerSecond,
            fragment.momentumKilogramMetersPerSecond);
        result.collocatedKineticEnergyJoules +=
            fragment.collocatedKineticEnergyJoules;

        auto& component =
            result.components[fragment.connectedComponentIndex];
        addVector(
            component.momentumKilogramMetersPerSecond,
            fragment.momentumKilogramMetersPerSecond);
        addVector(
            component.diagonalKineticEnergyByAxisJoules,
            fragment.diagonalKineticEnergyByAxisJoules);
        component.collocatedKineticEnergyJoules +=
            fragment.collocatedKineticEnergyJoules;
    }
    result.minimumFragmentStaggeringKineticEnergyJoules =
        minimumFragmentStaggering;
    result.staggeringKineticEnergyJoules =
        result.diagonalKineticEnergyJoules
        - result.collocatedKineticEnergyJoules;
    for (auto& component : result.components) {
        component.diagonalKineticEnergyJoules =
            component.diagonalKineticEnergyByAxisJoules.x
            + component.diagonalKineticEnergyByAxisJoules.y
            + component.diagonalKineticEnergyByAxisJoules.z;
        component.staggeringKineticEnergyJoules =
            component.diagonalKineticEnergyJoules
            - component.collocatedKineticEnergyJoules;
        if (!finiteVector(component.momentumKilogramMetersPerSecond)
            || !finiteVector(component.diagonalKineticEnergyByAxisJoules)
            || !std::isfinite(component.diagonalKineticEnergyJoules)
            || !std::isfinite(component.collocatedKineticEnergyJoules)
            || !std::isfinite(component.staggeringKineticEnergyJoules)
            || component.staggeringKineticEnergyJoules
                < -energyTolerance({
                    component.diagonalKineticEnergyJoules,
                    component.collocatedKineticEnergyJoules})) {
            throw std::invalid_argument(
                "opening velocity-state component energy is invalid");
        }
    }
    if (!finiteVector(result.diagonalMassByAxisKilograms)
        || !std::isfinite(result.physicalMassKilograms)
        || !(result.physicalMassKilograms > 0.0)
        || !finiteVector(result.momentumKilogramMetersPerSecond)
        || !finiteVector(result.diagonalKineticEnergyByAxisJoules)
        || !std::isfinite(result.diagonalKineticEnergyJoules)
        || !std::isfinite(result.collocatedKineticEnergyJoules)
        || !std::isfinite(result.staggeringKineticEnergyJoules)
        || result.staggeringKineticEnergyJoules
            < -energyTolerance({
                result.diagonalKineticEnergyJoules,
                result.collocatedKineticEnergyJoules})) {
        throw std::invalid_argument(
            "opening velocity-state aggregate energy is invalid");
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "opening velocity-state preflight storage changed");
    }
    result.fingerprint = stateFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(result);
    return result;
}

PlanarPressureRegionFragmentOpeningVelocityState captureState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
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
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        metric, grid, sweep, fragments, topology, baseMetric,
        openingDefinitions, openings, limits.metricLimits);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits.acceptedStateLimits);
    if (metric.dofs.size() > limits.maximumSamples
        || metric.fragments.size() > limits.maximumFragments
        || metric.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "opening accepted velocity-state entity limit exceeded");
    }
    const std::size_t ownedBytes = checkedAdd(
        checkedMultiply(
            metric.dofs.size(),
            sizeof(PlanarPressureRegionFragmentOpeningVelocityStateSample)),
        checkedAdd(
            checkedMultiply(
                metric.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateFragment)),
            checkedMultiply(
                metric.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityStateComponent))));
    const std::size_t workingBytes = checkedMultiply(
        checkedMultiply(metric.dofs.size(), 3), sizeof(double));
    if (ownedBytes > limits.maximumOwnedBytes
        || workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening accepted velocity-state mapping storage limit exceeded");
    }
    std::vector<double> normal(metric.dofs.size(), 0.0);
    std::vector<double> material(metric.dofs.size(), 0.0);
    std::vector<double> relative(metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        const auto& link = topology.links.at(dof.sourceFaceLinkIndex);
        if (link.stableId != dof.sourceFaceLinkStableId) {
            throw std::logic_error(
                "opening accepted velocity-state link binding is invalid");
        }
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            if (link.kind
                != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
                throw std::logic_error(
                    "opening accepted velocity-state grid DOF is invalid");
            }
            relative[dof.dofIndex] = acceptedState
                .orientedTopologyLinkVelocityMetersPerSecond.at(
                    link.linkIndex);
            normal[dof.dofIndex] = relative[dof.dofIndex];
            continue;
        }
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            throw std::logic_error(
                "opening accepted velocity-state wall DOF is invalid");
        }
        const double minusMaterial = volumeRates
            .fragments.at(link.minusFragmentIndex)
            .upperBoundaryVelocityMetersPerSecond;
        const double plusMaterial = volumeRates
            .fragments.at(link.plusFragmentIndex)
            .lowerBoundaryVelocityMetersPerSecond;
        if (minusMaterial != plusMaterial) {
            throw std::invalid_argument(
                "opening accepted velocity-state wall motion disagrees");
        }
        material[dof.dofIndex] = minusMaterial;
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                OpeningPatch) {
            const auto& patch = openings.patches.at(
                dof.sourceOpeningPatchIndex);
            const auto& flux = acceptedState.openingFlux.patches.at(
                dof.sourceOpeningPatchIndex);
            if (patch.patchStableId != dof.sourceOpeningPatchStableId
                || flux.patchStableId != patch.patchStableId
                || flux.sourceFaceLinkStableId != link.stableId) {
                throw std::logic_error(
                    "opening accepted velocity-state aperture binding is invalid");
            }
            relative[dof.dofIndex] =
                flux.relativeNormalVelocityMetersPerSecond;
        }
        normal[dof.dofIndex] = material[dof.dofIndex]
            + relative[dof.dofIndex];
    }
    return buildState(
        metric, normal, material, relative,
        acceptedState.settings.projection.densityKgPerCubicMeter,
        true, acceptedState.fingerprint, volumeRates.fingerprint, limits);
}

} // namespace

PlanarPressureRegionFragmentOpeningVelocityState
buildPlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const std::span<const double> normalVelocityMetersPerSecond,
    const std::span<const double> materialNormalVelocityMetersPerSecond,
    const std::span<const double> relativeNormalVelocityMetersPerSecond,
    const double densityKgPerCubicMeter,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    return buildState(
        metric, normalVelocityMetersPerSecond,
        materialNormalVelocityMetersPerSecond,
        relativeNormalVelocityMetersPerSecond, densityKgPerCubicMeter,
        false, 0, 0, limits);
}

PlanarPressureRegionFragmentOpeningVelocityState
capturePlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
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
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    return captureState(
        acceptedState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, baseMetric, metric, limits);
}

void validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(
    const PlanarPressureRegionFragmentOpeningVelocityState& state) {
    if (state.version
            != planarPressureRegionFragmentOpeningVelocityStateVersion
        || state.fingerprint == 0 || state.sourceMetricFingerprint == 0
        || state.mappedFromAcceptedEndpoint
            != (state.sourceAcceptedStateFingerprint != 0)
        || state.mappedFromAcceptedEndpoint
            != (state.sourceVolumeRateFingerprint != 0)
        || !std::isfinite(state.densityKgPerCubicMeter)
        || !(state.densityKgPerCubicMeter > 0.0)
        || state.samples.empty() || state.fragments.empty()
        || state.components.empty()
        || state.sharedRegionGridSampleCount
                + state.solidWallTraceSampleCount
                + state.openingPatchSampleCount
            != state.samples.size()
        || !finiteVector(state.diagonalMassByAxisKilograms)
        || !std::isfinite(state.physicalMassKilograms)
        || !(state.physicalMassKilograms > 0.0)
        || !finiteVector(state.momentumKilogramMetersPerSecond)
        || !finiteVector(state.diagonalKineticEnergyByAxisJoules)
        || !std::isfinite(state.diagonalKineticEnergyJoules)
        || !std::isfinite(state.collocatedKineticEnergyJoules)
        || !std::isfinite(state.staggeringKineticEnergyJoules)
        || !std::isfinite(
            state.maximumAbsoluteNormalVelocityMetersPerSecond)
        || !std::isfinite(state.maximumCollocatedSpeedMetersPerSecond)
        || !std::isfinite(
            state.maximumAbsoluteVelocityCompositionResidualMetersPerSecond)
        || !std::isfinite(
            state.minimumFragmentStaggeringKineticEnergyJoules)
        || state.ownedStorageBytes != ownedStorageBytes(state)
        || state.workingStorageBytes == 0
        || state.fingerprint != stateFingerprint(state)) {
        throw std::invalid_argument(
            "opening velocity-state integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityState& state,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(state);
    if (state.samples.size() > limits.maximumSamples
        || state.fragments.size() > limits.maximumFragments
        || state.components.size() > limits.maximumComponents
        || state.ownedStorageBytes > limits.maximumOwnedBytes
        || state.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening velocity-state validation limit exceeded");
    }
    const std::size_t validationWorkingBytes = checkedMultiply(
        checkedMultiply(state.samples.size(), 3), sizeof(double));
    if (validationWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening velocity-state validation storage limit exceeded");
    }
    std::vector<double> normal;
    std::vector<double> material;
    std::vector<double> relative;
    normal.reserve(state.samples.size());
    material.reserve(state.samples.size());
    relative.reserve(state.samples.size());
    for (const auto& sample : state.samples) {
        normal.push_back(sample.normalVelocityMetersPerSecond);
        material.push_back(sample.materialNormalVelocityMetersPerSecond);
        relative.push_back(sample.relativeNormalVelocityMetersPerSecond);
    }
    if (state != buildState(
            metric, normal, material, relative,
            state.densityKgPerCubicMeter,
            state.mappedFromAcceptedEndpoint,
            state.sourceAcceptedStateFingerprint,
            state.sourceVolumeRateFingerprint, limits)) {
        throw std::invalid_argument(
            "opening velocity state is corrupted");
    }
}

void validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
    const PlanarPressureRegionFragmentOpeningVelocityState& state,
    const PlanarPressureRegionFragmentOpeningAcceptedState& acceptedState,
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
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PlanarPressureRegionFragmentOpeningVelocityStateLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(state);
    if (!state.mappedFromAcceptedEndpoint
        || state != captureState(
            acceptedState, pressureOperator, basePressureOperator, grid,
            sweep, fragments, topology, volumeRates, openingDefinitions,
            openings, resistanceDefinitions, baseMetric, metric, limits)) {
        throw std::invalid_argument(
            "opening accepted velocity state is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
