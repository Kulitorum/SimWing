#include "fluid/planar_region_fragment_velocity_metric.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::uint64_t dofIdentityDomain = 0x5357'5246'564d'4554ULL;

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
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits) {
    if (limits.maximumDofs == 0 || limits.maximumFragments == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0) {
        throw std::invalid_argument(
            "planar regional velocity-metric limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional velocity-metric storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional velocity-metric storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t dofCount,
                              const std::size_t fragmentCount,
                              const std::size_t componentCount) {
    return checkedSum(
        checkedSum(
            checkedProduct(
                dofCount,
                sizeof(PlanarPressureRegionFragmentVelocityDof)),
            checkedProduct(
                fragmentCount,
                sizeof(
                    PlanarPressureRegionFragmentVelocityMetricFragment))),
        checkedProduct(
            componentCount,
            sizeof(PlanarPressureRegionFragmentVelocityMetricComponent)));
}

double vectorCoordinate(const Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "planar regional velocity-metric axis is invalid");
}

double& vectorCoordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "planar regional velocity-metric axis is invalid");
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

double closureTolerance(const double scale) {
    return 4096.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale);
}

double halfDistance(
    const PlanarPressureRegionFragmentFaceLink& link,
    const PlanarPressureRegionFragmentSet& fragments,
    const std::size_t fragmentIndex) {
    const auto& fragment = fragments.fragments.at(fragmentIndex);
    if (link.axis != fragments.axis) {
        return 0.5 * vectorCoordinate(fragments.spacingMeters, link.axis);
    }
    if (fragmentIndex == link.minusFragmentIndex) {
        return fragment.unwrappedUpperCoordinateMeters
            - fragment.unwrappedAxisCentroidMeters;
    }
    if (fragmentIndex == link.plusFragmentIndex) {
        return fragment.unwrappedAxisCentroidMeters
            - fragment.unwrappedLowerCoordinateMeters;
    }
    throw std::invalid_argument(
        "planar regional velocity-metric link misses its fragment");
}

std::uint64_t dofStableId(
    const std::uint64_t sourceLinkStableId,
    const PlanarPressureRegionFragmentVelocityDofKind kind) {
    Fingerprint fingerprint;
    fingerprint.integer(dofIdentityDomain);
    fingerprint.integer(sourceLinkStableId);
    fingerprint.enumeration(kind);
    return fingerprint.value();
}

std::uint64_t metricFingerprint(
    const PlanarPressureRegionFragmentVelocityMetric& metric) {
    Fingerprint fingerprint;
    fingerprint.integer(metric.version);
    fingerprint.integer(metric.sourceFragmentFingerprint);
    fingerprint.integer(metric.sourceTopologyFingerprint);
    fingerprint.enumeration(metric.profileAxis);
    fingerprint.integer(static_cast<std::uint64_t>(metric.dofs.size()));
    for (const auto& dof : metric.dofs) {
        fingerprint.integer(static_cast<std::uint64_t>(dof.dofIndex));
        fingerprint.integer(dof.stableId);
        fingerprint.enumeration(dof.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            dof.sourceFaceLinkIndex));
        fingerprint.integer(dof.sourceFaceLinkStableId);
        fingerprint.enumeration(dof.axis);
        fingerprint.integer(dof.surfaceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            dof.ownerFragmentIndex));
        fingerprint.integer(dof.ownerFragmentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            dof.oppositeFragmentIndex));
        fingerprint.integer(dof.oppositeFragmentStableId);
        fingerprint.integer(static_cast<std::uint64_t>(dof.componentIndex));
        fingerprint.integer(dof.regionStableId);
        for (const double value : {
                 dof.areaSquareMeters,
                 dof.ownerHalfDistanceMeters,
                 dof.oppositeHalfDistanceMeters,
                 dof.ownerDualVolumeCubicMeters,
                 dof.oppositeDualVolumeCubicMeters,
                 dof.dualVolumeCubicMeters}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(metric.fragments.size()));
    for (const auto& fragment : metric.fragments) {
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.fragmentIndex));
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            fragment.velocityDofIncidenceCount));
        for (const double value : {
                 fragment.dualVolumeByAxisCubicMeters.x,
                 fragment.dualVolumeByAxisCubicMeters.y,
                 fragment.dualVolumeByAxisCubicMeters.z,
                 fragment.sourceVolumeCubicMeters,
                 fragment.volumeClosureResidualByAxisCubicMeters.x,
                 fragment.volumeClosureResidualByAxisCubicMeters.y,
                 fragment.volumeClosureResidualByAxisCubicMeters.z}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(metric.components.size()));
    for (const auto& component : metric.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.stableId);
        fingerprint.integer(component.regionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.sharedGridDofCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.pressureLayerTraceDofCount));
        for (const double value : {
                 component.dualVolumeByAxisCubicMeters.x,
                 component.dualVolumeByAxisCubicMeters.y,
                 component.dualVolumeByAxisCubicMeters.z,
                 component.sourceVolumeCubicMeters,
                 component.volumeClosureResidualByAxisCubicMeters.x,
                 component.volumeClosureResidualByAxisCubicMeters.y,
                 component.volumeClosureResidualByAxisCubicMeters.z}) {
            fingerprint.real(value);
        }
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        metric.sharedRegionGridDofCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        metric.pressureLayerTraceDofCount));
    for (const double value : {
             metric.sharedRegionGridDualVolumeCubicMeters,
             metric.pressureLayerTraceDualVolumeCubicMeters,
             metric.totalDualVolumeCubicMeters,
             metric.dualVolumeByAxisCubicMeters.x,
             metric.dualVolumeByAxisCubicMeters.y,
             metric.dualVolumeByAxisCubicMeters.z,
             metric.domainVolumeClosureResidualByAxisCubicMeters.x,
             metric.domainVolumeClosureResidualByAxisCubicMeters.y,
             metric.domainVolumeClosureResidualByAxisCubicMeters.z,
             metric.maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
             metric.maximumAbsoluteComponentVolumeClosureResidualCubicMeters}) {
        fingerprint.real(value);
    }
    fingerprint.integer(static_cast<std::uint64_t>(metric.ownedStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentVelocityMetric buildMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentTopology(
        topology, grid, sweep, fragments, limits.topologyLimits);
    const std::size_t expectedDofs = checkedSum(
        topology.sameRegionGridLinkCount,
        checkedProduct(topology.pressureLayerWallLinkCount, 2));
    if (expectedDofs > limits.maximumDofs
        || fragments.fragments.size() > limits.maximumFragments
        || topology.components.size() > limits.maximumComponents) {
        throw std::length_error(
            "planar regional velocity-metric count limit exceeded");
    }

    PlanarPressureRegionFragmentVelocityMetric result;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.profileAxis = topology.profileAxis;
    result.ownedStorageBytes = ownedStorageBytes(
        expectedDofs, fragments.fragments.size(),
        topology.components.size());
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "planar regional velocity-metric storage limit exceeded");
    }
    result.dofs.reserve(expectedDofs);
    result.fragments.resize(fragments.fragments.size());
    for (std::size_t index = 0; index < fragments.fragments.size(); ++index) {
        const auto& source = fragments.fragments[index];
        auto& fragment = result.fragments[index];
        fragment.fragmentIndex = index;
        fragment.stableId = source.stableId;
        fragment.regionStableId = source.regionStableId;
        fragment.componentIndex = topology.fragments[index].componentIndex;
        fragment.sourceVolumeCubicMeters = source.volumeCubicMeters;
    }
    result.components.reserve(topology.components.size());
    for (const auto& source : topology.components) {
        result.components.push_back({
            source.componentIndex, source.stableId, source.regionStableId,
            0, 0, {}, source.volumeCubicMeters});
    }

    std::set<std::uint64_t> stableIds;
    const auto append = [&](
        const PlanarPressureRegionFragmentFaceLink& link,
        const PlanarPressureRegionFragmentVelocityDofKind kind,
        const std::size_t ownerFragmentIndex,
        const std::size_t oppositeFragmentIndex,
        const double ownerHalfDistanceMeters,
        const double oppositeHalfDistanceMeters,
        const std::size_t componentIndex,
        const std::uint64_t regionStableId) {
        PlanarPressureRegionFragmentVelocityDof dof;
        dof.dofIndex = result.dofs.size();
        dof.stableId = dofStableId(link.stableId, kind);
        dof.kind = kind;
        dof.sourceFaceLinkIndex = link.linkIndex;
        dof.sourceFaceLinkStableId = link.stableId;
        dof.axis = link.axis;
        dof.surfaceStableId = link.surfaceStableId;
        dof.ownerFragmentIndex = ownerFragmentIndex;
        dof.ownerFragmentStableId =
            fragments.fragments[ownerFragmentIndex].stableId;
        dof.oppositeFragmentIndex = oppositeFragmentIndex;
        dof.oppositeFragmentStableId =
            fragments.fragments[oppositeFragmentIndex].stableId;
        dof.componentIndex = componentIndex;
        dof.regionStableId = regionStableId;
        dof.areaSquareMeters = link.areaSquareMeters;
        dof.ownerHalfDistanceMeters = ownerHalfDistanceMeters;
        dof.oppositeHalfDistanceMeters = oppositeHalfDistanceMeters;
        dof.ownerDualVolumeCubicMeters =
            link.areaSquareMeters * ownerHalfDistanceMeters;
        dof.oppositeDualVolumeCubicMeters =
            link.areaSquareMeters * oppositeHalfDistanceMeters;
        dof.dualVolumeCubicMeters =
            dof.ownerDualVolumeCubicMeters
            + dof.oppositeDualVolumeCubicMeters;
        if (!stableIds.insert(dof.stableId).second
            || !std::isfinite(dof.dualVolumeCubicMeters)
            || !(dof.dualVolumeCubicMeters > 0.0)) {
            throw std::invalid_argument(
                "planar regional velocity-metric DOF is invalid");
        }
        result.dofs.push_back(dof);
    };

    for (const auto& link : topology.links) {
        const double minusHalfDistance = halfDistance(
            link, fragments, link.minusFragmentIndex);
        const double plusHalfDistance = halfDistance(
            link, fragments, link.plusFragmentIndex);
        const double reconstructedCenterDistance =
            minusHalfDistance + plusHalfDistance;
        if (!std::isfinite(minusHalfDistance)
            || !std::isfinite(plusHalfDistance)
            || !(minusHalfDistance > 0.0)
            || !(plusHalfDistance > 0.0)
            || std::abs(
                   reconstructedCenterDistance
                   - link.centerDistanceMeters)
                > closureTolerance(link.centerDistanceMeters)) {
            throw std::invalid_argument(
                "planar regional velocity-metric half distances are "
                "inconsistent");
        }
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            append(
                link,
                PlanarPressureRegionFragmentVelocityDofKind::SharedRegionGrid,
                link.minusFragmentIndex, link.plusFragmentIndex,
                minusHalfDistance, plusHalfDistance,
                link.minusComponentIndex, link.minusRegionStableId);
        } else {
            append(
                link,
                PlanarPressureRegionFragmentVelocityDofKind::
                    PressureLayerMinusTrace,
                link.minusFragmentIndex, link.plusFragmentIndex,
                minusHalfDistance, 0.0,
                link.minusComponentIndex, link.minusRegionStableId);
            append(
                link,
                PlanarPressureRegionFragmentVelocityDofKind::
                    PressureLayerPlusTrace,
                link.plusFragmentIndex, link.minusFragmentIndex,
                plusHalfDistance, 0.0,
                link.plusComponentIndex, link.plusRegionStableId);
        }
    }
    if (result.dofs.size() != expectedDofs) {
        throw std::invalid_argument(
            "planar regional velocity-metric DOF count is incomplete");
    }

    for (const auto& dof : result.dofs) {
        auto& owner = result.fragments[dof.ownerFragmentIndex];
        ++owner.velocityDofIncidenceCount;
        vectorCoordinate(owner.dualVolumeByAxisCubicMeters, dof.axis) +=
            dof.ownerDualVolumeCubicMeters;
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            auto& opposite = result.fragments[dof.oppositeFragmentIndex];
            ++opposite.velocityDofIncidenceCount;
            vectorCoordinate(
                opposite.dualVolumeByAxisCubicMeters, dof.axis) +=
                dof.oppositeDualVolumeCubicMeters;
            ++result.sharedRegionGridDofCount;
            result.sharedRegionGridDualVolumeCubicMeters +=
                dof.dualVolumeCubicMeters;
            ++result.components[dof.componentIndex].sharedGridDofCount;
        } else {
            ++result.pressureLayerTraceDofCount;
            result.pressureLayerTraceDualVolumeCubicMeters +=
                dof.dualVolumeCubicMeters;
            ++result.components[
                dof.componentIndex].pressureLayerTraceDofCount;
        }
        vectorCoordinate(
            result.components[dof.componentIndex]
                .dualVolumeByAxisCubicMeters,
            dof.axis) += dof.dualVolumeCubicMeters;
        vectorCoordinate(result.dualVolumeByAxisCubicMeters, dof.axis) +=
            dof.dualVolumeCubicMeters;
    }
    result.totalDualVolumeCubicMeters =
        result.sharedRegionGridDualVolumeCubicMeters
        + result.pressureLayerTraceDualVolumeCubicMeters;

    for (auto& fragment : result.fragments) {
        fragment.volumeClosureResidualByAxisCubicMeters = {
            fragment.dualVolumeByAxisCubicMeters.x
                - fragment.sourceVolumeCubicMeters,
            fragment.dualVolumeByAxisCubicMeters.y
                - fragment.sourceVolumeCubicMeters,
            fragment.dualVolumeByAxisCubicMeters.z
                - fragment.sourceVolumeCubicMeters,
        };
        const double residual = maximumAbsoluteComponent(
            fragment.volumeClosureResidualByAxisCubicMeters);
        result.maximumAbsoluteFragmentVolumeClosureResidualCubicMeters =
            std::max(
                result
                    .maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
                residual);
        if (fragment.velocityDofIncidenceCount != 6
            || residual > closureTolerance(fragment.sourceVolumeCubicMeters)) {
            throw std::invalid_argument(
                "planar regional velocity-metric fragment closure failed");
        }
    }
    for (auto& component : result.components) {
        component.volumeClosureResidualByAxisCubicMeters = {
            component.dualVolumeByAxisCubicMeters.x
                - component.sourceVolumeCubicMeters,
            component.dualVolumeByAxisCubicMeters.y
                - component.sourceVolumeCubicMeters,
            component.dualVolumeByAxisCubicMeters.z
                - component.sourceVolumeCubicMeters,
        };
        const double residual = maximumAbsoluteComponent(
            component.volumeClosureResidualByAxisCubicMeters);
        result.maximumAbsoluteComponentVolumeClosureResidualCubicMeters =
            std::max(
                result
                    .maximumAbsoluteComponentVolumeClosureResidualCubicMeters,
                residual);
        if (residual > closureTolerance(component.sourceVolumeCubicMeters)) {
            throw std::invalid_argument(
                "planar regional velocity-metric component closure failed");
        }
    }
    const double domainVolume = grid.cellVolumeCubicMeters()
        * static_cast<double>(grid.cellCount());
    result.domainVolumeClosureResidualByAxisCubicMeters = {
        result.dualVolumeByAxisCubicMeters.x - domainVolume,
        result.dualVolumeByAxisCubicMeters.y - domainVolume,
        result.dualVolumeByAxisCubicMeters.z - domainVolume,
    };
    if (maximumAbsoluteComponent(
            result.domainVolumeClosureResidualByAxisCubicMeters)
        > closureTolerance(domainVolume)) {
        throw std::invalid_argument(
            "planar regional velocity-metric domain closure failed");
    }
    result.fingerprint = metricFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentVelocityMetric
buildPlanarPressureRegionFragmentVelocityMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits) {
    return buildMetric(grid, sweep, fragments, topology, limits);
}

void validatePlanarPressureRegionFragmentVelocityMetric(
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits) {
    if (metric != buildMetric(grid, sweep, fragments, topology, limits)) {
        throw std::invalid_argument(
            "planar regional velocity metric is corrupted");
    }
}

} // namespace simwing::fsi::fluid
