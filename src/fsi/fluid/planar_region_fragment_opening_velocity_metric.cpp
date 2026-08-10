#include "fluid/planar_region_fragment_opening_velocity_metric.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;
constexpr std::size_t missingIndex = std::numeric_limits<std::size_t>::max();

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

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "opening velocity-metric storage overflows");
    }
    return first + second;
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "opening velocity-metric storage overflows");
    }
    return first * second;
}

void validateLimits(
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits) {
    if (limits.maximumDofs == 0 || limits.maximumFragments == 0
        || limits.maximumComponents == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "opening velocity-metric limits are invalid");
    }
}

double& coordinate(Vector3& value, const GridFaceAxis axis) {
    switch (axis) {
    case GridFaceAxis::X: return value.x;
    case GridFaceAxis::Y: return value.y;
    case GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "opening velocity-metric axis is invalid");
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({std::abs(value.x), std::abs(value.y),
                     std::abs(value.z)});
}

double closureTolerance(const double scale) {
    return 8192.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, std::abs(scale));
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t dofStableId(
    const PlanarPressureRegionFragmentOpeningVelocityDofKind kind,
    const std::uint64_t sourceStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint32_t{0x4f564d31U});
    fingerprint.enumeration(kind);
    fingerprint.integer(sourceStableId);
    return fingerprint.value();
}

std::size_t ownedStorageBytes(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric) {
    return checkedAdd(
        checkedMultiply(
            metric.dofs.size(),
            sizeof(PlanarPressureRegionFragmentOpeningVelocityDof)),
        checkedAdd(
            checkedMultiply(
                metric.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityMetricFragment)),
            checkedMultiply(
                metric.components.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityMetricComponent))));
}

std::uint64_t metricFingerprint(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric) {
    Fingerprint fingerprint;
    fingerprint.integer(metric.version);
    fingerprint.integer(metric.sourceBaseMetricFingerprint);
    fingerprint.integer(metric.sourceOpeningFingerprint);
    fingerprint.integer(metric.sourceFragmentFingerprint);
    fingerprint.integer(metric.sourceTopologyFingerprint);
    fingerprint.enumeration(metric.profileAxis);
    fingerprint.integer(static_cast<std::uint64_t>(metric.dofs.size()));
    for (const auto& dof : metric.dofs) {
        for (const std::size_t value : {
                 dof.dofIndex,
                 dof.sourceFaceLinkIndex,
                 dof.sourceOpeningPatchIndex,
                 dof.ownerFragmentIndex,
                 dof.oppositeFragmentIndex,
                 dof.ownerBaseComponentIndex,
                 dof.oppositeBaseComponentIndex,
                 dof.connectedComponentIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(dof.stableId);
        fingerprint.enumeration(dof.kind);
        fingerprint.integer(dof.sourceFaceLinkStableId);
        fingerprint.integer(dof.sourceOpeningPatchStableId);
        fingerprint.enumeration(dof.axis);
        fingerprint.integer(dof.surfaceStableId);
        fingerprint.integer(dof.ownerFragmentStableId);
        fingerprint.integer(dof.oppositeFragmentStableId);
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
        for (const std::size_t value : {
                 fragment.fragmentIndex,
                 fragment.baseComponentIndex,
                 fragment.connectedComponentIndex,
                 fragment.velocityDofIncidenceCount}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(fragment.stableId);
        fingerprint.integer(fragment.regionStableId);
        fingerprintVector(
            fingerprint, fragment.dualVolumeByAxisCubicMeters);
        fingerprint.real(fragment.sourceVolumeCubicMeters);
        fingerprintVector(
            fingerprint, fragment.volumeClosureResidualByAxisCubicMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(metric.components.size()));
    for (const auto& component : metric.components) {
        for (const std::size_t value : {
                 component.componentIndex,
                 component.baseComponentCount,
                 component.fragmentCount,
                 component.velocityDofIncidenceCount}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(component.stableId);
        fingerprintVector(
            fingerprint, component.dualVolumeByAxisCubicMeters);
        fingerprint.real(component.sourceVolumeCubicMeters);
        fingerprintVector(
            fingerprint, component.volumeClosureResidualByAxisCubicMeters);
    }
    for (const std::size_t value : {
             metric.sharedRegionGridDofCount,
             metric.solidWallTraceDofCount,
             metric.openingPatchDofCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    for (const double value : {
             metric.sharedRegionGridDualVolumeCubicMeters,
             metric.solidWallTraceDualVolumeCubicMeters,
             metric.openingPatchDualVolumeCubicMeters,
             metric.totalDualVolumeCubicMeters,
             metric.totalPressureWallAreaSquareMeters,
             metric.totalSolidWallAreaSquareMeters,
             metric.totalOpeningAreaSquareMeters,
             metric.wallAreaPartitionResidualSquareMeters}) {
        fingerprint.real(value);
    }
    fingerprintVector(fingerprint, metric.dualVolumeByAxisCubicMeters);
    fingerprintVector(
        fingerprint, metric.domainVolumeClosureResidualByAxisCubicMeters);
    fingerprint.real(
        metric.maximumAbsoluteFragmentVolumeClosureResidualCubicMeters);
    fingerprint.real(
        metric.maximumAbsoluteComponentVolumeClosureResidualCubicMeters);
    fingerprint.integer(static_cast<std::uint64_t>(metric.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(metric.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentOpeningVelocityMetric buildMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentVelocityMetric(
        baseMetric, grid, sweep, fragments, topology,
        limits.baseMetricLimits);
    validatePlanarPressureRegionFragmentOpenings(
        openings, grid, sweep, fragments, topology,
        openingDefinitions, limits.openingLimits);
    if (fragments.fragments.size() > limits.maximumFragments
        || openings.connectedComponents.size() > limits.maximumComponents) {
        throw std::length_error(
            "opening velocity-metric entity limit exceeded");
    }
    const std::size_t fullyOpenWallCount = static_cast<std::size_t>(
        std::ranges::count_if(
            openings.partitions,
            [](const auto& partition) {
                return partition.solidAreaSquareMeters == 0.0;
            }));
    if (fullyOpenWallCount > topology.pressureLayerWallLinkCount) {
        throw std::logic_error(
            "opening velocity-metric fully open wall count is invalid");
    }
    const std::size_t solidTraceDofCount = checkedMultiply(
        topology.pressureLayerWallLinkCount - fullyOpenWallCount, 2);
    const std::size_t dofCount = checkedAdd(
        checkedAdd(
            topology.sameRegionGridLinkCount,
            solidTraceDofCount),
        openings.patches.size());
    if (dofCount > limits.maximumDofs) {
        throw std::length_error(
            "opening velocity-metric DOF limit exceeded");
    }
    const std::size_t maximumOwnedBytes = checkedAdd(
        checkedMultiply(
            dofCount,
            sizeof(PlanarPressureRegionFragmentOpeningVelocityDof)),
        checkedAdd(
            checkedMultiply(
                fragments.fragments.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityMetricFragment)),
            checkedMultiply(
                openings.connectedComponents.size(),
                sizeof(
                    PlanarPressureRegionFragmentOpeningVelocityMetricComponent))));
    if (maximumOwnedBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "opening velocity-metric owned storage limit exceeded");
    }
    const std::size_t linkMapCount = checkedMultiply(
        topology.links.size(), 4);
    const std::size_t workingBytes = checkedMultiply(
        checkedAdd(linkMapCount, openings.patches.size()),
        sizeof(std::size_t));
    if (workingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening velocity-metric working storage limit exceeded");
    }

    std::vector<std::size_t> sharedByLink(
        topology.links.size(), missingIndex);
    std::vector<std::size_t> minusTraceByLink(
        topology.links.size(), missingIndex);
    std::vector<std::size_t> plusTraceByLink(
        topology.links.size(), missingIndex);
    for (const auto& dof : baseMetric.dofs) {
        if (dof.sourceFaceLinkIndex >= topology.links.size()) {
            throw std::logic_error(
                "opening velocity-metric base DOF link is invalid");
        }
        auto* target = &sharedByLink[dof.sourceFaceLinkIndex];
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                PressureLayerMinusTrace) {
            target = &minusTraceByLink[dof.sourceFaceLinkIndex];
        } else if (dof.kind
                   == PlanarPressureRegionFragmentVelocityDofKind::
                       PressureLayerPlusTrace) {
            target = &plusTraceByLink[dof.sourceFaceLinkIndex];
        }
        if (*target != missingIndex) {
            throw std::logic_error(
                "opening velocity-metric base DOF mapping is duplicated");
        }
        *target = dof.dofIndex;
    }
    std::vector<std::size_t> partitionByLink(
        topology.links.size(), missingIndex);
    for (const auto& partition : openings.partitions) {
        if (partition.sourceFaceLinkIndex >= topology.links.size()
            || partitionByLink[partition.sourceFaceLinkIndex]
                != missingIndex) {
            throw std::logic_error(
                "opening velocity-metric wall partition mapping is invalid");
        }
        partitionByLink[partition.sourceFaceLinkIndex] =
            partition.partitionIndex;
    }
    std::vector<std::size_t> patchOrder(openings.patches.size());
    std::iota(patchOrder.begin(), patchOrder.end(), 0);
    std::ranges::sort(patchOrder, [&](const std::size_t first,
                                     const std::size_t second) {
        const auto& lhs = openings.patches[first];
        const auto& rhs = openings.patches[second];
        return std::pair{lhs.sourceFaceLinkIndex, lhs.patchStableId}
            < std::pair{rhs.sourceFaceLinkIndex, rhs.patchStableId};
    });

    PlanarPressureRegionFragmentOpeningVelocityMetric result;
    result.sourceBaseMetricFingerprint = baseMetric.fingerprint;
    result.sourceOpeningFingerprint = openings.fingerprint;
    result.sourceFragmentFingerprint = fragments.fingerprint;
    result.sourceTopologyFingerprint = topology.fingerprint;
    result.profileAxis = topology.profileAxis;
    result.workingStorageBytes = workingBytes;
    result.dofs.reserve(dofCount);
    result.fragments.resize(fragments.fragments.size());
    for (std::size_t index = 0; index < fragments.fragments.size(); ++index) {
        const auto& source = fragments.fragments[index];
        const std::size_t baseComponentIndex =
            topology.fragments[index].componentIndex;
        if (baseComponentIndex >= openings.baseComponents.size()) {
            throw std::logic_error(
                "opening velocity-metric base component is missing");
        }
        result.fragments[index] = {
            index,
            source.stableId,
            source.regionStableId,
            baseComponentIndex,
            openings.baseComponents[baseComponentIndex]
                .connectedComponentIndex,
            0,
            {},
            source.volumeCubicMeters,
            {},
        };
    }
    result.components.reserve(openings.connectedComponents.size());
    for (const auto& source : openings.connectedComponents) {
        result.components.push_back({
            source.componentIndex,
            source.stableId,
            source.baseComponentCount,
            source.fragmentCount,
            0,
            {},
            source.volumeCubicMeters,
            {},
        });
    }

    std::set<std::uint64_t> stableIds;
    const auto append = [&result, &stableIds](
        const PlanarPressureRegionFragmentOpeningVelocityDofKind kind,
        const PlanarPressureRegionFragmentFaceLink& link,
        const PlanarPressureRegionFragmentVelocityDof& minusBaseDof,
        const PlanarPressureRegionFragmentVelocityDof& plusBaseDof,
        const std::size_t sourceOpeningPatchIndex,
        const std::uint64_t sourceOpeningPatchStableId,
        const double areaSquareMeters,
        const bool hasOpposite) {
        const std::size_t ownerFragmentIndex =
            kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        SolidWallPlusTrace
                ? link.plusFragmentIndex
                : link.minusFragmentIndex;
        const std::size_t oppositeFragmentIndex =
            ownerFragmentIndex == link.minusFragmentIndex
            ? link.plusFragmentIndex : link.minusFragmentIndex;
        const auto& ownerBase =
            ownerFragmentIndex == link.minusFragmentIndex
            ? minusBaseDof : plusBaseDof;
        const auto& oppositeBase =
            ownerFragmentIndex == link.minusFragmentIndex
            ? plusBaseDof : minusBaseDof;
        const auto& ownerFragment = result.fragments[ownerFragmentIndex];
        const auto& oppositeFragment =
            result.fragments[oppositeFragmentIndex];
        if (hasOpposite
            && ownerFragment.connectedComponentIndex
                != oppositeFragment.connectedComponentIndex) {
            throw std::logic_error(
                "opening velocity-metric shared DOF crosses connected components");
        }
        const std::uint64_t sourceStableId =
            sourceOpeningPatchStableId != 0
            ? sourceOpeningPatchStableId : link.stableId;
        PlanarPressureRegionFragmentOpeningVelocityDof dof;
        dof.dofIndex = result.dofs.size();
        dof.stableId = dofStableId(kind, sourceStableId);
        dof.kind = kind;
        dof.sourceFaceLinkIndex = link.linkIndex;
        dof.sourceFaceLinkStableId = link.stableId;
        dof.sourceOpeningPatchIndex = sourceOpeningPatchIndex;
        dof.sourceOpeningPatchStableId = sourceOpeningPatchStableId;
        dof.axis = link.axis;
        dof.surfaceStableId = link.surfaceStableId;
        dof.ownerFragmentIndex = ownerFragmentIndex;
        dof.ownerFragmentStableId = ownerFragment.stableId;
        dof.oppositeFragmentIndex = oppositeFragmentIndex;
        dof.oppositeFragmentStableId = oppositeFragment.stableId;
        dof.ownerBaseComponentIndex = ownerFragment.baseComponentIndex;
        dof.oppositeBaseComponentIndex =
            oppositeFragment.baseComponentIndex;
        dof.connectedComponentIndex =
            ownerFragment.connectedComponentIndex;
        dof.areaSquareMeters = areaSquareMeters;
        dof.ownerHalfDistanceMeters = ownerBase.ownerHalfDistanceMeters;
        dof.oppositeHalfDistanceMeters = hasOpposite
            ? (kind
                       == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                           SharedRegionGrid
                   ? ownerBase.oppositeHalfDistanceMeters
                   : oppositeBase.ownerHalfDistanceMeters)
            : 0.0;
        dof.ownerDualVolumeCubicMeters =
            areaSquareMeters * dof.ownerHalfDistanceMeters;
        dof.oppositeDualVolumeCubicMeters =
            areaSquareMeters * dof.oppositeHalfDistanceMeters;
        dof.dualVolumeCubicMeters = dof.ownerDualVolumeCubicMeters
            + dof.oppositeDualVolumeCubicMeters;
        if (!stableIds.insert(dof.stableId).second
            || !std::isfinite(dof.areaSquareMeters)
            || !(dof.areaSquareMeters > 0.0)
            || !std::isfinite(dof.dualVolumeCubicMeters)
            || !(dof.dualVolumeCubicMeters > 0.0)) {
            throw std::invalid_argument(
                "opening velocity-metric DOF is invalid");
        }
        result.dofs.push_back(dof);
    };

    std::size_t patchOffset = 0;
    for (const auto& link : topology.links) {
        const std::size_t sharedIndex = sharedByLink[link.linkIndex];
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            if (sharedIndex == missingIndex
                || minusTraceByLink[link.linkIndex] != missingIndex
                || plusTraceByLink[link.linkIndex] != missingIndex) {
                throw std::logic_error(
                    "opening velocity-metric shared base mapping is incomplete");
            }
            const auto& base = baseMetric.dofs[sharedIndex];
            append(
                PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid,
                link, base, base, 0, 0, link.areaSquareMeters, true);
            continue;
        }
        const std::size_t minusIndex = minusTraceByLink[link.linkIndex];
        const std::size_t plusIndex = plusTraceByLink[link.linkIndex];
        if (sharedIndex != missingIndex || minusIndex == missingIndex
            || plusIndex == missingIndex) {
            throw std::logic_error(
                "opening velocity-metric wall base mapping is incomplete");
        }
        const auto& minusBase = baseMetric.dofs[minusIndex];
        const auto& plusBase = baseMetric.dofs[plusIndex];
        result.totalPressureWallAreaSquareMeters += link.areaSquareMeters;
        double solidArea = link.areaSquareMeters;
        const std::size_t partitionIndex =
            partitionByLink[link.linkIndex];
        if (partitionIndex != missingIndex) {
            solidArea = openings.partitions[partitionIndex]
                .solidAreaSquareMeters;
        }
        result.totalSolidWallAreaSquareMeters += solidArea;
        if (solidArea > 0.0) {
            append(
                PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SolidWallMinusTrace,
                link, minusBase, plusBase, 0, 0, solidArea, false);
            append(
                PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SolidWallPlusTrace,
                link, minusBase, plusBase, 0, 0, solidArea, false);
        }
        while (patchOffset < patchOrder.size()
               && openings.patches[patchOrder[patchOffset]]
                      .sourceFaceLinkIndex < link.linkIndex) {
            throw std::logic_error(
                "opening velocity-metric patch ordering is incomplete");
        }
        while (patchOffset < patchOrder.size()
               && openings.patches[patchOrder[patchOffset]]
                      .sourceFaceLinkIndex == link.linkIndex) {
            const auto& patch = openings.patches[patchOrder[patchOffset]];
            if (patch.sourceFaceLinkStableId != link.stableId
                || patch.minusFragmentIndex != link.minusFragmentIndex
                || patch.plusFragmentIndex != link.plusFragmentIndex
                || patch.centerDistanceMeters != link.centerDistanceMeters) {
                throw std::logic_error(
                    "opening velocity-metric aperture geometry is inconsistent");
            }
            append(
                PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch,
                link, minusBase, plusBase, patch.patchIndex,
                patch.patchStableId, patch.areaSquareMeters, true);
            result.totalOpeningAreaSquareMeters += patch.areaSquareMeters;
            ++patchOffset;
        }
    }
    if (patchOffset != patchOrder.size()) {
        throw std::logic_error(
            "opening velocity-metric aperture coverage is incomplete");
    }
    result.wallAreaPartitionResidualSquareMeters =
        result.totalPressureWallAreaSquareMeters
        - (result.totalSolidWallAreaSquareMeters
           + result.totalOpeningAreaSquareMeters);

    for (const auto& dof : result.dofs) {
        ++result.fragments[dof.ownerFragmentIndex]
              .velocityDofIncidenceCount;
        ++result.components[dof.connectedComponentIndex]
              .velocityDofIncidenceCount;
        coordinate(
            result.fragments[dof.ownerFragmentIndex]
                .dualVolumeByAxisCubicMeters,
            dof.axis) += dof.ownerDualVolumeCubicMeters;
        coordinate(
            result.components[dof.connectedComponentIndex]
                .dualVolumeByAxisCubicMeters,
            dof.axis) += dof.ownerDualVolumeCubicMeters;
        if (dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid
            || dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch) {
            ++result.fragments[dof.oppositeFragmentIndex]
                  .velocityDofIncidenceCount;
            ++result.components[dof.connectedComponentIndex]
                  .velocityDofIncidenceCount;
            coordinate(
                result.fragments[dof.oppositeFragmentIndex]
                    .dualVolumeByAxisCubicMeters,
                dof.axis) += dof.oppositeDualVolumeCubicMeters;
            coordinate(
                result.components[dof.connectedComponentIndex]
                    .dualVolumeByAxisCubicMeters,
                dof.axis) += dof.oppositeDualVolumeCubicMeters;
        }
        coordinate(result.dualVolumeByAxisCubicMeters, dof.axis) +=
            dof.dualVolumeCubicMeters;
        result.totalDualVolumeCubicMeters += dof.dualVolumeCubicMeters;
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            ++result.sharedRegionGridDofCount;
            result.sharedRegionGridDualVolumeCubicMeters +=
                dof.dualVolumeCubicMeters;
        } else if (dof.kind
                   == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                       OpeningPatch) {
            ++result.openingPatchDofCount;
            result.openingPatchDualVolumeCubicMeters +=
                dof.dualVolumeCubicMeters;
        } else {
            ++result.solidWallTraceDofCount;
            result.solidWallTraceDualVolumeCubicMeters +=
                dof.dualVolumeCubicMeters;
        }
    }

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
        if (residual > closureTolerance(fragment.sourceVolumeCubicMeters)) {
            throw std::invalid_argument(
                "opening velocity-metric fragment volume closure failed");
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
                "opening velocity-metric component volume closure failed");
        }
    }
    const double domainVolume = grid.cellVolumeCubicMeters()
        * static_cast<double>(grid.cellCount());
    result.domainVolumeClosureResidualByAxisCubicMeters = {
        result.dualVolumeByAxisCubicMeters.x - domainVolume,
        result.dualVolumeByAxisCubicMeters.y - domainVolume,
        result.dualVolumeByAxisCubicMeters.z - domainVolume,
    };
    const double areaTolerance = closureTolerance(
        result.totalPressureWallAreaSquareMeters);
    if (!finiteVector(result.dualVolumeByAxisCubicMeters)
        || !finiteVector(result.domainVolumeClosureResidualByAxisCubicMeters)
        || maximumAbsoluteComponent(
               result.domainVolumeClosureResidualByAxisCubicMeters)
            > closureTolerance(domainVolume)
        || std::abs(result.wallAreaPartitionResidualSquareMeters)
            > areaTolerance
        || result.totalOpeningAreaSquareMeters
            != openings.totalOpeningAreaSquareMeters) {
        throw std::invalid_argument(
            "opening velocity-metric domain closure failed");
    }
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.dofs.size() != dofCount
        || result.ownedStorageBytes != maximumOwnedBytes) {
        throw std::logic_error(
            "opening velocity-metric preflight storage changed");
    }
    result.fingerprint = metricFingerprint(result);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentOpeningVelocityMetric
buildPlanarPressureRegionFragmentOpeningVelocityMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits) {
    return buildMetric(
        grid, sweep, fragments, topology, baseMetric,
        openingDefinitions, openings, limits);
}

void validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric) {
    if (metric.version
            != planarPressureRegionFragmentOpeningVelocityMetricVersion
        || metric.fingerprint == 0
        || metric.sourceBaseMetricFingerprint == 0
        || metric.sourceOpeningFingerprint == 0
        || metric.sourceFragmentFingerprint == 0
        || metric.sourceTopologyFingerprint == 0
        || metric.dofs.empty() || metric.fragments.empty()
        || metric.components.empty()
        || !finiteVector(metric.dualVolumeByAxisCubicMeters)
        || !finiteVector(metric.domainVolumeClosureResidualByAxisCubicMeters)
        || !std::isfinite(metric.totalDualVolumeCubicMeters)
        || !(metric.totalDualVolumeCubicMeters > 0.0)
        || !std::isfinite(metric.totalPressureWallAreaSquareMeters)
        || !(metric.totalPressureWallAreaSquareMeters > 0.0)
        || !std::isfinite(metric.totalSolidWallAreaSquareMeters)
        || metric.totalSolidWallAreaSquareMeters < 0.0
        || !std::isfinite(metric.totalOpeningAreaSquareMeters)
        || metric.totalOpeningAreaSquareMeters < 0.0
        || !std::isfinite(metric.wallAreaPartitionResidualSquareMeters)
        || !std::isfinite(
            metric.maximumAbsoluteFragmentVolumeClosureResidualCubicMeters)
        || !std::isfinite(
            metric.maximumAbsoluteComponentVolumeClosureResidualCubicMeters)
        || metric.sharedRegionGridDofCount
                + metric.solidWallTraceDofCount
                + metric.openingPatchDofCount
            != metric.dofs.size()
        || metric.ownedStorageBytes != ownedStorageBytes(metric)
        || metric.workingStorageBytes == 0
        || metric.fingerprint != metricFingerprint(metric)) {
        throw std::invalid_argument(
            "opening velocity-metric integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentOpeningVelocityMetric(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    const std::span<
        const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(metric);
    if (metric.dofs.size() > limits.maximumDofs
        || metric.fragments.size() > limits.maximumFragments
        || metric.components.size() > limits.maximumComponents
        || metric.ownedStorageBytes > limits.maximumOwnedBytes
        || metric.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "opening velocity-metric validation limit exceeded");
    }
    if (metric != buildMetric(
            grid, sweep, fragments, topology, baseMetric,
            openingDefinitions, openings, limits)) {
        throw std::invalid_argument(
            "opening velocity metric is foreign to its source");
    }
}

} // namespace simwing::fsi::fluid
