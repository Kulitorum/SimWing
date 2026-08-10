#pragma once

#include "fluid/planar_region_fragment_opening.h"
#include "fluid/planar_region_fragment_velocity_metric.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentOpeningVelocityMetricVersion = 1;

enum class PlanarPressureRegionFragmentOpeningVelocityDofKind
    : std::uint8_t {
    SharedRegionGrid = 1,
    SolidWallMinusTrace = 2,
    SolidWallPlusTrace = 3,
    OpeningPatch = 4,
};

struct PlanarPressureRegionFragmentOpeningVelocityDof {
    std::size_t dofIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentOpeningVelocityDofKind kind =
        PlanarPressureRegionFragmentOpeningVelocityDofKind::
            SharedRegionGrid;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    std::size_t sourceOpeningPatchIndex = 0;
    std::uint64_t sourceOpeningPatchStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t surfaceStableId = 0;
    std::size_t ownerFragmentIndex = 0;
    std::uint64_t ownerFragmentStableId = 0;
    std::size_t oppositeFragmentIndex = 0;
    std::uint64_t oppositeFragmentStableId = 0;
    std::size_t ownerBaseComponentIndex = 0;
    std::size_t oppositeBaseComponentIndex = 0;
    std::size_t connectedComponentIndex = 0;
    double areaSquareMeters = 0.0;
    double ownerHalfDistanceMeters = 0.0;
    double oppositeHalfDistanceMeters = 0.0;
    double ownerDualVolumeCubicMeters = 0.0;
    double oppositeDualVolumeCubicMeters = 0.0;
    double dualVolumeCubicMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityDof&) const =
        default;
};

struct PlanarPressureRegionFragmentOpeningVelocityMetricFragment {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t baseComponentIndex = 0;
    std::size_t connectedComponentIndex = 0;
    std::size_t velocityDofIncidenceCount = 0;
    Vector3 dualVolumeByAxisCubicMeters;
    double sourceVolumeCubicMeters = 0.0;
    Vector3 volumeClosureResidualByAxisCubicMeters;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityMetricFragment&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningVelocityMetricComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t baseComponentCount = 0;
    std::size_t fragmentCount = 0;
    std::size_t velocityDofIncidenceCount = 0;
    Vector3 dualVolumeByAxisCubicMeters;
    double sourceVolumeCubicMeters = 0.0;
    Vector3 volumeClosureResidualByAxisCubicMeters;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityMetricComponent&)
        const = default;
};

struct PlanarPressureRegionFragmentOpeningVelocityMetricLimits {
    PlanarPressureRegionFragmentVelocityMetricLimits baseMetricLimits;
    PlanarPressureRegionFragmentOpeningLimits openingLimits;
    std::size_t maximumDofs = 140'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 6144ULL * 1024ULL * 1024ULL;
    std::size_t maximumWorkingBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Opening-aware diagonal face geometry for regional momentum ownership. A
// same-region Cartesian face retains its shared dual volume. Every pressure-
// layer wall half-volume is instead partitioned by physical area: retained
// solid area owns an independent one-sided trace, while each aperture patch
// owns one shared negative-to-positive degree with distinct minus and plus
// half-volumes. Thus no inertia is duplicated or lost when a wall is partly
// or fully open.
//
// Fragment and opening-connected-component incidences close to their physical
// volume independently on X/Y/Z. Stable degree identity survives a topology-
// stable metric change. This immutable product owns geometry only: it assigns
// no wall/aperture velocity, density, momentum, transport, pressure step, or
// production-worker state.
struct PlanarPressureRegionFragmentOpeningVelocityMetric {
    std::uint32_t version =
        planarPressureRegionFragmentOpeningVelocityMetricVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceBaseMetricFingerprint = 0;
    std::uint64_t sourceOpeningFingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityDof> dofs;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityMetricFragment>
        fragments;
    std::vector<PlanarPressureRegionFragmentOpeningVelocityMetricComponent>
        components;
    std::size_t sharedRegionGridDofCount = 0;
    std::size_t solidWallTraceDofCount = 0;
    std::size_t openingPatchDofCount = 0;
    double sharedRegionGridDualVolumeCubicMeters = 0.0;
    double solidWallTraceDualVolumeCubicMeters = 0.0;
    double openingPatchDualVolumeCubicMeters = 0.0;
    double totalDualVolumeCubicMeters = 0.0;
    double totalPressureWallAreaSquareMeters = 0.0;
    double totalSolidWallAreaSquareMeters = 0.0;
    double totalOpeningAreaSquareMeters = 0.0;
    double wallAreaPartitionResidualSquareMeters = 0.0;
    Vector3 dualVolumeByAxisCubicMeters;
    Vector3 domainVolumeClosureResidualByAxisCubicMeters;
    double maximumAbsoluteFragmentVolumeClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteComponentVolumeClosureResidualCubicMeters = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t workingStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentOpeningVelocityMetric&) const =
        default;
};

[[nodiscard]] PlanarPressureRegionFragmentOpeningVelocityMetric
buildPlanarPressureRegionFragmentOpeningVelocityMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits =
        {});

void validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric);

void validatePlanarPressureRegionFragmentOpeningVelocityMetric(
    const PlanarPressureRegionFragmentOpeningVelocityMetric& metric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetric& baseMetric,
    std::span<const PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const PlanarPressureRegionFragmentOpeningSet& openings,
    const PlanarPressureRegionFragmentOpeningVelocityMetricLimits& limits =
        {});

} // namespace simwing::fsi::fluid
