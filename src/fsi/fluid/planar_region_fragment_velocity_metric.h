#pragma once

#include "fluid/planar_region_fragment_topology.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t
    planarPressureRegionFragmentVelocityMetricVersion = 1;

enum class PlanarPressureRegionFragmentVelocityDofKind : std::uint8_t {
    SharedRegionGrid = 1,
    PressureLayerMinusTrace = 2,
    PressureLayerPlusTrace = 3,
};

struct PlanarPressureRegionFragmentVelocityDof {
    std::size_t dofIndex = 0;
    std::uint64_t stableId = 0;
    PlanarPressureRegionFragmentVelocityDofKind kind =
        PlanarPressureRegionFragmentVelocityDofKind::SharedRegionGrid;
    std::size_t sourceFaceLinkIndex = 0;
    std::uint64_t sourceFaceLinkStableId = 0;
    GridFaceAxis axis = GridFaceAxis::X;
    std::uint64_t surfaceStableId = 0;
    std::size_t ownerFragmentIndex = 0;
    std::uint64_t ownerFragmentStableId = 0;
    std::size_t oppositeFragmentIndex = 0;
    std::uint64_t oppositeFragmentStableId = 0;
    std::size_t componentIndex = 0;
    std::uint64_t regionStableId = 0;
    double areaSquareMeters = 0.0;
    double ownerHalfDistanceMeters = 0.0;
    double oppositeHalfDistanceMeters = 0.0;
    double ownerDualVolumeCubicMeters = 0.0;
    double oppositeDualVolumeCubicMeters = 0.0;
    double dualVolumeCubicMeters = 0.0;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityDof&) const = default;
};

struct PlanarPressureRegionFragmentVelocityMetricFragment {
    std::size_t fragmentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t velocityDofIncidenceCount = 0;
    Vector3 dualVolumeByAxisCubicMeters;
    double sourceVolumeCubicMeters = 0.0;
    Vector3 volumeClosureResidualByAxisCubicMeters;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityMetricFragment&) const =
        default;
};

struct PlanarPressureRegionFragmentVelocityMetricComponent {
    std::size_t componentIndex = 0;
    std::uint64_t stableId = 0;
    std::uint64_t regionStableId = 0;
    std::size_t sharedGridDofCount = 0;
    std::size_t pressureLayerTraceDofCount = 0;
    Vector3 dualVolumeByAxisCubicMeters;
    double sourceVolumeCubicMeters = 0.0;
    Vector3 volumeClosureResidualByAxisCubicMeters;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityMetricComponent&) const =
        default;
};

struct PlanarPressureRegionFragmentVelocityMetricLimits {
    PlanarPressureRegionFragmentTopologyLimits topologyLimits;
    std::size_t maximumDofs = 120'000'000;
    std::size_t maximumFragments = 20'000'000;
    std::size_t maximumComponents = 20'000'000;
    std::size_t maximumOwnedBytes = 4096ULL * 1024ULL * 1024ULL;
};

// Diagonal orthogonal face metric for the regional fragment graph. A
// same-region Cartesian face owns one shared normal-velocity DOF with dual
// volume area*(minus half distance + plus half distance). A pressure-layer
// wall instead owns two independent one-sided trace DOFs, each with only its
// adjacent area*half-distance volume. No velocity or momentum is therefore
// averaged across fabric.
//
// The dual-volume incidences close each fragment and pressure component to its
// physical volume independently on X/Y/Z. This immutable product defines
// inertia geometry only; it owns no velocity values, density, kinetic-energy
// acceptance, wall prescription, advection, or production state.
struct PlanarPressureRegionFragmentVelocityMetric {
    std::uint32_t version =
        planarPressureRegionFragmentVelocityMetricVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t sourceFragmentFingerprint = 0;
    std::uint64_t sourceTopologyFingerprint = 0;
    GridFaceAxis profileAxis = GridFaceAxis::X;
    std::vector<PlanarPressureRegionFragmentVelocityDof> dofs;
    std::vector<PlanarPressureRegionFragmentVelocityMetricFragment> fragments;
    std::vector<PlanarPressureRegionFragmentVelocityMetricComponent> components;
    std::size_t sharedRegionGridDofCount = 0;
    std::size_t pressureLayerTraceDofCount = 0;
    double sharedRegionGridDualVolumeCubicMeters = 0.0;
    double pressureLayerTraceDualVolumeCubicMeters = 0.0;
    double totalDualVolumeCubicMeters = 0.0;
    Vector3 dualVolumeByAxisCubicMeters;
    Vector3 domainVolumeClosureResidualByAxisCubicMeters;
    double maximumAbsoluteFragmentVolumeClosureResidualCubicMeters = 0.0;
    double maximumAbsoluteComponentVolumeClosureResidualCubicMeters = 0.0;
    std::size_t ownedStorageBytes = 0;

    bool operator==(
        const PlanarPressureRegionFragmentVelocityMetric&) const = default;
};

[[nodiscard]] PlanarPressureRegionFragmentVelocityMetric
buildPlanarPressureRegionFragmentVelocityMetric(
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits = {});

void validatePlanarPressureRegionFragmentVelocityMetric(
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const PeriodicCartesianGrid& grid,
    const PlanarPressureRegionSweepLedger& sweep,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentTopology& topology,
    const PlanarPressureRegionFragmentVelocityMetricLimits& limits = {});

} // namespace simwing::fsi::fluid
