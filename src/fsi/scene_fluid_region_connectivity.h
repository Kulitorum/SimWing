#pragma once

#include "scene_fluid_region_continuity.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionConnectivityVersion = 1;
inline constexpr std::uint32_t
    sceneFluidRegionComponentContinuityVersion = 1;

struct SceneFluidRegionConnectivityLimits {
    std::size_t maximumRegions = 1'000'000;
    std::size_t maximumOpenings = 1'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumConnectivityBytes = 256ULL * 1024ULL * 1024ULL;
};

struct SceneFluidConnectedRegion {
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    RegionKind kind = RegionKind::Cell;
    std::size_t componentIndex = 0;

    bool operator==(const SceneFluidConnectedRegion&) const = default;
};

struct SceneFluidConnectedOpening {
    std::size_t openingIndex = 0;
    StableId openingId = invalidStableId;
    OpeningRole role = OpeningRole::Intake;
    std::size_t negativeSideRegionIndex = 0;
    std::size_t positiveSideRegionIndex = 0;
    std::size_t componentIndex = 0;

    bool operator==(const SceneFluidConnectedOpening&) const = default;
};

struct SceneFluidRegionComponent {
    std::size_t componentIndex = 0;
    StableId gaugeRegionId = invalidStableId;
    bool containsOutside = false;
    std::size_t firstRegionMember = 0;
    std::size_t regionCount = 0;
    std::size_t firstOpeningMember = 0;
    std::size_t openingCount = 0;

    bool operator==(const SceneFluidRegionComponent&) const = default;
};

// Authored openings are undirected pressure-connectivity edges even though
// their transport sign remains negative-side to positive-side. Components
// and members are canonicalized by stable ID; the smallest region ID owns the
// component's future pressure gauge. This product changes no grid topology.
struct SceneFluidRegionConnectivity {
    std::uint32_t version = sceneFluidRegionConnectivityVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::size_t ownedStorageBytes = 0;
    std::vector<SceneFluidConnectedRegion> regions;
    std::vector<SceneFluidConnectedOpening> openings;
    std::vector<SceneFluidRegionComponent> components;
    std::vector<std::size_t> componentRegionIndices;
    std::vector<std::size_t> componentOpeningIndices;

    bool operator==(const SceneFluidRegionConnectivity&) const = default;
};

[[nodiscard]] SceneFluidRegionConnectivity
buildSceneFluidRegionConnectivity(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidRegionConnectivityLimits& limits = {});

void validateSceneFluidRegionConnectivity(
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidSurfaceDefinition& surface);

struct SceneFluidRegionComponentContinuity {
    std::size_t componentIndex = 0;
    StableId gaugeRegionId = invalidStableId;
    bool containsOutside = false;
    std::size_t regionCount = 0;
    std::size_t openingCount = 0;
    double previousOpeningSourceResidualCubicMetersPerSecond = 0.0;
    double currentOpeningSourceResidualCubicMetersPerSecond = 0.0;
    double openingSourceToleranceCubicMetersPerSecond = 0.0;
    double geometryVolumeChangeCubicMeters = 0.0;
    double integratedOutwardRelativeVolumeCubicMeters = 0.0;
    double continuityResidualCubicMeters = 0.0;
    double toleranceCubicMeters = 0.0;
    bool compatible = false;

    bool operator==(
        const SceneFluidRegionComponentContinuity&) const = default;
};

struct SceneFluidRegionComponentContinuitySet {
    std::uint32_t version =
        sceneFluidRegionComponentContinuityVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t surfaceDefinitionFingerprint = 0;
    std::uint64_t regionConnectivityFingerprint = 0;
    std::uint64_t regionContinuityFingerprint = 0;
    std::uint64_t previousOpeningFluxFingerprint = 0;
    std::uint64_t currentOpeningFluxFingerprint = 0;
    std::uint64_t previousAcceptedStepCount = 0;
    std::uint64_t currentAcceptedStepCount = 0;
    double durationSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t failedComponentCount = 0;
    double maximumAbsoluteOpeningSourceResidualCubicMetersPerSecond = 0.0;
    double maximumAbsoluteContinuityResidualCubicMeters = 0.0;
    double globalContinuityResidualCubicMeters = 0.0;
    bool allComponentsCompatible = false;
    std::vector<SceneFluidRegionComponentContinuity> components;

    bool operator==(
        const SceneFluidRegionComponentContinuitySet&) const = default;
};

[[nodiscard]] SceneFluidRegionComponentContinuitySet
auditSceneFluidRegionComponentContinuity(
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux,
    const SceneFluidRegionConnectivityLimits& limits = {});

void validateSceneFluidRegionComponentContinuity(
    const SceneFluidRegionComponentContinuitySet& components,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidRegionContinuitySet& continuity,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidCellVolumeSet& previousVolumes,
    const SceneFluidCellVolumeSet& currentVolumes,
    const SceneFluidOpeningFluxSet& previousFlux,
    const SceneFluidOpeningFluxSet& currentFlux);

} // namespace simwing::fsi
