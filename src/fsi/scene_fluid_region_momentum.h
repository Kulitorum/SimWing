#pragma once

#include "scene_fluid_pressure_projection.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidRegionMomentumVersion = 2;

struct SceneFluidRegionMomentumLimits {
    std::size_t maximumControlVolumes = 50'000'000;
    std::size_t maximumMomentumBytes =
        4ULL * 1024ULL * 1024ULL * 1024ULL;
};

struct SceneFluidRegionMomentumControlVolume {
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t cellIndex = 0;
    std::size_t regionIndex = 0;
    StableId regionId = invalidStableId;
    std::size_t componentIndex = 0;
    double volumeCubicMeters = 0.0;
    fluid::Vector3 velocityMetersPerSecond;
    fluid::Vector3 momentumKilogramMetersPerSecond;
    std::array<double, 3> sampledFaceAreaSquareMeters{};
    std::array<std::size_t, 3> sampledLinkCounts{};

    bool operator==(
        const SceneFluidRegionMomentumControlVolume&) const = default;
};

struct SceneFluidRegionMomentumDiagnostics {
    std::size_t controlVolumeCount = 0;
    std::size_t linkCount = 0;
    std::size_t openingLinkCount = 0;
    std::size_t embeddedOpeningLinkCount = 0;
    std::size_t normalEquationControlCount = 0;
    std::size_t sampledComponentCount = 0;
    std::size_t fallbackComponentCount = 0;
    fluid::Vector3 totalMomentumKilogramMetersPerSecond;
    double kineticEnergyJoules = 0.0;
    double maximumAbsoluteVelocityMetersPerSecond = 0.0;
    double maximumLinkNormalVelocityResidualMetersPerSecond = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidRegionMomentumDiagnostics&) const = default;
};

// Accepted pressure-link flow reconstructed as one collocated momentum vector
// per positive cell/region control volume. Controls incident only to Cartesian
// links retain the exact component-wise area average. A control incident to an
// embedded link instead uses a bounded three-dimensional normal-equation solve
// over all its link normals, retaining the cell-centred MAC predictor in the
// unconstrained nullspace. The sampled component arrays continue to describe
// Cartesian-face coverage only. This is immutable momentum ownership for later
// conservative transport; it does not itself advance, project, or impose a
// wall model.
struct SceneFluidRegionMomentumState {
    std::uint32_t version = sceneFluidRegionMomentumVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t pressureProjectionFingerprint = 0;
    std::uint64_t pressureControlVolumeFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t fallbackVelocityFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t ownedStorageBytes = 0;
    SceneFluidRegionMomentumDiagnostics diagnostics;
    std::vector<SceneFluidRegionMomentumControlVolume> controlVolumes;

    bool operator==(const SceneFluidRegionMomentumState&) const = default;
};

[[nodiscard]] SceneFluidRegionMomentumState
reconstructSceneFluidRegionMomentumState(
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond,
    const SceneFluidRegionMomentumLimits& limits = {});

void validateSceneFluidRegionMomentumStateIntegrity(
    const SceneFluidRegionMomentumState& momentum);

// Validates the persisted state against its accepted pressure epoch without
// requiring the transient MAC fallback field that was used to construct it.
// The fallback field remains cryptographically bound by its stored
// fingerprint and the state's own integrity fingerprint.
void validateSceneFluidRegionMomentumStateBinding(
    const SceneFluidRegionMomentumState& momentum,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection);

void validateSceneFluidRegionMomentumState(
    const SceneFluidRegionMomentumState& momentum,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidPressureProjection& projection,
    const fluid::MacVelocityField& fallbackVelocityMetersPerSecond);

} // namespace simwing::fsi
