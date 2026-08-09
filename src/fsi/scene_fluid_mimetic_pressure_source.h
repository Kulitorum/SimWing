#pragma once

#include "scene_fluid_mimetic_control_cell.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t sceneFluidMimeticPressureSourceVersion = 1;

struct SceneFluidMimeticPressureSourceSettings {
    double densityKgPerCubicMeter = 1.225;
    double timeStepSeconds = 1.0 / 60.0;

    bool operator==(
        const SceneFluidMimeticPressureSourceSettings&) const = default;
};

struct SceneFluidMimeticPressureSourceLimits {
    std::size_t maximumControlCells = 50'000'000;
    std::size_t maximumComponents = 1'000'000;
    std::size_t maximumOwnedBytes = 2048ULL * 1024ULL * 1024ULL;
};

struct SceneFluidMimeticPressureControlSource {
    std::size_t controlCellIndex = 0;
    std::size_t controlVolumeIndex = 0;
    std::uint64_t stableId = 0;
    std::size_t componentIndex = 0;
    double geometryVolumeChangeRateCubicMetersPerSecond = 0.0;
    double predictedNetOutwardVolumeFlowRateCubicMetersPerSecond = 0.0;
    double predictedContinuityResidualCubicMetersPerSecond = 0.0;
    double integratedSourcePascalsMeters = 0.0;

    bool operator==(
        const SceneFluidMimeticPressureControlSource&) const = default;
};

// Immutable physical-unit bridge into the mimetic trace solve. Each control
// receives -(rho/dt)*(dV/dt + predicted net outward volume flow), matching the
// production projection's integrated Pa*m sign convention. This product does
// not sample link flows or choose moving-volume rates; those remain explicit
// upstream owners.
struct SceneFluidMimeticPressureSourceSet {
    std::uint32_t version = sceneFluidMimeticPressureSourceVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    SceneFluidMimeticPressureSourceSettings settings;
    std::size_t ownedStorageBytes = 0;
    std::size_t componentCount = 0;
    double maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond = 0.0;
    double maximumAbsolutePredictedNetOutwardVolumeRateCubicMetersPerSecond =
        0.0;
    double maximumAbsolutePredictedContinuityResidualCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteComponentIntegratedSourcePascalsMeters = 0.0;
    std::vector<SceneFluidMimeticPressureControlSource> controls;
    std::vector<double> componentContinuityResidualsCubicMetersPerSecond;
    std::vector<double> componentIntegratedSourcesPascalsMeters;

    bool operator==(
        const SceneFluidMimeticPressureSourceSet&) const = default;
};

[[nodiscard]] SceneFluidMimeticPressureSourceSet
buildSceneFluidMimeticPressureSources(
    const SceneFluidMimeticControlCellSet& controlCells,
    std::span<const double>
        predictedNetOutwardVolumeFlowRatesCubicMetersPerSecond,
    const SceneFluidMimeticPressureSourceSettings& settings = {},
    const SceneFluidMimeticPressureSourceLimits& limits = {});

[[nodiscard]] SceneFluidMimeticPressureSourceSet
buildSceneFluidMimeticPressureSources(
    const SceneFluidMimeticControlCellSet& controlCells,
    std::span<const double>
        predictedNetOutwardVolumeFlowRatesCubicMetersPerSecond,
    std::span<const double>
        geometryVolumeChangeRatesCubicMetersPerSecond,
    const SceneFluidMimeticPressureSourceSettings& settings = {},
    const SceneFluidMimeticPressureSourceLimits& limits = {});

void validateSceneFluidMimeticPressureSourceIntegrity(
    const SceneFluidMimeticPressureSourceSet& sources);

void validateSceneFluidMimeticPressureSources(
    const SceneFluidMimeticPressureSourceSet& sources,
    const SceneFluidMimeticControlCellSet& controlCells);

[[nodiscard]] std::vector<double>
sceneFluidMimeticIntegratedCellSources(
    const SceneFluidMimeticPressureSourceSet& sources);

} // namespace simwing::fsi
