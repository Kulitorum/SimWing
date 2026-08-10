#pragma once

#include "scene_fluid_mimetic_pressure_audit.h"
#include "scene_fluid_opening_patch.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticCorrectedTraceFlowVersion = 1;
inline constexpr std::uint32_t
    sceneFluidMimeticMacVelocityCollapseVersion = 1;

struct SceneFluidMimeticCorrectedTrace {
    std::size_t sharedTraceOrdinal = 0;
    std::size_t traceIndex = 0;
    std::size_t sourceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::uint64_t sourceStableId = 0;
    std::size_t componentIndex = 0;
    std::size_t minusControlCellIndex = 0;
    std::size_t plusControlCellIndex = 0;
    double predictedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;
    double pressureCorrectionVolumeFlowRateCubicMetersPerSecond = 0.0;
    double correctedRelativeVolumeFlowRateCubicMetersPerSecond = 0.0;

    bool operator==(
        const SceneFluidMimeticCorrectedTrace&) const = default;
};

// Accepted pressure is reconstructed through the same condensed/full-system
// boundary used by the atomic solve. The pressure-induced outward half-face
// flux is converted from Pa*m to m^3/s by dt/rho and added to the corresponding
// predicted shared flow. Material-wall traces remain exactly local and are not
// published as transport degrees of freedom.
struct SceneFluidMimeticCorrectedTraceFlow {
    std::uint32_t version =
        sceneFluidMimeticCorrectedTraceFlowVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t pressureAuditFingerprint = 0;
    std::uint64_t mimeticControlCellFingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t condensedTraceSystemFingerprint = 0;
    std::uint64_t predictedTraceFlowFingerprint = 0;
    std::uint64_t pressureSourceFingerprint = 0;
    std::uint64_t pressureStateFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double densityKgPerCubicMeter = 0.0;
    double timeStepSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t cartesianTraceCount = 0;
    std::size_t authoredOpeningTraceCount = 0;
    double maximumAbsolutePressureCorrectionCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectedFlowCubicMetersPerSecond = 0.0;
    double maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond =
        0.0;
    double correctedContinuityToleranceCubicMetersPerSecond = 0.0;
    bool finite = false;
    bool accepted = false;
    std::vector<SceneFluidMimeticCorrectedTrace> traces;
    std::vector<double>
        correctedContinuityResidualsCubicMetersPerSecond;

    bool operator==(
        const SceneFluidMimeticCorrectedTraceFlow&) const = default;
};

[[nodiscard]] SceneFluidMimeticCorrectedTraceFlow
correctSceneFluidMimeticTraceFlows(
    const SceneFluidMimeticPressureAuditEndpoint& endpoint);

void validateSceneFluidMimeticCorrectedTraceFlowIntegrity(
    const SceneFluidMimeticCorrectedTraceFlow& flow);

struct SceneFluidMimeticMacVelocityCollapseDiagnostics {
    std::uint32_t version =
        sceneFluidMimeticMacVelocityCollapseVersion;
    std::uint64_t correctedTraceFlowFingerprint = 0;
    std::uint64_t pressureFaceLinkFingerprint = 0;
    std::uint64_t openingPatchFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t faceCount = 0;
    std::size_t cartesianTraceCount = 0;
    std::size_t authoredOpeningTraceCount = 0;
    std::size_t embeddedOpeningTraceCount = 0;
    std::size_t multiTraceFaceCount = 0;
    double maximumAbsoluteVelocityMetersPerSecond = 0.0;
    double maximumSubfaceVelocityDeviationMetersPerSecond = 0.0;
    double maximumVolumeFlowClosureCubicMetersPerSecond = 0.0;
    bool finite = false;

    bool operator==(
        const SceneFluidMimeticMacVelocityCollapseDiagnostics&) const =
        default;
};

struct SceneFluidMimeticMacVelocityCollapse {
    explicit SceneFluidMimeticMacVelocityCollapse(
        const fluid::PeriodicCartesianGrid& grid)
        : velocityMetersPerSecond(grid) {}

    fluid::MacVelocityField velocityMetersPerSecond;
    SceneFluidMimeticMacVelocityCollapseDiagnostics diagnostics;

    bool operator==(
        const SceneFluidMimeticMacVelocityCollapse&) const = default;
};

// Area-collapses corrected Cartesian subface traces onto one bulk normal MAC
// velocity. Face-owned opening traces restore accepted cap sweep before the
// collapse. Cell-owned opening traces have no unique Cartesian face and remain
// explicit in diagnostics instead of being smeared into the bulk field.
[[nodiscard]] SceneFluidMimeticMacVelocityCollapse
collapseSceneFluidMimeticCorrectedMacVelocity(
    const SceneFluidMimeticCorrectedTraceFlow& correctedFlow,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const fluid::PeriodicCartesianGrid& grid);

} // namespace simwing::fsi
