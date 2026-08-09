#pragma once

#include "fluid/mimetic_wall_condensation.h"
#include "scene_fluid_mimetic_trace_system.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t
    sceneFluidMimeticCondensedTraceSystemVersion = 1;
inline constexpr std::size_t invalidSceneFluidMimeticCondensedTraceIndex =
    std::numeric_limits<std::size_t>::max();

struct SceneFluidMimeticCondensedTraceSystemLimits {
    std::size_t maximumReducedTraces = 200'000'000;
    std::size_t maximumLocalCondensations = 50'000'000;
    std::size_t maximumLocalCondensationBytes =
        4096ULL * 1024ULL * 1024ULL;
    std::size_t maximumOwnedBytes = 8192ULL * 1024ULL * 1024ULL;
    fluid::MimeticWallCondensationSettings local;
};

struct SceneFluidMimeticCondensedTrace {
    std::size_t traceIndex = 0;
    std::size_t fullTraceIndex = 0;
    std::uint64_t stableId = 0;
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::size_t componentIndex = 0;
    double operatorDiagonal = 0.0;
    bool isGauge = false;

    bool operator==(
        const SceneFluidMimeticCondensedTrace&) const = default;
};

// Immutable global Schur topology over only shared Cartesian and
// authored-opening traces. Material-wall traces remain addressable in the
// bound full system for exact RHS condensation and post-solve reconstruction,
// but are absent from this reduced field. Every component owns a deterministic
// shared gauge. This audit product does not replace the production graph path.
struct SceneFluidMimeticCondensedTraceSystem {
    std::uint32_t version =
        sceneFluidMimeticCondensedTraceSystemVersion;
    std::uint64_t fingerprint = 0;
    std::uint64_t fullTraceSystemFingerprint = 0;
    std::uint64_t structureDefinitionFingerprint = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::size_t ownedStorageBytes = 0;
    std::size_t localCondensationStorageBytes = 0;
    std::size_t componentCount = 0;
    std::size_t eliminatedMaterialWallTraceCount = 0;
    double minimumPositiveOperatorDiagonal = 0.0;
    double maximumOperatorDiagonal = 0.0;
    std::vector<SceneFluidMimeticCondensedTrace> traces;
    std::vector<std::size_t> fullTraceReducedIndices;
    std::vector<std::size_t> componentGaugeTraceIndices;
    std::vector<fluid::MimeticWallCondensation> localCondensations;

    bool operator==(
        const SceneFluidMimeticCondensedTraceSystem&) const = default;
};

[[nodiscard]] SceneFluidMimeticCondensedTraceSystem
buildSceneFluidMimeticCondensedTraceSystem(
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystemLimits& limits = {});

void validateSceneFluidMimeticCondensedTraceSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem);

[[nodiscard]] std::vector<double>
applySceneFluidMimeticCondensedTraceOperator(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    std::span<const double> reducedTraceScalars);

// Condenses a full-system trace RHS to the shared trace field.
[[nodiscard]] std::vector<double>
condenseSceneFluidMimeticTraceRightHandSide(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    std::span<const double> fullTraceRightHandSide);

// Reconstructs the full trace field after the reduced system has been solved.
[[nodiscard]] std::vector<double> reconstructSceneFluidMimeticFullTraces(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    std::span<const double> fullTraceRightHandSide,
    std::span<const double> reducedTraceScalars);

} // namespace simwing::fsi
