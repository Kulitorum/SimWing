#pragma once

#include "fluid/moving_interface.h"
#include "transfer.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t uniformFluidStructureBridgeVersion = 1;

struct UniformFluidStructureBridgeSettings {
    ConservativeTransferSettings transfer;
    double maximumPressureTractionDeviationPascals = 1.0e-10;
    double absoluteAreaToleranceSquareMeters = 1.0e-12;
    double relativeAreaTolerance = 1.0e-12;
    double absoluteForceToleranceNewtons = 1.0e-10;
    double relativeForceTolerance = 1.0e-12;
    double absolutePowerToleranceWatts = 1.0e-10;
    double relativePowerTolerance = 1.0e-12;
};

// Independent source and target ledgers. The bridge accepts only a pressure
// field that is uniform on its face-aligned fluid surface, because the current
// moving-interface diagnostic does not retain face-to-triangle geometry.
struct UniformFluidStructureBridgeDiagnostics {
    std::uint32_t version = uniformFluidStructureBridgeVersion;
    std::uint64_t fluidSurfaceStableId = 0;
    std::size_t fluidFaceCount = 0;
    std::size_t structureTriangleCount = 0;
    double fluidAreaSquareMeters = 0.0;
    double structureAreaSquareMeters = 0.0;
    double areaResidualSquareMeters = 0.0;
    double maximumPressureTractionDeviationPascals = 0.0;
    StructureVector3 uniformPressureTractionPascals;
    StructureVector3 fluidPressureForceNewtons;
    StructureVector3 structureSurfaceForceNewtons;
    StructureVector3 forceResidualNewtons;
    double forceResidualNormNewtons = 0.0;
    double fluidPressurePowerWatts = 0.0;
    double structureSurfacePowerWatts = 0.0;
    double powerResidualWatts = 0.0;
    bool finite = true;

    bool operator==(
        const UniformFluidStructureBridgeDiagnostics&) const = default;
};

class UniformFluidStructureTransferResult final {
public:
    [[nodiscard]] const ConservativeTransferResult&
    transferResult() const noexcept;
    [[nodiscard]] const UniformFluidStructureBridgeDiagnostics&
    diagnostics() const noexcept;

    bool operator==(
        const UniformFluidStructureTransferResult&) const = default;

private:
    friend class UniformFluidStructureBridge;
    UniformFluidStructureTransferResult(
        ConservativeTransferResult transferResult,
        UniformFluidStructureBridgeDiagnostics diagnostics);

    ConservativeTransferResult transferResult_;
    UniformFluidStructureBridgeDiagnostics diagnostics_;
};

// Stable-ID bridge for the deliberately narrow verification case: one fluid
// surface maps to one structural coupling surface with uniform world-space
// pressure traction. Any nonuniform reconstruction, area mismatch, force
// mismatch, or interface-power mismatch is rejected before structural state
// can be mutated.
class UniformFluidStructureBridge final {
public:
    UniformFluidStructureBridge(
        const Structure& target,
        std::uint64_t fluidSurfaceStableId,
        std::vector<CouplingSurfaceNodeDefinition> nodes,
        std::vector<CouplingSurfaceTriangleDefinition> triangles);

    [[nodiscard]] std::uint64_t fluidSurfaceStableId() const noexcept;
    [[nodiscard]] const ConservativeSurfaceTransfer& transfer() const noexcept;

    [[nodiscard]] UniformFluidStructureTransferResult evaluate(
        const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
        std::span<const CouplingNodeKinematics> nodeKinematics,
        const UniformFluidStructureBridgeSettings& settings = {}) const;

private:
    std::uint64_t fluidSurfaceStableId_ = 0;
    ConservativeSurfaceTransfer transfer_;
};

} // namespace simwing::fsi
