#pragma once

#include "scene_fluid_surface.h"
#include "transfer.h"

#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

// Binds one authoritative scene fluid surface to the existing conservative
// structural transfer topology. Region/material/opening ownership remains in
// SceneFluidSurfaceDefinition; this adapter neither classifies a grid nor
// invents traction.
class SceneFluidSurfaceTransfer final {
public:
    SceneFluidSurfaceTransfer(
        const SceneFluidSurfaceDefinition& surface,
        const SceneStructureMappings& structureMappings,
        const Structure& target);

    [[nodiscard]] std::uint64_t surfaceDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t couplingSurfaceFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t targetDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::span<const CouplingSurfaceNodeDefinition>
    nodes() const noexcept;
    [[nodiscard]] std::span<const CouplingSurfaceTriangleDefinition>
    triangles() const noexcept;

    [[nodiscard]] std::vector<CouplingNodeKinematics> kinematics(
        const SceneFluidSurfaceState& state) const;

    [[nodiscard]] ConservativeTransferResult evaluate(
        const SceneFluidSurfaceState& state,
        std::span<const CouplingTriangleTraction> triangleTractions,
        const ConservativeTransferSettings& settings = {}) const;

    [[nodiscard]] ConservativeTransferResult evaluateQuadrature(
        const SceneFluidSurfaceState& state,
        std::span<const CouplingTriangleTractionQuadrature> quadrature,
        const ConservativeTransferSettings& settings = {}) const;

    void addLoadsTo(Structure& target,
                    const ConservativeTransferResult& result) const;

private:
    struct Topology {
        std::vector<CouplingSurfaceNodeDefinition> nodes;
        std::vector<CouplingSurfaceTriangleDefinition> triangles;
    };

    SceneFluidSurfaceTransfer(std::uint64_t surfaceDefinitionFingerprint,
                              const Structure& target,
                              Topology topology);
    static Topology makeTopology(
        const SceneFluidSurfaceDefinition& surface,
        const SceneStructureMappings& structureMappings,
        const Structure& target);

    std::uint64_t surfaceDefinitionFingerprint_ = 0;
    ConservativeSurfaceTransfer transfer_;
};

} // namespace simwing::fsi
