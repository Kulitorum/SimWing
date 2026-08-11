#pragma once

#include "amr_external_flow.h"
#include "scene.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simwing::fsi::amr {

inline constexpr std::uint32_t amrStaticWingInterfaceVersion = 1;

struct StaticWingInterfaceCell {
    // Area-weighted average of n*n^T for all positive-area material patches
    // owned by this cell. Components are xx, xy, xz, yy, yz, zz.
    std::array<double, 6> normalProjectionTensor{};
    double surfaceAreaSquareMeters = 0.0;

    [[nodiscard]] bool active() const noexcept {
        return surfaceAreaSquareMeters > 0.0;
    }
};

struct StaticWingInterfaceLevel {
    fluid::GridCellCounts cellCounts;
    std::vector<StaticWingInterfaceCell> cells;
};

struct StaticWingInterfaceVertex {
    StableId id = invalidStableId;
    Vec3 positionMeters;
};

struct StaticWingInterfaceTriangle {
    StableId id = invalidStableId;
    std::array<std::size_t, 3> vertexIndices{};
    StableId negativeSideRegionId = invalidStableId;
    StableId positiveSideRegionId = invalidStableId;
    SurfaceRole role = SurfaceRole::Skin;
};

struct StaticWingInterfaceDiagnostics {
    std::uint32_t version = amrStaticWingInterfaceVersion;
    std::size_t vertexCount = 0;
    std::size_t triangleCount = 0;
    std::size_t coarseCandidateCount = 0;
    std::size_t refinedCandidateCount = 0;
    std::size_t coarseCutCellCount = 0;
    std::size_t refinedCutCellCount = 0;
    std::size_t activeCompositeCutCellCount = 0;
    double surfaceAreaSquareMeters = 0.0;
    double refinedPatchSurfaceAreaSquareMeters = 0.0;
    Vec3 lowerMeters;
    Vec3 upperMeters;
    bool entirelyInsideWindTunnel = false;
    bool finite = false;
    bool accepted = false;
};

// Immutable scene-v2 material-surface binding for the fixed-wing AMR
// checkpoint. Exact triangle/cell clipping supplies positive-area cut cells;
// no AMReX embedded-boundary, pressure jump, leakage, or region topology is
// invented here. The dense per-level tensors are a deliberately labelled
// direct-forcing approximation used by the first one-way momentum experiment.
class StaticWingInterface final {
public:
    StaticWingInterface(
        const Scene& scene,
        const WindTunnelGridSettings& gridSettings = {});

    [[nodiscard]] const StaticWingInterfaceDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] const StaticWingInterfaceLevel&
    level(std::size_t levelIndex) const;
    [[nodiscard]] const std::vector<StaticWingInterfaceVertex>&
    vertices() const noexcept;
    [[nodiscard]] const std::vector<StaticWingInterfaceTriangle>&
    triangles() const noexcept;

private:
    StaticWingInterfaceDiagnostics diagnostics_;
    std::array<StaticWingInterfaceLevel, 2> levels_;
    std::vector<StaticWingInterfaceVertex> vertices_;
    std::vector<StaticWingInterfaceTriangle> triangles_;
};

} // namespace simwing::fsi::amr
