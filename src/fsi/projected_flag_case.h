#pragma once

#include "face_resolved_bridge.h"
#include "fluid/moving_interface.h"
#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <vector>

namespace simwing::fsi {

inline constexpr char projectedGustFlagCaseChecksum[] =
    "sha256:simwing-fixed-reference-projected-gust-flag-v1";
inline constexpr char projectedGustFlagCaseSolverId[] =
    "simwing-fsi-fixed-reference-projected-gust-flag-worker-v1";
inline constexpr std::size_t projectedFlagTilesPerSide = 4;
inline constexpr std::size_t projectedFlagNodesPerSide =
    projectedFlagTilesPerSide + 1;
inline constexpr std::uint64_t projectedFlagSurfaceStableId = 41'000;

// First fluid-to-flexible-fabric canonical. A finite, grid-aligned reference
// panel diverts an accelerating periodic cross-flow through an incompressible
// MAC projection. The complete fluid constraint reaction (adjacent pressure
// plus the direct constrained-face impulse) is conservatively mapped to a
// edge-clamped XPBD membrane. Two node rows encode both position and slope at
// that clamp, instead of leaving an unintended rigid hinge mode. Structural
// displacement is deliberately not
// fed back to the fixed CFD reference surface: this is an honest one-way
// fixed-reference slice, not a moving cut-cell or two-way energy claim.
class ProjectedGustFlagCase final {
public:
    ProjectedGustFlagCase();

    ProjectedGustFlagCase(const ProjectedGustFlagCase&) = delete;
    ProjectedGustFlagCase& operator=(const ProjectedGustFlagCase&) = delete;
    ProjectedGustFlagCase(ProjectedGustFlagCase&&) = delete;
    ProjectedGustFlagCase& operator=(ProjectedGustFlagCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::MovingInterfaceProjectionDiagnostics&
    fluidDiagnostics() const noexcept;
    [[nodiscard]] const PlanarFaceResolvedBridgeDiagnostics&
    transferDiagnostics() const noexcept;
    [[nodiscard]] double gustSpeedMetersPerSecond() const noexcept;
    [[nodiscard]] double maximumNormalDisplacementMeters() const;
    [[nodiscard]] double maximumFreeEdgeDisplacementMeters() const;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::FaceAlignedMovingInterface interface_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    Structure structure_;
    fluid::MovingInterfaceProjectionSettings fluidSettings_;
    PlanarFaceResolvedFluidStructureBridge bridge_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    std::vector<CouplingNodeKinematics> referenceKinematics_;
    fluid::MovingInterfaceProjectionDiagnostics fluidDiagnostics_;
    PlanarFaceResolvedBridgeDiagnostics transferDiagnostics_;
    std::vector<StructureVector3> lastNodalForcesNewtons_;
    double gustSpeedMetersPerSecond_ = 0.0;
};

} // namespace simwing::fsi
