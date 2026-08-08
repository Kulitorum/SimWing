#pragma once

#include "face_resolved_bridge.h"
#include "fluid/moving_interface.h"
#include "structure.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace simwing::fsi {

inline constexpr char ramAirCellCaseChecksum[] =
    "sha256:simwing-fixed-reference-open-ram-air-cell-v1";
inline constexpr char ramAirCellCaseSolverId[] =
    "simwing-fsi-fixed-reference-open-ram-air-cell-worker-v1";
inline constexpr std::size_t ramAirCellTilesPerEdge = 4;
inline constexpr std::size_t ramAirCellNodesPerEdge =
    ramAirCellTilesPerEdge + 1;
inline constexpr std::size_t ramAirCellPanelCount = 5;

struct RamAirCellDiagnostics {
    std::uint32_t version = 1;
    std::vector<PlanarFaceResolvedBridgeDiagnostics> panels;
    StructureVector3 fluidPressureForceNewtons;
    StructureVector3 fluidReactionForceNewtons;
    StructureVector3 transferredForceNewtons;
    StructureVector3 forceResidualNewtons;
    StructureVector3 fluidReactionMomentNewtonMeters;
    StructureVector3 transferredMomentNewtonMeters;
    StructureVector3 momentResidualNewtonMeters;
    double maximumPanelForceResidualNewtons = 0.0;
    double maximumPanelMomentResidualNewtonMeters = 0.0;
    double fluidDivergenceL2PerSecond = 0.0;
    bool finite = true;

    bool operator==(const RamAirCellDiagnostics&) const = default;
};

// A five-panel rectangular cell with an open leading face. A persistent
// periodic MAC field is accelerated through and around the fixed reference
// cavity; complete constraint reactions on the back, side, top, and bottom
// panels are independently mapped to one shared-node XPBD shell. The mouth's
// first two perimeter rows encode clamp position and slope. Displaced fabric
// is deliberately not returned to the CFD topology, so this remains a one-way
// fixed-reference inflation/deformation gate rather than moving-cut-cell FSI.
class RamAirCellCase final {
public:
    RamAirCellCase();
    ~RamAirCellCase();

    RamAirCellCase(const RamAirCellCase&) = delete;
    RamAirCellCase& operator=(const RamAirCellCase&) = delete;
    RamAirCellCase(RamAirCellCase&&) = delete;
    RamAirCellCase& operator=(RamAirCellCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::MovingInterfaceProjectionDiagnostics&
    fluidDiagnostics() const noexcept;
    [[nodiscard]] const RamAirCellDiagnostics& diagnostics() const noexcept;
    [[nodiscard]] double gustSpeedMetersPerSecond() const noexcept;
    [[nodiscard]] double openingMeanVelocityMetersPerSecond() const;
    [[nodiscard]] double openingRmsVelocityMetersPerSecond() const;
    [[nodiscard]] double maximumDisplacementMeters() const;
    [[nodiscard]] double maximumOutwardInflationMeters() const;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::FaceAlignedMovingInterface interface_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    Structure structure_;
    fluid::MovingInterfaceProjectionSettings fluidSettings_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    std::vector<std::unique_ptr<PlanarFaceResolvedFluidStructureBridge>>
        bridges_;
    std::vector<std::vector<CouplingNodeKinematics>> referenceKinematics_;
    fluid::MovingInterfaceProjectionDiagnostics fluidDiagnostics_;
    RamAirCellDiagnostics diagnostics_;
    std::vector<StructureVector3> lastNodalForcesNewtons_;
    double gustSpeedMetersPerSecond_ = 0.0;
};

} // namespace simwing::fsi
