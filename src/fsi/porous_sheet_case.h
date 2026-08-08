#pragma once

#include "coupling.h"
#include "face_resolved_bridge.h"
#include "fluid/porous_interface.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

struct CoupledPorousSheetCheckpointCodecAccess;

inline constexpr char coupledPorousSheetCaseChecksum[] =
    "sha256:simwing-coupled-porous-sheet-case-v3";
inline constexpr char coupledPorousSheetCaseSolverId[] =
    "simwing-fsi-coupled-porous-sheet-worker-v3";
inline constexpr std::uint32_t coupledPorousSheetDiagnosticsVersion = 3;
inline constexpr std::uint32_t coupledPorousSheetCheckpointVersion = 3;
inline constexpr std::uint64_t coupledPorousSheetCaseFingerprint =
    0x5e8d7a31c4b2960fULL;
inline constexpr std::size_t coupledPorousSheetInitialFaceCoordinate = 3;
inline constexpr std::size_t coupledPorousSheetPumpFaceCoordinate = 7;
inline constexpr std::uint64_t
    coupledPorousSheetMaximumOrdinaryRebaseCount =
        coupledPorousSheetPumpFaceCoordinate
        - coupledPorousSheetInitialFaceCoordinate - 1;

struct CoupledPorousSheetStepDiagnostics {
    std::uint32_t version = coupledPorousSheetDiagnosticsVersion;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double sheetPositionBeforeMeters = 0.0;
    double sheetPositionAtConstitutiveTimeMeters = 0.0;
    double sheetPositionAfterMeters = 0.0;
    double fluidVelocityBeforeMetersPerSecond = 0.0;
    double fluidVelocityAtConstitutiveTimeMetersPerSecond = 0.0;
    double fluidVelocityAfterMetersPerSecond = 0.0;
    double sheetVelocityBeforeMetersPerSecond = 0.0;
    double sheetVelocityAtConstitutiveTimeMetersPerSecond = 0.0;
    double sheetVelocityAfterMetersPerSecond = 0.0;
    StructureVector3 actualFluidImpulseNewtonSeconds;
    StructureVector3 pressureJumpImpulseOnFluidNewtonSeconds;
    StructureVector3 actualSheetImpulseNewtonSeconds;
    StructureVector3 porousImpulseOnSheetNewtonSeconds;
    StructureVector3 pumpImpulseNewtonSeconds;
    StructureVector3 momentumResidualNewtonSeconds;
    double momentumResidualNormNewtonSeconds = 0.0;
    double fluidKineticEnergyChangeJoules = 0.0;
    double sheetKineticEnergyChangeJoules = 0.0;
    double pumpWorkJoules = 0.0;
    double porousDissipatedEnergyJoules = 0.0;
    double energyResidualJoules = 0.0;
    double maximumFluidUniformityResidualMetersPerSecond = 0.0;
    double maximumSheetRigidMotionResidualMeters = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    std::size_t porousFaceCoordinate =
        coupledPorousSheetInitialFaceCoordinate;
    bool topologyRebasedThisStep = false;
    fluid::PorousProjectionDiagnostics fluidProjection;
    fluid::PorousSurfaceTractionDiagnostics porousTraction;
    PorousFaceResolvedBridgeDiagnostics bridge;
    TimeIntegratedTransferDiagnostics transfer;
    bool finite = true;
    bool accepted = false;

    bool operator==(
        const CoupledPorousSheetStepDiagnostics&) const = default;
};

struct CoupledPorousSheetCheckpoint {
    std::uint32_t version = coupledPorousSheetCheckpointVersion;
    std::uint64_t caseFingerprint = coupledPorousSheetCaseFingerprint;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    std::size_t porousFaceCoordinate =
        coupledPorousSheetInitialFaceCoordinate;

private:
    friend class CoupledPorousSheetCase;
    friend struct CoupledPorousSheetCheckpointCodecAccess;
    struct Detail;
    std::shared_ptr<const Detail> detail;
};

// A one-degree-of-freedom coupled verification oracle. A
// periodic pump drives flow through one planar Darcy sheet. The midpoint
// porous solve owns fluid pressure work and dissipation; the stable-ID bridge
// transfers only the equal-and-opposite sheet reaction into XPBD. The scalar
// midpoint relation is analytic for this linear material, so fluid, sheet,
// pump, and dissipation ledgers can be checked independently. The sheet may
// rebind to the next MAC face across several dual-cell boundaries while
// preserving physical position and accepted fluid state. It rejects the later
// collision with the pump surface. General strong coupling remains outside
// this case.
class CoupledPorousSheetCase final {
public:
    CoupledPorousSheetCase();

    CoupledPorousSheetCase(const CoupledPorousSheetCase&) = delete;
    CoupledPorousSheetCase& operator=(const CoupledPorousSheetCase&) = delete;
    CoupledPorousSheetCase(CoupledPorousSheetCase&&) = delete;
    CoupledPorousSheetCase& operator=(CoupledPorousSheetCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] CoupledPorousSheetCheckpoint checkpoint() const;
    void restore(const CoupledPorousSheetCheckpoint& checkpoint);

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] const CoupledPorousSheetStepDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] std::size_t porousFaceCoordinate() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::PorousProjectionSettings projectionSettings_;
    Structure structure_;
    PlanarFaceResolvedFluidStructureBridge bridge_;
    ConservativeMacroStepCoupling coupling_;
    viewer::StructureFrameMapping frameMapping_;
    StructureStepSettings stepSettings_;
    CoupledPorousSheetStepDiagnostics diagnostics_;
    std::uint64_t topologyRebaseCount_ = 0;
    std::size_t porousFaceCoordinate_ =
        coupledPorousSheetInitialFaceCoordinate;
};

} // namespace simwing::fsi
