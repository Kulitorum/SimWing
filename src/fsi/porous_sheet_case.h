#pragma once

#include "coupling.h"
#include "face_resolved_bridge.h"
#include "fluid/porous_interface.h"
#include "fluid/porous_topology.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

struct CoupledPorousSheetCheckpointCodecAccess;

enum class CoupledPorousSheetMotionDirection : std::int8_t {
    Negative = -1,
    Positive = 1,
};

inline constexpr char coupledPorousSheetCaseChecksum[] =
    "sha256:simwing-coupled-porous-sheet-case-v5";
inline constexpr char coupledPorousSheetCaseSolverId[] =
    "simwing-fsi-coupled-porous-sheet-worker-v5";
inline constexpr char coupledPorousSheetReverseCaseChecksum[] =
    "sha256:simwing-coupled-porous-sheet-reverse-case-v1";
inline constexpr char coupledPorousSheetReverseCaseSolverId[] =
    "simwing-fsi-coupled-porous-sheet-reverse-worker-v1";
inline constexpr std::uint32_t coupledPorousSheetDiagnosticsVersion = 5;
inline constexpr std::uint32_t coupledPorousSheetCheckpointVersion = 5;
inline constexpr std::uint64_t coupledPorousSheetCaseFingerprint =
    0xc4196d3b85e72a0fULL;
inline constexpr std::uint64_t coupledPorousSheetReverseCaseFingerprint =
    0x7d2a61f9b408ce35ULL;
inline constexpr std::size_t coupledPorousSheetGridFaceCount = 8;
inline constexpr std::size_t coupledPorousSheetInitialFaceCoordinate = 3;
inline constexpr std::size_t coupledPorousSheetPumpFaceCoordinate = 2;
inline constexpr std::size_t coupledPorousSheetReversePumpFaceCoordinate = 4;
inline constexpr std::uint64_t
    coupledPorousSheetMaximumOrdinaryRebaseCount = 6;
[[nodiscard]] inline constexpr fluid::MovingPorousFaceTopology
coupledPorousSheetTopologyAfterRebases(
    const CoupledPorousSheetMotionDirection direction,
    const std::uint64_t rebaseCount) {
    const std::int64_t unwrappedFaceCoordinate =
        static_cast<std::int64_t>(
            coupledPorousSheetInitialFaceCoordinate)
        + static_cast<std::int64_t>(direction)
            * static_cast<std::int64_t>(rebaseCount);
    const std::int64_t faceCount = static_cast<std::int64_t>(
        coupledPorousSheetGridFaceCount);
    const std::int64_t periodicImage = unwrappedFaceCoordinate >= 0
        ? unwrappedFaceCoordinate / faceCount
        : -((-unwrappedFaceCoordinate + faceCount - 1) / faceCount);
    return {
        fluid::movingPorousFaceTopologyVersion,
        fluid::GridFaceAxis::X,
        static_cast<std::size_t>(
            unwrappedFaceCoordinate - periodicImage * faceCount),
        periodicImage,
    };
}
[[nodiscard]] inline constexpr fluid::MovingPorousFaceTopology
coupledPorousSheetTopologyAfterRebases(const std::uint64_t rebaseCount) {
    return coupledPorousSheetTopologyAfterRebases(
        CoupledPorousSheetMotionDirection::Positive, rebaseCount);
}
inline constexpr fluid::MovingPorousFaceTopology
    coupledPorousSheetInitialTopology =
        coupledPorousSheetTopologyAfterRebases(0);
inline constexpr fluid::MovingPorousFaceTopology
    coupledPorousSheetTerminalSafeTopology =
        coupledPorousSheetTopologyAfterRebases(
            coupledPorousSheetMaximumOrdinaryRebaseCount);
inline constexpr fluid::MovingPorousFaceTopology
    coupledPorousSheetPumpCollisionTopology =
        coupledPorousSheetTopologyAfterRebases(
            coupledPorousSheetMaximumOrdinaryRebaseCount + 1);
inline constexpr fluid::MovingPorousFaceTopology
    coupledPorousSheetReverseTerminalSafeTopology =
        coupledPorousSheetTopologyAfterRebases(
            CoupledPorousSheetMotionDirection::Negative,
            coupledPorousSheetMaximumOrdinaryRebaseCount);
inline constexpr fluid::MovingPorousFaceTopology
    coupledPorousSheetReversePumpCollisionTopology =
        coupledPorousSheetTopologyAfterRebases(
            CoupledPorousSheetMotionDirection::Negative,
            coupledPorousSheetMaximumOrdinaryRebaseCount + 1);
static_assert(
    coupledPorousSheetPumpCollisionTopology.faceCoordinate
            == coupledPorousSheetPumpFaceCoordinate
        && coupledPorousSheetPumpCollisionTopology.periodicImage == 1);
static_assert(
    coupledPorousSheetReversePumpCollisionTopology.faceCoordinate
            == coupledPorousSheetReversePumpFaceCoordinate
        && coupledPorousSheetReversePumpCollisionTopology.periodicImage
            == -1);

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
    fluid::MovingPorousFaceTopology porousTopology =
        coupledPorousSheetInitialTopology;
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
    fluid::MovingPorousFaceTopology porousTopology =
        coupledPorousSheetInitialTopology;

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
// rebind to the next MAC face across several dual-cell boundaries and one
// periodic wrap while preserving unwrapped physical position and accepted
// fluid state. It rejects the later collision with the next periodic image of
// the pump surface. A focused negative-direction instance reverses the pump,
// momentum, rebases, periodic image, provenance, and checkpoint identity
// through the same boundary. General strong coupling remains outside this case.
class CoupledPorousSheetCase final {
public:
    explicit CoupledPorousSheetCase(
        CoupledPorousSheetMotionDirection direction =
            CoupledPorousSheetMotionDirection::Positive);

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
    [[nodiscard]] CoupledPorousSheetMotionDirection
    motionDirection() const noexcept;
    [[nodiscard]] std::uint64_t caseFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] const fluid::MovingPorousFaceTopology&
    porousTopology() const noexcept;
    [[nodiscard]] std::size_t porousFaceCoordinate() const noexcept;

private:
    CoupledPorousSheetMotionDirection direction_ =
        CoupledPorousSheetMotionDirection::Positive;
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
    fluid::MovingPorousFaceTopology porousTopology_ =
        coupledPorousSheetInitialTopology;
};

} // namespace simwing::fsi
