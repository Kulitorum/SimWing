#pragma once

#include "fluid/evolution.h"
#include "pressure_jump_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

struct MovingPorousFlowCaseCheckpointCodecAccess;

inline constexpr char movingPorousFlowCaseChecksum[] =
    "sha256:simwing-moving-planar-porous-flow-case-v1";
inline constexpr char movingPorousFlowCaseSolverId[] =
    "simwing-fsi-moving-planar-porous-flow-worker-v1";
inline constexpr std::uint32_t movingPorousFlowCaseCheckpointVersion = 1;
inline constexpr std::uint64_t movingPorousFlowCaseDefinitionFingerprint =
    0x6d70666c6f773031ULL;
inline constexpr double movingPorousFlowInitialSheetPositionMeters = 3.48;
inline constexpr double movingPorousFlowSheetVelocityMetersPerSecond = 0.4;
inline constexpr std::size_t movingPorousFlowInitialFaceCoordinate = 3;

// Immutable accepted worker state. Public orchestration metadata is mirrored
// inside the private owning payload so a partially edited checkpoint cannot be
// restored. Persistent encoding is deliberately a separate later boundary.
struct MovingPorousFlowCaseCheckpoint {
    std::uint32_t version = movingPorousFlowCaseCheckpointVersion;
    std::uint64_t caseDefinitionFingerprint =
        movingPorousFlowCaseDefinitionFingerprint;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::size_t pressureJumpCount = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
    double sheetPositionMeters = 0.0;
    double sheetVelocityMetersPerSecond = 0.0;
    std::uint64_t topologyRebaseCount = 0;
    fluid::MovingPorousFaceTopology porousTopology;

private:
    friend class MovingPorousFlowCase;
    friend struct MovingPorousFlowCaseCheckpointCodecAccess;
    struct Detail;
    std::shared_ptr<const Detail> detail;
};

// Visible Qt-free oracle for one prescribed translating porous plane inside a
// complete periodic flow step. Each advance samples the sheet at both porous
// half-step midpoints, admits retained or adjacent topology, and publishes the
// final wrapped crossing plane plus its unwrapped topology diagnostics. This
// is prescribed planar source motion, not fluid/structure coupling.
class MovingPorousFlowCase final {
public:
    MovingPorousFlowCase();

    MovingPorousFlowCase(const MovingPorousFlowCase&) = delete;
    MovingPorousFlowCase& operator=(const MovingPorousFlowCase&) = delete;
    MovingPorousFlowCase(MovingPorousFlowCase&&) = delete;
    MovingPorousFlowCase& operator=(MovingPorousFlowCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] MovingPorousFlowCaseCheckpoint checkpoint() const;
    void restore(const MovingPorousFlowCaseCheckpoint& checkpoint);

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::SharpPressureJumpField&
    pressureJumps() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSspRk2Settings&
    stepSettings() const noexcept;
    [[nodiscard]] const fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] const fluid::MovingPorousFaceTopology&
    porousTopology() const noexcept;
    [[nodiscard]] double sheetPositionMeters() const noexcept;
    [[nodiscard]] double sheetVelocityMetersPerSecond() const noexcept;
    [[nodiscard]] std::uint64_t topologyRebaseCount() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::SharpPressureJumpField pressureJumps_;
    fluid::PorousIterationSettings porousIteration_;
    fluid::PeriodicFlowStrangSspRk2Settings flowSettings_;
    fluid::MovingPlanarPorousFlowStrangSspRk2Diagnostics diagnostics_;
    fluid::MovingPorousFaceTopology porousTopology_;
    double sheetPositionMeters_ =
        movingPorousFlowInitialSheetPositionMeters;
    double sheetVelocityMetersPerSecond_ =
        movingPorousFlowSheetVelocityMetersPerSecond;
    std::uint64_t topologyRebaseCount_ = 0;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
