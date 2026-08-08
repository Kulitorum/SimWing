#pragma once

#include "fluid/evolution.h"
#include "fluid_frame.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

inline constexpr char periodicFlowCaseChecksum[] =
    "sha256:simwing-periodic-taylor-green-flow-case-v1";
inline constexpr char periodicFlowCaseSolverId[] =
    "simwing-fsi-periodic-flow-worker-v1";
inline constexpr std::uint32_t periodicFlowCaseCheckpointVersion = 1;
inline constexpr std::uint64_t periodicFlowCaseDefinitionFingerprint =
    0x706572666c6f7731ull;

// Immutable in-memory state for one committed periodic-flow worker step. The
// public metadata is suitable for orchestration; the private payload prevents
// callers from changing numerical arrays independently of their identity.
struct PeriodicFlowCaseCheckpoint {
    std::uint32_t version = periodicFlowCaseCheckpointVersion;
    std::uint64_t caseDefinitionFingerprint =
        periodicFlowCaseDefinitionFingerprint;
    fluid::GridCellCounts cellCounts;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    std::size_t scalarSampleCount = 0;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;

private:
    friend class PeriodicFlowCase;
    struct Detail;
    std::shared_ptr<const Detail> detail;
};

// A visible, deterministic CFD verification case. It advances a smooth
// Galilean-shifted Taylor-Green field through the bounded periodic
// Strang/SSPRK2 subcycler and publishes only accepted cell-centred snapshots.
// It is deliberately not a canopy, cut-cell, or aerodynamic-truth case.
class PeriodicFlowCase final {
public:
    PeriodicFlowCase();

    PeriodicFlowCase(const PeriodicFlowCase&) = delete;
    PeriodicFlowCase& operator=(const PeriodicFlowCase&) = delete;
    PeriodicFlowCase(PeriodicFlowCase&&) = delete;
    PeriodicFlowCase& operator=(PeriodicFlowCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] PeriodicFlowCaseCheckpoint checkpoint() const;
    void restore(const PeriodicFlowCaseCheckpoint& checkpoint);

    [[nodiscard]] const fluid::PeriodicCartesianGrid& grid() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& velocity() const noexcept;
    [[nodiscard]] const fluid::CellScalarField& pressure() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSubcyclingSettings&
    stepSettings() const noexcept;
    [[nodiscard]] const fluid::PeriodicFlowStrangSubcyclingDiagnostics&
    diagnostics() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;

private:
    fluid::PeriodicCartesianGrid grid_;
    fluid::MacVelocityField velocity_;
    fluid::CellScalarField pressure_;
    fluid::PeriodicFlowStrangSubcyclingSettings stepSettings_;
    fluid::PeriodicFlowStrangSubcyclingDiagnostics diagnostics_;
    std::uint64_t acceptedStepCount_ = 0;
    double simulationTimeSeconds_ = 0.0;
};

} // namespace simwing::fsi
