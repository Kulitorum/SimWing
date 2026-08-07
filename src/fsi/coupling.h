#pragma once

#include "transfer.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t interfaceImpulseExchangeVersion = 1;

struct CouplingNodeImpulse {
    std::uint64_t stableId = 0;
    std::size_t structureNode = 0;
    StructureVector3 impulseNewtonSeconds;

    bool operator==(const CouplingNodeImpulse&) const = default;
};

// Instantaneous surface and nodal ledgers are trapezoidally integrated in
// macro-step-local time. The two sides remain independent so these residuals
// diagnose the actual temporal transfer as well as the spatial transfer.
struct TimeIntegratedTransferDiagnostics {
    std::size_t sampleCount = 0;
    std::size_t intervalCount = 0;
    double durationSeconds = 0.0;
    StructureVector3 momentReferenceMeters;
    StructureVector3 integratedSurfaceImpulseNewtonSeconds;
    StructureVector3 transferredNodalImpulseNewtonSeconds;
    StructureVector3 impulseResidualNewtonSeconds;
    double impulseResidualNormNewtonSeconds = 0.0;
    StructureVector3 integratedSurfaceAngularImpulseNewtonMeterSeconds;
    StructureVector3 transferredNodalAngularImpulseNewtonMeterSeconds;
    StructureVector3 angularImpulseResidualNewtonMeterSeconds;
    double angularImpulseResidualNormNewtonMeterSeconds = 0.0;
    double integratedSurfaceWorkJoules = 0.0;
    double transferredNodalWorkJoules = 0.0;
    double workResidualJoules = 0.0;
    bool finite = true;

    bool operator==(const TimeIntegratedTransferDiagnostics&) const = default;
};

class TimeIntegratedTransferResult final {
public:
    [[nodiscard]] std::uint32_t version() const noexcept;
    [[nodiscard]] std::uint64_t surfaceFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t targetDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::span<const CouplingNodeImpulse>
    nodeImpulses() const noexcept;
    [[nodiscard]] const TimeIntegratedTransferDiagnostics&
    diagnostics() const noexcept;

    bool operator==(const TimeIntegratedTransferResult&) const = default;

private:
    friend class ConservativeMacroStepCoupling;
    TimeIntegratedTransferResult() = default;

    std::uint32_t version_ = interfaceImpulseExchangeVersion;
    std::uint64_t surfaceFingerprint_ = 0;
    std::uint64_t targetDefinitionFingerprint_ = 0;
    std::vector<CouplingNodeImpulse> nodeImpulses_;
    TimeIntegratedTransferDiagnostics diagnostics_;
};

// Integrates immutable instantaneous transfer samples over one fluid
// macro-step. Sample times are offsets: the first must be exactly zero and the
// last is the macro-step duration. Structural acceptance applies the resulting
// impulse as one equivalent average nodal load across Structure's internal
// substeps. It does not decide strong-coupling convergence.
class ConservativeMacroStepCoupling final {
public:
    explicit ConservativeMacroStepCoupling(
        const ConservativeSurfaceTransfer& transfer);

    [[nodiscard]] std::uint64_t surfaceFingerprint() const noexcept;
    [[nodiscard]] std::uint64_t targetDefinitionFingerprint() const noexcept;

    [[nodiscard]] TimeIntegratedTransferResult integrate(
        std::span<const double> sampleOffsetsSeconds,
        std::span<const ConservativeTransferResult> samples) const;

    // Validates the complete exchange before mutation. On any structural-step
    // failure, the exact checkpoint from before load application is restored.
    [[nodiscard]] StructureDiagnostics advanceStructure(
        Structure& target,
        const TimeIntegratedTransferResult& transfer,
        const StructureStepSettings& settings) const;

private:
    std::uint64_t surfaceFingerprint_ = 0;
    std::uint64_t targetDefinitionFingerprint_ = 0;
    std::vector<CouplingSurfaceNodeDefinition> nodes_;
};

} // namespace simwing::fsi
