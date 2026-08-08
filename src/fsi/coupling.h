#pragma once

#include "transfer.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi {

inline constexpr std::uint32_t interfaceImpulseExchangeVersion = 1;
inline constexpr std::uint32_t aitkenRelaxationCheckpointVersion = 1;
inline constexpr std::uint32_t strongCouplingIterationCheckpointVersion = 1;

struct AitkenRelaxationSettings {
    double initialRelaxation = 0.5;
    double minimumRelaxation = 0.05;
    double maximumRelaxation = 1.0;

    bool operator==(const AitkenRelaxationSettings&) const = default;
};

struct AitkenRelaxationDiagnostics {
    std::uint64_t completedIterationCount = 0;
    double relaxation = 0.0;
    double residualL2 = 0.0;
    double residualChangeL2 = 0.0;
    double denominator = 0.0;
    bool usedDynamicRelaxation = false;
    bool relaxationWasClipped = false;
    bool finite = true;

    bool operator==(const AitkenRelaxationDiagnostics&) const = default;
};

struct AitkenRelaxationCheckpoint {
    std::uint32_t version = aitkenRelaxationCheckpointVersion;
    std::uint64_t interfaceDefinitionFingerprint = 0;
    std::size_t valueCount = 0;
    AitkenRelaxationSettings settings;
    std::uint64_t completedIterationCount = 0;
    double relaxation = 0.0;
    std::vector<double> previousResidual;

    bool operator==(const AitkenRelaxationCheckpoint&) const = default;
};

// Stateful dynamic relaxation for one consistently scaled interface vector.
// The nonzero definition fingerprint binds the caller's component ordering and
// scaling contract, not merely its value count.
// The caller supplies the current relaxed iterate and the newly solved
// unrelaxed candidate. The first update uses the configured fixed factor;
// later updates use the vector Aitken delta-squared factor, clipped to explicit
// bounds. Checkpoints make macro-step iteration rollback exact. Invalid inputs
// and incompatible restores mutate neither the state nor the output vector.
class AitkenInterfaceRelaxation final {
public:
    explicit AitkenInterfaceRelaxation(
        std::uint64_t interfaceDefinitionFingerprint,
        std::size_t valueCount,
        const AitkenRelaxationSettings& settings = {});

    [[nodiscard]] std::uint64_t
    interfaceDefinitionFingerprint() const noexcept;
    [[nodiscard]] std::size_t valueCount() const noexcept;
    [[nodiscard]] const AitkenRelaxationSettings& settings() const noexcept;
    [[nodiscard]] std::uint64_t completedIterationCount() const noexcept;
    [[nodiscard]] double relaxation() const noexcept;
    [[nodiscard]] AitkenRelaxationCheckpoint checkpoint() const;

    void restore(const AitkenRelaxationCheckpoint& checkpoint);
    void reset() noexcept;

    [[nodiscard]] AitkenRelaxationDiagnostics relax(
        std::span<const double> current,
        std::span<const double> candidate,
        std::vector<double>& relaxed);

private:
    std::uint64_t interfaceDefinitionFingerprint_ = 0;
    std::size_t valueCount_ = 0;
    AitkenRelaxationSettings settings_;
    std::uint64_t completedIterationCount_ = 0;
    double relaxation_ = 0.0;
    std::vector<double> previousResidual_;
};

struct CouplingResidualNorms {
    double displacementMetres = 0.0;
    double displacementReferenceMetres = 0.0;
    double velocityMetersPerSecond = 0.0;
    double velocityReferenceMetersPerSecond = 0.0;
    double tractionNewtons = 0.0;
    double tractionReferenceNewtons = 0.0;

    bool operator==(const CouplingResidualNorms&) const = default;
};

struct CouplingConvergenceSettings {
    std::uint64_t minimumIterations = 2;
    std::uint64_t maximumIterations = 25;
    double absoluteDisplacementToleranceMetres = 1.0e-6;
    double relativeDisplacementTolerance = 1.0e-3;
    double displacementReferenceFloorMetres = 1.0e-3;
    double absoluteVelocityToleranceMetersPerSecond = 1.0e-5;
    double relativeVelocityTolerance = 1.0e-3;
    double velocityReferenceFloorMetersPerSecond = 1.0e-2;
    double absoluteTractionToleranceNewtons = 1.0e-3;
    double relativeTractionTolerance = 1.0e-3;
    double tractionReferenceFloorNewtons = 1.0;

    bool operator==(const CouplingConvergenceSettings&) const = default;
};

struct CouplingConvergenceDecision {
    std::uint64_t iteration = 0;
    CouplingResidualNorms residuals;
    double relativeDisplacement = 0.0;
    double relativeVelocity = 0.0;
    double relativeTraction = 0.0;
    bool displacementConverged = false;
    bool velocityConverged = false;
    bool tractionConverged = false;
    bool minimumIterationsSatisfied = false;
    bool converged = false;
    bool iterationLimitReached = false;
    bool finite = true;

    bool operator==(const CouplingConvergenceDecision&) const = default;
};

// Evaluates an already reduced interface residual set. Every displacement,
// velocity, and traction channel must pass both its absolute and relative
// tolerance; no single convenient residual can accept the iteration. Relative
// references use explicit positive floors. Reaching the maximum iteration is
// reported independently from convergence so the caller can roll back and
// reduce the macro-step.
[[nodiscard]] CouplingConvergenceDecision evaluateCouplingConvergence(
    std::uint64_t iteration,
    const CouplingResidualNorms& residuals,
    const CouplingConvergenceSettings& settings = {});

enum class StrongCouplingIterationStatus : std::uint8_t {
    Iterating = 1,
    Converged = 2,
    Exhausted = 3,
};

struct StrongCouplingIterationResult {
    StrongCouplingIterationStatus status =
        StrongCouplingIterationStatus::Iterating;
    AitkenRelaxationDiagnostics relaxation;
    CouplingConvergenceDecision convergence;

    bool operator==(const StrongCouplingIterationResult&) const = default;
};

struct StrongCouplingIterationCheckpoint {
    std::uint32_t version = strongCouplingIterationCheckpointVersion;
    CouplingConvergenceSettings convergenceSettings;
    AitkenRelaxationCheckpoint relaxation;
    std::vector<double> currentInterface;
    StrongCouplingIterationStatus status =
        StrongCouplingIterationStatus::Iterating;
    AitkenRelaxationDiagnostics lastRelaxation;
    CouplingConvergenceDecision lastConvergence;

    bool operator==(
        const StrongCouplingIterationCheckpoint&) const = default;
};

// Owns one macro-step's relaxed interface iterate and terminal convergence
// state. Solver checkpoints remain in their solver owners; this coordinator's
// checkpoint composes the algorithm state they must save beside them. Advance
// is transactional and terminal after convergence or budget exhaustion.
class StrongCouplingIteration final {
public:
    StrongCouplingIteration(
        std::uint64_t interfaceDefinitionFingerprint,
        std::span<const double> initialInterface,
        const AitkenRelaxationSettings& relaxationSettings = {},
        const CouplingConvergenceSettings& convergenceSettings = {});

    [[nodiscard]] std::span<const double> currentInterface() const noexcept;
    [[nodiscard]] StrongCouplingIterationStatus status() const noexcept;
    [[nodiscard]] std::uint64_t completedIterationCount() const noexcept;
    [[nodiscard]] const AitkenRelaxationDiagnostics&
    lastRelaxation() const noexcept;
    [[nodiscard]] const CouplingConvergenceDecision&
    lastConvergence() const noexcept;
    [[nodiscard]] StrongCouplingIterationCheckpoint checkpoint() const;

    [[nodiscard]] StrongCouplingIterationResult advance(
        std::span<const double> unrelaxedCandidate,
        const CouplingResidualNorms& residuals);
    void restore(const StrongCouplingIterationCheckpoint& checkpoint);
    void reset(std::span<const double> initialInterface);

private:
    CouplingConvergenceSettings convergenceSettings_;
    AitkenInterfaceRelaxation relaxation_;
    std::vector<double> currentInterface_;
    StrongCouplingIterationStatus status_ =
        StrongCouplingIterationStatus::Iterating;
    AitkenRelaxationDiagnostics lastRelaxation_;
    CouplingConvergenceDecision lastConvergence_;
};

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

    // Reduces three topology-bound interface iterates and two corresponding
    // traction transfers into origin- and Galilean-invariant maximum nodal
    // residuals. Displacement/velocity references are measured as changes
    // from the saved macro-step baseline; traction references use physical
    // nodal-force magnitudes.
    [[nodiscard]] CouplingResidualNorms measureResiduals(
        std::span<const CouplingNodeKinematics> baselineKinematics,
        std::span<const CouplingNodeKinematics> previousKinematics,
        std::span<const CouplingNodeKinematics> currentKinematics,
        const ConservativeTransferResult& previousTraction,
        const ConservativeTransferResult& currentTraction) const;

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
