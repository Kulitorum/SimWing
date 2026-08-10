#pragma once

#include "fluid/grid.h"
#include "scene.h"
#include "structure.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace simwing::fsi {

inline constexpr const char* frozenScenePressureSolverId =
    "frozen-scene-mimetic-pressure-v1";

struct FrozenScenePressureCaseSettings {
    fluid::GridCellCounts cellCounts{2, 2, 2};
    double domainPaddingMeters = 0.5;
    bool useExplicitDomain = false;
    fluid::Vector3 lowerMeters;
    fluid::Vector3 upperMeters;
    fluid::Vector3 backgroundWindMetersPerSecond{-0.85, 0.0, 0.0};
    double diagnosticPerturbationSpeedMetersPerSecond = 1.0;
    double timeStepSeconds = 1.0 / 60.0;
    double densityKgPerCubicMeter = 1.225;
};

struct FrozenScenePressureCaseDiagnostics {
    std::size_t pressureControlCount = 0;
    std::size_t sharedTraceCount = 0;
    std::size_t pressureIterationCount = 0;
    double maximumAbsolutePressureDifferencePascals = 0.0;
    double maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
        0.0;
    double maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond =
        0.0;
    double correctedContinuityToleranceCubicMetersPerSecond = 0.0;
    double maximumCollapsedMacVelocityMetersPerSecond = 0.0;
    double maximumCollapsedSubfaceVelocityDeviationMetersPerSecond = 0.0;
    std::size_t embeddedOpeningTraceCount = 0;
    StructureVector3 pressureForceNewtons;
    StructureVector3 pressureMomentNewtonMeters;
    double transferForceResidualNewtons = 0.0;
    double transferMomentResidualNewtonMeters = 0.0;
    bool finite = false;

    bool operator==(
        const FrozenScenePressureCaseDiagnostics&) const = default;
};

// Loads one already authoritative scene-v2 payload, freezes its Structure
// geometry, and evaluates one deterministic mixed-hybrid pressure projection
// from a prescribed periodic MAC predictor. Repeated advance calls publish the
// same immutable physics at increasing diagnostic-frame times. This is a
// geometry/pressure integration probe, not time-resolved external flow, a wake,
// a lift polar, or a two-way FSI step.
class FrozenScenePressureCase final {
public:
    explicit FrozenScenePressureCase(
        Scene scene,
        const FrozenScenePressureCaseSettings& settings = {});
    ~FrozenScenePressureCase();

    FrozenScenePressureCase(const FrozenScenePressureCase&) = delete;
    FrozenScenePressureCase& operator=(
        const FrozenScenePressureCase&) = delete;
    FrozenScenePressureCase(FrozenScenePressureCase&&) noexcept;
    FrozenScenePressureCase& operator=(
        FrozenScenePressureCase&&) noexcept;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] const FrozenScenePressureCaseDiagnostics& diagnostics()
        const noexcept;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace simwing::fsi
