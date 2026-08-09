#pragma once

#include "fluid/evolution.h"
#include "scene_fluid_pressure_coupling.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char scenePressureCellCaseChecksum[] =
    "sha256:simwing-scene-pressure-feedback-cell-v9";
inline constexpr char scenePressureCellCaseSolverId[] =
    "simwing-fsi-scene-pressure-feedback-worker-v9";
inline constexpr std::uint32_t scenePressureCellCheckpointVersion = 9;

struct ScenePressureCellDiagnostics {
    SceneFluidPressureCouplingStepDiagnostics coupling;
    SceneFluidPressureMacVelocityCollapseDiagnostics macVelocity;
    fluid::PeriodicFlowStrangSspRk2Diagnostics bulkFlow;
    SceneFluidRegionTransportDiagnostics regionTransport;
    bool usesRegionTransport = false;
    double targetMeanWindMetersPerSecond = 0.0;
    double meanWindBeforePumpMetersPerSecond = 0.0;
    double flowPumpForceNewtons = 0.0;
    StructureVector3 pressureForceNewtons;
    StructureVector3 wallForceNewtons;
    StructureVector3 totalFluidForceNewtons;
    double maximumAbsolutePressurePascals = 0.0;
    double maximumDisplacementMeters = 0.0;
    bool finite = false;

    bool operator==(const ScenePressureCellDiagnostics&) const = default;
};

struct ScenePressureCellCheckpoint {
    std::uint32_t version = scenePressureCellCheckpointVersion;
    SceneFluidPressureCouplingCheckpoint coupling;
    std::optional<SceneFluidRegionMomentumState> regionMomentum;
};

// A visible analytic open cell driven by a prescribed periodic mean-flow
// pump. The accepted collapsed MAC state receives one existing symmetric
// viscosity/projected-nonlinear-advection/viscosity bulk step before scene-v2
// pressure and conservative tangential material-wall exchange are strongly
// coupled back to the same XPBD Structure.
// Immersed-boundary momentum remains an approximation: cut-region pressure is
// exact for this subset, but the bulk stages do not yet resolve distinct
// per-region subface advection. Wall viscosity is a local cut-region closure,
// not a general immersed-boundary boundary-layer model.
class ScenePressureCellCase final {
public:
    ScenePressureCellCase();

    ScenePressureCellCase(const ScenePressureCellCase&) = delete;
    ScenePressureCellCase& operator=(const ScenePressureCellCase&) = delete;
    ScenePressureCellCase(ScenePressureCellCase&&) = delete;
    ScenePressureCellCase& operator=(ScenePressureCellCase&&) = delete;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();

    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] const fluid::MacVelocityField& predictedVelocity()
        const noexcept;
    [[nodiscard]] const ScenePressureCellDiagnostics& diagnostics()
        const noexcept;
    [[nodiscard]] const SceneFluidRegionMomentumState* acceptedRegionMomentum()
        const noexcept;
    [[nodiscard]] ScenePressureCellCheckpoint checkpoint() const;
    void restore(const ScenePressureCellCheckpoint& checkpoint);

private:
    Scene scene_;
    SceneFluidSurfaceAssembly surface_;
    SceneStructureAssembly assembly_;
    Structure structure_;
    SceneFluidPressureCoupling coupling_;
    fluid::MacVelocityField predictedVelocity_;
    std::optional<SceneFluidRegionMomentumState> acceptedRegionMomentum_;
    viewer::StructureFrameMapping frameMapping_;
    ScenePressureCellDiagnostics diagnostics_;
};

} // namespace simwing::fsi
