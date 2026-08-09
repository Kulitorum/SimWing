#pragma once

#include "scene_fluid_pressure_coupling.h"
#include "structure_frame.h"
#include "viewer_protocol.h"

#include <cstdint>

namespace simwing::fsi {

inline constexpr char scenePressureCellCaseChecksum[] =
    "sha256:simwing-scene-pressure-feedback-cell-v2";
inline constexpr char scenePressureCellCaseSolverId[] =
    "simwing-fsi-scene-pressure-feedback-worker-v2";
inline constexpr std::uint32_t scenePressureCellCheckpointVersion = 2;

struct ScenePressureCellDiagnostics {
    SceneFluidPressureCouplingStepDiagnostics coupling;
    SceneFluidPressureMacVelocityCollapseDiagnostics macVelocity;
    double actuatorForceNewtons = 0.0;
    StructureVector3 pressureForceNewtons;
    double maximumAbsolutePressurePascals = 0.0;
    double maximumDisplacementMeters = 0.0;
    bool finite = false;

    bool operator==(const ScenePressureCellDiagnostics&) const = default;
};

struct ScenePressureCellCheckpoint {
    std::uint32_t version = scenePressureCellCheckpointVersion;
    SceneFluidPressureCouplingCheckpoint coupling;
};

// A visible analytic open cell driven by a bounded sinusoidal apex force. Its
// geometry is rebuilt by the scene-v2 cut-cell path and its pressure load is
// strongly coupled back to the same XPBD Structure. This is a diagnostic
// pressure-feedback canonical, not a complete momentum-evolving CFD model.
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
    [[nodiscard]] ScenePressureCellCheckpoint checkpoint() const;
    void restore(const ScenePressureCellCheckpoint& checkpoint);

private:
    Scene scene_;
    SceneFluidSurfaceAssembly surface_;
    SceneStructureAssembly assembly_;
    Structure structure_;
    SceneFluidPressureCoupling coupling_;
    fluid::MacVelocityField predictedVelocity_;
    viewer::StructureFrameMapping frameMapping_;
    std::size_t apexNode_ = 0;
    ScenePressureCellDiagnostics diagnostics_;
};

} // namespace simwing::fsi
