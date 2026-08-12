#pragma once

#include "scene.h"
#include "scene_structure.h"
#include "structure.h"
#include "surface_aerodynamics.h"
#include "viewer_protocol.h"

#include <cstddef>
#include <memory>

namespace simwing::fsi {

inline constexpr char surfaceFlightSolverId[] =
    "simwing-xpbd-surface-flight-v1";

struct SurfaceFlightCaseSettings {
    SurfaceFlightCaseSettings();

    SceneStructureSettings sceneStructure;
    StructureStepSettings structureStep;
    SurfaceAerodynamicsSettings aerodynamics;
    double fabricStiffnessScale = 1.0;
    // The first real-wing milestone applies the reduced-order net surface
    // force by fabric mass. This preserves the authored XPBD load path while
    // avoiding unresolved per-triangle shock modes in the very fine export.
    // Raw pressure traction remains fully diagnosed for the next fidelity step.
    bool massWeightedFlightPreviewLoads = true;
};

struct SurfaceFlightCaseDiagnostics {
    SurfaceAerodynamicsDiagnostics aerodynamics;
    ConservativeTransferDiagnostics transfer;
    StructureDiagnostics structure;
    double maximumSurfaceLoadAccelerationMetersPerSecondSquared = 0.0;
    bool finite = false;
};

// Interactive real-wing path: authored scene-v2 topology is assembled once,
// aerodynamic and pneumatic traction is evaluated on the accepted moving
// triangles, transferred conservatively, and committed with one XPBD step.
// The candidate gas state is committed only after the structure accepts.
class SurfaceFlightCase final {
public:
    explicit SurfaceFlightCase(
        Scene scene,
        const SurfaceFlightCaseSettings& settings = {});
    ~SurfaceFlightCase();

    SurfaceFlightCase(const SurfaceFlightCase&) = delete;
    SurfaceFlightCase& operator=(const SurfaceFlightCase&) = delete;
    SurfaceFlightCase(SurfaceFlightCase&&) noexcept;
    SurfaceFlightCase& operator=(SurfaceFlightCase&&) noexcept;

    [[nodiscard]] viewer::TraceHeader traceHeader() const;
    [[nodiscard]] viewer::DiagnosticFrame advance();
    [[nodiscard]] const Structure& structure() const noexcept;
    [[nodiscard]] const StructureStepSettings& stepSettings() const noexcept;
    [[nodiscard]] std::uint64_t acceptedStepCount() const noexcept;
    [[nodiscard]] double simulationTimeSeconds() const noexcept;
    [[nodiscard]] const SurfaceFlightCaseDiagnostics& diagnostics() const
        noexcept;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace simwing::fsi
