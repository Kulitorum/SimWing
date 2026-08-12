#include "surface_flight_case.h"

#include "scene_fluid_surface_transfer.h"
#include "structure_frame.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

[[nodiscard]] SceneStructureAssembly requireStructureAssembly(
    const Scene& scene,
    const SceneStructureSettings& settings) {
    SceneStructureAssembly assembly = assembleSceneStructure(scene, {}, settings);
    if (!assembly.ok()) {
        std::string message = "surface flight could not assemble scene structure";
        if (!assembly.diagnostics.empty()) {
            message += ": " + assembly.diagnostics.front().message;
        }
        throw std::invalid_argument(message);
    }
    return assembly;
}

[[nodiscard]] Scene prepareFlightScene(
    Scene scene,
    const double stiffnessScale) {
    if (!std::isfinite(stiffnessScale) || stiffnessScale < 1.0) {
        throw std::invalid_argument(
            "surface flight fabric stiffness scale must be finite and at least one");
    }
    for (FabricMaterial& material : scene.fabricMaterials) {
        material.warpStiffnessNewtonsPerMeter *= stiffnessScale;
        material.weftStiffnessNewtonsPerMeter *= stiffnessScale;
        material.shearStiffnessNewtonsPerMeter *= stiffnessScale;
        material.bendingStiffnessNewtonMeters *= stiffnessScale;
    }
    for (SeamMaterial& material : scene.seamMaterials) {
        material.axialStiffnessNewtons *= stiffnessScale;
    }
    return scene;
}

[[nodiscard]] SceneFluidSurfaceAssembly requireFluidSurface(
    const Scene& scene) {
    SceneFluidSurfaceAssembly assembly = assembleSceneFluidSurface(scene);
    if (!assembly.ok()) {
        std::string message = "surface flight could not assemble fluid surface";
        if (!assembly.diagnostics.empty()) {
            message += ": " + assembly.diagnostics.front().message;
        }
        throw std::invalid_argument(message);
    }
    return assembly;
}

[[nodiscard]] viewer::Vec3d toViewer(const StructureVector3& value) {
    return {value.x, value.y, value.z};
}

void addGlobalScalar(viewer::DiagnosticFrame& frame,
                     std::string name,
                     std::string unit,
                     const double value) {
    frame.scalarFields.push_back(
        {std::move(name), std::move(unit), viewer::FieldAssociation::Global,
         {value}});
}

void addGlobalVector(viewer::DiagnosticFrame& frame,
                     std::string name,
                     std::string unit,
                     const StructureVector3& value) {
    frame.vectorFields.push_back(
        {std::move(name), std::move(unit), viewer::FieldAssociation::Global,
         {toViewer(value)}});
}

[[nodiscard]] double totalAirMass(const SurfaceAerodynamicsState& state) {
    double total = 0.0;
    for (const SurfacePneumaticCellState& cell : state.cells) {
        total += cell.airMassKilograms;
    }
    return total;
}

} // namespace

SurfaceFlightCaseSettings::SurfaceFlightCaseSettings() {
    // A flight preview reports line residual but does not certify it. Keep the
    // numerical finite bound remote from expected motion; predictive cases use
    // the strict assembly default. The residual remains visible in every frame.
    sceneStructure.suspensionSolverIterations = 48;
    sceneStructure.suspensionMaximumLineResidualMeters = 0.10;
    structureStep.timeStepSeconds = 1.0 / 120.0;
    structureStep.substeps = 1;
    structureStep.constraintIterations = 48;
    structureStep.cableConstraintSweepPairs = 3;
    structureStep.gravityMetersPerSecondSquared = {0.0, 0.0, -9.80665};
    // The exact diagnostic mesh contains tiny elements and is not yet the
    // purpose-built interactive structural mesh. Preview-only damping bounds
    // those unresolved high-frequency modes; it is reported as model policy,
    // not a calibrated aerodynamic drag term.
    structureStep.velocityDampingPerSecond = 8.0;
    // Reproducibility is per explicit worker count. Twelve is a fixed product
    // setting rather than a machine-dependent hardware-concurrency guess.
    structureStep.workerThreads = 12;
    aerodynamics.timeStepSeconds = structureStep.timeStepSeconds;
}

struct SurfaceFlightCase::Implementation {
    Implementation(Scene input, const SurfaceFlightCaseSettings& requested)
        : settings(requested),
          scene(prepareFlightScene(
              std::move(input), settings.fabricStiffnessScale)),
          structureAssembly(requireStructureAssembly(
              scene, settings.sceneStructure)),
          surfaceAssembly(requireFluidSurface(scene)),
          structure(structureAssembly.definition),
          transfer(surfaceAssembly.definition,
                   structureAssembly.mappings,
                   structure),
          frameMapping(viewer::makeStructureFrameMapping(
              scene, structureAssembly, structure)),
          aerodynamics(surfaceAssembly.definition, settings.aerodynamics),
          aerodynamicState(aerodynamics.initialState(
              captureSceneFluidSurfaceState(
                  surfaceAssembly.definition,
                  structureAssembly.mappings,
                  structure))) {
        if (settings.structureStep.timeStepSeconds
                != settings.aerodynamics.timeStepSeconds) {
            throw std::invalid_argument(
                "surface flight structure and aerodynamic time steps must match");
        }
    }

    SurfaceFlightCaseSettings settings;
    Scene scene;
    SceneStructureAssembly structureAssembly;
    SceneFluidSurfaceAssembly surfaceAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;
    viewer::StructureFrameMapping frameMapping;
    SurfaceAerodynamicsModel aerodynamics;
    SurfaceAerodynamicsState aerodynamicState;
    SurfaceFlightCaseDiagnostics diagnostics;
};

SurfaceFlightCase::SurfaceFlightCase(
    Scene scene,
    const SurfaceFlightCaseSettings& settings)
    : implementation_(
          std::make_unique<Implementation>(std::move(scene), settings)) {}

SurfaceFlightCase::~SurfaceFlightCase() = default;
SurfaceFlightCase::SurfaceFlightCase(SurfaceFlightCase&&) noexcept = default;
SurfaceFlightCase& SurfaceFlightCase::operator=(SurfaceFlightCase&&) noexcept =
    default;

viewer::TraceHeader SurfaceFlightCase::traceHeader() const {
    return {implementation_->scene.metadata.designChecksum,
            surfaceFlightSolverId};
}

viewer::DiagnosticFrame SurfaceFlightCase::advance() {
    Implementation& impl = *implementation_;
    const SceneFluidSurfaceState surfaceState = captureSceneFluidSurfaceState(
        impl.surfaceAssembly.definition,
        impl.structureAssembly.mappings,
        impl.structure);
    const SurfaceAerodynamicsCandidate candidate = impl.aerodynamics.advance(
        impl.aerodynamicState, surfaceState);
    const ConservativeTransferResult transferResult = impl.transfer.evaluate(
        surfaceState, candidate.triangleTractions);
    if (!transferResult.diagnostics().finite) {
        throw std::runtime_error(
            "surface flight conservative transfer was non-finite");
    }
    std::vector<StructureVector3> appliedSurfaceForces(
        impl.structure.definition().nodes.size());
    if (impl.settings.massWeightedFlightPreviewLoads) {
        double surfaceMass = 0.0;
        for (std::size_t node = 0;
             node < impl.structureAssembly.mappings.nodeVertexIds.size();
             ++node) {
            surfaceMass += impl.structure.definition().nodes[node].massKg;
        }
        if (!(surfaceMass > 0.0)) {
            throw std::logic_error(
                "surface flight has no dynamic fabric mass for preview loads");
        }
        const StructureVector3 previewForce{
            candidate.diagnostics.aerodynamicForceNewtons.x
                + candidate.diagnostics.pressureForceNewtons.x,
            candidate.diagnostics.aerodynamicForceNewtons.y
                + candidate.diagnostics.pressureForceNewtons.y,
            candidate.diagnostics.aerodynamicForceNewtons.z
                + candidate.diagnostics.pressureForceNewtons.z};
        for (std::size_t node = 0;
             node < impl.structureAssembly.mappings.nodeVertexIds.size();
             ++node) {
            const double fraction = impl.structure.definition().nodes[node].massKg
                / surfaceMass;
            appliedSurfaceForces[node] = {
                previewForce.x * fraction,
                previewForce.y * fraction,
                previewForce.z * fraction};
        }
    } else {
        for (const CouplingNodeLoad& load : transferResult.nodeLoads()) {
            appliedSurfaceForces[load.structureNode] = load.forceNewtons;
        }
    }
    double maximumSurfaceLoadAcceleration = 0.0;
    for (std::size_t node = 0; node < appliedSurfaceForces.size(); ++node) {
        const StructureVector3& load = appliedSurfaceForces[node];
        const double mass = impl.structure.definition().nodes[node].massKg;
        const double force = std::sqrt(
            load.x * load.x + load.y * load.y + load.z * load.z);
        if (mass > 0.0) {
            maximumSurfaceLoadAcceleration = std::max(
                maximumSurfaceLoadAcceleration, force / mass);
        }
    }

    const StructureCheckpoint before = impl.structure.checkpoint();
    try {
        for (std::size_t node = 0; node < appliedSurfaceForces.size(); ++node) {
            const StructureVector3& load = appliedSurfaceForces[node];
            impl.structure.addExternalForce(node, load);
        }
        StructureStepSettings activeStep = impl.settings.structureStep;
        const auto targetWind =
            impl.settings.aerodynamics.targetWindMetersPerSecond;
        const double targetWindSpeed = std::sqrt(
            targetWind.x * targetWind.x + targetWind.y * targetWind.y
            + targetWind.z * targetWind.z);
        const auto activeWind = candidate.diagnostics.windMetersPerSecond;
        const double activeWindSpeed = std::sqrt(
            activeWind.x * activeWind.x + activeWind.y * activeWind.y
            + activeWind.z * activeWind.z);
        const double windActivation = targetWindSpeed > 1.0e-12
            ? std::clamp(activeWindSpeed / targetWindSpeed, 0.0, 1.0)
            : 1.0;
        // Lift starts quadratically with the wind ramp. Bring weight in one
        // power more gently so the authored inflated rest shape is not asked
        // to carry the pilot before aerodynamic support exists.
        const double weightActivation = windActivation * windActivation
            * windActivation;
        activeStep.gravityMetersPerSecondSquared.x *= weightActivation;
        activeStep.gravityMetersPerSecondSquared.y *= weightActivation;
        activeStep.gravityMetersPerSecondSquared.z *= weightActivation;
        const StructureDiagnostics structureDiagnostics =
            impl.structure.step(activeStep);
        if (!structureDiagnostics.finite) {
            throw std::runtime_error(
                "surface flight XPBD step was non-finite");
        }

        viewer::StructureFrameContext context;
        context.sceneChecksum = impl.scene.metadata.designChecksum;
        context.solverCommit = surfaceFlightSolverId;
        context.timeStepSeconds = impl.settings.structureStep.timeStepSeconds;
        context.couplingIteration = 0;
        context.couplingResiduals.tractionNewtons =
            transferResult.diagnostics().forceResidualNormNewtons;
        context.couplingResiduals.structure =
            structureDiagnostics.maximumMembraneResidual;
        context.couplingResiduals.interfacePowerWatts =
            transferResult.diagnostics().powerResidualWatts;
        context.conservation.fluidMassKilograms = totalAirMass(
            candidate.nextState);
        context.conservation.interfaceForceResidualNewtons = toViewer(
            transferResult.diagnostics().forceResidualNewtons);
        context.conservation.interfaceMomentResidualNewtonMetres = toViewer(
            transferResult.diagnostics().momentResidualNewtonMeters);
        context.conservation.interfacePowerResidualWatts =
            transferResult.diagnostics().powerResidualWatts;

        viewer::DiagnosticFrame frame = viewer::buildStructureFrame(
            impl.structure, impl.frameMapping, context);
        if (frame.triangles.size() != candidate.triangleTractions.size()) {
            throw std::logic_error(
                "surface flight triangle diagnostics lost topology alignment");
        }
        std::vector<double> pressureJump(frame.triangles.size());
        std::vector<double> externalTraction(frame.triangles.size());
        for (std::size_t index = 0; index < frame.triangles.size(); ++index) {
            if (frame.triangles[index].stableId
                    != candidate.triangleTractions[index].stableId) {
                throw std::logic_error(
                    "surface flight triangle diagnostics changed stable order");
            }
            pressureJump[index] = candidate.trianglePressureJumpPascals[index];
            externalTraction[index] =
                candidate.triangleExternalTractionPascals[index];
        }
        frame.scalarFields.insert(
            frame.scalarFields.begin(),
            {"surface_aero.external_traction", "Pa",
             viewer::FieldAssociation::Triangle,
             std::move(externalTraction)});
        frame.scalarFields.insert(
            frame.scalarFields.begin(),
            {"surface_aero.pressure_jump", "Pa",
             viewer::FieldAssociation::Triangle,
             std::move(pressureJump)});

        std::vector<viewer::Vec3d> appliedForces(frame.vertices.size());
        for (std::size_t node = 0;
             node < appliedSurfaceForces.size() && node < appliedForces.size();
             ++node) {
            appliedForces[node] = toViewer(appliedSurfaceForces[node]);
        }
        frame.vectorFields.insert(
            frame.vectorFields.begin(),
            {"surface_aero.applied_force", "N",
             viewer::FieldAssociation::Vertex,
             std::move(appliedForces)});

        const auto& aero = candidate.diagnostics;
        addGlobalScalar(frame, "surface_aero.dynamic_pressure", "Pa",
                        aero.dynamicPressurePascals);
        addGlobalScalar(frame, "surface_aero.planform_area", "m^2",
                        aero.planformAreaSquareMeters);
        addGlobalScalar(frame, "surface_aero.angle_of_attack", "rad",
                        aero.angleOfAttackRadians);
        addGlobalScalar(frame, "surface_aero.lift_coefficient", "1",
                        aero.liftCoefficient);
        addGlobalScalar(frame, "surface_aero.drag_coefficient", "1",
                        aero.dragCoefficient);
        addGlobalScalar(frame, "surface_aero.lift", "N", aero.liftNewtons);
        addGlobalScalar(frame, "surface_aero.drag", "N", aero.dragNewtons);
        addGlobalScalar(frame, "surface_aero.minimum_cell_pressure", "Pa",
                        aero.minimumCellGaugePressurePascals);
        addGlobalScalar(frame, "surface_aero.maximum_cell_pressure", "Pa",
                        aero.maximumCellGaugePressurePascals);
        addGlobalScalar(frame, "surface_aero.minimum_cell_volume_fraction", "1",
                        aero.minimumCellVolumeFraction);
        addGlobalScalar(frame, "surface_aero.opening_mass_flow", "kg/s",
                        aero.totalOpeningMassFlowKilogramsPerSecond);
        addGlobalScalar(frame, "surface_aero.maximum_load_acceleration", "m/s^2",
                        maximumSurfaceLoadAcceleration);
        addGlobalScalar(frame, "surface_aero.mass_weighted_preview_load", "1",
                        impl.settings.massWeightedFlightPreviewLoads ? 1.0 : 0.0);
        addGlobalVector(frame, "surface_aero.wind", "m/s",
                        aero.windMetersPerSecond);
        addGlobalVector(frame, "surface_aero.relative_wind", "m/s",
                        aero.relativeWindMetersPerSecond);
        addGlobalVector(frame, "surface_aero.polar_force", "N",
                        aero.aerodynamicForceNewtons);
        addGlobalVector(frame, "surface_aero.pressure_force", "N",
                        aero.pressureForceNewtons);

        viewer::ProtocolError frameError;
        if (!viewer::validateFrame(frame, &frameError)) {
            throw std::runtime_error(
                "surface flight diagnostic frame is invalid: "
                + frameError.message);
        }

        impl.aerodynamicState = candidate.nextState;
        impl.diagnostics.aerodynamics = candidate.diagnostics;
        impl.diagnostics.transfer = transferResult.diagnostics();
        impl.diagnostics.structure = structureDiagnostics;
        impl.diagnostics.maximumSurfaceLoadAccelerationMetersPerSecondSquared =
            maximumSurfaceLoadAcceleration;
        impl.diagnostics.finite = candidate.diagnostics.finite
            && transferResult.diagnostics().finite
            && structureDiagnostics.finite;
        return frame;
    } catch (...) {
        impl.structure.restore(before);
        throw;
    }
}

const Structure& SurfaceFlightCase::structure() const noexcept {
    return implementation_->structure;
}

const StructureStepSettings& SurfaceFlightCase::stepSettings() const noexcept {
    return implementation_->settings.structureStep;
}

std::uint64_t SurfaceFlightCase::acceptedStepCount() const noexcept {
    return implementation_->structure.checkpoint().acceptedStepCount;
}

double SurfaceFlightCase::simulationTimeSeconds() const noexcept {
    return implementation_->structure.checkpoint().simulationTimeSeconds;
}

const SurfaceFlightCaseDiagnostics& SurfaceFlightCase::diagnostics() const
    noexcept {
    return implementation_->diagnostics;
}

} // namespace simwing::fsi
