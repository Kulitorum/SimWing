#include "canonical_case.h"

#include <cmath>
#include <stdexcept>

namespace simwing::fsi {
namespace {

StructureDefinition makeDefinition() {
    StructureDefinition definition;
    definition.nodes = {
        StructureNodeDefinition{{0.0, 0.0, 0.0}, 0.0, true},
        StructureNodeDefinition{{1.0, 0.0, 0.0}, 0.25, false},
        StructureNodeDefinition{{0.0, 1.0, 0.0}, 0.25, false},
    };
    definition.triangles.push_back({{0, 1, 2}});
    definition.constraints = {
        {StructureConstraintKind::Distance, 0, 1, 1.0, 1.0e-5},
        {StructureConstraintKind::Cable, 0, 2, 1.05, 1.0e-4},
        {StructureConstraintKind::SuspensionTie, 1, 2, std::sqrt(2.0),
         1.0e-5},
    };

    StructureMembraneDefinition membrane;
    membrane.triangle = 0;
    membrane.materialCoordinates = {
        StructureVector2{0.0, 0.0},
        StructureVector2{1.0, 0.0},
        StructureVector2{0.0, 1.0},
    };
    membrane.material.warpStiffnessNewtonsPerMeter = 800.0;
    membrane.material.weftStiffnessNewtonsPerMeter = 500.0;
    membrane.material.couplingStiffnessNewtonsPerMeter = 100.0;
    membrane.material.shearStiffnessNewtonsPerMeter = 180.0;
    membrane.material.dampingSeconds = 0.01;
    membrane.material.compressionStiffnessRatio = 1.0;
    membrane.role = StructureMaterialRole::Bulk;
    definition.membranes.push_back(membrane);
    return definition;
}

viewer::StructureFrameMappingDefinition makeFrameMapping() {
    viewer::StructureFrameMappingDefinition mapping;
    mapping.vertexStableIds = {1001, 1002, 1003};
    mapping.triangles = {{2001, 1, 2}};
    mapping.lines = {
        {3001, static_cast<std::uint32_t>(StructureConstraintKind::Distance)},
        {3002, static_cast<std::uint32_t>(StructureConstraintKind::Cable)},
        {3003,
         static_cast<std::uint32_t>(StructureConstraintKind::SuspensionTie)},
    };
    return mapping;
}

StructureStepSettings makeStepSettings() {
    StructureStepSettings settings;
    settings.timeStepSeconds = 1.0 / 120.0;
    settings.substeps = 2;
    settings.constraintIterations = 12;
    settings.cableConstraintSweepPairs = 0;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.05;
    settings.workerThreads = 0;
    return settings;
}

} // namespace

CanonicalStructuralCase::CanonicalStructuralCase()
    : structure_(makeDefinition()),
      frameMapping_(structure_, makeFrameMapping()),
      stepSettings_(makeStepSettings()) {}

viewer::TraceHeader CanonicalStructuralCase::traceHeader() const {
    return {canonicalCaseChecksum, canonicalCaseSolverId};
}

viewer::DiagnosticFrame CanonicalStructuralCase::advance() {
    // Constant, source-known loads make this a deterministic structural
    // integration case. They enter through public nodal load paths only.
    structure_.addExternalForce(1, {0.40, 0.05, -0.20});
    structure_.addExternalForce(2, {0.20, -0.05, -0.10});
    const StructureDiagnostics diagnostics = structure_.step(stepSettings_);
    if (!diagnostics.finite) {
        throw std::runtime_error(
            "canonical structural step produced non-finite diagnostics");
    }

    viewer::StructureFrameContext context;
    context.sceneChecksum = canonicalCaseChecksum;
    context.solverCommit = canonicalCaseSolverId;
    context.timeStepSeconds = stepSettings_.timeStepSeconds;
    context.couplingIteration = 0;
    // Keep dimensioned constraint errors in their named structure fields.
    // The generic coupling residual is dimensionless, so only the membrane
    // solver's normalized residual belongs here.
    context.couplingResiduals.structure =
        diagnostics.maximumMembraneResidual;
    return viewer::buildStructureFrame(structure_, frameMapping_, context);
}

const Structure& CanonicalStructuralCase::structure() const noexcept {
    return structure_;
}

const StructureStepSettings& CanonicalStructuralCase::stepSettings() const
    noexcept {
    return stepSettings_;
}

} // namespace simwing::fsi
