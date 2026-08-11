#pragma once

#include "structure.h"

#include <cstddef>

namespace simwing::fsi {

// Bounds for a transactional zero-load probe. Linear and angular velocity
// have separate limits because an XPBD step can return near its initial
// position while retaining momentum that would move the following step.
struct StructureRestAuditLimits {
    double maximumLinearDisplacementMeters = 1.0e-8;
    double maximumLinearSpeedMetersPerSecond = 1.0e-6;
    double maximumAngularDisplacementRadians = 1.0e-8;
    double maximumAngularSpeedRadiansPerSecond = 1.0e-6;
};

struct StructureRestAuditDiagnostics {
    bool stationary = false;
    std::size_t nodeCount = 0;
    std::size_t maximumDisplacementNode = 0;
    double maximumNodeDisplacementMeters = 0.0;
    double rmsNodeDisplacementMeters = 0.0;
    double maximumNodeSpeedMetersPerSecond = 0.0;
    double payloadDisplacementMeters = 0.0;
    double payloadRotationRadians = 0.0;
    double maximumHarnessDisplacementMeters = 0.0;
    double maximumPayloadOrHarnessSpeedMetersPerSecond = 0.0;
    double payloadAngularSpeedRadiansPerSecond = 0.0;
    std::size_t maximumInitialMembraneEdgeMismatchIndex = 0;
    double maximumInitialMembraneEdgeMismatchMeters = 0.0;
    std::uint64_t maximumInitialSuspensionExtensionSegmentStableId = 0;
    double maximumInitialSuspensionExtensionMeters = 0.0;
    double maximumSuspensionResidualMeters = 0.0;
};

// Probes whether the committed Structure state is stationary without gravity
// or external loads. The caller's complete checkpoint is restored on both
// success and failure. Nonzero pending loads, initial velocity, or a nonzero
// gravity/damping reference are rejected because they are not a rest audit.
[[nodiscard]] StructureRestAuditDiagnostics auditStructureRestState(
    Structure& structure,
    const StructureStepSettings& stepSettings,
    const StructureRestAuditLimits& limits = {});

} // namespace simwing::fsi
